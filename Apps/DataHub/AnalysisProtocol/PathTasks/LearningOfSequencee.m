function LearningOfSequencee(Data, StartTrial, EndTrial,ShowFigureFlag,OutputFlag);

ECode;
DataPath= Data.Path;
FileName=Data.FileName;
TaskType=Data.TaskType;
TaskCode=Data.TaskCode;



%for compatibility of the old events 

codecheck = FindOutTime(Data.EventChannel,Data.EventTimeChannel,FIRSTTGT);
codecheck=codecheck(~isnan(codecheck ));

if isempty(codecheck)
    disp("Old events, auto transform segment events");
    
    CheckCode=[TGTCD;TGTOFF];
    AbnormalTrials=AbnormalTrialCheck(Data,CheckCode);
    
    OldCode = 1;
    SegEvent=[SHOWFIXCD,TGTCD];%Event Flow of interest
    
   
   

else
    
%{\
%Abnormal check for equality
CheckCode=[[FIRSTTGT:FIFTHTGT];[TGTOFF,TGTOFF,TGTOFF,TGTOFF,TGTOFF]];
AbnormalTrials=AbnormalTrialCheck(Data,CheckCode);
%}
 OldCode =0;
 
 SegEvent=[SHOWFIXCD,FIRSTTGT:FIFTHTGT,FINALFIX];%Event Flow of interest
 
 
 
 
end

 PlotInterval = [-500,1000];
 PlotEyeTrace = 1;
 IntervalForFP =[100,500];
 IntervalForTgt=[-600,0];
 PlotError = 0;
 
 ReChunk=0;%Temporally put here;

 



%SegEventEnd=[SHOWFIXCD,FIXOFF,TGTOFF,RWDOFFCD ,NaN,EXTRAFP];%Event end
%SegEventEnd=[TGTOFF];%Event end
SegEventEnd=NaN*ones(1,length(SegEvent));
SegEventStr={'FixOn','1TGT','2TGT','3TGT','4TGT','5TGT','LastFix'};
MarkerSizeEvent = 20*ones(1,length(SegEventStr));
MarkerSizeEvent([1,end]) = 10;
MarkerColorEvent = {'#44AA99','#FE6100','#FE6100','#FE6100','#FE6100','#FE6100','#44AA99'};
DistratorColorEvent = {'#44AA99','#648FFF','#648FFF','#648FFF','#648FFF','#648FFF','#44AA99'};

PlotStartOffset=1000*ones(1,length(SegEvent));%in ms
PlotStopOffset=NaN*ones(1,length(SegEvent));%in ms

%{
FRISTTGT = 171;
SECONDTGT = 172;
THIRDTGT = 173;
FOURTHTGT = 174;

FINALFIX = 179;
FAILFINALFIX = 180;
%}


%Parameter to define the sliding window to measure firing rates
StepSize=5;%in ms
BinWidth=15;%in ms
%Color Definition for each scenes

Color=[lines(6);prism(6);pink(6)];
%Color(2,:)=copper(1)+0.5;
Color=[Color;hsv(10)];
Color(18,:)=[0.1,0.457,0.194];
Color=[Color;Color*0.5];
Color=[Color;Color*0.8];
%Color Definition for each fractals
%{
Color=[lines(6);prism(6) ];
Color(2,:)=copper(1)+0.5;
Color=[Color;Color*0.5];
Color=[Color;Color*0.8];
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Unpack all the data base and task settings



Selection=StartTrial:EndTrial;
%{
%Abnormal check %Especially For Robin081120FORAGEXPLO4.rst
CheckCode=SCENE_ON;%This code should appear once in each trial
AbnormalTrials=AbnormalTrialCheck(Data,CheckCode);
%}
Selection=Selection(~ismember(Selection,AbnormalTrials));

TrialNum=length(Selection);

%FixationPointIndex=Data.FixPointIndex(Selection,:);
%FixLoc=Data.FixLoc(Selection,:,:);


TargetIndex=Data.TargetIndex(Selection,2:end,1);
TargetName=Data.TargetIndex(Selection,2:end,2);
%TargetLoc=Data.TargetLoc(Selection,:,:);

TargetNum=Data.TargetNum(Selection,2:end);
TargetLoc=Data.TargetLoc(Selection,2:end,:);
TargetPolarLoc= Data.TargetPolarLoc(Selection,2:end,:);

%SuccessfulMat = TransformResult(TargetNum);%The last one is not accurate

%CorrectTrialIndex = sum(SuccessfulMat,2)==size(TargetNum,2);




%Sequence information from the library
ParaLib=Data.ParaLib;
if ismember('HyperSetNO,ManuelFlag,LocSeqCode',keys(ParaLib))
    Seq_Info=ParaLib('HyperSetNO,ManuelFlag,LocSeqCode');
    Seq_Info=Seq_Info(Selection,:);
    SeqNo=Seq_Info(:,2);
end

if ismember('LocationSeqTgt,LocationSeqDis',keys(ParaLib))
    Seq_InfoLoc=ParaLib('LocationSeqTgt,LocationSeqDis');
    Seq_InfoLoc=ReOrganizeData(Seq_InfoLoc(Selection,:));
end


%Reconstruct sequence code from the target location
SeqType = floor(SeqCode(squeeze(TargetPolarLoc(:,:,1))));
SeqAll = [SeqNo,SeqType(:,end)];
SeqAllUnique = unique(SeqAll,'rows','stable');
SeqAllUnique =uint64(SeqAllUnique(~isnan(SeqAllUnique(:,2)),:));

%For convenience, represent sequence with NO
SeqNoUnique = SeqAllUnique(:,1);
SeqCodeUnique = SeqAllUnique(:,2);

UniqueSeq = SeqNoUnique;%Use sequence NO for the seq loops
SeqCodeUse = SeqNo;



%Event channel

EventChannel=Data.EventChannel(Selection,:);
EventTimeChannel=Data.EventTimeChannel(Selection,:);
EventBin=Data.EventBin(Selection,:);
%GetEyeData
EyeData=Data.EyeDataRaw;
EyeChannel_X_L=EyeData(Selection,1);
EyeChannel_Y_L=EyeData(Selection,2);
%EyeChannel_X_R=EyeData(Selection,3);
%EyeChannel_Y_R=EyeData(Selection,4);




EyeBinWidth=Data.EyeBinWidth;

%Recalibrate eyes according to the extra FP
EyeCalInfo=Data.EyeCalInfo;%Each row: X_Gain,X_Bias; Y_Gain; Y_bias; X=(raw+X_bias)*X_Gain
EyeCalInfo_L=EyeCalInfo(1,:);


%Eye_Eccentricity=cellfun(@(x,y) sqrt(x.^2+y.^2),EyeChannel_X_L,EyeChannel_Y_L,'uniform',0);

%Now load the spike channel

SpikeChannelWhole=Data.SpikeChannel;

for index = 1:length(SpikeChannelWhole)
    SpikeChannel = SpikeChannelWhole(index);
    SpikeTime=squeeze(Data.SpikeTimeData(index,:,Selection))';
    
    
    


SpikeBin=squeeze(Data.SpikeBinData(:,:,Selection))';

%
%It took a long time to chunk data according to different event marker.
%To save time, save the results to the data folder
%OutputFileName=sprintf('%s_%s_N%s_C%s.mat',MonkeyName,string(RecordDateOriginal),string(NeuronNum),string(ChannelNum));

path=DataPath;

cd(path);


%Get File Information
num = regexp(FileName, '\d+', 'match');
%RecordDateOriginal=cell2mat(num(1));

NeuronNum=cell2mat(num(2));
SpikeChannel=0;
ChannelNum=SpikeChannel;

SpikeFileName=strcat(FileName(1:end-4),'_Seq_Spike_N',string(NeuronNum),'_C',string(ChannelNum),'.mat');

if exist(SpikeFileName) & ~ReChunk
    load(SpikeFileName);
    disp('preanalysized chuncked data aligned on tgt-on has been loaded to the workspace') ;
else
   disp('Chuncking the data, please wait.') 
   

%Segregate spike data from fractal onset to fractal offset

SegmentData=EventSegTool(SegEvent,SegEventEnd,PlotStartOffset,PlotStopOffset,EventChannel,EventTimeChannel,SpikeTime,StepSize,BinWidth);

save(SpikeFileName,'SegmentData');
disp('Seqed spikes sucessfully exported'); 
load(SpikeFileName);
disp('preanalysized chuncked data aligned on tgt-on has been loaded to the workspace') ;



    
end

%Retrive data from the field
%Total average for each event
MeanFR=SegmentData.MeanFR; %[Mean,sem]%Averall mean and sem
MeanFRRaw=SegmentData.MeanFRRaw;%=AverageFiringRate_Event;

%Relative time bin for each event
RelativeTimeMarker=SegmentData.RelativeTimeMarker;%=EventTimeMarker;
%Absolute time marker for each event
AbosoluteTimeMarker=SegmentData.AbsoluteTimeMarker;



%TgtOnTime_test=FindOutTime(EventChannel,EventTimeChannel,FRISTTGT);

%FPOnTime_test=FindOutTime(EventChannel,EventTimeChannel,SHOWFIXCD);
%tmp=AbosoluteTimeMarker{2}(:,1);


%Transform the data if it is the old code

if OldCode
    MeanPSTH=SegmentData.MeanPSTHEach;%=[PSTH_Event_Time',PSTH_Event'];
    MeanPSTH_Time=MeanPSTH(:,1);    
    MeanPSTH_FR=MeanPSTH(:,2);
    
    MeanPSTH_Time ={ MeanPSTH_Time{1}{1},MeanPSTH_Time{2}{1:end},MeanPSTH_Time{1}{2}};
    MeanPSTH_FR = { MeanPSTH_FR{1}{1},MeanPSTH_FR{2}{1:end},MeanPSTH_FR{1}{2}};
    
    EndCheck =ReproduceFromEvent(EventChannel,SHOWFIXCD);
    EndCheck = double((cellfun(@numel,EndCheck)==2));
    EndCheck(EndCheck==0)=NaN;
    
    
    
else
    %PSTH for each event
    MeanPSTH=SegmentData.MeanPSTH;%=[PSTH_Event_Time',PSTH_Event'];
    MeanPSTH_Time=MeanPSTH(:,1);    
    MeanPSTH_FR=MeanPSTH(:,2);

    EndCheck =ReproduceFromEvent(EventChannel,FINALFIX);
    EndCheck = double(~cellfun(@isempty,EndCheck));
    EndCheck(EndCheck==0)=NaN;

end


RewardTime=FindOutTime(EventChannel,EventTimeChannel,RWDCD);
TgtOffTime=FindOutTime(EventChannel,EventTimeChannel,TGTOFF);
FixOffTime=FindOutTime(EventChannel,EventTimeChannel,FIXOFF);





%TimeMarker=SegmentData.RelativeTimeMarker{1}; %[Mean,sem]%Averall mean and sem



%Trial use: here only use successful trials
%TrialsForUse =

%Correct Rate
TotalTrialNum = size(EventChannel,1);
CorrectTrialIndicator = cell2mat(ReproduceFromEvent(EventChannel,[ERROR,CORRECTCD]));
CorrectTrialNum = sum(CorrectTrialIndicator ==CORRECTCD);
ErrorTrialNum = sum(CorrectTrialIndicator ==ERROR);


CorrectTrialMat = repmat([CorrectTrialIndicator==CORRECTCD]',1,7);


SuccessfulMat=ReorganizeEye( ReproduceFromEvent(EventChannel,[MATCHOBJ,NONMATCHOBJ]));
SuccessfulMat(SuccessfulMat==MATCHOBJ)=1;
SuccessfulMat(SuccessfulMat==NONMATCHOBJ)=0;
SuccessfulMatEach = [ones(size(SuccessfulMat,1),1),SuccessfulMat,EndCheck'];

SuccessfulMatFull = SuccessfulMatEach ;

%Align Eye Data with events
EyeForAlign.Bin =EyeBinWidth ;
EyeForAlign.X =EyeChannel_X_L;
EyeForAlign.Y =EyeChannel_Y_L;
EyeForAlign.MaxTime = cellfun(@numel,EyeChannel_X_L);

Eye_Aligned = AlignEyeDataWithEvents(SegEvent,SegEventEnd,PlotStartOffset,PlotStopOffset,EventChannel,EventTimeChannel,EyeForAlign);

Eye_Aligned_X=Eye_Aligned.X;
Eye_Aligned_Y=Eye_Aligned.Y;
Eye_Aligned_Time = Eye_Aligned.Time;

if OldCode
   
    
    Eye_Aligned_Time ={ Eye_Aligned_Time{1,1},Eye_Aligned_Time{2,1:end},Eye_Aligned_Time{1,2}};
    Eye_Aligned_X = { Eye_Aligned_X{1,1},Eye_Aligned_X{2,1:end},Eye_Aligned_X{1,2}};
    Eye_Aligned_Y = { Eye_Aligned_Y{1,1},Eye_Aligned_Y{2,1:end},Eye_Aligned_Y{1,2}};
    
   
    
    
    
end


RewardTimeSeq=[];
LastTgtOffTimeSeq=[];




%Loop over the sequence 
for i = 1:length(UniqueSeq)
%for i = 1:1
    sel = SeqCodeUse == UniqueSeq(i);
    %sel = SeqCodeUse == 1006;
    SucIndex = SuccessfulMatFull(sel,:);

    
    
    TargetLoc_curr_X = nanmean(TargetLoc(sel,:,1),1);
    TargetLoc_curr_Y = nanmean(TargetLoc(sel,:,2),1);
    curr_X(i,:) = [0,0,TargetLoc_curr_X,0,0];
    curr_Y(i,:) = [0,0,TargetLoc_curr_Y,0,0];
    
    max_fr_seq(i) = 0; 
    
    TotalTrialsSeq(i) = sum(sel);
    CorrectTrialSeq(i) = sum(CorrectTrialIndicator(sel) ==CORRECTCD);
    for j = 1:length(SegEventStr)
        %Neural responses
        
            %SuccessfulMatcurr=SuccessfulMatFull(:,j);
            
            %%%%%%%%%%%  Select all correct hypersets   %%%%%%%%%%
            %{
            SuccessfulMatcurr=SuccessfulMatFull(:,j);
            
            select = ~isnan(SuccessfulMatFull(:,j));
            
            correctSelector = (CorrectTrialIndicator==CORRECTCD)' & sel;%correct trial index for the current scene
            
            select_idx = correctSelector(select == 1);
            
            error_curr = (CorrectTrialIndicator==ERROR)' & sel;  
            select_error = error_curr(select==1);
            
            select_eye_correct = correctSelector;
            %}
            
            %{
            PSTH_Time_curr = nanmean(MeanPSTH_Time{j}(sel(SuccessfulMatcurr == 1),:),1);
            
             
            PSTH_curr = MeanPSTH_FR{j}(sel(SuccessfulMatcurr == 1) ,:);
            %}
            %%%%%%%%%%%% If select all correct sets in the hyperset %%%%%%%
            SuccessfulMatcurr=SuccessfulMatFull(:,j);
            select = ~isnan(SuccessfulMatFull(:,j));
            correctSelector = SuccessfulMatcurr==1  & sel;%correct trial index for the current scene
            
           % select_idx = correctSelector(select == 1 &  correctSelector==1);%Modify to commondate for the only correct selection
           select_idx = correctSelector(select == 1 );
            error_curr =SuccessfulMatcurr==0 & sel;  
           % select_error = error_curr(select==1 & error_curr==1);%Modify to commondate for the only correct selection
            select_error = error_curr(select==1);
            select_eye_correct = correctSelector ;
            
            
            
            
            
            
            PSTH_Time_curr = nanmean(MeanPSTH_Time{j}(select_idx,:),1);
           
  
            
             
            PSTH_curr = MeanPSTH_FR{j}(select_idx ,:);
            
            
            %A{j}=PSTH_curr;
            
             
             
             PSTH_curr_error = MeanPSTH_FR{j}(select_error ,:);
             PSTH_Time_curr_error =  nanmean(MeanPSTH_Time{j}(select_error ,:),1);
        
        
        PSTH_Mean_Time{i,j} = PSTH_Time_curr;
        if size(PSTH_curr,1)<3
            PSTH_curr = NaN*PSTH_curr;
        end
        PSTH_Mean_Event{i,j} = nanmean(PSTH_curr,1);
        PSTH_SEM_Event{i,j} = nanstd(PSTH_curr,[],1)./sqrt(sum(~isnan(PSTH_curr),1));
        
        if size(PSTH_curr_error,1)<3
            PSTH_curr_error = NaN*PSTH_curr_error;
        end
        PSTH_Mean_Time_error{i,j} = PSTH_Time_curr_error;
        PSTH_Mean_Event_error{i,j} = nanmean(PSTH_curr_error ,1);
        PSTH_SEM_Event_error{i,j} = nanstd(PSTH_curr_error ,[],1)./sqrt(sum(~isnan(PSTH_curr_error ),1));
        
        
        sel_tmp =PSTH_Mean_Time{i,j}>=PlotInterval(1) & PSTH_Mean_Time{i,j}<=PlotInterval(2);
        if max(PSTH_Mean_Event{i,j}(sel_tmp)+PSTH_SEM_Event{i,j}(sel_tmp)) > max_fr_seq(i) 
            max_fr_seq(i) = max(PSTH_Mean_Event{i,j}(sel_tmp)+PSTH_SEM_Event{i,j}(sel_tmp));
            
        end
        %Eye trace
        
        
        
        
            
        curr_Eye_x{i,j} = Eye_Aligned_X{j}(select_eye_correct,:);
        curr_Eye_y{i,j} = Eye_Aligned_Y{j}(select_eye_correct,:);
        
        curr_Eye_x_error{i,j} = Eye_Aligned_X{j}(sel&(SuccessfulMatcurr == 0),:);
        curr_Eye_y_error{i,j} = Eye_Aligned_Y{j}(sel&(SuccessfulMatcurr == 0),:);
       
        
        curr_Eye_time{i,j}=Eye_Aligned_Time{j}(select_eye_correct,:);
        
        
        if OldCode
            if j==1
                MarkerTime = FindOutTime(EventChannel,EventTimeChannel,SegEvent(1));
                MarkerTime=MarkerTime(:,1);
            elseif j== length(SegEventStr)
                MarkerTime = FindOutTime(EventChannel,EventTimeChannel,SegEvent(1));
                MarkerTime=MarkerTime(:,2);
            else
                MarkerTime = FindOutTime(EventChannel,EventTimeChannel,SegEvent(2));
                MarkerTime=MarkerTime(:,j-1);
   
            end
             
            
        else
            MarkerTime = FindOutTime(EventChannel,EventTimeChannel,SegEvent(j));
        end
       
        
        
        if j >1
         %Reward Time shown
         
         if j-1 <size(RewardTime,2)
             
        
         
            RewardTime_tmp=RewardTime(:,j-1)-MarkerTime;

            RewardTimeSeq(i,j)=nanmean(RewardTime_tmp(sel&(SuccessfulMatcurr == 1)));
         else
             RewardTimeSeq(i,j)=NaN;
             
         end
         
         if j >2
         TgtOffTime_tmp=TgtOffTime(:,j-2)-MarkerTime;
         
         LastTgtOffTimeSeq(i,j)=nanmean(TgtOffTime_tmp(sel&(SuccessfulMatcurr == 1)));
         end
         if j==length(SegEventStr)
             LastFixOffTime_tmp=FixOffTime(:,2)-MarkerTime;
             LastFixOffTimeSeq(i,j)=nanmean(FixOffTime_tmp(sel&(SuccessfulMatcurr == 1)));
         
         end
         
       %{  
       Marker_tmp = AbosoluteTimeMarker{j}(sel(SuccessfulMatcurr == 1),:);
       
        %Time of reward on and tgt off
        
       RewardTime_tmp =RewardTime(sel&(SuccessfulMatcurr == 1),j-1)-Marker_tmp(:,1);
        if j <6
            TgtOffTime_tmp=TgtOffTime(sel&(SuccessfulMatcurr == 1),j-1)-Marker_tmp(:,1);
        end
         %}
       else
            %Add fix off time
            
            FixOffTime_tmp=FixOffTime(:,1)-MarkerTime;
            
            FixOffTimeSeq(i,j)=nanmean(FixOffTime_tmp(sel&(SuccessfulMatcurr == 1)));
       end
        
        
    end
    
    
    
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% For figure plotting%%%%%%%%%%%%%%%%%
if ShowFigureFlag

    Plot_First = 20; %First 10 trials 
    Plot_Last = 20;%Last 10 trials
    
  %Loop over the sequence type
  StartNO = 130;
 if  PlotEyeTrace 
     Fig1Size = [50,100, 1200,800];
     ylevel2 = 0.15;
     ysize = 0.3;
 else
      Fig1Size = [50,100, 1200,400];
      ylevel2 = 0.3;
      ysize = 0.4;
 end
  %PlotEyeInterval=[-200,0];
  
      
  for i =1 :length(UniqueSeq)
      if min(UniqueSeq)>=1000
            figtitlestr{i}=sprintf("Scene Path Task, HyperSet:%d,TotalTrial:%d,Correct:%d",UniqueSeq(i),TotalTrialsSeq(i),CorrectTrialSeq(i));
            figtitlestr_save{i} = sprintf("Scene_S%d_CT%d",UniqueSeq(i),CorrectTrialSeq(i));
      else
          figtitlestr{i}=sprintf("Fractal Path Task, HyperSet:%d,TotalTrial:%d,Correct:%d",UniqueSeq(i),TotalTrialsSeq(i),CorrectTrialSeq(i));
           figtitlestr_save{i} = sprintf("Frac_S%d_CT%d",UniqueSeq(i),CorrectTrialSeq(i));
          
      end
      Fig1Size(2)=Fig1Size(2)+(i-1)*100;
      fig{i}=PrepareFigure(StartNO+i-1,'w',Fig1Size,'Name',figtitlestr{i});
      
      
      for j = 1:length(SegEventStr)
          %First plot the target location
          %title(SegEventStr{j});
          axes('Position',[(j-1)*0.136,0.75,0.2,0.2]);
          rc_w = 40;rc_h = 40;
          
          ctr_x =curr_X(i,2*(j-1)+1:2*j);
          ctr_y =curr_Y(i,2*(j-1)+1:2*j);
          
          rectangle('Position',[0 0 rc_w rc_h],'Facecolor','#dcded5')
          hold on
          
          plot(ctr_x(1)+rc_w/2,ctr_y(1)+rc_h/2,'o','MarkerSize',MarkerSizeEvent(j),'MarkerFaceColor',MarkerColorEvent{j},'color',MarkerColorEvent{j});
          plot(ctr_x(2)+rc_w/2,ctr_y(2)+rc_h/2,'o','MarkerSize',MarkerSizeEvent(j),'MarkerFaceColor',DistratorColorEvent{j},'color',DistratorColorEvent{j});
          axis equal;
          axis off;
          %Add a vector poting from the last location to the current
          %location
          if j>1
              last_x = curr_X(i,2*(j-2)+1);
              last_y = curr_Y(i,2*(j-2)+1);
              curr_x = ctr_x(1);
              curr_y = ctr_y(1);
              diff_x = curr_x - last_x;
              diff_y = curr_y - last_y;
              
              
              quiver(last_x+rc_w/2,last_y+rc_h/2,diff_x,diff_y,'LineWidth',3,'MaxHeadSize',5,'Color','r');
              
          end
        if PlotEyeTrace 
          %Eye traces
          axes('Position',[(j-1)*0.136,0.5,0.2,0.2]);
          
          if j==1|j==length(SegEventStr)
              PlotEyeInterval=IntervalForFP;
              
          else
               PlotEyeInterval=IntervalForTgt;
               
              
          end
          
          %time_tmp=-PlotStartOffset:-PlotStartOffset+size(curr_Eye_x{i,j},2);
          
          time_tmp = curr_Eye_time{i,j}(1,:);
          
          sel = time_tmp>=PlotEyeInterval(1) & time_tmp<=PlotEyeInterval(2); 
           rectangle('Position',[0 0 rc_w rc_h],'Facecolor','#dcded5');
           
        
             %{
              if ~isempty(curr_Eye_x_error{i,j})
                plot((curr_Eye_x_error{i,j}(:,sel)+rc_w/2)',(curr_Eye_y_error{i,j}(:,sel)+rc_h/2)','-r','color','#211902','LineWidth',3);
              end
             %}
             hold on
             if ~isempty(curr_Eye_x{i,j})
                 tmp_x = curr_Eye_x{i,j}(:,sel);
                 num_x = size(tmp_x);
                 num_x = num_x(1);

                plot((curr_Eye_x{i,j}(1:min(Plot_First,num_x),sel)+rc_w/2)',(curr_Eye_y{i,j}(1:min(Plot_First,num_x),sel)+rc_h/2)','-r','color',MarkerColorEvent{j},'LineWidth',2);
               
             end
          xlim([-20+rc_w/2,20+rc_w/2]);
          ylim([-20+rc_h/2,20+rc_h/2]);
        %  box on;
          
          axis equal;
          axis off;
          
        end %iF PLOT EYE TRACE
          
          
          
          %Last few trials
          axes('Position',[0.04+(j-1)*0.13,ylevel2,0.2,0.2]);
           if j==1|j==length(SegEventStr)
              PlotEyeInterval=IntervalForFP;
              
          else
               PlotEyeInterval=IntervalForTgt;
               
              
          end
          
          %time_tmp=-PlotStartOffset:-PlotStartOffset+size(curr_Eye_x{i,j},2);
          
          time_tmp = curr_Eye_time{i,j}(1,:);
          
          sel = time_tmp>=PlotEyeInterval(1) & time_tmp<=PlotEyeInterval(2); 
           rectangle('Position',[0 0 rc_w rc_h],'Facecolor','#dcded5');
           
        
             %{
              if ~isempty(curr_Eye_x_error{i,j})
                plot((curr_Eye_x_error{i,j}(:,sel)+rc_w/2)',(curr_Eye_y_error{i,j}(:,sel)+rc_h/2)','-r','color','#211902','LineWidth',3);
              end
             %}
             hold on
             if ~isempty(curr_Eye_x{i,j})
                 tmp_x = curr_Eye_x{i,j}(:,sel);
                 num_x = size(tmp_x);
                 num_x = num_x(1);

                plot((curr_Eye_x{i,j}(max(num_x - Plot_Last+1,1):end,sel)+rc_w/2)',(curr_Eye_y{i,j}(max(num_x - Plot_Last+1,1):end,sel)+rc_h/2)','-r','color',MarkerColorEvent{j},'LineWidth',2);
               
             end
          xlim([-20+rc_w/2,20+rc_w/2]);
          ylim([-20+rc_h/2,20+rc_h/2]);
        %  box on;
          
          axis equal;
          axis off;
          
      end
 
          

          %{
          sel = PSTH_Mean_Time{i,j}>=PlotInterval(1) & PSTH_Mean_Time{i,j}<=PlotInterval(2);  
          if PlotError
            shadedErrorBar( PSTH_Mean_Time_error{i,j}(sel),PSTH_Mean_Event_error{i,j}(sel), PSTH_SEM_Event_error{i,j}(sel),'lineprops',{'#211902'},'transparent',1,'patchSaturation',0.3);
            %Error trials
          end
            
                 hold on 
                 shadedErrorBar(PSTH_Mean_Time{i,j}(sel),PSTH_Mean_Event{i,j}(sel),PSTH_SEM_Event{i,j}(sel),'lineprops',{Color(i,:)},'transparent',1,'patchSaturation',0.3);
           
           %Correct trials
        %Add markers
        
        plot([0,0],[0,max_fr_seq(i)],'--k');
          %}
        %{
        %Add reward time
        if j>1
            plot([RewardTimeSeq(i,j),RewardTimeSeq(i,j)],[0,max_fr_seq(i)],'--k','color',[0.9290, 0.6940, 0.1250]);
           
        else
            %Add fix off time
            plot([FixOffTimeSeq(i,j),FixOffTimeSeq(i,j)],[0,max_fr_seq(i)],'--k','color',[0, 0.4470, 0.7410]);
        end
        %Add pre-target off time
        if j>2
            plot([LastTgtOffTimeSeq(i,j),LastTgtOffTimeSeq(i,j)],[0,max_fr_seq(i)],'--k','color','r');
        end
        if j == length(SegEventStr)
           plot([LastFixOffTimeSeq(i,j),LastFixOffTimeSeq(i,j)],[0,max_fr_seq(i)],'--k','color',[0, 0.4470, 0.7410]);   
        end
        
          
          set(findall(gca, 'Type', 'Line'),'LineWidth',3);
          
          ylim([0,max(max_fr_seq(i),1)]);
          
          xlabel(SegEventStr{j});
          ylabel('FR(Hz)');
          set(gca,'LineWidth',3,'FontSize',10,'FontWeight','Bold');
          if j>1
            h=gca;
            h.YAxis.Visible = 'off'; % remove y-axis
            hold on;
          end
        %}
      end
      
      
      
  end
    


    
    
    
%end %end of the if plotfigure




%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Output to
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%files%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Setup output directory
Workingdirectory=pwd;


%Set up output path according to the system
OperationSystem = computer;%Get the operation system information

if strcmp(OperationSystem(1:3),"PCW")  
    %PC
     MarkerFolder='DataHub';
elseif strcmp(OperationSystem(1:3),"MAC")  
    %Mac
     MarkerFolder='DataAnalysis';
 end
Flag=strfind(Workingdirectory,MarkerFolder);
BasicDirectory=Workingdirectory(1:Flag+length(MarkerFolder));

%%%%%%%%%%%Export data%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%OutputFolerName='ForagingTask';
OutPath=strcat(BasicDirectory,'Results',BasicDirectory(end));%,OutputFolerName);
if ~exist(OutPath)
    mkdir(OutPath)
end
cd(OutPath);
FileName=Data.FileName;
FileName=FileName(1:end-4);


%Get File Information
num = regexp(FileName, '\d+', 'match');
RecordDateOriginal=cell2mat(num(1));

NeuronNum=cell2mat(num(2));
ChannelNum=SpikeChannel;


%To get the monkey name
RemainFiles=erase(FileName,RecordDateOriginal);
Delimiter=find(isstrprop(RemainFiles,'upper')==1);
MonkeyName=RemainFiles(Delimiter(1):Delimiter(2)-1);


%Output file name
OutputFileName=sprintf('%s_%s_N%s_C%s.mat',MonkeyName,string(RecordDateOriginal),string(NeuronNum),string(ChannelNum));
OutputFigureName=sprintf('%s_%s_N%s_C%s',MonkeyName,string(RecordDateOriginal),string(NeuronNum),string(ChannelNum));

%Store figures in common folders for each neuron
%{
if TaskCode==FRACTAL_PATH_TASK
    ProtocolMarker='FractalPath';
elseif  TaskCode==SCENE_PATH_TASK
    ProtocolMarker='ScenePath';
end
%}
ProtocolMarker = 'SeqPathTask';
FolderType = 'Separate';
FigurePath=strcat(BasicDirectory,'Figures',BasicDirectory(end),ProtocolMarker,BasicDirectory(end),FolderType,BasicDirectory(end));%,OutputFolerName);
if ~exist(FigurePath)
    mkdir(FigurePath)
end

if ShowFigureFlag & OutputFlag
    for i=1:length(fig)
        cd(FigurePath);
        OutputFigureNameFig=[];
        titlestr=figtitlestr_save{i};
        
        FolderName =OutputFigureName;
       
         if ~exist(FolderName,'dir')
          mkdir(FolderName);
         end
          cd(FolderName);
    
        OutputFigureNameFig=strcat(titlestr,'.jpg');
        saveas(fig{i},OutputFigureNameFig);
        
        
    end
        
    
    disp('Figures have been exported to the common folder');
    
end




end %end of the spike channel loop

end
function code = SeqCode(loc);
TrialIndex = loc(:,1);
LocLib = [270,0,90,180,315,45,135,225];
CodeLib = [1:length(LocLib)];

code = loc;
for i = LocLib
    code(abs(code- i)<0.01) = CodeLib(LocLib == i);
    
end

code_end = zeros(size(code,1),1);
for i = 1:size(code,2)
    code_end = code_end + code(:,i)*10^(size(code,2)-i);
end

code = [TrialIndex,code,code_end];


end
