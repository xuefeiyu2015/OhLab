function Learning_NeuralResponse(Data, StartTrial, EndTrial,ShowFigureFlag,OutputFlag);


ECode;
DataPath= Data.Path;
FileName=Data.FileName;
TaskType=Data.TaskType;
TaskCode=Data.TaskCode;



%Neural responses change during learning

codecheck = FindOutTime(Data.EventChannel,Data.EventTimeChannel,FIRSTTGT);
codecheck=codecheck(~isnan(codecheck ));

if isempty(codecheck)
    disp("Old events, auto transform segment events");
    
    CheckCode=[TGTCD,SCENE_ON;TGTOFF,SCENE_OFF];
    AbnormalTrials=AbnormalTrialCheck(Data,CheckCode);
    
    CheckCodeExtra=[SCENE_ON;SCENE_OFF];
    AbnormalTrialsExtra=AbnormalTrialCheck(Data,CheckCodeExtra);
    AbnormalTrials=[AbnormalTrials,AbnormalTrialsExtra];
    
    
    OldCode = 1;
    SegEvent=[SHOWFIXCD,TGTCD];%Event Flow of interest
    
   
   

else
    
%{\
%Abnormal check for equality
CheckCode=[[FIRSTTGT:FIFTHTGT,SCENE_ON];[TGTOFF,TGTOFF,TGTOFF,TGTOFF,TGTOFF,SCENE_OFF]];
AbnormalTrials=AbnormalTrialCheck(Data,CheckCode);
CheckCodeExtra=[SCENE_ON;SCENE_OFF];
AbnormalTrialsExtra=AbnormalTrialCheck(Data,CheckCodeExtra);
AbnormalTrials=[AbnormalTrials,AbnormalTrialsExtra];

%}
 OldCode =0;
 
 SegEvent=[SHOWFIXCD,FIRSTTGT:FIFTHTGT,FINALFIX];%Event Flow of interest
 
 
 
 
end



 PlotInterval = repmat([-500,1000],length(SegEvent),1);
 PlotInterval(end,:)=[-500,1500];%Extra length for the final reward
 PlotIntervalExtra=[-500,1000;-500,1000];%For the scene on/off
 %PlotEyeTrace = 1;
 %IntervalForFP =[100,500];
 %IntervalForTgt=[-200,0];
 %
 
 NTrialLearning = 10;%First N correct set trials, last N correct set trials

 PlotDiagram = 1;%Plot the original target location
 PlotAll= 1;%Plot all the sequence responses on one figure
 PlotError = 1;%Plot the error trials



%SegEventEnd=[SHOWFIXCD,FIXOFF,TGTOFF,RWDOFFCD ,NaN,EXTRAFP];%Event end
%SegEventEnd=[TGTOFF];%Event end
SegEventEnd=NaN*ones(1,length(SegEvent));
%SegEventEnd(1) =SCENE_OFF;

SegEventStr={'FixOn','1TGT','2TGT','3TGT','4TGT','5TGT','LastFix'};
% 
MarkerSizeEvent = 20*ones(1,length(SegEventStr));
MarkerSizeEvent([1,end]) = 10;

MarkerColorEvent = {'#44AA99','#FE6100','#FE6100','#FE6100','#FE6100','#FE6100','#44AA99'};
DistratorColorEvent = {'#44AA99','#648FFF','#648FFF','#648FFF','#648FFF','#648FFF','#44AA99'};
%}
PlotStartOffset=500*ones(1,length(SegEvent));%in ms
PlotStopOffset=NaN*ones(1,length(SegEvent));%in ms

%Add the scene responses for the scene based learning task
SegEventExtra = [SCENE_ON,SCENE_OFF];
SegEventEndExtra=NaN*ones(1,length(SegEventExtra));
SegEventStrExtra = {'SceneOn','SceneOff'};

PlotStartOffsetExtra=500*ones(1,length(SegEventExtra));%in ms
PlotStopOffsetExtra=NaN*ones(1,length(SegEventExtra));%in ms


%PlotStopOffset(1) = 0;


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

Color_Early = Color*0.3;
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

SequenceLength = sum(~isnan(TargetNum),2);
UnqiueSeqLength = unique(SequenceLength);

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
ChannelIDWhole=Data.ChannelID;
ClusterIDWhole=Data.ClusterID;

for index = 1:length(SpikeChannelWhole)
    SpikeChannel = SpikeChannelWhole(index);
    SpikeTime=squeeze(Data.SpikeTimeData(index,:,Selection))';
    
    ChannelID = ChannelIDWhole(index);
    ClusterID = ClusterIDWhole(index);
    

%SpikeBin=squeeze(Data.SpikeBinData(:,:,Selection))';


%
%It took a long time to chunk data according to different event marker.
%To save time, save the results to the data folder
%OutputFileName=sprintf('%s_%s_N%s_C%s.mat',MonkeyName,string(RecordDateOriginal),string(NeuronNum),string(ChannelNum));

path=DataPath;

cd(path);

ReChunk=0;%Temporally put here;
ReChunkExtra = 0;
%Get File Information
num = regexp(FileName, '\d+', 'match');
%RecordDateOriginal=cell2mat(num(1));

NeuronNum=cell2mat(num(2));
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
disp('preanalysized saccades has been loaded to the workspace') ;
   

end

%Retrive data from the field
%Total average for each event
MeanFR=SegmentData.MeanFR; %[Mean,sem]%Averall mean and sem
MeanFRRaw=SegmentData.MeanFRRaw;%=AverageFiringRate_Event;

%Relative time bin for each event
RelativeTimeMarker=SegmentData.RelativeTimeMarker;%=EventTimeMarker;
%Absolute time marker for each event
AbosoluteTimeMarker=SegmentData.AbsoluteTimeMarker;




%{
TgtOnTime_test=FindOutTime(EventChannel,EventTimeChannel,FRISTTGT);

FPOnTime_test=FindOutTime(EventChannel,EventTimeChannel,SHOWFIXCD);
tmp=AbosoluteTimeMarker{2}(:,1);

%}
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

%Add the scene on/off aligned responses to the task
SpikeFileNameExtra=strcat(FileName(1:end-4),'_Seq_Spike_N',string(NeuronNum),'_C',string(ChannelNum),'SCENEONOFFs','.mat');
if exist(SpikeFileNameExtra) & ~ReChunkExtra
    load(SpikeFileNameExtra);
    disp('preanalysized chuncked data aligned on scene-on/off has been loaded to the workspace') ;
else
   disp('Chuncking the extra data, please wait.') 

%Segregate spike data from fractal onset to fractal offset

SegmentDataExtra=EventSegTool(SegEventExtra,SegEventEndExtra,PlotStartOffsetExtra,PlotStopOffsetExtra,EventChannel,EventTimeChannel,SpikeTime,StepSize,BinWidth);

save(SpikeFileNameExtra,'SegmentDataExtra');
disp('Seqed spikes on scene on/off sucessfully exported'); 
load(SpikeFileNameExtra);
disp('preanalysized chuncked data aligned on scene on/off has been loaded to the workspace') ;
   

end
%Insert the aligned scene on/off responses to the whole responses 
MeanPSTHExtra=SegmentDataExtra.MeanPSTH;%=[PSTH_Event_Time',PSTH_Event'];
MeanPSTH_TimeExtra=MeanPSTHExtra(:,1);    
MeanPSTH_FRExtra=MeanPSTHExtra(:,2);

MeanPSTH_Time = {MeanPSTH_TimeExtra{1},MeanPSTH_Time{1:length(SegEvent)},MeanPSTH_TimeExtra{2}}';

MeanPSTH_FR = {MeanPSTH_FRExtra{1},MeanPSTH_FR{1:length(SegEvent)},MeanPSTH_FRExtra{2}}';


SegEventStrAll ={SegEventStrExtra{1},SegEventStr{1:length(SegEvent)},SegEventStrExtra{2}};
PlotIntervalAll = [PlotIntervalExtra(1,:);PlotInterval;PlotIntervalExtra(2,:)];

SegEventAll = [SegEventExtra(1),SegEvent,SegEventExtra(2)];

%{
RewardTime=FindOutTime(EventChannel,EventTimeChannel,RWDCD);
TgtOffTime=FindOutTime(EventChannel,EventTimeChannel,TGTOFF);
FixOffTime=FindOutTime(EventChannel,EventTimeChannel,FIXOFF);
%}
 

%TimeMarker=SegmentData.RelativeTimeMarker{1}; %[Mean,sem]%Averall mean and sem



%Trial use: here only use successful trials
%TrialsForUse =

%Correct Rate
TotalTrialNum = size(EventChannel,1);
CorrectTrialIndicator = cell2mat(ReproduceFromEvent(EventChannel,[ERROR,CORRECTCD]));
CorrectTrialNum = sum(CorrectTrialIndicator ==CORRECTCD);
ErrorTrialNum = sum(CorrectTrialIndicator ==ERROR);


SuccessfulMat=ReorganizeEye( ReproduceFromEvent(EventChannel,[MATCHOBJ,NONMATCHOBJ]));
SuccessfulMat(SuccessfulMat==MATCHOBJ)=1;
SuccessfulMat(SuccessfulMat==NONMATCHOBJ)=0;
SuccessfulMatFull = [ones(size(SuccessfulMat,1),1),SuccessfulMat,EndCheck'];

SuccessfulMatExtra =[ones(size(SuccessfulMat,1),1),(CorrectTrialIndicator ==CORRECTCD)'];%Only include successful scene off

SuccessfulMatFull =[SuccessfulMatExtra(:,1),SuccessfulMatFull,SuccessfulMatExtra(:,2)];%Only include successful scene off




RewardTimeSeq=[];
LastTgtOffTimeSeq=[];

%Loop over the sequence 
for i = 1:length(UniqueSeq)
%for i = 1:1
    sel = SeqCodeUse == UniqueSeq(i);
    %sel = SeqCodeUse == 1006;
    %
   % SucIndex = SuccessfulMatFull(sel,:);
    
    TargetLoc_curr_X = nanmean(TargetLoc(sel,:,1),1);
    TargetLoc_curr_Y = nanmean(TargetLoc(sel,:,2),1);
    curr_X(i,:) = [0,0,TargetLoc_curr_X,0,0];
    curr_Y(i,:) = [0,0,TargetLoc_curr_Y,0,0];
    
   
   % TotalTrialsSeq(i) = sum(sel);
   % CorrectTrialSeq(i) = sum(CorrectTrialIndicator(sel) ==CORRECTCD);
    %}
   max_fr_seq(i) = 0; 
    
    
    for j = 1:length(SegEventStrAll)
        %Neural responses
        
            %SuccessfulMatcurr=SuccessfulMatFull(:,j);
            
            %%%%%%%%%%%% If select all correct sets in the hyperset %%%%%%%
            %
            SuccessfulMatcurr=SuccessfulMatFull(:,j);
            select = ~isnan(SuccessfulMatFull(:,j));
            correctSelector = SuccessfulMatcurr==1 & sel;%correct trial index for the current scene
            
            select_idx = correctSelector(select == 1);
            error_curr =SuccessfulMatcurr==0 & sel;  
            select_error = error_curr(select==1);
            
           % select_eye_correct = correctSelector ;
     
            PSTH_Time_curr = nanmean(MeanPSTH_Time{j}(select_idx,:),1);         
            PSTH_curr = MeanPSTH_FR{j}(select_idx ,:);
 
             
            PSTH_curr_error = MeanPSTH_FR{j}(select_error ,:);
            PSTH_Time_curr_error =  nanmean(MeanPSTH_Time{j}(select_error ,:),1);
            
            %Separate the correct trials into the first N and last N
            %NTrialLearning
            if size(PSTH_curr,1)<NTrialLearning
                
                PSTH_curr_Early = NaN*PSTH_curr;
                PSTH_curr_Late = NaN*PSTH_curr;
            else
                PSTH_curr_Early = PSTH_curr(1:NTrialLearning,:);
                PSTH_curr_Late = PSTH_curr(end-NTrialLearning+1:end,:);
                
            end
     
            PSTH_Mean_Time{i,j} = PSTH_Time_curr;
        
        PSTH_Mean_Event_Early{i,j} = nanmean(PSTH_curr_Early,1);
        PSTH_SEM_Event_Early{i,j} = nanstd(PSTH_curr_Early,[],1)./sqrt(sum(~isnan(PSTH_curr_Early),1));
        
        PSTH_Mean_Event_Late{i,j} = nanmean(PSTH_curr_Late,1);
        PSTH_SEM_Event_Late{i,j} = nanstd(PSTH_curr_Late,[],1)./sqrt(sum(~isnan(PSTH_curr_Late),1));
        
        PSTH_Mean_Event{i,j} = nanmean(PSTH_curr,1);
        PSTH_SEM_Event{i,j} = nanstd(PSTH_curr,[],1)./sqrt(sum(~isnan(PSTH_curr),1));
        
        
        
        if size(PSTH_curr_error,1)<3
            PSTH_curr_error = NaN*PSTH_curr_error;
        end
        PSTH_Mean_Time_error{i,j} = PSTH_Time_curr_error;
        PSTH_Mean_Event_error{i,j} = nanmean(PSTH_curr_error ,1);
        PSTH_SEM_Event_error{i,j} = nanstd(PSTH_curr_error ,[],1)./sqrt(sum(~isnan(PSTH_curr_error ),1));
        
        
        sel_tmp =PSTH_Mean_Time{i,j}>=PlotIntervalAll(j,1) & PSTH_Mean_Time{i,j}<=PlotIntervalAll(j,2);
        if max(max(PSTH_Mean_Event_Early{i,j}(sel_tmp)+PSTH_SEM_Event_Early{i,j}(sel_tmp)),max(PSTH_Mean_Event_Late{i,j}(sel_tmp)+PSTH_SEM_Event_Late{i,j}(sel_tmp))) > max_fr_seq(i) 
            max_fr_seq(i) = max(max(PSTH_Mean_Event_Early{i,j}(sel_tmp)+PSTH_SEM_Event_Early{i,j}(sel_tmp)),max(PSTH_Mean_Event_Late{i,j}(sel_tmp)+PSTH_SEM_Event_Late{i,j}(sel_tmp)));
            
        end
       
        
        
        
    end
    
    
    
end
 
 
    


if ShowFigureFlag
    
  %Loop over the sequence type
  StartNO = 140;
  
 
  %PlotEyeInterval=[-200,0];
  Fig1Size = [50,100, 1400,1000];
  
  Fig2Size = [30,100, 1200,1200];
  
  %Try to plot all sequence on one or two figures
  TotalSeq = length(UniqueSeq);
  NumOfSeqInFig = 6;
  NumOfFig = max(ceil(TotalSeq/NumOfSeqInFig),1);
  
  
  

  for ii = 1:NumOfFig
      if min(UniqueSeq)>=1000
            figtitlestr{ii}=sprintf("Scene Path Task, HyperSet:%d,%d,%d,%d,%d,%d",UniqueSeq((ii-1)*NumOfSeqInFig+1:min(ii*NumOfSeqInFig,TotalSeq)));
            figtitlestr_save{ii} = sprintf("LearningScene_S%d%d%d%d%d%d_S%d_C_%d_Cl%d",UniqueSeq((ii-1)*NumOfSeqInFig+1:min(ii*NumOfSeqInFig,TotalSeq)),SpikeChannel,ChannelID,ClusterID);
            TaskMarker = 1;%Scene Path
      else
          figtitlestr{ii}=sprintf("Fractal Path Task, HyperSet:%d,%d,%d,%d,%d,%d",UniqueSeq((ii-1)*NumOfSeqInFig+1:min(ii*NumOfSeqInFig,TotalSeq)));
           figtitlestr_save{ii} = sprintf("LearningFrac_S%d%d%d%d%d%d_S%d_C_%d_Cl%ds",UniqueSeq((ii-1)*NumOfSeqInFig+1:min(ii*NumOfSeqInFig,TotalSeq)),SpikeChannel,ChannelID,ClusterID);
           TaskMarker = 2;%Fractal Path
          
      end
      
      Fig1Size(2)=Fig1Size(2)+(ii-1)*100;
      fig{ii}=PrepareFigure(StartNO+ii-1,'w',Fig1Size,'Name',figtitlestr{ii});
      
      if PlotDiagram
          %Plot all the seq diagram on one figure
          fig{NumOfFig+ii}=PrepareFigure(StartNO+NumOfFig+ii-1,'w',Fig2Size,'Name',strcat(figtitlestr{ii},'diagram'));
          
      end
      
      if PlotAll
          %Plot all the sequence responses on one figure
          fig{NumOfFig*(1+PlotDiagram)+ii}=PrepareFigure(StartNO+NumOfFig*(1+PlotDiagram)+ii-1,'w',Fig1Size,'Name',strcat(figtitlestr{ii},'all'));
      end
      
      
      %Loop over each hyperset and each events
      pidx = 1;
      for i = (ii-1)*NumOfSeqInFig+1:min(ii*NumOfSeqInFig,length(UniqueSeq))
          for j = 1:length(SegEventStrAll)
              %First plot the diagram if nessary
              if PlotDiagram
                  set(groot,'CurrentFigure',fig{NumOfFig+ii});
                  
                %No need to plot the scene onset or offset event
                
                    
                if j<= length(SegEventStrAll)-1 && j>=2
                    jcorr = j-1;
                  
                  %{
                    if ii == 1
                        subplot(NumOfSeqInFig,length(SegEventStr)-2,pidx);
                    else
                        subplot(TotalSeq - NumOfSeqInFig,length(SegEventStr)-2,pidx);
                    end
                    %}
                    if ii == 1
                        deidx = i;
                    else
                        deidx = i-NumOfSeqInFig;
                    end
                    axes('Position',[(jcorr-1)*0.136,0.8-(deidx-1)*0.155,0.15,0.15]);
                
                 rc_w = 40;rc_h = 40;
          
                ctr_x =curr_X(i,2*(jcorr-1)+1:2*jcorr);
                ctr_y =curr_Y(i,2*(jcorr-1)+1:2*jcorr);
          
                rectangle('Position',[0 0 rc_w rc_h],'Facecolor','#dcded5')
                hold on
          
                plot(ctr_x(1)+rc_w/2,ctr_y(1)+rc_h/2,'o','MarkerSize',MarkerSizeEvent(jcorr),'MarkerFaceColor',MarkerColorEvent{jcorr},'color',MarkerColorEvent{jcorr});
                plot(ctr_x(2)+rc_w/2,ctr_y(2)+rc_h/2,'o','MarkerSize',MarkerSizeEvent(jcorr),'MarkerFaceColor',DistratorColorEvent{jcorr},'color',DistratorColorEvent{jcorr});
                axis equal;
                axis off;
                %Add a vector poting from the last location to the current
                %location
                    if jcorr>1
                        last_x = curr_X(i,2*(jcorr-2)+1);
                        last_y = curr_Y(i,2*(jcorr-2)+1);
                        curr_x = ctr_x(1);
                        curr_y = ctr_y(1);
                        diff_x = curr_x - last_x;
                        diff_y = curr_y - last_y;
              
              
                        quiver(last_x+rc_w/2,last_y+rc_h/2,diff_x,diff_y,'LineWidth',3,'MaxHeadSize',5,'Color','r');
              
                    end % end of plotting quiver
                    
                    if mod(pidx,length(SegEventStr)-2) == 1
                        axes('Position',[0.01,0.8-(deidx-1)*0.155,0.05,0.05]);
                    
                        %annotation('textbox', [0.15, 0.85-(deidx-1)*0.155, 0.2, 0.04], 'string', sprintf('Seq%d',UniqueSeq(i)),'FontSize',20,'FontWeight','Bold','LineStyle','none','FitBoxToText','on');
                        ht = text(0.1,0.5,sprintf('Seq%d',UniqueSeq(i)));
                        set(ht,'Rotation',90);
                        set(ht,'FontSize',18);
                        axis off;
                    end
                   %{  
                    if i == min(ii*NumOfSeqInFig,length(UniqueSeq)) && j == length(SegEventStr)
                        if  NumOfFig > 1   
                            data_msg=sprintf(' %s, SpikeChannel: %d,RecordChannel: %d, ClusterID: % d (%d) ',FileName,SpikeChannel,ChannelID,ClusterID,ii);
                        else
                            data_msg=sprintf(' %s, SpikeChannel: %d,RecordChannel: %d, ClusterID: % d ',FileName,SpikeChannel,ChannelID,ClusterID);
     
                        end
                            annotation('textbox', [0.15, 0.95, 0.7, 0.04], 'string', data_msg,'FontSize',20,'FontWeight','Bold','LineStyle','none','FitBoxToText','on');
                    end
                %}
                
                    
                    
                end %End of the if j in 2 to end-1
                
                
                
                
                
                  
                  
                  
              end %End of the if plot diagram
              
              %%Plot the neural responses for early trials and late trials
              
              sel = PSTH_Mean_Time{i,j}>=PlotIntervalAll(j,1) & PSTH_Mean_Time{i,j}<=PlotIntervalAll(j,2); 
              
              
              set(groot,'CurrentFigure',fig{ii});
              if ii == 1
                  subplot(NumOfSeqInFig,length(SegEventStrAll),pidx);
              else
                  subplot(TotalSeq - NumOfSeqInFig,length(SegEventStrAll),pidx);
              end
              shadedErrorBar(PSTH_Mean_Time{i,j}(sel),PSTH_Mean_Event_Early{i,j}(sel),PSTH_SEM_Event_Early{i,j}(sel),'lineprops',{Color_Early(i,:)},'transparent',1,'patchSaturation',0.3);
              hold on
              shadedErrorBar(PSTH_Mean_Time{i,j}(sel),PSTH_Mean_Event_Late{i,j}(sel),PSTH_SEM_Event_Late{i,j}(sel),'lineprops',{Color(i,:)},'transparent',1,'patchSaturation',0.1);
              set(findall(gca, 'Type', 'Line'),'LineWidth',3);
              
              
              %ylim([0,max(max_fr_seq(i),1)]);
              ylim([0,max(max(max_fr_seq),1)]);
              %plot([0,0],[0,max_fr_seq(i)],'--k');
              plot([0,0],[0,max(max_fr_seq)],'--k');
              
              
            if j>1
                h=gca;
                h.YAxis.Visible = 'off'; % remove y-axis
                hold on;
            end
             set(gca,'LineWidth',3,'FontSize',10,'FontWeight','Bold');
            if mod(pidx,length(SegEventStrAll)) == 1
                ylabel(sprintf('Seq%d',UniqueSeq(i)),'FontSize',15);
            end
            if ceil(pidx/length(SegEventStrAll)) == NumOfSeqInFig && ii==1
                xlabel(SegEventStrAll{j},'FontSize',15); 
            elseif (ceil(pidx/length(SegEventStr)) == (length(UniqueSeq) - NumOfSeqInFig)) && ii==2
                xlabel(SegEventStrAll{j},'FontSize',15); 
                    
            end
            
            %%Plot the average responses across the whole task periods
            if PlotAll
                set(groot,'CurrentFigure',fig{NumOfFig*(1+PlotDiagram)+ii});
              if ii == 1
                  subplot(NumOfSeqInFig,length(SegEventStrAll),pidx);
              else
                  subplot(TotalSeq - NumOfSeqInFig,length(SegEventStrAll),pidx);
              end
              if PlotError 
                shadedErrorBar(PSTH_Mean_Time{i,j}(sel),PSTH_Mean_Event_error{i,j}(sel),PSTH_SEM_Event_error{i,j}(sel),'lineprops',{'#211902'},'transparent',1,'patchSaturation',0.3);
                hold on
              end
              
              shadedErrorBar(PSTH_Mean_Time{i,j}(sel),PSTH_Mean_Event{i,j}(sel),PSTH_SEM_Event{i,j}(sel),'lineprops',{Color(i,:)},'transparent',1,'patchSaturation',0.1);
              
              set(findall(gca, 'Type', 'Line'),'LineWidth',3);
              
              
              %ylim([0,max(max_fr_seq(i),1)]);
              ylim([0,max(max(max_fr_seq),1)]);
              %plot([0,0],[0,max_fr_seq(i)],'--k');
              plot([0,0],[0,max(max_fr_seq)],'--k');
              
              
                if j>1
                    h=gca;
                    h.YAxis.Visible = 'off'; % remove y-axis
                    hold on;
                end
                set(gca,'LineWidth',3,'FontSize',10,'FontWeight','Bold');
                if mod(pidx,length(SegEventStrAll)) == 1
                    ylabel(sprintf('Seq%d',UniqueSeq(i)),'FontSize',15);
                end
                if ceil(pidx/length(SegEventStrAll)) == NumOfSeqInFig && ii==1
                    xlabel(SegEventStrAll{j},'FontSize',15); 
                elseif (ceil(pidx/length(SegEventStrAll)) == (length(UniqueSeq) - NumOfSeqInFig)) && ii==2
                    xlabel(SegEventStrAll{j},'FontSize',15); 
                    
                end
                
           
                
            end
            
           
            
             %ylabel('UniqueSeq(i)');
              pidx=pidx+1;
          end
          
      end
      
      
      
    %Add data information
    set(groot,'CurrentFigure',fig{ii});
        if  NumOfFig > 1   
            data_msg=sprintf(' %s, SpikeChannel: %d,RecordChannel: %d, ClusterID: % d (%d) Early v.s. Late',FileName,SpikeChannel,ChannelID,ClusterID,ii);
        else
            data_msg=sprintf(' %s, SpikeChannel: %d,RecordChannel: %d, ClusterID: % d Early v.s. Late',FileName,SpikeChannel,ChannelID,ClusterID);
     
        end
        annotation('textbox', [0.15, 0.95, 0.7, 0.04], 'string', data_msg,'FontSize',20,'FontWeight','Bold','LineStyle','none','FitBoxToText','on');
        
    if PlotDiagram
        set(groot,'CurrentFigure',fig{NumOfFig+ii});
        if  NumOfFig > 1   
            data_msg=sprintf(' %s, SpikeChannel: %d,RecordChannel: %d, ClusterID: % d (%d) Diagram',FileName,SpikeChannel,ChannelID,ClusterID,ii);
        else
            data_msg=sprintf(' %s, SpikeChannel: %d,RecordChannel: %d, ClusterID: % d Diagram',FileName,SpikeChannel,ChannelID,ClusterID);
     
        end
        annotation('textbox', [0.15, 0.95, 0.7, 0.04], 'string', data_msg,'FontSize',20,'FontWeight','Bold','LineStyle','none','FitBoxToText','on');
      
    end
    
    if PlotAll
        set(groot,'CurrentFigure',fig{NumOfFig*(1+PlotDiagram)+ii});
        if  NumOfFig > 1   
            data_msg=sprintf(' %s, SpikeChannel: %d,RecordChannel: %d, ClusterID: % d (%d) Correct v.s. Error',FileName,SpikeChannel,ChannelID,ClusterID,ii);
        else
            data_msg=sprintf(' %s, SpikeChannel: %d,RecordChannel: %d, ClusterID: % d Correct v.s. Error',FileName,SpikeChannel,ChannelID,ClusterID);
     
        end
        annotation('textbox', [0.15, 0.95, 0.7, 0.04], 'string', data_msg,'FontSize',20,'FontWeight','Bold','LineStyle','none','FitBoxToText','on');
      
    end
    
    
    
    
      
  end %end of the numOfFig loop
  
 
      
  

    
    
    
end %end of the if plotfigure




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
%OutputFileName=sprintf('%s_%s_N%s_C%s.mat',MonkeyName,string(RecordDateOriginal),string(NeuronNum),string(ChannelNum));
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
OutputFolerName = 'AllInOneFig';
if TaskMarker == 1
    FolderType = 'ScenePath';
else
    FolderType = 'FractalPath';
end
FigurePath=strcat(BasicDirectory,'Figures',BasicDirectory(end),ProtocolMarker,BasicDirectory(end),OutputFolerName,BasicDirectory(end),FolderType,BasicDirectory(end));
if ~exist(FigurePath)
    mkdir(FigurePath)
end

if ShowFigureFlag & OutputFlag
    for i=1:NumOfFig
        
        %OutputFigureNameFig=[];
        
        titlestr=figtitlestr_save{i};
        
       % FolderName =OutputFigureName;
       %Into a common folder
       %Into 3 folders: Diagram; Learning; Correct
       %Diagram
       if PlotDiagram
           cd(FigurePath);
           FolderName ='Diagram';
           if ~exist(FolderName,'dir')
            mkdir(FolderName);
           end
            cd(FolderName);
            % OutputFigureNameFig=strcat(titlestr,'.pdf');
            OutputFigureNameFig=strcat(titlestr,'.jpg');
            saveas(fig{i+NumOfFig},OutputFigureNameFig);
             %print(fig{i+NumOfFig}, '-dpdf', OutputFigureNameFig,'-bestfit'); 
            % print(fig{i+NumOfFig}, '-dpdf', OutputFigureNameFig); 
           
       end
       %Correct
       if PlotAll
           cd(FigurePath);
           FolderName ='Correct';
           if ~exist(FolderName,'dir')
            mkdir(FolderName);
           end
            cd(FolderName);
            %OutputFigureNameFig=strcat(titlestr,'.pdf');
            OutputFigureNameFig=strcat(titlestr,'.jpg');
             saveas(fig{i+NumOfFig*(PlotDiagram+1)},OutputFigureNameFig);
             %print(fig{i+NumOfFig*(PlotDiagram+1)}, '-dpdf', OutputFigureNameFig,'-bestfit'); 
       end
       %Learning
       cd(FigurePath);
       FolderName ='Learning';
       if ~exist(FolderName,'dir')
          mkdir(FolderName);
       end
       cd(FolderName);
    
        %OutputFigureNameFig=strcat(titlestr,'.pdf');
        OutputFigureNameFig=strcat(titlestr,'.jpg');
        saveas(fig{i},OutputFigureNameFig);
        %print(fig{i}, '-dpdf', OutputFigureNameFig,'-bestfit'); 
        %print(fig{i}, '-dpdf', OutputFigureNameFig); 
        
        
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



