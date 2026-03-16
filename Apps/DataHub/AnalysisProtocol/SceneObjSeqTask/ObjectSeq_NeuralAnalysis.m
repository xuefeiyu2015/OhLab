function ObjectSeq_NeuralAnalysis(Data, StartTrial, EndTrial,ShowFigureFlag,OutputFlag);

ECode;


CheckCode=[TGTCD;TGTOFF];
AbnormalTrials=AbnormalTrialCheck(Data,CheckCode);
SegEvent=[SHOWFIXCD,TGTCD];%Event Flow of interest
    

 PlotInterval = [-500,1000];
 
 IntervalForFP =[100,500];
 IntervalForTgt=[-200,0];

 

%SegEventEnd=[SHOWFIXCD,FIXOFF,TGTOFF,RWDOFFCD ,NaN,EXTRAFP];%Event end
%SegEventEnd=[TGTOFF];%Event end
SegEventEnd=NaN*ones(1,length(SegEvent));
SegEventStr={'FixOn','1TGT','2TGT','3TGT','4TGT','LastFix'};
MarkerSizeEvent = 20*ones(1,length(SegEventStr));
MarkerSizeEvent([1,end]) = 10;
MarkerColorEvent = {'#44AA99','#FE6100','#FE6100','#FE6100','#FE6100','#44AA99'};
DistratorColorEvent = {'#44AA99','#648FFF','#648FFF','#648FFF','#648FFF','#44AA99'};

PlotStartOffset=500*ones(1,length(SegEvent));%in ms
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
ECode;
FileName=Data.FileName;
TaskType=Data.TaskType;
TaskCode=Data.TaskCode;

Selection=StartTrial:EndTrial;
%{
%Abnormal check %Especially For Robin081120FORAGEXPLO4.rst
CheckCode=SCENE_ON;%This code should appear once in each trial
AbnormalTrials=AbnormalTrialCheck(Data,CheckCode);

Selection=Selection(~ismember(Selection,AbnormalTrials));
%}
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

SceneIndex = Data.SceneIndex(Selection,2);
SceneName = Data.SceneIndex(Selection,3);

%Sequence information from the library
ParaLib=Data.ParaLib;
DistractorType = ParaLib('DistractorType');
DistractorType = DistractorType(Selection,2);
UniqueDistractor = unique(DistractorType)

UniqueScene = unique(SceneName);

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

SpikeChannel=Data.SpikeChannel;
SpikeTime=squeeze(Data.SpikeTimeData(:,:,Selection))';
%SpikeBin=squeeze(Data.SpikeBinData(:,:,Selection))';


%Segregate spike data from fractal onset to fractal offset

SegmentData=EventSegTool(SegEvent,SegEventEnd,PlotStartOffset,PlotStopOffset,EventChannel,EventTimeChannel,SpikeTime,StepSize,BinWidth);

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
%Re-organize into a line form

    MeanPSTH=SegmentData.MeanPSTHEach;%=[PSTH_Event_Time',PSTH_Event'];
    MeanPSTH_Time=MeanPSTH(:,1);    
    MeanPSTH_FR=MeanPSTH(:,2);
    
    MeanPSTH_Time ={ MeanPSTH_Time{1}{1},MeanPSTH_Time{2}{1:end},MeanPSTH_Time{1}{end}};
    MeanPSTH_FR = { MeanPSTH_FR{1}{1},MeanPSTH_FR{2}{1:end},MeanPSTH_FR{1}{end}};
    
    EndCheck =ReproduceFromEvent(EventChannel,SHOWFIXCD);
    EndCheck = double((cellfun(@numel,EndCheck)==2));
    EndCheck(EndCheck==0)=NaN;
    
    
    

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


SuccessfulMat=ReorganizeEye( ReproduceFromEvent(EventChannel,[MATCHOBJ,NONMATCHOBJ]));
SuccessfulMat(SuccessfulMat==MATCHOBJ)=1;
SuccessfulMat(SuccessfulMat==NONMATCHOBJ)=0;
SuccessfulMatFull = [ones(size(SuccessfulMat,1),1),SuccessfulMat,EndCheck'];


%Align Eye Data with events
EyeForAlign.Bin =EyeBinWidth ;
EyeForAlign.X =EyeChannel_X_L;
EyeForAlign.Y =EyeChannel_Y_L;
EyeForAlign.MaxTime = cellfun(@numel,EyeChannel_X_L);

Eye_Aligned = AlignEyeDataWithEvents(SegEvent,SegEventEnd,PlotStartOffset,PlotStopOffset,EventChannel,EventTimeChannel,EyeForAlign);

Eye_Aligned_X=Eye_Aligned.X;
Eye_Aligned_Y=Eye_Aligned.Y;
Eye_Aligned_Time = Eye_Aligned.Time;


   
Eye_Aligned_Time ={ Eye_Aligned_Time{1,1},Eye_Aligned_Time{2,1:end},Eye_Aligned_Time{1,end}};
Eye_Aligned_X = { Eye_Aligned_X{1,1},Eye_Aligned_X{2,1:end},Eye_Aligned_X{1,end}};
Eye_Aligned_Y = { Eye_Aligned_Y{1,1},Eye_Aligned_Y{2,1:end},Eye_Aligned_Y{1,end}};
    
     

RewardTimeSeq=[];
LastTgtOffTimeSeq=[];

for m = 1 : length(UniqueDistractor)
%Loop over the Scenes 
for i = 1:length(UniqueScene)
    sel = SceneName == UniqueScene(i) & DistractorType ==UniqueDistractor(m);
    SucIndex = SuccessfulMatFull(sel,:); 
    max_fr_seq(i) = 0; 
    
    TotalTrialsSeq(i) = sum(sel);
    CorrectTrialSeq(i) = sum(CorrectTrialIndicator(sel) ==CORRECTCD);
    for j = 1:length(SegEventStr)
        %Neural responses
        
            SuccessfulMatcurr=SuccessfulMatFull(:,j);
            
            PSTH_Time_curr = nanmean(MeanPSTH_Time{j}(sel(SuccessfulMatcurr == 1),:),1);
           
            PSTH_curr = MeanPSTH_FR{j}(sel(SuccessfulMatcurr == 1) ,:);
            
             PSTH_curr_error = MeanPSTH_FR{j}(sel(SuccessfulMatcurr == 0) ,:);
             PSTH_Time_curr_error =  nanmean(MeanPSTH_Time{j}(sel(SuccessfulMatcurr == 0),:),1);
        
        
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
        %{
        %Eye trace
            
        curr_Eye_x{i,j} = Eye_Aligned_X{j}(sel&(SuccessfulMatcurr == 1),:);
        curr_Eye_y{i,j} = Eye_Aligned_Y{j}(sel&(SuccessfulMatcurr == 1),:);
        
        curr_Eye_x_error{i,j} = Eye_Aligned_X{j}(sel&(SuccessfulMatcurr == 0),:);
        curr_Eye_y_error{i,j} = Eye_Aligned_Y{j}(sel&(SuccessfulMatcurr == 0),:);
       
        curr_Eye_time{i,j}=nanmean(Eye_Aligned_Time{j}(sel&(SuccessfulMatcurr == 1),:),1);
        
         %}
        
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
            
      
        
        
        if j >1
         %Reward Time shown
         
         
         
         RewardTime_tmp=RewardTime(:,j-1)-MarkerTime;

         RewardTimeSeq(i,j)=nanmean(RewardTime_tmp(sel&(SuccessfulMatcurr == 1)));
         
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
end %end of the loop over the unique distractor

if ShowFigureFlag
    
  %Loop over the sequence type
  StartNO = 130;
      Fig1Size = [50,100, 1200,400];
      ylevel2 = 0.3;
      ysize = 0.4;

  %PlotEyeInterval=[-200,0];
  
      
  for i =1 :length(UniqueScene)
      
      figtitlestr{i}=sprintf("Scene ObjectSeq Task, HyperSet:%d,TotalTrial:%d,Correct:%d",UniqueScene(i),TotalTrialsSeq(i),CorrectTrialSeq(i));
      
      Fig1Size(2)=Fig1Size(2)+(i-1)*100;
      fig{i}=PrepareFigure(StartNO+i-1,'w',Fig1Size,'Name',figtitlestr{i});
      
      
      for j = 2:length(SegEventStr)
          %First plot the target location
          %title(SegEventStr{j});
          %{
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
          
          time_tmp=-PlotStartOffset:-PlotStartOffset+size(curr_Eye_x{i,j},2);
          sel = time_tmp>=PlotEyeInterval(1) & time_tmp<=PlotEyeInterval(2); 
           rectangle('Position',[0 0 rc_w rc_h],'Facecolor','#dcded5');
           
        
             %{
              if ~isempty(curr_Eye_x_error{i,j})
                plot((curr_Eye_x_error{i,j}(:,sel)+rc_w/2)',(curr_Eye_y_error{i,j}(:,sel)+rc_h/2)','-r','color','#211902','LineWidth',3);
              end
             %}
             hold on
             if ~isempty(curr_Eye_x{i,j})
                plot((curr_Eye_x{i,j}(:,sel)+rc_w/2)',(curr_Eye_y{i,j}(:,sel)+rc_h/2)','-r','color',MarkerColorEvent{j},'LineWidth',2);
               
             end
          xlim([-20+rc_w/2,20+rc_w/2]);
          ylim([-20+rc_h/2,20+rc_h/2]);
        %  box on;
          
          axis equal;
          axis off;
          
          end
          %}
          
          
          %Neural Responses
          axes('Position',[0.04+(j-2)*0.13,ylevel2,0.12,ysize]);
          sel = PSTH_Mean_Time{i,j}>=PlotInterval(1) & PSTH_Mean_Time{i,j}<=PlotInterval(2);             
          shadedErrorBar( PSTH_Mean_Time_error{i,j}(sel),PSTH_Mean_Event_error{i,j}(sel), PSTH_SEM_Event_error{i,j}(sel),'lineprops',{'#211902'},'transparent',1,'patchSaturation',0.3);
            %Error trials
          hold on 
           shadedErrorBar(PSTH_Mean_Time{i,j}(sel),PSTH_Mean_Event{i,j}(sel),PSTH_SEM_Event{i,j}(sel),'lineprops',{Color(i,:)},'transparent',1,'patchSaturation',0.3);
        %Correct trials
        %Add markers
        
        plot([0,0],[0,max_fr_seq(i)],'--k');
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
          ylim([0,max_fr_seq(i)]);
          xlabel(SegEventStr{j});
          ylabel('FR(Hz)');
          set(gca,'LineWidth',3,'FontSize',10,'FontWeight','Bold');
          if j>1
            h=gca;
            h.YAxis.Visible = 'off'; % remove y-axis
            hold on;
          end
          
      end
      
      
      
  end
    


    
    
    
end


end






