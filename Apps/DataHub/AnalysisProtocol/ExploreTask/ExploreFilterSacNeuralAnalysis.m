function ExploreFilterSacNeuralAnalysis(Data, StartTrial, EndTrial,ShowFigureFlag,OutputFlag);
%Function to Analysis the neural responses to saccade and gaze during the
%exploration period
%Using the traditional thresholding method to filter saccades
ECode;
TaskType=Data.TaskType;
DataPath= Data.Path;
FileName=Data.FileName;

Selection=StartTrial:EndTrial;%Select the trials according to the pannel
%Abnormal check %Especially For Robin081120FORAGEXPLO4.rst
CheckCode=SCENE_ON;%This code should appear once in each trial
AbnormalTrials=AbnormalTrialCheck(Data,CheckCode);
Selection=Selection(~ismember(Selection,AbnormalTrials));
TrialNum=length(Selection);
%User defined parameter for the neural/behavioral analysis
TargetOnOffset=100;%in ms, the time selection bin after target onset; 
StartDelay=-50;%in ms, the time selection bin after fixation offset;

%For Plot of PSTH aligned on Each fixation on and each fixation Off
PlotFixationOn_StartOffset=100;
PlotFixationOn_StopOffset=300;%Duration of plotting

PlotFixationOff_StartOffset=300;
PlotFixationOff_StopOffset=100;%Duration of plotting

PlotSaccadeOn_StartOffset=100;
PlotSaccadeOn_StopOffset=300;

PlotSaccadeOff_StartOffset=300;
PlotSaccadeOff_StopOffset=100;

%Criterias for selecting saccades
V_Threshold=40;
ContinueBin=5;
ISI_Threshold=20;


%For segmentation of  the neural responses 
SegEvent=[SHOWFIXCD];%Event Flow of interest
SegEventEnd=[TGTCD];%Event end
SegEventStr={'Explore'};

PlotStartOffset=[-StartDelay];%in ms
PlotStopOffset=[TargetOnOffset];%in ms

EventTypeFlag=[2];%1:Number same with the trial number, 2:Number same with the target number

%Spon Event
SponEvent=[SCENE_OFF];
SponPreTimeBuffer=-500;
SponEventEnd=NaN;
SponStartOffset=0;
SponStopOffset=1000;

%Parameter to define the sliding window to measure firing rates
StepSize=10;%in ms
BinWidth=15;%in ms

%Plot Raster Flag
PlotRasterFlag=0;


%ResponseLatencyThreshold=3;

%Some user defined parameter for visualization
%Color Definition for each scenes
Color={'r','g','b','m',[0.4,0.4,0],[0.28,0.78,0.89],'k',[0.6,0.1,1],[0.6,0.1,0.4],[0.3,0.5,1],[0.6,0.5,0.1],[0,0.8,0.4]};
ColorForEachMode=[0.83 0.14 0.14
             1.00 0.54 0.00
             0.47 0.25 0.80
             0.25 0.80 0.54];
AssociateBlock=[10,11,0,1];
Blockstr={'FixLoc,FixSceneSeq','FixLoc,RandSceneSeq','RandLoc,FixSceneSeq','RandLoc,RandSceneSeq'};


%BlockSepTrial=unique(SceneLib(:,1));
BlockType=Data.BlockType(:,:);

BlockType=BlockType(Selection,2);

%Get File Information
num = regexp(Data.FileName, '\d+', 'match');
RecordDateOriginal=cell2mat(num(1));
RecordDate=datetime(RecordDateOriginal,'InputFormat','MMddyy');

%Correct early block type
%Before 03,18,2020, I only use 0 to indicate fix scene sequence and 1 to
%indicate random scene sequence, both with fixed location sequence
if RecordDate<datestr('03,17,2020')
    
   BlockType(BlockType==0)=10;
    BlockType(BlockType==1)=11;
   
    
end

%{
BlockType(BlockType==10)=11;
BlockType(BlockType==0)=1;
%}


%BlockType(BlockType(:,3)>0,3)=999;



%GetEyeData
EyeData=Data.EyeDataRaw;
EyeChannel_X_L=EyeData(Selection,1);
EyeChannel_Y_L=EyeData(Selection,2);
%EyeChannel_X_R=EyeData(Selection,3);
%EyeChannel_Y_R=EyeData(Selection,4);

EyeBinWidth=Data.EyeBinWidth;
EyeCalInfo=Data.EyeCalInfo;%Each row: X_Gain,X_Bias; Y_Gain; Y_bias; X=(raw+X_bias)*X_Gain
EyeCalInfo_L=EyeCalInfo(1,:);

%Load event data
EventChannel=Data.EventChannel(Selection,:);
EventTimeChannel=Data.EventTimeChannel(Selection,:);
EventBin=Data.EventBin(Selection,:);
%Get the information about the extra FP
%
%Recalibrate eyes according to the extra FP
ParaLib=Data.ParaLib;

if(ismember('ExtraFixPointID,Style,Size,Loc_X,Loc_Y',keys(ParaLib)))
ExtraFP_Info=ParaLib('ExtraFixPointID,Style,Size,Loc_X,Loc_Y');
ExtraFP_Loc=ExtraFP_Info(Selection,end-1:end);
Ref_Loc=ExtraFP_Loc;

%Get the information about the extra FP
%ExtraFPOnset=FindOutTime(EventChannel,EventTimeChannel,EXTRAFP);
ExtraEyeIn=FindOutTime(EventChannel,EventBin,EXTRAFIXSTARTCD);
ExtraEyeOff=FindOutTime(EventChannel,EventBin,EXTRAFIXOFF);

EyeInBin_X=SelectEyeInterval(EyeChannel_X_L,ExtraEyeIn,ExtraEyeOff);

EyeInBin_Y=SelectEyeInterval(EyeChannel_Y_L,ExtraEyeIn,ExtraEyeOff);


else
    
FP_Loc=Data.FixLoc(Selection,:);
Ref_Loc=FP_Loc;
%Get the information about the extra FP
%ExtraFPOnset=FindOutTime(EventChannel,EventTimeChannel,EXTRAFP);
EyeIn=FindOutTime(EventChannel,EventBin,FIXSTARTCD );
EyeOff=FindOutTime(EventChannel,EventBin,FIXOFF);

EyeInBin_X=SelectEyeInterval(EyeChannel_X_L,EyeIn,EyeOff);

EyeInBin_Y=SelectEyeInterval(EyeChannel_Y_L,EyeIn,EyeOff);
end
EyeInBin_X=EyeInBin_X{1};
EyeInBin_Y=EyeInBin_Y{1};
EyeDataForCal={EyeInBin_X,EyeInBin_Y};

EyeDataRaw={EyeChannel_X_L,EyeChannel_Y_L};
CalPlotCheck=0;%Examine the calibration results
CaledEye=EyeCalibration(EyeDataRaw,EyeDataForCal,EyeCalInfo_L,Ref_Loc,CalPlotCheck);


EyeChannel_X_L=CaledEye.EyeTrace{1};
EyeChannel_Y_L=CaledEye.EyeTrace{2};

%}

%Load the spike channel
SpikeChannel=Data.SpikeChannel;
SpikeTime=squeeze(Data.SpikeTimeData(:,:,Selection))';


%Get scene and fractals information
SceneName=Data.SceneIndex(Selection,3);
TargetName=Data.TargetIndex(Selection,2:end,2);
TargetIndex=Data.TargetIndex(Selection,2:end,1);
TargetLoc=Data.TargetLoc(Selection,2:end,:);

%TargetNum=Data.TargetNum(Selection,2:end);
TargetNum=sum(~isnan(TargetName),2);
TargetSize=Data.TargetSize(Selection,2:end);
TargetPolarAngle=Data.TargetPolarLoc(Selection,2:end,1);
TargetPolarEccentricity=Data.TargetPolarLoc(Selection,2:end,2);

ChoiceOutcome=Data.ChoiceOutcome(Selection,2:end);

%Select the trials at lease has one case of advance eye in

AdvanceInTargetNum=sum(EventChannel==ADVEYE_IN,2);
%AdvanceInTargetTrial=AdvanceInTargetNum>0;
AdvanceInTargetTrial=AdvanceInTargetNum==TargetNum;

strAdvanceIn=fprintf('There are %d Advance-Eye-In trials\n',sum(AdvanceInTargetTrial));


%Select Searching time
 FixOffTime=FindOutTime(EventChannel,EventTimeChannel,FIXOFF);
 %TargetOnTime=FixOffTime+300;
 AdvanceEyeInTime=FindOutTime(EventChannel,EventTimeChannel,ADVEYE_IN );
 
 TargetOnTime=FindOutTime(EventChannel,EventTimeChannel,TGTCD);
 OutScreenTime=FindOutTime(EventChannel,EventTimeChannel,EYEOUTSCREEN);
 
 FPOnTime=FindOutTime(EventChannel,EventTimeChannel,SHOWFIXCD);
 %screen out the searching trials
 FixOffTime=FixOffTime(AdvanceInTargetTrial,:);
 OutScreenTime=OutScreenTime(AdvanceInTargetTrial,:);
 FPOnTime=FPOnTime(AdvanceInTargetTrial,:);
 
 FixOffTime=ScreenOutOverScreen(FixOffTime,OutScreenTime);
 
 
 TargetOnTime=TargetOnTime(AdvanceInTargetTrial,:);
 %Remove the col
 NaNCol_Check=sum(isnan(TargetOnTime),1);

 TargetOnTime(:,NaNCol_Check==size(TargetOnTime,1))=[];
 
 SceneName=SceneName(AdvanceInTargetTrial,:);
 
 EyeChannel_X_L=EyeChannel_X_L(AdvanceInTargetTrial);
 EyeChannel_Y_L=EyeChannel_Y_L(AdvanceInTargetTrial);
 
 TargetLoc=TargetLoc(AdvanceInTargetTrial,:,:);
 TargetPolarAngle=TargetPolarAngle(AdvanceInTargetTrial,:);
 
 TargetNum=TargetNum(AdvanceInTargetTrial);
 BlockType=BlockType(AdvanceInTargetTrial);
 
 [UniqueBlockType,orderindex,reprojectorderindex]=unique(BlockType,'row');
NumberOfBlocks=max(reprojectorderindex);
 
 %If the target onset time and fixofftime dimension are not aligned 
 if size(TargetOnTime,2)~=size(FixOffTime,2)
     FixOffTime=FixOffTime(:,1:size(TargetOnTime,2));
     disp('Warning!FixOffTime Corrected!!')
 end
     
 
SearchingTime= (TargetOnTime-FixOffTime)/1000;%Change into s for presentation


%Segregate spike data from Fixation onset to target onset
SegmentData=EventSegTool(SegEvent,SegEventEnd,PlotStartOffset,PlotStopOffset,EventChannel,EventTimeChannel,SpikeTime,StepSize,BinWidth);

%Retrive data from the field
%Total average for each event
MeanFR=SegmentData.MeanFR; %[Mean,sem]%Averall mean and sem
MeanFRRaw=SegmentData.MeanFRRaw;%=AverageFiringRate_Event;

%PSTH(Spike Density Function) for each event
MeanPSTH=SegmentData.MeanPSTH;%=[PSTH_Event_Time',PSTH_Event'];
MeanPSTH_Time=MeanPSTH{1};
MeanPSTH_FR=MeanPSTH{2};

%Raster and raw spike count PSTH
AlignedRaster=SegmentData.RasterAfterAlign{1};
PSTH_Event_SpikeCount=SegmentData.PSTH_Event_SpikeCount;


%Relative time bin for each event
RelativeTimeMarker=SegmentData.RelativeTimeMarker{1};%=EventTimeMarker;
%Absolute time marker for each event
AbosoluteTimeMarker=SegmentData.AbsoluteTimeMarker{1};


%Reorganize parameters into line library
%Rearrange target into lines
%Target
TargetIndex_Line=TargetIndex;
TargetIndex_Line=reshape(TargetIndex_Line,numel(TargetIndex_Line),1);

NoNNaNSelector=~isnan(TargetIndex_Line);
TargetIndex_Line=TargetIndex_Line(~isnan(TargetIndex_Line));

TargetName_Line=TargetName;
TargetName_Line=reshape(TargetName_Line,numel(TargetName_Line),1);
TargetName_Line=TargetName_Line(~isnan(TargetName_Line));

%TargetNumLine
TargetNum_Line=repmat(TargetNum,1,max(TargetNum));
TargetNum_Line=reshape(TargetNum_Line,numel(TargetNum_Line),1);
TargetNum_Line=TargetNum_Line(NoNNaNSelector);

%Rearrange block type into lines
BlockType_Line=repmat(BlockType,1,max(TargetNum));
BlockType_Line=reshape(BlockType_Line,numel(BlockType_Line),1);
BlockType_Line=BlockType_Line(NoNNaNSelector);

%Rearrange searching time into lines
SearchingTime_Line=reshape(SearchingTime,numel(SearchingTime),1);
SearchingTime_Line=SearchingTime_Line(NoNNaNSelector);

%Rearrange_trialIndex into lines
TrialIndex=[1:size(SearchingTime,1)]';
TrialIndex_Line=repmat(TrialIndex,1,max(TargetNum));
TrialIndex_Line=reshape(TrialIndex_Line,numel(TrialIndex_Line),1);
TrialIndex_Line=TrialIndex_Line(NoNNaNSelector);

%FR for spon
EventTimeChannelForSpon=EventTimeChannel-SponPreTimeBuffer;

SegmentDataSpon=EventSegTool(SponEvent,SponEventEnd,SponStartOffset,SponStopOffset,EventChannel,EventTimeChannelForSpon,SpikeTime,StepSize,BinWidth);

Mean_Spon=SegmentDataSpon.MeanFR;
Mean_Spon_ForEachTrial=SegmentDataSpon.MeanFRRaw{1};

Mean_Spon_ForEachTrial_Line=repmat(Mean_Spon_ForEachTrial,1,max(TargetNum));
Mean_Spon_ForEachTrial_Line=reshape(Mean_Spon_ForEachTrial_Line,numel(Mean_Spon_ForEachTrial_Line),1);
Mean_Spon_ForEachTrial_Line=Mean_Spon_ForEachTrial_Line(NoNNaNSelector);


%Select out eye trace for each target
SegEventTime=FindOutTime(EventChannel,EventTimeChannel,SegEvent(1));%For this protocol, there is only one event
SegEventEndTime=FindOutTime(EventChannel,EventTimeChannel,SegEventEnd(1));%For this protocol, there is only one event


StartTimeBin=floor((SegEventTime-PlotStartOffset(1))/EyeBinWidth);         
EndTimeBin=ceil((SegEventEndTime+PlotStopOffset(1))/EyeBinWidth);
         
Eye_Explore_X=SelectEyeInterval(EyeChannel_X_L,StartTimeBin,EndTimeBin)';


Eye_Explore_Y=SelectEyeInterval(EyeChannel_Y_L,StartTimeBin,EndTimeBin)';



%Eye_Explore_Line=ReorganizeData(Eye_Explore_X,Eye_Explore_Y);
%Eye_Explore_Line=Eye_Explore_Line(NoNNaNSelector,:);
EyeChannel_X_Explore=ReorganizeData(Eye_Explore_X);
EyeChannel_Y_Explore=ReorganizeData(Eye_Explore_Y);

EyeChannel_X_Explore=EyeChannel_X_Explore(NoNNaNSelector,:);
EyeChannel_Y_Explore=EyeChannel_Y_Explore(NoNNaNSelector,:);
%Using thresholding to select saccades

Saccades=SelectSaccade(EyeChannel_X_Explore,EyeChannel_Y_Explore,EyeBinWidth,V_Threshold,ContinueBin,ISI_Threshold);

SaccadeAngle=arrayfun(@(x)  x.SaccadeAngle,Saccades,'UniformOutput' ,0)';
SaccadeAmplitude=arrayfun(@(x)  x.SaccadeAmplitude,Saccades,'UniformOutput' ,0)';

SaccadeStartTime=arrayfun(@(x)  x.SaccadeStartTime,Saccades,'UniformOutput' ,0)';
SaccadeEndTime=arrayfun(@(x)  x.SaccadeEndTime,Saccades,'UniformOutput' ,0)';

SaccadeStartPoint=arrayfun(@(x)  x.SaccadeStartPoint,Saccades,'UniformOutput' ,0)';
SaccadeEndPoint=arrayfun(@(x)  x.SaccadeEndPoint,Saccades,'UniformOutput' ,0)';

%Screen out the valid saccades and saccade start time after fixoff

AmplitudeThreshold=5;
TimeThreshold_ori=FixOffTime-StartTimeBin; 
TimeThreshold=reshape(TimeThreshold_ori,numel(TimeThreshold_ori),1)-50;
TimeThreshold=TimeThreshold(NoNNaNSelector);
TimeThreshold=num2cell(TimeThreshold);
%TimeThresholdEnd=num2cell(MaxTimeLength-PlotSacPostBuffer);
%{
SaccadeAngle_ITI=fillinNaNs(cellfun(@(x,y,z,m) z(((x>AmplitudeThreshold)&(y<=m)))',SaccadeAmplitude,SaccadeEndTime,SaccadeAngle,TimeThresholdEnd,'UniformOutput' ,0));
SaccadeEndPoint_ITI=(cellfun(@(x,y,z,m) z(((x>AmplitudeThreshold)&(y<=m)),:)',SaccadeAmplitude,SaccadeEndTime,SaccadeEndPoint,TimeThresholdEnd,'UniformOutput' ,0));
SaccadeStartPoint_ITI=(cellfun(@(x,y,z,m) z(((x>AmplitudeThreshold)&(y<=m)),:)',SaccadeAmplitude,SaccadeEndTime,SaccadeStartPoint,TimeThresholdEnd,'UniformOutput' ,0));
%}

SaccadeAngle_Explore=fillinNaNs(cellfun(@(x,y,z,m) z(((x>AmplitudeThreshold)&(y>m)))',SaccadeAmplitude,SaccadeStartTime,SaccadeAngle,TimeThreshold,'UniformOutput' ,0));
SaccadeEndPoint_Explore=(cellfun(@(x,y,z,m) z(((x>AmplitudeThreshold)&(y>m)),:)',SaccadeAmplitude,SaccadeStartTime,SaccadeEndPoint,TimeThreshold,'UniformOutput' ,0));
SaccadeStartPoint_Explore=(cellfun(@(x,y,z,m) z(((x>AmplitudeThreshold)&(y>m)),:)',SaccadeAmplitude,SaccadeStartTime,SaccadeStartPoint,TimeThreshold,'UniformOutput' ,0));
%}


SaccadeNumInEachTrial=cellfun(@(x) size(x,2),SaccadeEndPoint_Explore);

SaccadeEndPoint_Explore_X=NaN*ones(length(SaccadeNumInEachTrial),max(SaccadeNumInEachTrial));
SaccadeEndPoint_Explore_Y=NaN*ones(length(SaccadeNumInEachTrial),max(SaccadeNumInEachTrial));

SaccadeStartPoint_Explore_X=NaN*ones(length(SaccadeNumInEachTrial),max(SaccadeNumInEachTrial));
SaccadeStartPoint_Explore_Y=NaN*ones(length(SaccadeNumInEachTrial),max(SaccadeNumInEachTrial));


SaccadeEndPoint_Explore_X(SaccadeNumInEachTrial>0,:)=fillinNaNs(cellfun(@(x) x(1,:),SaccadeEndPoint_Explore(SaccadeNumInEachTrial>0,:),'UniformOutput' ,0));
SaccadeEndPoint_Explore_Y(SaccadeNumInEachTrial>0,:)=fillinNaNs(cellfun(@(x) x(2,:),SaccadeEndPoint_Explore(SaccadeNumInEachTrial>0,:),'UniformOutput' ,0));

SaccadeStartPoint_Explore_X(SaccadeNumInEachTrial>0,:)=fillinNaNs(cellfun(@(x) x(1,:),SaccadeStartPoint_Explore(SaccadeNumInEachTrial>0,:),'UniformOutput' ,0));
SaccadeStartPoint_Explore_Y(SaccadeNumInEachTrial>0,:)=fillinNaNs(cellfun(@(x) x(2,:),SaccadeStartPoint_Explore(SaccadeNumInEachTrial>0,:),'UniformOutput' ,0));

%{
SaccadeStartTime_Explore=fillinNaNs(cellfun(@(x,y,z,m) z(((x>AmplitudeThreshold)&(y<=m)),:)',SaccadeAmplitude,SaccadeStartTime,SaccadeStartTime,TimeThresholdEnd,'UniformOutput' ,0));
SaccadeEndTime_Explore=fillinNaNs(cellfun(@(x,y,z,m) z(((x>AmplitudeThreshold)&(y<=m)),:)',SaccadeAmplitude,SaccadeStartTime,SaccadeEndTime,TimeThresholdEnd,'UniformOutput' ,0));


SaccadeAmplitude_Explore=fillinNaNs(cellfun(@(x,y,z,m) z(((x>AmplitudeThreshold)&(y<=m)),:)',SaccadeAmplitude,SaccadeStartTime,SaccadeAmplitude,TimeThresholdEnd,'UniformOutput' ,0));
%}
%

SaccadeStartTime_Explore=fillinNaNs(cellfun(@(x,y,z,m) z(((x>AmplitudeThreshold)&(y>m)),:)',SaccadeAmplitude,SaccadeStartTime,SaccadeStartTime,TimeThreshold,'UniformOutput' ,0));
SaccadeEndTime_Explore=fillinNaNs(cellfun(@(x,y,z,m) z(((x>AmplitudeThreshold)&(y>m)),:)',SaccadeAmplitude,SaccadeStartTime,SaccadeEndTime,TimeThreshold,'UniformOutput' ,0));


SaccadeAmplitude_Explore=fillinNaNs(cellfun(@(x,y,z,m) z(((x>AmplitudeThreshold)&(y>m)),:)',SaccadeAmplitude,SaccadeStartTime,SaccadeAmplitude,TimeThreshold,'UniformOutput' ,0));

%BeforeExplore Fixation and Saccade
%FixTime_BeforeExplore_Start



% Loop over the blocks

for i=1:NumberOfBlocks
    
      currBlocks=(BlockType_Line==UniqueBlockType(i));%& ismember(TargetNum_Line,[1]);%& ismember(SceneName_Line,SceneName_Unique(1));%
      FixationTimeStart_currBlocks=SaccadeEndTime_Explore(currBlocks,1:end-1);
      FixationTimeEnd_currBlocks=SaccadeStartTime_Explore(currBlocks,2:end);
      
      SaccadeTime_Explore_Start_currBlocks=SaccadeStartTime_Explore(currBlocks,:);
      SaccadeTime_Explore_End_currBlocks=SaccadeEndTime_Explore(currBlocks,:);
      
      %Remove extra NaN Fixations 
      NaNCol=sum(isnan(FixationTimeEnd_currBlocks),1);
      RemoveCol=NaNCol==size(FixationTimeEnd_currBlocks,1);
      
      FixationTimeStart_currBlocks(:,RemoveCol)=[];
      FixationTimeEnd_currBlocks(:,RemoveCol)=[];
      
      FixationDuration_currBlocks=FixationTimeEnd_currBlocks-FixationTimeStart_currBlocks;
      FixationDuration_currBlocks=reshape(FixationDuration_currBlocks,numel(FixationDuration_currBlocks),1);
      
      %Align PSTH according to the fixation onset time
      SpikeTime_currBlocks=AlignedRaster(currBlocks,:);
      SpikeStart_currBlocks=-PlotStartOffset*ones(sum(currBlocks),1);
      SpikeEnd_currBlocks=RelativeTimeMarker(currBlocks,2);
      %PlotStopOffset*ones(sum(currBlocks),1)+RelativeTimeMarker(currBlocks,2);
      SpikeTime_currBlocks_full=[SpikeTime_currBlocks, SpikeStart_currBlocks,SpikeEnd_currBlocks];
       AlignedSpike=SpikeAlignTool(FixationTimeStart_currBlocks, SpikeTime_currBlocks_full,StepSize,BinWidth);
      
     PSTH_Fixation_Explore=AlignedSpike.PSTH_Aligned;
     TimeSequence_Fixation_Explore=AlignedSpike.PSTH_Time_Aligned;
     SpikeTime_Fixation_Explore_currBlocks=AlignedSpike.SpikeTime_Aligned_Whole;
  
    
    %Select PSTH interval
     SelectPSTHinterval_Fixation_Explore=[-PlotFixationOn_StartOffset,PlotFixationOn_StopOffset];
    %SelectPSTHinterval_Fixation_Explore=[-nanmean(nanmean(FixationTimeEnd_currBlocks-FixationTimeStart_currBlocks)),100];
    [PSTH_Fixation_Explore,TimeSequence_Fixation_Explore]=SelectPSTHInterval(PSTH_Fixation_Explore,TimeSequence_Fixation_Explore,SelectPSTHinterval_Fixation_Explore);
    
     %Select spike time interval
     SpikeRaster_FixationOn_Explore_currBlocks=SelectSpikeTimeInterval(SpikeTime_Fixation_Explore_currBlocks,-PlotFixationOn_StartOffset*ones(size(SpikeTime_Fixation_Explore_currBlocks,1),1),FixationDuration_currBlocks+PlotFixationOn_StopOffset);

    %Remove NaNs
    NaN_Renove=isnan(FixationDuration_currBlocks);
    PSTH_Fixation_Explore(NaN_Renove,:)=[];
    TimeSequence_Fixation_Explore(NaN_Renove,:)=[];
    SpikeRaster_FixationOn_Explore_currBlocks(NaN_Renove,:)=[];
    %FixationDuration_currBlocks(NaN_Renove,:)=[];
    
    
    
    
    
     %Align PSTH according to the fixation offset time
      AlignedSpike=SpikeAlignTool(FixationTimeEnd_currBlocks, SpikeTime_currBlocks_full,StepSize,BinWidth);
      
     PSTH_FixationEnd_Explore=AlignedSpike.PSTH_Aligned;
     TimeSequence_FixationEnd_Explore=AlignedSpike.PSTH_Time_Aligned;
     SpikeTime_FixationEnd_Explore_currBlocks=AlignedSpike.SpikeTime_Aligned_Whole;
  
    
    %Select PSTH interval
     SelectPSTHinterval_FixationEnd_Explore=[-PlotFixationOff_StartOffset,PlotFixationOff_StopOffset];
    %SelectPSTHinterval_Fixation_Explore=[-nanmean(nanmean(FixationTimeEnd_currBlocks-FixationTimeStart_currBlocks)),100];
    [PSTH_FixationEnd_Explore,TimeSequence_FixationEnd_Explore]=SelectPSTHInterval(PSTH_FixationEnd_Explore,TimeSequence_FixationEnd_Explore,SelectPSTHinterval_FixationEnd_Explore);
    
     %Select spike time interval
     SpikeRaster_FixationEnd_Explore_currBlocks=SelectSpikeTimeInterval(SpikeTime_FixationEnd_Explore_currBlocks,-FixationDuration_currBlocks-PlotFixationOff_StartOffset,PlotFixationOff_StopOffset*ones(size(SpikeTime_FixationEnd_Explore_currBlocks,1),1));

    %Remove NaNs
   
    PSTH_FixationEnd_Explore(NaN_Renove,:)=[];
    TimeSequence_FixationEnd_Explore(NaN_Renove,:)=[];
    SpikeRaster_FixationEnd_Explore_currBlocks(NaN_Renove,:)=[];
  
    FixationDuration_currBlocks(NaN_Renove)=[];
    
    %Align PSTH according to the saccade onset time
    %Remove extra NaN Fixations 
      NaNCol=sum(isnan(SaccadeTime_Explore_Start_currBlocks),1);
      RemoveCol=NaNCol==size(SaccadeTime_Explore_Start_currBlocks,1);
      
      SaccadeTime_Explore_Start_currBlocks(:,RemoveCol)=[];
      SaccadeTime_Explore_End_currBlocks(:,RemoveCol)=[];
      
      SaccadeDuration_currBlocks=SaccadeTime_Explore_Start_currBlocks-SaccadeTime_Explore_End_currBlocks;
      SaccadeDuration_currBlocks=reshape(SaccadeDuration_currBlocks,numel(SaccadeDuration_currBlocks),1);
      
      
       AlignedSpike=SpikeAlignTool(SaccadeTime_Explore_Start_currBlocks, SpikeTime_currBlocks_full,StepSize,BinWidth);
      
     PSTH_Saccade_Explore=AlignedSpike.PSTH_Aligned;
     TimeSequence_Saccade_Explore=AlignedSpike.PSTH_Time_Aligned;
     SpikeTime_Saccade_Explore_currBlocks=AlignedSpike.SpikeTime_Aligned_Whole;
  
    
    %Select PSTH interval
     SelectPSTHinterval_Saccade_Explore=[-PlotSaccadeOn_StartOffset,PlotSaccadeOn_StopOffset];
    %SelectPSTHinterval_Fixation_Explore=[-nanmean(nanmean(FixationTimeEnd_currBlocks-FixationTimeStart_currBlocks)),100];
    [PSTH_Saccade_Explore,TimeSequence_Saccade_Explore]=SelectPSTHInterval(PSTH_Saccade_Explore,TimeSequence_Saccade_Explore,SelectPSTHinterval_Saccade_Explore);
    
     %Select spike time interval
     SpikeRaster_SaccadeOn_Explore_currBlocks=SelectSpikeTimeInterval(SpikeTime_Saccade_Explore_currBlocks,-PlotSaccadeOn_StartOffset*ones(size(SpikeTime_Saccade_Explore_currBlocks,1),1),SaccadeDuration_currBlocks+PlotSaccadeOn_StopOffset);

    %Remove NaNs
    NaN_Renove=isnan(SaccadeDuration_currBlocks);
    PSTH_Saccade_Explore(NaN_Renove,:)=[];
    TimeSequence_Saccade_Explore(NaN_Renove,:)=[];
    SpikeRaster_SaccadeOn_Explore_currBlocks(NaN_Renove,:)=[];
    %FixationDuration_currBlocks(NaN_Renove,:)=[];
    
    
    
    
    
     %Align PSTH according to the saccade offset time
      AlignedSpike=SpikeAlignTool(SaccadeTime_Explore_End_currBlocks, SpikeTime_currBlocks_full,StepSize,BinWidth);
      
     PSTH_SaccadeEnd_Explore=AlignedSpike.PSTH_Aligned;
     TimeSequence_SaccadeEnd_Explore=AlignedSpike.PSTH_Time_Aligned;
     SpikeTime_SaccadeEnd_Explore_currBlocks=AlignedSpike.SpikeTime_Aligned_Whole;
  
    
    %Select PSTH interval
     SelectPSTHinterval_SaccadeEnd_Explore=[-PlotSaccadeOff_StartOffset,PlotSaccadeOff_StopOffset];
    %SelectPSTHinterval_Fixation_Explore=[-nanmean(nanmean(FixationTimeEnd_currBlocks-FixationTimeStart_currBlocks)),100];
    [PSTH_SaccadeEnd_Explore,TimeSequence_SaccadeEnd_Explore]=SelectPSTHInterval(PSTH_SaccadeEnd_Explore,TimeSequence_SaccadeEnd_Explore,SelectPSTHinterval_SaccadeEnd_Explore);
    
     %Select spike time interval
     SpikeRaster_SaccadeEnd_Explore_currBlocks=SelectSpikeTimeInterval(SpikeTime_SaccadeEnd_Explore_currBlocks,-SaccadeDuration_currBlocks-PlotSaccadeOff_StartOffset,PlotSaccadeOff_StopOffset*ones(size(SpikeTime_SaccadeEnd_Explore_currBlocks,1),1));

    %Remove NaNs
   
    PSTH_SaccadeEnd_Explore(NaN_Renove,:)=[];
    TimeSequence_SaccadeEnd_Explore(NaN_Renove,:)=[];
    SpikeRaster_SaccadeEnd_Explore_currBlocks(NaN_Renove,:)=[];
  
    SaccadeDuration_currBlocks(NaN_Renove)=[];
    
    
    
    
    
    
    %{
    
    %Align PSTH according to the before explore fixation onset time
   
    FixationOnset_BeforeExplore_currBlock=FixTime_BeforeExplore_Start(currBlocks,:);
    FixationEnd_BeforeExplore_currBlock=FixTime_BeforeExplore_End(currBlocks,:);
    %Remove extra NaN Fixations 
      NaNCol=sum(isnan(FixationOnset_BeforeExplore_currBlock),1);
      RemoveCol=NaNCol==size(FixationOnset_BeforeExplore_currBlock,1);
      
      FixationOnset_BeforeExplore_currBlock(:,RemoveCol)=[];
      FixationEnd_BeforeExplore_currBlock(:,RemoveCol)=[];
      
      FixationDuration_BeforeExplore_currBlocks=FixationOnset_BeforeExplore_currBlock-FixationEnd_BeforeExplore_currBlock;
      FixationDuration_BeforeExplore_currBlocks=reshape( FixationDuration_BeforeExplore_currBlocks,numel( FixationDuration_BeforeExplore_currBlocks),1);
    
   
   % Align PSTH according to the before explore fixation onset time
    
    
       AlignedSpike=SpikeAlignTool(FixationOnset_BeforeExplore_currBlock, SpikeTime_currBlocks_full,StepSize,BinWidth);
    
     
      
     PSTH_FixationOnset_BeforeExplore=AlignedSpike.PSTH_Aligned;
     TimeSequence_FixationOnset_BeforeExplore=AlignedSpike.PSTH_Time_Aligned;
     SpikeTime_FixationOnset_BeforeExplore_currBlocks=AlignedSpike.SpikeTime_Aligned_Whole;
  
    
    %Select PSTH interval
     SelectPSTHinterval_FixationOn_BeforeExplore=[-PlotFixationOn_StartOffset,PlotFixationOn_StopOffset];
    %SelectPSTHinterval_Fixation_Explore=[-nanmean(nanmean(FixationTimeEnd_currBlocks-FixationTimeStart_currBlocks)),100];
    
    [PSTH_FixationOn_BeforeExplore,TimeSequence_FixationOn_BeforeExplore]=SelectPSTHInterval(PSTH_FixationOnset_BeforeExplore,TimeSequence_FixationOnset_BeforeExplore,SelectPSTHinterval_FixationOn_BeforeExplore);
    
     %Select spike time interval
     SpikeRaster_FixationOn_BeforeExplore_currBlocks=SelectSpikeTimeInterval(SpikeTime_FixationOnset_BeforeExplore_currBlocks,-PlotFixationOn_StartOffset*ones(size(SpikeTime_FixationOnset_BeforeExplore_currBlocks,1),1),  FixationDuration_BeforeExplore_currBlocks+PlotFixationOn_StopOffset);
%

%Remove NaNs
    NaN_Renove=isnan( FixationDuration_BeforeExplore_currBlocks);
    PSTH_FixationOn_BeforeExplore(NaN_Renove,:)=[];
    TimeSequence_FixationOn_BeforeExplore(NaN_Renove,:)=[];
     SpikeRaster_FixationOn_BeforeExplore_currBlocks(NaN_Renove,:)=[];





     %Align PSTH according to the before exploration fixation offset time
      AlignedSpike=SpikeAlignTool(FixationEnd_BeforeExplore_currBlock, SpikeTime_currBlocks_full,StepSize,BinWidth);
      
     PSTH_FixationEnd_BeforeExplore=AlignedSpike.PSTH_Aligned;
     TimeSequence_FixationEnd_BeforeExplore=AlignedSpike.PSTH_Time_Aligned;
     SpikeTime_FixationEnd_BeforeExplore_currBlocks=AlignedSpike.SpikeTime_Aligned_Whole;
  
    
    %Select PSTH interval
     SelectPSTHinterval_FixationEnd_BeforeExplore=[-PlotFixationOff_StartOffset,PlotFixationOff_StopOffset];
    %SelectPSTHinterval_Fixation_Explore=[-nanmean(nanmean(FixationTimeEnd_currBlocks-FixationTimeStart_currBlocks)),100];
     [PSTH_FixationEnd_BeforeExplore,TimeSequence_FixationEnd_BeforeExplore]=SelectPSTHInterval(PSTH_FixationEnd_BeforeExplore,TimeSequence_FixationEnd_BeforeExplore,SelectPSTHinterval_FixationEnd_BeforeExplore);
    
   
     %Select spike time interval
      SpikeRaster_FixationEnd_BeforeExplore_currBlocks=SelectSpikeTimeInterval(SpikeTime_FixationEnd_BeforeExplore_currBlocks,-FixationDuration_BeforeExplore_currBlocks-PlotFixationOff_StartOffset,PlotFixationOff_StopOffset*ones(size(SpikeTime_FixationEnd_BeforeExplore_currBlocks,1),1));

      %Remove NaNs
    
    PSTH_FixationEnd_BeforeExplore(NaN_Renove,:)=[];
    TimeSequence_FixationEnd_BeforeExplore(NaN_Renove,:)=[];
     SpikeRaster_FixationEnd_BeforeExplore_currBlocks(NaN_Renove,:)=[];
     FixationDuration_BeforeExplore_currBlocks(NaN_Renove)=[];
     %}
       %Store useful information for plotting 
       %First Gaze on
    Time_GazeOn_Explore_Block{i}=nanmean(TimeSequence_Fixation_Explore,1);
    PSTH_GazeOn_Explore_mean_Block{i}=nanmean(PSTH_Fixation_Explore,1);
    PSTH_GazeOn_Explore_sem_Block{i}=nanstd(PSTH_Fixation_Explore,[],1)./sqrt(sum(~isnan(PSTH_Fixation_Explore),1));
     
    Raster_GazeOn_Explore_Block{i}=SpikeRaster_FixationOn_Explore_currBlocks;
    GazeDuration_Explore_Block{i}=FixationDuration_currBlocks;
      
      
       %Then Gaze off
    Time_GazeOff_Explore_Block{i}=nanmean(TimeSequence_FixationEnd_Explore,1);
    PSTH_GazeOff_Explore_mean_Block{i}=nanmean(PSTH_FixationEnd_Explore,1);
    PSTH_GazeOff_Explore_sem_Block{i}=nanstd(PSTH_FixationEnd_Explore,[],1)./sqrt(sum(~isnan(PSTH_FixationEnd_Explore),1));
     
    Raster_GazeOff_Explore_Block{i}=SpikeRaster_FixationEnd_Explore_currBlocks;
    
    %{
    %Before exploration gaze onset
    Time_GazeOn_BeforeExplore_Block{i}=nanmean(TimeSequence_FixationOn_BeforeExplore,1);
    PSTH_GazeOn_BeforeExplore_mean_Block{i}=nanmean(PSTH_FixationOn_BeforeExplore,1);
    PSTH_GazeOn_BeforeExplore_sem_Block{i}=nanstd(PSTH_FixationOn_BeforeExplore,[],1)./sqrt(sum(~isnan(PSTH_FixationOn_BeforeExplore),1));
     
    Raster_GazeOn_BeforeExplore_Block{i}=SpikeRaster_FixationOn_BeforeExplore_currBlocks;
    GazeDuration_BeforeExplore_Block{i}= FixationDuration_BeforeExplore_currBlocks;
    
     %Before exploration gaze offset
    Time_GazeOff_BeforeExplore_Block{i}=nanmean(TimeSequence_FixationEnd_BeforeExplore,1);
    PSTH_GazeOff_BeforeExplore_mean_Block{i}=nanmean(PSTH_FixationEnd_BeforeExplore,1);
    PSTH_GazeOff_BeforeExplore_sem_Block{i}=nanstd(PSTH_FixationEnd_BeforeExplore,[],1)./sqrt(sum(~isnan(PSTH_FixationEnd_BeforeExplore),1));
     
    Raster_GazeOff_BeforeExplore_Block{i}=SpikeRaster_FixationEnd_BeforeExplore_currBlocks;
   %}
    %Explore saccade onset
    
    Time_SacOn_Explore_Block{i}=nanmean(TimeSequence_Saccade_Explore,1);
    PSTH_SacOn_Explore_mean_Block{i}=nanmean(PSTH_Saccade_Explore,1);
    PSTH_SacOn_Explore_sem_Block{i}=nanstd(PSTH_Saccade_Explore,[],1)./sqrt(sum(~isnan(PSTH_Saccade_Explore),1));
     
    Raster_SacOn_Explore_Block{i}=SpikeRaster_SaccadeOn_Explore_currBlocks;
    SacDuration_Explore_Block{i}=SaccadeDuration_currBlocks;
      
      
       %Then saccade off
    Time_SacOff_Explore_Block{i}=nanmean(TimeSequence_SaccadeEnd_Explore,1);
    PSTH_SacOff_Explore_mean_Block{i}=nanmean(PSTH_SaccadeEnd_Explore,1);
    PSTH_SacOff_Explore_sem_Block{i}=nanstd(PSTH_SaccadeEnd_Explore,[],1)./sqrt(sum(~isnan(PSTH_SaccadeEnd_Explore),1));
     
    Raster_SacOff_Explore_Block{i}=SpikeRaster_SaccadeEnd_Explore_currBlocks;
    
    %{
    %Re-order the raster according to the fixation duration
   [ FixDur_Order,OrderIdx]=sort(FixationDuration__BeforeExplore_currBlocks);
    SpikeRaster_FixationEnd_BeforeExplore_currBlocks_Ordered= SpikeRaster_FixationEnd_BeforeExplore_currBlocks(OrderIdx,:);
    
    subplot(ceil(NumberOfBlocks/2),2,i)
    SpikeRaterPlot( SpikeRaster_FixationEnd_BeforeExplore_currBlocks_Ordered,1,1,ColorForEachMode(AssociateBlock==UniqueBlockType(i),:),1);
   % subplot(ceil(NumberOfBlocks/2)*2,2,i)
  %  plot(nanmean(TimeSequence_FixationEnd_Explore),nanmean(PSTH_FixationEnd_Explore));
   % hold on
   % errorbar(nanmean(TimeSequence_Fixation_Explore),nanmean(PSTH_Fixation_Explore),nanstd(PSTH_Fixation_Explore)/sqrt(sum(currBlocks)));
   
   
   
    
    
   %  xlim( SelectPSTHinterval_Fixation_Explore);
   %xlim([-nanmean(nanmean(FixationTimeEnd_currBlocks-FixationTimeStart_currBlocks)),100]);
    %}
end


if ShowFigureFlag
    FigureStartNum=130;
    FigureIndex=1;
    LineWidth=2;
    Barlength=1;
    if PlotRasterFlag
    %Figure1, Raster of gaze on/off response across blocks
    figtitlestr{FigureIndex}='RasterGazeOnOff_Explore';
    fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[50,100, 800,800],'Name',figtitlestr{FigureIndex});
    StartLocOn=1;
    StartLocOff=1;
    
    for i=1:NumberOfBlocks
        currBlock=UniqueBlockType(i);
        duration=GazeDuration_Explore_Block{i};
        %Re-order the raster according to the fixation duration
       [ FixDur_Order,OrderIdx]=sort(duration,'descend');
   
        subplot(1,2,1); %First, raster of gaze on
        
        raster_on_reorder=Raster_GazeOn_Explore_Block{i}(OrderIdx,:);
        endmarker=FixDur_Order;
        NumTrials=length(endmarker);
        
        %Event-on line
            plot([0,0],[StartLocOn-1,StartLocOn+NumTrials],'--k','Color','#929498','LineWidth',1);
            hold on
            %Spike raster
            SpikeRaterPlot(raster_on_reorder,StartLocOn,Barlength,ColorForEachMode(AssociateBlock==currBlock,:),LineWidth);
            
            %Add event-end marker
            y_tmp=[StartLocOn:StartLocOn+NumTrials-1]';
            x_tmp=FixDur_Order;
            plot(x_tmp,y_tmp,'>','MarkerEdgeColor','#707276','MarkerFaceColor','#707276','MarkerSize',2);
            xlim([-PlotStartOffset,600]);
            xlabel('Time from gaze onset');
            ylabel('Trial #');
            title('Gaze On during explore');
            set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
          %  legend(Blockstr{AssociateBlock==currBlock});
        
            StartLocOn=StartLocOn+NumTrials;
            
           subplot(1,2,2); %Then, raster of gaze off
        
        raster_off_reorder=Raster_GazeOff_Explore_Block{i}(OrderIdx,:);
        
        
        %Event-on line
            plot([0,0],[StartLocOff-1,StartLocOff+NumTrials],'--k','Color','#929498','LineWidth',1);
            hold on
            %Spike raster
            SpikeRaterPlot(raster_off_reorder,StartLocOff,Barlength,ColorForEachMode(AssociateBlock==currBlock,:),LineWidth);
            
            %Add event-end marker
            y_tmp=[StartLocOff:StartLocOff+NumTrials-1]';
            x_tmp=-FixDur_Order;
            plot(x_tmp,y_tmp,'>','MarkerEdgeColor','#707276','MarkerFaceColor','#707276','MarkerSize',2);
            xlim([-600,PlotStopOffset]);
            xlabel('Time from gaze offset');
            ylabel('Trial #');
            title('Gaze Off during explore');
            set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
          %  legend(Blockstr{AssociateBlock==currBlock});
        
            StartLocOff=StartLocOff+NumTrials;
               
            
    end     
        
   
     FigureStartNum=FigureStartNum+1;
     FigureIndex=FigureIndex+1;
        
    end %End of the if plot raster flag
    
    
  %Figure2, SpikeDensity of gaze on/off response across blocks
    figtitlestr{FigureIndex}='SDFGazeOnOff_Explore';
    fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[70,100, 800,400],'Name',figtitlestr{FigureIndex});
     
        
         
        for i=1:NumberOfBlocks
            
            %Spike Density  PSTH
           
            currBlock=UniqueBlockType(i);
            BlockOrder(i)=find(AssociateBlock==UniqueBlockType(i));
            %For Gaze onset
            subplot(1,2,1);
            
            
            curr_Response=PSTH_GazeOn_Explore_mean_Block{i};
            curr_Response_err=PSTH_GazeOn_Explore_sem_Block{i};
            
            hold on
            shadedErrorBar(Time_GazeOn_Explore_Block{i},curr_Response,curr_Response_err,'lineprops',{ColorForEachMode(AssociateBlock==currBlock,:)},'transparent',1,'patchSaturation',0.5);%PSTH for fix
            set(findall(gca, 'Type', 'Line'),'LineWidth',3);
            hold on 
           
            
            xlim([-PlotFixationOn_StartOffset,PlotFixationOn_StopOffset]);
            xlabel('Time from gaze onset');
            ylabel('Firing Rate(Hz)');
           
             title('Gaze Onset');
            set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
            
            
            
            %For Gaze offset
            subplot(1,2,2);
            
            %Event-on line
            curr_Response=PSTH_GazeOff_Explore_mean_Block{i};
            curr_Response_err=PSTH_GazeOff_Explore_sem_Block{i};
          
            hold on
            shadedErrorBar(Time_GazeOff_Explore_Block{i},curr_Response,curr_Response_err,'lineprops',{ColorForEachMode(AssociateBlock==currBlock,:)},'transparent',1,'patchSaturation',0.5);%PSTH for fix
            set(findall(gca, 'Type', 'Line'),'LineWidth',3);
            hold on 
           
            
            xlim([-PlotFixationOff_StartOffset,PlotFixationOff_StopOffset]);
            xlabel('Time from gaze offset');
            title('Gaze Offset');
            set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
            
            %}
            
        end
        % ylim([0,max(  Interval_max)]);
         subplot(1,2,1);
         %Event-on line
         yrange=ylim;
          plot([0,0],[yrange(1),yrange(2)],'--k','Color','#929498','LineWidth',1);
          subplot(1,2,2);
         %Event-on line
         yrange=ylim;
            plot([0,0],[yrange(1),yrange(2)],'--k','Color','#929498','LineWidth',1);
            legend(Blockstr{ BlockOrder});
    
    FigureIndex=FigureIndex+1;
    FigureStartNum=FigureStartNum+1; 
    
    %{
    %Figure3, SpikeDensity of before explore gaze on/off response across blocks
    figtitlestr{FigureIndex}='SDFGazeOnOff_BeforeExplore';
    fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[90,100, 800,400],'Name',figtitlestr{FigureIndex});
     
        
         
        for i=1:NumberOfBlocks
            
            %Spike Density  PSTH
           
            currBlock=UniqueBlockType(i);
            BlockOrder(i)=find(AssociateBlock==UniqueBlockType(i));
            %For Gaze onset
            subplot(1,2,1);
            
            
            curr_Response=PSTH_GazeOn_BeforeExplore_mean_Block{i};
            curr_Response_err=PSTH_GazeOn_BeforeExplore_sem_Block{i};
            
            hold on
            shadedErrorBar(Time_GazeOn_BeforeExplore_Block{i},curr_Response,curr_Response_err,'lineprops',{ColorForEachMode(AssociateBlock==currBlock,:)},'transparent',1,'patchSaturation',0.5);
            set(findall(gca, 'Type', 'Line'),'LineWidth',3);
            hold on 
           
            
            xlim([-PlotFixationOn_StartOffset,PlotFixationOn_StopOffset]);
            xlabel('Time from gaze onset');
            ylabel('Firing Rate(Hz)');
           
             title('Gaze Onset');
            set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
            
            
            
            %For Gaze offset
            subplot(1,2,2);
            
            %Event-on line
            curr_Response=PSTH_GazeOff_BeforeExplore_mean_Block{i};
            curr_Response_err=PSTH_GazeOff_BeforeExplore_sem_Block{i};
          
            hold on
            shadedErrorBar(Time_GazeOff_BeforeExplore_Block{i},curr_Response,curr_Response_err,'lineprops',{ColorForEachMode(AssociateBlock==currBlock,:)},'transparent',1,'patchSaturation',0.5);%PSTH for fix
            set(findall(gca, 'Type', 'Line'),'LineWidth',3);
            hold on 
           
            
            xlim([-PlotFixationOff_StartOffset,PlotFixationOff_StopOffset]);
            xlabel('Time from gaze offset');
            title('Gaze Offset');
            set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
            
            
            
        end
        % ylim([0,max(  Interval_max)]);
         subplot(1,2,1);
         %Event-on line
         yrange=ylim;
          plot([0,0],[yrange(1),yrange(2)],'--k','Color','#929498','LineWidth',1);
          subplot(1,2,2);
         %Event-on line
         yrange=ylim;
            plot([0,0],[yrange(1),yrange(2)],'--k','Color','#929498','LineWidth',1);
            legend(Blockstr{ BlockOrder});
    
    FigureIndex=FigureIndex+1;
    FigureStartNum=FigureStartNum+1; 
    %}
    
   %Figure4, SpikeDensity of sac on/off response across blocks
    figtitlestr{FigureIndex}='SDFSaccadeOnOff_Explore';
    fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[110,100, 800,400],'Name',figtitlestr{FigureIndex});
     
        
         
        for i=1:NumberOfBlocks
            
            %Spike Density  PSTH
           
            currBlock=UniqueBlockType(i);
            BlockOrder(i)=find(AssociateBlock==UniqueBlockType(i));
            %For Gaze onset
            subplot(1,2,1);
            
            
            curr_Response=PSTH_SacOn_Explore_mean_Block{i};
            curr_Response_err=PSTH_SacOn_Explore_sem_Block{i};
            
            hold on
            shadedErrorBar(Time_SacOn_Explore_Block{i},curr_Response,curr_Response_err,'lineprops',{ColorForEachMode(AssociateBlock==currBlock,:)},'transparent',1,'patchSaturation',0.5);%PSTH for fix
            set(findall(gca, 'Type', 'Line'),'LineWidth',3);
            hold on 
           
            
            xlim([-PlotSaccadeOn_StartOffset,PlotSaccadeOn_StopOffset]);
            xlabel('Time from sac onset');
            ylabel('Firing Rate(Hz)');
           
             title('Saccade Onset');
            set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
            
            
            
            %For Saccade offset
            subplot(1,2,2);
            
            %Event-on line
            curr_Response=PSTH_SacOff_Explore_mean_Block{i};
            curr_Response_err=PSTH_SacOff_Explore_sem_Block{i};
          
            hold on
            shadedErrorBar(Time_SacOff_Explore_Block{i},curr_Response,curr_Response_err,'lineprops',{ColorForEachMode(AssociateBlock==currBlock,:)},'transparent',1,'patchSaturation',0.5);%PSTH for fix
            set(findall(gca, 'Type', 'Line'),'LineWidth',3);
            hold on 
           
            
            xlim([-PlotSaccadeOff_StartOffset,PlotSaccadeOff_StopOffset]);
            xlabel('Time from saccade offset');
            title('Saccade Offset');
            set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
            
            %}
            
        end
        % ylim([0,max(  Interval_max)]);
         subplot(1,2,1);
         %Event-on line
         yrange=ylim;
          plot([0,0],[yrange(1),yrange(2)],'--k','Color','#929498','LineWidth',1);
          subplot(1,2,2);
         %Event-on line
         yrange=ylim;
            plot([0,0],[yrange(1),yrange(2)],'--k','Color','#929498','LineWidth',1);
            legend(Blockstr{ BlockOrder});
    
    FigureIndex=FigureIndex+1;
    FigureStartNum=FigureStartNum+1;   
    

    
    
    
    
end





end
function organized=ReorganizeEye(X,Y);
organized=[];
for i=1:length(X)
    X_curr=X{i};
    Y_curr=Y{i};
    organized=[organized;cellfun(@(x,y) [x;y],X_curr,Y_curr,'uniform',0)];
    
  
    
end


end
function organized=ReorganizeData(data);
organized=[];
for i=1:length(data)
    data_curr=data{i};
    maxnumel=max(cellfun(@numel,data_curr));
    organized{i}=cell2mat(cellfun(@(x) [x,NaN*ones(1,maxnumel-length(x))],data_curr,'uniform',0));
    
end

 maxnumel=max(cellfun(@(x) size(x,2),organized));
  organized=cell2mat(cellfun(@(x) [x,NaN*ones(size(x,1),maxnumel-size(x,2))],organized,'uniform',0)');


end
