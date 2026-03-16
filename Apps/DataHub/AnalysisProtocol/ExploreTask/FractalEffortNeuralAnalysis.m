function FractalEffortNeuralAnalysis(Data, StartTrial, EndTrial,ShowFigureFlag,OutputFlag);
%Protocol to analysis the achieving fractal effect and fractal order
%Xuefei Yu 09142020
ECode;
TaskType=Data.TaskType;
TaskCode=Data.TaskCode;
DataPath= Data.Path;
FileName=Data.FileName;

Selection=StartTrial:EndTrial;%Select the trials according to the pannel
%Abnormal check %Especially For Robin081120FORAGEXPLO4.rst
CheckCode=SCENE_ON;%This code should appear once in each trial
AbnormalTrials=AbnormalTrialCheck(Data,CheckCode);
Selection=Selection(~ismember(Selection,AbnormalTrials));
TrialNum=length(Selection);
%User defined parameter for the neural/behavioral analysis
%TargetOnOffset=100;%in ms, the time selection bin after target onset; 
%StartDelay=-50;%in ms, the time selection bin after fixation offset;

%For Plot of PSTH aligned on Each fixation on and each fixation Off
PlotFixationOn_StartOffset=100;
PlotFixationOn_StopOffset=300;%Duration of plotting

PlotFixationOff_StartOffset=300;
PlotFixationOff_StopOffset=100;%Duration of plotting

PlotSaccadeOn_StartOffset=100;
PlotSaccadeOn_StopOffset=300;

PlotSaccadeOff_StartOffset=300;
PlotSaccadeOff_StopOffset=100;


%For segmentation of  the neural responses 
SegEvent=[FIXOFF,TGTCD];%Event Flow of interest
SegEventEnd=[TGTCD,TGTOFF ];%Event end
SegEventStr={'Explore','Fractal'};

PlotStartOffset=[700,600];%in ms
%PlotStartOffset=[1700];%in ms
PlotStopOffset=[0,0];%in ms

EventTypeFlag=[2,2];%1:Number same with the trial number, 2:Number same with the target number

%Spon Event
SponEvent=[SCENE_OFF];
SponPreTimeBuffer=-500;
SponEventEnd=NaN;
SponStartOffset=0;
SponStopOffset=1000;

%Parameter to define the sliding window to measure firing rates
StepSize=5;%in ms
BinWidth=15;%in ms

%Plot Raster Flag
PlotRasterFlag=0;

%Criterias for selecting saccades
V_Threshold=40;
ContinueBin=5;
ISI_Threshold=20;

%ResponseLatencyThreshold=3;

%Some user defined parameter for visualization
%Color Definition for each saccade nums
Color={'#FB5441','#4173fb','#FBBB41','#0EC335',[0.28,0.78,0.89],'k',[0.6,0.1,1],[0.6,0.1,0.4],[0.3,0.5,1],[0.6,0.5,0.1],[0,0.8,0.4],'r','g','b','m'};
ColorForEachMode=[0.83 0.14 0.14
             1.00 0.54 0.00
             0.47 0.25 0.80
             0.25 0.80 0.54];
AssociateBlock=[10,11,0,1];
Blockstr={'FixLoc,FixSceneSeq','FixLoc,RandSceneSeq','RandLoc,FixSceneSeq','RandLoc,RandSceneSeq'};
%Get File Information
num = regexp(Data.FileName, '\d+', 'match');
RecordDateOriginal=cell2mat(num(1));
RecordDate=datetime(RecordDateOriginal,'InputFormat','MMddyy');
%To get the monkey name
RemainFiles=erase(Data.FileName,RecordDateOriginal);
Delimiter=find(isstrprop(RemainFiles,'upper')==1);
MonkeyName=RemainFiles(Delimiter(1):Delimiter(2)-1);


%BlockSepTrial=unique(SceneLib(:,1));
BlockType=Data.BlockType(:,:);
if RecordDate<datestr('04,03,2021') && TaskCode==FORAGING_EXPLORE_TASK_ITL
    fix_tmp=BlockType(:,3)>0;
    BlockType(fix_tmp,2)=BlockType(fix_tmp,2)-10;
    
end
BlockType=BlockType(Selection,2);


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
if RecordDate<datestr('10,28,2020') && strcmp(MonkeyName,'Robin')
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
end
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
     
 
SearchingTime= (TargetOnTime-FixOffTime);%/1000;%Change into s for presentation


%Segregate spike data from Fixation onset to target onset
SegmentData=EventSegTool(SegEvent,SegEventEnd,PlotStartOffset,PlotStopOffset,EventChannel,EventTimeChannel,SpikeTime,StepSize,BinWidth);

%Retrive data from the field
%Total average for each event
MeanFR=SegmentData.MeanFR; %[Mean,sem]%Averall mean and sem
MeanFRRaw=SegmentData.MeanFRRaw;%=AverageFiringRate_Event;

%PSTH(Spike Density Function) for each event
MeanPSTH=SegmentData.MeanPSTH;%=[PSTH_Event_Time',PSTH_Event'];
MeanPSTH_Time=MeanPSTH{1,1};
MeanPSTH_FR=MeanPSTH{1,2};

%Raster and raw spike count PSTH
AlignedRaster=SegmentData.RasterAfterAlign{1};
PSTH_Event_SpikeCount=SegmentData.PSTH_Event_SpikeCount{1};


%Relative time bin for each event
RelativeTimeMarker=SegmentData.RelativeTimeMarker{1};%=EventTimeMarker;
%Absolute time marker for each event
AbosoluteTimeMarker=SegmentData.AbsoluteTimeMarker{1};

%Get the PSTH after object onset
MeanPSTH_Time_Fractal=MeanPSTH{2,1};
MeanPSTH_FR_Fractal=MeanPSTH{2,2};


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

GoodObj_Line=mod(TargetIndex_Line,2)==1;
BadObj_Line=mod(TargetIndex_Line,2)==0;

TargetPolarAngle_Line=TargetPolarAngle;
TargetPolarAngle_Line=reshape(TargetPolarAngle,numel(TargetPolarAngle),1);
TargetPolarAngle_Line=TargetPolarAngle_Line(NoNNaNSelector);


%TargetOrderLine
TargetOrder=1:max(TargetNum);
TotalNum=length(TargetNum);
TargetOrder_Mat=repmat(TargetOrder,TotalNum,1);
TargetOrder_Line=reshape(TargetOrder_Mat,numel(TargetOrder_Mat),1);
TargetOrder_Line=TargetOrder_Line(NoNNaNSelector);


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

AmplitudeThreshold=1;
TimeThreshold_ori=FixOffTime-StartTimeBin; 
TimeThreshold=reshape(TimeThreshold_ori,numel(TimeThreshold_ori),1);
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

EyeChannel_Ecc_Explore=sqrt(EyeChannel_X_Explore.^2+EyeChannel_Y_Explore.^2);
EyeChannel_V_Explore=abs(diff(EyeChannel_Ecc_Explore,1,2));

SaccadeSep=[1,2,3];
ObjOrderSep=[1,2];
% Loop over the blocks
max_Val=0;
min_Val=999;
for i=1:NumberOfBlocks
    
      currBlocks=(BlockType_Line==UniqueBlockType(i));%&(SaccadeNumInEachTrial==1);%&(SearchingTime_Line<500)%& ismember(TargetNum_Line,[1])& SearchingTime_Line>1000;%& ismember(SceneName_Line,SceneName_Unique(1));%
      
    
     %Separate into different groups
     %%
%I. Sort neural responses according to the Rearching time

%GroupStr={'Highest','Higher Medium','Mid Medium','Lowest'};


Spontmp=Mean_Spon_ForEachTrial_Line(currBlocks);
RearchingTime_currBlock=SearchingTime_Line( currBlocks);
PSTH_currBlock=MeanPSTH_FR(currBlocks,:);%;-Spontmp;
PSTH_Time_currBlock=MeanPSTH_Time(currBlocks,:);
MeanPSTH_Time_currBlock{i}=nanmean(PSTH_Time_currBlock,1);

%RT_First_Target=SacEnd_First_Target;
%SepGroupAmp=3;
%SepGroupAmpNum= SepGroupAmp;
%SepGroupAmp={[1200,600,300]};
SepGroupAmp={[1100]};
%GroupStr={'Highest','Higher Medium','Mid Medium','Lowest'};
GroupStr={'Long search (>500)','Short search (<500)'};
SepGroupAmpNum=length(SepGroupAmp{1})+1;


%After Fractal Onset   
PSTH_Fractal_currBlock=MeanPSTH_FR_Fractal(currBlocks,:);%-Spontmp;
PSTH_Time__Fractal_currBlock=MeanPSTH_Time_Fractal(currBlocks,:);

PSTH_group_temp=SortDataAccordingToGroup(PSTH_Fractal_currBlock,RearchingTime_currBlock,SepGroupAmp,'descend');


PSTH_RTGroup_Fractal_mean{i}=PSTH_group_temp.mean(1:SepGroupAmpNum,:);
PSTH_RTGroup_Fractal_sem{i}=PSTH_group_temp.sem(1:SepGroupAmpNum,:);


PSTH_RTGroup_Fractal_Time{i}=nanmean(PSTH_Time__Fractal_currBlock,1);

%Good Object

%After Fractal Onset   
SelGood=currBlocks&GoodObj_Line;
SelBad=currBlocks&BadObj_Line;
SponGoodtmp=Mean_Spon_ForEachTrial_Line(SelGood);
SponBadtmp=Mean_Spon_ForEachTrial_Line(SelBad);

PSTH_GoodFractal_currBlock=MeanPSTH_FR_Fractal(SelGood,:);%-SponGoodtmp;
PSTH_Time_GoodFractal_currBlock=MeanPSTH_Time_Fractal(SelGood,:);

PSTH_group_temp=SortDataAccordingToGroup(PSTH_GoodFractal_currBlock,SearchingTime_Line(SelGood),SepGroupAmp,'descend');


PSTH_RTGroup_GoodFractal_mean{i}=PSTH_group_temp.mean(1:SepGroupAmpNum,:);
PSTH_RTGroup_GoodFractal_sem{i}=PSTH_group_temp.sem(1:SepGroupAmpNum,:);


PSTH_RTGroup_GoodFractal_Time{i}=nanmean(PSTH_Time_GoodFractal_currBlock,1);
%Bad Object

%After Fractal Onset   
PSTH_BadFractal_currBlock=MeanPSTH_FR_Fractal(SelBad,:);%-SponBadtmp;
PSTH_Time_BadFractal_currBlock=MeanPSTH_Time_Fractal(SelBad,:);

PSTH_group_temp=SortDataAccordingToGroup(PSTH_BadFractal_currBlock,SearchingTime_Line(SelBad),SepGroupAmp,'descend');


PSTH_RTGroup_BadFractal_mean{i}=PSTH_group_temp.mean(1:SepGroupAmpNum,:);
PSTH_RTGroup_BadFractal_sem{i}=PSTH_group_temp.sem(1:SepGroupAmpNum,:);


PSTH_RTGroup_BadFractal_Time{i}=nanmean(PSTH_Time_BadFractal_currBlock,1);
%}
     
    %Separate into different fractal orders
  for j=1:length(ObjOrderSep)+1
          if j==length(ObjOrderSep)+1
             sel=currBlocks & TargetOrder_Line>ObjOrderSep(end);
             order_string{i,j}=sprintf('Fractal Order> %d (N=%d)',ObjOrderSep(end),sum(sel));
          else
            sel=currBlocks & TargetOrder_Line==ObjOrderSep(j);
            order_string{i,j}=sprintf('Fractal Order %d (N=%d)',ObjOrderSep(j),sum(sel));
          end 
      
      %PSTH after researching
      MeanPSTH_Time_Fractal_tmp=MeanPSTH_Time_Fractal(sel,:);
      MeanPSTH_PSTH_Fractal_tmp=MeanPSTH_FR_Fractal(sel,:);%-Mean_Spon_ForEachTrial_Line(sel);
      
      PSTH_currBlock_OrderSep_Fractal_mean{i,j}=nanmean(MeanPSTH_PSTH_Fractal_tmp,1);
      PSTH_currBlock_OrderSep_Fractal_sem{i,j}=nanstd(MeanPSTH_PSTH_Fractal_tmp,[],1)./sqrt(sum(~isnan(MeanPSTH_PSTH_Fractal_tmp),1));
     PSTH_time_currBlock_OrderSep_Fractal_mean{i,j}=nanmean(MeanPSTH_Time_Fractal_tmp,1);
     
      
     %Loop over object types
     GoodSel=sel&GoodObj_Line;
     %PSTH of Good after researching
     Spontmp_GoodSel=Mean_Spon_ForEachTrial_Line(GoodSel);
      MeanPSTH_Time_Fractal_tmp=MeanPSTH_Time_Fractal( GoodSel,:);
      MeanPSTH_PSTH_Fractal_tmp=MeanPSTH_FR_Fractal( GoodSel,:);%-Spontmp_GoodSel;
      
      PSTH_currBlock_OrderSep_GoodFractal_mean{i,j}=nanmean(MeanPSTH_PSTH_Fractal_tmp,1);
      PSTH_currBlock_OrderSep_GoodFractal_sem{i,j}=nanstd(MeanPSTH_PSTH_Fractal_tmp,[],1)./sqrt(sum(~isnan(MeanPSTH_PSTH_Fractal_tmp),1));
     PSTH_time_currBlock_OrderSep_GoodFractal_mean{i,j}=nanmean(MeanPSTH_Time_Fractal_tmp,1);
     
     BadSel=sel&BadObj_Line;
     Spontmp_BadSel=Mean_Spon_ForEachTrial_Line(BadSel);
     %PSTH of Good after researching
      MeanPSTH_Time_Fractal_tmp=MeanPSTH_Time_Fractal( BadSel,:);
      MeanPSTH_PSTH_Fractal_tmp=MeanPSTH_FR_Fractal( BadSel,:);%-Spontmp_BadSel;
      
      PSTH_currBlock_OrderSep_BadFractal_mean{i,j}=nanmean(MeanPSTH_PSTH_Fractal_tmp,1);
      PSTH_currBlock_OrderSep_BadFractal_sem{i,j}=nanstd(MeanPSTH_PSTH_Fractal_tmp,[],1)./sqrt(sum(~isnan(MeanPSTH_PSTH_Fractal_tmp),1));
     PSTH_time_currBlock_OrderSep_BadFractal_mean{i,j}=nanmean(MeanPSTH_Time_Fractal_tmp,1);
     
      if j==length(ObjOrderSep)+1
             
             order_stringGood{i,j}=sprintf('Fractal Order> %d (N=%d)',ObjOrderSep(end),sum(GoodSel));
             order_stringBad{i,j}=sprintf('Fractal Order> %d (N=%d)',ObjOrderSep(end),sum(BadSel));
      else
             order_stringGood{i,j}=sprintf('Fractal Order %d (N=%d)',ObjOrderSep(j),sum(GoodSel));
            order_stringBad{i,j}=sprintf('Fractal Order %d (N=%d)',ObjOrderSep(j),sum(BadSel));
      end
     
     
  end
%{  
  %Move to an independent protocol
  %Measure the location tuning across conditions
  TargetPolarAngle_curr=TargetPolarAngle_Line(currBlocks);
  TargetPolarAngle_curr_unique{i}=unique(TargetPolarAngle_curr);
  for j=1:length(TargetPolarAngle_curr_unique{i})
      ang=TargetPolarAngle_curr_unique{i}(j);
      sel_ang=TargetPolarAngle_curr==ang;
      if sum(sel_ang)>3 %At lease three
           MeanPSTH_currBlock_FractalLocation{i,j}=nanmean(PSTH_Fractal_currBlock(sel_ang,:),1);
           SemPSTH_currBlock_FractalLocation{i,j}=nanstd(PSTH_Fractal_currBlock(sel_ang,:),[],1)./sqrt(sum(~isnan(PSTH_Fractal_currBlock(sel_ang,:)),1));
           
            MeanPSTH_currBlock_BeforeFractalLocation{i,j}=nanmean( PSTH_currBlock(sel_ang,:),1);
            SemPSTH_currBlock_BeforeFractalLocation{i,j}=nanstd( PSTH_currBlock(sel_ang,:),[],1)./sqrt(sum(~isnan( PSTH_currBlock(sel_ang,:)),1));
            
         
      end
      
      
  end
  
 %}
  
end

%PlotStartOffset
if ShowFigureFlag
    FigureStartNum=140;
    FigureIndex=1;
     %Figure6, Response after Fractal Onset for different subsequent
    %object orders for each block
    figtitlestr{FigureIndex}='PSTH_Fractal_DifferentFractalOrder_Blocks';
    fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[50,100, 800,800],'Name',figtitlestr{FigureIndex});
    for i=1:NumberOfBlocks
    
  
        subplot(ceil(NumberOfBlocks/2),2,i)
       
        for j=1:(length(ObjOrderSep)+1)
         
           % Eye_V_Explore_currBlock_SacSep{i,j}
            currBlock=UniqueBlockType(i);
            BlockOrder=find(AssociateBlock==currBlock);
             shadedErrorBar(PSTH_time_currBlock_OrderSep_Fractal_mean{i,j},PSTH_currBlock_OrderSep_Fractal_mean{i,j},PSTH_currBlock_OrderSep_Fractal_sem{i,j},'lineprops',{Color{j}},'transparent',1,'patchSaturation',0.5);%PSTH for fix
             set(findall(gca, 'Type', 'Line'),'LineWidth',3);
             hold on
              xlim([-PlotStartOffset(2),500]);
              %{
              v_tmp=Eye_V_Explore_currBlock_OrderSep{i,j};
              x_tmp=repmat((-PlotStartOffset(1):-PlotStartOffset(1)+size(v_tmp,2)-1),size(v_tmp,1),1);
            % plot(x_tmp',yrange(2)*1.2+v_tmp','-k','color','#9D9EA1');
            h=plot(x_tmp',max_Val*1.2+v_tmp','-k','color',ColorForEachMode(AssociateBlock==currBlock,:));
            %Some legend off
            for ii=1:length(h)
            set( get( get( h(ii), 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );
            end
              %}
 
        end
       
         legend(order_string{i,:});
        
        
       %  ylim([min_Val*0.8,max_Val*1.2+2]);   
         %  ylim([0,150]); 
       % xlim([-PlotStartOffset,RearchingTime_SacSep{i,j}]);
     %  xlim([-PlotStartOffset,RearchingTime_SacSep{1,1}]);
        xlabel('Time from FractalOn(ms)');
        ylabel('Firing Rate(Hz)');
       title(Blockstr{BlockOrder});
         set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
               
            
    end     
        
   
     FigureStartNum=FigureStartNum+1;
     FigureIndex=FigureIndex+1;
     
     
     %Figure2, Response after Different Fractal Type Onset for different subsequent
    %object orders for each block
    figtitlestr{FigureIndex}='PSTH_GooDbADFractal_DifferentFractalOrder_Blocks';
    fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[70,100, 800,800],'Name',figtitlestr{FigureIndex});
    for i=1:NumberOfBlocks
    
  %Good Object
        subplot(NumberOfBlocks,2,(i-1)*2+1)
       
        for j=1:(length(ObjOrderSep)+1)
         
           % Eye_V_Explore_currBlock_SacSep{i,j}
            currBlock=UniqueBlockType(i);
            BlockOrder=find(AssociateBlock==currBlock);
             shadedErrorBar(PSTH_time_currBlock_OrderSep_GoodFractal_mean{i,j},PSTH_currBlock_OrderSep_GoodFractal_mean{i,j},PSTH_currBlock_OrderSep_GoodFractal_sem{i,j},'lineprops',{Color{j}},'transparent',1,'patchSaturation',0.5);%PSTH for fix
             set(findall(gca, 'Type', 'Line'),'LineWidth',3);
             hold on
              xlim([-PlotStartOffset(2),500]);
             
             
 
        end
       
         legend(order_stringGood{i,:});
        
       % ylim([0,150]);
       %  ylim([min_Val*0.8,max_Val*1.2+2]);   
         %  ylim([0,150]); 
       % xlim([-PlotStartOffset,RearchingTime_SacSep{i,j}]);
     %  xlim([-PlotStartOffset,RearchingTime_SacSep{1,1}]);
        xlabel('Time from Good FractalOn(ms)');
        ylabel('Firing Rate(Hz) over spon');
       title(Blockstr{BlockOrder});
         set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
   %Bad Object   
        subplot(NumberOfBlocks,2,(i)*2)
       
        for j=1:(length(ObjOrderSep)+1)
         
           % Eye_V_Explore_currBlock_SacSep{i,j}
            currBlock=UniqueBlockType(i);
            BlockOrder=find(AssociateBlock==currBlock);
             shadedErrorBar(PSTH_time_currBlock_OrderSep_BadFractal_mean{i,j},PSTH_currBlock_OrderSep_BadFractal_mean{i,j},PSTH_currBlock_OrderSep_BadFractal_sem{i,j},'lineprops',{Color{j}},'transparent',1,'patchSaturation',0.5);%PSTH for fix
             set(findall(gca, 'Type', 'Line'),'LineWidth',3);
             hold on
              xlim([-PlotStartOffset(2),600]);
             
             
 
        end
       
         legend(order_stringBad{i,:});
        
        
       %  ylim([min_Val*0.8,max_Val*1.2+2]);   
         %  ylim([0,150]); 
       % xlim([-PlotStartOffset,RearchingTime_SacSep{i,j}]);
     %  xlim([-PlotStartOffset,RearchingTime_SacSep{1,1}]);
        xlabel('Time from Bad FractalOn(ms)');
        ylabel('Firing Rate(Hz)over spon');
       title(Blockstr{BlockOrder});
         set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
        
           % ylim([0,150]); 
    end     
        
   
     FigureStartNum=FigureStartNum+1;
     FigureIndex=FigureIndex+1;
   %Figure3 Response for fractals for different s
    %Searching groups for each mode
    figtitlestr{FigureIndex}='PSTH_Fractal_DifferentRT_Blocks';
    fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[130,100, 800,800],'Name',figtitlestr{FigureIndex});
%    tiledlayout(ceil(NumberOfBlocks/2),2);

     for i=1:NumberOfBlocks
     
         
   subplot(ceil(NumberOfBlocks/2),2,i);
        
     %   ax(i)=nexttile;
       
         for j=1:(SepGroupAmpNum)
           % Eye_V_Explore_currBlock_SacSep{i,j}
            currBlock=UniqueBlockType(i);
            BlockOrder(i)=find(AssociateBlock==currBlock);
            % shadedErrorBar(PSTH_RTGroup_Fractal_Time{i},PSTH_RTGroup_Fractal_mean{i}(j,:),PSTH_RTGroup_Fractal_sem{i}(j,:),'lineprops',{ColorForEachMode(AssociateBlock==currBlock,:)},'transparent',1,'patchSaturation',0.5);%PSTH for fix
              shadedErrorBar(PSTH_RTGroup_Fractal_Time{i},PSTH_RTGroup_Fractal_mean{i}(j,:),PSTH_RTGroup_Fractal_sem{i}(j,:),'lineprops',{Color{j}},'transparent',1,'patchSaturation',0.5);%PSTH for fix
            
            set(findall(gca, 'Type', 'Line'),'LineWidth',3);
             hold on
              xlim([-PlotStartOffset(2),600]);
             %{ 
              v_tmp=Eye_V_Explore_currBlock_SacSep{i,j};
              x_tmp=repmat((-PlotStartOffset:-PlotStartOffset+size(v_tmp,2)-1),size(v_tmp,1),1);
            % plot(x_tmp',yrange(2)*1.2+v_tmp','-k','color','#9D9EA1');
            h=plot(x_tmp',max_Val*1.2+v_tmp','-k','color',ColorForEachMode(AssociateBlock==currBlock,:));
            %Some legend off
            for ii=1:length(h)
            set( get( get( h(ii), 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );
            end
              %}
 
        end
        %if j==1
        
         legend(GroupStr)
        %end
        
        % ylim([min_Val*0.8,max_Val*1.2+2]);   
          %  ylim([0,150]); 
       % xlim([-PlotStartOffset,RearchingTime_SacSep{i,j}]);
     %  xlim([-PlotStartOffset,RearchingTime_SacSep{1,1}]);
        xlabel('Time from Fractal On(ms)');
        ylabel('Firing Rate(Hz)over spon');
        title(Blockstr{BlockOrder(i)});
         set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
               
            
     end 
   %linkaxes(ax,'y');
        
   
     FigureStartNum=FigureStartNum+1;
     FigureIndex=FigureIndex+1;   
     
   %Figure6 Response for different type of fractals for different s
    %Searching groups for each mode
    figtitlestr{FigureIndex}='PSTH_GoodBadFractal_DifferentRT_Blocks';
    fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[150,100, 800,800],'Name',figtitlestr{FigureIndex});
  %  tiledlayout(NumberOfBlocks,2);
     for i=1:NumberOfBlocks
    
         
   %Good fractal 
        
        subplot(NumberOfBlocks,2,(i-1)*2+1)
       %  ax((i-1)*2+1) = nexttile;
       
         for j=1:(SepGroupAmpNum)
           % Eye_V_Explore_currBlock_SacSep{i,j}
            currBlock=UniqueBlockType(i);
            BlockOrder(i)=find(AssociateBlock==currBlock);
             shadedErrorBar(PSTH_RTGroup_GoodFractal_Time{i},PSTH_RTGroup_GoodFractal_mean{i}(j,:),PSTH_RTGroup_GoodFractal_sem{i}(j,:),'lineprops',{Color{j}},'transparent',1,'patchSaturation',0.5);%PSTH for fix
             set(findall(gca, 'Type', 'Line'),'LineWidth',3);
             hold on
              xlim([-PlotStartOffset(2),600]);
             %{ 
              v_tmp=Eye_V_Explore_currBlock_SacSep{i,j};
              x_tmp=repmat((-PlotStartOffset:-PlotStartOffset+size(v_tmp,2)-1),size(v_tmp,1),1);
            % plot(x_tmp',yrange(2)*1.2+v_tmp','-k','color','#9D9EA1');
            h=plot(x_tmp',max_Val*1.2+v_tmp','-k','color',ColorForEachMode(AssociateBlock==currBlock,:));
            %Some legend off
            for ii=1:length(h)
            set( get( get( h(ii), 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );
            end
              %}
 
        end
      %  if j==1
         legend(GroupStr);
         
       % end
        
        % ylim([min_Val*0.8,max_Val*1.2+2]);   
          %  ylim([0,150]); 
       % xlim([-PlotStartOffset,RearchingTime_SacSep{i,j}]);
     %  xlim([-PlotStartOffset,RearchingTime_SacSep{1,1}]);
        xlabel('Time from Good Fractal On(ms)');
        ylabel('Firing Rate(Hz)over spon');
       % title(GroupStr{j});
        title(Blockstr{BlockOrder(i)});
         set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
         
         
        %Bad Object 
      subplot(NumberOfBlocks,2,i*2)
     %   ax(i*2) = nexttile;
       for j=1:(SepGroupAmpNum)
           % Eye_V_Explore_currBlock_SacSep{i,j}
            currBlock=UniqueBlockType(i);
            BlockOrder(i)=find(AssociateBlock==currBlock);
             shadedErrorBar(PSTH_RTGroup_BadFractal_Time{i},PSTH_RTGroup_BadFractal_mean{i}(j,:),PSTH_RTGroup_BadFractal_sem{i}(j,:),'lineprops',{Color{j}},'transparent',1,'patchSaturation',0.5);%PSTH for fix
             set(findall(gca, 'Type', 'Line'),'LineWidth',3);
             hold on
              xlim([-PlotStartOffset(2),600]);
             %{ 
              v_tmp=Eye_V_Explore_currBlock_SacSep{i,j};
              x_tmp=repmat((-PlotStartOffset:-PlotStartOffset+size(v_tmp,2)-1),size(v_tmp,1),1);
            % plot(x_tmp',yrange(2)*1.2+v_tmp','-k','color','#9D9EA1');
            h=plot(x_tmp',max_Val*1.2+v_tmp','-k','color',ColorForEachMode(AssociateBlock==currBlock,:));
            %Some legend off
            for ii=1:length(h)
            set( get( get( h(ii), 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );
            end
              %}
 
        end
        %if j==1
         legend(GroupStr);
        % legend(Blockstr{BlockOrder});
        %end
        
        % ylim([min_Val*0.8,max_Val*1.2+2]);   
          %  ylim([0,150]); 
       % xlim([-PlotStartOffset,RearchingTime_SacSep{i,j}]);
     %  xlim([-PlotStartOffset,RearchingTime_SacSep{1,1}]);
        xlabel('Time from Bad Fractal On(ms)');
        ylabel('Firing Rate(Hz)over spon');
      %  title(GroupStr{j});
        title(Blockstr{BlockOrder(i)});
         set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
                         
            
     end     
      
   %linkaxes(ax,'y');
   
     FigureStartNum=FigureStartNum+1;
     FigureIndex=FigureIndex+1;  
 %{    
       %Figure7 Response for different fractals in each location  for each mode
    figtitlestr{FigureIndex}='PSTH_FractalLocation_Different_Blocks';
    fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[170,100, 800,800],'Name',figtitlestr{FigureIndex});
    tiledlayout(NumberOfBlocks,2);
     for i=1:NumberOfBlocks
    
         
  
       % subplot(NumberOfBlocks,2,(i-1)*2+1)
        ax((i-1)*2+1) = nexttile;
       
       
           % Eye_V_Explore_currBlock_SacSep{i,j}
           
            currBlock=UniqueBlockType(i);
            BlockOrder(i)=find(AssociateBlock==currBlock);
            for j=1:length( TargetPolarAngle_curr_unique{i})
             shadedErrorBar(PSTH_RTGroup_GoodFractal_Time{i}, MeanPSTH_currBlock_FractalLocation{i,j}, SemPSTH_currBlock_FractalLocation{i,j},'lineprops',{Color{j}},'transparent',1,'patchSaturation',0.5);%PSTH for fix
             set(findall(gca, 'Type', 'Line'),'LineWidth',3);
             hold on
              xlim([-PlotStartOffset(2),600]);
             %{ 
              v_tmp=Eye_V_Explore_currBlock_SacSep{i,j};
              x_tmp=repmat((-PlotStartOffset:-PlotStartOffset+size(v_tmp,2)-1),size(v_tmp,1),1);
            % plot(x_tmp',yrange(2)*1.2+v_tmp','-k','color','#9D9EA1');
            h=plot(x_tmp',max_Val*1.2+v_tmp','-k','color',ColorForEachMode(AssociateBlock==currBlock,:));
            %Some legend off
            for ii=1:length(h)
            set( get( get( h(ii), 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );
            end
              %}
 
       
      %  if j==1
       
            end
              legend( string(TargetPolarAngle_curr_unique{i}));
       % end
        
        % ylim([min_Val*0.8,max_Val*1.2+2]);   
          %  ylim([0,150]); 
       % xlim([-PlotStartOffset,RearchingTime_SacSep{i,j}]);
     %  xlim([-PlotStartOffset,RearchingTime_SacSep{1,1}]);
        xlabel('Time from Fractal On(ms)');
        ylabel('Firing Rate(Hz)over spon');
       % title(GroupStr{j});
        title(Blockstr{BlockOrder(i)});
         set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
         
        ax(i*2) = nexttile;
        
          for j=1:length( TargetPolarAngle_curr_unique{i})
             shadedErrorBar(MeanPSTH_Time_currBlock{i},MeanPSTH_currBlock_BeforeFractalLocation{i,j},SemPSTH_currBlock_BeforeFractalLocation{i,j},'lineprops',{Color{j}},'transparent',1,'patchSaturation',0.5);%PSTH for fix
             set(findall(gca, 'Type', 'Line'),'LineWidth',3);
             hold on
              xlim([-PlotStartOffset(1),500]);
          
         end
              legend( string(TargetPolarAngle_curr_unique{i}));
       
        xlabel('Time from Fix off(ms)');
        ylabel('Firing Rate(Hz)over spon');
       % title(GroupStr{j});
        title(Blockstr{BlockOrder(i)});
         set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
         
         
       
     end     
      
   linkaxes(ax,'y');
   
     FigureStartNum=FigureStartNum+1;
     FigureIndex=FigureIndex+1;  
     
     %}
  
   
  
     
     
end   
    
%}

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
