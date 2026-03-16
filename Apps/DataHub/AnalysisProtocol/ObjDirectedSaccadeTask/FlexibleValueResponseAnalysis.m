function FlexibleValueResponseAnalysis(Data, StartTrial, EndTrial,ShowFigureFlag,OutputFlag);
%Function to analysis flexible value task
%Xuefei Yu 01152022

ECode;
SegEvent=[TGTCD];%Event Flow of interest
%SegEventEnd=[SHOWFIXCD,FIXOFF,TGTOFF,RWDOFFCD ,NaN,EXTRAFP];%Event end
SegEventEnd=[TGTOFF];%Event end
SegEventStr={'TargetOn'};
AmpMatchFlag=0;%1:Amplitude match on ;0:Amplitude Match off


PlotStartOffset=[1000];%in ms
PlotStopOffset=[8000];%in ms

ExamplePlot=0;%PLots for presentations, turn off for regular use

ResponseLatencyThreshold=3;

%Parameter to define the sliding window to measure firing rates
StepSize=5;%in ms
BinWidth=20;%in ms


%Criterias for sexwlecting saccades
V_Threshold=40;
ContinueBin=5;
ISI_Threshold=20;

%Time bin pre saccade and time bin offset after saccade onset
PlotSacPreBuffer=300;%in ms
PlotSacPostBuffer=500;%in ms

%PlotFracOffPreBuffer=500;%in ms
%PlotFracOffPostBuffer=2000;%in ms


%Color Definition for each fractals
%{
Color=[lines(6);prism(6) ];
Color(2,:)=copper(1)+0.5;
Color=[Color;Color*0.5];
%}
Color={'#648FFF','#FE6100','#009E73','#88CCEE','#CC6677','#785EF0'};
Color_RGB=[100,143,255;
           254,97,0;
           0,158,115;
          136,204,238;
          204,102,119;
          120,94,240]/255;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Unpack all the data base and task settings
ECode;
FileName=Data.FileName;
TaskType=Data.TaskType;
TaskCode=Data.TaskCode;
TrialNum=EndTrial-StartTrial+1;
Selection=StartTrial:EndTrial;


%Abnormal check
CheckCodePair=[TGTCD,TGTOFF];%This code should be the same as target on for good trials
AbnormalTrials=AbnormalTrialCheck(Data,CheckCodePair);
Selection=Selection(~ismember(Selection,AbnormalTrials));
    





%SpikeChannel=Data.SpikeChannel;



if strcmp(FileName,'Robin040321FRACSAC1.rst') & SpikeChannel==22
    Selection=Selection(~ismember(Selection,[14,40]));
end

%FixationPointIndex=Data.FixPointIndex(Selection,:);
%FixLoc=Data.FixLoc(Selection,:,:);

try
tmp=Data.ParaLib('FlipGoodBad');
FlipGoodBad=tmp(Selection,2);
catch
   FlipGoodBad=0*ones(sum(Selection),1);
   disp('Auto no flipping sign');
end


%Target Information
TargetIndex=Data.TargetIndex(Selection,2);%For this task, there is only target name, but no target index
try
TargetName=Data.TargetIndex(Selection,3);%For this task, there is only target name, but no target index
catch
    TargetName=TargetIndex;
    disp('No targetname');
end
%TargetLoc=Data.TargetLoc(Selection,:,:);

TargetNum=Data.TargetNum(Selection,2:end);

%Rearrange target into lines
TargetIndex_Line=TargetIndex;

TargetIndex_Line=reshape(TargetIndex_Line,numel(TargetIndex_Line),1);

NoNNaNSelector=~isnan(TargetIndex_Line);

TargetIndex_Line=TargetIndex_Line(~isnan(TargetIndex_Line));


TargetName_Line=TargetName;

TargetName_Line=reshape(TargetName_Line,numel(TargetName_Line),1);

TargetName_Line=TargetName_Line(~isnan(TargetName_Line));
%Analysis Good versus Bad

GoodObjectIndex_Raw=TargetIndex_Line<=4;
%%
%Flip the sign of good and bad object according to flipping index

GoodObjectIndex_Raw(FlipGoodBad==1)=~GoodObjectIndex_Raw(FlipGoodBad==1);

%% Target location

TargetLoc=Data.TargetLoc(Selection,2:end);
TargetPolar=Data.TargetPolarLoc(Selection,2:end);
TargetAngle = TargetPolar(:,1);
TargetEcc =  TargetPolar(:,2);

UniqueAngle = unique(TargetAngle);
UniqueEcc = unique(TargetEcc);





%%
%GetEyeData
EyeData=Data.EyeDataRaw;
EyeChannel_X_L=EyeData(Selection,1);
EyeChannel_Y_L=EyeData(Selection,2);
%EyeChannel_X_R=EyeData(Selection,3);
%EyeChannel_Y_R=EyeData(Selection,4);

Eye_Eccentricity=cellfun(@(x,y) sqrt(x.^2+y.^2),EyeChannel_X_L,EyeChannel_Y_L,'uniform',0);

EyeBinWidth=Data.EyeBinWidth;


EventChannel=Data.EventChannel(Selection,:);
EventTimeChannel=Data.EventTimeChannel(Selection,:);
EventBin=Data.EventBin(Selection,:);




%%Reproduce goodbad index from event channel
%ObjectIndex=
ObjectIndex=ReproduceFromEvent(EventChannel,[GOODOBJ  BADOBJ]);

GoodObjectIndex=(ObjectIndex==GOODOBJ)';


%Load the spike channel
SpikeChannelWhole=Data.SpikeChannel;
SpikeTimeWhole=Data.SpikeTimeData(:,:,Selection);
ChannelIDWhole=Data.ChannelID;
ClusterIDWhole=Data.ClusterID;

MultiChannel=0;
if length(SpikeChannelWhole)>1
    MultiChannel=1;
    
end
for spk=1:length(SpikeChannelWhole)
    SpikeChannel=SpikeChannelWhole(spk);
    SpikeTime=squeeze(SpikeTimeWhole(spk,:,:))';
    ChannelID=ChannelIDWhole(spk);
    ClusterID=ClusterIDWhole(spk);

%sSpikeTime=squeeze(Data.SpikeTimeData(:,:,Selection))';
%SpikeBin=squeeze(Data.SpikeBinData(:,:,Selection))';


%Select the period between target onset and target off set
FractalOnsetBin = FindOutTime(EventChannel,EventBin,TGTCD );
FractalOffsetBin = FindOutTime(EventChannel,EventBin,TGTOFF);
FractalOnsetTime = FindOutTime(EventChannel,EventTimeChannel,TGTCD );
FractalOffsetTime = FindOutTime(EventChannel,EventTimeChannel,TGTOFF );

EyeChannel_X_Selection= SelectEyeInterval(EyeChannel_X_L,FractalOnsetBin ,FractalOffsetBin )';
EyeChannel_Y_Selection= SelectEyeInterval(EyeChannel_Y_L,FractalOnsetBin ,FractalOffsetBin )';

EyeChannel_X_Selection=ReorganizeEye(EyeChannel_X_Selection);
EyeChannel_Y_Selection=ReorganizeEye(EyeChannel_Y_Selection);

EyeChannel_X_Selection=EyeChannel_X_Selection(NoNNaNSelector,:);
EyeChannel_Y_Selection=EyeChannel_Y_Selection(NoNNaNSelector,:);





%Analysis saccades for each trial
Saccades=SelectSaccade(EyeChannel_X_Selection,EyeChannel_Y_Selection,EyeBinWidth,V_Threshold,ContinueBin,ISI_Threshold);
SaccadeNumber=arrayfun(@(x)  x.NumOfSaccade,Saccades)';
SaccadeAngle=arrayfun(@(x)  x.SaccadeAngle,Saccades,'UniformOutput' ,0)';
SaccadeAmplitude=arrayfun(@(x)  x.SaccadeAmplitude,Saccades,'UniformOutput' ,0)';

SaccadeStartTime=arrayfun(@(x)  x.SaccadeStartTime,Saccades,'UniformOutput' ,0)';
SaccadeEndTime=arrayfun(@(x)  x.SaccadeEndTime,Saccades,'UniformOutput' ,0)';

SaccadeStartPoint=arrayfun(@(x)  x.SaccadeStartPoint,Saccades,'UniformOutput' ,0)';
SaccadeEndPoint=arrayfun(@(x)  x.SaccadeEndPoint,Saccades,'UniformOutput' ,0)';

%Screen out the valid saccades and saccade start time above 5;

AmplitudeThreshold=3;

TimeThreshold=50;

SaccadeAngle_Selection=cellfun(@(x,y,z) z(((x>AmplitudeThreshold)&(y>TimeThreshold))),SaccadeAmplitude,SaccadeStartTime,SaccadeAngle,'UniformOutput' ,0);
SaccadeEndPoint_Selection=cellfun(@(x,y,z) z(((x>AmplitudeThreshold)&(y>TimeThreshold)),:),SaccadeAmplitude,SaccadeStartTime,SaccadeEndPoint,'UniformOutput' ,0);
SaccadeStartPoint_Selection=cellfun(@(x,y,z) z(((x>AmplitudeThreshold)&(y>TimeThreshold)),:),SaccadeAmplitude,SaccadeStartTime,SaccadeStartPoint,'UniformOutput' ,0);

SaccadeStartTime_Selection=cellfun(@(x,y,z) z(((x>AmplitudeThreshold)&(y>TimeThreshold)),:),SaccadeAmplitude,SaccadeStartTime,SaccadeStartTime,'UniformOutput' ,0);
SaccadeEndTime_Selection=cellfun(@(x,y,z) z(((x>AmplitudeThreshold)&(y>TimeThreshold)),:),SaccadeAmplitude,SaccadeStartTime,SaccadeEndTime,'UniformOutput' ,0);


SaccadeAmplitude_Selection=cellfun(@(x,y,z) z(((x>AmplitudeThreshold)&(y>TimeThreshold)),:),SaccadeAmplitude,SaccadeStartTime,SaccadeAmplitude,'UniformOutput' ,0);
%Only look at the first saccade
SaccadeNumInEachTrial=cellfun(@numel,SaccadeAngle_Selection);


SaccadeAngle_First=NaN*ones(1,length(SaccadeNumInEachTrial));
SaccadeEndPoint_First=NaN*ones(length(SaccadeNumInEachTrial),2);
RT_First_Target=NaN*ones(1,length(SaccadeNumInEachTrial));
SacEnd_First_Target=NaN*ones(1,length(SaccadeNumInEachTrial));

SaccadeStartPoint_First=NaN*ones(length(SaccadeNumInEachTrial),2);

SaccadeAmplitude_First=NaN*ones(1,length(SaccadeNumInEachTrial));

SaccadeAngle_First(SaccadeNumInEachTrial>0)=cell2mat(cellfun(@(x) x(1,:),SaccadeAngle_Selection(SaccadeNumInEachTrial>0),'UniformOutput' ,0));
SaccadeEndPoint_First(SaccadeNumInEachTrial>0,:)=cell2mat(cellfun(@(x) x(1,:),SaccadeEndPoint_Selection(SaccadeNumInEachTrial>0,:),'UniformOutput' ,0));
RT_First_Target(SaccadeNumInEachTrial>0)=cell2mat(cellfun(@(x) x(1,:),SaccadeStartTime_Selection(SaccadeNumInEachTrial>0,:),'UniformOutput' ,0));
SacEnd_First_Target(SaccadeNumInEachTrial>0)=cell2mat(cellfun(@(x) x(1,:),SaccadeEndTime_Selection(SaccadeNumInEachTrial>0,:),'UniformOutput' ,0));


SaccadeStartPoint_First(SaccadeNumInEachTrial>0,:)=cell2mat(cellfun(@(x) x(1,:),SaccadeStartPoint_Selection(SaccadeNumInEachTrial>0,:),'UniformOutput' ,0));

SaccadeAmplitude_First(SaccadeNumInEachTrial>0)=cell2mat(cellfun(@(x) x(1,:),SaccadeAmplitude_Selection(SaccadeNumInEachTrial>0,:),'UniformOutput' ,0));




%Align neural responses according to the saccade onset time
%

%Organize the spike time using the same way
%Align spike time with target onset

AlignedSpike=SpikeAlignTool(FractalOnsetTime,SpikeTime,StepSize,BinWidth);

SpikeTimeAligned_TgtOn=AlignedSpike.SpikeTime_Aligned_Whole;


SpikeTimeAligned_TgtOn=SpikeTimeAligned_TgtOn(NoNNaNSelector,:);

%Align again on saccade onset
AlignedSpike=SpikeAlignTool(RT_First_Target',SpikeTimeAligned_TgtOn,StepSize,BinWidth);

PSTH_SacOn=AlignedSpike.PSTH_Aligned;
TimeSequence_SacOn=AlignedSpike.PSTH_Time_Aligned;
%Select interval
SelectPSTHinterval_SacOn=[-PlotSacPreBuffer,PlotSacPostBuffer];
[PSTH_SacOn,TimeSequence_SacOn]=SelectPSTHInterval(PSTH_SacOn,TimeSequence_SacOn,SelectPSTHinterval_SacOn);



TimeSequence_SacOn_mean=nanmean(TimeSequence_SacOn,1);
PSTH_SacOn_mean=nanmean(PSTH_SacOn,1);
PSTH_SacOn_sem=nanstd(PSTH_SacOn,[],1)./sqrt(sum(~isnan(PSTH_SacOn),1));

%Good versus Bad
PSTH_SacOn_Good_mean=nanmean(PSTH_SacOn(GoodObjectIndex,:),1);
PSTH_SacOn_Good_sem=nanstd(PSTH_SacOn(GoodObjectIndex,:),[],1)./sqrt(sum(~isnan(PSTH_SacOn(GoodObjectIndex,1)),1));

PSTH_SacOn_Bad_mean=nanmean(PSTH_SacOn(~GoodObjectIndex,:),1);
PSTH_SacOn_Bad_sem=nanstd(PSTH_SacOn(~GoodObjectIndex,:),[],1)./sqrt(sum(~isnan(PSTH_SacOn(~GoodObjectIndex,1)),1));


%Align Eye trace with saccade onset
%Select the eye eccentricity between target onset and target offset
Eye_Eccentricity_Tgt=SelectEyeInterval(Eye_Eccentricity,FractalOnsetBin,FractalOffsetBin)';
Eye_Eccentricity_Tgt=ReorganizeEye(Eye_Eccentricity_Tgt);
Eye_Eccentricity_Tgt=Eye_Eccentricity_Tgt(NoNNaNSelector,:);

%Select the interval around saccade onset
EyeTime=0:EyeBinWidth:size(Eye_Eccentricity_Tgt,2)-1;

EyeEcc_Aligned=AlignEyeData(Eye_Eccentricity_Tgt,EyeTime,RT_First_Target,[-PlotSacPreBuffer,PlotSacPostBuffer]);



EyeEcc_Aligned_mean=nanmean(EyeEcc_Aligned,1);
EyeEcc_Aligned_sem=nanstd(EyeEcc_Aligned,[],1)./sqrt(sum(~isnan(EyeEcc_Aligned),1));

%Good versus Bad

EyeEcc_Aligned_Good_mean=nanmean(EyeEcc_Aligned(GoodObjectIndex,:),1);
EyeEcc_Aligned_Good_sem=nanstd(EyeEcc_Aligned(GoodObjectIndex,:),[],1)./sqrt(sum(~isnan(EyeEcc_Aligned(GoodObjectIndex,:)),1));

EyeEcc_Aligned_Bad_mean=nanmean(EyeEcc_Aligned(~GoodObjectIndex,:),1);
EyeEcc_Aligned_Bad_sem=nanstd(EyeEcc_Aligned(~GoodObjectIndex,:),[],1)./sqrt(sum(~isnan(EyeEcc_Aligned(~GoodObjectIndex,:)),1));

%%
%%%%%%%%%%%%%%%%%%%%%%%Align data on object onset
%Segregate spike data from fractal onset to fractal offset

SegmentData=EventSegTool(SegEvent,SegEventEnd,PlotStartOffset,PlotStopOffset,EventChannel,EventTimeChannel,SpikeTime,StepSize,BinWidth);

%Retrive data from the field

%Relative time bin for each event
RelativeTimeMarker=SegmentData.RelativeTimeMarker{1};%=EventTimeMarker;
%Absolute time marker for each event
AbosoluteTimeMarker=SegmentData.AbsoluteTimeMarker{1};

ObjOnTimeMarker=nanmean(RelativeTimeMarker);

RwdTimeAbs=FindOutTime(EventChannel,EventBin,RWDCD);
RwdTimeRel=RwdTimeAbs-AbosoluteTimeMarker(:,1);

%Algin marker to sac-on & sac-off
RwdTimeRelSacEnd=RwdTimeRel-SacEnd_First_Target';
TgtOffTimeRelSacEnd=RelativeTimeMarker(:,2)-SacEnd_First_Target';

RwdTimeRelSacEndMean=nanmean(RwdTimeRelSacEnd);
TgtOffTimeRelSacEndMean=nanmean(TgtOffTimeRelSacEnd);

%PSTH for each event
MeanPSTH=SegmentData.MeanPSTH;%=[PSTH_Event_Time',PSTH_Event'];
MeanPSTH_Time=MeanPSTH{1};
MeanPSTH_FR_FractalOn=MeanPSTH{2};

PSTH_Time_FraOn_All=nanmean(MeanPSTH_Time,1);

%Get mean and sem of neural responses
MeanPSTH_FR_FractalOn_mean=nanmean(MeanPSTH_FR_FractalOn,1);
MeanPSTH_FR_FractalOn_sem=nanstd(MeanPSTH_FR_FractalOn,[],1)./sqrt(sum(~isnan(MeanPSTH_FR_FractalOn),1));
%Good vs.Bad
PSTH_FR_Good=MeanPSTH_FR_FractalOn(GoodObjectIndex,:);
PSTH_FR_Bad=MeanPSTH_FR_FractalOn(~GoodObjectIndex,:);

PSTH_FR_Good_mean=nanmean(PSTH_FR_Good,1);
PSTH_FR_Bad_mean=nanmean(PSTH_FR_Bad,1);

PSTH_FR_Good_sem=nanstd(PSTH_FR_Good,[],1)./sqrt(sum(~isnan(PSTH_FR_Good),1));
PSTH_FR_Bad_sem=nanstd(PSTH_FR_Bad,[],1)./sqrt(sum(~isnan(PSTH_FR_Bad),1));

%{
%Get mean and sem of eye eccentricity
%Get eye trance distance during the interval
StartTimeBin=FindOutTime(EventChannel,EventTimeChannel,SegEvent)-PlotStartOffset;
StartTimeBin(StartTimeBin<1)=1;
EndTimeBin=FindOutTime(EventChannel,EventTimeChannel,SegEventEnd )+PlotStopOffset;
  
Eye_Eccentricity_FractalOn=SelectEyeInterval(Eye_Eccentricity,StartTimeBin,EndTimeBin)';
Eye_Eccentricity_FractalOn=ReorganizeEye(Eye_Eccentricity_FractalOn);
Eye_Eccentricity_FractalOn=Eye_Eccentricity_FractalOn(NoNNaNSelector,:);





Eye_Eccentricity_FractalOn_mean=nanmean(Eye_Eccentricity_FractalOn,1);
Eye_Eccentricity_FractalOn_sem=nanstd(Eye_Eccentricity_FractalOn,1)./sqrt(sum(~isnan(Eye_Eccentricity_FractalOn),1));

%Eye eccentricity Good v.s. Bad
Eye_Eccentricity_FractalOn_Good=Eye_Eccentricity_FractalOn(GoodObjectIndex,:);
Eye_Eccentricity_FractalOn_Bad=Eye_Eccentricity_FractalOn(~GoodObjectIndex,:);


Eye_Ecc_Good_mean=nanmean(Eye_Eccentricity_FractalOn_Good,1);
Eye_Ecc_Bad_mean=nanmean(Eye_Eccentricity_FractalOn_Bad,1);

Eye_Ecc_Good_sem=nanstd(Eye_Eccentricity_FractalOn_Good,[],1)./sqrt(sum(~isnan(Eye_Eccentricity_FractalOn_Good),1));
Eye_Ecc_Bad_sem=nanstd(Eye_Eccentricity_FractalOn_Bad,[],1)./sqrt(sum(~isnan(Eye_Eccentricity_FractalOn_Bad),1));
%}
%{
%%
% %%%%%%%%%%%%%%%%%%%%%%%Align data on saccade offset

AlignedSpike=SpikeAlignTool(SacEnd_First_Target',SpikeTimeAligned_TgtOn,StepSize,BinWidth);

PSTH_SacOff=AlignedSpike.PSTH_Aligned;
TimeSequence_SacOff=AlignedSpike.PSTH_Time_Aligned;

%Select interval
SelectPSTHinterval_SacOff=[-PlotSacPreBuffer,PlotSacPostBuffer];
[PSTH_SacOff,TimeSequence_SacOff]=SelectPSTHInterval(PSTH_SacOff,TimeSequence_SacOff,SelectPSTHinterval_SacOff);



TimeSequence_SacOff_mean=nanmean(TimeSequence_SacOff,1);
PSTH_SacOff_mean=nanmean(PSTH_SacOff,1);
PSTH_SacOff_sem=nanstd(PSTH_SacOff,[],1)./sqrt(sum(~isnan(PSTH_SacOff),1));

%Good versus Bad
PSTH_SacOff_Good_mean=nanmean(PSTH_SacOff(GoodObjectIndex,:),1);
PSTH_SacOff_Good_sem=nanstd(PSTH_SacOff(GoodObjectIndex,:),[],1)./sqrt(sum(~isnan(PSTH_SacOff(GoodObjectIndex,1)),1));

PSTH_SacOff_Bad_mean=nanmean(PSTH_SacOff(~GoodObjectIndex,:),1);
PSTH_SacOff_Bad_sem=nanstd(PSTH_SacOff(~GoodObjectIndex,:),[],1)./sqrt(sum(~isnan(PSTH_SacOff(~GoodObjectIndex,1)),1));

%Align Eye trace with saccade offset

EyeEcc_Aligned_SacOff=AlignEyeData(Eye_Eccentricity_Tgt,EyeTime,SacEnd_First_Target,[-PlotSacPreBuffer,PlotSacPostBuffer]);

EyeEcc_Aligned_SacOff_mean=nanmean(EyeEcc_Aligned_SacOff,1);
EyeEcc_Aligned_SacOff_sem=nanstd(EyeEcc_Aligned_SacOff,[],1)./sqrt(sum(~isnan(EyeEcc_Aligned_SacOff),1));

%Good versus Bad

EyeEcc_Aligned_SacOff_Good_mean=nanmean(EyeEcc_Aligned_SacOff(GoodObjectIndex,:),1);
EyeEcc_Aligned_SacOff_Good_sem=nanstd(EyeEcc_Aligned_SacOff(GoodObjectIndex,:),[],1)./sqrt(sum(~isnan(EyeEcc_Aligned_SacOff(GoodObjectIndex,:)),1));

EyeEcc_Aligned_SacOff_Bad_mean=nanmean(EyeEcc_Aligned_SacOff(~GoodObjectIndex,:),1);
EyeEcc_Aligned_SacOff_Bad_sem=nanstd(EyeEcc_Aligned_SacOff(~GoodObjectIndex,:),[],1)./sqrt(sum(~isnan(EyeEcc_Aligned_SacOff(~GoodObjectIndex,:)),1));

%%
%}
%{
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Align eye trace on target offset
EyeTimeII=[-PlotStartOffset:-PlotStartOffset+size(Eye_Eccentricity_FractalOn,2)-1];

EyeEcc_Aligned_FraOff=AlignEyeData(Eye_Eccentricity_FractalOn,EyeTimeII,RelativeTimeMarker(:,2),[-PlotFracOffPreBuffer,PlotFracOffPostBuffer]);

%Good versus Bad

EyeEcc_Aligned_FracOff_Good_mean=nanmean(EyeEcc_Aligned_FraOff(GoodObjectIndex,:),1);
EyeEcc_Aligned_FracOff_Good_sem=nanstd(EyeEcc_Aligned_FraOff(GoodObjectIndex,:),[],1)./sqrt(sum(~isnan(EyeEcc_Aligned_FraOff(GoodObjectIndex,:)),1));

EyeEcc_Aligned_FracOff_Bad_mean=nanmean(EyeEcc_Aligned_FraOff(~GoodObjectIndex,:),1);
EyeEcc_Aligned_FracOff_Bad_sem=nanstd(EyeEcc_Aligned_FraOff(~GoodObjectIndex,:),[],1)./sqrt(sum(~isnan(EyeEcc_Aligned_FraOff(~GoodObjectIndex,:)),1));


%Get the neural data aligned on frac-off

SegEventOff=[TGTOFF];%Event Flow of interest
%SegEventEnd=[SHOWFIXCD,FIXOFF,TGTOFF,RWDOFFCD ,NaN,EXTRAFP];%Event end
SegEventEndOff=[NaN];%Event end
SegEventStrOff={'TargetOff'};


SegmentData=EventSegTool(SegEventOff,SegEventEndOff,PlotStartOffset,PlotStopOffset,EventChannel,EventTimeChannel,SpikeTime,StepSize,BinWidth);

%Retrive data from the field

%Relative time bin for each event
RelativeTimeMarker=SegmentData.RelativeTimeMarker{1};%=EventTimeMarker;
%Absolute time marker for each event
AbosoluteTimeMarker=SegmentData.AbsoluteTimeMarker{1};

ObjOffTimeMarker=nanmean(RelativeTimeMarker);


%PSTH for each event
MeanPSTH=SegmentData.MeanPSTH;%=[PSTH_Event_Time',PSTH_Event'];
MeanPSTH_TimeOff=MeanPSTH{1};
PSTH_FracOff=MeanPSTH{2};

PSTH_FracOff_mean=nanmean(PSTH_FracOff,1);

PSTH_FracOff_Good_mean=nanmean(PSTH_FracOff(GoodObjectIndex,:),1);
PSTH_FracOff_Good_sem=nanstd(PSTH_FracOff(GoodObjectIndex,:),[],1)./sqrt(sum(~isnan(PSTH_FracOff(GoodObjectIndex,1)),1));

PSTH_FracOff_Bad_mean=nanmean(PSTH_FracOff(~GoodObjectIndex,:),1);
PSTH_FracOff_Bad_sem=nanstd(PSTH_FracOff(~GoodObjectIndex,:),[],1)./sqrt(sum(~isnan(PSTH_FracOff(~GoodObjectIndex,1)),1));

TimeSequence_FracOff_mean=nanmean(MeanPSTH_TimeOff,1);
%}



%%
%Analyze how neural responses changed with the eye location after saccades
%Screen out the valid saccades and saccade start time above 5;
%{
AmplitudeThresholdAfter=3;

TimeThresholdAfter=num2cell(SacEnd_First_Target+50)';

SaccadeAngle_SelectionAfter=cellfun(@(x,y,z,m) z(((x>AmplitudeThresholdAfter)&(y>m))),SaccadeAmplitude,SaccadeStartTime,SaccadeAngle,TimeThresholdAfter,'UniformOutput' ,0);
SaccadeEndPoint_SelectionAfter=cellfun(@(x,y,z,m) z(((x>AmplitudeThresholdAfter)&(y>m)),:),SaccadeAmplitude,SaccadeStartTime,SaccadeEndPoint,TimeThresholdAfter,'UniformOutput' ,0);
SaccadeStartPoint_SelectionAfter=cellfun(@(x,y,z,m) z(((x>AmplitudeThresholdAfter)&(y>m)),:),SaccadeAmplitude,SaccadeStartTime,SaccadeStartPoint,TimeThresholdAfter,'UniformOutput' ,0);

SaccadeStartTime_SelectionAfter=cellfun(@(x,y,z,m) z(((x>AmplitudeThresholdAfter)&(y>m)),:),SaccadeAmplitude,SaccadeStartTime,SaccadeStartTime,TimeThresholdAfter,'UniformOutput' ,0);
SaccadeEndTime_SelectionAfter=cellfun(@(x,y,z,m) z(((x>AmplitudeThresholdAfter)&(y>m)),:),SaccadeAmplitude,SaccadeStartTime,SaccadeEndTime,TimeThresholdAfter,'UniformOutput' ,0);


SaccadeAmplitude_SelectionAfter=cellfun(@(x,y,z,m) z(((x>AmplitudeThresholdAfter)&(y>m)),:),SaccadeAmplitude,SaccadeStartTime,SaccadeAmplitude,TimeThresholdAfter,'UniformOutput' ,0);
%Only look at the first saccade
SaccadeNumInEachTrialAfter=cellfun(@numel,SaccadeAngle_SelectionAfter);
SaccadeNumInEachTrialAfter(isnan(SacEnd_First_Target))=0;

SaccadeAngle_After_First=NaN*ones(1,length(SaccadeNumInEachTrialAfter));

SaccadeEndPoint_After_First=NaN*ones(length(SaccadeNumInEachTrialAfter),2);

RT_First_Target_After=NaN*ones(1,length(SaccadeNumInEachTrialAfter));
SaccadeStartPoint_After_First=NaN*ones(length(SaccadeNumInEachTrialAfter),2);

SaccadeAmplitude_After_First=NaN*ones(1,length(SaccadeNumInEachTrialAfter));

SaccadeAngle_After_First(SaccadeNumInEachTrialAfter>0)=cell2mat(cellfun(@(x) x(1,:),SaccadeAngle_SelectionAfter(SaccadeNumInEachTrialAfter>0),'UniformOutput' ,0));
SaccadeEndPoint_After_First(SaccadeNumInEachTrialAfter>0,:)=cell2mat(cellfun(@(x) x(1,:),SaccadeEndPoint_SelectionAfter(SaccadeNumInEachTrialAfter>0,:),'UniformOutput' ,0));
RT_First_Target_After(SaccadeNumInEachTrialAfter>0)=cell2mat(cellfun(@(x) x(1,:),SaccadeStartTime_SelectionAfter(SaccadeNumInEachTrialAfter>0,:),'UniformOutput' ,0));
SacEnd_First_Target(SaccadeNumInEachTrialAfter>0)=cell2mat(cellfun(@(x) x(1,:),SaccadeEndTime_SelectionAfter(SaccadeNumInEachTrialAfter>0,:),'UniformOutput' ,0));


SaccadeStartPoint_After_First(SaccadeNumInEachTrialAfter>0,:)=cell2mat(cellfun(@(x) x(1,:),SaccadeStartPoint_SelectionAfter(SaccadeNumInEachTrialAfter>0,:),'UniformOutput' ,0));

SaccadeAmplitude_After_First(SaccadeNumInEachTrialAfter>0)=cell2mat(cellfun(@(x) x(1,:),SaccadeAmplitude_SelectionAfter(SaccadeNumInEachTrialAfter>0,:),'UniformOutput' ,0));
%}
%{
%Align data on the return saccade onset
%Align again on return saccade onset

AlignedSpike=SpikeAlignTool(RT_First_Target_After',SpikeTimeAligned_TgtOn,StepSize,BinWidth);

PSTH_SacOnReturn=AlignedSpike.PSTH_Aligned;
TimeSequence_SacOnReturn=AlignedSpike.PSTH_Time_Aligned;


%Select interval
SelectPSTHinterval_SacOnReturn=[-PlotSacPreBuffer,PlotSacPostBuffer];
[PSTH_SacOnReturn,TimeSequence_SacOnReturn]=SelectPSTHInterval(PSTH_SacOnReturn,TimeSequence_SacOnReturn,SelectPSTHinterval_SacOnReturn);


TimeSequence_SacOnReturn_mean=nanmean(TimeSequence_SacOnReturn,1);
PSTH_SacOnReturn_mean=nanmean(PSTH_SacOnReturn,1);
PSTH_SacOnReturn_sem=nanstd(PSTH_SacOnReturn,[],1)./sqrt(sum(~isnan(PSTH_SacOnReturn),1));

%Good versus Bad
PSTH_SacOnReturn_Good_mean=nanmean(PSTH_SacOnReturn(GoodObjectIndex,:),1);
PSTH_SacOnReturn_Good_sem=nanstd(PSTH_SacOnReturn(GoodObjectIndex,:),[],1)./sqrt(sum(~isnan(PSTH_SacOnReturn(GoodObjectIndex,1)),1));

PSTH_SacOnReturn_Bad_mean=nanmean(PSTH_SacOnReturn(~GoodObjectIndex,:),1);
PSTH_SacOnReturn_Bad_sem=nanstd(PSTH_SacOnReturn(~GoodObjectIndex,:),[],1)./sqrt(sum(~isnan(PSTH_SacOnReturn(~GoodObjectIndex,1)),1));


%Align Eye trace with Return saccade onset

%Select the interval around saccade onset

EyeEcc_AlignedAfter=AlignEyeData(Eye_Eccentricity_Tgt,EyeTime,RT_First_Target_After,[-PlotSacPreBuffer,PlotSacPostBuffer]);



EyeEcc_AlignedAfter_mean=nanmean(EyeEcc_AlignedAfter,1);
EyeEcc_AlignedAfter_sem=nanstd(EyeEcc_AlignedAfter,[],1)./sqrt(sum(~isnan(EyeEcc_AlignedAfter),1));

%Good versus Bad

EyeEcc_AlignedAfter_Good_mean=nanmean(EyeEcc_AlignedAfter(GoodObjectIndex,:),1);
EyeEcc_AlignedAfter_Good_sem=nanstd(EyeEcc_AlignedAfter(GoodObjectIndex,:),[],1)./sqrt(sum(~isnan(EyeEcc_AlignedAfter(GoodObjectIndex,:)),1));

EyeEcc_AlignedAfter_Bad_mean=nanmean(EyeEcc_AlignedAfter(~GoodObjectIndex,:),1);
EyeEcc_AlignedAfter_Bad_sem=nanstd(EyeEcc_AlignedAfter(~GoodObjectIndex,:),[],1)./sqrt(sum(~isnan(EyeEcc_AlignedAfter(~GoodObjectIndex,:)),1));

%}

%% Loop over the target directions
UniqueTarget = unique(TargetName);
for i = 1:length(UniqueAngle )
    sel_a = TargetAngle == UniqueAngle(i);% & TargetEcc == UniqueEcc(j);
   % for j = 1:length(UniqueEcc)
        %{
        sel = TargetAngle == UniqueAngle(i) & TargetEcc == UniqueEcc(j);
        
        %Align on target on 
        curr_GoodFracOn_Loc = MeanPSTH_FR_FractalOn(sel & GoodObjectIndex,:);
        Good_FracOn_Loc_Mean{i,j} = nanmean(curr_GoodFracOn_Loc,1);
        Good_FracOn_Loc_Sem{i,j} = nanstd(curr_GoodFracOn_Loc,[],1)./sqrt(sum(~isnan(curr_GoodFracOn_Loc),1));
        

        curr_BadFracOn_Loc = MeanPSTH_FR_FractalOn(sel & ~GoodObjectIndex,:);
        Bad_FracOn_Loc_Mean{i,j} = nanmean(curr_BadFracOn_Loc,1);
        Bad_FracOn_Loc_Sem{i,j} = nanstd(curr_BadFracOn_Loc,[],1)./sqrt(sum(~isnan(curr_BadFracOn_Loc),1));
        
        
        
        AngleEccorder{i,j} = [UniqueAngle(i),UniqueEcc(j)];
        
         %Align on sac on 
        curr_GoodSacOn_Loc = PSTH_SacOn(sel & GoodObjectIndex,:);
        Good_SacOn_Loc_Mean{i,j} = nanmean(curr_GoodSacOn_Loc,1);
        Good_SacOn_Loc_Sem{i,j} = nanstd(curr_GoodSacOn_Loc,[],1)./sqrt(sum(~isnan(curr_GoodSacOn_Loc),1));
        
        curr_BadSacOn_Loc = PSTH_SacOn(sel & ~GoodObjectIndex,:);
        Bad_SacOn_Loc_Mean{i,j} = nanmean(curr_BadSacOn_Loc,1);
        Bad_SacOn_Loc_Sem{i,j} = nanstd(curr_BadSacOn_Loc,[],1)./sqrt(sum(~isnan(curr_BadSacOn_Loc),1));
        %}
        %For each object, compare good vs bad
        for m = 1:length(UniqueTarget)
            sel = sel_a & TargetName == UniqueTarget(m);

            
           

            %Align on target on 
        curr_GoodFracOn_Loc = MeanPSTH_FR_FractalOn(sel & GoodObjectIndex,:);
        SameGood_FracOn_Loc_Mean{i,m} = nanmean(curr_GoodFracOn_Loc,1);
        SameGood_FracOn_Loc_Sem{i,m} = nanstd(curr_GoodFracOn_Loc,[],1)./sqrt(sum(~isnan(curr_GoodFracOn_Loc),1));
        

        curr_BadFracOn_Loc = MeanPSTH_FR_FractalOn(sel & ~GoodObjectIndex,:);
        SameBad_FracOn_Loc_Mean{i,m} = nanmean(curr_BadFracOn_Loc,1);
        SameBad_FracOn_Loc_Sem{i,m} = nanstd(curr_BadFracOn_Loc,[],1)./sqrt(sum(~isnan(curr_BadFracOn_Loc),1));
        
        
        
       % AngleEccorder{i,j} = [UniqueAngle(i),UniqueEcc(j)];
        
         %Align on sac on 
        curr_GoodSacOn_Loc = PSTH_SacOn(sel & GoodObjectIndex,:);
        SameGood_SacOn_Loc_Mean{i,m} = nanmean(curr_GoodSacOn_Loc,1);
        SameGood_SacOn_Loc_Sem{i,m} = nanstd(curr_GoodSacOn_Loc,[],1)./sqrt(sum(~isnan(curr_GoodSacOn_Loc),1));
        
        curr_BadSacOn_Loc = PSTH_SacOn(sel & ~GoodObjectIndex,:);
        SameBad_SacOn_Loc_Mean{i,m} = nanmean(curr_BadSacOn_Loc,1);
        SameBad_SacOn_Loc_Sem{i,m} = nanstd(curr_BadSacOn_Loc,[],1)./sqrt(sum(~isnan(curr_BadSacOn_Loc),1));

        AngleEccorder{i} = [UniqueAngle(i),UniqueEcc(1)];

        end

        
  %  end
end







%%
if ShowFigureFlag
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Figure 1: Average firing rate
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%for each location


figtitlestr{1}='PSTHAlignOnEachFracOnset_EachLocation';

fig{1}=PrepareFigure(110,'w',[50,100, 1200,800],'Name',figtitlestr{1});

for i = 1:length(UniqueAngle)
    for m = 1:length(UniqueTarget)

        %Align on target onset

        subplot(length(UniqueTarget),length(UniqueAngle)*2,(m-1)*length(UniqueAngle)*2 + (i-1)*2+1);

        shadedErrorBar(PSTH_Time_FraOn_All,SameGood_FracOn_Loc_Mean{i,m},SameGood_FracOn_Loc_Sem{i,m},'lineprops',{[1    0    0]},'transparent',1,'patchSaturation',0.3);

        hold on
        shadedErrorBar(PSTH_Time_FraOn_All,SameBad_FracOn_Loc_Mean{i,m},SameBad_FracOn_Loc_Sem{i,m},'lineprops',{[0    0    1]},'transparent',1,'patchSaturation',0.3);

       % yrange=ylim;
       % yrange(2)=max(MeanPSTH_FR_FractalOn_mean(PSTH_Time_FraOn_All<PlotSacPostBuffer & PSTH_Time_FraOn_All>-PlotSacPreBuffer))*1.05;
          xlim([-PlotSacPreBuffer, PlotSacPostBuffer]);%Keep in the same scale
       set(findall(gca, 'Type', 'Line'),'LineWidth',3);

        set(gca,'FontSize',15,'FontWeight','Bold');
        set(gca,'LineWidth',3);
        title(sprintf("Angle: %1.1d,Ecc: %1.1d",AngleEccorder{i}));
        
       xlabel(sprintf('Fractal%d Onset(ms)',UniqueTarget(m)),'FontSize',15,'FontWeight','Bold');
       if i == 1
            ylabel('Firing Rate(Hz)','FontSize',15,'FontWeight','Bold');
       end


       %Align on sac onset

        subplot(length(UniqueTarget),length(UniqueAngle)*2,(m-1)*length(UniqueAngle)*2 + (i-1)*2+2);

       shadedErrorBar(TimeSequence_SacOn_mean,SameGood_SacOn_Loc_Mean{i,m},SameGood_SacOn_Loc_Sem{i,m},'lineprops',{[1    0    0]},'transparent',1,'patchSaturation',0.3);

        hold on
        shadedErrorBar(TimeSequence_SacOn_mean,SameBad_SacOn_Loc_Mean{i,m},SameBad_SacOn_Loc_Sem{i,m},'lineprops',{[0    0    1]},'transparent',1,'patchSaturation',0.3);

       % yrange=ylim;
       % yrange(2)=max(MeanPSTH_FR_FractalOn_mean(PSTH_Time_FraOn_All<PlotSacPostBuffer & PSTH_Time_FraOn_All>-PlotSacPreBuffer))*1.05;
          xlim([-PlotSacPreBuffer, PlotSacPostBuffer]);%Keep in the same scale
       set(findall(gca, 'Type', 'Line'),'LineWidth',3);

        set(gca,'FontSize',15,'FontWeight','Bold');
        set(gca,'LineWidth',3);
        title(sprintf("Angle: %1.1d,Ecc: %1.1d",AngleEccorder{i}));
        xlabel('Sac Onset(ms)','FontSize',15,'FontWeight','Bold');
       if i == 1
            ylabel('Firing Rate(Hz)','FontSize',15,'FontWeight','Bold');
       end


        


    end
end




end %End of ShowFigureFlag



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

%Store figures in common folders
if TaskCode==OBJ_DIR_SAC_TASK  
ProtocolMarker='SatbleValue';
elseif  TaskCode==FLEXIBLE_VALUE_TASK 
    ProtocolMarker='FlexibleValue';
end
FigurePath=strcat(BasicDirectory,'Figures',BasicDirectory(end),ProtocolMarker,BasicDirectory(end));%,OutputFolerName);
if ~exist(FigurePath)
    mkdir(FigurePath)
end

if ShowFigureFlag & OutputFlag
    for i=1:length(fig)
        cd(FigurePath);
        OutputFigureNameFig=[];
        titlestr=figtitlestr{i};
        FolderName=titlestr;
         if ~exist(FolderName,'dir')
          mkdir(FolderName);
         end
          cd(FolderName);
    
        OutputFigureNameFig=strcat(OutputFigureName,'.jpg');
        saveas(fig{i},OutputFigureNameFig);
        
        
    end
        
    
    disp('Figures have been exported to the common folder');
    
end



%{


ExistFlag=0;
%Load the old file if exist
if exist(OutputFileName)

    load(OutputFileName);
    %Load OutputData into the memory
    if isfield(OutputFileName,'TargetSaccade')
    ExistFlag=1;
    end
end


%%%%%%%%%%%%%%%%%%%%OutputData For current analysis
OutputData_New=[];

OutputData_New.TargetSaccade.TrialType='GoodTrials';
%%%%%%%%%Save overall scene response characteristics%%%%%%%%%%%%%%%%%%%%%%%
%Organize the data
NameStr={'PSTH_SacOn_Mean_Time','PSTH_SacOn_FR','PSTH_SacOn_FR_Good','PSTH_SacOn_FR_Bad',...
    'PSTH_Eye_SacOn','PSTH_Eye_SacOn_Good','PSTH_Eye_SacOn_Bad',...
    'PSTH_SacOff_Mean_Time','PSTH_SacOff_FR','PSTH_SacOff_FR_Good','PSTH_SacOn_FR_Bad',...
    'PSTH_Eye_SacOff','PSTH_Eye_SacOff_Good','PSTH_Eye_SacOff_Bad',...
    'PSTH_RTGroup_FraOn_mean','PSTH_RTGroup_FraOn_sem','RT_Group','NumIn_RT_Group','PSTH_RTGroup_FraOn_Eye_mean','PSTH_RTGroup_FraOn_Eye_sem',...
    'PSTH_RTGroup_SacOn_mean','PSTH_RTGroup_SacOn_sem','PSTH_RTGroup_SacOn_Eye_mean','PSTH_RTGroup_SacOn_Eye_sem',...
    'PSTH_RTGroup_SacOff_mean','PSTH_RTGroup_SacOff_sem','PSTH_RTGroup_SacOff_Eye_mean','PSTH_RTGroup_SacOff_Eye_sem',...    
    'PSTH_AmpGroup_FraOn_mean','PSTH_AmpGroup_FraOn_sem','Amp_Group','NumIn_Amp_Group','PSTH_AmpGroup_FraOn_Eye_mean','PSTH_AmpGroup_FraOn_Eye_sem',...
    'PSTH_AmpGroup_SacOn_mean','PSTH_AnpGroup_SacOn_sem','PSTH_AmpGroup_SacOn_Eye_mean','PSTH_AmpGroup_SacOn_Eye_sem',...
    'PSTH_AmpGroup_SacOff_mean','PSTH_AnpGroup_SacOff_sem','PSTH_AmpGroup_SacOff_Eye_mean','PSTH_AmpGroup_SacOff_Eye_sem',...
    'PreSac_Amp_FR_Overall','p_Anova_PreSac_Amp','DurSac_Amp_FR_Overall','p_Anova_DurSac_Amp','PostSac_Amp_FR_Overall','p_Anova_PostSac_Amp',...
    'preSac_r_p','preSac_r_p_Good','preSac_r_p_Bad',...
    'durSac_r_p','durSac_r_p_Good','durSac_r_p_Bad',...
    'postSac_r_p','postSac_r_p_Good','postSac_r_p_Bad',...
    'SacVec','SacVectTuning_PreDurPost_mean','SacVectTuning_PreDurPost_sem','SacVectTuning_p','Pref_Vect','Pref_Vect_Length',...
    'PSTH_SacOnReturn_Mean_Time','PSTH_SacOnReturn_FR','PSTH_SacOnReturn_FR_Good','PSTH_SacOnReturn_FR_Bad',...
    'PSTH_Eye_SacOnReturn','PSTH_Eye_SacOnReturn_Good','PSTH_Eye_SacOnReturn_Bad'
    
    
};

DataLib={TimeSequence_SacOn_mean,[PSTH_SacOn_mean',PSTH_SacOn_sem'],[PSTH_SacOn_Good_mean',PSTH_SacOn_Good_sem'],[PSTH_SacOn_Bad_mean',PSTH_SacOn_Bad_sem'],...
[EyeEcc_Aligned_mean',EyeEcc_Aligned_sem'],[EyeEcc_Aligned_Good_mean',EyeEcc_Aligned_Good_sem'],[EyeEcc_Aligned_Bad_mean',EyeEcc_Aligned_Bad_sem'],...
TimeSequence_SacOff_mean,[PSTH_SacOff_mean',PSTH_SacOff_sem'],[PSTH_SacOff_Good_mean',PSTH_SacOff_Good_sem'],[PSTH_SacOff_Bad_mean',PSTH_SacOff_Bad_sem'],...
[EyeEcc_Aligned_SacOff_mean',EyeEcc_Aligned_SacOff_sem'],[EyeEcc_Aligned_SacOff_Good_mean',EyeEcc_Aligned_SacOff_Good_sem'],[EyeEcc_Aligned_SacOff_Bad_mean',EyeEcc_Aligned_SacOff_Bad_sem'],...
PSTH_FraOn_RTGroup_mean,PSTH_FraOn_RTGroup_sem,[RT_Group_mean',RT_Group_sem'],NumInGroupRT',EyeEcc_FraOn_RTGroup_mean,EyeEcc_FraOn_RTGroup_sem,...
PSTH_SacOn_RTGroup_mean,PSTH_SacOn_RTGroup_sem,EyeEcc_SacOn_RTGroup_mean,EyeEcc_SacOn_RTGroup_sem,...
PSTH_SacOff_RTGroup_mean,PSTH_SacOff_RTGroup_sem,EyeEcc_SacOff_RTGroup_mean,EyeEcc_SacOff_RTGroup_sem,...
PSTH_FraOn_AmpGroup_mean,PSTH_FraOn_AmpGroup_sem,[Amp_Group_mean',Amp_Group_sem'],NumInGroupAmp',EyeEcc_FracOn_AmpGroup_mean,EyeEcc_FracOn_AmpGroup_sem,...
PSTH_SacOn_AmpGroup_mean,PSTH_SacOn_AmpGroup_sem,EyeEcc_SacOn_AmpGroup_mean,EyeEcc_SacOn_AmpGroup_sem,...
PSTH_SacOff_AmpGroup_mean,PSTH_SacOff_AmpGroup_sem,EyeEcc_SacOff_AmpGroup_mean,EyeEcc_SacOff_AmpGroup_sem,...
[SpikeRate_PreSac_Amp_mean,SpikeRate_PreSac_Amp_sem],ANOVA_p_Amp_pre,[SpikeRate_DurSac_Amp_mean,SpikeRate_DurSac_Amp_sem],ANOVA_p_Amp_dur,[SpikeRate_PostSac_Amp_mean,SpikeRate_PostSac_Amp_sem],ANOVA_p_Amp_post,...
[r_preSac,p_preSac],[r_preSac_Good,p_preSac_Good],[r_preSac_Bad,p_preSac_Bad],...
[r_durSac,p_durSac],[r_durSac_Good,p_durSac_Good],[r_durSac_Bad,p_durSac_Bad],...
[r_postSac,p_postSac],[r_postSac_Good,p_postSac_Good],[r_postSac_Bad,p_postSac_Bad],...
SaccadeVector,FiringRate_PreDurPost_mean,FiringRate_PreDurPost_sem,SacVecTuning_p,PrefVector,PrefVector_Amp,...
TimeSequence_SacOnReturn_mean,[PSTH_SacOnReturn_mean',PSTH_SacOnReturn_sem'],[PSTH_SacOnReturn_Good_mean',PSTH_SacOnReturn_Good_sem'],[PSTH_SacOnReturn_Bad_mean',PSTH_SacOnReturn_Bad_sem'],...
[EyeEcc_AlignedAfter_mean',EyeEcc_AlignedAfter_sem'],[EyeEcc_AlignedAfter_Good_mean',EyeEcc_AlignedAfter_Good_sem'],[EyeEcc_AlignedAfter_Bad_mean',EyeEcc_AlignedAfter_Bad_sem']
};


DataStamp=containers.Map(NameStr,DataLib);





OutputData_New.TargetSaccade.Task=TaskType;
OutputData_New.TargetSaccade.TaskCode=TaskCode;
OutputData_New.TargetSaccade.DataStamp=DataStamp;


%Store figures
if ShowFigureFlag
   
    %OutputData_New.SceneTuning.Figures={fig99,fig100,fig101,fig102};
    %Output figures in a folder;
    FolderName=OutputFigureName;
    if ~exist(FolderName,'dir')
        mkdir(FolderName);
    end
    cd(FolderName);
    
    TotalFigure=length(fig);
    
    for j=1:TotalFigure
        OutputFigureName=strcat(OutputFigureName,figtitlestr{j},'.jpg');
        if ~isempty(fig{j})
        saveas(fig{j},OutputFigureName);
        end
        
    end

    disp('Figures have been exported to the neuron folder');

end

if ExistFlag
%Compare the old one with the new one
Task_Old=OutputData.TargetSaccade.TaskCode;

     if Task_Old~=TaskCode
        %Not the same task, add the new dataset into the old one
        OutputData.TargetSaccade(numel(OutputData.TargetSaccade)+1)=OutputData_New.TargetSaccade;
     else
         %If the same task,the same protocol,replace the old one with the new one
         
         OutputData.TargetSaccade=OutputData_New.TargetSaccade;
         
    
    
     end
else
    %Output the current file
   OutputData.TargetSaccade=OutputData_New.TargetSaccade; 
end



if OutputFlag
    cd(OutPath);
    save(OutputFileName,'OutputData');
    disp('Data Exported');
else
    disp('No data export');
end  

%}

end %End of the spikechannel

end

function organized=ReorganizeEye(data);
organized=[];
for i=1:length(data)
    data_curr=data{i};
    maxnumel=max(cellfun(@numel,data_curr));
    organized{i}=cell2mat(cellfun(@(x) [x,NaN*ones(1,maxnumel-length(x))],data_curr,'uniform',0));
    
end

 maxnumel=max(cellfun(@(x) size(x,2),organized));
  organized=cell2mat(cellfun(@(x) [x,NaN*ones(size(x,1),maxnumel-size(x,2))],organized,'uniform',0)');

end
function FR=SelectFringRateInterval(FR_Raw,Time,Interval);
for i=1:size(Interval,1)
    if sum(~isnan(Interval(i,:)))==2
        FR(i)=nanmean(FR_Raw(i,Time>=Interval(i,1) & Time<=Interval(i,2)));
    else
         FR(i)=NaN;
    end
    
    
    
end




end
function ObjectIndex=ReproduceFromEvent(event,code)
for i=1:size(event,1)
    
    ObjectIndex(i)=event(i,ismember(event(i,:),code));
    
  
    
end








end