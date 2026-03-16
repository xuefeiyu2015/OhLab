function DifferentPeriodCompare(Data, StartTrial, EndTrial,ShowFigureFlag,OutputFlag);
%Function to analysis neural responses and target property and saccades
%properties
%Xuefei Yu 05122020

ECode;
SegEvent=[SHOWFIXCD,TGTCD,RWDCD];%Event Flow of interest
%SegEventEnd=[SHOWFIXCD,FIXOFF,TGTOFF,RWDOFFCD ,NaN,EXTRAFP];%Event end
SegEventEnd=[NaN,TGTOFF,NaN];%Event end
SegEventStr={'FixOn','TargetOn','Reward'};
AmpMatchFlag=0;%1:Amplitude match on ;0:Amplitude Match off


PlotStartOffset=[200,200,200];%in ms
PlotStopOffset=[800,8000,400];%in ms

ExamplePlot=0;%PLots for presentations, turn off for regular use

ResponseLatencyThreshold=3;

%Parameter to define the sliding window to measure firing rates
StepSize=5;%in ms
BinWidth=15;%in ms


%Criterias for sexwlecting saccades
V_Threshold=40;
ContinueBin=5;
ISI_Threshold=20;

%Time bin pre saccade and time bin offset after saccade onset
PlotSacPreBuffer=500;%in ms
PlotSacPostBuffer=500;%in ms

PlotFracOffPreBuffer=500;%in ms
PlotFracOffPostBuffer=2000;%in ms


%Color Definition for each fractals
%{ee
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
    

%Remove choice trials
try
TaskType = Data.ParaLib('TaskType');
%Only maintain the sigle target trials
%SingleTarget =TaskType(TaskType(:,2) == 1,1) ;
SingleTarget =find(TaskType(:,2) == 1) ;
catch
    SingleTarget = 1:length(Selection);
end

Selection = Selection(ismember(Selection,SingleTarget));



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
%TargetPolar=Data.TargetPolarLoc(Selection,:,2:end);
if length(size(Data.TargetPolarLoc)) == 3

    TargetAngle = Data.TargetPolarLoc(Selection,2,1);

    TargetEcc =  Data.TargetPolarLoc(Selection,2,2);
else
    TargetAngle = Data.TargetPolarLoc(Selection,2);

    TargetEcc =  Data.TargetPolarLoc(Selection,3);

end

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
ObjectIndex=cell2mat(ReproduceFromEvent(EventChannel,[GOODOBJ  BADOBJ]));

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

%find mean and std for z-score
FR_Each = (sum(~isnan(SpikeTime),2)-2)./(SpikeTime(:,end)-SpikeTime(:,end-1))*1000;
FR_Mean = mean(FR_Each);
FR_Std = std(FR_Each);


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


%Raster
Raster_SacOn_whole = AlignedSpike.SpikeTime_Aligned_Whole;



%Select interval
SelectPSTHinterval_SacOn=[-PlotSacPreBuffer,PlotSacPostBuffer];
[PSTH_SacOn,TimeSequence_SacOn]=SelectPSTHInterval(PSTH_SacOn,TimeSequence_SacOn,SelectPSTHinterval_SacOn);



TimeSequence_SacOn_mean=nanmean(TimeSequence_SacOn,1);
PSTH_SacOn_mean=nanmean(PSTH_SacOn,1);
PSTH_SacOn_sem=nanstd(PSTH_SacOn,[],1)./sqrt(sum(~isnan(PSTH_SacOn),1));

%Change raster
Raster_SacOn= SelectSpikeTimeInterval(Raster_SacOn_whole,SelectPSTHinterval_SacOn(1)*ones(size(Raster_SacOn_whole,1),1),SelectPSTHinterval_SacOn(end)*ones(size(Raster_SacOn_whole,1),1));

%Good versus Bad
PSTH_SacOn_Good_mean=nanmean(PSTH_SacOn(GoodObjectIndex,:),1);
PSTH_SacOn_Good_sem=nanstd(PSTH_SacOn(GoodObjectIndex,:),[],1)./sqrt(sum(~isnan(PSTH_SacOn(GoodObjectIndex,1)),1));

PSTH_SacOn_Bad_mean=nanmean(PSTH_SacOn(~GoodObjectIndex,:),1);
PSTH_SacOn_Bad_sem=nanstd(PSTH_SacOn(~GoodObjectIndex,:),[],1)./sqrt(sum(~isnan(PSTH_SacOn(~GoodObjectIndex,1)),1));

%Good versus Bad raster
%Raster_SacOn_Good = Raster_SacOn(GoodObjectIndex,:);
%Raster_SacOn_Bad = Raster_SacOn(~GoodObjectIndex,:);

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

%ObjOnTimeMarker=nanmean(RelativeTimeMarker);

RwdTimeAbs=FindOutTime(EventChannel,EventBin,RWDCD);
RwdTimeRel=RwdTimeAbs-AbosoluteTimeMarker(:,1);

%Algin marker to sac-on & sac-off
RwdTimeRelSacEnd=RwdTimeRel-SacEnd_First_Target';
TgtOffTimeRelSacEnd=RelativeTimeMarker(:,2)-SacEnd_First_Target';

RwdTimeRelSacEndMean=nanmean(RwdTimeRelSacEnd);
TgtOffTimeRelSacEndMean=nanmean(TgtOffTimeRelSacEnd);

%PSTH for each event
MeanPSTH=SegmentData.MeanPSTH;%=[PSTH_Event_Time',PSTH_Event'];
MeanPSTH_Time=MeanPSTH(:,1);
MeanPSTH_FR=MeanPSTH(:,2);

%keyboard
%MeanPSTH_FR_FractalOn=MeanPSTH{:,2};

%Raster 
%Raster_TgtOn = SegmentData.RasterAfterAlign{1};


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

%Good versus Bad raster
%Raster_TgtOn_Good = Raster_TgtOn(GoodObjectIndex,:);
%Raster_TgtOn_Bad = Raster_TgtOn(~GoodObjectIndex,:);


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




%%
%Analyze how neural responses changed with the eye location after saccades
%Screen out the valid saccades and saccade start time above 5;
%
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



%% Loop over the target directions
for i = 1:length(UniqueAngle )
    for j = 1:length(UniqueEcc)
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

        %Average 
        %%
       
        TimeTgt = [100,400];
        Ave_FR_TgtOn_Good = nanmean(curr_GoodFracOn_Loc(:,PSTH_Time_FraOn_All >= TimeTgt(1) & PSTH_Time_FraOn_All <= TimeTgt(2)),2);
        Ave_FR_TgtOn_Bad = nanmean(curr_BadFracOn_Loc(:,PSTH_Time_FraOn_All >= TimeTgt(1) & PSTH_Time_FraOn_All <= TimeTgt(2)),2);

        TimeSac = [-200,0];
        N = 100;
        permuteN = 100;

        [TgtAUC(i,j),bestz,perm] = rocN(Ave_FR_TgtOn_Good,Ave_FR_TgtOn_Bad,N,permuteN);
        p_AUC_Tgt(i,j) = perm.pValue;



        Ave_FR_SacOn_Good = nanmean(curr_GoodSacOn_Loc(:,TimeSequence_SacOn_mean >= TimeSac(1) & TimeSequence_SacOn_mean <= TimeSac(2)),2);
        Ave_FR_SacOn_Bad = nanmean(curr_BadSacOn_Loc(:,TimeSequence_SacOn_mean >= TimeSac(1) & TimeSequence_SacOn_mean <= TimeSac(2)),2);
    
        [h,p_Tgt(i)]= ttest2(Ave_FR_TgtOn_Good,Ave_FR_TgtOn_Bad);
        [h,p_Sac(i)]= ttest2(Ave_FR_SacOn_Good,Ave_FR_SacOn_Bad);


        %% Ave
        %Average every 100ms
        index_i = 1;
        for ii = 0:100:300
            time_interval = PSTH_Time_FraOn_All >= ii & PSTH_Time_FraOn_All <= ii+100;
            Ave_FR_TgtOn_Good_Interval = nanmean(curr_GoodFracOn_Loc(:,time_interval),2);
            Ave_FR_TgtOn_Bad_Interval = nanmean(curr_BadFracOn_Loc(:,time_interval),2);

            N = 100;
        permuteN = 100;

        [TgtAUC_interval(i,index_i),bestz,perm] = rocN(Ave_FR_TgtOn_Good_Interval,Ave_FR_TgtOn_Bad_Interval ,N,permuteN);
        p_AUC_Tgt_interval(i,index_i)= perm.pValue;




            [h,p_interval(i,index_i )]=ttest2(Ave_FR_TgtOn_Good_Interval ,Ave_FR_TgtOn_Bad_Interval);
            index_i=index_i +1;

           
        end


        index_i = 1;
        for ii = 0:100:300
            time_interval = PSTH_Time_FraOn_All >= ii & PSTH_Time_FraOn_All <= ii+100;
            Ave_FR_TgtOn_Good_Interval = nanmean(curr_GoodFracOn_Loc(:,time_interval),2);
            Ave_FR_TgtOn_Bad_Interval = nanmean(curr_BadFracOn_Loc(:,time_interval),2);


            time_interval_spon = PSTH_Time_FraOn_All <=0 & PSTH_Time_FraOn_All >= -100;

            Ave_FR_TgtOn_Spon_Good = nanmean(curr_GoodFracOn_Loc(:,time_interval_spon),2);
            Ave_FR_TgtOn_Spon_Bad = nanmean(curr_BadFracOn_Loc(:,time_interval_spon),2);

            Spon = [Ave_FR_TgtOn_Spon_Good;Ave_FR_TgtOn_Spon_Bad];
            FR = [Ave_FR_TgtOn_Good_Interval;Ave_FR_TgtOn_Bad_Interval];

            [h,p_response_interval(i,index_i)] = ttest(FR,Spon);

            %{
            N = 100;
        permuteN = 100;

        [TgtAUC_interval(i,index_i),bestz,perm] = rocN(Ave_FR_TgtOn_Good_Interval,Ave_FR_TgtOn_Bad_Interval ,N,permuteN);
        p_AUC_Tgt_interval(i,index_i)= perm.pValue;




            [h,p_interval(i,index_i )]=ttest2(Ave_FR_TgtOn_Good_Interval ,Ave_FR_TgtOn_Bad_Interval);
            %}
            index_i=index_i +1;

           
        end
       

        TgtGood(i) = mean(Ave_FR_TgtOn_Good);
        TgtBad(i) = mean(Ave_FR_TgtOn_Bad);

        SacGood(i) = mean(Ave_FR_SacOn_Good);
        SacBad(i) = mean(Ave_FR_SacOn_Bad);

        %Raster
        Raster_TgtOn_Good{i}=Raster_TgtOn(sel & GoodObjectIndex,:);
        Raster_TgtOn_Bad{i}=Raster_TgtOn(sel & ~GoodObjectIndex,:);

        Raster_SacOn_Good{i}=Raster_SacOn(sel & GoodObjectIndex,:);
        Raster_SacOn_Bad{i}=Raster_SacOn(sel & ~GoodObjectIndex,:);




        
    end
end

%Use contralateral RF for export

RF = UniqueAngle > 0 & UniqueAngle <= 180;
if length(RF)>1
    RF_New = 0*ones(1,length(RF));
    first = find(RF==1,1,'first');
    RF_New(first) =1;
    RF=RF_New';

end


%Transfer to real RF
%{
if length(UniqueAngle) >=2
if (TgtGood(1)+TgtBad(1)) > (TgtGood(2)+TgtBad(2)) 
    RF = UniqueAngle==UniqueAngle(1);
else
    RF = UniqueAngle==UniqueAngle(2);
end
else
    RF = UniqueAngle==UniqueAngle(1);
end
%}
p_response_interval_RF=p_response_interval(find(RF==1),:);
p_response_interval_nonRF=p_response_interval(find(RF==0,1,'first'),:);

TgtAUC_interval_RF=TgtAUC_interval(find(RF==1),:);
TgtAUC_interval_nonRF=TgtAUC_interval(find(RF==0),:);

p_interval_RF = p_interval(find(RF==1));
p_interval_nonRF = p_interval(find(RF==0));



%Time for export
ExportTime = [-400,400];
SelectTime = PSTH_Time_FraOn_All>ExportTime(1) & PSTH_Time_FraOn_All<ExportTime(2);
Time_TgtOn = PSTH_Time_FraOn_All(SelectTime);

ExportTimeSac = [-200,400];
SelectTimeSac = TimeSequence_SacOn_mean>ExportTimeSac(1) & TimeSequence_SacOn_mean<ExportTimeSac(2);
Time_SacOn = TimeSequence_SacOn_mean(SelectTimeSac);




PSTH_TgtOn_Good_RF = Good_FracOn_Loc_Mean{find(RF==1),1}(SelectTime);
PSTH_TgtOn_Bad_RF = Bad_FracOn_Loc_Mean{find(RF==1),1}(SelectTime);

PSTH_TgtOn_Good_RF_z = (PSTH_TgtOn_Good_RF - FR_Mean)/FR_Std;
PSTH_TgtOn_Bad_RF_z = (PSTH_TgtOn_Bad_RF - FR_Mean)/FR_Std;




PSTH_SacOn_Good_RF = Good_SacOn_Loc_Mean{find(RF==1),1}(SelectTimeSac);
PSTH_SacOn_Bad_RF = Bad_SacOn_Loc_Mean{find(RF==1),1}(SelectTimeSac);

PSTH_SacOn_Good_RF_z = (PSTH_SacOn_Good_RF - FR_Mean)/FR_Std;
PSTH_SacOn_Bad_RF_z = (PSTH_SacOn_Bad_RF - FR_Mean)/FR_Std;

if length(UniqueAngle) == 1
    PSTH_TgtOn_Good_nonRF = NaN*Good_FracOn_Loc_Mean{1}(SelectTime);
    PSTH_TgtOn_Bad_nonRF = NaN*Bad_FracOn_Loc_Mean{1}(SelectTime);

    PSTH_SacOn_Good_nonRF = NaN * Good_SacOn_Loc_Mean{1}(SelectTimeSac);
    PSTH_SacOn_Bad_nonRF = NaN * Bad_SacOn_Loc_Mean{1}(SelectTimeSac);

    TgtGood_nonRF = NaN;
    TgtBad_nonRF = NaN;
    SacGood_nonRF = NaN;
    SacBad_nonRF = NaN;

    



else
    PSTH_TgtOn_Good_nonRF = Good_FracOn_Loc_Mean{find(RF==0,1,'first'),1}(SelectTime);
    PSTH_TgtOn_Bad_nonRF = Bad_FracOn_Loc_Mean{find(RF==0,1,'first'),1}(SelectTime);

    PSTH_SacOn_Good_nonRF = Good_SacOn_Loc_Mean{find(RF==0,1,'first'),1}(SelectTimeSac);
    PSTH_SacOn_Bad_nonRF = Bad_SacOn_Loc_Mean{find(RF==0,1,'first'),1}(SelectTimeSac);

    TgtGood_nonRF = TgtGood(find(RF==0,1,'first'));
    TgtBad_nonRF = TgtBad(find(RF==0,1,'first'));
    SacGood_nonRF = SacGood(find(RF==0,1,'first'));
    SacBad_nonRF = SacBad(find(RF==0,1,'first'));
end



PSTH_TgtOn_Good_nonRF_z = (PSTH_TgtOn_Good_nonRF - FR_Mean)/FR_Std;
PSTH_TgtOn_Bad_nonRF_z = (PSTH_TgtOn_Bad_nonRF- FR_Mean)/FR_Std;

PSTH_SacOn_Good_nonRF_z = (PSTH_SacOn_Good_nonRF- FR_Mean)/FR_Std;
PSTH_SacOn_Bad_nonRF_z = (PSTH_SacOn_Bad_nonRF- FR_Mean)/FR_Std;



TgtGood_RF = TgtGood(find(RF==1));
TgtBad_RF = TgtBad(find(RF==1));



SacGood_RF = SacGood(find(RF==1));
SacBad_RF = SacBad(find(RF==1));



p_Tgt_RF  = p_Tgt(find(RF==1));
p_Sac_RF  = p_Sac(find(RF==1));

p_Tgt_nonRF  = p_Tgt(find(RF==0,1,'first'));
p_Sac_nonRF  = p_Sac(find(RF==0,1,'first'));


TgtAUC_RF = TgtAUC(find(RF==1),1);
TgtAUC_nonRF = TgtAUC(find(RF==0,1,'first'),1);

p_AUC_Tgt_RF = p_AUC_Tgt(find(RF==1),1);
p_AUC_Tgt_nonRF = p_AUC_Tgt(find(RF==0,1,'first'),1);





if p_Tgt_RF < 0.05
    if TgtGood_RF > TgtBad_RF
        PSTH_TgtOn_Prefer_RF = PSTH_TgtOn_Good_RF;
        PSTH_TgtOn_nonPrefer_RF = PSTH_TgtOn_Bad_RF;

        PSTH_SacOn_Prefer_RF = PSTH_SacOn_Good_RF;
        PSTH_SacOn_nonPrefer_RF = PSTH_SacOn_Bad_RF;

        PSTH_TgtOn_Prefer_nonRF = PSTH_TgtOn_Good_nonRF;
        PSTH_TgtOn_nonPrefer_nonRF = PSTH_TgtOn_Bad_nonRF;

        PSTH_SacOn_Prefer_nonRF = PSTH_SacOn_Good_nonRF;
        PSTH_SacOn_nonPrefer_nonRF = PSTH_SacOn_Bad_nonRF;

    else
        PSTH_TgtOn_Prefer_RF = PSTH_TgtOn_Bad_RF;
        PSTH_TgtOn_nonPrefer_RF = PSTH_TgtOn_Good_RF;

        PSTH_SacOn_Prefer_RF = PSTH_SacOn_Bad_RF;
        PSTH_SacOn_nonPrefer_RF = PSTH_SacOn_Good_RF;

        PSTH_TgtOn_Prefer_nonRF = PSTH_TgtOn_Bad_nonRF;
        PSTH_TgtOn_nonPrefer_nonRF = PSTH_TgtOn_Good_nonRF;

        PSTH_SacOn_Prefer_nonRF = PSTH_SacOn_Bad_nonRF;
        PSTH_SacOn_nonPrefer_nonRF = PSTH_SacOn_Good_nonRF;
    end


elseif p_Sac_RF < 0.05
    if SacGood_RF > SacBad_RF
        PSTH_TgtOn_Prefer_RF = PSTH_TgtOn_Good_RF;
        PSTH_TgtOn_nonPrefer_RF = PSTH_TgtOn_Bad_RF;

        PSTH_SacOn_Prefer_RF = PSTH_SacOn_Good_RF;
        PSTH_SacOn_nonPrefer_RF = PSTH_SacOn_Bad_RF;

        PSTH_TgtOn_Prefer_nonRF = PSTH_TgtOn_Good_nonRF;
        PSTH_TgtOn_nonPrefer_nonRF = PSTH_TgtOn_Bad_nonRF;

        PSTH_SacOn_Prefer_nonRF = PSTH_SacOn_Good_nonRF;
        PSTH_SacOn_nonPrefer_nonRF = PSTH_SacOn_Bad_nonRF;


    else
        PSTH_TgtOn_Prefer_RF = PSTH_TgtOn_Bad_RF;
        PSTH_TgtOn_nonPrefer_RF = PSTH_TgtOn_Good_RF;

        PSTH_SacOn_Prefer_RF = PSTH_SacOn_Bad_RF;
        PSTH_SacOn_nonPrefer_RF = PSTH_SacOn_Good_RF;

        PSTH_TgtOn_Prefer_nonRF = PSTH_TgtOn_Bad_nonRF;
        PSTH_TgtOn_nonPrefer_nonRF = PSTH_TgtOn_Good_nonRF;

        PSTH_SacOn_Prefer_nonRF = PSTH_SacOn_Bad_nonRF;
        PSTH_SacOn_nonPrefer_nonRF = PSTH_SacOn_Good_nonRF;
    end
    
else
    if TgtGood_RF > TgtBad_RF
        PSTH_TgtOn_Prefer_RF = PSTH_TgtOn_Good_RF;
        PSTH_TgtOn_nonPrefer_RF = PSTH_TgtOn_Bad_RF;

        PSTH_SacOn_Prefer_RF = PSTH_SacOn_Good_RF;
        PSTH_SacOn_nonPrefer_RF = PSTH_SacOn_Bad_RF;

        PSTH_TgtOn_Prefer_nonRF = PSTH_TgtOn_Good_nonRF;
        PSTH_TgtOn_nonPrefer_nonRF = PSTH_TgtOn_Bad_nonRF;

        PSTH_SacOn_Prefer_nonRF = PSTH_SacOn_Good_nonRF;
        PSTH_SacOn_nonPrefer_nonRF = PSTH_SacOn_Bad_nonRF;

    else
        PSTH_TgtOn_Prefer_RF = PSTH_TgtOn_Bad_RF;
        PSTH_TgtOn_nonPrefer_RF = PSTH_TgtOn_Good_RF;

        PSTH_SacOn_Prefer_RF = PSTH_SacOn_Bad_RF;
        PSTH_SacOn_nonPrefer_RF = PSTH_SacOn_Good_RF;

        PSTH_TgtOn_Prefer_nonRF = PSTH_TgtOn_Bad_nonRF;
        PSTH_TgtOn_nonPrefer_nonRF = PSTH_TgtOn_Good_nonRF;

        PSTH_SacOn_Prefer_nonRF = PSTH_SacOn_Bad_nonRF;
        PSTH_SacOn_nonPrefer_nonRF = PSTH_SacOn_Good_nonRF;
    end

end




PSTH_TgtOn_Prefer_RF_z = (PSTH_TgtOn_Prefer_RF - FR_Mean)/FR_Std;

PSTH_TgtOn_nonPrefer_RF_z = (PSTH_TgtOn_nonPrefer_RF - FR_Mean)/FR_Std;

PSTH_SacOn_Prefer_RF_z = (PSTH_SacOn_Prefer_RF - FR_Mean)/FR_Std;
PSTH_SacOn_nonPrefer_RF_z = (PSTH_SacOn_nonPrefer_RF- FR_Mean)/FR_Std;

PSTH_TgtOn_Prefer_nonRF_z = (PSTH_TgtOn_Prefer_nonRF - FR_Mean)/FR_Std;
PSTH_TgtOn_nonPrefer_nonRF_z = (PSTH_TgtOn_nonPrefer_nonRF - FR_Mean)/FR_Std;

PSTH_SacOn_Prefer_nonRF_z = (PSTH_SacOn_Prefer_nonRF - FR_Mean)/FR_Std;
PSTH_SacOn_nonPrefer_nonRF_z = (PSTH_SacOn_nonPrefer_nonRF- FR_Mean)/FR_Std;









%%
if ShowFigureFlag
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Figure 1: Average firing rate
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%for each location
if length(UniqueAngle)>2

figtitlestr{1}='PSTHAlignOnFracOnset_EachLocation';

fig{1}=PrepareFigure(110,'w',[50,100, 1200,800],'Name',figtitlestr{1});

%
%Set a layout for each angle 
if length(UniqueAngle) <= 2
  angle = [90,270];
  subloc = [1,3];  
  layout = [1,3];
  center = 1;
    
else
angle = 0:45:359;
subloc = [2,1,4,7,8,9,6,3];
if length(UniqueEcc) == 1 || (length(UniqueEcc) == 2 && ismember(0,UniqueEcc))
    layout = [3,3];
    center = 5;
else
    totalpart = (length(UniqueEcc) - ismember(0,UniqueEcc))*length(angle)+ismember(0,UniqueEcc);
    layout = [ceil(sqrt(totalpart)),ceil(sqrt(totalpart))];
    center = (layout(1)+1)*2+1;
    
end
end
%
for i = 1:length(UniqueAngle )
    for j = 1:length(UniqueEcc)
        if UniqueEcc(j) == 0
            index = center;
        else
            [minValue,closestIndex] = min(UniqueAngle(i)-angle);
            index = subloc(closestIndex) ;
        end
        %subplot(layout(1),layout(2),index);
        subplot(max(ceil(length(UniqueAngle)/2),1),2,i);
        
        shadedErrorBar(PSTH_Time_FraOn_All,Good_FracOn_Loc_Mean{i,j},Good_FracOn_Loc_Sem{i,j},'lineprops',{[1    0    0]},'transparent',1,'patchSaturation',0.3);

        hold on
        shadedErrorBar(PSTH_Time_FraOn_All,Bad_FracOn_Loc_Mean{i,j},Bad_FracOn_Loc_Sem{i,j},'lineprops',{[0    0    1]},'transparent',1,'patchSaturation',0.3);

       % yrange=ylim;
       % yrange(2)=max(MeanPSTH_FR_FractalOn_mean(PSTH_Time_FraOn_All<PlotSacPostBuffer & PSTH_Time_FraOn_All>-PlotSacPreBuffer))*1.05;
          xlim([-PlotSacPreBuffer, PlotSacPostBuffer]);%Keep in the same scale
       set(findall(gca, 'Type', 'Line'),'LineWidth',3);

        set(gca,'FontSize',15,'FontWeight','Bold');
        set(gca,'LineWidth',3);
        title(sprintf("Angle: %1.1d,Ecc: %1.1d",AngleEccorder{i,j}));
        if index == center + layout(1)
            xlabel('Time from Fractal Onset(ms)','FontSize',15,'FontWeight','Bold');
            
        end
        
        ylabel('Firing Rate(Hz)','FontSize',15,'FontWeight','Bold');
        

       
        
        
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Figure 1: Average firing rate
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%for each location


figtitlestr{2}='PSTHAlignOnSacOnset_EachLocation';

fig{2}=PrepareFigure(111,'w',[70,100, 1200,800],'Name',figtitlestr{2});



for i = 1:length(UniqueAngle )
    for j = 1:length(UniqueEcc)
        if UniqueEcc(j) == 0
            index = center;
        else
            [minValue,closestIndex] = min(abs(UniqueAngle(i)-angle));
            index = subloc(closestIndex) ;
        end
        subplot(max(ceil(length(UniqueAngle))/2,1),2,i);
        
        shadedErrorBar(TimeSequence_SacOn_mean,Good_SacOn_Loc_Mean{i,j},Good_SacOn_Loc_Sem{i,j},'lineprops',{[1    0    0]},'transparent',1,'patchSaturation',0.3);

        hold on
        shadedErrorBar(TimeSequence_SacOn_mean,Bad_SacOn_Loc_Mean{i,j},Bad_SacOn_Loc_Sem{i,j},'lineprops',{[0    0    1]},'transparent',1,'patchSaturation',0.3);

       % yrange=ylim;
       % yrange(2)=max(MeanPSTH_FR_FractalOn_mean(PSTH_Time_FraOn_All<PlotSacPostBuffer & PSTH_Time_FraOn_All>-PlotSacPreBuffer))*1.05;
          xlim([-PlotSacPreBuffer, PlotSacPostBuffer]);%Keep in the same scale
       set(findall(gca, 'Type', 'Line'),'LineWidth',3);

        set(gca,'FontSize',15,'FontWeight','Bold');
        set(gca,'LineWidth',3);
        title(sprintf("Angle: %1.1d,Ecc: %1.1d",AngleEccorder{i,j}));
        if index == center + layout(1)
            xlabel('Time from Saccade Onset(ms)','FontSize',15,'FontWeight','Bold');
            
        end
        
        ylabel('Firing Rate(Hz)','FontSize',15,'FontWeight','Bold');
        

       
        
        
    end
end

%%%%%%%%%%%%%%%%%%%%Figure 3 Raster
figtitlestr{3}='Raster_EachLocation';

fig{3}=PrepareFigure(112,'w',[50,100, 1200,800],'Name',figtitlestr{3});

StartLoc = 1;
LineWidth = 2;
Barlength = 1;


for i = 1:length(UniqueAngle)
    subplot(2,length(UniqueAngle),(i-1)*2+1);
    SpikeRaterPlot(Raster_TgtOn_Good{i},StartLoc,Barlength,'r',LineWidth);
    
 hold on
 IncreaseBar = size(Raster_TgtOn_Good{i},1);
 
    SpikeRaterPlot(Raster_TgtOn_Bad{i},StartLoc+IncreaseBar,Barlength,'b',LineWidth);
    yrange = ylim;
 hold on
 plot([0,0],yrange,'--k','LineWidth',3);
 xlim([-PlotSacPreBuffer, PlotSacPostBuffer]);
 xlabel('Time from Fractal Onset(ms)','FontSize',15,'FontWeight','Bold');
  title(sprintf("Angle: %1.1d",AngleEccorder{i,1}));


     subplot(2,length(UniqueAngle),(i-1)*2+2);
    SpikeRaterPlot(Raster_SacOn_Good{i},StartLoc,Barlength,'r',LineWidth);
    
 hold on
 IncreaseBar = size(Raster_SacOn_Good{i},1);
 
    SpikeRaterPlot(Raster_SacOn_Bad{i},StartLoc+IncreaseBar,Barlength,'b',LineWidth);
    yrange = ylim;
 hold on
 plot([0,0],yrange,'--k','LineWidth',3);
xlim([-PlotSacPreBuffer, PlotSacPostBuffer])

xlabel('Time from Saccade Onset(ms)','FontSize',15,'FontWeight','Bold');
 title(sprintf("Angle: %1.1d",AngleEccorder{i,1}));


end

else
%Put all the figures into one for better visualization
figtitlestr{1}='PSTH_GoodBad_EachLocation';

fig{1}=PrepareFigure(110,'w',[50,100, 1200,800],'Name',figtitlestr{1});

StartLoc = 1;
LineWidth = 2;
Barlength = 1;

    for i = 1:length(UniqueAngle)
        subplot(4,length(UniqueAngle),i);
        %Raster

        SpikeRaterPlot(Raster_TgtOn_Good{i},StartLoc,Barlength,'r',LineWidth);
        hold on
        IncreaseBar = size(Raster_TgtOn_Good{i},1);
 
        SpikeRaterPlot(Raster_TgtOn_Bad{i},StartLoc+IncreaseBar,Barlength,'b',LineWidth);
        yrange = ylim;
        hold on
        plot([0,0],yrange,'--k','LineWidth',3);
        set(findall(gca, 'Type', 'Line'),'LineWidth',3);
     %   xlim([-PlotSacPreBuffer, PlotSacPostBuffer]);
     xlim([-200, 400]);%Keep in the same scale
        title(sprintf("Angle: %1.1d,Ecc: %1.1d",AngleEccorder{i,1}));

        subplot(4,length(UniqueAngle),length(UniqueAngle)+i);
        %PSTH align on targt onset
        shadedErrorBar(PSTH_Time_FraOn_All,Good_FracOn_Loc_Mean{i,1},Good_FracOn_Loc_Sem{i,1},'lineprops',{[1    0    0]},'transparent',1,'patchSaturation',0.3);
        hold on
        shadedErrorBar(PSTH_Time_FraOn_All,Bad_FracOn_Loc_Mean{i,1},Bad_FracOn_Loc_Sem{i,1},'lineprops',{[0    0    1]},'transparent',1,'patchSaturation',0.3);

         hold on
         yrange2 = ylim;
        plot([0,0],yrange2,'--k','LineWidth',3);
        

       % xlim([-PlotSacPreBuffer, PlotSacPostBuffer]);%Keep in the same scale
       xlim([-200, 400]);%Keep in the same scale
       set(findall(gca, 'Type', 'Line'),'LineWidth',3);

        set(gca,'FontSize',15,'FontWeight','Bold');
        set(gca,'LineWidth',3);
    

        xlabel('Time from Fractal Onset(ms)','FontSize',15,'FontWeight','Bold');
        
%Raster align on saccade on
        subplot(4,length(UniqueAngle),2*length(UniqueAngle)+i);
         SpikeRaterPlot(Raster_SacOn_Good{i},StartLoc,Barlength,'r',LineWidth);
        hold on
        IncreaseBar = size(Raster_SacOn_Good{i},1);
 
        SpikeRaterPlot(Raster_SacOn_Bad{i},StartLoc+IncreaseBar,Barlength,'b',LineWidth);
        yrange = ylim;
        hold on
        plot([0,0],yrange,'--k','LineWidth',3);
       % xlim([-PlotSacPreBuffer, PlotSacPostBuffer]);
       xlim([-200, 200]);
        set(findall(gca, 'Type', 'Line'),'LineWidth',3);
        title(sprintf("Angle: %1.1d,Ecc: %1.1d",AngleEccorder{i,1}));


        subplot(4,length(UniqueAngle),3*length(UniqueAngle)+i);
        %PSTH align on sacccade onset
        shadedErrorBar(TimeSequence_SacOn_mean,Good_SacOn_Loc_Mean{i,j},Good_SacOn_Loc_Sem{i,j},'lineprops',{[1    0    0]},'transparent',1,'patchSaturation',0.3);

        hold on
        shadedErrorBar(TimeSequence_SacOn_mean,Bad_SacOn_Loc_Mean{i,j},Bad_SacOn_Loc_Sem{i,j},'lineprops',{[0    0    1]},'transparent',1,'patchSaturation',0.3);
       % yrange = ylim;
         plot([0,0],yrange2,'--k','LineWidth',3);

      %  xlim([-PlotSacPreBuffer, PlotSacPostBuffer]);%Keep in the same scale
       xlim([-200, 200]);
       set(findall(gca, 'Type', 'Line'),'LineWidth',3);

        set(gca,'FontSize',15,'FontWeight','Bold');
        set(gca,'LineWidth',3);
    

        xlabel('Time from Saccade Onset(ms)','FontSize',15,'FontWeight','Bold');




    end


end

%{
%%%%%%%%%Align on Fractal Onset%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(2,4,1); 
%shadedErrorBar(PSTH_Time_FraOn_All,MeanPSTH_FR_FractalOn_mean,MeanPSTH_FR_FractalOn_sem,'lineprops',{[0.9290    0.6940    0.1250]},'transparent',1,'patchSaturation',0.3);

%hold on
%ylim([0,max(MeanPSTH_FR_FractalOn_mean(PSTH_Time_FraOn_All<PlotSacPostBuffer & PSTH_Time_FraOn_All>-PlotSacPreBuffer))*1.05]);
%yrange=ylim;
%yrange(2)=max(MeanPSTH_FR_FractalOn_mean(PSTH_Time_FraOn_All<PlotSacPostBuffer & PSTH_Time_FraOn_All>-PlotSacPreBuffer))*1.05;

%shadedErrorBar([-PlotStartOffset:-PlotStartOffset+length(Eye_Eccentricity_FractalOn_mean)-1],Eye_Eccentricity_FractalOn_mean+yrange(2),Eye_Eccentricity_FractalOn_sem,'lineprops',{[0.9290    0.6940    0.1250]},'transparent',1,'patchSaturation',0.3);


shadedErrorBar([-PlotStartOffset:-PlotStartOffset+length(Eye_Ecc_Good_mean)-1],Eye_Ecc_Good_mean,Eye_Ecc_Good_sem,'lineprops',{[1    0    0]},'transparent',1,'patchSaturation',0.3);
hold on
shadedErrorBar([-PlotStartOffset:-PlotStartOffset+length(Eye_Ecc_Bad_mean)-1],Eye_Ecc_Bad_mean,Eye_Ecc_Bad_sem,'lineprops',{[0    0    1]},'transparent',1,'patchSaturation',0.3);
%ylim([0,yrange(2)+16]);




xlim([-PlotSacPreBuffer, PlotSacPostBuffer]);%Keep in the same scale
%yrange=ylim;
%ylim([0,yrange(2)]);
yrange=ylim;
plot([0,0],yrange,'--k');
set(findall(gca, 'Type', 'Line'),'LineWidth',3);

set(gca,'FontSize',15,'FontWeight','Bold');
set(gca,'LineWidth',3);
xlabel('Time from Fractal Onset(ms)','FontSize',15,'FontWeight','Bold');
ylabel('EyeEccentricity','FontSize',15,'FontWeight','Bold');

%%%%%%%%%Align on Saccade Onset%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(2,4,2);

%{
shadedErrorBar(TimeSequence_SacOn_mean,PSTH_SacOn_mean,PSTH_SacOn_sem,'lineprops',{[0.9290    0.6940    0.1250]},'transparent',1,'patchSaturation',0.3);

hold on

yrange=ylim;

shadedErrorBar([-PlotSacPreBuffer:EyeBinWidth: PlotSacPostBuffer],EyeEcc_Aligned_mean+yrange(2),EyeEcc_Aligned_sem,'lineprops',{[0.9290    0.6940    0.1250]},'transparent',1,'patchSaturation',0.3);
%}
shadedErrorBar([-PlotSacPreBuffer:-PlotSacPreBuffer+length(EyeEcc_Aligned_Good_mean)-1],EyeEcc_Aligned_Good_mean,EyeEcc_Aligned_Good_sem,'lineprops',{[1    0    0]},'transparent',1,'patchSaturation',0.3);
hold on
shadedErrorBar([-PlotSacPreBuffer:-PlotSacPreBuffer+length(EyeEcc_Aligned_Bad_mean)-1],EyeEcc_Aligned_Bad_mean,EyeEcc_Aligned_Bad_sem,'lineprops',{[0    0    1]},'transparent',1,'patchSaturation',0.3);
%ylim([0,yrange(2)+16]);

%yrange=ylim;
plot([0,0],yrange,'--k');
set(findall(gca, 'Type', 'Line'),'LineWidth',3);

set(gca,'FontSize',15,'FontWeight','Bold');
set(gca,'LineWidth',3);
xlabel('Time from Saccade Onset','FontSize',15,'FontWeight','Bold');
ylabel('Eye Eccentricity','FontSize',15,'FontWeight','Bold');

%%%%%%%%%Align on Saccade Offset%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(2,4,3);

%{
shadedErrorBar(TimeSequence_SacOff_mean,PSTH_SacOff_mean,PSTH_SacOff_sem,'lineprops',{[0.9290    0.6940    0.1250]},'transparent',1,'patchSaturation',0.3);
hold on

yrange=ylim;
yrange(2)=max(PSTH_SacOff_mean(TimeSequence_SacOff_mean<PlotSacPostBuffer & TimeSequence_SacOff_mean>-PlotSacPreBuffer))*1.05;
%}
%shadedErrorBar([-PlotSacPreBuffer:-PlotSacPreBuffer+length(EyeEcc_Aligned_SacOff_mean)-1],EyeEcc_Aligned_SacOff_mean+yrange(2),EyeEcc_Aligned_SacOff_sem,'lineprops',{[0.9290    0.6940    0.1250]},'transparent',1,'patchSaturation',0.3);
%ylim([0,yrange(2)+16]);

shadedErrorBar([-PlotSacPreBuffer:-PlotSacPreBuffer+length(EyeEcc_Aligned_SacOff_Good_mean)-1],EyeEcc_Aligned_SacOff_Good_mean,EyeEcc_Aligned_SacOff_Good_sem,'lineprops',{[1    0    0]},'transparent',1,'patchSaturation',0.3);
hold on
shadedErrorBar([-PlotSacPreBuffer:-PlotSacPreBuffer+length(EyeEcc_Aligned_SacOff_Bad_mean)-1],EyeEcc_Aligned_SacOff_Bad_mean,EyeEcc_Aligned_SacOff_Bad_sem,'lineprops',{[0    0    1]},'transparent',1,'patchSaturation',0.3);

xlim([-PlotSacPreBuffer, PlotSacPostBuffer]);%Keep in the same scale
yrange=ylim;
plot([0,0],yrange,'--k');

%Add reward and tgt off marker
plot([RwdTimeRelSacEndMean,RwdTimeRelSacEndMean],yrange,'--b');
plot([TgtOffTimeRelSacEndMean,TgtOffTimeRelSacEndMean],yrange,'--m');




set(findall(gca, 'Type', 'Line'),'LineWidth',3);
set(gca,'FontSize',15,'FontWeight','Bold');
set(gca,'LineWidth',3);

xlabel('Time from Saccade Offset','FontSize',15,'FontWeight','Bold');
ylabel('Eye Eccentricity','FontSize',15,'FontWeight','Bold');



%%%%%%%%%Align on Fractal Offset%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(2,4,4);

%{
shadedErrorBar(TimeSequence_SacOff_mean,PSTH_SacOff_mean,PSTH_SacOff_sem,'lineprops',{[0.9290    0.6940    0.1250]},'transparent',1,'patchSaturation',0.3);
hold on

yrange=ylim;
yrange(2)=max(PSTH_SacOff_mean(TimeSequence_SacOff_mean<PlotSacPostBuffer & TimeSequence_SacOff_mean>-PlotSacPreBuffer))*1.05;
%}
%shadedErrorBar([-PlotSacPreBuffer:-PlotSacPreBuffer+length(EyeEcc_Aligned_SacOff_mean)-1],EyeEcc_Aligned_SacOff_mean+yrange(2),EyeEcc_Aligned_SacOff_sem,'lineprops',{[0.9290    0.6940    0.1250]},'transparent',1,'patchSaturation',0.3);
%ylim([0,yrange(2)+16]);

shadedErrorBar([-PlotFracOffPreBuffer:-PlotFracOffPreBuffer+length(EyeEcc_Aligned_FracOff_Good_mean)-1],EyeEcc_Aligned_FracOff_Good_mean,EyeEcc_Aligned_FracOff_Good_sem,'lineprops',{[1    0    0]},'transparent',1,'patchSaturation',0.3);
hold on
shadedErrorBar([-PlotFracOffPreBuffer:-PlotFracOffPreBuffer+length(EyeEcc_Aligned_FracOff_Good_mean)-1],EyeEcc_Aligned_FracOff_Bad_mean,EyeEcc_Aligned_FracOff_Bad_sem,'lineprops',{[0    0    1]},'transparent',1,'patchSaturation',0.3);

xlim([-PlotFracOffPreBuffer, PlotFracOffPostBuffer]);%Keep in the same scale
yrange=ylim;
plot([0,0],yrange,'--k');






set(findall(gca, 'Type', 'Line'),'LineWidth',3);
set(gca,'FontSize',15,'FontWeight','Bold');
set(gca,'LineWidth',3);

xlabel('Time from Frac Offset','FontSize',15,'FontWeight','Bold');
ylabel('Eye Eccentricity','FontSize',15,'FontWeight','Bold');




%%%%%%%%%Align on Fractal Onset: Good vs Bad %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(2,4,5);
%%
shadedErrorBar(PSTH_Time_FraOn_All,PSTH_FR_Good_mean,PSTH_FR_Good_sem,'lineprops',{[1    0    0]},'transparent',1,'patchSaturation',0.3);

hold on
shadedErrorBar(PSTH_Time_FraOn_All,PSTH_FR_Bad_mean,PSTH_FR_Bad_sem,'lineprops',{[0    0    1]},'transparent',1,'patchSaturation',0.3);

yrange=ylim;
yrange(2)=max(MeanPSTH_FR_FractalOn_mean(PSTH_Time_FraOn_All<PlotSacPostBuffer & PSTH_Time_FraOn_All>-PlotSacPreBuffer))*1.05;

%shadedErrorBar([-PlotStartOffset:-PlotStartOffset+length(Eye_Ecc_Good_mean)-1],Eye_Ecc_Good_mean+yrange(2),Eye_Ecc_Good_sem,'lineprops',{[1    0    0]},'transparent',1,'patchSaturation',0.3);
%shadedErrorBar([-PlotStartOffset:-PlotStartOffset+length(Eye_Ecc_Bad_mean)-1],Eye_Ecc_Bad_mean+yrange(2),Eye_Ecc_Bad_sem,'lineprops',{[0    0    1]},'transparent',1,'patchSaturation',0.3);
%ylim([0,yrange(2)+16]);

xlim([-PlotSacPreBuffer, PlotSacPostBuffer]);%Keep in the same scale
yrange=ylim;
plot([0,0],yrange,'--k');

set(findall(gca, 'Type', 'Line'),'LineWidth',3);
set(gca,'FontSize',15,'FontWeight','Bold');
set(gca,'LineWidth',3);

xlabel('Time from Fractal Onset(ms)','FontSize',15,'FontWeight','Bold');
ylabel('Spikes/s','FontSize',15,'FontWeight','Bold');



%%%%%%%%%Align on Saccade Onset: Good vs Bad %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(2,4,6);
shadedErrorBar(TimeSequence_SacOn_mean,PSTH_SacOn_Good_mean,PSTH_SacOn_Good_sem,'lineprops',{[1    0    0]},'transparent',1,'patchSaturation',0.3);

hold on
shadedErrorBar(TimeSequence_SacOn_mean,PSTH_SacOn_Bad_mean,PSTH_SacOn_Bad_sem,'lineprops',{[0    0    1]},'transparent',1,'patchSaturation',0.3);

yrange=ylim;
yrange(2)=max(PSTH_SacOn_mean(TimeSequence_SacOn_mean<PlotSacPostBuffer & TimeSequence_SacOn_mean>-PlotSacPreBuffer))*1.05;
%{
shadedErrorBar([-PlotSacPreBuffer:-PlotSacPreBuffer+length(EyeEcc_Aligned_Good_mean)-1],EyeEcc_Aligned_Good_mean+yrange(2),EyeEcc_Aligned_Good_sem,'lineprops',{[1    0    0]},'transparent',1,'patchSaturation',0.3);
shadedErrorBar([-PlotSacPreBuffer:-PlotSacPreBuffer+length(EyeEcc_Aligned_Bad_mean)-1],EyeEcc_Aligned_Bad_mean+yrange(2),EyeEcc_Aligned_Bad_sem,'lineprops',{[0    0    1]},'transparent',1,'patchSaturation',0.3);
ylim([0,yrange(2)+16]);
%}
xlim([-PlotSacPreBuffer, PlotSacPostBuffer]);%Keep in the same scale
yrange=ylim;
plot([0,0],yrange,'--k');

set(findall(gca, 'Type', 'Line'),'LineWidth',3);
set(gca,'FontSize',15,'FontWeight','Bold');
set(gca,'LineWidth',3);

xlabel('Time from Saccade Onset','FontSize',15,'FontWeight','Bold');
ylabel('Spikes/s','FontSize',15,'FontWeight','Bold');

%%%%%%%%%Align on Saccade Offset: Good vs Bad %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(2,4,7);
shadedErrorBar(TimeSequence_SacOff_mean,PSTH_SacOff_Good_mean,PSTH_SacOff_Good_sem,'lineprops',{[1    0    0]},'transparent',1,'patchSaturation',0.3);

hold on
shadedErrorBar(TimeSequence_SacOff_mean,PSTH_SacOff_Bad_mean,PSTH_SacOff_Bad_sem,'lineprops',{[0    0    1]},'transparent',1,'patchSaturation',0.3);

yrange=ylim;
yrange(2)=max(PSTH_SacOff_mean(TimeSequence_SacOff_mean<PlotSacPostBuffer & TimeSequence_SacOff_mean>-PlotSacPreBuffer))*1.05;
%{
shadedErrorBar([-PlotSacPreBuffer:-PlotSacPreBuffer+length(EyeEcc_Aligned_SacOff_Good_mean)-1],EyeEcc_Aligned_SacOff_Good_mean+yrange(2),EyeEcc_Aligned_SacOff_Good_sem,'lineprops',{[1    0    0]},'transparent',1,'patchSaturation',0.3);
shadedErrorBar([-PlotSacPreBuffer:-PlotSacPreBuffer+length(EyeEcc_Aligned_SacOff_Bad_mean)-1],EyeEcc_Aligned_SacOff_Bad_mean+yrange(2),EyeEcc_Aligned_SacOff_Bad_sem,'lineprops',{[0    0    1]},'transparent',1,'patchSaturation',0.3);
ylim([0,yrange(2)+16]);
%}
xlim([-PlotSacPreBuffer, PlotSacPostBuffer]);%Keep in the same scale
yrange=ylim;
plot([0,0],yrange,'--k');

%Add reward and tgt off marker
plot([RwdTimeRelSacEndMean,RwdTimeRelSacEndMean],yrange,'--b');
plot([TgtOffTimeRelSacEndMean,TgtOffTimeRelSacEndMean],yrange,'--m');


set(findall(gca, 'Type', 'Line'),'LineWidth',3);
set(gca,'FontSize',15,'FontWeight','Bold');
set(gca,'LineWidth',3);

xlabel('Time from Saccade Offset','FontSize',15,'FontWeight','Bold');
ylabel('Spikes/s','FontSize',15,'FontWeight','Bold');

%%%%%%%%%Align on Fractal Offset: Good vs Bad %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(2,4,8);
shadedErrorBar(TimeSequence_FracOff_mean,PSTH_FracOff_Good_mean,PSTH_FracOff_Good_sem,'lineprops',{[1    0    0]},'transparent',1,'patchSaturation',0.3);

hold on
shadedErrorBar(TimeSequence_FracOff_mean,PSTH_FracOff_Bad_mean,PSTH_FracOff_Bad_sem,'lineprops',{[0    0    1]},'transparent',1,'patchSaturation',0.3);

yrange=ylim;
yrange(2)=max(PSTH_FracOff_mean(TimeSequence_FracOff_mean<PlotFracOffPreBuffer & TimeSequence_FracOff_mean>-PlotFracOffPreBuffer))*1.05;
%{
shadedErrorBar([-PlotSacPreBuffer:-PlotSacPreBuffer+length(EyeEcc_Aligned_SacOff_Good_mean)-1],EyeEcc_Aligned_SacOff_Good_mean+yrange(2),EyeEcc_Aligned_SacOff_Good_sem,'lineprops',{[1    0    0]},'transparent',1,'patchSaturation',0.3);
shadedErrorBar([-PlotSacPreBuffer:-PlotSacPreBuffer+length(EyeEcc_Aligned_SacOff_Bad_mean)-1],EyeEcc_Aligned_SacOff_Bad_mean+yrange(2),EyeEcc_Aligned_SacOff_Bad_sem,'lineprops',{[0    0    1]},'transparent',1,'patchSaturation',0.3);
ylim([0,yrange(2)+16]);
%}
xlim([-PlotFracOffPreBuffer, PlotFracOffPostBuffer]);%Keep in the same scale
yrange=ylim;
plot([0,0],yrange,'--k');



set(findall(gca, 'Type', 'Line'),'LineWidth',3);
set(gca,'FontSize',15,'FontWeight','Bold');
set(gca,'LineWidth',3);

xlabel('Time from Fractal Offset','FontSize',15,'FontWeight','Bold');
ylabel('Spikes/s','FontSize',15,'FontWeight','Bold');



%Add title

title_msg=sprintf(' %s, SpikeChannel: %d, TrialNum: %d ',FileName,SpikeChannel,TrialNum);
a=annotation('textbox',[0.25,0.68,0.3,0.3],'String',title_msg,'FitBoxToText','on');
a.FontSize=25;
a.FontWeight='Bold';
a.LineStyle='none';

%}
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
%
%Store figures in common folders
if TaskCode==OBJ_DIR_SAC_TASK  
ProtocolMarker='SatbleValue';
%ProtocolMarker='NovelObject_Low';

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



%}


ExistFlag=0;
%Load the old file if exist
if TaskCode==OBJ_DIR_SAC_TASK  
    ProtocolMarker='StableValue';
elseif  TaskCode==FLEXIBLE_VALUE_TASK 
    ProtocolMarker='FlexibleValue';
end
if exist(OutputFileName)

    load(OutputFileName);
    %Load OutputData into the memory
    if isfield(OutputData,ProtocolMarker)
    ExistFlag=1;
    end
end


%%%%%%%%%%%%%%%%%%%%OutputData For current analysis
OutputData_New=[];
if TaskCode==OBJ_DIR_SAC_TASK 
    OutputData_New.StableValue.TrialType='GoodTrials';
elseif  TaskCode==FLEXIBLE_VALUE_TASK 
    OutputData_New.FlexibleValue.TrialType='GoodTrials';
end

%%%%%%%%%Save overall scene response characteristics%%%%%%%%%%%%%%%%%%%%%%%
%Organize the data
NameStr={'PSTH_TgtOn_Good_RF','PSTH_TgtOn_Bad_RF','PSTH_TgtOn_Good_nonRF','PSTH_TgtOn_Bad_nonRF',...
    'PSTH_TgtOn_Good_RF_z','PSTH_TgtOn_Bad_RF_z',...
    'PSTH_TgtOn_Good_nonRF_z','PSTH_TgtOn_Bad_nonRF_z',...
    'TgtGood_RF','TgtBad_RF','TgtGood_nonRF','TgtBad_nonRF',...
    'PSTH_SacOn_Good_RF','PSTH_SacOn_Bad_RF','PSTH_SacOn_Good_nonRF','PSTH_SacOn_Bad_nonRF',...
    'PSTH_SacOn_Good_RF_z','PSTH_SacOn_Bad_RF_z',...
    'PSTH_SacOn_Good_nonRF_z','PSTH_SacOn_Bad_nonRF_z',...
    'SacGood_RF','SacBad_RF','SacGood_nonRF','SacBad_nonRF',...    
    'PSTH_TgtOn_Prefer_RF','PSTH_TgtOn_nonPrefer_RF','PSTH_SacOn_Prefer_RF','PSTH_SacOn_nonPrefer_RF',...
    'PSTH_TgtOn_Prefer_nonRF','PSTH_TgtOn_nonPrefer_nonRF','PSTH_SacOn_Prefer_nonRF','PSTH_SacOn_nonPrefer_nonRF',...
    'Time_TgtOn','Time_SacOn','p_Tgt_RF','p_Sac_RF','p_Tgt_nonRF','p_Sac_nonRF'...
    'PSTH_TgtOn_Prefer_RF_z','PSTH_TgtOn_nonPrefer_RF_z','PSTH_SacOn_Prefer_RF_z','PSTH_SacOn_nonPrefer_RF_z',...
    'PSTH_TgtOn_Prefer_nonRF_z','PSTH_TgtOn_nonPrefer_nonRF_z','PSTH_SacOn_Prefer_nonRF_z','PSTH_SacOn_nonPrefer_nonRF_z',...
    'TgtAUC_RF','TgtAUC_nonRF','p_AUC_Tgt_RF','p_AUC_Tgt_nonRF ',...
    'p_response_interval_RF','p_response_interval_nonRF','TgtAUC_interval_RF','TgtAUC_interval_nonRF','p_interval_RF','p_interval_nonRF'

  

};

DataLib={PSTH_TgtOn_Good_RF,PSTH_TgtOn_Bad_RF,PSTH_TgtOn_Good_nonRF,PSTH_TgtOn_Bad_nonRF,...
    PSTH_TgtOn_Good_RF_z,PSTH_TgtOn_Bad_RF_z,...
    PSTH_TgtOn_Good_nonRF_z,PSTH_TgtOn_Bad_nonRF_z,...
    TgtGood_RF,TgtBad_RF,TgtGood_nonRF,TgtBad_nonRF...
    PSTH_SacOn_Good_RF,PSTH_SacOn_Bad_RF,PSTH_SacOn_Good_nonRF,PSTH_SacOn_Bad_nonRF,...
    PSTH_SacOn_Good_RF_z,PSTH_SacOn_Bad_RF_z,...
    PSTH_SacOn_Good_nonRF_z,PSTH_SacOn_Bad_nonRF_z,...
    SacGood_RF,SacBad_RF,SacGood_nonRF,SacBad_nonRF...
    PSTH_TgtOn_Prefer_RF,PSTH_TgtOn_nonPrefer_RF,PSTH_SacOn_Prefer_RF,PSTH_SacOn_nonPrefer_RF,...
    PSTH_TgtOn_Prefer_nonRF,PSTH_TgtOn_nonPrefer_nonRF,PSTH_SacOn_Prefer_nonRF,PSTH_SacOn_nonPrefer_nonRF,...
    Time_TgtOn,Time_SacOn,p_Tgt_RF,p_Sac_RF,p_Tgt_nonRF,p_Sac_nonRF...
    PSTH_TgtOn_Prefer_RF_z,PSTH_TgtOn_nonPrefer_RF_z,PSTH_SacOn_Prefer_RF_z,PSTH_SacOn_nonPrefer_RF_z,...
    PSTH_TgtOn_Prefer_nonRF_z,PSTH_TgtOn_nonPrefer_nonRF_z,PSTH_SacOn_Prefer_nonRF_z,PSTH_SacOn_nonPrefer_nonRF_z,...
    TgtAUC_RF,TgtAUC_nonRF,p_AUC_Tgt_RF,p_AUC_Tgt_nonRF,...
    p_response_interval_RF,p_response_interval_nonRF,TgtAUC_interval_RF,TgtAUC_interval_nonRF,p_interval_RF,p_interval_nonRF...
    

};


DataStamp=containers.Map(NameStr,DataLib);




if TaskCode==OBJ_DIR_SAC_TASK 
    OutputData_New.StableValue.Task=TaskType;
    OutputData_New.StableValue.TaskCode=TaskCode;
    OutputData_New.StableValue.DataStamp=DataStamp;


%Store figures
%{
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
%}
if ExistFlag
%Compare the old one with the new one
Task_Old=OutputData.StableValue.TaskCode;

     if Task_Old~=TaskCode
        %Not the same task, add the new dataset into the old one
        OutputData.StableValue(numel(OutputData.StableValue)+1)=OutputData_New.StableValue;
     else
         %If the same task,the same protocol,replace the old one with the new one
         
         OutputData.StableValue=OutputData_New.StableValue;
         
    
    
     end
else
    %Output the current file
   OutputData.StableValue=OutputData_New.StableValue; 
end

elseif  TaskCode==FLEXIBLE_VALUE_TASK 
        OutputData_New.FlexibleValue.Task=TaskType;
        OutputData_New.FlexibleValue.TaskCode=TaskCode;
        OutputData_New.FlexibleValue.DataStamp=DataStamp;


%Store figures
%{
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
%}
if ExistFlag
%Compare the old one with the new one
Task_Old=OutputData.FlexibleValue.TaskCode;

     if Task_Old~=TaskCode
        %Not the same task, add the new dataset into the old one
        OutputData.FlexibleValue(numel(OutputData.FlexibleValue)+1)=OutputData_New.FlexibleValue;
     else
         %If the same task,the same protocol,replace the old one with the new one
         
         OutputData.FlexibleValue=OutputData_New.FlexibleValue;
         
    
    
     end
else
    %Output the current file
   OutputData.FlexibleValue=OutputData_New.FlexibleValue; 
end


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