function  TargetSaccadeAnalysis(Data, StartTrial, EndTrial,ShowFigureFlag,OutputFlag);
%Function to analysis neural responses and target property and saccades
%properties
%Xuefei Yu 05122020
ECode;


SegEvent=[TGTCD];%Event Flow of interest
%SegEventEnd=[SHOWFIXCD,FIXOFF,TGTOFF,RWDOFFCD ,NaN,EXTRAFP];%Event end
SegEventEnd=[TGTOFF];%Event end
SegEventStr={'TargetOn'};
AmpMatchFlag=0;%1:Amplitude match on ;0:Amplitude Match off


PlotStartOffset=[200];%in ms
PlotStopOffset=[100];%in ms

ExamplePlot=0;%PLots for presentations, turn off for regular use

ResponseLatencyThreshold=3;

%Parameter to define the sliding window to measure firing rates
StepSize=5;%in ms
BinWidth=15;%in ms


%Criterias for selecting saccades
V_Threshold=40;
ContinueBin=5;
ISI_Threshold=20;

%Time bin pre saccade and time bin offset after saccade onset
PlotSacPreBuffer=200;%in ms
PlotSacPostBuffer=400;%in ms




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

Selection=StartTrial:EndTrial;

%Abnormal check %Especially For Robin081120FORAGEXPLO4.rst
CheckCode=SCENE_ON;%This code should appear once in each trial
AbnormalTrials=AbnormalTrialCheck(Data,CheckCode);

Selection=Selection(~ismember(Selection,AbnormalTrials));
TrialNum=length(Selection);


%FixationPointIndex=Data.FixPointIndex(Selection,:);
%FixLoc=Data.FixLoc(Selection,:,:);

%Target Information
TargetIndex=Data.TargetIndex(Selection,2:end,1);
TargetName=Data.TargetIndex(Selection,2:end,2);
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
GoodObjectIndex=mod(TargetIndex_Line,2)==1;


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


%Recalibrate eyes according to the extra FP
EyeCalInfo=Data.EyeCalInfo;%Each row: X_Gain,X_Bias; Y_Gain; Y_bias; X=(raw+X_bias)*X_Gain
EyeCalInfo_L=EyeCalInfo(1,:);

ParaLib=Data.ParaLib;
ExtraFP_Info=ParaLib('ExtraFixPointID,Style,Size,Loc_X,Loc_Y');
ExtraFP_Loc=ExtraFP_Info(Selection,end-1:end);


%Get the information about the extra FP
%ExtraFPOnset=FindOutTime(EventChannel,EventTimeChannel,EXTRAFP);
ExtraEyeIn=FindOutTime(EventChannel,EventBin,EXTRAFIXSTARTCD);
ExtraEyeOff=FindOutTime(EventChannel,EventBin,EXTRAFIXOFF);

EyeInBin_X=SelectEyeInterval(EyeChannel_X_L,ExtraEyeIn,ExtraEyeOff);
EyeInBin_X=EyeInBin_X{1};
EyeInBin_Y=SelectEyeInterval(EyeChannel_Y_L,ExtraEyeIn,ExtraEyeOff);
EyeInBin_Y=EyeInBin_Y{1};
EyeDataForCal={EyeInBin_X,EyeInBin_Y};

EyeDataRaw={EyeChannel_X_L,EyeChannel_Y_L};
CalPlotCheck=0;%Examine the calibration results
CaledEye=EyeCalibration(EyeDataRaw,EyeDataForCal,EyeCalInfo_L,ExtraFP_Loc,CalPlotCheck);
EyeChannel_X_L=CaledEye.EyeTrace{1};
EyeChannel_Y_L=CaledEye.EyeTrace{2};


SpikeChannel=Data.SpikeChannel;
SpikeTime=squeeze(Data.SpikeTimeData(:,:,Selection))';
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


%Get mean and sem of eye eccentricity
%Get eye trance distance during the interval
StartTimeBin=FindOutTime(EventChannel,EventTimeChannel,SegEvent)-PlotStartOffset;
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

EyeEcc_Aligned_FraOff=AlignEyeData(Eye_Eccentricity_FractalOn,EyeTimeII,RelativeTimeMarker(:,2),[-PlotSacPreBuffer,PlotSacPostBuffer]);


%%
%I. Sort neural responses according to the saccade latency
SepGroupRT=3;
GroupStr={'Highest','Medium','Lowest'};

%1. For neural responses aligned with fractal onset
%RT_First_Target=SacEnd_First_Target;

FacOn_temp=SortDataAccordingToGroup(MeanPSTH_FR_FractalOn,RT_First_Target,SepGroupRT,'descend');


PSTH_FraOn_RTGroup_mean=FacOn_temp.mean;
PSTH_FraOn_RTGroup_sem=FacOn_temp.sem;

RT_Group_mean=FacOn_temp.MarkerMean;
RT_Group_sem=FacOn_temp.MarkerSem;

NumInGroupRT=FacOn_temp.Num;
RT_SortIndex=FacOn_temp.SortIndex;





%Also need the eye trace

FacOn_Eyetemp=SortDataAccordingToGroup(Eye_Eccentricity_FractalOn,RT_First_Target,SepGroupRT,'descend');
EyeEcc_FraOn_RTGroup_mean=FacOn_Eyetemp.mean;
EyeEcc_FraOn_RTGroup_sem=FacOn_Eyetemp.sem;

%2. For neural responses aligned with Saccade onset

SacOn_temp=SortDataAccordingToGroup(PSTH_SacOn,RT_First_Target,SepGroupRT,'descend');

PSTH_SacOn_RTGroup_mean=SacOn_temp.mean;
PSTH_SacOn_RTGroup_sem=SacOn_temp.sem;

%Also need the eye trace

SaccadeOn_Eyetemp=SortDataAccordingToGroup(EyeEcc_Aligned,RT_First_Target,SepGroupRT,'descend');
EyeEcc_SacOn_RTGroup_mean=SaccadeOn_Eyetemp.mean;
EyeEcc_SacOn_RTGroup_sem=SaccadeOn_Eyetemp.sem;

%3. For neural responses aligned with Saccade Offset

SacOff_temp=SortDataAccordingToGroup(PSTH_SacOff,RT_First_Target,SepGroupRT,'descend');

PSTH_SacOff_RTGroup_mean=SacOff_temp.mean;
PSTH_SacOff_RTGroup_sem=SacOff_temp.sem;

%Also need the eye trace

SaccadeOff_Eyetemp=SortDataAccordingToGroup(EyeEcc_Aligned_SacOff,RT_First_Target,SepGroupRT,'descend');
EyeEcc_SacOff_RTGroup_mean=SaccadeOff_Eyetemp.mean;
EyeEcc_SacOff_RTGroup_sem=SaccadeOff_Eyetemp.sem;


%II. Sort neural responses according to the Saccade Amplitude
SepGroupAmp=SepGroupRT;
SepGroupAmpNum=3;
%SepGroupAmp={[12,10]};
%1. For neural responses aligned with fractal onset

FacOn_temp=SortDataAccordingToGroup(MeanPSTH_FR_FractalOn,SaccadeAmplitude_First,SepGroupAmp,'descend');


PSTH_FraOn_AmpGroup_mean=FacOn_temp.mean;
PSTH_FraOn_AmpGroup_sem=FacOn_temp.sem;

Amp_Group_mean=FacOn_temp.MarkerMean;
Amp_Group_sem=FacOn_temp.MarkerSem;

NumInGroupAmp=FacOn_temp.Num;
Amp_SortIndex=FacOn_temp.SortIndex;

%Also need the eye trace

FracOn_Eyetemp=SortDataAccordingToGroup(Eye_Eccentricity_FractalOn,SaccadeAmplitude_First,SepGroupAmp,'descend');
EyeEcc_FracOn_AmpGroup_mean=FracOn_Eyetemp.mean;
EyeEcc_FracOn_AmpGroup_sem=FracOn_Eyetemp.sem;

%2. For neural responses aligned with Saacde onset

SacOn_temp=SortDataAccordingToGroup(PSTH_SacOn,SaccadeAmplitude_First,SepGroupAmp,'descend');

PSTH_SacOn_AmpGroup_mean=SacOn_temp.mean;
PSTH_SacOn_AmpGroup_sem=SacOn_temp.sem;

%Also need the eye trace

SaccadeOn_Eyetemp=SortDataAccordingToGroup(EyeEcc_Aligned,SaccadeAmplitude_First,SepGroupAmp,'descend');
EyeEcc_SacOn_AmpGroup_mean=SaccadeOn_Eyetemp.mean;
EyeEcc_SacOn_AmpGroup_sem=SaccadeOn_Eyetemp.sem;


%3. For neural responses aligned with Saccade Offset

SacOff_temp=SortDataAccordingToGroup(PSTH_SacOff,SaccadeAmplitude_First,SepGroupAmp,'descend');

PSTH_SacOff_AmpGroup_mean=SacOff_temp.mean;
PSTH_SacOff_AmpGroup_sem=SacOff_temp.sem;

%Also need the eye trace

SaccadeOff_Eyetemp=SortDataAccordingToGroup(EyeEcc_Aligned_SacOff,SaccadeAmplitude_First,SepGroupAmp,'descend');
EyeEcc_SacOff_AmpGroup_mean=SaccadeOff_Eyetemp.mean;
EyeEcc_SacOff_AmpGroup_sem=SaccadeOff_Eyetemp.sem;




%%
%%%%%%%%%%Analysis the neural responses (aligned on saccade onset) and
%%%%%%%%%%saccade amplitude on case-by-case basis
%Use spike rate
%PreSaccade
PreSacDur=[-100,0];%in ms
SpikeRate_PreSac=nanmean(PSTH_SacOn(:,[TimeSequence_SacOn_mean<=PreSacDur(2) &TimeSequence_SacOn_mean>=PreSacDur(1)]),2);
%Group by amplitude( Sac On)
SacOn_tempPre=SortDataAccordingToGroup(SpikeRate_PreSac,SaccadeAmplitude_First,SepGroupAmp,'descend');
SpikeRate_PreSac_Amp_mean=SacOn_tempPre.mean;
SpikeRate_PreSac_Amp_sem=SacOn_tempPre.sem;

%For stats:

SpikeRate_PreSac_Amp_line=SacOn_tempPre.LineData;
SpikeRate_PreSac_Amp_lineIndicator=SacOn_tempPre.GrouopIndicator;

ANOVA_p_Amp_pre=anova1(SpikeRate_PreSac_Amp_line,SpikeRate_PreSac_Amp_lineIndicator,'off');


%During Saccade:
SacIn=[0*ones(length(SacEnd_First_Target),1),SacEnd_First_Target'-RT_First_Target'];
SpikeRate_DurSac=SelectFringRateInterval(PSTH_SacOn,TimeSequence_SacOn_mean,SacIn)';

%Group by amplitude( DurSac)
SacOn_tempDur=SortDataAccordingToGroup(SpikeRate_DurSac,SaccadeAmplitude_First,SepGroupAmp,'descend');
SpikeRate_DurSac_Amp_mean=SacOn_tempDur.mean;
SpikeRate_DurSac_Amp_sem=SacOn_tempDur.sem;
%For stats:

SpikeRate_DurSac_Amp_line=SacOn_tempDur.LineData;
SpikeRate_DurSac_Amp_lineIndicator=SacOn_tempDur.GrouopIndicator;

ANOVA_p_Amp_dur=anova1(SpikeRate_DurSac_Amp_line,SpikeRate_DurSac_Amp_lineIndicator,'off');


%After Saccade:
PostSacDur=[0,100];%in ms
SpikeRate_PostSac=nanmean(PSTH_SacOff(:,[TimeSequence_SacOff_mean<=PostSacDur(2) &TimeSequence_SacOff_mean>=PostSacDur(1)]),2);


%Group by amplitude( PostSac)
SacOn_tempPost=SortDataAccordingToGroup(SpikeRate_PostSac,SaccadeAmplitude_First,SepGroupAmp,'descend');
SpikeRate_PostSac_Amp_mean=SacOn_tempPost.mean;
SpikeRate_PostSac_Amp_sem=SacOn_tempPost.sem;
%For stats:

SpikeRate_PostSac_Amp_line=SacOn_tempPost.LineData;
SpikeRate_PostSac_Amp_lineIndicator=SacOn_tempPost.GrouopIndicator;

ANOVA_p_Amp_post=anova1(SpikeRate_PostSac_Amp_line,SpikeRate_PostSac_Amp_lineIndicator,'off');



%%Fitting
%PreSaccade:
 [B,BINT,R,RINT,STATS] = regress(SpikeRate_PreSac,[ones(length(SaccadeAmplitude_First),1),SaccadeAmplitude_First']);
 r_preSac=sqrt(STATS(1));
 p_preSac=STATS(3);
 %Line
 LinepreSac=B;
 %Good v.s. Bad
 [B,BINT,R,RINT,STATS] = regress(SpikeRate_PreSac(GoodObjectIndex),[ones(length(SaccadeAmplitude_First((GoodObjectIndex))),1),SaccadeAmplitude_First((GoodObjectIndex))']);
 r_preSac_Good=sqrt(STATS(1));
 p_preSac_Good=STATS(3);
 
 [B,BINT,R,RINT,STATS] = regress(SpikeRate_PreSac(~GoodObjectIndex),[ones(length(SaccadeAmplitude_First((~GoodObjectIndex))),1),SaccadeAmplitude_First((~GoodObjectIndex))']);
 r_preSac_Bad=sqrt(STATS(1));
 p_preSac_Bad=STATS(3);
 
 
 
 %DurSaccade:
 [B,BINT,R,RINT,STATS] = regress(SpikeRate_DurSac,[ones(length(SaccadeAmplitude_First),1),SaccadeAmplitude_First']);
 r_durSac=sqrt(STATS(1));
 p_durSac=STATS(3);
 %Line
 LinedurSac=B;
 
 %Good v.s. Bad
 [B,BINT,R,RINT,STATS] = regress(SpikeRate_DurSac(GoodObjectIndex),[ones(length(SaccadeAmplitude_First((GoodObjectIndex))),1),SaccadeAmplitude_First((GoodObjectIndex))']);
 r_durSac_Good=sqrt(STATS(1));
 p_durSac_Good=STATS(3);
 
 
 [B,BINT,R,RINT,STATS] = regress(SpikeRate_DurSac(~GoodObjectIndex),[ones(length(SaccadeAmplitude_First((~GoodObjectIndex))),1),SaccadeAmplitude_First((~GoodObjectIndex))']);
 r_durSac_Bad=sqrt(STATS(1));
 p_durSac_Bad=STATS(3);
 
 
 
 
 %PostSaccade:
 [B,BINT,R,RINT,STATS] = regress(SpikeRate_PostSac,[ones(length(SaccadeAmplitude_First),1),SaccadeAmplitude_First']);
 r_postSac=sqrt(STATS(1));
 p_postSac=STATS(3);
 %Line
 LinepostSac=B;
 %Good v.s. Bad
 [B,BINT,R,RINT,STATS] = regress(SpikeRate_PostSac(GoodObjectIndex),[ones(length(SaccadeAmplitude_First((GoodObjectIndex))),1),SaccadeAmplitude_First((GoodObjectIndex))']);
 r_postSac_Good=sqrt(STATS(1));
 p_postSac_Good=STATS(3);
 
 [B,BINT,R,RINT,STATS] = regress(SpikeRate_PostSac(~GoodObjectIndex),[ones(length(SaccadeAmplitude_First((~GoodObjectIndex))),1),SaccadeAmplitude_First((~GoodObjectIndex))']);
 r_postSac_Bad=sqrt(STATS(1));
 p_postSac_Bad=STATS(3);
 %Saccade Vector Tuning
 %Group responses according to saccade vectors
 SeparationPoint=[22.5:45:360];
 
 SaccadeVectorTuning=PoloarHistGroup(SaccadeAngle_First,[SpikeRate_PreSac,SpikeRate_DurSac,SpikeRate_PostSac],SeparationPoint);
 SaccadeVector=SaccadeVectorTuning.Center;
 FiringRate_PreDurPost_mean=SaccadeVectorTuning.DataMean;
 FiringRate_PreDurPost_sem=SaccadeVectorTuning.DataSEM;
 SacVecTuning_p=SaccadeVectorTuning.ANOVA_p;
 
 PrefVector=SaccadeVectorTuning.PrefVector;
PrefVector_Amp=SaccadeVectorTuning.PrefVector_Amp;
 
 
 
 %}

%{
 %%Use Spike Count
SpikeTimeSelected_SacOn=SpikeTimeSelected_TgtOn-RT_First_Target';%Align on saccade on

%Measure spike counts according to the pre-saccade time,Dur-Sac,
%Post-Saccade

%PreSaccade
PreSacDur=[-200,0];%in ms
SpikeCount_PreSac=sum(SpikeTimeSelected_SacOn>=PreSacDur(1) & SpikeTimeSelected_SacOn<=PreSacDur(2),2);
SpikeRate_PreSac=SpikeCount_PreSac/(diff(PreSacDur)*(10^-3));

%During Saccade:
SacIn=[0*ones(length(SacEnd_First_Target),1),SacEnd_First_Target'-RT_First_Target'];

SpikeCount_DurSac=sum(SpikeTimeSelected_SacOn>=SacIn(:,1) & SpikeTimeSelected_SacOn<=SacIn(:,2),2);
SpikeRate_DurSac=SpikeCount_DurSac./(SacIn(:,2)*10^-3);

%After Saccade:
SpikeTimeSelected_SacOff=SpikeTimeSelected_TgtOn-SacEnd_First_Target';%Align on saccade off
PostSacDur=[0,100];%in ms
SpikeCount_AfterSac=sum(SpikeTimeSelected_SacOff>=PostSacDur(1) & SpikeTimeSelected_SacOff<=PostSacDur(2),2);
SpikeRate_PostSac=SpikeCount_AfterSac./(diff(PostSacDur)*10^-3);
%}
 
 if AmpMatchFlag
 
 %Separate good and object with amplitude match
 
SacAmp_Good=SaccadeAmplitude_First(GoodObjectIndex);
SacAmp_Bad=SaccadeAmplitude_First(~GoodObjectIndex);

nanremover_Good=~isnan(SacAmp_Good);
nanremover_Bad=~isnan(SacAmp_Bad);

SacAmp_Good=SacAmp_Good(nanremover_Good);
SacAmp_Bad=SacAmp_Bad(nanremover_Bad);


SpikeRate_DurSac_Goodtemp=SpikeRate_DurSac(GoodObjectIndex);
SpikeRate_DurSac_Badtemp=SpikeRate_DurSac(~GoodObjectIndex);


SpikeRate_DurSac_Goodtemp=SpikeRate_DurSac_Goodtemp(nanremover_Good);
SpikeRate_DurSac_Badtemp=SpikeRate_DurSac_Badtemp(nanremover_Bad);

SpikeRate_PreSac_Goodtemp=SpikeRate_PreSac(GoodObjectIndex);
SpikeRate_PreSac_Badtemp=SpikeRate_PreSac(~GoodObjectIndex);


SpikeRate_PreSac_Goodtemp=SpikeRate_PreSac_Goodtemp(nanremover_Good);
SpikeRate_PreSac_Badtemp=SpikeRate_PreSac_Badtemp(nanremover_Bad);



SpikeRate_PostSac_Goodtemp=SpikeRate_PostSac(GoodObjectIndex);
SpikeRate_PostSac_Badtemp=SpikeRate_PostSac(~GoodObjectIndex);


SpikeRate_PostSac_Goodtemp=SpikeRate_PostSac_Goodtemp(nanremover_Good);
SpikeRate_PostSac_Badtemp=SpikeRate_PostSac_Badtemp(nanremover_Bad);


SacAmp_Good_mean=mean(SacAmp_Good);
SacAmp_Bad_mean=mean(SacAmp_Bad);

MaxLoopNum=1000;
LoopNum=0;
GroupNum=ceil(min(length(SacAmp_Good),length(SacAmp_Bad))*1/2);

while abs(SacAmp_Good_mean-SacAmp_Bad_mean)>0.01 && LoopNum<=MaxLoopNum
   Selection_Good=randperm(length(SacAmp_Good),GroupNum); 
   Selection_Bad=randperm(length(SacAmp_Bad),GroupNum); 
   
   SacAmp_Good_mean=mean(SacAmp_Good(Selection_Good));
   SacAmp_Bad_mean=mean(SacAmp_Bad(Selection_Bad));
   
  LoopNum=LoopNum+1;
  if LoopNum>MaxLoopNum
      disp('Amplitude match failed')
  end
    
end
%}
%{
Selection_Good=find(SacAmp_Good>10 & SacAmp_Good<20);
Selection_Bad=find(SacAmp_Bad>10 & SacAmp_Bad<20);

Selection_Good=Selection_Good(randperm(length(Selection_Good),GroupNum));
Selection_Bad=Selection_Bad(randperm(length(Selection_Bad),GroupNum));
%}

 SacAmp_Good_match=SacAmp_Good(Selection_Good);
  SacAmp_Bad_match=SacAmp_Bad(Selection_Bad);
  
  
  SpikeRate_PreSac_Goodmatch=SpikeRate_PreSac_Goodtemp(Selection_Good);
  SpikeRate_PreSac_Badmatch=SpikeRate_PreSac_Badtemp(Selection_Bad);
  
  SpikeRate_DurSac_Goodmatch=SpikeRate_DurSac_Goodtemp(Selection_Good);
  SpikeRate_DurSac_Badmatch=SpikeRate_DurSac_Badtemp(Selection_Bad);
  
   SpikeRate_PostSac_Goodmatch=SpikeRate_PostSac_Goodtemp(Selection_Good);
  SpikeRate_PostSac_Badmatch=SpikeRate_PostSac_Badtemp(Selection_Bad);
  
  SpikeRate_PreSac_Match_mean=[mean(SpikeRate_PreSac_Goodmatch),mean(SpikeRate_PreSac_Badmatch)];
  SpikeRate_PreSac_Match_sem=[std(SpikeRate_PreSac_Goodmatch)/sqrt(length(SpikeRate_PreSac_Goodmatch)),mean(SpikeRate_PreSac_Badmatch)/sqrt(length(SpikeRate_PreSac_Badmatch))];
  
   SpikeRate_DurSac_Match_mean=[mean(SpikeRate_DurSac_Goodmatch),mean(SpikeRate_DurSac_Badmatch)];
  SpikeRate_DurSac_Match_sem=[std(SpikeRate_DurSac_Goodmatch)/sqrt(length(SpikeRate_DurSac_Goodmatch)),mean(SpikeRate_DurSac_Badmatch)/sqrt(length(SpikeRate_DurSac_Badmatch))];
  
   SpikeRate_PostSac_Match_mean=[mean(SpikeRate_PostSac_Goodmatch),mean(SpikeRate_PostSac_Badmatch)];
  SpikeRate_PostSac_Match_sem=[std(SpikeRate_PostSac_Goodmatch)/sqrt(length(SpikeRate_PostSac_Goodmatch)),mean(SpikeRate_PostSac_Badmatch)/sqrt(length(SpikeRate_PostSac_Badmatch))];
  
  %{
   [h,p_Pre_AmpMatch]=ttest2(SpikeRate_PreSac_Goodmatch,SpikeRate_PreSac_Badmatch);
  [h,p_Dur_AmpMatch]=ttest2(SpikeRate_DurSac_Goodmatch,SpikeRate_DurSac_Badmatch);
  [h,p_Post_AmpMatch]=ttest2(SpikeRate_PostSac_Goodmatch,SpikeRate_PostSac_Badmatch);
  %}
  
  [p_Pre_AmpMatch]=ranksum(SpikeRate_PreSac_Goodmatch,SpikeRate_PreSac_Badmatch);
  [p_Dur_AmpMatch]=ranksum(SpikeRate_DurSac_Goodmatch,SpikeRate_DurSac_Badmatch);
  [p_Post_AmpMatch]=ranksum(SpikeRate_PostSac_Goodmatch,SpikeRate_PostSac_Badmatch);

 end

%Align Eye trace with saccade offset
%{
EyeEcc_Aligned_SacOff=AlignEyeData(Eye_Eccentricity_Tgt,EyeTime,SacEnd_First_Target,[-PlotSacPreBuffer,PlotSacPostBuffer]);

EyeEcc_Aligned_SacOff_mean=nanmean(EyeEcc_Aligned_SacOff,1);
EyeEcc_Aligned_SacOff_sem=nanstd(EyeEcc_Aligned_SacOff,[],1)./sqrt(sum(~isnan(EyeEcc_Aligned_SacOff),1));
%}
%PSTH_SacOff
%TimeSequence_SacOff]


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


%%
if ShowFigureFlag
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Figure 1: Average firing rate
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%for each fractal
figtitlestr{1}='PSTHAlignOnSacOnset';

fig{1}=PrepareFigure(110,'w',[50,100, 1200,800],'Name',figtitlestr{1});


%%%%%%%%%Align on Fractal Onset%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(2,3,1); 
shadedErrorBar(PSTH_Time_FraOn_All,MeanPSTH_FR_FractalOn_mean,MeanPSTH_FR_FractalOn_sem,'lineprops',{[0.9290    0.6940    0.1250]},'transparent',1,'patchSaturation',0.3);

hold on
%ylim([0,max(MeanPSTH_FR_FractalOn_mean(PSTH_Time_FraOn_All<PlotSacPostBuffer & PSTH_Time_FraOn_All>-PlotSacPreBuffer))*1.05]);
yrange=ylim;
yrange(2)=max(MeanPSTH_FR_FractalOn_mean(PSTH_Time_FraOn_All<PlotSacPostBuffer & PSTH_Time_FraOn_All>-PlotSacPreBuffer))*1.05;

shadedErrorBar([-PlotStartOffset:-PlotStartOffset+length(Eye_Eccentricity_FractalOn_mean)-1],Eye_Eccentricity_FractalOn_mean+yrange(2),Eye_Eccentricity_FractalOn_sem,'lineprops',{[0.9290    0.6940    0.1250]},'transparent',1,'patchSaturation',0.3);
xlim([-PlotSacPreBuffer, PlotSacPostBuffer]);%Keep in the same scale
yrange=ylim;
ylim([0,yrange(2)]);
yrange=ylim;
plot([0,0],yrange,'--k');
set(findall(gca, 'Type', 'Line'),'LineWidth',3);

set(gca,'FontSize',15,'FontWeight','Bold');
set(gca,'LineWidth',3);
xlabel('Time from Fractal Onset(ms)','FontSize',25,'FontWeight','Bold');
ylabel('Spikes/s','FontSize',25,'FontWeight','Bold');

%%%%%%%%%Align on Saccade Onset%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(2,3,2);


shadedErrorBar(TimeSequence_SacOn_mean,PSTH_SacOn_mean,PSTH_SacOn_sem,'lineprops',{[0.9290    0.6940    0.1250]},'transparent',1,'patchSaturation',0.3);

hold on

yrange=ylim;

shadedErrorBar([-PlotSacPreBuffer:EyeBinWidth: PlotSacPostBuffer],EyeEcc_Aligned_mean+yrange(2),EyeEcc_Aligned_sem,'lineprops',{[0.9290    0.6940    0.1250]},'transparent',1,'patchSaturation',0.3);
yrange=ylim;
plot([0,0],yrange,'--k');
set(findall(gca, 'Type', 'Line'),'LineWidth',3);

set(gca,'FontSize',15,'FontWeight','Bold');
set(gca,'LineWidth',3);
xlabel('Time from Saccade Onset','FontSize',25,'FontWeight','Bold');
ylabel('Spikes/s','FontSize',25,'FontWeight','Bold');

%%%%%%%%%Align on Saccade Offset%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(2,3,3);


shadedErrorBar(TimeSequence_SacOff_mean,PSTH_SacOff_mean,PSTH_SacOff_sem,'lineprops',{[0.9290    0.6940    0.1250]},'transparent',1,'patchSaturation',0.3);
hold on

yrange=ylim;
yrange(2)=max(PSTH_SacOff_mean(TimeSequence_SacOff_mean<PlotSacPostBuffer & TimeSequence_SacOff_mean>-PlotSacPreBuffer))*1.05;

shadedErrorBar([-PlotSacPreBuffer:-PlotSacPreBuffer+length(EyeEcc_Aligned_SacOff_mean)-1],EyeEcc_Aligned_SacOff_mean+yrange(2),EyeEcc_Aligned_SacOff_sem,'lineprops',{[0.9290    0.6940    0.1250]},'transparent',1,'patchSaturation',0.3);
ylim([0,yrange(2)+16]);

xlim([-PlotSacPreBuffer, PlotSacPostBuffer]);%Keep in the same scale
yrange=ylim;
plot([0,0],yrange,'--k');

set(findall(gca, 'Type', 'Line'),'LineWidth',3);
set(gca,'FontSize',15,'FontWeight','Bold');
set(gca,'LineWidth',3);

xlabel('Time from Saccade Offset','FontSize',25,'FontWeight','Bold');
ylabel('Spikes/s','FontSize',25,'FontWeight','Bold');







%%%%%%%%%Align on Fractal Onset: Good vs Bad %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(2,3,4);
%%
shadedErrorBar(PSTH_Time_FraOn_All,PSTH_FR_Good_mean,PSTH_FR_Good_sem,'lineprops',{[1    0    0]},'transparent',1,'patchSaturation',0.3);

hold on
shadedErrorBar(PSTH_Time_FraOn_All,PSTH_FR_Bad_mean,PSTH_FR_Bad_sem,'lineprops',{[0    0    1]},'transparent',1,'patchSaturation',0.3);

yrange=ylim;
yrange(2)=max(MeanPSTH_FR_FractalOn_mean(PSTH_Time_FraOn_All<PlotSacPostBuffer & PSTH_Time_FraOn_All>-PlotSacPreBuffer))*1.05;

shadedErrorBar([-PlotStartOffset:-PlotStartOffset+length(Eye_Ecc_Good_mean)-1],Eye_Ecc_Good_mean+yrange(2),Eye_Ecc_Good_sem,'lineprops',{[1    0    0]},'transparent',1,'patchSaturation',0.3);
shadedErrorBar([-PlotStartOffset:-PlotStartOffset+length(Eye_Ecc_Bad_mean)-1],Eye_Ecc_Bad_mean+yrange(2),Eye_Ecc_Bad_sem,'lineprops',{[0    0    1]},'transparent',1,'patchSaturation',0.3);
ylim([0,yrange(2)+16]);

xlim([-PlotSacPreBuffer, PlotSacPostBuffer]);%Keep in the same scale
yrange=ylim;
plot([0,0],yrange,'--k');

set(findall(gca, 'Type', 'Line'),'LineWidth',3);
set(gca,'FontSize',15,'FontWeight','Bold');
set(gca,'LineWidth',3);

xlabel('Time from Fractal Onset(ms)','FontSize',25,'FontWeight','Bold');
ylabel('Spikes/s','FontSize',25,'FontWeight','Bold');



%%%%%%%%%Align on Saccade Onset: Good vs Bad %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(2,3,5);
shadedErrorBar(TimeSequence_SacOn_mean,PSTH_SacOn_Good_mean,PSTH_SacOn_Good_sem,'lineprops',{[1    0    0]},'transparent',1,'patchSaturation',0.3);

hold on
shadedErrorBar(TimeSequence_SacOn_mean,PSTH_SacOn_Bad_mean,PSTH_SacOn_Bad_sem,'lineprops',{[0    0    1]},'transparent',1,'patchSaturation',0.3);

yrange=ylim;
yrange(2)=max(PSTH_SacOn_mean(TimeSequence_SacOn_mean<PlotSacPostBuffer & TimeSequence_SacOn_mean>-PlotSacPreBuffer))*1.05;

shadedErrorBar([-PlotSacPreBuffer:-PlotSacPreBuffer+length(EyeEcc_Aligned_Good_mean)-1],EyeEcc_Aligned_Good_mean+yrange(2),EyeEcc_Aligned_Good_sem,'lineprops',{[1    0    0]},'transparent',1,'patchSaturation',0.3);
shadedErrorBar([-PlotSacPreBuffer:-PlotSacPreBuffer+length(EyeEcc_Aligned_Bad_mean)-1],EyeEcc_Aligned_Bad_mean+yrange(2),EyeEcc_Aligned_Bad_sem,'lineprops',{[0    0    1]},'transparent',1,'patchSaturation',0.3);
ylim([0,yrange(2)+16]);

xlim([-PlotSacPreBuffer, PlotSacPostBuffer]);%Keep in the same scale
yrange=ylim;
plot([0,0],yrange,'--k');

set(findall(gca, 'Type', 'Line'),'LineWidth',3);
set(gca,'FontSize',15,'FontWeight','Bold');
set(gca,'LineWidth',3);

xlabel('Time from Saccade Onset','FontSize',25,'FontWeight','Bold');
ylabel('Spikes/s','FontSize',25,'FontWeight','Bold');

%%%%%%%%%Align on Saccade Offset: Good vs Bad %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(2,3,6);
shadedErrorBar(TimeSequence_SacOff_mean,PSTH_SacOff_Good_mean,PSTH_SacOff_Good_sem,'lineprops',{[1    0    0]},'transparent',1,'patchSaturation',0.3);

hold on
shadedErrorBar(TimeSequence_SacOff_mean,PSTH_SacOff_Bad_mean,PSTH_SacOff_Bad_sem,'lineprops',{[0    0    1]},'transparent',1,'patchSaturation',0.3);

yrange=ylim;
yrange(2)=max(PSTH_SacOff_mean(TimeSequence_SacOff_mean<PlotSacPostBuffer & TimeSequence_SacOff_mean>-PlotSacPreBuffer))*1.05;

shadedErrorBar([-PlotSacPreBuffer:-PlotSacPreBuffer+length(EyeEcc_Aligned_SacOff_Good_mean)-1],EyeEcc_Aligned_SacOff_Good_mean+yrange(2),EyeEcc_Aligned_SacOff_Good_sem,'lineprops',{[1    0    0]},'transparent',1,'patchSaturation',0.3);
shadedErrorBar([-PlotSacPreBuffer:-PlotSacPreBuffer+length(EyeEcc_Aligned_SacOff_Bad_mean)-1],EyeEcc_Aligned_SacOff_Bad_mean+yrange(2),EyeEcc_Aligned_SacOff_Bad_sem,'lineprops',{[0    0    1]},'transparent',1,'patchSaturation',0.3);
ylim([0,yrange(2)+16]);

xlim([-PlotSacPreBuffer, PlotSacPostBuffer]);%Keep in the same scale
yrange=ylim;
plot([0,0],yrange,'--k');

set(findall(gca, 'Type', 'Line'),'LineWidth',3);
set(gca,'FontSize',15,'FontWeight','Bold');
set(gca,'LineWidth',3);

xlabel('Time from Saccade Offset','FontSize',25,'FontWeight','Bold');
ylabel('Spikes/s','FontSize',25,'FontWeight','Bold');


%Add title

title_msg=sprintf(' %s, SpikeChannel: %d, TrialNum: %d ',FileName,SpikeChannel,TrialNum);
a=annotation('textbox',[0.25,0.68,0.3,0.3],'String',title_msg,'FitBoxToText','on');
a.FontSize=25;
a.FontWeight='Bold';
a.LineStyle='none';


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Figure 2: PSTH for different
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%RTs
figtitlestr{2}='PSTHGroupedByRT';

fig{2}=PrepareFigure(111,'w',[100,100, 1200,800],'Name',figtitlestr{2});
title('PSTH grouped by Saccade Reaction Time','FontSize',25,'FontWeight','Bold');
axis off;
%First for PSTH aligned on fractal onset
%Eye Trace
axes('Position',[0.05,0.7,0.25,0.15]);
for i=1:SepGroupRT 
    shadedErrorBar([-PlotStartOffset:-PlotStartOffset+length(EyeEcc_FraOn_RTGroup_mean(i,:))-1],EyeEcc_FraOn_RTGroup_mean(i,:),EyeEcc_FraOn_RTGroup_sem(i,:),'lineprops',{Color{i}},'transparent',1,'patchSaturation',0.3);  
    hold on

end
ylim([0,18]);
yrange=ylim;
h=plot([0,0],yrange,'--k');
xlim([-PlotSacPreBuffer, PlotSacPostBuffer]);%Keep in the same scale
set(findall(gca, 'Type', 'Line'),'LineWidth',3);
xticks([]);
xticklabels([]);
set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
ylabel('Eye Ecc(degree)','FontSize',20,'FontWeight','Bold');

%Neural firing
 axes('Position',[0.05,0.15,0.25,0.5]);
 for i=1:SepGroupRT   
    shadedErrorBar(PSTH_Time_FraOn_All,PSTH_FraOn_RTGroup_mean(i,:),PSTH_FraOn_RTGroup_sem(i,:),'lineprops',{Color{i}},'transparent',1,'patchSaturation',0.3);
    hold on
 end 
 yrange=ylim;
 h=plot([0,0],yrange,'--k');
  xlim([-PlotSacPreBuffer, PlotSacPostBuffer]);%Keep in the same scale
  yrange(1)=min(min(PSTH_FraOn_RTGroup_mean(1:end-1,PSTH_Time_FraOn_All>-PlotSacPreBuffer &PSTH_Time_FraOn_All<PlotSacPostBuffer )))*0.8;
  yrange(2)=max(max(PSTH_FraOn_RTGroup_mean(1:end-1,PSTH_Time_FraOn_All>-PlotSacPreBuffer &PSTH_Time_FraOn_All<PlotSacPostBuffer )))*1.2;
  ylim(yrange);
  
  set(findall(gca, 'Type', 'Line'),'LineWidth',3);   
  set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
  legend(GroupStr,'Location','northwest');
ylabel('Firing Rate (Spikes/s)','FontSize',20,'FontWeight','Bold');
xlabel('Time from Fractal Onset(ms)','FontSize',20,'FontWeight','Bold');



set(get(get(h,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');


%Second for PSTH aligned on saccade onset
%Eye Trace
axes('Position',[0.35,0.7,0.25,0.15]);
for i=1:SepGroupRT 
    shadedErrorBar([-PlotSacPreBuffer:EyeBinWidth: PlotSacPostBuffer],EyeEcc_SacOn_RTGroup_mean(i,:),EyeEcc_SacOn_RTGroup_sem(i,:),'lineprops',{Color{i}},'transparent',1,'patchSaturation',0.3);  
    hold on

end
ylim([0,18]);
yrange=ylim;
h=plot([0,0],yrange,'--k');
xlim([-PlotSacPreBuffer, PlotSacPostBuffer]);%Keep in the same scale
set(findall(gca, 'Type', 'Line'),'LineWidth',3);
xticks([]);
xticklabels([]);
set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
ylabel('Eye Ecc(degree)','FontSize',20,'FontWeight','Bold');

%Neural firing
 axes('Position',[0.35,0.15,0.25,0.5]);
 for i=1:SepGroupRT   
    shadedErrorBar(TimeSequence_SacOn_mean,PSTH_SacOn_RTGroup_mean(i,:),PSTH_SacOn_RTGroup_sem(i,:),'lineprops',{Color{i}},'transparent',1,'patchSaturation',0.3);
    hold on
 end 
 yrange=ylim;
 h=plot([0,0],yrange,'--k');
  xlim([-PlotSacPreBuffer, PlotSacPostBuffer]);%Keep in the same scale
  set(findall(gca, 'Type', 'Line'),'LineWidth',3);   
  set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
  legend(GroupStr,'Location','northwest');
ylabel('Firing Rate (Spikes/s)','FontSize',20,'FontWeight','Bold');
xlabel('Time from Saccade Onset(ms)','FontSize',20,'FontWeight','Bold');



set(get(get(h,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');

%Third, for PSTH aligned on saccade offset
%Eye Trace
axes('Position',[0.65,0.7,0.25,0.15]);
for i=1:SepGroupRT 
    shadedErrorBar([-PlotSacPreBuffer:EyeBinWidth: PlotSacPostBuffer],EyeEcc_SacOff_RTGroup_mean(i,:),EyeEcc_SacOff_RTGroup_sem(i,:),'lineprops',{Color{i}},'transparent',1,'patchSaturation',0.3);  
    hold on

end
ylim([0,18]);
yrange=ylim;
h=plot([0,0],yrange,'--k');
xlim([-PlotSacPreBuffer, PlotSacPostBuffer]);%Keep in the same scale
set(findall(gca, 'Type', 'Line'),'LineWidth',3);
xticks([]);
xticklabels([]);
set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
ylabel('Eye Ecc(degree)','FontSize',20,'FontWeight','Bold');

%Neural firing
 axes('Position',[0.65,0.15,0.25,0.5]);
 for i=1:SepGroupRT   
    shadedErrorBar(TimeSequence_SacOff_mean,PSTH_SacOff_RTGroup_mean(i,:),PSTH_SacOff_RTGroup_sem(i,:),'lineprops',{Color{i}},'transparent',1,'patchSaturation',0.3);
    hold on
 end 
 yrange=ylim;
 h=plot([0,0],yrange,'--k');
  xlim([-PlotSacPreBuffer, PlotSacPostBuffer]);%Keep in the same scale
  set(findall(gca, 'Type', 'Line'),'LineWidth',3);   
  set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
  legend(GroupStr,'Location','northwest');
ylabel('Firing Rate (Spikes/s)','FontSize',20,'FontWeight','Bold');
xlabel('Time from Saccade Offset(ms)','FontSize',20,'FontWeight','Bold');



set(get(get(h,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Figure 3: PSTH for different
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Saccade Amplitudes
figtitlestr{3}='PSTHGroupedByAmp';

fig{3}=PrepareFigure(112,'w',[150,100, 1200,800],'Name',figtitlestr{3});
title('PSTH grouped by Saccade Amplitude','FontSize',25,'FontWeight','Bold');
axis off;
%First for PSTH aligned on fractal onset
%Eye Trace
axes('Position',[0.05,0.7,0.25,0.15]);
for i=1:SepGroupAmpNum
    shadedErrorBar([-PlotStartOffset:-PlotStartOffset+length(EyeEcc_FracOn_AmpGroup_mean(i,:))-1],EyeEcc_FracOn_AmpGroup_mean(i,:),EyeEcc_FracOn_AmpGroup_sem(i,:),'lineprops',{Color{i}},'transparent',1,'patchSaturation',0.3);  
    hold on

end
ylim([0,18]);
yrange=ylim;

  
h=plot([0,0],yrange,'--k');
xlim([-PlotSacPreBuffer, PlotSacPostBuffer]);%Keep in the same scale
set(findall(gca, 'Type', 'Line'),'LineWidth',3);
xticks([]);
xticklabels([]);
set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
ylabel('Eye Ecc(degree)','FontSize',20,'FontWeight','Bold');

%Neural firing
axes('Position',[0.05,0.15,0.25,0.5]);
 for i=1:SepGroupAmpNum 
    shadedErrorBar(PSTH_Time_FraOn_All,PSTH_FraOn_AmpGroup_mean(i,:),PSTH_FraOn_AmpGroup_sem(i,:),'lineprops',{Color{i}},'transparent',1,'patchSaturation',0.3);
    hold on
 end 
 yrange=ylim;
 yrange(1)=min(min(PSTH_FraOn_RTGroup_mean(1:end-1,PSTH_Time_FraOn_All>-PlotSacPreBuffer &PSTH_Time_FraOn_All<PlotSacPostBuffer )))*0.8;
yrange(2)=max(max(PSTH_FraOn_RTGroup_mean(1:end-1,PSTH_Time_FraOn_All>-PlotSacPreBuffer &PSTH_Time_FraOn_All<PlotSacPostBuffer )))*1.2;
 ylim(yrange);
 
 h=plot([0,0],yrange,'--k');
  xlim([-PlotSacPreBuffer, PlotSacPostBuffer]);%Keep in the same scale
  set(findall(gca, 'Type', 'Line'),'LineWidth',3);   
  set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
  legend(GroupStr,'Location','northwest');
ylabel('Firing Rate (Spikes/s)','FontSize',20,'FontWeight','Bold');
xlabel('Time from Fractal Onset(ms)','FontSize',20,'FontWeight','Bold');



set(get(get(h,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');


%Second for PSTH aligned on saccade onset
%Eye Trace
axes('Position',[0.35,0.7,0.25,0.15]);
for i=1:SepGroupAmpNum
    shadedErrorBar([-PlotSacPreBuffer:EyeBinWidth: PlotSacPostBuffer],EyeEcc_SacOn_AmpGroup_mean(i,:),EyeEcc_SacOn_AmpGroup_sem(i,:),'lineprops',{Color{i}},'transparent',1,'patchSaturation',0.3);  
    hold on

end
ylim([0,18]);
yrange=ylim;
h=plot([0,0],yrange,'--k');
xlim([-PlotSacPreBuffer, PlotSacPostBuffer]);%Keep in the same scale
set(findall(gca, 'Type', 'Line'),'LineWidth',3);
xticks([]);
xticklabels([]);
set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
ylabel('Eye Ecc(degree)','FontSize',20,'FontWeight','Bold');

%Neural firing
 axes('Position',[0.35,0.15,0.25,0.5]);
 for i=1:SepGroupAmpNum
    shadedErrorBar(TimeSequence_SacOn_mean,PSTH_SacOn_AmpGroup_mean(i,:),PSTH_SacOn_AmpGroup_sem(i,:),'lineprops',{Color{i}},'transparent',1,'patchSaturation',0.3);
    hold on
 end 

 
 yrange=ylim;
 h=plot([0,0],yrange,'--k');
  xlim([-PlotSacPreBuffer, PlotSacPostBuffer]);%Keep in the same scale
  set(findall(gca, 'Type', 'Line'),'LineWidth',3);   
  set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
  legend(GroupStr,'Location','northwest');
ylabel('Firing Rate (Spikes/s)','FontSize',20,'FontWeight','Bold');
xlabel('Time from Saccade Onset(ms)','FontSize',20,'FontWeight','Bold');



set(get(get(h,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');


%Third, for PSTH aligned on saccade offset
%Eye Trace
axes('Position',[0.65,0.7,0.25,0.15]);
for i=1:SepGroupAmpNum
    shadedErrorBar([-PlotSacPreBuffer:EyeBinWidth: PlotSacPostBuffer],EyeEcc_SacOff_AmpGroup_mean(i,:),EyeEcc_SacOff_AmpGroup_sem(i,:),'lineprops',{Color{i}},'transparent',1,'patchSaturation',0.3);  
    hold on

end
ylim([0,18]);
yrange=ylim;
h=plot([0,0],yrange,'--k');
xlim([-PlotSacPreBuffer, PlotSacPostBuffer]);%Keep in the same scale
set(findall(gca, 'Type', 'Line'),'LineWidth',3);
xticks([]);
xticklabels([]);
set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
ylabel('Eye Ecc(degree)','FontSize',20,'FontWeight','Bold');

%Neural firing
 axes('Position',[0.65,0.15,0.25,0.5]);
 for i=1:SepGroupAmpNum  
    shadedErrorBar(TimeSequence_SacOff_mean,PSTH_SacOff_AmpGroup_mean(i,:),PSTH_SacOff_AmpGroup_sem(i,:),'lineprops',{Color{i}},'transparent',1,'patchSaturation',0.3);
    hold on
 end 
 yrange=ylim;
 h=plot([0,0],yrange,'--k');
  xlim([-PlotSacPreBuffer, PlotSacPostBuffer]);%Keep in the same scale
  set(findall(gca, 'Type', 'Line'),'LineWidth',3);   
  set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
  legend(GroupStr,'Location','northwest');
ylabel('Firing Rate (Spikes/s)','FontSize',20,'FontWeight','Bold');
xlabel('Time from Saccade Offset(ms)','FontSize',20,'FontWeight','Bold');



set(get(get(h,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');

%%
%%%%%%%%%%%%%Figure 3: Firing Rate as a function of spike amplitude
figtitlestr{4}='RFwithSacAmp';

fig{4}=PrepareFigure(113,'w',[175,200, 1200,600],'Name',figtitlestr{4});
subplot(2,3,1)
s_g=scatter(SaccadeAmplitude_First(GoodObjectIndex),SpikeRate_PreSac(GoodObjectIndex),40,'or','filled','MarkerEdgeColor','w');
hold on
s_b=scatter(SaccadeAmplitude_First(~GoodObjectIndex),SpikeRate_PreSac(~GoodObjectIndex),40,'ob','filled','MarkerEdgeColor','w');


%Add fitting lines
xrange=xlim;
x_temp=xrange(1):0.1:xrange(2);
y_line=LinepreSac(1)+LinepreSac(2)*x_temp;
plot(x_temp,y_line,'-','Color','#969391','LineWidth',3);
legend('Good','Bad','Fitting','Location','NorthWest');

set(gca,'FontSize',20,'FontWeight','Bold','LineWidth',3);
xlabel('Saccade Amplitude(deg)','FontSize',25,'FontWeight','Bold');
ylabel('PreSacFR(Hz)','FontSize',25,'FontWeight','Bold');

subplot(2,3,2)
scatter(SaccadeAmplitude_First(GoodObjectIndex),SpikeRate_DurSac(GoodObjectIndex),40,'or','filled','MarkerEdgeColor','w');
hold on
scatter(SaccadeAmplitude_First(~GoodObjectIndex),SpikeRate_DurSac(~GoodObjectIndex),40,'ob','filled','MarkerEdgeColor','w');
%Add fitting lines
xrange=xlim;
x_temp=xrange(1):0.1:xrange(2);
y_line=LinedurSac(1)+LinedurSac(2)*x_temp;
plot(x_temp,y_line,'-','Color','#969391','LineWidth',3);

set(gca,'FontSize',20,'FontWeight','Bold','LineWidth',3);
xlabel('Saccade Amplitude(deg)','FontSize',25,'FontWeight','Bold');
ylabel('DurSacFR(Hz)','FontSize',25,'FontWeight','Bold');

subplot(2,3,3)
scatter(SaccadeAmplitude_First(GoodObjectIndex),SpikeRate_PostSac(GoodObjectIndex),40,'or','filled','MarkerEdgeColor','w');
hold on
scatter(SaccadeAmplitude_First(~GoodObjectIndex),SpikeRate_PostSac(~GoodObjectIndex),40,'ob','filled','MarkerEdgeColor','w');

%Add fitting lines
xrange=xlim;
x_temp=xrange(1):0.1:xrange(2);
y_line=LinepostSac(1)+LinepostSac(2)*x_temp;
plot(x_temp,y_line,'-','Color','#969391','LineWidth',3);

set(gca,'FontSize',20,'FontWeight','Bold','LineWidth',3);
xlabel('Saccade Amplitude(deg)','FontSize',25,'FontWeight','Bold');
ylabel('PostSacFR(Hz)','FontSize',25,'FontWeight','Bold');


%Bar plot for Pearson's Correlation Coefficient (PreSaccade)
subplot(2,3,4);
g=bar([r_preSac,r_preSac_Good,r_preSac_Bad]);
g.FaceColor='Flat';
g.CData=[[150,147,145]/255;
          [1,0,0];
          [0,0,1]];
%Add significance level
XData=g.XData;
yData=g.YData;

pCorrcoef=[p_preSac,p_preSac_Good,p_preSac_Bad];
yrange=ylim;

text(XData(pCorrcoef<0.05)-0.045,yData(pCorrcoef<0.05)+yrange(2)*0.1,'*','FontSize',40,'Color','r')
text(XData(pCorrcoef<0.01)-0.045,yData(pCorrcoef<0.01)+yrange(2)*0.2,'*','FontSize',40,'Color','r')
text(XData(pCorrcoef<0.001)-0.045,yData(pCorrcoef<0.001)+yrange(2)*0.3,'*','FontSize',40,'Color','r')
ylim([yrange(1),yrange(2)*1.3]);
box off


xticklabels({'All','Good','Bad'});
set(gca,'FontSize',20,'FontWeight','Bold','LineWidth',3);
ylabel(['Peason''s Correlation',newline,'Coefficient'],'FontSize',25,'FontWeight','Bold');
xlabel('PreSac','FontSize',25,'FontWeight','Bold');

%Bar plot for Pearson's Correlation Coefficient (During Saccade)
subplot(2,3,5);
g=bar([r_durSac,r_durSac_Good,r_durSac_Bad]);
g.FaceColor='Flat';
g.CData=[[150,147,145]/255;
          [1,0,0];
          [0,0,1]];
%Add significance level
XData=g.XData;
yData=g.YData;

pCorrcoef=[p_durSac,p_durSac_Good,p_durSac_Bad];
yrange=ylim;

text(XData(pCorrcoef<0.05)-0.045,yData(pCorrcoef<0.05)+yrange(2)*0.1,'*','FontSize',40,'Color','r')
text(XData(pCorrcoef<0.01)-0.045,yData(pCorrcoef<0.01)+yrange(2)*0.2,'*','FontSize',40,'Color','r')
text(XData(pCorrcoef<0.001)-0.045,yData(pCorrcoef<0.001)+yrange(2)*0.3,'*','FontSize',40,'Color','r')
ylim([yrange(1),yrange(2)*1.3]);
box off


xticklabels({'All','Good','Bad'});
set(gca,'FontSize',20,'FontWeight','Bold','LineWidth',3);
%ylabel(['Peason''s Correlation',newline,'Coefficient'],'FontSize',25,'FontWeight','Bold');
xlabel('DurSac','FontSize',25,'FontWeight','Bold');


%Bar plot for Pearson's Correlation Coefficient (During Saccade)
subplot(2,3,6);
g=bar([r_postSac,r_postSac_Good,r_postSac_Bad]);
g.FaceColor='Flat';
g.CData=[[150,147,145]/255;
          [1,0,0];
          [0,0,1]];
%Add significance level
XData=g.XData;
yData=g.YData;

pCorrcoef=[p_postSac,p_postSac_Good,p_postSac_Bad];
yrange=ylim;

text(XData(pCorrcoef<0.05)-0.045,yData(pCorrcoef<0.05)+yrange(2)*0.1,'*','FontSize',40,'Color','r')
text(XData(pCorrcoef<0.01)-0.045,yData(pCorrcoef<0.01)+yrange(2)*0.2,'*','FontSize',40,'Color','r')
text(XData(pCorrcoef<0.001)-0.045,yData(pCorrcoef<0.001)+yrange(2)*0.3,'*','FontSize',40,'Color','r')
ylim([yrange(1),yrange(2)*1.3]);
box off


xticklabels({'All','Good','Bad'});
set(gca,'FontSize',20,'FontWeight','Bold','LineWidth',3);
%ylabel(['Peason''s Correlation',newline,'Coefficient'],'FontSize',25,'FontWeight','Bold');
xlabel('PostSac','FontSize',25,'FontWeight','Bold');


title_msg={'Correlation Coefficient on case-by-case basis (Pre,Dur,Post Sac)'};
a=annotation('textbox',[0.2,0.8,0.3,0.2],'String',title_msg,'FitBoxToText','on');
a.FontSize=25;
a.FontWeight='Bold';
a.LineStyle='none';




%}
if ExamplePlot==1
%% Example Eye trace aligned on different events
figtitlestr{5}='Example Eye movements';

fig{5}=PrepareFigure(114,'w',[185,200, 1200,400],'Name',figtitlestr{5});

ExampleTrial=[1:20];

%Align on target onset
subplot(2,2,1);
X=repmat([-PlotStartOffset:-PlotStartOffset+length(Eye_Eccentricity_FractalOn(1,:))-1],length(ExampleTrial),1);
plot(X',Eye_Eccentricity_FractalOn(ExampleTrial,:)','-k')

yticks([0,15]);
box off
set(gca,'LineWidth',3,'FontSize',25,'FontWeight','Bold');
xlabel('Time from fractal onset (ms)');
xlim([-100,500])

%Align on Saccade onset
subplot(2,2,2);

X2=repmat([-PlotStartOffset:-PlotStartOffset+length(EyeEcc_Aligned(1,:))-1],length(ExampleTrial),1);
plot(X2',EyeEcc_Aligned(ExampleTrial,:)','-k')

yticks([0,15]);
box off
set(gca,'LineWidth',3,'FontSize',25,'FontWeight','Bold');
xlabel('Time from Saccade onset (ms)');

%Align on Saccade offset
subplot(2,2,3);

X2=repmat([-PlotStartOffset:-PlotStartOffset+length(EyeEcc_Aligned_SacOff(1,:))-1],length(ExampleTrial),1);
plot(X2',EyeEcc_Aligned_SacOff(ExampleTrial,:)','-k')

yticks([0,15]);
box off
set(gca,'LineWidth',3,'FontSize',25,'FontWeight','Bold');
xlabel('Time from Saccade offset (ms)');

%Align on Fractal offset

subplot(2,2,4);
X=repmat([-PlotStartOffset:-PlotStartOffset+length(EyeEcc_Aligned_FraOff(1,:))-1],length(ExampleTrial),1);
plot(X',EyeEcc_Aligned_FraOff(ExampleTrial,:)','-k')

yticks([0,15]);
box off
set(gca,'LineWidth',3,'FontSize',25,'FontWeight','Bold');
xlabel('Time from fractal offset (ms)');


%%%%%%%%%%%%%%% Group separation  %%%%%%%%%%%%%%%%%%%%%


figtitlestr{6}='RTandAmplitudeSeparations';

fig{6}=PrepareFigure(115,'w',[195,200, 1200,400],'Name',figtitlestr{6});
%Saccade Latency
subplot(1,3,1);
b=bar([1:length(RT_Group_mean)],[RT_Group_mean]);
b.FaceColor='flat';
b.EdgeColor='flat';
b.CData=Color_RGB(1:length(RT_Group_mean),:);
hold on
errorbar([1:length(RT_Group_mean)],[RT_Group_mean],[RT_Group_sem],'.k','LineWidth',3);
set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
xticks([1:length(RT_Group_mean)]);
xticklabels(GroupStr);
ylabel('Saccade Latency (ms)','FontSize',25);
box off

%Saccade Amplitude
subplot(1,3,2);
b=bar([1:length(Amp_Group_mean)],[Amp_Group_mean]);
b.FaceColor='flat';
b.EdgeColor='flat';
b.CData=Color_RGB(1:length(Amp_Group_mean),:);
hold on
errorbar([1:length(Amp_Group_mean)],[Amp_Group_mean],[Amp_Group_sem],'.k','LineWidth',3);
set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
xticks([1:length(Amp_Group_mean)]);
xticklabels(GroupStr);
ylabel('Saccade Amplitude (deg)','FontSize',25);
box off

end %%End of the plot example




%% Amplitude group stats
figtitlestr{7}='Group by amplitude';

fig{7}=PrepareFigure(115,'w',[195,200, 1200,400],'Name',figtitlestr{7});

SpikeRate_PreSac_Amp_mean=SpikeRate_PreSac_Amp_mean(1:SepGroupRT);
SpikeRate_PreSac_Amp_sem=SpikeRate_PreSac_Amp_sem(1:SepGroupRT);

SpikeRate_DurSac_Amp_mean=SpikeRate_DurSac_Amp_mean(1:SepGroupRT);
SpikeRate_DurSac_Amp_sem=SpikeRate_DurSac_Amp_sem(1:SepGroupRT);

SpikeRate_PostSac_Amp_mean=SpikeRate_PostSac_Amp_mean(1:SepGroupRT);
SpikeRate_PostSac_Amp_sem=SpikeRate_PostSac_Amp_sem(1:SepGroupRT);


subplot(1,3,1)


b=bar([1:length(SpikeRate_PreSac_Amp_mean)],SpikeRate_PreSac_Amp_mean);
b.FaceColor='flat';
b.EdgeColor='flat';
b.CData=Color_RGB(1:length(SpikeRate_PreSac_Amp_mean),:);
hold on
errorbar([1:length(SpikeRate_PreSac_Amp_mean)],SpikeRate_PreSac_Amp_mean,SpikeRate_PreSac_Amp_sem,'.k','LineWidth',3);
set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
xticks([1:length(SpikeRate_PreSac_Amp_mean)]);
xticklabels(GroupStr);
ylabel('Firing Rate (Hz)','FontSize',25);
box off

%Add significance level
XData=b.XData;
yData=b.YData;

yrange=ylim;


if ANOVA_p_Amp_pre<0.05 & ANOVA_p_Amp_pre>=0.01
    text(XData(2),yrange(2)*1.1,'*','FontSize',40,'Color','r');
end
    
if ANOVA_p_Amp_pre<0.01 & ANOVA_p_Amp_pre>=0.001
     text(XData(2)-0.3,yrange(2)*1.1,'*','FontSize',40,'Color','r');
     text(XData(2)+0.3,yrange(2)*1.1,'*','FontSize',40,'Color','r');
     
    
end
if ANOVA_p_Amp_pre<0.001 
     text(XData(2)-0.3,yrange(2)*1.1,'*','FontSize',40,'Color','r');
     text(XData(2),yrange(2)*1.1,'*','FontSize',40,'Color','r');
     text(XData(2)+0.3,yrange(2)*1.1,'*','FontSize',40,'Color','r');
    
end
       
        
        
   

ylim([yrange(1),yrange(2)*1.3]);
box off
title('Pre-Sac','FontSize',25,'FontWeight','Bold');



%During Saccade
subplot(1,3,2)


b=bar([1:length(SpikeRate_DurSac_Amp_mean)],SpikeRate_DurSac_Amp_mean);
b.FaceColor='flat';
b.EdgeColor='flat';
b.CData=Color_RGB(1:length(SpikeRate_DurSac_Amp_mean),:);
hold on
errorbar([1:length(SpikeRate_DurSac_Amp_mean)],SpikeRate_DurSac_Amp_mean,SpikeRate_DurSac_Amp_sem,'.k','LineWidth',3);
set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
xticks([1:length(SpikeRate_DurSac_Amp_mean)]);
xticklabels(GroupStr);
ylabel('Firing Rate (Hz)','FontSize',25);
box off

%Add significance level
XData=b.XData;
yData=b.YData;

yrange=ylim;

if ANOVA_p_Amp_dur<0.05 & ANOVA_p_Amp_dur>=0.01
    text(XData(2),yrange(2)*1.1,'*','FontSize',40,'Color','r');
end
    
if ANOVA_p_Amp_dur<0.01 & ANOVA_p_Amp_dur>=0.001
     text(XData(2)-0.3,yrange(2)*1.1,'*','FontSize',40,'Color','r');
     text(XData(2)+0.3,yrange(2)*1.1,'*','FontSize',40,'Color','r');
     
    
end
if ANOVA_p_Amp_dur<0.001 
     text(XData(2)-0.3,yrange(2)*1.1,'*','FontSize',40,'Color','r');
     text(XData(2),yrange(2)*1.1,'*','FontSize',40,'Color','r');
     text(XData(2)+0.3,yrange(2)*1.1,'*','FontSize',40,'Color','r');
    
end

ylim([yrange(1),yrange(2)*1.3]);
box off
title('Dur-Sac','FontSize',25,'FontWeight','Bold');


subplot(1,3,3)


b=bar([1:length(SpikeRate_PostSac_Amp_mean)],SpikeRate_PostSac_Amp_mean);
b.FaceColor='flat';
b.EdgeColor='flat';
b.CData=Color_RGB(1:length(SpikeRate_PostSac_Amp_mean),:);
hold on
errorbar([1:length(SpikeRate_PostSac_Amp_mean)],SpikeRate_PostSac_Amp_mean,SpikeRate_PostSac_Amp_sem,'.k','LineWidth',3);
set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
xticks([1:length(SpikeRate_PostSac_Amp_mean)]);
xticklabels(GroupStr);
ylabel('Firing Rate (Hz)','FontSize',25);
box off

%Add significance level
XData=b.XData;
yData=b.YData;

yrange=ylim;

if ANOVA_p_Amp_post<0.05 & ANOVA_p_Amp_post>=0.01
    text(XData(2),yrange(2)*1.1,'*','FontSize',40,'Color','r');
end
    
if ANOVA_p_Amp_post<0.01 & ANOVA_p_Amp_post>=0.001
     text(XData(2)-0.3,yrange(2)*1.1,'*','FontSize',40,'Color','r');
     text(XData(2)+0.3,yrange(2)*1.1,'*','FontSize',40,'Color','r');
     
    
end
if ANOVA_p_Amp_post<0.001 
     text(XData(2)-0.3,yrange(2)*1.1,'*','FontSize',40,'Color','r');
     text(XData(2),yrange(2)*1.1,'*','FontSize',40,'Color','r');
     text(XData(2)+0.3,yrange(2)*1.1,'*','FontSize',40,'Color','r');
    
end

ylim([yrange(1),yrange(2)*1.3]);
box off
title('Post-Sac','FontSize',25,'FontWeight','Bold');

if ExamplePlot==1
%%
%Polar Plot for saccade angle and amplitude

figtitlestr{8}='Example Eye movements';

fig{8}=PrepareFigure(116,'w',[200,200, 600,400],'Name',figtitlestr{8});

p= polarplot((SaccadeAngle_First)*pi/180',SaccadeAmplitude_First,'o');
p.MarkerFaceColor='#EC8A26';
p.MarkerEdgeColor='w';
p.MarkerSize=10;
set(gca,'ThetaDir' , 'counterclockwise');
set(gca,'ThetaZeroLocation','top')
set(gca,'FontSize',25,'FontWeight','Bold','LineWidth',3);

thetaticks([0:45:359]);
thetaticklabels([0:45:359]);
title('Saccade Amplitude c.s. Saccade Angle','FontSize',25,'FontWeight','Bold');
end



%%
%Polar Plot for saccade vector tuning

figtitlestr{9}='SaccadeVectorTuning_Average';

fig{9}=PrepareFigure(117,'w',[220,200, 1200,400],'Name',figtitlestr{9});


subtitlestr={'Pre-Sac Tuning','Dur-Sac Tuning','Post-Sac Tuning'};

Selection_temp=~isnan(SaccadeVector);
SaccadeVector=SaccadeVector(Selection_temp);
FiringRate_PreDurPost_mean=FiringRate_PreDurPost_mean(Selection_temp,:);
FiringRate_PreDurPost_sem=FiringRate_PreDurPost_sem(Selection_temp,:);


for j=1:3
%Pre-saccade , During, Post saccade tuning
subplot(1,3,j);
%{
p= polarplot((SaccadeAngle_First)*pi/180',SpikeRate_PreSac,'or');
p.MarkerFaceColor='r';
p.MarkerEdgeColor='r';
%}

%hold on
polarwitherrorbar([SaccadeVector/180*pi SaccadeVector(1)/180*pi],[FiringRate_PreDurPost_mean(:,j)' FiringRate_PreDurPost_mean(1,j)],[FiringRate_PreDurPost_sem(:,j)' FiringRate_PreDurPost_sem(1,j)],...
    '-b',3);
hold on
h_p = polarplot([PrefVector(j) PrefVector(j)]/180*pi,[0 PrefVector_Amp(j)],'r');
    set(h_p,'LineWidth',2.5);

set(gca,'ThetaDir' , 'counterclockwise');
set(gca,'ThetaZeroLocation','top')
set(gca,'FontSize',25,'FontWeight','Bold','LineWidth',3);

thetaticks(sort(SaccadeVector));
thetaticklabels(sort(SaccadeVector));
title(subtitlestr{j},'FontSize',25,'FontWeight','Bold');
pstr=sprintf('p= %2.3f',SacVecTuning_p(j));
a=annotation('textbox',[0.15+(j-1)/3.5,0.06,0.1,0.1] ,'String',pstr,'FitBoxToText','on','LineStyle','none');
a.FontSize=30;
a.FontWeight='Bold';
if SacVecTuning_p(j)<0.05
    a.Color='r';
end

rlim([0,max(max(FiringRate_PreDurPost_mean))*1.2]);

 
 

end




if AmpMatchFlag

%Amplitude Match plot

figtitlestr{10}='Amplitude Match';

fig{10}=PrepareFigure(120,'w',[260,200, 1200,400],'Name',figtitlestr{10});


%Pre-saccade , During, Post saccade responses with amplitude match
subplot(1,3,1);

b=bar([1:length(SpikeRate_PreSac_Match_mean)],SpikeRate_PreSac_Match_mean);
b.FaceColor='flat';
b.EdgeColor='flat';
b.CData=[1,0,0;0,0,1];
hold on
errorbar([1:length(SpikeRate_PreSac_Match_mean)],SpikeRate_PreSac_Match_mean,SpikeRate_PreSac_Match_sem,'.k','LineWidth',3);
set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
xticks([1:length(SpikeRate_PreSac_Match_mean)]);
xticklabels({'Good','Bad'});

%Add significant level between good and bad

yrange=ylim;
if p_Pre_AmpMatch<0.05
ctr2 = bsxfun(@plus, b.XData, [b.XOffset]');
hold on
plot(ctr2(1:2), [1 1]*yrange(2)*1.35, '-k', 'LineWidth',3,'Color','r')
if p_Pre_AmpMatch<0.05 & p_Pre_AmpMatch>=0.01
   plot(mean(ctr2(1:2)), yrange(2)*1.4, '*k','MarkerSize',15,'Color','r','LineWidth',3);
elseif p_Pre_AmpMatch<0.01 & p_Pre_AmpMatch>=0.001
   plot(mean(ctr2(1:2))*0.9, yrange(2)*1.4, '*k','MarkerSize',15,'Color','r','LineWidth',3);
   plot(mean(ctr2(1:2))*1.1, yrange(2)*1.4, '*k','MarkerSize',15,'Color','r','LineWidth',3);
elseif p_Pre_AmpMatch<0.001 
   plot(mean(ctr2(1:2))*0.9, yrange(2)*1.4, '*k','MarkerSize',15,'Color','r','LineWidth',3);
   plot(mean(ctr2(1:2))*1.1, yrange(2)*1.4, '*k','MarkerSize',15,'Color','r','LineWidth',3);
   plot(mean(ctr2(1:2))*1.3, yrange(2)*1.4, '*k','MarkerSize',15,'Color','r','LineWidth',3);
end
end


ylabel('Firing Rate (Hz)','FontSize',25);
box off
title('Pre-Sac','FontSize',25);

subplot(1,3,2);

b=bar([1:length(SpikeRate_DurSac_Match_mean)],SpikeRate_DurSac_Match_mean);
b.FaceColor='flat';
b.EdgeColor='flat';
b.CData=[1,0,0;0,0,1];
hold on
errorbar([1:length(SpikeRate_DurSac_Match_mean)],SpikeRate_DurSac_Match_mean,SpikeRate_DurSac_Match_sem,'.k','LineWidth',3);
set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
xticks([1:length(SpikeRate_DurSac_Match_mean)]);
xticklabels({'Good','Bad'});

%Add significant level between good and bad

yrange=ylim;
if p_Dur_AmpMatch<0.05
ctr2 = bsxfun(@plus, b.XData, [b.XOffset]');
hold on
plot(ctr2(1:2), [1 1]*yrange(2)*1.35, '-k', 'LineWidth',3,'Color','r')
if p_Dur_AmpMatch<0.05 & p_Dur_AmpMatch>=0.01
   plot(mean(ctr2(1:2)), yrange(2)*1.4, '*k','MarkerSize',15,'Color','r','LineWidth',3);
elseif p_Dur_AmpMatch<0.01 & p_Dur_AmpMatch>=0.001
   plot(mean(ctr2(1:2))*0.9, yrange(2)*1.4, '*k','MarkerSize',15,'Color','r','LineWidth',3);
   plot(mean(ctr2(1:2))*1.1, yrange(2)*1.4, '*k','MarkerSize',15,'Color','r','LineWidth',3);
elseif p_Dur_AmpMatch<0.001 
   plot(mean(ctr2(1:2))*0.9, yrange(2)*1.4, '*k','MarkerSize',15,'Color','r','LineWidth',3);
   plot(mean(ctr2(1:2))*1.1, yrange(2)*1.4, '*k','MarkerSize',15,'Color','r','LineWidth',3);
   plot(mean(ctr2(1:2))*1.3, yrange(2)*1.4, '*k','MarkerSize',15,'Color','r','LineWidth',3);
end
end

ylabel('Firing Rate (Hz)','FontSize',25);
box off
title('Dur-Sac','FontSize',25);

subplot(1,3,3);

b=bar([1:length(SpikeRate_PostSac_Match_mean)],SpikeRate_PostSac_Match_mean);
b.FaceColor='flat';
b.EdgeColor='flat';
b.CData=[1,0,0;0,0,1];
hold on
errorbar([1:length(SpikeRate_PostSac_Match_mean)],SpikeRate_PostSac_Match_mean,SpikeRate_PostSac_Match_sem,'.k','LineWidth',3);
set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
xticks([1:length(SpikeRate_PostSac_Match_mean)]);
xticklabels({'Good','Bad'});

%Add significant level between good and bad

yrange=ylim;
if p_Post_AmpMatch<0.05
ctr2 = bsxfun(@plus, b.XData, [b.XOffset]');
hold on
plot(ctr2(1:2), [1 1]*yrange(2)*1.35, '-k', 'LineWidth',3,'Color','r')
if p_Post_AmpMatch<0.05 & p_Post_AmpMatch>=0.01
   plot(mean(ctr2(1:2)), yrange(2)*1.4, '*k','MarkerSize',15,'Color','r','LineWidth',3);
elseif p_Post_AmpMatch<0.01 & p_Post_AmpMatch>=0.001
   plot(mean(ctr2(1:2))*0.9, yrange(2)*1.4, '*k','MarkerSize',15,'Color','r','LineWidth',3);
   plot(mean(ctr2(1:2))*1.1, yrange(2)*1.4, '*k','MarkerSize',15,'Color','r','LineWidth',3);
elseif p_Post_AmpMatch<0.001 
   plot(mean(ctr2(1:2))*0.9, yrange(2)*1.4, '*k','MarkerSize',15,'Color','r','LineWidth',3);
   plot(mean(ctr2(1:2))*1.1, yrange(2)*1.4, '*k','MarkerSize',15,'Color','r','LineWidth',3);
   plot(mean(ctr2(1:2))*1.3, yrange(2)*1.4, '*k','MarkerSize',15,'Color','r','LineWidth',3);
end
end

ylabel('Firing Rate (Hz)','FontSize',25);
box off
title('Post-Sac','FontSize',25);

end %End of the if amplitude match

%PSTH aligned on Returned Saccade
figtitlestr{11}='ReturnSaccadePSTH';

fig{11}=PrepareFigure(121,'w',[280,200, 1200,400],'Name',figtitlestr{11});

%%%%%%%%%Align on Returned Saccade Onset%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(1,2,1);

shadedErrorBar(TimeSequence_SacOnReturn_mean,PSTH_SacOnReturn_mean,PSTH_SacOnReturn_sem,'lineprops',{[0.9290    0.6940    0.1250]},'transparent',1,'patchSaturation',0.3);

hold on

yrange=ylim;

shadedErrorBar([-PlotSacPreBuffer:EyeBinWidth: PlotSacPostBuffer],EyeEcc_AlignedAfter_mean+yrange(2)*1.5,EyeEcc_AlignedAfter_sem,'lineprops',{[0.9290    0.6940    0.1250]},'transparent',1,'patchSaturation',0.3);
yrange=ylim;
plot([0,0],yrange,'--k');
set(findall(gca, 'Type', 'Line'),'LineWidth',3);

set(gca,'FontSize',15,'FontWeight','Bold');
set(gca,'LineWidth',3);
xlabel('Time from Return Saccade Onset','FontSize',25,'FontWeight','Bold');
ylabel('Spikes/s','FontSize',25,'FontWeight','Bold');


%%%%%%%%%Align on Return Saccade Onset: Good vs Bad %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(1,2,2);
shadedErrorBar(TimeSequence_SacOnReturn_mean,PSTH_SacOnReturn_Good_mean,PSTH_SacOnReturn_Good_sem,'lineprops',{[1    0    0]},'transparent',1,'patchSaturation',0.3);

hold on
shadedErrorBar(TimeSequence_SacOnReturn_mean,PSTH_SacOnReturn_Bad_mean,PSTH_SacOnReturn_Bad_sem,'lineprops',{[0    0    1]},'transparent',1,'patchSaturation',0.3);

yrange=ylim;
yrange(2)=max(PSTH_SacOnReturn_mean(TimeSequence_SacOnReturn_mean<PlotSacPostBuffer & TimeSequence_SacOnReturn_mean>-PlotSacPreBuffer))*1.05;

shadedErrorBar([-PlotSacPreBuffer:-PlotSacPreBuffer+length(EyeEcc_AlignedAfter_Good_mean)-1],EyeEcc_AlignedAfter_Good_mean+yrange(2)*1.5,EyeEcc_AlignedAfter_Good_sem,'lineprops',{[1    0    0]},'transparent',1,'patchSaturation',0.3);
shadedErrorBar([-PlotSacPreBuffer:-PlotSacPreBuffer+length(EyeEcc_AlignedAfter_Bad_mean)-1],EyeEcc_AlignedAfter_Bad_mean+yrange(2)*1.5,EyeEcc_AlignedAfter_Bad_sem,'lineprops',{[0    0    1]},'transparent',1,'patchSaturation',0.3);
ylim([0,yrange(2)*1.5+20]);

xlim([-PlotSacPreBuffer, PlotSacPostBuffer]);%Keep in the same scale
yrange=ylim;
plot([0,0],yrange,'--k');

set(findall(gca, 'Type', 'Line'),'LineWidth',3);
set(gca,'FontSize',15,'FontWeight','Bold');
set(gca,'LineWidth',3);

xlabel('Time from Return Saccade Onset','FontSize',25,'FontWeight','Bold');
ylabel('Spikes/s','FontSize',25,'FontWeight','Bold');




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
