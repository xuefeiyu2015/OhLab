function OneDRAnalysis(Data, StartTrial, EndTrial,ShowFigureFlag,OutputFlag);

%Function to analysis neural responses during one direciton reward task
%Xuefei Yu 01172023
%{
%Section to export file for students
export_path = '/Users/xuefeiyu/Documents/XuefeiFile/WorkRelated/Data/FEF_Inactivation_OHLab';
% Prepare data for export
InjectionCondition = zeros(Data.TrialNum,1);%Before injection
InjectionCondition(161:end) = NaN; %After injection


Data.InjectionCondition = InjectionCondition;

output_file = fullfile(export_path, 'Adams110725FEFInactivation_OneDR.mat');
save(output_file, 'Data');
%disp('file exported')
%}

%keyboard


ECode;
% Lab colour convention, location value: high = purple, low = green
[ColorHigh, ColorLow] = ValueColors('Location');
SegEvent=[TGTCD];%Event Flow of interest
%SegEventEnd=[SHOWFIXCD,FIXOFF,TGTOFF,RWDOFFCD ,NaN,EXTRAFP];%Event end
SegEventEnd=[TGTOFF];%Event end
SegEventStr={'TargetOn'};
AmpMatchFlag=0;%1:Amplitude match on ;0:Amplitude Match off
FileName_raw = Data.FileName; 

if FileName_raw(1) == 'A'
    Monkey_curr = 2;%Adams
else
    Monkey_curr = 1;%Robin
end


PlotStartOffset=[1000];%in ms
PlotStopOffset=[8000];%in ms

ExamplePlot=0;%PLots for presentations, turn off for regular use

ResponseLatencyThreshold=3;

%Parameter to define the sliding window to measure firing rates
StepSize=10;%in ms
BinWidth=10;%in ms


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

%TargetNum=Data.TargetNum(Selection,2:end);
%{
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
%}
%% Target location

% Target position. 2-D layout (Robin): [ID x y] / [ID angle ecc].
% 3-D layout (Adams): trials x [ID target1 target2] x [x y] / [angle ecc];
% indexing (:,2:end) on it mixes target 2 into the columns (ecc became the
% angle of target 2), so read target 1 explicitly.
if ndims(Data.TargetPolarLoc) == 3
    TargetPolar = [Data.TargetPolarLoc(Selection,2,1), Data.TargetPolarLoc(Selection,2,2)];
else
    TargetPolar = Data.TargetPolarLoc(Selection,2:end);
end
if ndims(Data.TargetLoc) == 3
    TargetLoc = [Data.TargetLoc(Selection,2,1), Data.TargetLoc(Selection,2,2)];
else
    TargetLoc = Data.TargetLoc(Selection,2:end);
end
TargetAngle = TargetPolar(:,1);
TargetEcc =  TargetPolar(:,2);

UniqueAngle = uniquetol(TargetAngle,0.1);
UniqueEcc = uniquetol(TargetEcc,0.1);





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
try
ObjectIndex=cell2mat(ReproduceFromEvent(EventChannel,[GOODOBJ  BADOBJ]));

GoodObjectIndex=(ObjectIndex==GOODOBJ)';
catch
    keyboard
end

% Drop the first n trials after each value switch so adaptation trials are
% not included in Good/Bad averages.
% A switch (value reversal) is where the high-value side changes:
% high-value side = target angle on high-value trials, the opposite angle on
% low-value trials. This works for both TargetIndex layouts; the previous
% diff(TargetIndex(:,2)) rule found no switch in the 3-D layout (Adams),
% where that column is all zeros.
nSwitchExclude = 5;

GoodSide = round(TargetAngle(:));
GoodSide(~GoodObjectIndex(:)) = mod(GoodSide(~GoodObjectIndex(:)) + 180, 360);
SwitchIdx = find(diff(GoodSide)~=0)+1;
DropMask = false(numel(GoodSide),1);
for idx = 1:length(SwitchIdx)
    DropMask(SwitchIdx(idx):min(SwitchIdx(idx)+nSwitchExclude-1,numel(GoodSide))) = true;
end
KeepMask = ~DropMask;
if any(DropMask)
    fprintf('Excluding %d trials (first %d after each of %d value switches)\n', ...
        sum(DropMask), nSwitchExclude, numel(SwitchIdx));
end

% All-trial copies, used only to draw the trial-order rasters with the
% excluded switch trials in grey. All analyses below use the kept trials.
AllTrials.Selection        = Selection;
AllTrials.KeepMask         = KeepMask;
AllTrials.IsGood           = GoodObjectIndex(:);
AllTrials.TargetAngle      = TargetAngle(:);
AllTrials.EyeX             = EyeChannel_X_L;
AllTrials.EyeY             = EyeChannel_Y_L;
AllTrials.EventChannel     = EventChannel;
AllTrials.EventTimeChannel = EventTimeChannel;
AllTrials.EventBin         = EventBin;

% Apply exclusion mask to all trial-wise variables
Selection = Selection(KeepMask);
GoodObjectIndex = GoodObjectIndex(KeepMask);
TargetLoc = TargetLoc(KeepMask,:);
TargetPolar = TargetPolar(KeepMask,:);
TargetAngle = TargetAngle(KeepMask);
TargetEcc = TargetEcc(KeepMask);
EyeChannel_X_L = EyeChannel_X_L(KeepMask);
EyeChannel_Y_L = EyeChannel_Y_L(KeepMask);
Eye_Eccentricity = Eye_Eccentricity(KeepMask);
EventChannel = EventChannel(KeepMask,:);
EventTimeChannel = EventTimeChannel(KeepMask,:);
EventBin = EventBin(KeepMask,:);
UniqueAngle = uniquetol(TargetAngle,0.1);
UniqueEcc = uniquetol(TargetEcc,0.1);


%Select the period between target onset and target off set
FractalOnsetBin = FindOutTime(EventChannel,EventBin,TGTCD );


FixOffsetBin = FindOutTime(EventChannel,EventBin,FIXOFF);

AverageDelay  = nanmean(FixOffsetBin - FractalOnsetBin);
disp('Overlap delay: ms')
disp(AverageDelay )



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
% Get the mean firing rate
FR_All = (sum(~isnan(SpikeTime),2)-2)./SpikeTime(end)*1000;
FR_Mean = mean(FR_All);
FR_Std = std(FR_All);

%Select the period between target onset and target off set
FractalOnsetBin = FindOutTime(EventChannel,EventBin,TGTCD );
FractalOffsetBin = FindOutTime(EventChannel,EventBin,FIXOFF)+1000;
FractalOnsetTime = FindOutTime(EventChannel,EventTimeChannel,TGTCD );
%FractalOffsetTime = FindOutTime(EventChannel,EventTimeChannel,TGTOFF );

EyeChannel_X_Selection= SelectEyeInterval(EyeChannel_X_L,FractalOnsetBin ,FractalOffsetBin )';
EyeChannel_Y_Selection= SelectEyeInterval(EyeChannel_Y_L,FractalOnsetBin ,FractalOffsetBin )';

EyeChannel_X_Selection=ReorganizeEye(EyeChannel_X_Selection);
EyeChannel_Y_Selection=ReorganizeEye(EyeChannel_Y_Selection);

%EyeChannel_X_Selection=EyeChannel_X_Selection(NoNNaNSelector,:);
%EyeChannel_Y_Selection=EyeChannel_Y_Selection(NoNNaNSelector,:);





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


%SpikeTimeAligned_TgtOn=SpikeTimeAligned_TgtOn(NoNNaNSelector,:);

%Align again on saccade onset
AlignedSpike=SpikeAlignTool(RT_First_Target',SpikeTimeAligned_TgtOn,StepSize,BinWidth);

PSTH_SacOn=AlignedSpike.PSTH_Aligned;
TimeSequence_SacOn=AlignedSpike.PSTH_Time_Aligned;
%Select interval
SelectPSTHinterval_SacOn=[-PlotSacPreBuffer,PlotSacPostBuffer];
[PSTH_SacOn,TimeSequence_SacOn]=SelectPSTHInterval(PSTH_SacOn,TimeSequence_SacOn,SelectPSTHinterval_SacOn);

%Raster
Raster_SacOn_whole = AlignedSpike.SpikeTime_Aligned_Whole;

%Change raster
Raster_SacOn= SelectSpikeTimeInterval(Raster_SacOn_whole,SelectPSTHinterval_SacOn(1)*ones(size(Raster_SacOn_whole,1),1),SelectPSTHinterval_SacOn(end)*ones(size(Raster_SacOn_whole,1),1));


PSTH_SacOn_Good = PSTH_SacOn(GoodObjectIndex,:);
PSTH_SacOn_Bad = PSTH_SacOn(~GoodObjectIndex,:);



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
%Eye_Eccentricity_Tgt=Eye_Eccentricity_Tgt(NoNNaNSelector,:);

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

%SpikeCount
PSTH_TgtOn_SpikeCount = SegmentData.PSTH_Event_SpikeCount{1}/(StepSize*10^-3);


%%Raster
%Raster 
Raster_TgtOn = SegmentData.RasterAfterAlign{1};

% Display only: rasters for ALL trials (incl. excluded switch trials), so
% the trial-order rasters can show them in grey. Only the excluded trials
% are aligned here (same alignment as above); kept trials reuse the rasters
% already computed, then both are merged back into trial order.
if ShowFigureFlag
    Excl = ~AllTrials.KeepMask;
    ExclTrials = AllTrials;
    for fld = {'EyeX', 'EyeY'}
        ExclTrials.(fld{1}) = AllTrials.(fld{1})(Excl);
    end
    for fld = {'EventChannel', 'EventTimeChannel', 'EventBin'}
        ExclTrials.(fld{1}) = AllTrials.(fld{1})(Excl, :);
    end
    SpikeTimeAll = squeeze(Data.SpikeTimeData(spk,:,AllTrials.Selection(Excl)))';
    RasterCfg = struct('SegEvent', SegEvent, 'SegEventEnd', SegEventEnd, ...
        'PlotStartOffset', PlotStartOffset, 'PlotStopOffset', PlotStopOffset, ...
        'StepSize', StepSize, 'BinWidth', BinWidth, 'EyeBinWidth', EyeBinWidth, ...
        'V_Threshold', V_Threshold, 'ContinueBin', ContinueBin, ...
        'ISI_Threshold', ISI_Threshold, 'AmplitudeThreshold', AmplitudeThreshold, ...
        'TimeThreshold', TimeThreshold, 'SacInterval', SelectPSTHinterval_SacOn);
    if any(Excl)
        if sum(Excl) == 1
            SpikeTimeAll = SpikeTimeAll(:)';           % squeeze of one trial
        end
        [ExclRaster_TgtOn, ExclRaster_SacOn] = BuildAllTrialRasters(ExclTrials, SpikeTimeAll, RasterCfg);
    else
        ExclRaster_TgtOn = []; ExclRaster_SacOn = [];
    end
    RasterAll_TgtOn = MergeTrialRows(Raster_TgtOn, ExclRaster_TgtOn, AllTrials.KeepMask);
    RasterAll_SacOn = MergeTrialRows(Raster_SacOn, ExclRaster_SacOn, AllTrials.KeepMask);
end

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
StartTimeBin(StartTimeBin<1)=1;
EndTimeBin=FindOutTime(EventChannel,EventTimeChannel,SegEventEnd )+PlotStopOffset;
  
Eye_Eccentricity_FractalOn=SelectEyeInterval(Eye_Eccentricity,StartTimeBin,EndTimeBin)';
Eye_Eccentricity_FractalOn=ReorganizeEye(Eye_Eccentricity_FractalOn);
%Eye_Eccentricity_FractalOn=Eye_Eccentricity_FractalOn(NoNNaNSelector,:);




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


%{
%%
%{
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

%}


%%
%{
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
%{
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
%}

%% Loop over the target directions
for i = 1:length(UniqueAngle )
    %for j = 1:length(UniqueEcc)
        sel = TargetAngle == UniqueAngle(i) ;
        
        %Align on target on 
        currFracOn_Loc = MeanPSTH_FR_FractalOn(sel,:);
        curr_GoodFracOn_Loc = MeanPSTH_FR_FractalOn(sel & GoodObjectIndex,:);
        Good_FracOn_Loc_Mean{i} = nanmean(curr_GoodFracOn_Loc,1);
        Good_FracOn_Loc_Sem{i} = nanstd(curr_GoodFracOn_Loc,[],1)./sqrt(sum(~isnan(curr_GoodFracOn_Loc),1));
        
        curr_BadFracOn_Loc = MeanPSTH_FR_FractalOn(sel & ~GoodObjectIndex,:);
        Bad_FracOn_Loc_Mean{i} = nanmean(curr_BadFracOn_Loc,1);
        Bad_FracOn_Loc_Sem{i} = nanstd(curr_BadFracOn_Loc,[],1)./sqrt(sum(~isnan(curr_BadFracOn_Loc),1));
        
        
        
        AngleEccorder{i} = [UniqueAngle(i),UniqueEcc(1)];
        
         %Align on sac on 
        curr_SacOn_Loc = PSTH_SacOn(sel,:);
        curr_GoodSacOn_Loc = PSTH_SacOn(sel & GoodObjectIndex,:);
        Good_SacOn_Loc_Mean{i} = nanmean(curr_GoodSacOn_Loc,1);
        Good_SacOn_Loc_Sem{i} = nanstd(curr_GoodSacOn_Loc,[],1)./sqrt(sum(~isnan(curr_GoodSacOn_Loc),1));
        
        curr_BadSacOn_Loc = PSTH_SacOn(sel & ~GoodObjectIndex,:);
        Bad_SacOn_Loc_Mean{i} = nanmean(curr_BadSacOn_Loc,1);
        Bad_SacOn_Loc_Sem{i} = nanstd(curr_BadSacOn_Loc,[],1)./sqrt(sum(~isnan(curr_BadSacOn_Loc),1));

        %%
        %{
        TimePreTgt = [-200,-80];
        Ave_FR_preTgtOn_Good = nanmean(curr_GoodFracOn_Loc(:,PSTH_Time_FraOn_All >= TimePreTgt(1) & PSTH_Time_FraOn_All <= TimePreTgt(2)),2);
        Ave_FR_preTgtOn_Bad = nanmean(curr_BadFracOn_Loc(:,PSTH_Time_FraOn_All >= TimePreTgt(1) & PSTH_Time_FraOn_All <= TimePreTgt(2)),2);
        %}
        TimeTgt = [100,200];
        Ave_FR_TgtOn_Good = nanmean(curr_GoodFracOn_Loc(:,PSTH_Time_FraOn_All >= TimeTgt(1) & PSTH_Time_FraOn_All <= TimeTgt(2)),2);
        Ave_FR_TgtOn_Bad = nanmean(curr_BadFracOn_Loc(:,PSTH_Time_FraOn_All >= TimeTgt(1) & PSTH_Time_FraOn_All <= TimeTgt(2)),2);

        TimeSac = [-200,0];

        Ave_FR_SacOn_Good = nanmean(curr_GoodSacOn_Loc(:,TimeSequence_SacOn_mean >= TimeSac(1) & TimeSequence_SacOn_mean <= TimeSac(2)),2);
        Ave_FR_SacOn_Bad = nanmean(curr_BadSacOn_Loc(:,TimeSequence_SacOn_mean >= TimeSac(1) & TimeSequence_SacOn_mean <= TimeSac(2)),2);

       % [h,p_PreTgt(i)] = ttest2(Ave_FR_preTgtOn_Good,Ave_FR_preTgtOn_Bad);
       % [h,p_Tgt(i)]= ttest2(Ave_FR_TgtOn_Good,Ave_FR_TgtOn_Bad);
       % [h,p_Sac(i)]= ttest2(Ave_FR_SacOn_Good,Ave_FR_SacOn_Bad);
%{
        preTgtGood(i) = mean(Ave_FR_preTgtOn_Good);
        preTgtBad(i) = mean(Ave_FR_preTgtOn_Bad);

        TgtGood(i) = mean(Ave_FR_TgtOn_Good);
        TgtBad(i) = mean(Ave_FR_TgtOn_Bad);

        SacGood(i) = mean(Ave_FR_SacOn_Good);
        SacBad(i) = mean(Ave_FR_SacOn_Bad);
%}
        %}
        %For visual responses
        VisTgt = [60,200];
        Ave_FR_TgtOn_Vis = nanmean(currFracOn_Loc(:,PSTH_Time_FraOn_All >= VisTgt(1) & PSTH_Time_FraOn_All <= VisTgt(2)),2);
        SponTime = [-500,-400];
        

        Ave_FR_TgtOn_Spon = nanmean(currFracOn_Loc(:,PSTH_Time_FraOn_All >= SponTime(1) & PSTH_Time_FraOn_All <= SponTime(2)),2);

        Vis_Dff(i) = mean(Ave_FR_TgtOn_Vis)-mean(Ave_FR_TgtOn_Spon)';

        [h,p_Vis(i)]= ttest(Ave_FR_TgtOn_Vis ,Ave_FR_TgtOn_Spon);
        Vis_Dff(i) = mean(Ave_FR_TgtOn_Vis)-mean(Ave_FR_TgtOn_Spon)';
%{
        Ave_FR_SacOn_Good = nanmean(curr_GoodSacOn_Loc(:,TimeSequence_SacOn_mean >= TimeSac(1) & TimeSequence_SacOn_mean <= TimeSac(2)),2);
        Ave_FR_SacOn_Bad = nanmean(curr_BadSacOn_Loc(:,TimeSequence_SacOn_mean >= TimeSac(1) & TimeSequence_SacOn_mean <= TimeSac(2)),2);
    
        [h,p_Tgt(i)]= ttest2(Ave_FR_TgtOn_Good,Ave_FR_TgtOn_Bad);
        [h,p_Sac(i)]= ttest2(Ave_FR_SacOn_Good,Ave_FR_SacOn_Bad);
%}
         %For saccade responses 

        SacTgt = [-200,20];
        Ave_FR_SacOn = nanmean(curr_SacOn_Loc(:,TimeSequence_SacOn_mean >= VisTgt(1) & TimeSequence_SacOn_mean <= VisTgt(2)),2);
        
       

        [h,p_pureSac(i)]= ttest(Ave_FR_SacOn,Ave_FR_TgtOn_Spon);
        pureSac_Diff(i) = mean(Ave_FR_SacOn)-mean(Ave_FR_TgtOn_Spon)';

%
        N = 100;
        permuteN = 1000;

        [TgtAUC(i),bestz,perm] = rocN(Ave_FR_TgtOn_Good,Ave_FR_TgtOn_Bad,N,permuteN);
        p_AUC_Tgt(i) = perm.pValue;
        p_AUC_Tgt_ranksum(i) = ranksum(Ave_FR_TgtOn_Good, Ave_FR_TgtOn_Bad);


        [SacAUC(i),bestz,perm] = rocN(Ave_FR_SacOn_Good,Ave_FR_SacOn_Bad,N,permuteN);
        p_AUC_Sac(i) = perm.pValue;
%}
         %Raster
        Raster_TgtOn_Good{i}=Raster_TgtOn(sel & GoodObjectIndex,:);
        Raster_TgtOn_Bad{i}=Raster_TgtOn(sel & ~GoodObjectIndex,:);
        Raster_TgtOn_Ordered{i} = Raster_TgtOn(sel,:);      % all trials in original order
        ValueLabel_TgtOn{i}     = GoodObjectIndex(sel);      % 1=high value, 0=low value

        Raster_SacOn_Good{i}=Raster_SacOn(sel & GoodObjectIndex,:);
        Raster_SacOn_Bad{i}=Raster_SacOn(sel & ~GoodObjectIndex,:);
        Raster_SacOn_Ordered{i} = Raster_SacOn(sel,:);
        ValueLabel_SacOn{i}     = GoodObjectIndex(sel);
        if ShowFigureFlag
            % All trials at this location in trial order, incl. excluded
            selAll = abs(AllTrials.TargetAngle - UniqueAngle(i)) < 0.1;
            RasterAll_TgtOn_Ordered{i} = RasterAll_TgtOn(selAll,:);
            RasterAll_SacOn_Ordered{i} = RasterAll_SacOn(selAll,:);
            ValueLabel_All{i}          = AllTrials.IsGood(selAll);
            Kept_All{i}                = AllTrials.KeepMask(selAll);
        end

        %Time sequence comparison
        % for tgt onset
        %N=100;
        %permuteN=1000;
        for m = 1:length(PSTH_Time_FraOn_All)
            good_tgt=curr_GoodFracOn_Loc(:,m);
            bad_tgt =curr_BadFracOn_Loc(:,m);
            [h,p_tgttime_value_loc(i,m)]=ttest2(good_tgt,bad_tgt);

            h_tgttime_value_loc(i,m) = p_tgttime_value_loc(i,m)<0.05;

            goodbaddiff = mean(good_tgt)-mean(bad_tgt);
            goodpreferindex = goodbaddiff>0;
            
            h_tgttime_value_loc_good(i,m) = goodpreferindex & h_tgttime_value_loc(i,m);
            h_tgttime_value_loc_bad(i,m) = ~goodpreferindex & h_tgttime_value_loc(i,m);

            [AUC_tgttime_loc(i,m),~,perm] = rocN(good_tgt,bad_tgt,N,permuteN);
            AUC_p_time(i,m) = perm.pValue;
           
            valid_g = ~isnan(good_tgt);
       valid_b =  ~isnan(bad_tgt);
       if sum(valid_g)>3 & sum(valid_b) >3
           AUC_p_time_ranksum(i,m) = ranksum(good_tgt(valid_g),bad_tgt(valid_b));
       else

            AUC_p_time_ranksum(i,m) = NaN;
       end
            


        end

        %for saccade onset
        for m = 1:length(TimeSequence_SacOn_mean)
            good_sac=curr_GoodSacOn_Loc(:,m);
            bad_sac =curr_BadSacOn_Loc(:,m);
            [h,p_sactime_value_loc(i,m)]=ttest2(good_sac,bad_sac);

            h_sactime_value_loc(i,m) = p_sactime_value_loc(i,m)<0.05;

            goodbaddiff = mean(good_sac)-mean(bad_sac);
            goodpreferindex = goodbaddiff>0;
            h_sactime_value_loc_good(i,m) = goodpreferindex & h_sactime_value_loc(i,m);
            h_sactime_value_loc_bad(i,m) = ~goodpreferindex & h_sactime_value_loc(i,m);
            AUC_sactime_loc(i,m) = rocN(good_sac,bad_sac);


        end

%{
        % Peak AUC around target onset (value-related)
        sigma_ms = 30;
        win = 50;
        NROC = N;
        NpermutePeak = 2000;
        auc_sm_whole = GaussSmooth(PSTH_Time_FraOn_All, AUC_tgttime_loc(i,:), sigma_ms);
        %Choose a time interval to search for the peak
        PeakSearch = [-500,400];
        auc_sm = auc_sm_whole(PSTH_Time_FraOn_All>=PeakSearch(1) & PSTH_Time_FraOn_All<=PeakSearch(2));
        auc_seach_time = PSTH_Time_FraOn_All(PSTH_Time_FraOn_All>=PeakSearch(1) & PSTH_Time_FraOn_All<=PeakSearch(2));
        % Use max absolute deviation from 0.5 so neurons preferring "bad"
        % (AUC < 0.5) are correctly identified — not just max(AUC).
        [~, peakIdx] = max(abs(auc_sm - 0.5));
        %peakAUC_sm(i,1) = auc_sm(peakIdx);
        peakTime_AUC(i,1) = auc_seach_time(peakIdx);
        sel_win = PSTH_Time_FraOn_All >= peakTime_AUC(i)-win & PSTH_Time_FraOn_All <= peakTime_AUC(i)+win;
        %peakAUC_avg(i,1) = mean(AUC_tgttime_loc(i,sel_win),'omitnan');
        %TrialIndexValue = GoodObjectIndex(sel);
        %p_peakAUC(i,1) = permutationROCpeak(TrialIndexValue,currFracOn_Loc,PSTH_Time_FraOn_All,peakTime_AUC(i),peakAUC_avg(i), win, NROC, NpermutePeak);
%}
        % Latency of ROC (first of 3 consecutive significant bins in the same
        % direction as the peak — handles both good-preferring (AUC>0.5) and
        % bad-preferring (AUC<0.5) neurons)
        NContinuous = 3;
        alphaROC = 0.05;
        sigMask = p_tgttime_value_loc(i,:) < alphaROC;
        if TgtAUC(i) > 0.5
            dirMask = AUC_tgttime_loc(i,:) > 0.5; % good-preferring neuron
        else
            dirMask = AUC_tgttime_loc(i,:) < 0.5; % bad-preferring neuron
        end
        LatencyTimeAUC(i) = NaN;
        for idxL = 1:(length(PSTH_Time_FraOn_All)-NContinuous+1)
            segSig = sigMask(idxL:idxL+NContinuous-1);
            segDir = dirMask(idxL:idxL+NContinuous-1);
            if all(segSig) && all(segDir)
                LatencyTimeAUC(i) = PSTH_Time_FraOn_All(idxL);
                break;
            end
        end

            sigMask_ranksum = AUC_p_time_ranksum(i,:) < alphaROC;
            
            LatencyTimeAUC_ranksum(i) = NaN;

            for idxL = 1:(length(PSTH_Time_FraOn_All)-NContinuous+1)
                segSig = sigMask_ranksum(idxL:idxL+NContinuous-1);
                segDir = dirMask(idxL:idxL+NContinuous-1);
                if all(segSig) && all(segDir)
                    LatencyTimeAUC_ranksum(i) = PSTH_Time_FraOn_All(idxL);
                    break;
                end
            end


        
   % end
end %End of the target direction loop

%RF Fix as contralateral
%{
meanAll = TgtGood + TgtBad;
max_idx = meanAll == max(meanAll);
RF = UniqueAngle(max_idx);
%}
%{
if Monkey_curr == 1 %Robin
    max_idx = UniqueAngle >= 0 & UniqueAngle < 180;
else
    max_idx = UniqueAngle >= 180 & UniqueAngle <360;

end
RF = UniqueAngle(max_idx);
nonRF = ~max_idx;
%}

if length(UniqueAngle) >=2
if   Vis_Dff(1) > Vis_Dff(2)% (TgtGood(1)+TgtBad(1)) > (TgtGood(2)+TgtBad(2)) 
    max_idx  = UniqueAngle==UniqueAngle(1);
else
    max_idx  = UniqueAngle==UniqueAngle(2);
end
else
    max_idx  = UniqueAngle==UniqueAngle(1);
end

RF = UniqueAngle(max_idx);
nonRF = ~max_idx;

%{
TgtGood_RF = TgtGood(max_idx);
TgtBad_RF = TgtBad(max_idx);

p_Tgt_RF = p_Tgt(max_idx);

PreTgtGood_RF = preTgtGood(max_idx);
PreTgtBad_RF = preTgtBad(max_idx);

p_preTgt_RF = p_PreTgt(max_idx);

SacGood_RF = SacGood(max_idx);
SacBad_RF = SacBad(max_idx);

p_Sac_RF = p_Sac(max_idx);
%}
TgtAUC_RF = TgtAUC(max_idx);
p_AUC_Tgt_RF =p_AUC_Tgt(max_idx);
p_AUC_Tgt_RF_ranksum = p_AUC_Tgt_ranksum(max_idx);

SacAUC_RF = SacAUC(max_idx);
p_AUC_Sac_RF =p_AUC_Sac(max_idx);
%}
%peakTime_AUC_RF = peakTime_AUC(max_idx);
%peakAUC_avg_RF = peakAUC_avg(max_idx);
%p_peakAUC_RF = p_peakAUC(max_idx);
LatencyTimeAUC_RF = LatencyTimeAUC(max_idx);
LatencyTimeAUC_ranksum_RF = LatencyTimeAUC_ranksum(max_idx);

%{
TgtGood_nonRF = TgtGood(nonRF);
TgtBad_nonRF = TgtBad(nonRF);

p_Tgt_nonRF = p_Tgt(nonRF);

PreTgtGood_nonRF = preTgtGood(nonRF);
PreTgtBad_nonRF = preTgtBad(nonRF);

p_preTgt_nonRF = p_PreTgt(nonRF);

SacGood_nonRF = SacGood(nonRF);
SacBad_nonRF = SacBad(nonRF);

p_Sac_nonRF = p_Sac(nonRF);
%}
TgtAUC_nonRF = TgtAUC(nonRF);
p_AUC_Tgt_nonRF =p_AUC_Tgt(nonRF);
p_AUC_Tgt_nonRF_ranksum = p_AUC_Tgt_ranksum(nonRF);

SacAUC_nonRF = SacAUC(nonRF);
p_AUC_Sac_nonRF =p_AUC_Sac(nonRF);

p_Vis_RF = p_Vis(max_idx);
VisDiff_RF = Vis_Dff(max_idx);


p_pureSac_RF =p_pureSac(max_idx);
pureSacDiff_RF =pureSac_Diff(max_idx);

%{
p_Tgt_nonRF  = p_Tgt(nonRF);
    p_Sac_nonRF  = p_Sac(nonRF);

    TgtAUC_nonRF = TgtAUC(nonRF);
    p_AUC_Tgt_nonRF = p_AUC_Tgt(nonRF);
%}
p_Vis_nonRF = p_Vis(nonRF);
VisDiff_nonRF = Vis_Dff(nonRF);

p_pureSac_nonRF =p_pureSac(nonRF);
pureSacDiff_nonRF =pureSac_Diff(nonRF);

%peakTime_AUC_nonRF = peakTime_AUC(nonRF);
%peakAUC_avg_nonRF = peakAUC_avg(nonRF);
%p_peakAUC_nonRF = p_peakAUC(nonRF);
LatencyTimeAUC_nonRF = LatencyTimeAUC(nonRF);
LatencyTimeAUC_ranksum_nonRF = LatencyTimeAUC_ranksum(nonRF);






ExportInterval = [-600,400];
sel = PSTH_Time_FraOn_All >= ExportInterval(1) & PSTH_Time_FraOn_All <= ExportInterval(2); 
SelectTime=sel;
Time_TgtOn=PSTH_Time_FraOn_All(sel);

PSTH_Good_TgtOn_RF = Good_FracOn_Loc_Mean{max_idx}(:,sel);
PSTH_Bad_TgtOn_RF = Bad_FracOn_Loc_Mean{max_idx}(:,sel);

PSTH_Good_TgtOn_nonRF = Good_FracOn_Loc_Mean{~max_idx}(:,sel);
PSTH_Bad_TgtOn_nonRF = Bad_FracOn_Loc_Mean{~max_idx}(:,sel);

PSTH_Good_TgtOn_RF_z = (PSTH_Good_TgtOn_RF - FR_Mean)/FR_Std;
PSTH_Bad_TgtOn_RF_z = (PSTH_Bad_TgtOn_RF - FR_Mean)/FR_Std;  

PSTH_Good_TgtOn_nonRF_z = (PSTH_Good_TgtOn_nonRF - FR_Mean)/FR_Std;
PSTH_Bad_TgtOn_nonRF_z = (PSTH_Bad_TgtOn_nonRF - FR_Mean)/FR_Std;


hPSTH_TgtOn_Value_RF = h_tgttime_value_loc(max_idx,sel);

hPSTH_TgtOn_Good_RF = h_tgttime_value_loc_good(max_idx,sel);
hPSTH_TgtOn_Bad_RF = h_tgttime_value_loc_bad(max_idx,sel);

hPSTH_TgtOn_Value_nonRF = h_tgttime_value_loc(nonRF,sel);

hPSTH_TgtOn_Good_nonRF = h_tgttime_value_loc_good(nonRF,sel);
hPSTH_TgtOn_Bad_nonRF = h_tgttime_value_loc_bad(nonRF,sel);

 AUCPSTH_TgtOn_Value_RF = AUC_tgttime_loc(max_idx,sel);
 AUCPSTH_TgtOn_Value_nonRF = AUC_tgttime_loc(nonRF,sel);
 




ExportInterval_sac = [-200,400];
sel = TimeSequence_SacOn_mean >= ExportInterval_sac(1) & TimeSequence_SacOn_mean <= ExportInterval_sac(2); 
Time_SacOn=TimeSequence_SacOn_mean(sel); 

PSTH_SacGood_RF = Good_SacOn_Loc_Mean{max_idx}(:,sel);
PSTH_SacBad_RF = Bad_SacOn_Loc_Mean{max_idx}(:,sel);


PSTH_SacGood_nonRF = Good_SacOn_Loc_Mean{~max_idx}(:,sel);
PSTH_SacBad_nonRF = Bad_SacOn_Loc_Mean{~max_idx}(:,sel);

hPSTH_SacOn_Value_RF = h_sactime_value_loc(max_idx,sel);

hPSTH_SacOn_Good_RF = h_sactime_value_loc_good(max_idx,sel);
hPSTH_SacOn_Bad_RF = h_sactime_value_loc_bad(max_idx,sel);

 
    hPSTH_SacOn_Value_nonRF = h_sactime_value_loc(nonRF,sel);

    hPSTH_SacOn_Good_nonRF = h_sactime_value_loc_good(nonRF,sel);
    hPSTH_SacOn_Bad_nonRF = h_sactime_value_loc_bad(nonRF,sel);

    AUCPSTH_SacOn_Value_RF = AUC_sactime_loc(max_idx,sel);

    
    AUCPSTH_SacOn_Value_nonRF = AUC_sactime_loc(nonRF,sel);


%z-score
PSTH_SacGood_RF_z = (PSTH_SacGood_RF - FR_Mean)/FR_Std;
PSTH_SacBad_RF_z = (PSTH_SacBad_RF - FR_Mean)/FR_Std;

PSTH_SacGood_nonRF_z = (PSTH_SacGood_nonRF - FR_Mean)/FR_Std;
PSTH_SacBad_nonRF_z = (PSTH_SacBad_nonRF - FR_Mean)/FR_Std;
%{
TgtGood_RF_z = (TgtGood_RF - FR_Mean)/FR_Std;
PreTgtGood_RF_z = (PreTgtGood_RF - FR_Mean)/FR_Std;

TgtBad_RF_z = (TgtBad_RF - FR_Mean)/FR_Std;
PreTgtBad_RF_z = (PreTgtBad_RF - FR_Mean)/FR_Std;

%}

%%
        if ShowFigureFlag
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Figure 1: Average firing rate
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%for each location

%{
figtitlestr{1}='PSTHEachLocation';

fig{1}=PrepareFigure(110,'w',[50,100, 1200,800],'Name',figtitlestr{1});

for i = 1:length(UniqueAngle)
    subplot(length(UniqueAngle),2,(i-1)*length(UniqueAngle)+1);%Target onset

    shadedErrorBar(PSTH_Time_FraOn_All,Good_FracOn_Loc_Mean{i},Good_FracOn_Loc_Sem{i},'lineprops',{ColorHigh},'transparent',1,'patchSaturation',0.3);

    hold on
    shadedErrorBar(PSTH_Time_FraOn_All,Bad_FracOn_Loc_Mean{i},Bad_FracOn_Loc_Sem{i},'lineprops',{ColorLow},'transparent',1,'patchSaturation',0.3);

       xlim([-PlotSacPreBuffer, PlotSacPostBuffer]);%Keep in the same scale
       set(findall(gca, 'Type', 'Line'),'LineWidth',3);

        set(gca,'FontSize',15,'FontWeight','Bold');
        set(gca,'LineWidth',3);
        title(sprintf("Angle: %.0f,Ecc: %1.1f",AngleEccorder{i}));
        if i==1
            xlabel('Target Onset(ms)','FontSize',15,'FontWeight','Bold');
            
        end
        
        ylabel('Firing Rate(Hz)','FontSize',15,'FontWeight','Bold');

         subplot(length(UniqueAngle),2,(i-1)*length(UniqueAngle)+2);%Saccade onset

          shadedErrorBar(TimeSequence_SacOn_mean,Good_SacOn_Loc_Mean{i},Good_SacOn_Loc_Sem{i},'lineprops',{ColorHigh},'transparent',1,'patchSaturation',0.3);

        hold on
          shadedErrorBar(TimeSequence_SacOn_mean,Bad_SacOn_Loc_Mean{i},Bad_SacOn_Loc_Sem{i},'lineprops',{ColorLow},'transparent',1,'patchSaturation',0.3);

       % yrange=ylim;
       % yrange(2)=max(MeanPSTH_FR_FractalOn_mean(PSTH_Time_FraOn_All<PlotSacPostBuffer & PSTH_Time_FraOn_All>-PlotSacPreBuffer))*1.05;
          xlim([-PlotSacPreBuffer, PlotSacPostBuffer]);%Keep in the same scale
       set(findall(gca, 'Type', 'Line'),'LineWidth',3);

        set(gca,'FontSize',15,'FontWeight','Bold');
        set(gca,'LineWidth',3);
        title(sprintf("Angle: %.0f,Ecc: %1.1f",AngleEccorder{i}));
        if i==1
            xlabel('Time from Saccade Onset(ms)','FontSize',15,'FontWeight','Bold');
            
        end
        
        ylabel('Firing Rate(Hz)','FontSize',15,'FontWeight','Bold');
        


end %End of the uniqueAngle


figtitlestr{2}='PSTHEachLocationWithRaster';

fig{2}=PrepareFigure(111,'w',[50,100, 1200,800],'Name',figtitlestr{2});

StartLoc = 1;
LineWidth = 2;
Barlength = 1;

for i = 1:length(UniqueAngle)
    subplot(2,length(UniqueAngle),(i-1)*2+1);
    SpikeRaterPlot(Raster_TgtOn_Good{i},StartLoc,Barlength,ColorHigh,LineWidth);
    
 hold on
 IncreaseBar = size(Raster_TgtOn_Good{i},1);
 
    SpikeRaterPlot(Raster_TgtOn_Bad{i},StartLoc+IncreaseBar,Barlength,ColorLow,LineWidth);
    yrange = ylim;
 hold on
 plot([0,0],yrange,'--k','LineWidth',3);
 xlim([-PlotSacPreBuffer, PlotSacPostBuffer]);
 xlabel('Time from Fractal Onset(ms)','FontSize',15,'FontWeight','Bold');
  title(sprintf("Angle: %1.1d",AngleEccorder{i}));
 

    subplot(2,length(UniqueAngle),(i-1)*2+2);
    SpikeRaterPlot(Raster_SacOn_Good{i},StartLoc,Barlength,ColorHigh,LineWidth);
    
 hold on
 IncreaseBar = size(Raster_SacOn_Good{i},1);
 
    SpikeRaterPlot(Raster_SacOn_Bad{i},StartLoc+IncreaseBar,Barlength,ColorLow,LineWidth);
    yrange = ylim;
 hold on
 plot([0,0],yrange,'--k','LineWidth',3);
xlim([-PlotSacPreBuffer, PlotSacPostBuffer])

xlabel('Time from Saccade Onset(ms)','FontSize',15,'FontWeight','Bold');


    

end
%}
%Put all the figures into one for better visualization
figtitlestr{1}='PSTH_GoodBad_EachLocation';

fig{1}=PrepareFigure(110,'w',[50,100, 1200,800],'Name',figtitlestr{1});

StartLoc = 1;
LineWidth = 2;
Barlength = 1;

    for i = 1:length(UniqueAngle)
        subplot(4,length(UniqueAngle),i);
        % Raster of ALL trials in trial order: purple = high value, green = low
        % value, grey = first nSwitchExclude trials after a value switch
        % (excluded from all analyses)
        nT_tgt = size(RasterAll_TgtOn_Ordered{i}, 1);
        Xred=[]; Yred=[]; Xblue=[]; Yblue=[]; Xgrey=[]; Ygrey=[];
        for t = 1:nT_tgt
            row = RasterAll_TgtOn_Ordered{i}(t,:);
            y0 = StartLoc+t-1; y1 = y0+Barlength;
            Xt = repmat(row,2,1);
            Yt = repmat([y0;y1],1,size(row,2));
            if ~Kept_All{i}(t)
                Xgrey=[Xgrey,Xt]; Ygrey=[Ygrey,Yt];
            elseif ValueLabel_All{i}(t)
                Xred=[Xred,Xt]; Yred=[Yred,Yt];
            else
                Xblue=[Xblue,Xt]; Yblue=[Yblue,Yt];
            end
        end
        hold on;
        if ~isempty(Xgrey), plot(Xgrey, Ygrey, 'Color', [0.7 0.7 0.7], 'LineWidth', LineWidth); end
        if ~isempty(Xred),  plot(Xred,  Yred,  'Color', ColorHigh, 'LineWidth', LineWidth); end
        if ~isempty(Xblue), plot(Xblue, Yblue, 'Color', ColorLow, 'LineWidth', LineWidth); end
        ylim([Barlength, StartLoc-1+nT_tgt+Barlength*2]);
        yrange = ylim;
        hold on
        plot([0,0],yrange,'--k','LineWidth',3);
        set(findall(gca, 'Type', 'Line'),'LineWidth',3);
     %   xlim([-PlotSacPreBuffer, PlotSacPostBuffer]);
     xlim([-200, 400]);%Keep in the same scale
     
        title(sprintf("Angle: %.0f,Ecc: %.1f",AngleEccorder{i}));

        subplot(4,length(UniqueAngle),length(UniqueAngle)+i);
        %PSTH align on targt onset
        shadedErrorBar(PSTH_Time_FraOn_All,Good_FracOn_Loc_Mean{i},Good_FracOn_Loc_Sem{i},'lineprops',{ColorHigh},'transparent',1,'patchSaturation',0.3);
        hold on
        shadedErrorBar(PSTH_Time_FraOn_All,Bad_FracOn_Loc_Mean{i},Bad_FracOn_Loc_Sem{i},'lineprops',{ColorLow},'transparent',1,'patchSaturation',0.3);

         hold on
         yrange2 = ylim;
        plot([0,0],yrange2,'--k','LineWidth',3);
        

       % xlim([-PlotSacPreBuffer, PlotSacPostBuffer]);%Keep in the same scale
       xlim([-200, 400]);%Keep in the same scale
       set(findall(gca, 'Type', 'Line'),'LineWidth',3);

        set(gca,'FontSize',15,'FontWeight','Bold');
        set(gca,'LineWidth',3);
    

        xlabel('Time from Fractal Onset(ms)','FontSize',15,'FontWeight','Bold');
        ylim([0,200]);
%Raster align on saccade on
        subplot(4,length(UniqueAngle),2*length(UniqueAngle)+i);
        % Raster of ALL trials in trial order: purple = high value, green = low
        % value, grey = first nSwitchExclude trials after a value switch
        % (excluded from all analyses)
        nT_sac = size(RasterAll_SacOn_Ordered{i}, 1);
        Xred=[]; Yred=[]; Xblue=[]; Yblue=[]; Xgrey=[]; Ygrey=[];
        for t = 1:nT_sac
            row = RasterAll_SacOn_Ordered{i}(t,:);
            y0 = StartLoc+t-1; y1 = y0+Barlength;
            Xt = repmat(row,2,1);
            Yt = repmat([y0;y1],1,size(row,2));
            if ~Kept_All{i}(t)
                Xgrey=[Xgrey,Xt]; Ygrey=[Ygrey,Yt];
            elseif ValueLabel_All{i}(t)
                Xred=[Xred,Xt]; Yred=[Yred,Yt];
            else
                Xblue=[Xblue,Xt]; Yblue=[Yblue,Yt];
            end
        end
        hold on;
        if ~isempty(Xgrey), plot(Xgrey, Ygrey, 'Color', [0.7 0.7 0.7], 'LineWidth', LineWidth); end
        if ~isempty(Xred),  plot(Xred,  Yred,  'Color', ColorHigh, 'LineWidth', LineWidth); end
        if ~isempty(Xblue), plot(Xblue, Yblue, 'Color', ColorLow, 'LineWidth', LineWidth); end
        ylim([Barlength, StartLoc-1+nT_sac+Barlength*2]);
        yrange = ylim;
        hold on
        plot([0,0],yrange,'--k','LineWidth',3);
       % xlim([-PlotSacPreBuffer, PlotSacPostBuffer]);
       xlim([-200, 200]);
        set(findall(gca, 'Type', 'Line'),'LineWidth',3);
        title(sprintf("Angle: %.0f,Ecc: %.1f",AngleEccorder{i}));


        subplot(4,length(UniqueAngle),3*length(UniqueAngle)+i);
        %PSTH align on sacccade onset
        shadedErrorBar(TimeSequence_SacOn_mean,Good_SacOn_Loc_Mean{i},Good_SacOn_Loc_Sem{i},'lineprops',{ColorHigh},'transparent',1,'patchSaturation',0.3);

        hold on
        shadedErrorBar(TimeSequence_SacOn_mean,Bad_SacOn_Loc_Mean{i},Bad_SacOn_Loc_Sem{i},'lineprops',{ColorLow},'transparent',1,'patchSaturation',0.3);
       % yrange = ylim;
         plot([0,0],yrange2,'--k','LineWidth',3);

      %  xlim([-PlotSacPreBuffer, PlotSacPostBuffer]);%Keep in the same scale
       xlim([-200, 200]);
       set(findall(gca, 'Type', 'Line'),'LineWidth',3);

        set(gca,'FontSize',15,'FontWeight','Bold');
        set(gca,'LineWidth',3);
    

        xlabel('Time from Saccade Onset(ms)','FontSize',15,'FontWeight','Bold');
       
      
    end

     info = sprintf('%s_C%d_p(RFVis):%1.3f;p(NonRFVis):%1.3f;p_pureSac(RF):%1.3f;p_pureSac(NonRF):%1.3f;',FileName_raw,SpikeChannel,p_Vis_RF,p_Vis_nonRF,p_pureSac_RF,p_pureSac_nonRF);
     info = sprintf('%s\nRaster: purple = high, green = low, grey = first %d trials after each value switch (excluded, %d trials)', ...
         info, nSwitchExclude, sum(~AllTrials.KeepMask));

      %  showingtitle=sprintf('%s_%s_N%s_C%s',MonkeyName,string(RecordDateOriginal),string(NeuronNum),string(ChannelNum));
      annotation('textbox', [0.1 0.9 0.8 0.1], 'String', info, 'EdgeColor', 'none', ...
           'HorizontalAlignment', 'center', 'FontSize', 10, 'Interpreter', 'none');

      set(gcf, 'Renderer', 'painters')
%{
% Export fig{1} as vectorized EPS to Desktop
DesktopPath = fullfile(getenv('HOME'), 'Desktop');
EpsFileName = fullfile(DesktopPath, ['RasterPSTH_onedrFEF_loconly.eps']);
print(fig{1}, EpsFileName, '-depsc', '-painters');
disp(['Exported EPS to: ' EpsFileName]);

%}
%keyboard
%{
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

for i = 1:length(UniqueAngle )
    for j = 1:length(UniqueEcc)
        if UniqueEcc(j) == 0
            index = center;
        else
            [minValue,closestIndex] = min(abs(UniqueAngle(i)-angle));
            index = subloc(closestIndex) ;
        end
        subplot(layout(1),layout(2),index);
        
        shadedErrorBar(PSTH_Time_FraOn_All,Good_FracOn_Loc_Mean{i,j},Good_FracOn_Loc_Sem{i,j},'lineprops',{ColorHigh},'transparent',1,'patchSaturation',0.3);

        hold on
        shadedErrorBar(PSTH_Time_FraOn_All,Bad_FracOn_Loc_Mean{i,j},Bad_FracOn_Loc_Sem{i,j},'lineprops',{ColorLow},'transparent',1,'patchSaturation',0.3);

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
        subplot(layout(1),layout(2),index);
        
        shadedErrorBar(TimeSequence_SacOn_mean,Good_SacOn_Loc_Mean{i,j},Good_SacOn_Loc_Sem{i,j},'lineprops',{ColorHigh},'transparent',1,'patchSaturation',0.3);

        hold on
        shadedErrorBar(TimeSequence_SacOn_mean,Bad_SacOn_Loc_Mean{i,j},Bad_SacOn_Loc_Sem{i,j},'lineprops',{ColorLow},'transparent',1,'patchSaturation',0.3);

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
%}

        end %End of the if plotfigure

%% Output to files
%Setup output directory
Workingdirectory=pwd;
MarkerFolder='DataAnalysis';
%MarkerFolder='DataHub';
Flag=strfind(Workingdirectory,MarkerFolder);
BasicDirectory=Workingdirectory(1:Flag+length(MarkerFolder));


%%%%%%%%%%%Export data%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Set up output path
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
OutputFigureName=sprintf('%s_%s_N%s_C%s_%CID_%NID',MonkeyName,string(RecordDateOriginal),string(NeuronNum),string(ChannelNum),string(ChannelID),string(ClusterID));



if TaskCode == ONEDR_TASK  
    ProtocolMarker='OneDR';
%ProtocolMarker='NovelObject_Low';

elseif  TaskCode==OBJ_LOCATIONVALUE_TASK
    ProtocolMarker='ObjLocationValue';
end


ExistFlag=0;
%Load the old file if exist
if exist(OutputFileName)

    load(OutputFileName);
    %Load OutputData into the memory
    if isfield(OutputData,ProtocolMarker)
    ExistFlag=1;
    end
end




%{
%FileName=strcat(FileName,'_RearchingTime.mat');
if SaveFlag
save(FileName,'OutputData');
disp('Seaching Time Data Exported');
else
    disp('No data export');
end
%}
%OutputData For current analysis
OutputData_New=[];
if TaskCode == ONEDR_TASK
    OutputData_New.OneDR.TrialType='GoodTrials';
elseif  TaskCode == OBJ_LOCATIONVALUE_TASK
    OutputData_New.OBJLOCVAL.TrialType='GoodTrials';
end






%%%%%%%%%Save overall scene response characteristics%%%%%%%%%%%%%%%%%%%%%%%
%Organize the data

% Per-trial decoder variables (all conditions, raw firing rate, no z-score)
sel_RF_trials = ismember(TargetAngle, UniqueAngle(max_idx)); % logical: which trials are at RF
psth_targeton_raw = PSTH_TgtOn_SpikeCount(:, SelectTime);           % all trials x T, raw FR
target_index      = double(GoodObjectIndex);                  % 1=high value, 0=low value, all trials
RF_index          = double(sel_RF_trials);                    % 1=into RF, 0=out of RF, all trials
n_high_RF         = sum( GoodObjectIndex &  sel_RF_trials);
n_low_RF          = sum(~GoodObjectIndex &  sel_RF_trials);
n_high_nonRF      = sum( GoodObjectIndex & ~sel_RF_trials);
n_low_nonRF       = sum(~GoodObjectIndex & ~sel_RF_trials);


NameStr={
    'PSTH_Good_TgtOn_RF','PSTH_Bad_TgtOn_RF','PSTH_Good_TgtOn_nonRF','PSTH_Bad_TgtOn_nonRF',...
    'PSTH_SacGood_nonRF','PSTH_SacBad_nonRF','PSTH_SacGood_RF','PSTH_SacBad_RF',...
    'Time_TgtOn','Time_SacOn',...
    'PSTH_Good_TgtOn_RF_z','PSTH_Bad_TgtOn_RF_z','PSTH_Good_TgtOn_nonRF_z','PSTH_Bad_TgtOn_nonRF_z',...
    'PSTH_SacGood_RF_z','PSTH_SacBad_RF_z','PSTH_SacGood_nonRF_z','PSTH_SacBad_nonRF_z',...  
    'p_Vis_RF','VisDiff_RF','p_Vis_nonRF','VisDiff_nonRF','p_pureSac_RF','p_pureSac_nonRF','pureSacDiff_RF','pureSacDiff_nonRF',...
    'hPSTH_TgtOn_Value_RF','hPSTH_TgtOn_Good_RF','hPSTH_TgtOn_Bad_RF','hPSTH_SacOn_Value_RF','hPSTH_SacOn_Good_RF','hPSTH_SacOn_Bad_RF',...
    'hPSTH_TgtOn_Value_nonRF','hPSTH_TgtOn_Good_nonRF','hPSTH_TgtOn_Bad_nonRF','hPSTH_SacOn_Value_nonRF','hPSTH_SacOn_Good_nonRF','hPSTH_SacOn_Bad_nonRF',...
    'AUCPSTH_TgtOn_Value_RF','AUCPSTH_TgtOn_Value_nonRF','AUCPSTH_SacOn_Value_RF','AUCPSTH_SacOn_Value_nonRF',...
    'UniqueAngle',...
     'LatencyTimeAUC_RF','LatencyTimeAUC_nonRF',...
    'psth_targeton_raw','target_index','RF_index',...
    'n_high_RF','n_low_RF','n_high_nonRF','n_low_nonRF',...
    'TgtAUC_RF','p_AUC_Tgt_RF','TgtAUC_nonRF','p_AUC_Tgt_nonRF',...
    'p_AUC_Tgt_RF_ranksum','p_AUC_Tgt_nonRF_ranksum','LatencyTimeAUC_ranksum_RF','LatencyTimeAUC_ranksum_nonRF'


    };

DataLib={
    PSTH_Good_TgtOn_RF,PSTH_Bad_TgtOn_RF,PSTH_Good_TgtOn_nonRF,PSTH_Bad_TgtOn_nonRF,...
    PSTH_SacGood_nonRF,PSTH_SacBad_nonRF,PSTH_SacGood_RF,PSTH_SacBad_RF,...
    Time_TgtOn,Time_SacOn,...
    PSTH_Good_TgtOn_RF_z,PSTH_Bad_TgtOn_RF_z,PSTH_Good_TgtOn_nonRF_z,PSTH_Bad_TgtOn_nonRF_z,...
    PSTH_SacGood_RF_z,PSTH_SacBad_RF_z,PSTH_SacGood_nonRF_z ,PSTH_SacBad_nonRF_z,...
    p_Vis_RF,VisDiff_RF,p_Vis_nonRF,VisDiff_nonRF,p_pureSac_RF,p_pureSac_nonRF,pureSacDiff_RF,pureSacDiff_nonRF,...
    hPSTH_TgtOn_Value_RF,hPSTH_TgtOn_Good_RF,hPSTH_TgtOn_Bad_RF,hPSTH_SacOn_Value_RF,hPSTH_SacOn_Good_RF,hPSTH_SacOn_Bad_RF,...
    hPSTH_TgtOn_Value_nonRF,hPSTH_TgtOn_Good_nonRF,hPSTH_TgtOn_Bad_nonRF,hPSTH_SacOn_Value_nonRF,hPSTH_SacOn_Good_nonRF,hPSTH_SacOn_Bad_nonRF,...
    AUCPSTH_TgtOn_Value_RF,AUCPSTH_TgtOn_Value_nonRF,AUCPSTH_SacOn_Value_RF,AUCPSTH_SacOn_Value_nonRF,...
    UniqueAngle,...
    LatencyTimeAUC_RF,LatencyTimeAUC_nonRF,...
    psth_targeton_raw,target_index,RF_index,...
    n_high_RF,n_low_RF,n_high_nonRF,n_low_nonRF,...
    TgtAUC_RF,p_AUC_Tgt_RF,TgtAUC_nonRF,p_AUC_Tgt_nonRF,...
    p_AUC_Tgt_RF_ranksum,p_AUC_Tgt_nonRF_ranksum,LatencyTimeAUC_ranksum_RF,LatencyTimeAUC_ranksum_nonRF

    };
%{
TgtGood_nonRF = TgtGood(nonRF);
TgtBad_nonRF = TgtBad(nonRF);

p_Tgt_nonRF = p_Tgt(nonRF);

PreTgtGood_nonRF = preTgtGood(nonRF);
PreTgtBad_nonRF = preTgtBad(nonRF);

p_preTgt_nonRF = p_PreTgt(nonRF);

SacGood_nonRF = SacGood(nonRF);
SacBad_nonRF = SacBad(nonRF);

p_Sac_nonRF = p_Sac(nonRF);

TgtAUC_nonRF = TgtAUC(nonRF);
p_AUC_Tgt_nonRF =p_AUC_Tgt(nonRF);

SacAUC_nonRF = SacAUC(nonRF);
p_AUC_Sac_nonRF =p_AUC_Sac(nonRF);
%}

%{
TgtGood_RF = TgtGood(max_idx);
TgtBad_RF = TgtBad(max_idx);

p_Tgt_RF = p_Tgt(max_idx);

PreTgtGood_RF = preTgtGood(max_idx);
PreTgtBad_RF = preTgtBad(max_idx);

p_preTgt_RF = p_PreTgt(max_idx);

SacGood_RF = SacGood(max_idx);
SacBad_RF = SacBad(max_idx);

p_Sac_RF = p_Sac(max_idx);


PSTH_Good_TgtOn_RF = Good_FracOn_Loc_Mean{max_idx};
PSTH_Bad_TgtOn_RF = Bad_FracOn_Loc_Mean{max_idx};

PSTH_SacGood_RF = Good_SacOn_Loc_Mean{max_idx};
PSTH_SacBad_RF = Bad_SacOn_Loc_Mean{max_idx};

PSTH_Good_TgtOn_nonRF = Good_FracOn_Loc_Mean{~max_idx};
PSTH_Bad_TgtOn_nonRF = Bad_FracOn_Loc_Mean{~max_idx};

PSTH_SacGood_nonRF = Good_SacOn_Loc_Mean{~max_idx};
PSTH_SacBad_nonRF = Bad_SacOn_Loc_Mean{~max_idx};

%}

DataStamp=containers.Map(NameStr,DataLib);



if TaskCode==ONEDR_TASK

OutputData_New.OneDR.Task=TaskType;
OutputData_New.OneDR.TaskCode=TaskCode;
OutputData_New.OneDR.DataStamp=DataStamp;

if ExistFlag
%Compare the old one with the new one
Task_Old=OutputData.OneDR.TaskCode;

     if Task_Old~=TaskCode
        %Not the same task, add the new dataset into the old one
        OutputData.OneDR(numel(OutputData.OneDR)+1)=OutputData_New.OneDR;
     else
         %If the same task, replace the old one with the new one
         OutputData.OneDR=OutputData_New.OneDR;
    
    
     end
else
    %Output the current file
   OutputData.OneDR=OutputData_New.OneDR; 
end


else
    OutputData_New.OBJLOCVAL.Task=TaskType;
    OutputData_New.OBJLOCVAL.TaskCode=TaskCode;
    OutputData_New.OBJLOCVAL.DataStamp=DataStamp;

    if ExistFlag
    %Compare the old one with the new one
        Task_Old=OutputData.OBJLOCVAL.TaskCode;

        if Task_Old~=TaskCode
        %Not the same task, add the new dataset into the old one
            OutputData.OBJLOCVAL(numel(OutputData.OBJLOCVAL)+1)=OutputData_New.OBJLOCVAL;
        else
         %If the same task, replace the old one with the new one
            OutputData.OBJLOCVAL=OutputData_New.OBJLOCVAL;
    
    
        end
    else
    %Output the current file
        OutputData.OBJLOCVAL=OutputData_New.OBJLOCVAL; 
    end


end


%Store figures in common folders
if TaskCode == ONEDR_TASK  
    ProtocolMarker='OneDR';
%ProtocolMarker='NovelObject_Low';

elseif  TaskCode==OBJ_LOCATIONVALUE_TASK
    ProtocolMarker='ObjLocationValue';
end


%ProtocolMarker='OneDR';


FigurePath=strcat(BasicDirectory,'Figures',BasicDirectory(end),ProtocolMarker,BasicDirectory(end));%,OutputFolerName);
if ~exist(FigurePath)
    mkdir(FigurePath)
end

if ShowFigureFlag & OutputFlag
    for i=1:length(fig)
        cd(FigurePath);
        OutputFigureNameFig=[];
        titlestr=figtitlestr{i};
        FolderName=strcat(titlestr,'_FEF');
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
    
    for jj=1:TotalFigure
        if ~isempty(fig{jj})
        
         OutputFigureName_Final=strcat(OutputFigureName,figtitlestr{jj},'.jpg');
         saveas(fig{jj},OutputFigureName_Final);
        end
        
            
        
    end

    disp('Figures have been exported to the neuron folder');

end

%}








if OutputFlag
    cd(OutPath);
    save(OutputFileName,'OutputData');
    disp('Data Exported');
else
    disp('No data export');
end  



    end %End of the spike channel loop
end %End of the protocol

function All = MergeTrialRows(Kept, Excl, KeepMask)
%MERGETRIALROWS Put kept-trial and excluded-trial raster rows back into
%trial order; rows NaN-padded to a common width. Display only.
nCol = max(size(Kept, 2), size(Excl, 2));
All = nan(numel(KeepMask), nCol);
All(KeepMask, 1:size(Kept, 2)) = Kept;
if ~isempty(Excl)
    All(~KeepMask, 1:size(Excl, 2)) = Excl;
end
end

function [RasterTgt, RasterSac] = BuildAllTrialRasters(AllTrials, SpikeTime, Cfg)
%BUILDALLTRIALRASTERS Target-onset and saccade-onset aligned rasters for a
%set of trials (here: the excluded switch trials), using the same alignment
%as the main OneDRAnalysis pipeline. Display only.
%
% Inputs:
%   AllTrials : struct with EyeX, EyeY, EventChannel, EventTimeChannel,
%               EventBin for the trials to align
%   SpikeTime : [nTrial x nSpike] spike times of one channel, all trials
%   Cfg       : alignment / saccade parameters copied from OneDRAnalysis
% Outputs:
%   RasterTgt : spike times aligned to target onset (EventSegTool)
%   RasterSac : spike times aligned to first target saccade onset,
%               limited to Cfg.SacInterval
ECode;
EC = AllTrials.EventChannel; ET = AllTrials.EventTimeChannel; EB = AllTrials.EventBin;

% First target saccade onset (same rule as the main pipeline)
OnBin  = FindOutTime(EC, EB, TGTCD);
OffBin = FindOutTime(EC, EB, FIXOFF) + 1000;
OnTime = FindOutTime(EC, ET, TGTCD);
X = ReorganizeEye(SelectEyeInterval(AllTrials.EyeX, OnBin, OffBin)');
Y = ReorganizeEye(SelectEyeInterval(AllTrials.EyeY, OnBin, OffBin)');
Sac = SelectSaccade(X, Y, Cfg.EyeBinWidth, Cfg.V_Threshold, Cfg.ContinueBin, Cfg.ISI_Threshold);
RT = nan(numel(Sac), 1);
for t = 1:numel(Sac)
    if Sac(t).NumOfSaccade > 0
        k = find(Sac(t).SaccadeAmplitude(:) > Cfg.AmplitudeThreshold & ...
                 Sac(t).SaccadeStartTime(:) > Cfg.TimeThreshold, 1, 'first');
        if ~isempty(k)
            RT(t) = Sac(t).SaccadeStartTime(k);
        end
    end
end

% Saccade-onset aligned raster
A1 = SpikeAlignTool(OnTime, SpikeTime, Cfg.StepSize, Cfg.BinWidth);
A2 = SpikeAlignTool(RT, A1.SpikeTime_Aligned_Whole, Cfg.StepSize, Cfg.BinWidth);
Raw = A2.SpikeTime_Aligned_Whole;
n = size(Raw, 1);
RasterSac = SelectSpikeTimeInterval(Raw, Cfg.SacInterval(1)*ones(n,1), Cfg.SacInterval(end)*ones(n,1));

% Target-onset aligned raster
Seg = EventSegTool(Cfg.SegEvent, Cfg.SegEventEnd, Cfg.PlotStartOffset, Cfg.PlotStopOffset, ...
    EC, ET, SpikeTime, Cfg.StepSize, Cfg.BinWidth);
RasterTgt = Seg.RasterAfterAlign{1};
end

