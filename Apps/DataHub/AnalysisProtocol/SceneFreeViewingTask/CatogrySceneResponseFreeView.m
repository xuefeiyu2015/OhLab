function CatogrySceneResponseFreeView(Data, StartTrial, EndTrial,ShowFigureFlag,OutputFlag);

%%%%%%%%%%%%%%%%%User controller starts here
ECode;
legendOn=1;

%Event of interest
%{
SegEvent=[SCENE_ON,SHOWFIXCD,TGTCD,RWDCD ,SCENE_OFF,SCENE_OFF];%Event Flow of interest
%SegEventEnd=[SHOWFIXCD,FIXOFF,TGTOFF,RWDOFFCD ,NaN,EXTRAFP];%Event end
SegEventEnd=[SHOWFIXCD,FIXOFF,TGTOFF,RWDOFFCD ,NaN,NaN];%Event end
SegEventStr={'SceneOn','FixOn','TgtOn','RewardOn','SceneOff','Spon'};

PlotStartOffset=[100,100,100,100,100,-1000];%in ms
PlotStopOffset=[100,100,100,100,200,2000];%in ms
%}

SegEvent=[SCENE_ON,SCENE_OFF];%Event Flow of interest
%SegEventEnd=[SHOWFIXCD,FIXOFF,TGTOFF,RWDOFFCD ,NaN,EXTRAFP];%Event end
%SegEventEnd=[SHOWFIXCD];%Event end
%SegEventEnd=[SCENE_OFF,NaN ];%Event end
SegEventEnd=[NaN,NaN ];%Event end
SegEventStr={'SceneOn','SceneOff'};




GlobalPlot=0;%Flag to show the plot for each scene and the average for all scenes

PlotStartOffset=[1000,500];%in ms
PlotStopOffset=[2000,2000];%in ms

plotrange=[-PlotStartOffset(1),2000];

plotrangeoff=[500,1000+PlotStopOffset(1)];

ResponseLatencyThreshold=3;

%Parameter to define the sliding window to measure firing rates
StepSize=10;%in ms
BinWidth=20;%in ms
%{
%Parameter seting for interval selecting plot period
StartOffset=100;%in ms
StopOffset=50;%in ms
%For scene offset
StartOffset2=500;%in ms
StopOffset2=50;%in ms
IntervalAfterMarkerStart=1000;%in ms;%This for align marker
IntervalAfterMarkerStartOff=500;%in ms;%For scene offset
%}
%Color Definition for each scenes

Color=[lines(6);prism(6);pink(6)];
%Color(2,:)=copper(1)+0.5;
Color(9,:)=[0.28,0.78,0.89];
Color=[Color;hsv(10)];
Color(18,:)=[0.1,0.457,0.194];
ColorForScene=Color;


ColorForEachMode=[0.83 0.14 0.14
             1.00 0.54 0.00
             0.47 0.25 0.80
             0.25 0.80 0.54];
AssociateBlock=[10,11,0,1];
Blockstr={'FixLoc,FixSceneSeq','FixLoc,RandSceneSeq','RandLoc,FixSceneSeq','RandLoc,RandSceneSeq'};

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
%
if  strcmp(Data.FileName,'Robin110820SCENEPV1.rst')
    EndTrial=min(Selection(end),178);
    Selection=Selection(3):EndTrial;
    
end
%}
%{
if  strcmp(Data.FileName,'Robin013121SCENEPV1.rst')
    RemoveTrial=[10,251];
    Selection=Selection(~ismember(Selection,RemoveTrial));
end
%}
%
if  strcmp(Data.FileName,'Robin112620SCENEPV1.rst')
    RemoveTrial=[121,186];
    Selection=Selection(~ismember(Selection,RemoveTrial));
end
%}
%SceneType
FixLocScene=1010:1017;
RandLocScene=1000:1007;
DuelLocScene=[1054:1057,1060:1063];
NovelSceneInterval=[2000,3999];
NaturalSceneInterval=[2000,2999];
MapSceneInterval=[3000,3999];



TrialNum=length(Selection);



%FixationPointIndex=Data.FixPointIndex(Selection,:);
%FixLoc=Data.FixLoc(Selection,:,:);


%TargetIndex=Data.TargetIndex(Selection,:,:);
%TargetLoc=Data.TargetLoc(Selection,:,:);

%TargetNum=Data.TargetNum(Selection,:);
%TargetSize=Data.TargetSize(Selection,:);
%TargetPolarLoc=Data.TargetPolarLoc(Selection,:,:);
%ChoiceOutcome=Data.ChoiceOutcome(Selection,:);


%GoodObjectLib=Data.GoodObjectLib;
%BadObjectLib=Data.BadObjectLib;
%SceneLib=Data.SceneLib;

%GoodObject=arrayfun(@(x)( ismember(x,GoodObjectLib(:,3))),TargetIndex(:,2:end,2));%Use real name
%GoodObjectID=TargetIndex(:,2:end,2);

%BadObject=arrayfun(@(x)( ismember(x,BadObjectLib(:,3))),TargetIndex(:,2:end,2));%Use real name
%BadObjectID=TargetIndex(:,2:end,2);
%GoodBadIndex=GoodObject*1+BadObject*2;%1:Good Object;2:Bad Object

SceneIndexType=Data.SceneIndex;
if size(size(SceneIndexType),2)==3
    SceneIndex=Data.SceneIndex(Selection,2:end,1);
    SceneName=Data.SceneIndex(Selection,2:end,2);
else
    SceneIndex=Data.SceneIndex(Selection,2);
    SceneName=Data.SceneIndex(Selection,3);
    
end 


%SceneIndex=SceneIndex(:,[6:8]);%Scene Name
%EyeChannel_X=Data.EyeChannel_X_L(Selection,:); %Only import the first two channels
%EyeChannel_Y=Data.EyeChannel_Y_L(Selection,:);
%EyeBinWidth=Data.EyeBinWidth;
%EyeTime=Data.EyeTime(Selection,:);

EventChannel=Data.EventChannel(Selection,:);
EventTimeChannel=Data.EventTimeChannel(Selection,:);
EventBin=Data.EventBin(Selection,:);

%SpikeChannel=Data.SpikeChannel;
%SpikeTime=squeeze(Data.SpikeTimeData(:,:,Selection))';
%SpikeBin=squeeze(Data.SpikeBinData(:,:,Selection))';

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

%Segregate spike data from scene onset to fixation onset

SegmentData=EventSegTool(SegEvent,SegEventEnd,PlotStartOffset,PlotStopOffset,EventChannel,EventTimeChannel,SpikeTime,StepSize,BinWidth);

%Retrive data from the field
%Total average for each event
MeanFR=SegmentData.MeanFR; %[Mean,sem]%Averall mean and sem
MeanFRRaw=SegmentData.MeanFRRaw{1};%=AverageFiringRate_Event;

%PSTH for each event
MeanPSTH=SegmentData.MeanPSTH;%=[PSTH_Event_Time',PSTH_Event'];
MeanPSTH_Time=MeanPSTH{1,1};
MeanPSTH_FR=MeanPSTH{1,2};

MeanPSTH_SceneOff_Time=MeanPSTH{2,1};
MeanPSTH_SceneOff_FR=MeanPSTH{2,2};


%Relative time bin for each event
RelativeTimeMarker=SegmentData.RelativeTimeMarker{1};%=EventTimeMarker;
%Absolute time marker for each event
AbosoluteTimeMarker=SegmentData.AbsoluteTimeMarker{1};

RelativeTimeMarkerOff=SegmentData.RelativeTimeMarker{2};%=EventTimeMarker;
%Absolute time marker for each event
AbosoluteTimeMarkerOff=SegmentData.AbsoluteTimeMarker{2};

%Mean time
SceneOnMarker=nanmean(RelativeTimeMarker);

%Reward time
RwdTimeAbs=FindOutTime(EventChannel,EventBin,RWDCD);
RwdTimeRel=RwdTimeAbs-AbosoluteTimeMarker(:,1);
RwdTimeRelMean=nanmean(RwdTimeRel);

%{
%As there is only on scene onset for each trial, we don't need the
following 

%Average firing rate for each event for each incidance in each trial
MeanFREach=SegmentData.MeanFREach;%=AverageFiringRate_Event_Each;

%Average PSTH for each event for each incidance in each trial
MeanPSTHEach=SegmentData.MeanPSTHEach;%=[PSTH_Event_Each_Time',PSTH_Event_Each'];

%Event time marker for each event for each incidance in each trial
RelativeTimeMarkerEach=SegmentData.RelativeTimeMarkerEach;%=EventTimeMarkerEach;
AbosoluteTimeMarkerEach=SegmentData.AbosoluteTimeMarkerEach;%=AbosoluteTimeMarkerEach;
%}
%
%Get the average firing rate 500ms before scene onset as the baseline
SponEvent=[SCENE_ON];
PreTimeBuffer=500;
SponEventEnd=NaN;
StartOffset=0;
StopOffset=PreTimeBuffer;

EventTimeChannelForSpon=EventTimeChannel-PreTimeBuffer;


SegmentDataSpon=EventSegTool(SponEvent,SponEventEnd,StartOffset,StopOffset,EventChannel,EventTimeChannelForSpon,SpikeTime,StepSize,BinWidth);

Mean_Spon=SegmentDataSpon.MeanFR;
Mean_Spon_ForEachTrial=SegmentDataSpon.MeanFRRaw{1};


%MeanPSTH_FR=MeanPSTH_FR-Mean_Spon_ForEachTrial;%Minus spon
%MeanFRRaw=MeanFRRaw-Mean_Spon_ForEachTrial;%Minus spon


%Find out when the firing rate start to show scene response

[p_All,h_All]=ttest2_group(MeanPSTH_FR,Mean_Spon_ForEachTrial,'both');

p_All(1:ceil(PlotStartOffset/StepSize))=NaN;
h_All(1:ceil(PlotStartOffset/StepSize))=NaN;
h_All(h_All==0)=NaN;%For plot

%PSTH mean and sem

MeanPSTH_FR_mean=nanmean(MeanPSTH_FR);%Minus spon
MeanPSTH_FR_sem=nanstd(MeanPSTH_FR)./sqrt(sum(~isnan(MeanPSTH_FR),1));

%Find out response latency

ResponseLatency=FindOutLatency(nanmean(MeanPSTH_Time),h_All,ResponseLatencyThreshold);

if isempty(ResponseLatency)
    TotalResponseDuration=NaN;
    PropResponseDuration=NaN;
    
    
    FR_Sig_Whole=NaN;
    PeakFiringTime=NaN;
    FirstIntervalDuration=NaN;
    ResponseLatency=NaN;
    PeakFiring=NaN;

    
    
else
%Calculate response duration--the total time within particular period
%during the whole response window
ResponseDuration=nanmean(MeanPSTH_Time);
ResponseDuration=ResponseDuration(h_All==1);
ResponseDuration=ResponseDuration(ResponseDuration<=min(RelativeTimeMarker(:,2)) & ResponseDuration>=ResponseLatency);
TotalResponseDuration=length(ResponseDuration)*StepSize;
TotalEventDuration=min(RelativeTimeMarker(:,2));
PropResponseDuration=TotalResponseDuration/TotalEventDuration*100;




%First significant interval 
 %Peak Firing Rate in the first significant interval
Significant_TimePoint_temp=ResponseDuration(ResponseDuration>=ResponseLatency);
FirstIntervalEnd_Whole=find(diff(Significant_TimePoint_temp)>StepSize,1,'first');

    if isempty(FirstIntervalEnd_Whole) || FirstIntervalEnd_Whole>length(Significant_TimePoint_temp)
        FirstIntervalEnd_Whole=length(Significant_TimePoint_temp);
        
    end
    
    
FirstSigInterval_Whole=[1:FirstIntervalEnd_Whole];


FR_Sig_Whole=MeanPSTH_FR_mean(nanmean(MeanPSTH_Time)>=ResponseLatency);
    
[PeakFiring,loc]=max(FR_Sig_Whole(FirstSigInterval_Whole));
    
PeakFiringTime=Significant_TimePoint_temp(loc);

%{    
%First Interval Duration
FirstIntervalDuration=diff(Significant_TimePoint_temp([FirstSigInterval_Whole(1),FirstSigInterval_Whole(end)]));

%}
end 

%Rearrange scene into lines
SceneIndex_Line=SceneIndex;

SceneIndex_Line=reshape(SceneIndex_Line,numel(SceneIndex_Line),1);

NoNaNSelector=~isnan(SceneIndex_Line);

SceneIndex_Line=SceneIndex_Line(~isnan(SceneIndex_Line));

%Rearrange scene name into lines
SceneName_Line=SceneName;

SceneName_Line=reshape(SceneName_Line,numel(SceneName_Line),1);
SceneName_Line=SceneName_Line(NoNaNSelector);
%{
disp(sprintf('Number of scene measured: %d',length(SceneIndex_Line)))
disp(sprintf('Number of unique scene measured: %d',length(unique(SceneIndex_Line))))
%}


%ScrambleSceneMarker=FIRST_SCRAMBLE_SCENE-550:FIRST_SCRAMBLE_SCENE-550+11;
%GreySceneMarker=FIRST_GREY_SCENE-550:FIRST_GREY_SCENE-550+11;
ScrambleSceneMarker=FIRST_SCRAMBLE_SCENE:FIRST_SCRAMBLE_SCENE+11;
GreySceneMarker=FIRST_GREY_SCENE:FIRST_GREY_SCENE+11;

%Separate grey and scramble scene from normal scenes
%NormalSceneSelector=SceneName_Line>0;

%GreySceneSelector=ismember(SceneIndex_Line,GreySceneMarker) & ~NormalSceneSelector;
%ScrambleSceneSelector=ismember(SceneIndex_Line,ScrambleSceneMarker) & ~NormalSceneSelector;

GreySceneSelector=ismember(SceneIndex_Line,GreySceneMarker);
ScrambleSceneSelector=ismember(SceneIndex_Line,ScrambleSceneMarker);


NormalSceneSelector=~(GreySceneSelector|ScrambleSceneSelector);

GreySceneName_Line=SceneName_Line(GreySceneSelector);
ScrambleSceneName_Line=SceneName_Line(ScrambleSceneSelector);



NormalSceneName_Line=SceneName_Line(NormalSceneSelector);

%Scene type selector

FixLocSceneSelector=ismember(SceneName_Line,FixLocScene);
FixLocSceneScrambleSelector=ismember(abs(SceneName_Line),FixLocScene) & ScrambleSceneSelector;

RandLocSceneSelector=ismember(SceneName_Line,RandLocScene);
RandLocSceneScrambleSelector=ismember(abs(SceneName_Line),RandLocScene) & ScrambleSceneSelector;

DuelLocSceneSelector=ismember(SceneName_Line,DuelLocScene);
DuelLocSceneScrambleSelector=ismember(abs(SceneName_Line),DuelLocScene) & ScrambleSceneSelector;




NovelSceneSelector=(SceneName_Line>=NovelSceneInterval(1) & SceneName_Line<=NovelSceneInterval(2) );
NovelSceneScrambleSelector=(abs(SceneName_Line)>=NovelSceneInterval(1) & abs(SceneName_Line)<=NovelSceneInterval(2) ) & ScrambleSceneSelector;


NaturalSceneSelector=(SceneName_Line>=NaturalSceneInterval(1) & SceneName_Line<=NaturalSceneInterval(2) );
NaturalSceneScrambleSelector=(abs(SceneName_Line)>=NaturalSceneInterval(1) & abs(SceneName_Line)<=NaturalSceneInterval(2) ) & ScrambleSceneSelector;


MapSceneSelector=(SceneName_Line>=MapSceneInterval(1) & SceneName_Line<=MapSceneInterval(2) );
MapSceneScrambleSelector=(abs(SceneName_Line)>=MapSceneInterval(1) & abs(SceneName_Line)<=MapSceneInterval(2) ) & ScrambleSceneSelector;


FamiliarSceneSelector=(~NovelSceneSelector) & NormalSceneSelector ;
FamiliarSceneScrambleSelector=(~NovelSceneScrambleSelector) & ScrambleSceneSelector;




%Correctness check
Overlapcheck=sum((~(GreySceneSelector|ScrambleSceneSelector))&NormalSceneSelector)==sum(NormalSceneSelector);
if ~Overlapcheck
    disp('Sth wrong with the scene definition, please check the program!!');
    msgbox('Sth wrong with the scene definition, please check the program!!','Error','error');
    return
    
end

%Average according to the selected intervals
AveInterval=[0,1000];
PSTH_Time=nanmean(MeanPSTH_Time,1);

SelInterval=PSTH_Time>AveInterval(1) & PSTH_Time<AveInterval(2);
MeanFRSceneOn=nanmean(MeanPSTH_FR(:,SelInterval),2);

SponAveInterval=[-500,0];
SelIntervalSpon=PSTH_Time>SponAveInterval(1) & PSTH_Time<SponAveInterval(2);
MeanFRSceneOnSpon=nanmean(MeanPSTH_FR(:,SelIntervalSpon),2);



%Mean and PSTH for each scene type
NormalScene_Mean_EachTrial=MeanFRSceneOn(NormalSceneSelector);
NormalScene_Mean=nanmean(NormalScene_Mean_EachTrial);
NormalScene_Sem=nanstd(NormalScene_Mean_EachTrial)/sqrt(sum(NormalSceneSelector));

ScrambleScene_Mean_EachTrial=MeanFRSceneOn(ScrambleSceneSelector);
ScrambleScene_Mean=nanmean(ScrambleScene_Mean_EachTrial);
ScrambleScene_Sem=nanstd(ScrambleScene_Mean_EachTrial)/sqrt(sum(ScrambleSceneSelector));

GreyScene_Mean_EachTrial=MeanFRSceneOn(GreySceneSelector);
GreyScene_Mean=nanmean(GreyScene_Mean_EachTrial);
GreyScene_Sem=nanstd(GreyScene_Mean_EachTrial)/sqrt(sum(GreySceneSelector));

PreSceneSpon_Mean=nanmean(MeanFRSceneOnSpon);
PreSceneSpon_Sem=nanstd(MeanFRSceneOnSpon)/sqrt(sum(~isnan(MeanFRSceneOnSpon)));

%Statistics between Normal and Scramble and Grey Average
[h,p_SceneScram]=ttest2(NormalScene_Mean_EachTrial,ScrambleScene_Mean_EachTrial);
[h,p_SceneGrey]=ttest2(NormalScene_Mean_EachTrial,GreyScene_Mean_EachTrial);

%Statistics between Scene and Spon

[h,p_SceneSpon]=ttest2(NormalScene_Mean_EachTrial,MeanFRSceneOnSpon);
[h,p_ScrambleSpon]=ttest2(ScrambleScene_Mean_EachTrial,MeanFRSceneOnSpon);
[h,p_GreySpon]=ttest2(GreyScene_Mean_EachTrial,MeanFRSceneOnSpon);

%PSTH



PSTH_NormalScene_Mean_EachTrial=MeanPSTH_FR(NormalSceneSelector,:);
PSTH_NormalScene_Mean=nanmean(PSTH_NormalScene_Mean_EachTrial,1);
PSTH_NormalScene_Sem=nanstd(PSTH_NormalScene_Mean_EachTrial,[],1)./sqrt(sum(~isnan(PSTH_NormalScene_Mean_EachTrial),1));

PSTH_ScrambleScene_Mean_EachTrial=MeanPSTH_FR(ScrambleSceneSelector,:);
PSTH_ScrambleScene_Mean=nanmean(PSTH_ScrambleScene_Mean_EachTrial,1);
PSTH_ScrambleScene_Sem=nanstd(PSTH_ScrambleScene_Mean_EachTrial,[],1)./sqrt(sum(~isnan(PSTH_ScrambleScene_Mean_EachTrial),1));

PSTH_GreyScene_Mean_EachTrial=MeanPSTH_FR(GreySceneSelector,:);
PSTH_GreyScene_Mean=nanmean(PSTH_GreyScene_Mean_EachTrial,1);
PSTH_GreyScene_Sem=nanstd(PSTH_GreyScene_Mean_EachTrial,[],1)./sqrt(sum(~isnan(PSTH_GreyScene_Mean_EachTrial),1));



FixLocScene_Mean_EachTrial=MeanFRSceneOn(FixLocSceneSelector);
FixLocScene_Mean=nanmean(FixLocScene_Mean_EachTrial);
FixLocScene_Sem=nanstd(FixLocScene_Mean_EachTrial)/sqrt(sum(FixLocScene_Mean_EachTrial));

FixLocScene_Scramble_Mean_EachTrial=MeanFRSceneOn(FixLocSceneScrambleSelector);
FixLocScene_Scramble_Mean=nanmean(FixLocScene_Scramble_Mean_EachTrial);
FixLocScene_Scramble_Sem=nanstd(FixLocScene_Scramble_Mean_EachTrial)/sqrt(sum(FixLocSceneScrambleSelector));


RandLocScene_Mean_EachTrial=MeanFRSceneOn(RandLocSceneSelector);
RandLocScene_Mean=nanmean(RandLocScene_Mean_EachTrial);
RandLocScene_Sem=nanstd(RandLocScene_Mean_EachTrial)/sqrt(sum(RandLocScene_Mean_EachTrial));

RandLocScene_Scramble_Mean_EachTrial=MeanFRSceneOn(RandLocSceneScrambleSelector);
RandLocScene_Scramble_Mean=nanmean(RandLocScene_Scramble_Mean_EachTrial);
RandLocScene_Scramble_Sem=nanstd(RandLocScene_Scramble_Mean_EachTrial)/sqrt(sum(RandLocSceneScrambleSelector));

DuelScene_Mean_EachTrial=MeanFRSceneOn(DuelLocSceneSelector);
DuelScene_Mean=nanmean(DuelScene_Mean_EachTrial);
DuelScene_Sem=nanstd(DuelScene_Mean_EachTrial)/sqrt(sum(DuelScene_Mean_EachTrial));

DuelScene_Scramble_Mean_EachTrial=MeanFRSceneOn(DuelLocSceneScrambleSelector);
DuelScene_Scramble_Mean=nanmean(DuelScene_Scramble_Mean_EachTrial);
DuelScene_Scramble_Sem=nanstd(DuelScene_Scramble_Mean_EachTrial)/sqrt(sum(DuelScene_Scramble_Mean_EachTrial));


%PSTH



PSTH_FixLocScene_Mean_EachTrial=MeanPSTH_FR(FixLocSceneSelector,:);
PSTH_FixLocScene_Mean=nanmean(PSTH_FixLocScene_Mean_EachTrial,1);
PSTH_FixLocScene_Sem=nanstd(PSTH_FixLocScene_Mean_EachTrial,[],1)./sqrt(sum(~isnan(PSTH_FixLocScene_Mean_EachTrial),1));

PSTH_RandLocScene_Mean_EachTrial=MeanPSTH_FR(RandLocSceneSelector,:);
PSTH_RandLocScene_Mean=nanmean(PSTH_RandLocScene_Mean_EachTrial,1);
PSTH_RandLocScene_Sem=nanstd(PSTH_RandLocScene_Mean_EachTrial,[],1)./sqrt(sum(~isnan(PSTH_RandLocScene_Mean_EachTrial),1));

PSTH_DuelScene_Mean_EachTrial=MeanPSTH_FR(DuelLocSceneSelector,:);
PSTH_DuelScene_Mean=nanmean(PSTH_DuelScene_Mean_EachTrial,1);
PSTH_DuelScene_Sem=nanstd(PSTH_DuelScene_Mean_EachTrial,[],1)./sqrt(sum(~isnan(PSTH_DuelScene_Mean_EachTrial),1));




PSTH_FixLocScene_Scramble_Mean_EachTrial=MeanPSTH_FR(FixLocSceneScrambleSelector,:);
PSTH_FixLocScene_Scramble_Mean=nanmean(PSTH_FixLocScene_Scramble_Mean_EachTrial,1);
PSTH_FixLocScene_Scramble_Sem=nanstd(PSTH_FixLocScene_Scramble_Mean_EachTrial,[],1)./sqrt(sum(~isnan(PSTH_FixLocScene_Scramble_Mean_EachTrial),1));

PSTH_RandLocScene_Scramble_Mean_EachTrial=MeanPSTH_FR(RandLocSceneScrambleSelector,:);
PSTH_RandLocScene_Scramble_Mean=nanmean(PSTH_RandLocScene_Scramble_Mean_EachTrial,1);
PSTH_RandLocScene_Scramble_Sem=nanstd(PSTH_RandLocScene_Scramble_Mean_EachTrial,[],1)./sqrt(sum(~isnan(PSTH_RandLocScene_Scramble_Mean_EachTrial),1));

PSTH_DuelScene_Scramble_Mean_EachTrial=MeanPSTH_FR(DuelLocSceneScrambleSelector,:);
PSTH_DuelScene_Scramble_Mean=nanmean(PSTH_DuelScene_Scramble_Mean_EachTrial,1);
PSTH_DuelScene_Scramble_Sem=nanstd(PSTH_DuelScene_Scramble_Mean_EachTrial,[],1)./sqrt(sum(~isnan(PSTH_DuelScene_Scramble_Mean_EachTrial),1));






%Statistics between Fix and Rand
[h,p_FixRand]=ttest2(FixLocScene_Mean_EachTrial,RandLocScene_Mean_EachTrial);

[h,p_FixRand_Scram]=ttest2(FixLocScene_Scramble_Mean_EachTrial,RandLocScene_Scramble_Mean_EachTrial);
%Statistics between Fix and Duel
[h,p_FixDuel]=ttest2(FixLocScene_Mean_EachTrial,DuelScene_Mean_EachTrial);

[h,p_FixDuel_Scram]=ttest2(FixLocScene_Scramble_Mean_EachTrial,DuelScene_Scramble_Mean_EachTrial);

%Statistics between Rand and Duel
[h,p_RandDuel]=ttest2(RandLocScene_Mean_EachTrial,DuelScene_Mean_EachTrial);

[h,p_RandDuel_Scram]=ttest2(RandLocScene_Scramble_Mean_EachTrial,DuelScene_Scramble_Mean_EachTrial);


%Familiar v.s. Novel
FamiliarScene_Mean_EachTrial=MeanFRSceneOn(FamiliarSceneSelector);
FamiliarScene_Mean=nanmean(FamiliarScene_Mean_EachTrial);
FamiliarScene_Sem=nanstd(FamiliarScene_Mean_EachTrial)/sqrt(sum(FamiliarSceneSelector));

FamiliarScene_Scramble_Mean_EachTrial=MeanFRSceneOn(FamiliarSceneScrambleSelector);
FamiliarScene_Scramble_Mean=nanmean(FamiliarScene_Scramble_Mean_EachTrial);
FamiliarScene_Scramble_Sem=nanstd(FamiliarScene_Scramble_Mean_EachTrial)/sqrt(sum(FamiliarSceneScrambleSelector));


NovelScene_Mean_EachTrial=MeanFRSceneOn(NovelSceneSelector);
NovelScene_Mean=nanmean(NovelScene_Mean_EachTrial);
NovelScene_Sem=nanstd(NovelScene_Mean_EachTrial)/sqrt(sum(NovelSceneSelector));

NovelScene_Scramble_Mean_EachTrial=MeanFRSceneOn(NovelSceneScrambleSelector);
NovelScene_Scramble_Mean=nanmean(NovelScene_Scramble_Mean_EachTrial);
NovelScene_Scramble_Sem=nanstd(NovelScene_Scramble_Mean_EachTrial)/sqrt(sum(NovelSceneScrambleSelector));

%Statistics between Familiar and Novel 
[h,p_FamiliarNovel]=ttest2(FamiliarScene_Mean_EachTrial,NovelScene_Mean_EachTrial);

[h,p_FamiliarNovel_Scram]=ttest2(FamiliarScene_Scramble_Mean_EachTrial,NovelScene_Scramble_Mean_EachTrial);

%PSTH



PSTH_FamiliarScene_Mean_EachTrial=MeanPSTH_FR(FamiliarSceneSelector,:);
PSTH_FamiliarScene_Mean=nanmean(PSTH_FamiliarScene_Mean_EachTrial,1);
PSTH_FamiliarScene_Sem=nanstd(PSTH_FamiliarScene_Mean_EachTrial,[],1)./sqrt(sum(~isnan(PSTH_FamiliarScene_Mean_EachTrial),1));

PSTH_NovelScene_Mean_EachTrial=MeanPSTH_FR(NovelSceneSelector,:);
PSTH_NovelScene_Mean=nanmean(PSTH_NovelScene_Mean_EachTrial,1);
PSTH_NovelScene_Sem=nanstd(PSTH_NovelScene_Mean_EachTrial,[],1)./sqrt(sum(~isnan(PSTH_NovelScene_Mean_EachTrial),1));

PSTH_FamiliarScene_Scramble_Mean_EachTrial=MeanPSTH_FR(FamiliarSceneScrambleSelector,:);
PSTH_FamiliarScene_Scramble_Mean=nanmean(PSTH_FamiliarScene_Scramble_Mean_EachTrial,1);
PSTH_FamiliarScene_Scramble_Sem=nanstd(PSTH_FamiliarScene_Scramble_Mean_EachTrial,[],1)./sqrt(sum(~isnan(PSTH_FamiliarScene_Scramble_Mean_EachTrial),1));

PSTH_NovelScene_Scramble_Mean_EachTrial=MeanPSTH_FR(NovelSceneScrambleSelector,:);
PSTH_NovelScene_Scramble_Mean=nanmean(PSTH_NovelScene_Scramble_Mean_EachTrial,1);
PSTH_NovelScene_Scramble_Sem=nanstd(PSTH_NovelScene_Scramble_Mean_EachTrial,[],1)./sqrt(sum(~isnan(PSTH_NovelScene_Scramble_Mean_EachTrial),1));




NaturalScene_Mean_EachTrial=MeanFRSceneOn(NaturalSceneSelector);
NaturalScene_Mean=nanmean(NaturalScene_Mean_EachTrial);
NaturalScene_Sem=nanstd(NaturalScene_Mean_EachTrial)/sqrt(sum(NaturalSceneSelector));

NaturalScene_Scramble_Mean_EachTrial=MeanFRSceneOn(NaturalSceneScrambleSelector);
NaturalScene_Scramble_Mean=nanmean(NaturalScene_Scramble_Mean_EachTrial);
NaturalScene_Scramble_Sem=nanstd(NaturalScene_Scramble_Mean_EachTrial)/sqrt(sum(NaturalSceneScrambleSelector));


MapScene_Mean_EachTrial=MeanFRSceneOn(MapSceneSelector);
MapScene_Mean=nanmean(MapScene_Mean_EachTrial);
MapScene_Sem=nanstd(MapScene_Mean_EachTrial)/sqrt(sum(MapSceneSelector));

MapScene_Scramble_Mean_EachTrial=MeanFRSceneOn(MapSceneScrambleSelector);
MapScene_Scramble_Mean=nanmean(MapScene_Scramble_Mean_EachTrial);
MapScene_Scramble_Sem=nanstd(MapScene_Scramble_Mean_EachTrial)/sqrt(sum(MapSceneScrambleSelector));


%PSTH



PSTH_NaturalScene_Mean_EachTrial=MeanPSTH_FR(NaturalSceneSelector,:);
PSTH_NaturalScene_Mean=nanmean(PSTH_NaturalScene_Mean_EachTrial,1);
PSTH_NaturalScene_Sem=nanstd(PSTH_NaturalScene_Mean_EachTrial,[],1)./sqrt(sum(~isnan(PSTH_NaturalScene_Mean_EachTrial),1));

PSTH_MapScene_Mean_EachTrial=MeanPSTH_FR(MapSceneSelector,:);
PSTH_MapScene_Mean=nanmean(PSTH_MapScene_Mean_EachTrial,1);
PSTH_MapScene_Sem=nanstd(PSTH_MapScene_Mean_EachTrial,[],1)./sqrt(sum(~isnan(PSTH_MapScene_Mean_EachTrial),1));

PSTH_NaturalScene_Scramble_Mean_EachTrial=MeanPSTH_FR(NaturalSceneScrambleSelector,:);
PSTH_NaturalScene_Scramble_Mean=nanmean(PSTH_NaturalScene_Scramble_Mean_EachTrial,1);
PSTH_NaturalScene_Scramble_Sem=nanstd(PSTH_NaturalScene_Scramble_Mean_EachTrial,[],1)./sqrt(sum(~isnan(PSTH_NaturalScene_Scramble_Mean_EachTrial),1));

PSTH_MapScene_Scramble_Mean_EachTrial=MeanPSTH_FR(MapSceneScrambleSelector,:);
PSTH_MapScene_Scramble_Mean=nanmean(PSTH_MapScene_Scramble_Mean_EachTrial,1);
PSTH_MapScene_Scramble_Sem=nanstd(PSTH_MapScene_Scramble_Mean_EachTrial,[],1)./sqrt(sum(~isnan(PSTH_MapScene_Scramble_Mean_EachTrial),1));



%Statistics between Natural and Map Scene
[h,p_NaturalMap]=ttest2(NaturalScene_Mean_EachTrial,MapScene_Mean_EachTrial);

[h,p_NaturalMap_Scram]=ttest2(NaturalScene_Scramble_Mean_EachTrial,MapScene_Scramble_Mean_EachTrial);

%Find PSTH for each scene

%#################Separate trials according to normal scene name
%UniqueScene=unique(SceneIndex_Line);
UniqueScene=unique(SceneName_Line);
NumOfCondition=length(UniqueScene);

%Calculate and plot the average firing rate according to scenes

for i=1:NumOfCondition
    currentScene=UniqueScene(i);
    
  %  trialindex=SceneIndex_Line==currentScene;
  trialindex=SceneName_Line==currentScene;
  %{
    %For average Scene Tuning
    FR_Scene=MeanFRRaw(trialindex');
  
    FR_Scene_mean(i)=mean(FR_Scene);
    FR_Scene_sem(i)=std(FR_Scene)/sqrt(sum(trialindex));
    FR_Scene_SSE(i)=sum((FR_Scene-mean(FR_Scene)).^2); 
    
    %Measure the mean firing rate of the spon for each scene
    Spon_ForEachScene=Mean_Spon_ForEachTrial(trialindex);
    Spon_ForEachScene_Mean(i)=mean(Mean_Spon_ForEachTrial(trialindex));
    Spon_ForEachScene_SEM(i)=std(Mean_Spon_ForEachTrial(trialindex))/sqrt(sum(trialindex));
    
     FR_Spon_SSE(i)=sum((Spon_ForEachScene-mean(Spon_ForEachScene)).^2); 
    
    %p value from the spon
    [h_SceneFromSpon(i),p_SceneFromSpon(i)]=ttest2(currentScene,Spon_ForEachScene,'tail','right');
    %}
    
    
    
    
   %For scene-on  PSTH
    tmpmean=mean(MeanPSTH_FR(trialindex,:),1);
    tmpsem=std(MeanPSTH_FR(trialindex,:),[],1)./sqrt(sum(~isnan(MeanPSTH_FR(trialindex,:)),1));
    
    %Few number correction
    numCol=sum(~isnan(MeanPSTH_FR(trialindex,:)),1);
    removeCol=numCol<3;
    
    tmpmean(removeCol)=NaN;
    tmpsem(removeCol)=NaN;
    SpikeFiringRateSceneOn_mean(i,:)=tmpmean;
    SpikeFiringRateSceneOn_sem(i,:)=tmpsem;
   
    RelativeTimeMarkerForeachScene(i)=min(RelativeTimeMarker(trialindex,2));
    
    
    SceneStr{i}=sprintf('Scene%d',UniqueScene(i));
    
    %SSE dynamics
    FR_SSE_Scene_Time(i,:)=sum((MeanPSTH_FR(trialindex,:)-mean(MeanPSTH_FR(trialindex,:),1)).^2,1); 
   
    
    
    %Find out response latency, response duration, peak firing rate, peak
    %firing time;
    
    [p_EachScene(i,:),h_EachScene(i,:)]=ttest2_group(MeanPSTH_FR(trialindex,:),Mean_Spon_ForEachTrial,'right');
    
    p_EachScene(i,1:ceil(PlotStartOffset/StepSize))=NaN;
    h_EachScene(i,1:ceil(PlotStartOffset/StepSize))=NaN;
    h_EachScene(i,h_EachScene(i,:)==0)=NaN;
    
    %Response latency
   
   
    Latency=FindOutLatency(nanmean(MeanPSTH_Time),h_EachScene(i,:),ResponseLatencyThreshold);
    if ~isempty(Latency)
        ResponseLatencyForEachScene(i)=Latency;
        
    else
       ResponseLatencyForEachScene(i)=NaN; 
    end
   
    
    %Response Duration   
    ResponseDurationScene=nanmean(MeanPSTH_Time);
    
    ResponseDurationScene=ResponseDurationScene(h_EachScene(i,:)==1);
    ResponseDurationScene=ResponseDurationScene(ResponseDurationScene<=min(RelativeTimeMarker(:,2)) & ResponseDurationScene>=ResponseLatencyForEachScene(i));
     
    TotalResponseDurationScene(i)=length(ResponseDurationScene)*StepSize;
    TotalEventDuration=min(RelativeTimeMarker(:,2));
    PropResponseDurationScene(i)=TotalResponseDurationScene(i)/TotalEventDuration*100;
    
    %Peak Firing Rate in the first significant interval
    Significant_TimePoint=ResponseDurationScene(ResponseDurationScene>=ResponseLatencyForEachScene(i));
    FirstIntervalEnd=find(diff(Significant_TimePoint)>StepSize,1,'first');
    if isempty(FirstIntervalEnd) || FirstIntervalEnd>length(Significant_TimePoint)
        FirstIntervalEnd=length(Significant_TimePoint);
        
    end
    
    
    FirstSigInterval=[1:FirstIntervalEnd];
    
    FR_Sig=SpikeFiringRateSceneOn_mean(i,nanmean(MeanPSTH_Time)>=ResponseLatencyForEachScene(i));
    
    if ~isempty(FR_Sig) && ~isempty(FirstSigInterval)  
        
    [PeakFiringScene(i),loc]=max(FR_Sig(FirstSigInterval));
       
    PeakFiringTimeScene(i)=Significant_TimePoint(loc);
    
    %First Interval Duration
    FirstIntervalDurationScene(i)=diff(Significant_TimePoint([FirstSigInterval(1),FirstSigInterval(end)]));
    else
        PeakFiringScene(i)=NaN;
        PeakFiringTimeScene(i)=NaN;
        FirstIntervalDurationScene(i)=NaN;
            
    end
    
    
    %PeakFiringTimeScene(i)=ResponseDurationScene(find(ResponseDurationScene==PeakFiringScene(i),1,'first'));
    

    %SpikeFiringRateSceneOff_mean(i,:)=mean(SpikeFiringRateOff(trialindex,:),1);
    %SpikeFiringRateSceneOff_sem(i,:)=std(SpikeFiringRateOff(trialindex,:),[],1)/sqrt(sum(trialindex));
    
    

end
%Add spon
%SceneStr{i+1}='SponAll';

%Separate PSTH for each scene for each scene category
%Fix Loc Scene
UniqueFixLocSceneSel=ismember(UniqueScene,FixLocScene);
UniqueFixLocScene=UniqueScene(UniqueFixLocSceneSel);

PSTH_FixLocSceneSep_Mean=SpikeFiringRateSceneOn_mean(UniqueFixLocSceneSel,:);
PSTH_FixLocSceneSep_Sem=SpikeFiringRateSceneOn_sem(UniqueFixLocSceneSel,:);

%Rand Loc Scene
UniqueRandLocSceneSel=ismember(UniqueScene,RandLocScene);
UniqueRandLocScene=UniqueScene(UniqueRandLocSceneSel);

PSTH_RandLocSceneSep_Mean=SpikeFiringRateSceneOn_mean(UniqueRandLocSceneSel,:);
PSTH_RandLocSceneSep_Sem=SpikeFiringRateSceneOn_sem(UniqueRandLocSceneSel,:);

%Duel Scene
UniqueDuelSceneSel=ismember(UniqueScene,DuelLocScene);
UniqueDuelScene=UniqueScene(UniqueDuelSceneSel);

PSTH_DuelSceneSep_Mean=SpikeFiringRateSceneOn_mean(UniqueDuelSceneSel,:);
PSTH_DuelSceneSep_Sem=SpikeFiringRateSceneOn_sem(UniqueDuelSceneSel,:);

%Natural Scene
UniqueNaturalSceneSel=ismember(UniqueScene,[NaturalSceneInterval(1):NaturalSceneInterval(2)]);
UniqueNaturalScene=UniqueScene(UniqueNaturalSceneSel);

PSTH_NaturalSceneSep_Mean=SpikeFiringRateSceneOn_mean(UniqueNaturalSceneSel,:);
PSTH_NaturalSceneSep_Sem=SpikeFiringRateSceneOn_sem(UniqueNaturalSceneSel,:);

%Map Scene
UniqueMapSceneSel=ismember(UniqueScene,[MapSceneInterval(1):MapSceneInterval(2)]);
UniqueMapScene=UniqueScene(UniqueMapSceneSel);

PSTH_MapSceneSep_Mean=SpikeFiringRateSceneOn_mean(UniqueMapSceneSel,:);
PSTH_MapSceneSep_Sem=SpikeFiringRateSceneOn_sem(UniqueMapSceneSel,:);




%%
%If show figures 
if ShowFigureFlag
   
FigureStartNum=120;
    FigureIndex=1;
%%%%%%%%%%%%%%%%%%Figure 1: Barplot for comparison of different scene types

figtitlestr{FigureIndex}='Bar_plot_for_different_scene_types';
fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[50,100, 1200,600],'Name',figtitlestr{FigureIndex});


%1.Compare Normal Scene response with Scramble Scene and Grey Scenes
%axes('Position',[0.12,0.5,0.8,0.4])
subplot(2,2,1);
g=bar([NormalScene_Mean,ScrambleScene_Mean,GreyScene_Mean,PreSceneSpon_Mean],'LineWidth',3);



g.FaceColor = 'flat';
g.EdgeColor = 'flat';
ColorChosen= Color(1:4,:);

g.CData = ColorChosen;

xlim([0.2,4.8]);


XData=g.XData;
yData=g.YData;

XData=XData(1:end-1);
yData=yData(1:end-1);

%{
%Add significance level
text(XData(h_SceneFromSpon==1)-0.045,yData(h_SceneFromSpon==1)+max(FR_Scene_mean_Sorted+FR_Scene_sem(Index_Sort))*0.1,'*','FontSize',40)
text(XData(p_SceneFromSpon<0.01)-0.045,yData(p_SceneFromSpon<0.01)+max(FR_Scene_mean_Sorted+FR_Scene_sem(Index_Sort))*0.2,'*','FontSize',40)
text(XData(p_SceneFromSpon<0.001)-0.045,yData(p_SceneFromSpon<0.001)+max(FR_Scene_mean_Sorted+FR_Scene_sem(Index_Sort))*0.3,'*','FontSize',40)

ylim([0,max(FR_Scene_mean_Sorted+FR_Scene_sem(Index_Sort))*1.5]);
%}
hold on
errorbar([NormalScene_Mean,ScrambleScene_Mean,GreyScene_Mean,PreSceneSpon_Mean],[NormalScene_Sem,ScrambleScene_Sem,GreyScene_Sem,PreSceneSpon_Sem],'.k','LineWidth',3);

set(gca,'LineWidth',3)
xticks(1:4);
xticklabels({'Scene','Scramble','Grey','PreSceneSpon'});
set(gca,'FontSize',15,'FontWeight','Bold')

ylabel(['\fontsize{20}Firing Rate ',newline, '\fontsize{20}(Spikes/s)']);
box off
h=gca; h.XAxis.TickLength = [0 0];

%Statistics
yrange=ylim;
if p_SceneScram<0.05
ctr2 = bsxfun(@plus, g.XData, [g.XOffset]');
hold on
plot(ctr2(1:2), [1 1]*yrange(2)*1.35, '-k', 'LineWidth',3,'Color','r')
if p_SceneScram<0.05 & p_SceneScram>=0.01
   plot(mean(ctr2(1:2)), yrange(2)*1.4, '*k','MarkerSize',15,'Color','r','LineWidth',3);
elseif p_SceneScram<0.01 & p_SceneScram>=0.001
   plot(mean(ctr2(1:2))*0.9, yrange(2)*1.4, '*k','MarkerSize',15,'Color','r','LineWidth',3);
   plot(mean(ctr2(1:2))*1.1, yrange(2)*1.4, '*k','MarkerSize',15,'Color','r','LineWidth',3);
elseif p_SceneScram<0.001 
   plot(mean(ctr2(1:2))*0.9, yrange(2)*1.4, '*k','MarkerSize',15,'Color','r','LineWidth',3);
   plot(mean(ctr2(1:2))*1.1, yrange(2)*1.4, '*k','MarkerSize',15,'Color','r','LineWidth',3);
   plot(mean(ctr2(1:2))*1.3, yrange(2)*1.4, '*k','MarkerSize',15,'Color','r','LineWidth',3);
end

end

if p_SceneGrey<0.05
yrange=yrange*1.1;
ctr2 = bsxfun(@plus, g.XData, [g.XOffset]');
hold on
plot(ctr2([1,3]), [1 1]*yrange(2)*1.35, '-k', 'LineWidth',3,'Color','r')
if p_SceneGrey<0.05 & p_SceneGrey>=0.01
   plot(mean(ctr2([1,3])), yrange(2)*1.4, '*k','MarkerSize',15,'Color','r','LineWidth',3);
elseif p_SceneGrey<0.01 & p_SceneGrey>=0.001
   plot(mean(ctr2([1,3]))*0.9, yrange(2)*1.4, '*k','MarkerSize',15,'Color','r','LineWidth',3);
   plot(mean(ctr2([1,3]))*1.1, yrange(2)*1.4, '*k','MarkerSize',15,'Color','r','LineWidth',3);
elseif p_SceneGrey<0.001 
   plot(mean(ctr2([1,3]))*0.9, yrange(2)*1.4, '*k','MarkerSize',15,'Color','r','LineWidth',3);
   plot(mean(ctr2([1,3]))*1.1, yrange(2)*1.4, '*k','MarkerSize',15,'Color','r','LineWidth',3);
   plot(mean(ctr2([1,3]))*1.3, yrange(2)*1.4, '*k','MarkerSize',15,'Color','r','LineWidth',3);
end

end
%Comparison between scene and spon
if p_SceneSpon<0.05
yrange=yrange+1.5;
ctr2 = bsxfun(@plus, g.XData, [g.XOffset]');
hold on
plot(ctr2([1,4]), [1 1]*yrange(2)*1.35, '-k', 'LineWidth',3,'Color','r')
if p_SceneSpon<0.05 & p_SceneSpon>=0.01
   plot(mean(ctr2([1,4])), yrange(2)*1.4, '*k','MarkerSize',15,'Color','r','LineWidth',3);
elseif p_SceneSpon<0.01 & p_SceneSpon>=0.001
   plot(mean(ctr2([1,4]))*0.9, yrange(2)*1.4, '*k','MarkerSize',15,'Color','r','LineWidth',3);
   plot(mean(ctr2([1,4]))*1.1, yrange(2)*1.4, '*k','MarkerSize',15,'Color','r','LineWidth',3);
elseif p_SceneSpon<0.001 
   plot(mean(ctr2([1,4]))*0.9, yrange(2)*1.4, '*k','MarkerSize',15,'Color','r','LineWidth',3);
   plot(mean(ctr2([1,4]))*1.1, yrange(2)*1.4, '*k','MarkerSize',15,'Color','r','LineWidth',3);
   plot(mean(ctr2([1,4]))*1.3, yrange(2)*1.4, '*k','MarkerSize',15,'Color','r','LineWidth',3);
end

end

if p_ScrambleSpon<0.05
yrange=yrange+1.5;
ctr2 = bsxfun(@plus, g.XData, [g.XOffset]');
hold on
plot(ctr2([2,4]), [1 1]*yrange(2)*1.35, '-k', 'LineWidth',3,'Color','r')
if p_ScrambleSpon<0.05 & p_ScrambleSpon>=0.01
   plot(mean(ctr2([2,4])), yrange(2)*1.4, '*k','MarkerSize',15,'Color','r','LineWidth',3);
elseif p_ScrambleSpon<0.01 & p_ScrambleSpon>=0.001
   plot(mean(ctr2([2,4]))*0.9, yrange(2)*1.4, '*k','MarkerSize',15,'Color','r','LineWidth',3);
   plot(mean(ctr2([2,4]))*1.1, yrange(2)*1.4, '*k','MarkerSize',15,'Color','r','LineWidth',3);
elseif p_ScrambleSpon<0.001 
   plot(mean(ctr2([2,4]))*0.9, yrange(2)*1.4, '*k','MarkerSize',15,'Color','r','LineWidth',3);
   plot(mean(ctr2([2,4]))*1.1, yrange(2)*1.4, '*k','MarkerSize',15,'Color','r','LineWidth',3);
   plot(mean(ctr2([2,4]))*1.3, yrange(2)*1.4, '*k','MarkerSize',15,'Color','r','LineWidth',3);
end

end

if p_GreySpon<0.05
yrange=yrange+1.5;
ctr2 = bsxfun(@plus, g.XData, [g.XOffset]');
hold on
plot(ctr2([3,4]), [1 1]*yrange(2)*1.35, '-k', 'LineWidth',3,'Color','r')
if p_GreySpon<0.05 & p_GreySpon>=0.01
   plot(mean(ctr2([3,4])), yrange(2)*1.4, '*k','MarkerSize',15,'Color','r','LineWidth',3);
elseif p_GreySpon<0.01 & p_GreySpon>=0.001
   plot(mean(ctr2([3,4]))*0.9, yrange(2)*1.4, '*k','MarkerSize',15,'Color','r','LineWidth',3);
   plot(mean(ctr2([3,4]))*1.0, yrange(2)*1.4, '*k','MarkerSize',15,'Color','r','LineWidth',3);
elseif p_GreySpon<0.001 
   plot(mean(ctr2([3,4]))*0.9, yrange(2)*1.4, '*k','MarkerSize',15,'Color','r','LineWidth',3);
   plot(mean(ctr2([3,4]))*1.0, yrange(2)*1.4, '*k','MarkerSize',15,'Color','r','LineWidth',3);
   plot(mean(ctr2([3,4]))*1.1, yrange(2)*1.4, '*k','MarkerSize',15,'Color','r','LineWidth',3);
end

end




%2.Compare Familiar Scene, Novel Scene response with their scrambles
%axes('Position',[0.12,0.5,0.8,0.4])
subplot(2,2,2);


g=bar([FamiliarScene_Mean,FamiliarScene_Scramble_Mean,NovelScene_Mean,NovelScene_Scramble_Mean],'LineWidth',3);
g.FaceColor = 'flat';
g.EdgeColor = 'flat';
ColorChosen= Color(1:4,:);

g.CData = ColorChosen;

xlim([0.2,4.8]);


XData=g.XData;
yData=g.YData;

XData=XData(1:end-1);
yData=yData(1:end-1);


hold on
errorbar([FamiliarScene_Mean,FamiliarScene_Scramble_Mean,NovelScene_Mean,NovelScene_Scramble_Mean],[FamiliarScene_Sem,FamiliarScene_Scramble_Sem,NovelScene_Sem,NovelScene_Scramble_Sem],'.k','LineWidth',3);

set(gca,'LineWidth',3)
xticks(1:4);
xticklabels({'Familiar','FScram','Novel','NScram'});
set(gca,'FontSize',15,'FontWeight','Bold')

ylabel(['\fontsize{20}Firing Rate ',newline, '\fontsize{20}(Spikes/s)']);
box off
h=gca; h.XAxis.TickLength = [0 0];




%Statistics


yrange=ylim;
if p_FamiliarNovel<0.05
ctr2 = bsxfun(@plus, g.XData, [g.XOffset]');
hold on
plot(ctr2([1,3]), [1 1]*yrange(2)*1.35, '-k', 'LineWidth',3,'Color','r')
if p_FamiliarNovel<0.05 & p_FamiliarNovel>=0.01
   plot(mean(ctr2([1,3])), yrange(2)*1.4, '*k','MarkerSize',15,'Color','r','LineWidth',3);
elseif p_FamiliarNovel<0.01 & p_FamiliarNovel>=0.001
   plot(mean(ctr2([1,3]))*0.9, yrange(2)*1.4, '*k','MarkerSize',15,'Color','r','LineWidth',3);
   plot(mean(ctr2([1,3]))*1.1, yrange(2)*1.4, '*k','MarkerSize',15,'Color','r','LineWidth',3);
elseif p_FamiliarNovel<0.001 
   plot(mean(ctr2([1,3]))*0.9, yrange(2)*1.4, '*k','MarkerSize',15,'Color','r','LineWidth',3);
   plot(mean(ctr2([1,3]))*1.1, yrange(2)*1.4, '*k','MarkerSize',15,'Color','r','LineWidth',3);
   plot(mean(ctr2([1,3]))*1.3, yrange(2)*1.4, '*k','MarkerSize',15,'Color','r','LineWidth',3);
end

end




if p_FamiliarNovel_Scram<0.05
yrange=yrange*1.2;
ctr2 = bsxfun(@plus, g.XData, [g.XOffset]');
hold on
plot(ctr2([2,4]), [1 1]*yrange(2)*1.35, '-k', 'LineWidth',3,'Color','r')
if p_FamiliarNovel_Scram<0.05 & p_FamiliarNovel_Scram>=0.01
   plot(mean(ctr2([2,4])), yrange(2)*1.4, '*k','MarkerSize',15,'Color','r','LineWidth',3);
elseif p_FamiliarNovel_Scram<0.01 & p_FamiliarNovel_Scram>=0.001
   plot(mean(ctr2([2,4]))*0.9, yrange(2)*1.4, '*k','MarkerSize',15,'Color','r','LineWidth',3);
   plot(mean(ctr2([2,4]))*1.1, yrange(2)*1.4, '*k','MarkerSize',15,'Color','r','LineWidth',3);
elseif p_FamiliarNovel_Scram<0.001 
   plot(mean(ctr2([2,4]))*0.9, yrange(2)*1.4, '*k','MarkerSize',15,'Color','r','LineWidth',3);
   plot(mean(ctr2([2,4]))*1.1, yrange(2)*1.4, '*k','MarkerSize',15,'Color','r','LineWidth',3);
   plot(mean(ctr2([2,4]))*1.3, yrange(2)*1.4, '*k','MarkerSize',15,'Color','r','LineWidth',3);
end

end



%[h,p_FixRand]=ttest2(FixLocScene_Mean_EachTrial,RandLocScene_Mean_EachTrial);

%[h,p_FixRand_Scram]=ttest2(FixLocScene_Scramble_Mean_EachTrial,RandLocScene_Scramble_Mean_EachTrial);

%3.Compare FixLocScene, RandLoc Scene response with their scrambles
%axes('Position',[0.12,0.5,0.8,0.4])
subplot(2,2,3);


g=bar([FixLocScene_Mean,RandLocScene_Mean,DuelScene_Mean,FixLocScene_Scramble_Mean,RandLocScene_Scramble_Mean,DuelScene_Scramble_Mean],'LineWidth',3);
g.FaceColor = 'flat';
g.EdgeColor = 'flat';
ColorChosen= Color(1:6,:);

g.CData = ColorChosen;

xlim([0.2,6.8]);


XData=g.XData;
yData=g.YData;

XData=XData(1:end-1);
yData=yData(1:end-1);


hold on
errorbar([FixLocScene_Mean,RandLocScene_Mean,DuelScene_Mean,FixLocScene_Scramble_Mean,RandLocScene_Scramble_Mean,DuelScene_Scramble_Mean],[FixLocScene_Sem,RandLocScene_Sem,DuelScene_Sem,FixLocScene_Scramble_Sem,RandLocScene_Scramble_Sem,DuelScene_Scramble_Sem],'.k','LineWidth',3);

set(gca,'LineWidth',3)
xticks(1:6);
xticklabels({'FixL','RandL','Duel','FixLS','RandLS','DuelS'});
set(gca,'FontSize',15,'FontWeight','Bold')

ylabel(['\fontsize{20}Firing Rate ',newline, '\fontsize{20}(Spikes/s)']);
box off
h=gca; h.XAxis.TickLength = [0 0];



%Statistics

yrange=ylim;
if p_FixRand<0.05
ctr2 = bsxfun(@plus, g.XData, [g.XOffset]');
hold on
plot(ctr2([1,2]), [1 1]*yrange(2)*1.35, '-k', 'LineWidth',3,'Color','r')
if p_FixRand<0.05 & p_FixRand>=0.01
   plot(mean(ctr2([1,2])), yrange(2)*1.4, '*k','MarkerSize',15,'Color','r','LineWidth',3);
elseif p_FixRand<0.01 & p_FixRand>=0.001
   plot(mean(ctr2([1,2]))*0.9, yrange(2)*1.4, '*k','MarkerSize',15,'Color','r','LineWidth',3);
   plot(mean(ctr2([1,2]))*1.1, yrange(2)*1.4, '*k','MarkerSize',15,'Color','r','LineWidth',3);
elseif p_FixRand<0.001 
   plot(mean(ctr2([1,2]))*0.9, yrange(2)*1.4, '*k','MarkerSize',15,'Color','r','LineWidth',3);
   plot(mean(ctr2([1,2]))*1.1, yrange(2)*1.4, '*k','MarkerSize',15,'Color','r','LineWidth',3);
   plot(mean(ctr2([1,2]))*1.3, yrange(2)*1.4, '*k','MarkerSize',15,'Color','r','LineWidth',3);
end

end


if p_FixDuel<0.05
yrange=yrange*1.1;
ctr2 = bsxfun(@plus, g.XData, [g.XOffset]');
hold on
plot(ctr2([1,3]), [1 1]*yrange(2)*1.35, '-k', 'LineWidth',3,'Color','r')
if p_FixDuel<0.05 & p_FixDuel>=0.01
   plot(mean(ctr2([1,3])), yrange(2)*1.4, '*k','MarkerSize',15,'Color','r','LineWidth',3);
elseif p_FixDuel<0.01 & p_FixDuel>=0.001
   plot(mean(ctr2([1,3]))*0.9, yrange(2)*1.4, '*k','MarkerSize',15,'Color','r','LineWidth',3);
   plot(mean(ctr2([1,3]))*1.1, yrange(2)*1.4, '*k','MarkerSize',15,'Color','r','LineWidth',3);
elseif p_FixDuel<0.001 
   plot(mean(ctr2([1,3]))*0.9, yrange(2)*1.4, '*k','MarkerSize',15,'Color','r','LineWidth',3);
   plot(mean(ctr2([1,3]))*1.1, yrange(2)*1.4, '*k','MarkerSize',15,'Color','r','LineWidth',3);
   plot(mean(ctr2([1,3]))*1.3, yrange(2)*1.4, '*k','MarkerSize',15,'Color','r','LineWidth',3);
end

end

if p_RandDuel<0.05
yrange=yrange*1.1;
ctr2 = bsxfun(@plus, g.XData, [g.XOffset]');
hold on
plot(ctr2([2,3]), [1 1]*yrange(2)*1.35, '-k', 'LineWidth',3,'Color','r')
if p_RandDuel<0.05 & p_RandDuel>=0.01
   plot(mean(ctr2([2,3])), yrange(2)*1.4, '*k','MarkerSize',15,'Color','r','LineWidth',3);
elseif p_RandDuel<0.01 & p_RandDuel>=0.001
   plot(mean(ctr2([2,3]))*0.9, yrange(2)*1.4, '*k','MarkerSize',15,'Color','r','LineWidth',3);
   plot(mean(ctr2([2,3]))*1.1, yrange(2)*1.4, '*k','MarkerSize',15,'Color','r','LineWidth',3);
elseif p_RandDuel<0.001 
   plot(mean(ctr2([2,3]))*0.9, yrange(2)*1.4, '*k','MarkerSize',15,'Color','r','LineWidth',3);
   plot(mean(ctr2([2,3]))*1.1, yrange(2)*1.4, '*k','MarkerSize',15,'Color','r','LineWidth',3);
   plot(mean(ctr2([2,3]))*1.3, yrange(2)*1.4, '*k','MarkerSize',15,'Color','r','LineWidth',3);
end

end



if p_FixRand_Scram<0.05
yrange=yrange*1.1;
ctr2 = bsxfun(@plus, g.XData, [g.XOffset]');
hold on
plot(ctr2([4,5]), [1 1]*yrange(2)*1.35, '-k', 'LineWidth',3,'Color','r')
if p_FixRand_Scram<0.05 & p_FixRand_Scram>=0.01
   plot(mean(ctr2([4,5])), yrange(2)*1.4, '*k','MarkerSize',15,'Color','r','LineWidth',3);
elseif p_FixRand_Scram<0.01 & p_FixRand_Scram>=0.001
   plot(mean(ctr2([4,5]))*0.9, yrange(2)*1.4, '*k','MarkerSize',15,'Color','r','LineWidth',3);
   plot(mean(ctr2([4,5]))*1.1, yrange(2)*1.4, '*k','MarkerSize',15,'Color','r','LineWidth',3);
elseif p_FixRand_Scram<0.001 
   plot(mean(ctr2([4,5]))*0.9, yrange(2)*1.4, '*k','MarkerSize',15,'Color','r','LineWidth',3);
   plot(mean(ctr2([4,5]))*1.1, yrange(2)*1.4, '*k','MarkerSize',15,'Color','r','LineWidth',3);
   plot(mean(ctr2([4,5]))*1.3, yrange(2)*1.4, '*k','MarkerSize',15,'Color','r','LineWidth',3);
end

end

if p_FixDuel_Scram<0.05
yrange=yrange*1.1;
ctr2 = bsxfun(@plus, g.XData, [g.XOffset]');
hold on
plot(ctr2([4,6]), [1 1]*yrange(2)*1.35, '-k', 'LineWidth',3,'Color','r')
if p_FixDuel_Scram<0.05 & p_FixDuel_Scram>=0.01
   plot(mean(ctr2([4,6])), yrange(2)*1.4, '*k','MarkerSize',15,'Color','r','LineWidth',3);
elseif p_FixDuel_Scram<0.01 & p_FixDuel_Scram>=0.001
   plot(mean(ctr2([4,6]))*0.9, yrange(2)*1.4, '*k','MarkerSize',15,'Color','r','LineWidth',3);
   plot(mean(ctr2([4,6]))*1.1, yrange(2)*1.4, '*k','MarkerSize',15,'Color','r','LineWidth',3);
elseif p_FixDuel_Scram<0.001 
   plot(mean(ctr2([4,6]))*0.9, yrange(2)*1.4, '*k','MarkerSize',15,'Color','r','LineWidth',3);
   plot(mean(ctr2([4,6]))*1.1, yrange(2)*1.4, '*k','MarkerSize',15,'Color','r','LineWidth',3);
   plot(mean(ctr2([4,6]))*1.3, yrange(2)*1.4, '*k','MarkerSize',15,'Color','r','LineWidth',3);
end

end

if p_RandDuel_Scram<0.05
yrange=yrange*1.1;
ctr2 = bsxfun(@plus, g.XData, [g.XOffset]');
hold on
plot(ctr2([5,6]), [1 1]*yrange(2)*1.35, '-k', 'LineWidth',3,'Color','r')
if p_RandDuel_Scram<0.05 & p_RandDuel_Scram>=0.01
   plot(mean(ctr2([5,6])), yrange(2)*1.4, '*k','MarkerSize',15,'Color','r','LineWidth',3);
elseif p_RandDuel_Scram<0.01 & p_RandDuel_Scram>=0.001
   plot(mean(ctr2([5,6]))*0.9, yrange(2)*1.4, '*k','MarkerSize',15,'Color','r','LineWidth',3);
   plot(mean(ctr2([5,6]))*1.1, yrange(2)*1.4, '*k','MarkerSize',15,'Color','r','LineWidth',3);
elseif p_RandDuel_Scram<0.001 
   plot(mean(ctr2([5,6]))*0.9, yrange(2)*1.4, '*k','MarkerSize',15,'Color','r','LineWidth',3);
   plot(mean(ctr2([5,6]))*1.1, yrange(2)*1.4, '*k','MarkerSize',15,'Color','r','LineWidth',3);
   plot(mean(ctr2([5,6]))*1.3, yrange(2)*1.4, '*k','MarkerSize',15,'Color','r','LineWidth',3);
end

end

%{
[h,p_FixDuel]=ttest2(FixLocScene_Mean_EachTrial,DuelScene_Mean_EachTrial);

[h,p_FixDuel_Scram]=ttest2(FixLocScene_Scramble_Mean_EachTrial,DuelScene_Scramble_Mean_EachTrial);

%Statistics between Rand and Duel
[h,p_RandDuel]=ttest2(RandLocScene_Mean_EachTrial,DuelScene_Mean_EachTrial);

[h,p_RandDuel_Scram]=ttest2(RandLocScene_Scramble_Mean_EachTrial,DuelScene_Scramble_Mean_EachTrial);

%}


%4.Compare MapScene, Natural Scene response with their scrambles
%axes('Position',[0.12,0.5,0.8,0.4])
subplot(2,2,4);

g=bar([NaturalScene_Mean,NaturalScene_Scramble_Mean,MapScene_Mean,MapScene_Scramble_Mean],'LineWidth',3);
g.FaceColor = 'flat';
g.EdgeColor = 'flat';
ColorChosen= Color(1:4,:);

g.CData = ColorChosen;

xlim([0.2,4.8]);


XData=g.XData;
yData=g.YData;

XData=XData(1:end-1);
yData=yData(1:end-1);


hold on

errorbar([NaturalScene_Mean,NaturalScene_Scramble_Mean,MapScene_Mean,MapScene_Scramble_Mean],[NaturalScene_Sem,NaturalScene_Scramble_Sem,MapScene_Sem,MapScene_Scramble_Sem],'.k','LineWidth',3);

set(gca,'LineWidth',3)
xticks(1:4);
xticklabels({'Natural','NScram','Map','MScram'});
set(gca,'FontSize',15,'FontWeight','Bold')

ylabel(['\fontsize{20}Firing Rate ',newline, '\fontsize{20}(Spikes/s)']);
box off
h=gca; h.XAxis.TickLength = [0 0];

%Statistics


yrange=ylim;
if p_NaturalMap<0.05
ctr2 = bsxfun(@plus, g.XData, [g.XOffset]');
hold on
plot(ctr2([1,3]), [1 1]*yrange(2)*1.35, '-k', 'LineWidth',3,'Color','r')
if p_NaturalMap<0.05 & p_NaturalMap>=0.01
   plot(mean(ctr2([1,3])), yrange(2)*1.4, '*k','MarkerSize',15,'Color','r','LineWidth',3);
elseif p_NaturalMap<0.01 & p_NaturalMap>=0.001
   plot(mean(ctr2([1,3]))*0.9, yrange(2)*1.4, '*k','MarkerSize',15,'Color','r','LineWidth',3);
   plot(mean(ctr2([1,3]))*1.1, yrange(2)*1.4, '*k','MarkerSize',15,'Color','r','LineWidth',3);
elseif p_NaturalMap<0.001 
   plot(mean(ctr2([1,3]))*0.9, yrange(2)*1.4, '*k','MarkerSize',15,'Color','r','LineWidth',3);
   plot(mean(ctr2([1,3]))*1.1, yrange(2)*1.4, '*k','MarkerSize',15,'Color','r','LineWidth',3);
   plot(mean(ctr2([1,3]))*1.3, yrange(2)*1.4, '*k','MarkerSize',15,'Color','r','LineWidth',3);
end

end




if p_NaturalMap_Scram<0.05
yrange=yrange*1.5;
ctr2 = bsxfun(@plus, g.XData, [g.XOffset]');
hold on
plot(ctr2([2,4]), [1 1]*yrange(2)*1.35, '-k', 'LineWidth',3,'Color','r')
if p_NaturalMap_Scram<0.05 & p_NaturalMap_Scram>=0.01
   plot(mean(ctr2([2,4])), yrange(2)*1.4, '*k','MarkerSize',15,'Color','r','LineWidth',3);
elseif p_NaturalMap_Scram<0.01 & p_NaturalMap_Scram>=0.001
   plot(mean(ctr2([2,4]))*0.9, yrange(2)*1.4, '*k','MarkerSize',15,'Color','r','LineWidth',3);
   plot(mean(ctr2([2,4]))*1.1, yrange(2)*1.4, '*k','MarkerSize',15,'Color','r','LineWidth',3);
elseif p_NaturalMap_Scram<0.001 
   plot(mean(ctr2([2,4]))*0.9, yrange(2)*1.4, '*k','MarkerSize',15,'Color','r','LineWidth',3);
   plot(mean(ctr2([2,4]))*1.1, yrange(2)*1.4, '*k','MarkerSize',15,'Color','r','LineWidth',3);
   plot(mean(ctr2([2,4]))*1.3, yrange(2)*1.4, '*k','MarkerSize',15,'Color','r','LineWidth',3);
end

end
 data_msg=strcat(FileName,' SpikeChannel:',num2str(SpikeChannel),';ChannelID:',num2str(ChannelID),';ClusterID:',num2str(ClusterID),';TrialNum: ',num2str(TrialNum));

    annotation('textbox', [0.17, 0.95, 0.7, 0.04], 'string', data_msg,'FontSize',20,'FontWeight','Bold','LineStyle','none')
     

FigureStartNum=FigureStartNum+1;
 FigureIndex=FigureIndex+1;  

%%%%%%%%%%%%%%%%%%Figure 2: PSTH of different scene types 


figtitlestr{FigureIndex}='PSTH_between_Screen_Scramble_Grey';
fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[70,100, 1500,600],'Name',figtitlestr{FigureIndex});



subplot(2,4,1);
%Scene v.s. Scramble v.s. Grey
shadedErrorBar(PSTH_Time,PSTH_NormalScene_Mean,PSTH_NormalScene_Sem,'lineprops',{[1,0,0]},'transparent',1,'patchSaturation',0.3)
%plot(TimeSequence(i,:),SpikeFiringRateSceneOn_mean(i,:),'Color',Color(i,:))
hold on
shadedErrorBar(PSTH_Time,PSTH_ScrambleScene_Mean,PSTH_ScrambleScene_Sem,'lineprops',{[0,1,0]},'transparent',1,'patchSaturation',0.3)
shadedErrorBar(PSTH_Time,PSTH_GreyScene_Mean,PSTH_GreyScene_Sem,'lineprops',{[0.2,0.2,0.2]},'transparent',1,'patchSaturation',0.3)


yrange=ylim;
h1=plot([SceneOnMarker(1),SceneOnMarker(1)],[yrange(1),yrange(2)],'--k','LineWidth',1.5);
h2=plot([SceneOnMarker(2),SceneOnMarker(2)],[yrange(1),yrange(2)],'--k','LineWidth',1.5);

 set( get( get( h1, 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );
 set( get( get( h2, 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );

 set(findall(gca, 'Type', 'Line'),'LineWidth',3); 

xlabel('Time from scene onset');
ylabel('Firing Rate(Hz)');
legend('Scene','Scramble','Grey');
set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
xlim(plotrange);

subplot(2,4,2);
%Fix Loc v.s. Rand Loc
shadedErrorBar(PSTH_Time,PSTH_FixLocScene_Mean,PSTH_FixLocScene_Sem,'lineprops',{[1,0,0]},'transparent',1,'patchSaturation',0.3)
%plot(TimeSequence(i,:),SpikeFiringRateSceneOn_mean(i,:),'Color',Color(i,:))
hold on
shadedErrorBar(PSTH_Time,PSTH_RandLocScene_Mean,PSTH_RandLocScene_Sem,'lineprops',{[0,0,1]},'transparent',1,'patchSaturation',0.3)
shadedErrorBar(PSTH_Time,PSTH_DuelScene_Mean,PSTH_DuelScene_Sem,'lineprops',{[0,1,0]},'transparent',1,'patchSaturation',0.3)

yrange=ylim;
h1=plot([SceneOnMarker(1),SceneOnMarker(1)],[yrange(1),yrange(2)],'--k','LineWidth',1.5);
h2=plot([SceneOnMarker(2),SceneOnMarker(2)],[yrange(1),yrange(2)],'--k','LineWidth',1.5);

 set( get( get( h1, 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );
 set( get( get( h2, 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );
 

set(findall(gca, 'Type', 'Line'),'LineWidth',3); 
xlabel('Time from scene onset');
ylabel('Firing Rate(Hz)');
legend('FixLoc','RandLoc','Duel','FLScram','RLScram','DuelScram');
set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
xlim(plotrange);


subplot(2,4,3);
%Familiar v.s. Novel
shadedErrorBar(PSTH_Time,PSTH_FamiliarScene_Mean,PSTH_FamiliarScene_Sem,'lineprops',{[1,0,0]},'transparent',1,'patchSaturation',0.3)
%plot(TimeSequence(i,:),SpikeFiringRateSceneOn_mean(i,:),'Color',Color(i,:))
hold on
shadedErrorBar(PSTH_Time,PSTH_NovelScene_Mean,PSTH_NovelScene_Sem,'lineprops',{[0,0,1]},'transparent',1,'patchSaturation',0.3)

yrange=ylim;
h1=plot([SceneOnMarker(1),SceneOnMarker(1)],[yrange(1),yrange(2)],'--k','LineWidth',1.5);
h2=plot([SceneOnMarker(2),SceneOnMarker(2)],[yrange(1),yrange(2)],'--k','LineWidth',1.5);

 set( get( get( h1, 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );
 set( get( get( h2, 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );

set(findall(gca, 'Type', 'Line'),'LineWidth',3); 
xlabel('Time from scene onset');
ylabel('Firing Rate(Hz)');
legend('Familiar','Novel');
set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
xlim(plotrange);


subplot(2,4,4);
%Map v.s. Natural
shadedErrorBar(PSTH_Time,PSTH_NaturalScene_Mean,PSTH_NaturalScene_Sem,'lineprops',{[1,0,0]},'transparent',1,'patchSaturation',0.3);
%plot(TimeSequence(i,:),SpikeFiringRateSceneOn_mean(i,:),'Color',Color(i,:))
hold on
shadedErrorBar(PSTH_Time,PSTH_MapScene_Mean,PSTH_MapScene_Sem,'lineprops',{[0,0,1]},'transparent',1,'patchSaturation',0.3);

yrange=ylim;
h1=plot([SceneOnMarker(1),SceneOnMarker(1)],[yrange(1),yrange(2)],'--k','LineWidth',1.5);
h2=plot([SceneOnMarker(2),SceneOnMarker(2)],[yrange(1),yrange(2)],'--k','LineWidth',1.5);

 set( get( get( h1, 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );
 set( get( get( h2, 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );
 

set(findall(gca, 'Type', 'Line'),'LineWidth',3); 
xlabel('Time from scene onset');
ylabel('Firing Rate(Hz)');
legend('Natural','Map');
set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');

xlim(plotrange);



subplot(2,4,6);
%Fix Loc v.s. Rand Loc(Scramble)
shadedErrorBar(PSTH_Time,PSTH_FixLocScene_Scramble_Mean,PSTH_FixLocScene_Scramble_Sem,'lineprops',{[0.5,0,0]},'transparent',1,'patchSaturation',0.3);
hold on
shadedErrorBar(PSTH_Time,PSTH_RandLocScene_Scramble_Mean,PSTH_RandLocScene_Scramble_Sem,'lineprops',{[0,0,0.5]},'transparent',1,'patchSaturation',0.3);
shadedErrorBar(PSTH_Time,PSTH_DuelScene_Scramble_Mean,PSTH_DuelScene_Scramble_Sem,'lineprops',{[0,0.5,0]},'transparent',1,'patchSaturation',0.3);

yrange=ylim;
h1=plot([SceneOnMarker(1),SceneOnMarker(1)],[yrange(1),yrange(2)],'--k','LineWidth',1.5);
h2=plot([SceneOnMarker(2),SceneOnMarker(2)],[yrange(1),yrange(2)],'--k','LineWidth',1.5);

 set( get( get( h1, 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );
 set( get( get( h2, 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );
 
set(findall(gca, 'Type', 'Line'),'LineWidth',3); 
xlabel('Time from scene onset');
ylabel('Firing Rate(Hz)');
legend('FLScram','RLScram','DuelScram');
set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
xlim(plotrange);

subplot(2,4,7);

%Familiar v.s. Novel(Scramble)


shadedErrorBar(PSTH_Time,PSTH_FamiliarScene_Scramble_Mean,PSTH_FamiliarScene_Scramble_Sem,'lineprops',{[0.5,0,0]},'transparent',1,'patchSaturation',0.3)
hold on
shadedErrorBar(PSTH_Time,PSTH_NovelScene_Scramble_Mean,PSTH_NovelScene_Scramble_Sem,'lineprops',{[0,0,0.5]},'transparent',1,'patchSaturation',0.3)

yrange=ylim;
h1=plot([SceneOnMarker(1),SceneOnMarker(1)],[yrange(1),yrange(2)],'--k','LineWidth',1.5);
h2=plot([SceneOnMarker(2),SceneOnMarker(2)],[yrange(1),yrange(2)],'--k','LineWidth',1.5);

 set( get( get( h1, 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );
 set( get( get( h2, 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );
 
set(findall(gca, 'Type', 'Line'),'LineWidth',3); 
xlabel('Time from scene onset');
ylabel('Firing Rate(Hz)');
legend('FScram','NScram');
set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
xlim(plotrange);


subplot(2,4,8);
%Map v.s. Natural(Scramble)
shadedErrorBar(PSTH_Time,PSTH_NaturalScene_Scramble_Mean,PSTH_NaturalScene_Scramble_Sem,'lineprops',{[0.5,0,0]},'transparent',1,'patchSaturation',0.3);
hold on
shadedErrorBar(PSTH_Time,PSTH_MapScene_Scramble_Mean,PSTH_MapScene_Scramble_Sem,'lineprops',{[0,0,0.5]},'transparent',1,'patchSaturation',0.3);

yrange=ylim;
h1=plot([SceneOnMarker(1),SceneOnMarker(1)],[yrange(1),yrange(2)],'--k','LineWidth',1.5);
h2=plot([SceneOnMarker(2),SceneOnMarker(2)],[yrange(1),yrange(2)],'--k','LineWidth',1.5);

 set( get( get( h1, 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );
 set( get( get( h2, 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );
 
set(findall(gca, 'Type', 'Line'),'LineWidth',3); 
xlabel('Time from scene onset');
ylabel('Firing Rate(Hz)');
legend('NScram','MScram');
set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');

xlim(plotrange);
 data_msg=strcat(FileName,' SpikeChannel:',num2str(SpikeChannel),';ChannelID:',num2str(ChannelID),';ClusterID:',num2str(ClusterID),';TrialNum: ',num2str(TrialNum));

    annotation('textbox', [0.17, 0.95, 0.7, 0.04], 'string', data_msg,'FontSize',20,'FontWeight','Bold','LineStyle','none')
     


FigureStartNum=FigureStartNum+1;
 FigureIndex=FigureIndex+1;  
 
 
 %%%%%%%%%%%%%%%%%%Figure 3: PSTH of each scene in each type


figtitlestr{FigureIndex}='PSTH_InsideEachCategory';
fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[90,100, 1500,600],'Name',figtitlestr{FigureIndex});

%plotrange=[-PlotStartOffset(1),1500];
if isempty(UniqueDuelScene)
    totalCatNum=4;
else
    totalCatNum=5;
end

subplot(2,ceil(totalCatNum/2),1);
%Fix Location Scenes
for i=1:length(UniqueFixLocScene)
    shadedErrorBar(PSTH_Time,PSTH_FixLocSceneSep_Mean(i,:),PSTH_FixLocSceneSep_Sem(i,:),'lineprops',{[ColorForScene(i,:)]},'transparent',1,'patchSaturation',0.3);
    hold on
    
end
set(findall(gca, 'Type', 'Line'),'LineWidth',3); 
yrange=ylim;
h1=plot([SceneOnMarker(1),SceneOnMarker(1)],[yrange(1),yrange(2)],'--k','LineWidth',1.5);
hold on
h2=plot([SceneOnMarker(2),SceneOnMarker(2)],[yrange(1),yrange(2)],'--k','LineWidth',1.5);

 set( get( get( h1, 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );
 set( get( get( h2, 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );
 
xlabel('Time from scene onset');
ylabel('Firing Rate(Hz)');
title('FixLoc Scenes');
if legendOn
legend(num2str(UniqueFixLocScene));
end
set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
xlim(plotrange);



subplot(2,ceil(totalCatNum/2),2);
%Rand Location Scenes
for i=1:length(UniqueRandLocScene)
    shadedErrorBar(PSTH_Time,PSTH_RandLocSceneSep_Mean(i,:),PSTH_RandLocSceneSep_Sem(i,:),'lineprops',{[ColorForScene(i,:)]},'transparent',1,'patchSaturation',0.3);
    hold on
    
end
set(findall(gca, 'Type', 'Line'),'LineWidth',3); 
yrange=ylim;
h1=plot([SceneOnMarker(1),SceneOnMarker(1)],[yrange(1),yrange(2)],'--k','LineWidth',1.5);
hold on
h2=plot([SceneOnMarker(2),SceneOnMarker(2)],[yrange(1),yrange(2)],'--k','LineWidth',1.5);

 set( get( get( h1, 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );
 hold on
 set( get( get( h2, 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );
 
xlabel('Time from scene onset');
ylabel('Firing Rate(Hz)');
title('RandLoc Scenes');
if legendOn
legend(num2str(UniqueRandLocScene));
end
set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
xlim(plotrange);


nextlabel=3;   

if totalCatNum==5
    subplot(2,ceil(totalCatNum/2),3);
%Duel Location Scenes
for i=1:length(UniqueDuelScene)
    shadedErrorBar(PSTH_Time,PSTH_DuelSceneSep_Mean(i,:),PSTH_DuelSceneSep_Sem(i,:),'lineprops',{[ColorForScene(i,:)]},'transparent',1,'patchSaturation',0.3);
    hold on
    
end
yrange=ylim;
h1=plot([SceneOnMarker(1),SceneOnMarker(1)],[yrange(1),yrange(2)],'--k','LineWidth',1.5);
h2=plot([SceneOnMarker(2),SceneOnMarker(2)],[yrange(1),yrange(2)],'--k','LineWidth',1.5);

 set( get( get( h1, 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );
 set( get( get( h2, 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );
 
set(findall(gca, 'Type', 'Line'),'LineWidth',3); 
xlabel('Time from scene onset');
ylabel('Firing Rate(Hz)');
title('Duel Scenes');
if legendOn
legend(num2str(UniqueDuelScene));
end
set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
xlim(plotrange);
   
 nextlabel=4;   
end


 subplot(2,ceil(totalCatNum/2),nextlabel);
%Natural Location Scenes
for i=1:length(UniqueNaturalScene)
    shadedErrorBar(PSTH_Time,PSTH_NaturalSceneSep_Mean(i,:),PSTH_NaturalSceneSep_Sem(i,:),'lineprops',{[ColorForScene(i,:)]},'transparent',1,'patchSaturation',0.3);
    hold on
    
end
set(findall(gca, 'Type', 'Line'),'LineWidth',3); 

yrange=ylim;
h1=plot([SceneOnMarker(1),SceneOnMarker(1)],[yrange(1),yrange(2)],'--k','LineWidth',1.5);
h2=plot([SceneOnMarker(2),SceneOnMarker(2)],[yrange(1),yrange(2)],'--k','LineWidth',1.5);

 set( get( get( h1, 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );
 set( get( get( h2, 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );
 
xlabel('Time from scene onset');
ylabel('Firing Rate(Hz)');
title('Natural Scenes');
if legendOn
legend(num2str(UniqueNaturalScene));
end
set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
xlim(plotrange);

nextlabel=nextlabel+1;


subplot(2,ceil(totalCatNum/2),nextlabel);
%Map Location Scenes
for i=1:length(UniqueMapScene)
    shadedErrorBar(PSTH_Time,PSTH_MapSceneSep_Mean(i,:),PSTH_MapSceneSep_Sem(i,:),'lineprops',{[ColorForScene(i,:)]},'transparent',1,'patchSaturation',0.3);
    hold on
    
end
set(findall(gca, 'Type', 'Line'),'LineWidth',3); 

yrange=ylim;
h1=plot([SceneOnMarker(1),SceneOnMarker(1)],[yrange(1),yrange(2)],'--k','LineWidth',1.5);
h2=plot([SceneOnMarker(2),SceneOnMarker(2)],[yrange(1),yrange(2)],'--k','LineWidth',1.5);

 set( get( get( h1, 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );
 set( get( get( h2, 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );
 
xlabel('Time from scene onset');
ylabel('Firing Rate(Hz)');
title('Map Scenes');
if legendOn
legend(num2str(UniqueMapScene));
end
set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
xlim(plotrange);

 data_msg=strcat(FileName,' SpikeChannel:',num2str(SpikeChannel),';ChannelID:',num2str(ChannelID),';ClusterID:',num2str(ClusterID),';TrialNum: ',num2str(TrialNum));

    annotation('textbox', [0.17, 0.95, 0.7, 0.04], 'string', data_msg,'FontSize',20,'FontWeight','Bold','LineStyle','none')

 
 
 FigureStartNum=FigureStartNum+1;
 FigureIndex=FigureIndex+1;  
 
 
 
 %%%%%%%%%%%%%%%%%%Figure 4: same as figure 2, but focus on the scene off
 %%%%%%%%%%%%%%%%%%and reward period


figtitlestr{FigureIndex}='PSTHSceneOff';
fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[110,100, 1500,600],'Name',figtitlestr{FigureIndex});



subplot(2,4,1);
%Scene v.s. Scramble v.s. Grey
shadedErrorBar(PSTH_Time,PSTH_NormalScene_Mean,PSTH_NormalScene_Sem,'lineprops',{[1,0,0]},'transparent',1,'patchSaturation',0.3)
%plot(TimeSequence(i,:),SpikeFiringRateSceneOn_mean(i,:),'Color',Color(i,:))
hold on
shadedErrorBar(PSTH_Time,PSTH_ScrambleScene_Mean,PSTH_ScrambleScene_Sem,'lineprops',{[0,1,0]},'transparent',1,'patchSaturation',0.3)
shadedErrorBar(PSTH_Time,PSTH_GreyScene_Mean,PSTH_GreyScene_Sem,'lineprops',{[0.2,0.2,0.2]},'transparent',1,'patchSaturation',0.3)


yrange=ylim;
h1=plot([SceneOnMarker(1),SceneOnMarker(1)],[yrange(1),yrange(2)],'--k','LineWidth',1.5);
h2=plot([SceneOnMarker(2),SceneOnMarker(2)],[yrange(1),yrange(2)],'--k','LineWidth',1.5);

h3=plot([RwdTimeRelMean,RwdTimeRelMean],[yrange(1),yrange(2)],'--b','LineWidth',1.5);




 set( get( get( h1, 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );
 set( get( get( h2, 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );
 set( get( get( h3, 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );

 set(findall(gca, 'Type', 'Line'),'LineWidth',3); 

xlabel('Time from scene onset');
ylabel('Firing Rate(Hz)');
legend('Scene','Scramble','Grey');
set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');

xlim(plotrangeoff);

subplot(2,4,2);
%Fix Loc v.s. Rand Loc
shadedErrorBar(PSTH_Time,PSTH_FixLocScene_Mean,PSTH_FixLocScene_Sem,'lineprops',{[1,0,0]},'transparent',1,'patchSaturation',0.3)
%plot(TimeSequence(i,:),SpikeFiringRateSceneOn_mean(i,:),'Color',Color(i,:))
hold on
shadedErrorBar(PSTH_Time,PSTH_RandLocScene_Mean,PSTH_RandLocScene_Sem,'lineprops',{[0,0,1]},'transparent',1,'patchSaturation',0.3)
shadedErrorBar(PSTH_Time,PSTH_DuelScene_Mean,PSTH_DuelScene_Sem,'lineprops',{[0,1,0]},'transparent',1,'patchSaturation',0.3)

yrange=ylim;
h1=plot([SceneOnMarker(1),SceneOnMarker(1)],[yrange(1),yrange(2)],'--k','LineWidth',1.5);
h2=plot([SceneOnMarker(2),SceneOnMarker(2)],[yrange(1),yrange(2)],'--k','LineWidth',1.5);
h3=plot([RwdTimeRelMean,RwdTimeRelMean],[yrange(1),yrange(2)],'--b','LineWidth',1.5);

 set( get( get( h1, 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );
 set( get( get( h2, 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );
 set( get( get( h3, 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );
 

set(findall(gca, 'Type', 'Line'),'LineWidth',3); 
xlabel('Time from scene onset');
ylabel('Firing Rate(Hz)');
legend('FixLoc','RandLoc','Duel','FLScram','RLScram','DuelScram');
set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
xlim(plotrangeoff);


subplot(2,4,3);
%Familiar v.s. Novel
shadedErrorBar(PSTH_Time,PSTH_FamiliarScene_Mean,PSTH_FamiliarScene_Sem,'lineprops',{[1,0,0]},'transparent',1,'patchSaturation',0.3)
%plot(TimeSequence(i,:),SpikeFiringRateSceneOn_mean(i,:),'Color',Color(i,:))
hold on
shadedErrorBar(PSTH_Time,PSTH_NovelScene_Mean,PSTH_NovelScene_Sem,'lineprops',{[0,0,1]},'transparent',1,'patchSaturation',0.3)

yrange=ylim;
h1=plot([SceneOnMarker(1),SceneOnMarker(1)],[yrange(1),yrange(2)],'--k','LineWidth',1.5);
h2=plot([SceneOnMarker(2),SceneOnMarker(2)],[yrange(1),yrange(2)],'--k','LineWidth',1.5);
h3=plot([RwdTimeRelMean,RwdTimeRelMean],[yrange(1),yrange(2)],'--b','LineWidth',1.5);

 set( get( get( h1, 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );
 set( get( get( h2, 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );
  set( get( get( h3, 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );
 

set(findall(gca, 'Type', 'Line'),'LineWidth',3); 
xlabel('Time from scene onset');
ylabel('Firing Rate(Hz)');
legend('Familiar','Novel');
set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
xlim(plotrangeoff);


subplot(2,4,4);
%Map v.s. Natural
shadedErrorBar(PSTH_Time,PSTH_NaturalScene_Mean,PSTH_NaturalScene_Sem,'lineprops',{[1,0,0]},'transparent',1,'patchSaturation',0.3);
%plot(TimeSequence(i,:),SpikeFiringRateSceneOn_mean(i,:),'Color',Color(i,:))
hold on
shadedErrorBar(PSTH_Time,PSTH_MapScene_Mean,PSTH_MapScene_Sem,'lineprops',{[0,0,1]},'transparent',1,'patchSaturation',0.3);

yrange=ylim;
h1=plot([SceneOnMarker(1),SceneOnMarker(1)],[yrange(1),yrange(2)],'--k','LineWidth',1.5);
h2=plot([SceneOnMarker(2),SceneOnMarker(2)],[yrange(1),yrange(2)],'--k','LineWidth',1.5);
h3=plot([RwdTimeRelMean,RwdTimeRelMean],[yrange(1),yrange(2)],'--b','LineWidth',1.5);


 set( get( get( h1, 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );
 set( get( get( h2, 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );
  set( get( get( h3, 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );
 
 

set(findall(gca, 'Type', 'Line'),'LineWidth',3); 
xlabel('Time from scene onset');
ylabel('Firing Rate(Hz)');
legend('Natural','Map');
set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');

xlim(plotrangeoff);



subplot(2,4,6);
%Fix Loc v.s. Rand Loc(Scramble)
shadedErrorBar(PSTH_Time,PSTH_FixLocScene_Scramble_Mean,PSTH_FixLocScene_Scramble_Sem,'lineprops',{[0.5,0,0]},'transparent',1,'patchSaturation',0.3);
hold on
shadedErrorBar(PSTH_Time,PSTH_RandLocScene_Scramble_Mean,PSTH_RandLocScene_Scramble_Sem,'lineprops',{[0,0,0.5]},'transparent',1,'patchSaturation',0.3);
shadedErrorBar(PSTH_Time,PSTH_DuelScene_Scramble_Mean,PSTH_DuelScene_Scramble_Sem,'lineprops',{[0,0.5,0]},'transparent',1,'patchSaturation',0.3);

yrange=ylim;
h1=plot([SceneOnMarker(1),SceneOnMarker(1)],[yrange(1),yrange(2)],'--k','LineWidth',1.5);
h2=plot([SceneOnMarker(2),SceneOnMarker(2)],[yrange(1),yrange(2)],'--k','LineWidth',1.5);
h3=plot([RwdTimeRelMean,RwdTimeRelMean],[yrange(1),yrange(2)],'--b','LineWidth',1.5);


 set( get( get( h1, 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );
 set( get( get( h2, 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );
  set( get( get( h3, 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );
 

 
set(findall(gca, 'Type', 'Line'),'LineWidth',3); 
xlabel('Time from scene onset');
ylabel('Firing Rate(Hz)');
legend('FLScram','RLScram','DuelScram');
set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
xlim(plotrangeoff);

subplot(2,4,7);

%Familiar v.s. Novel(Scramble)


shadedErrorBar(PSTH_Time,PSTH_FamiliarScene_Scramble_Mean,PSTH_FamiliarScene_Scramble_Sem,'lineprops',{[0.5,0,0]},'transparent',1,'patchSaturation',0.3)
hold on
shadedErrorBar(PSTH_Time,PSTH_NovelScene_Scramble_Mean,PSTH_NovelScene_Scramble_Sem,'lineprops',{[0,0,0.5]},'transparent',1,'patchSaturation',0.3)

yrange=ylim;
h1=plot([SceneOnMarker(1),SceneOnMarker(1)],[yrange(1),yrange(2)],'--k','LineWidth',1.5);
h2=plot([SceneOnMarker(2),SceneOnMarker(2)],[yrange(1),yrange(2)],'--k','LineWidth',1.5);
h3=plot([RwdTimeRelMean,RwdTimeRelMean],[yrange(1),yrange(2)],'--b','LineWidth',1.5);

 set( get( get( h1, 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );
 set( get( get( h2, 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );
  set( get( get( h3, 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );
 
 
set(findall(gca, 'Type', 'Line'),'LineWidth',3); 
xlabel('Time from scene onset');
ylabel('Firing Rate(Hz)');
legend('FScram','NScram');
set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
xlim(plotrangeoff);


subplot(2,4,8);
%Map v.s. Natural(Scramble)
shadedErrorBar(PSTH_Time,PSTH_NaturalScene_Scramble_Mean,PSTH_NaturalScene_Scramble_Sem,'lineprops',{[0.5,0,0]},'transparent',1,'patchSaturation',0.3);
hold on
shadedErrorBar(PSTH_Time,PSTH_MapScene_Scramble_Mean,PSTH_MapScene_Scramble_Sem,'lineprops',{[0,0,0.5]},'transparent',1,'patchSaturation',0.3);

yrange=ylim;
h1=plot([SceneOnMarker(1),SceneOnMarker(1)],[yrange(1),yrange(2)],'--k','LineWidth',1.5);
h2=plot([SceneOnMarker(2),SceneOnMarker(2)],[yrange(1),yrange(2)],'--k','LineWidth',1.5);
h3=plot([RwdTimeRelMean,RwdTimeRelMean],[yrange(1),yrange(2)],'--b','LineWidth',1.5);

 set( get( get( h1, 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );
 set( get( get( h2, 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );
  set( get( get( h3, 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );
 
set(findall(gca, 'Type', 'Line'),'LineWidth',3); 
xlabel('Time from scene onset');
ylabel('Firing Rate(Hz)');
legend('NScram','MScram');
set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');

xlim(plotrangeoff);
 data_msg=strcat(FileName,' SpikeChannel:',num2str(SpikeChannel),';ChannelID:',num2str(ChannelID),';ClusterID:',num2str(ClusterID),';TrialNum: ',num2str(TrialNum));

    annotation('textbox', [0.17, 0.95, 0.7, 0.04], 'string', data_msg,'FontSize',20,'FontWeight','Bold','LineStyle','none')
     


FigureStartNum=FigureStartNum+1;
 FigureIndex=FigureIndex+1;  
 
 
 %RwdTimeRelMean
 
%{ 
 %Separate PSTH for each scene for each scene category
%Fix Loc Scene
UniqueFixLocSceneSel=ismember(UniqueScene,FixLocScene);
UniqueFixLocScene=UniqueScene(UniqueFixLocSceneSel);

PSTH_FixLocSceneSep_Mean=SpikeFiringRateSceneOn_mean(UniqueFixLocSceneSel,:);
PSTH_FixLocSceneSep_Sem=SpikeFiringRateSceneOn_sem(UniqueFixLocSceneSel,:);

%Rand Loc Scene
UniqueRandLocSceneSel=ismember(UniqueScene,RandLocScene);
UniqueRandLocScene=UniqueScene(UniqueRandLocSceneSel);

PSTH_RandLocSceneSep_Mean=SpikeFiringRateSceneOn_mean(UniqueRandLocSceneSel,:);
PSTH_RandLocSceneSep_Sem=SpikeFiringRateSceneOn_sem(UniqueRandLocSceneSel,:);

%Duel Scene
UniqueDuelSceneSel=ismember(UniqueScene,DuelLocScene);
UniqueDuelScene=UniqueScene(UniqueDuelSceneSel);

PSTH_DuelSceneSep_Mean=SpikeFiringRateSceneOn_mean(UniqueDuelSceneSel,:);
PSTH_DuelSceneSep_Sem=SpikeFiringRateSceneOn_sem(UniqueDuelSceneSel,:);

%Natural Scene
UniqueNaturalSceneSel=ismember(UniqueScene,[NaturalSceneInterval(1):NaturalSceneInterval(2)]);
UniqueNaturalScene=UniqueScene(UniqueNaturalSceneSel);

%Map Scene
UniqueMapSceneSel=ismember(UniqueScene,[MapSceneInterval(1):MapSceneInterval(2)]);
UniqueMapScene=UniqueScene(UniqueMapSceneSel);
%}
%{
%%
%%%%%%%%%%%%%%%%%%Figure 7: PSTH of different scene types for each scene 
figtitlestr{7}='PSTH_between_Screen_Stable_FlexibleEach';

fig{7}=PrepareFigure(105,'w',[90,100, 1000,600],'Name',figtitlestr{7});
%{
PSTH_SceneOnPureContextScene=[SpikeFiringRateSceneOn_mean(ismember(UniqueScene,FixLocScene),:);...
    SpikeFiringRateSceneOn_mean(ismember(UniqueScene,RandLocScene),:);];

PSTH_SceneOnPureContextSceneScramble=[SpikeFiringRateSceneOn_mean(ismember(UniqueScene,-FixLocScene),:);...
    SpikeFiringRateSceneOn_mean(ismember(UniqueScene,-RandLocScene),:);];
%}
FullContextScene=[1000:1007,1010:1017];
UniqueContextScenes=[UniqueScene(ismember(UniqueScene,RandLocScene));UniqueScene(ismember(UniqueScene,FixLocScene))];
 %First half
 NumOfContextScenes=size(PSTH_SceneOnPureContextScene,1);
 
 SceneOnMarker=nanmean(RelativeTimeMarker,1);
for i=1:ceil(NumOfContextScenes/2)
    
    currScenePlot=UniqueContextScenes(i);
   
  subplot(2,2,1);
  shadedErrorBar(PSTH_Time,PSTH_SceneOnPureContextScene(i,:),PSTH_SceneOnPureContextScene_Sem(i,:),'lineprops',{ColorForScene(find(currScenePlot==FullContextScene),:)},'transparent',1,'patchSaturation',0.5);
  hold on
  %
  %Add time marker
  if i==ceil(length(NumOfContextScenes)/2)
      hold on
 % plot([SceneOnMarker(1), SceneOnMarker(1)],[0,max(max(PSTH_SceneOnPureContextScene+PSTH_SceneOnPureContextScene_Sem))],'-k','Color','#707276','LineWidth',1.5);
 % plot([SceneOnMarker(2), SceneOnMarker(2)],[0,max(max(PSTH_SceneOnPureContextScene+PSTH_SceneOnPureContextScene_Sem))],'-k','Color','#707276','LineWidth',1.5);
 
   h1=plot([SceneOnMarker(1), SceneOnMarker(1)],[0,max(max(PSTH_SceneOnPureContextScene+PSTH_SceneOnPureContextScene_Sem))],'--k','Color','#707276','LineWidth',1.5);
  h2=plot([SceneOnMarker(2), SceneOnMarker(2)],[0,max(max(PSTH_SceneOnPureContextScene+PSTH_SceneOnPureContextScene_Sem))],'--k','Color','#707276','LineWidth',1.5);
  
  set( get( get( h1, 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );
 set( get( get( h2, 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );
 
  end
  
  if i==ceil(NumOfContextScenes/2)
      %Add the Grey Scene
       shadedErrorBar(PSTH_Time,PSTH_GreyScene_Mean,PSTH_GreyScene_Sem,'lineprops',{'#94969A'},'transparent',1,'patchSaturation',0.5);

      
  end
     

  

  set(findall(gca, 'Type', 'Line'),'LineWidth',3);
  xlabel('Time from scene onset(ms)');
  ylabel('Firing Rate(Hz)');
  

  title('Full Context SceneOnset');
  set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
   xlim([-500,500+SceneOnMarker(2)]);
   
  subplot(2,2,2);
  if ~isempty(PSTH_SceneOnPureContextSceneScramble)
   shadedErrorBar(PSTH_Time,PSTH_SceneOnPureContextSceneScramble(i,:),PSTH_SceneOnPureContextSceneScramble_Sem(i,:),'lineprops',{ColorForScene(find(currScenePlot==FullContextScene),:)},'transparent',1,'patchSaturation',0.5);
  hold on
   %Add time marker
  h1=plot([SceneOnMarker(1), SceneOnMarker(1)],[0,max(max(PSTH_SceneOnPureContextScene+PSTH_SceneOnPureContextScene_Sem))],'--k','Color','#707276','LineWidth',1.5);
  h2=plot([SceneOnMarker(2), SceneOnMarker(2)],[0,max(max(PSTH_SceneOnPureContextScene+PSTH_SceneOnPureContextScene_Sem))],'--k','Color','#707276','LineWidth',1.5);
  
  set( get( get( h1, 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );
 set( get( get( h2, 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );
  
  
    xlim([-500,500+SceneOnMarker(2)]);
    hold on
    if i==NumOfContextScenes/2
      %Add the Grey Scene
       shadedErrorBar(PSTH_Time,PSTH_GreyScene_Mean,PSTH_GreyScene_Sem,'lineprops',{'#94969A'},'transparent',1,'patchSaturation',0.5);

      
    end
     
   set(findall(gca, 'Type', 'Line'),'LineWidth',3);
    xlabel('Time from scene onset(ms)');
  ylabel('Firing Rate(Hz)');

   title('Scramble Scene Onset');
   set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
     hold on
  end
    
end
ld=legend([string(UniqueContextScenes(1:ceil(length(UniqueContextScenes)/2)))]);
%ld=legend([string(UniqueScene);'SSI']);
ld.Box='off';
ld.Position=[0.9   0.7    0.0854    0.0880];

%Second half
for i=ceil(NumOfContextScenes/2)+1:NumOfContextScenes
   currScenePlot=UniqueContextScenes(i);
   
  subplot(2,2,3);
  shadedErrorBar(PSTH_Time,PSTH_SceneOnPureContextScene(i,:),PSTH_SceneOnPureContextScene_Sem(i,:),'lineprops',{ColorForScene(find(currScenePlot==FullContextScene),:)},'transparent',1,'patchSaturation',0.5);
  hold on
  %
  %Add time marker
  if i==NumOfContextScenes
      hold on
  %Add time marker
  h1=plot([SceneOnMarker(1), SceneOnMarker(1)],[0,max(max(PSTH_SceneOnPureContextScene+PSTH_SceneOnPureContextScene_Sem))],'--k','Color','#707276','LineWidth',1.5);
  h2=plot([SceneOnMarker(2), SceneOnMarker(2)],[0,max(max(PSTH_SceneOnPureContextScene+PSTH_SceneOnPureContextScene_Sem))],'--k','Color','#707276','LineWidth',1.5);
  
   
  set( get( get( h1, 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );
 set( get( get( h2, 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );
  end
  
  if i==NumOfContextScenes
      %Add the Grey Scene
       shadedErrorBar(PSTH_Time,PSTH_GreyScene_Mean,PSTH_GreyScene_Sem,'lineprops',{'#94969A'},'transparent',1,'patchSaturation',0.5);

      
  end
  xlim([-500,500+SceneOnMarker(2)]);
     
  set(findall(gca, 'Type', 'Line'),'LineWidth',3);
  xlabel('Time from scene onset(ms)');
  ylabel('Firing Rate(Hz)');
 % xlim([-PlotSceneOnMarkerStart,PlotSceneOnMarkerStop+SceneOnMarker(2)]);

  title('Full Context Scene');
  set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
     
     
  subplot(2,2,4);
  if ~isempty(PSTH_SceneOnPureContextSceneScramble)
    shadedErrorBar(PSTH_Time,PSTH_SceneOnPureContextSceneScramble(i,:),PSTH_SceneOnPureContextSceneScramble_Sem(i,:),'lineprops',{ColorForScene(find(currScenePlot==FullContextScene),:)},'transparent',1,'patchSaturation',0.5);
%{
    
   %Add time marker
  if i==length(UniqueScenes)
    plot([SceneOffMarker(1), SceneOffMarker(1)],[0,max(cellfun(@(x,y) max(x+y),SceneOff_Resp_Mean_Each,SceneOff_Resp_Sem_Each))],'--k','LineWidth',1.5);
  plot([SceneOffMarker(2), SceneOffMarker(2)],[0,max(cellfun(@(x,y) max(x+y),SceneOff_Resp_Mean_Each,SceneOff_Resp_Sem_Each))],'--k','LineWidth',1.5);
  end
    %}
    if i==NumOfContextScenes
      hold on
  %Add time marker
  h1=plot([SceneOnMarker(1), SceneOnMarker(1)],[0,max(max(PSTH_SceneOnPureContextScene+PSTH_SceneOnPureContextScene_Sem))],'--k','Color','#707276','LineWidth',1.5);
  h2=plot([SceneOnMarker(2), SceneOnMarker(2)],[0,max(max(PSTH_SceneOnPureContextScene+PSTH_SceneOnPureContextScene_Sem))],'--k','Color','#707276','LineWidth',1.5);
  
   
  set( get( get( h1, 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );
 set( get( get( h2, 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );
    end
  
    hold on
    if i==NumOfContextScenes
      %Add the Grey Scene
       shadedErrorBar(PSTH_Time,PSTH_GreyScene_Mean,PSTH_GreyScene_Sem,'lineprops',{'#94969A'},'transparent',1,'patchSaturation',0.5);

      
    end
  xlim([-500,500+SceneOnMarker(2)]);
   set(findall(gca, 'Type', 'Line'),'LineWidth',3);
    xlabel('Time from scene onset(ms)');
  ylabel('Firing Rate(Hz)');
  title('Scramble Scene Onset');
%  xlim([-PlotSceneOffMarkerStart,PlotSceneOffMarkerStop+SceneOffMarker(2)]);
 %  title('SceneOffset');
   set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
     hold on
  end

end
  ld=legend([string(UniqueContextScenes(ceil(length(UniqueContextScenes)/2)+1:end));'GreyScene']);
%ld=legend([string(UniqueScene);'SSI']);
ld.Box='off';
ld.Position=[0.92   0.5    0.0854    0.0880];
%ld=legend([string(UniqueScene);'SSI']);
 %Add data information
  %  data_msg=sprintf(' %s, SpikeChannel: %d, TrialNum: %d ',FileName,SpikeChannel,TrialNum);
%    data_msg=strcat(FileName,' SpikeChannel:',num2str(SpikeChannel),';ChannelID:',num2str(ChannelID),';ClusterID:',num2str(ClusterID),';TrialNum: ',num2str(TrialNum),celltype);

 %   annotation('textbox', [0.17, 0.95, 0.7, 0.04], 'string', data_msg,'FontSize',20,'FontWeight','Bold','LineStyle','none')
    

%%
figtitlestr{8}='PSTH_between_Screen_Stable_FlexibleEach';

fig{8}=PrepareFigure(106,'w',[110,100, 1000,600],'Name',figtitlestr{8});

%{
pSceneDiffWithTime_Stable(RemovalFlagOn)=NaN;
rsquare_psth_Stable(RemovalFlagOn)=NaN;

h_sig_ssi_withtime_Stable=floor(pSceneDiffWithTime_Stable<0.05);
h_sig_ssi_withtime_Stable(h_sig_ssi_withtime_Stable==0)=NaN;
%}

%Context Scene
subplot(2,1,1)
plot(PSTH_Time,rsquare_psth_Stable,'-k','Color',ColorForEachMode(AssociateBlock==11,:),'LineWidth',3);
hold on
plot(PSTH_Time,rsquare_psth_Flexible,'-k','Color',ColorForEachMode(AssociateBlock==1,:),'LineWidth',3);

title('Scene Selectivity Dynamics for stable scenes and flexible scenes');
xlabel('Time from scene onset');
ylabel('SSI');
set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');

 xlim([-500,500+SceneOnMarker(2)]);

%Add significance
sig=area(PSTH_Time,h_sig_ssi_withtime_Stable*0.01,'FaceColor',ColorForEachMode(AssociateBlock==11,:),'EdgeColor','none','FaceAlpha',0.5,'EdgeAlpha',0.5);
%Add significance
sig=area(PSTH_Time,-h_sig_ssi_withtime_Flexible*0.01,'FaceColor',ColorForEachMode(AssociateBlock==1,:),'EdgeColor','none','FaceAlpha',0.5,'EdgeAlpha',0.5);
box off




title('Scene Selectivity Dynamics for stable and flexible scenes');
xlabel('Time from scene onset');
ylabel('SSI');
set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
box off

subplot(2,1,2)
plot(PSTH_Time,rsquare_psth_StableScramble,'-k','Color',ColorForEachMode(AssociateBlock==11,:),'LineWidth',3);
hold on
plot(PSTH_Time,rsquare_psth_FlexibleScramble,'-k','Color',ColorForEachMode(AssociateBlock==1,:),'LineWidth',3);
xlim([-500,500+SceneOnMarker(2)]);

title('Scene Selectivity Dynamics for scramble scenes');
xlabel('Time from scene onset');
ylabel('SSI');
set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');

%Add significance
sig=area(PSTH_Time,h_sig_ssi_withtime_StableScramble*0.01,'FaceColor',ColorForEachMode(AssociateBlock==11,:),'EdgeColor','none','FaceAlpha',0.5,'EdgeAlpha',0.5);
%Add significance
sig=area(PSTH_Time,-h_sig_ssi_withtime_FlexibleScramble*0.01,'FaceColor',ColorForEachMode(AssociateBlock==1,:),'EdgeColor','none','FaceAlpha',0.5,'EdgeAlpha',0.5);
box off

%data_msg=strcat(FileName,' SpikeChannel:',num2str(SpikeChannel),';ChannelID:',num2str(ChannelID),';ClusterID:',num2str(ClusterID),';TrialNum: ',num2str(TrialNum),celltype);

% annotation('textbox', [0.1, 0.96, 0.8, 0.04], 'string', data_msg,'FontSize',15,'FontWeight','Bold','LineStyle','none');
 

 %%
 if ~isempty(ContextScene)
 figtitlestr{9}='PearsonCorrelationBetweenScramNorm';
 fig{9}=PrepareFigure(107,'w',[120,100, 1000,600],'Name',figtitlestr{9});
 MarkerSize=150;
 scatter(FRSceneMean_ContextStableForComparison,FRSceneMean_ScrambleStableForComparison,MarkerSize*ones(1,size(FRSceneMean_ContextStableForComparison,1)),'o','MarkerFaceColor',ColorForEachMode(AssociateBlock==11,:),'MarkerEdgeColor','w');
 hold on
  scatter(FRSceneMean_ContextFlexibleForComparison,FRSceneMean_ScrambleFlexibleForComparison,MarkerSize*ones(1,size(FRSceneMean_ContextFlexibleForComparison,1)),'o','MarkerFaceColor',ColorForEachMode(AssociateBlock==1,:),'MarkerEdgeColor','w');
 axis equal;
 %Add diagnal line
 x=min([FRSceneMean_ContextStableForComparison,FRSceneMean_ScrambleStableForComparison,FRSceneMean_ContextFlexibleForComparison,FRSceneMean_ScrambleFlexibleForComparison])-0.1:max([FRSceneMean_ContextStableForComparison,FRSceneMean_ScrambleStableForComparison,FRSceneMean_ContextFlexibleForComparison,FRSceneMean_ScrambleFlexibleForComparison]+1);
 
 plot(x,x,'--k','LineWidth',2);
 xlim([min(x),max(x)]);
 ylim([min(x),max(x)])
 
 set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
box off;
xlabel('Over spon response of stable&flexible scene');
ylabel('Over spon response of scramble scene')
 
 %FRSceneMean_Stable,FRSceneMean_Flexible,FRSceneMean_StableScramble,FRSceneMean_FlexibleScramble
 end
%} 
end %end of if show pictures


%Output to files
%Setup output directory
Workingdirectory=pwd;
MarkerFolder='DataAnalysis';
%MarkerFolder='DataHub';
Flag=strfind(Workingdirectory,MarkerFolder);
BasicDirectory=Workingdirectory(1:Flag+length(MarkerFolder));

%{
%Raw responses for output 

FRSceneMean_Stable=FR_Scene_mean(ismember(UniqueScene,FixLocScene))-Spon_ForEachScene_Mean(ismember(UniqueScene,FixLocScene));
FRSceneMean_Flexible=FR_Scene_mean(ismember(UniqueScene,RandLocScene))-Spon_ForEachScene_Mean(ismember(UniqueScene,RandLocScene));

FRSceneMean_StableScramble=FR_Scene_mean(ismember(UniqueScene,StabelSceneScramble))-Spon_ForEachScene_Mean(ismember(UniqueScene,StabelSceneScramble));
FRSceneMean_FlexibleScramble=FR_Scene_mean(ismember(UniqueScene,FlexibleSceneScramble))-Spon_ForEachScene_Mean(ismember(UniqueScene,FlexibleSceneScramble));

[mat_ContextVSScramble]=corrcoef([FRSceneMean_Stable,FRSceneMean_Flexible],[FRSceneMean_StableScramble,FRSceneMean_FlexibleScramble]);
p_ContextVSScramble=mat_ContextVSScramble(1,2);
%}

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

%Store figures in common folders
ProtocolMarker='ScenePassiveView';
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
    if isfield(OutputFileName,'SceneTuningPassive')
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

OutputData_New.SceneTuningPassive.TrialType='GoodTrials';

%%%%%%%%%Save overall scene response characteristics%%%%%%%%%%%%%%%%%%%%%%%
%Organize the data
NameStr={'TrialNum','StabelScene','FlexibleScene','FRSceneMean_Stable','FRSceneMean_Flexible','FRSceneMean_StableScramble','FRSceneMean_FlexibleScramble','p_ContextVSScramble','r_ContextVSScramble',...
'R2_Whole','R2_Scramble','R2_Stable','R2_Flexible','p_SceneOnAnova','p_SceneOnAnovaScramble','p_SceneOnAnovaStable','p_SceneOnAnovaFlexible',...
'PropOfSig_Stable','PropOfSig_Flexible'
};

DataLib={TrialNum,OutputStableScenes,OutputFlexibleScenes,FRSceneMean_Stable,FRSceneMean_Flexible,FRSceneMean_StableScramble,FRSceneMean_FlexibleScramble,p_ContextVSScramble,r_ContextVSScramble,...
  R2_Whole, R2_Scramble,R2_Stable,R2_Flexible,p_SceneOnAnova, p_SceneOnAnovaScramble,p_SceneOnAnovaStable,p_SceneOnAnovaFlexible,...
  PropOfSig_Stable,PropOfSig_Flexible

};


DataStamp=containers.Map(NameStr,DataLib);





OutputData_New.SceneTuningPassive.Task=TaskType;
OutputData_New.SceneTuningPassiveg.TaskCode=TaskCode;
OutputData_New.SceneTuningPassive.DataStamp=DataStamp;
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





if ExistFlag
%Compare the old one with the new one
Task_Old=OutputData.SceneTuningPassive.TaskCode;

     if Task_Old~=TaskCode
        %Not the same task, add the new dataset into the old one
        OutputData(numel(OutputData)+1)=OutputData_New;
     else
         %If the same task, replace the old one with the new one
         OutputData=OutputData_New; 
    
    
     end
else
    %Output the current file
   OutputData=OutputData_New; 
end



if OutputFlag
    cd(OutPath);
    save(OutputFileName,'OutputData');
    disp('Data Exported');
else
    disp('No data export');
end  
%}








end



end