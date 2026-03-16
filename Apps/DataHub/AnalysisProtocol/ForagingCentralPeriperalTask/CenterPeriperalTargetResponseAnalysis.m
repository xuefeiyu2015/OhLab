function  CenterPeriperalTargetResponseAnalysis(Data, StartTrial, EndTrial,ShowFigureFlag,OutputFlag);
%Prototol to analyse the response three kinds of periperal
%target(Good/Bad/Cue) and the response difference between central gaze and
%periperal gaze
%Xuefei Yu 11012020


ECode;


SegEvent=[TGTCD,CUEON ];%Event Flow of interest
%SegEventEnd=[SHOWFIXCD,FIXOFF,TGTOFF,RWDOFFCD ,NaN,EXTRAFP];%Event end
SegEventEnd=[TGTOFF,TGTCD];%Event end
SegEventStr={'TargetOn','CueOn'};

AverageInterval_Central=[0,300];%Interval for polar plot
AverageInterval_Peripheral=[0,300];%Interval for polar plot

PlotStartOffset=[500,500];%in ms
PlotStopOffset=[100,100];%in ms

ResponseLatencyThreshold=3;

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
ColorForPeripheralMode={'#EA1A41','#648FFF','#36AF64'};
ColorForCentralMode={'#EA1A41','#648FFF'};
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
%Abnormal check %Especially For Robin081120FORAGEXPLO4.rst
CheckCode=SCENE_ON;%This code should appear once in each trial
AbnormalTrials=AbnormalTrialCheck(Data,CheckCode);

Selection=Selection(~ismember(Selection,AbnormalTrials));
TrialNum=length(Selection);

%FixationPointIndex=Data.FixPointIndex(Selection,:);
%FixLoc=Data.FixLoc(Selection,:,:);


TargetIndex=Data.TargetIndex(Selection,2:end,1);
TargetName=Data.TargetIndex(Selection,2:end,2);
%TargetLoc=Data.TargetLoc(Selection,:,:);

TargetNum=Data.TargetNum(Selection,2:end);

SceneIndex=Data.SceneIndex(Selection,2);
SceneName=Data.SceneIndex(Selection,3);
TargetAngle=Data.TargetPolarLoc(Selection,2:end,1);





%Event channel

EventChannel=Data.EventChannel(Selection,:);
EventTimeChannel=Data.EventTimeChannel(Selection,:);
EventBin=Data.EventBin(Selection,:);

%Get the cueonindex data
ParaLib=Data.ParaLib;

if(ismember('TgtCueOnFlag',keys(ParaLib)))
CueOnFlagInfo=ParaLib('TgtCueOnFlag');
CueOnFlagInfo=ReOrganizeData(CueOnFlagInfo);

CueOnFlag=CueOnFlagInfo(Selection,2:end);




else
    disp('TgtCueOnFlag information, change the analysis program');
    return 
    

end






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

%Skip the eye calibration for the central periperal task
%{
%Recalibrate eyes according to the extra FP
try
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
    
FP_Loc=Data.FixLoc(Selection,2:end,:);
Ref_Loc=squeeze(FP_Loc(:,1,:));
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

catch
    disp('Offline Eye calibration Failed');
end
%}
Eye_Eccentricity=cellfun(@(x,y) sqrt(x.^2+y.^2),EyeChannel_X_L,EyeChannel_Y_L,'uniform',0);



%EventBin=Data.EventBin(Selection,:);

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


TimeMarker=SegmentData.RelativeTimeMarker; %[Mean,sem]%Averall mean and sem

%PSTH for each event
MeanPSTH=SegmentData.MeanPSTH;%=[PSTH_Event_Time',PSTH_Event'];
MeanPSTH_Time=MeanPSTH(:,1);
MeanPSTH_FR=MeanPSTH(:,2);

%{
%Get the average firing rate 200ms before fix off as the baseline
SponEvent=[FIXOFF];
PreTimeBuffer=200;
SponEventEnd=NaN;
StartOffset=0;
StopOffset=PreTimeBuffer;

EventTimeChannelForSpon=EventTimeChannel-PreTimeBuffer;


SegmentDataSpon=EventSegTool(SponEvent,SponEventEnd,StartOffset,StopOffset,EventChannel,EventTimeChannelForSpon,SpikeTime,StepSize,BinWidth);

Mean_Spon=SegmentDataSpon.MeanFR;
Mean_Spon_ForEachTrial=SegmentDataSpon.MeanFRRaw{1};
%}
%Not consider about the eye temporally 
%{
%Get eye trance distance during the interval
StartTimeBin=FindOutTime(EventChannel,EventTimeChannel,SegEvent)-PlotStartOffset;
EndTimeBin=FindOutTime(EventChannel,EventTimeChannel,SegEventEnd )+PlotStopOffset;
  
Eye_Eccentricity_Tgt=SelectEyeInterval(Eye_Eccentricity,StartTimeBin,EndTimeBin)';
Eye_Eccentricity_Tgt=ReorganizeEye(Eye_Eccentricity_Tgt);
%}




%Rearrange target into lines
TargetIndex_Line=TargetIndex;

TargetIndex_Line=reshape(TargetIndex_Line,numel(TargetIndex_Line),1);

NoNNaNSelector=~isnan(TargetIndex_Line);

TargetIndex_Line=TargetIndex_Line(~isnan(TargetIndex_Line));


TargetName_Line=TargetName;

TargetName_Line=reshape(TargetName_Line,numel(TargetName_Line),1);

TargetName_Line=TargetName_Line(~isnan(TargetName_Line));

GoodObjectIndex=mod(TargetIndex_Line,2)==1;

%Rearrange cueonindex into lines
CueOnFlag_Line=CueOnFlag;
CueOnFlag_Line=reshape(CueOnFlag,numel(CueOnFlag),1);
CueOnFlag_Line=CueOnFlag_Line(NoNNaNSelector);

%Rearrange scenename into lines
SceneName_Line=repmat(SceneName,1,size(TargetName,2));
SceneName_Line=reshape(SceneName_Line,numel(SceneName_Line),1);
SceneName_Line=SceneName_Line(NoNNaNSelector);

%Rearrange target location into lines
TargetAngle_Line=TargetAngle;
TargetAngle_Line=reshape(TargetAngle_Line,numel(TargetAngle_Line),1);

TargetAngle_Line=TargetAngle_Line(~isnan(TargetAngle));


%{
%Select eye trace
Eye_Eccentricity_Tgt=Eye_Eccentricity_Tgt(NoNNaNSelector,:);

%Analysis Good versus Bad
GoodObjectIndex=mod(TargetIndex_Line,2)==1;

%Eye eccentricity Good v.s. Bad
Eye_Eccentricity_Tgt_Good=Eye_Eccentricity_Tgt(GoodObjectIndex,:);
Eye_Eccentricity_Tgt_Bad=Eye_Eccentricity_Tgt(~GoodObjectIndex,:);


Eye_Ecc_Good_mean=nanmean(Eye_Eccentricity_Tgt_Good);
Eye_Ecc_Bad_mean=nanmean(Eye_Eccentricity_Tgt_Bad);

Eye_Ecc_Good_sem=nanstd(Eye_Eccentricity_Tgt_Good)./sqrt(sum(~isnan(Eye_Eccentricity_Tgt_Good),1));
Eye_Ecc_Bad_sem=nanstd(Eye_Eccentricity_Tgt_Bad)./sqrt(sum(~isnan(Eye_Eccentricity_Tgt_Bad),1));

Eye_Ecc_All_mean=nanmean(Eye_Eccentricity_Tgt);
Eye_Ecc_All_sem=nanstd(Eye_Eccentricity_Tgt)./sqrt(sum(~isnan(Eye_Eccentricity_Tgt),1));
%}

%All conditions
GoodPeripheralCondition=CueOnFlag_Line==0 & GoodObjectIndex;
BadPeripheralCondition=CueOnFlag_Line==0 & ~GoodObjectIndex;

GoodCentralCondition= GoodObjectIndex & CueOnFlag_Line==1;
BadCentralCondition= ~GoodObjectIndex & CueOnFlag_Line==1;


%PSTH for all these conditions
PSTH_GoodPeripheral=MeanPSTH_FR{1}(GoodPeripheralCondition,:);
PSTH_BadPeripheral=MeanPSTH_FR{1}(BadPeripheralCondition,:);
PSTH_CuePeripheral=MeanPSTH_FR{2};

PSTH_Time_GoodPeripheral=nanmean(MeanPSTH_Time{1}(GoodPeripheralCondition,:),1);
PSTH_Time_BadPeripheral=nanmean(MeanPSTH_Time{1}(BadPeripheralCondition,:),1);
PSTH_Time_CuePeripheral=nanmean(MeanPSTH_Time{2},1);


PSTH_GoodCentral=MeanPSTH_FR{1}(GoodCentralCondition,:);
PSTH_BadCentral=MeanPSTH_FR{1}(BadCentralCondition,:);

PSTH_Time_GoodCentral=nanmean(MeanPSTH_Time{1}(GoodCentralCondition,:),1);
PSTH_Time_BadCentral=nanmean(MeanPSTH_Time{1}(BadCentralCondition,:),1);



%PSTH Average for good and bad and cue for central and periperal
PSTH_Mean_GoodPeripheral=nanmean(PSTH_GoodPeripheral,1);
PSTH_Mean_BadPeripheral=nanmean(PSTH_BadPeripheral,1);
PSTH_Mean_CuePeripheral=nanmean(PSTH_CuePeripheral,1);

PSTH_Sem_GoodPeripheral=nanstd(PSTH_GoodPeripheral,[],1)./sqrt(sum(~isnan(PSTH_GoodPeripheral),1));
PSTH_Sem_BadPeripheral=nanstd(PSTH_BadPeripheral,[],1)./sqrt(sum(~isnan(PSTH_BadPeripheral),1));
PSTH_Sem_CuePeripheral=nanstd(PSTH_CuePeripheral,[],1)./sqrt(sum(~isnan(PSTH_CuePeripheral),1));

PSTH_Mean_GoodCentral=nanmean(PSTH_GoodCentral,1);
PSTH_Mean_BadCentral=nanmean(PSTH_BadCentral,1);


PSTH_Sem_GoodCentral=nanstd(PSTH_GoodCentral,[],1)./sqrt(sum(~isnan(PSTH_GoodCentral),1));
PSTH_Sem_BadCentral=nanstd(PSTH_BadCentral,[],1)./sqrt(sum(~isnan(PSTH_BadCentral),1));


%Select interval for an average for both good, bad, cue, central and peripheral

Interval_Sel_Central=  (PSTH_Time_GoodCentral>=AverageInterval_Central(1) & PSTH_Time_GoodCentral<=AverageInterval_Central(2)) ;
Interval_Sel_Peripheral=  (PSTH_Time_GoodPeripheral>=AverageInterval_Peripheral(1) & PSTH_Time_GoodCentral<=AverageInterval_Peripheral(2)) ;

PSTH_Mean_GoodPeripheral_Interval=nanmean(PSTH_GoodPeripheral(:,Interval_Sel_Peripheral),2);
PSTH_Mean_BadPeripheral_Interval=nanmean(PSTH_BadPeripheral(:,Interval_Sel_Peripheral),2);
PSTH_Mean_CuePeripheral_Interval=nanmean(PSTH_CuePeripheral(:,Interval_Sel_Peripheral),2);


PSTH_Mean_GoodCentral_Interval=nanmean(PSTH_GoodCentral(:,Interval_Sel_Central),2);
PSTH_Mean_BadCentral_Interval=nanmean(PSTH_BadCentral(:,Interval_Sel_Central),2);


TargetAngle_GoodPeripheral=TargetAngle_Line(GoodPeripheralCondition);
TargetAngle_BadPeripheral=TargetAngle_Line(BadPeripheralCondition);
TargetAngle_CuePeripheral=TargetAngle_Line(CueOnFlag_Line==1);

TargetAngle_GoodCentral=TargetAngle_Line(GoodCentralCondition);
TargetAngle_BadCentral=TargetAngle_Line(BadCentralCondition);

%Transfer into location tuning
SeparationPoint=[22.5:45:360];
 
%Good Peripheral
 LocationTuning=PoloarHistGroup(TargetAngle_GoodPeripheral,PSTH_Mean_GoodPeripheral_Interval,SeparationPoint);
 LocationVector_GoodPeripheral=LocationTuning.Center;
 FR_Mean_GoodPeripheral=LocationTuning.DataMean;
 FR_Sem_GoodPeripheral=LocationTuning.DataSEM;
 GoodPeripheral_Tuning_p=LocationTuning.ANOVA_p;
 
 PrefLoc_GoodPeripheral=LocationTuning.PrefVector;
PrefLoc_Amp_GoodPeripheral=LocationTuning.PrefVector_Amp;

%Bad Peripheral
LocationTuning=PoloarHistGroup(TargetAngle_BadPeripheral,PSTH_Mean_BadPeripheral_Interval,SeparationPoint);
 LocationVector_BadPeripheral=LocationTuning.Center;
 FR_Mean_BadPeripheral=LocationTuning.DataMean;
 FR_Sem_BadPeripheral=LocationTuning.DataSEM;
 BadPeripheral_Tuning_p=LocationTuning.ANOVA_p;
 
 PrefLoc_BadPeripheral=LocationTuning.PrefVector;
PrefLoc_Amp_BadPeripheral=LocationTuning.PrefVector_Amp;

%CuePeripheral
LocationTuning=PoloarHistGroup(TargetAngle_CuePeripheral,PSTH_Mean_CuePeripheral_Interval,SeparationPoint);
 LocationVector_CuePeripheral=LocationTuning.Center;
 FR_Mean_CuePeripheral=LocationTuning.DataMean;
 FR_Sem_CuePeripheral=LocationTuning.DataSEM;
 CuePeripheral_Tuning_p=LocationTuning.ANOVA_p;
 
 PrefLoc_CuePeripheral=LocationTuning.PrefVector;
PrefLoc_Amp_CuePeripheral=LocationTuning.PrefVector_Amp;

%Good Central
 LocationTuning=PoloarHistGroup(TargetAngle_GoodCentral,PSTH_Mean_GoodCentral_Interval,SeparationPoint);
 LocationVector_GoodCentral=LocationTuning.Center;
 FR_Mean_GoodCentral=LocationTuning.DataMean;
 FR_Sem_GoodCentral=LocationTuning.DataSEM;
 GoodCentral_Tuning_p=LocationTuning.ANOVA_p;
 
 PrefLoc_GoodCentral=LocationTuning.PrefVector;
PrefLoc_Amp_GoodCentral=LocationTuning.PrefVector_Amp;

%Bad Central
LocationTuning=PoloarHistGroup(TargetAngle_BadCentral,PSTH_Mean_BadCentral_Interval,SeparationPoint);
 LocationVector_BadCentral=LocationTuning.Center;
 FR_Mean_BadCentral=LocationTuning.DataMean;
 FR_Sem_BadCentral=LocationTuning.DataSEM;
 BadCentral_Tuning_p=LocationTuning.ANOVA_p;
 
 PrefLoc_BadCentral=LocationTuning.PrefVector;
PrefLoc_Amp_BadCentral=LocationTuning.PrefVector_Amp;

%Put them together for plotting

FR_Mean_Peripheral=[FR_Mean_GoodPeripheral,FR_Mean_BadPeripheral,FR_Mean_CuePeripheral];
FR_Sem_Peripheral=[FR_Sem_GoodPeripheral,FR_Sem_BadPeripheral,FR_Sem_CuePeripheral];
Peripheral_Tuning_p=[GoodPeripheral_Tuning_p,BadPeripheral_Tuning_p,CuePeripheral_Tuning_p];
LocationVector_Peripheral=[LocationVector_GoodPeripheral',LocationVector_BadPeripheral',LocationVector_CuePeripheral'];
PrefLoc_Peripheral=[PrefLoc_GoodPeripheral,PrefLoc_BadPeripheral,PrefLoc_CuePeripheral];
PrefLoc_Amp_Peripheral=[PrefLoc_Amp_GoodPeripheral,PrefLoc_Amp_BadPeripheral,PrefLoc_Amp_CuePeripheral];

FR_Mean_Central=[FR_Mean_GoodCentral,FR_Mean_BadCentral];
FR_Sem_Central=[FR_Sem_GoodCentral,FR_Sem_BadCentral];
Central_Tuning_p=[GoodCentral_Tuning_p,BadCentral_Tuning_p];
LocationVector_Central=[LocationVector_GoodCentral',LocationVector_BadCentral'];
PrefLoc_Central=[PrefLoc_GoodCentral,PrefLoc_BadCentral];
PrefLoc_Amp_Central=[PrefLoc_Amp_GoodCentral,PrefLoc_Amp_BadCentral];

%Loop over objects on each scenes

UniqueSceneName=unique(SceneName_Line);

NumberOfScenes=length(UniqueSceneName);

%%
%%%%%%%%%%%%%Loop over the scenes
for i=1:NumberOfScenes
    currScene=UniqueSceneName(i);
    currGoodPeripheral=(SceneName_Line==currScene)&GoodPeripheralCondition;
    currBadPeripheral=(SceneName_Line==currScene)&BadPeripheralCondition;
    sceneSel=SceneName_Line==currScene;
    currCuePeripheral=sceneSel(CueOnFlag_Line==1);
    
    currGoodCentral=(SceneName_Line==currScene)&GoodCentralCondition;
    currBadCentral=(SceneName_Line==currScene)&BadCentralCondition;
    
    
   %PSTH for all these conditions
PSTH_GoodPeripheral=MeanPSTH_FR{1}(currGoodPeripheral,:);
PSTH_BadPeripheral=MeanPSTH_FR{1}(currBadPeripheral,:);
PSTH_CuePeripheral=MeanPSTH_FR{2}(currCuePeripheral,:);

PSTH_Time_GoodPeripheral_currScene{i}=nanmean(MeanPSTH_Time{1}(currGoodPeripheral,:),1);
PSTH_Time_BadPeripheral_currScene{i}=nanmean(MeanPSTH_Time{1}(currBadPeripheral,:),1);
PSTH_Time_CuePeripheral_currScene{i}=nanmean(MeanPSTH_Time{2}(currCuePeripheral,:),1);


PSTH_GoodCentral=MeanPSTH_FR{1}(currGoodCentral,:);
PSTH_BadCentral=MeanPSTH_FR{1}(currBadCentral,:);

PSTH_Time_GoodCentral_currScene{i}=nanmean(MeanPSTH_Time{1}(currGoodCentral,:),1);
PSTH_Time_BadCentral_currScene{i}=nanmean(MeanPSTH_Time{1}(currBadCentral,:),1);



%PSTH Average for good and bad and cue for central and periperal
PSTH_Mean_GoodPeripheral_currScene{i}=nanmean(PSTH_GoodPeripheral,1);
PSTH_Mean_BadPeripheral_currScene{i}=nanmean(PSTH_BadPeripheral,1);
PSTH_Mean_CuePeripheral_currScene{i}=nanmean(PSTH_CuePeripheral,1);

PSTH_Sem_GoodPeripheral_currScene{i}=nanstd(PSTH_GoodPeripheral,[],1)./sqrt(sum(~isnan(PSTH_GoodPeripheral),1));
PSTH_Sem_BadPeripheral_currScene{i}=nanstd(PSTH_BadPeripheral,[],1)./sqrt(sum(~isnan(PSTH_BadPeripheral),1));
PSTH_Sem_CuePeripheral_currScene{i}=nanstd(PSTH_CuePeripheral,[],1)./sqrt(sum(~isnan(PSTH_CuePeripheral),1));

PSTH_Mean_GoodCentral_currScene{i}=nanmean(PSTH_GoodCentral,1);
PSTH_Mean_BadCentral_currScene{i}=nanmean(PSTH_BadCentral,1);


PSTH_Sem_GoodCentral_currScene{i}=nanstd(PSTH_GoodCentral,[],1)./sqrt(sum(~isnan(PSTH_GoodCentral),1));
PSTH_Sem_BadCentral_currScene{i}=nanstd(PSTH_BadCentral,[],1)./sqrt(sum(~isnan(PSTH_BadCentral),1));

    
 
    SceneStr{i}=sprintf('Scene%d',currScene);
    
    
end

%Compare the location selectivity between the central condition and the
%peripheral condition

UniqueTargetAngle=uniquetol(TargetAngle_Line,0.001);
NumUniqueAngle=length(UniqueTargetAngle);

for i=1:NumUniqueAngle
    
    currAngle=UniqueTargetAngle(i);
    currGoodPeripheral=(TargetAngle_Line==currAngle)&GoodPeripheralCondition;
    currBadPeripheral=(TargetAngle_Line==currAngle)&BadPeripheralCondition;
    
    AngleSel=TargetAngle_Line==currAngle;
    
    currCuePeripheral=AngleSel(CueOnFlag_Line==1);
    
    currGoodCentral=(TargetAngle_Line==currAngle)&GoodCentralCondition;
    currBadCentral=(TargetAngle_Line==currAngle)&BadCentralCondition;
    
     %PSTH for all these conditions
PSTH_GoodPeripheral=MeanPSTH_FR{1}(currGoodPeripheral,:);
PSTH_BadPeripheral=MeanPSTH_FR{1}(currBadPeripheral,:);
PSTH_CuePeripheral=MeanPSTH_FR{2}(currCuePeripheral,:);

PSTH_Time_GoodPeripheral_currAngle{i}=nanmean(MeanPSTH_Time{1}(currGoodPeripheral,:),1);
PSTH_Time_BadPeripheral_currAngle{i}=nanmean(MeanPSTH_Time{1}(currBadPeripheral,:),1);
PSTH_Time_CuePeripheral_currAngle{i}=nanmean(MeanPSTH_Time{2}(currCuePeripheral,:),1);


PSTH_GoodCentral=MeanPSTH_FR{1}(currGoodCentral,:);
PSTH_BadCentral=MeanPSTH_FR{1}(currBadCentral,:);

PSTH_Time_GoodCentral_currAngle{i}=nanmean(MeanPSTH_Time{1}(currGoodCentral,:),1);
PSTH_Time_BadCentral_currAngle{i}=nanmean(MeanPSTH_Time{1}(currBadCentral,:),1);



%PSTH Average for good and bad and cue for central and periperal
PSTH_Mean_GoodPeripheral_currAngle{i}=nanmean(PSTH_GoodPeripheral,1);
PSTH_Mean_BadPeripheral_currAngle{i}=nanmean(PSTH_BadPeripheral,1);
PSTH_Mean_CuePeripheral_currAngle{i}=nanmean(PSTH_CuePeripheral,1);

PSTH_Sem_GoodPeripheral_currAngle{i}=nanstd(PSTH_GoodPeripheral,[],1)./sqrt(sum(~isnan(PSTH_GoodPeripheral),1));
PSTH_Sem_BadPeripheral_currAngle{i}=nanstd(PSTH_BadPeripheral,[],1)./sqrt(sum(~isnan(PSTH_BadPeripheral),1));
PSTH_Sem_CuePeripheral_currAngle{i}=nanstd(PSTH_CuePeripheral,[],1)./sqrt(sum(~isnan(PSTH_CuePeripheral),1));

PSTH_Mean_GoodCentral_currAngle{i}=nanmean(PSTH_GoodCentral,1);
PSTH_Mean_BadCentral_currAngle{i}=nanmean(PSTH_BadCentral,1);


PSTH_Sem_GoodCentral_currAngle{i}=nanstd(PSTH_GoodCentral,[],1)./sqrt(sum(~isnan(PSTH_GoodCentral),1));
PSTH_Sem_BadCentral_currAngle{i}=nanstd(PSTH_BadCentral,[],1)./sqrt(sum(~isnan(PSTH_BadCentral),1));


  
    
    
    
    
end



%}

%%

if ShowFigureFlag

%%
%%%%%%%%%%%%%%%%%%%%%%%% Figure 1 PSTH for good, bad, cue, peripheral and
%%%%%%%%%%%%%%%%%%%%%%%% central
%%%%%%%%%%%%%%%%%%%%%%%% bad%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figtitlestr{1}='PSTH_CenteralPeripheralGoodBadCue';

fig{1}=PrepareFigure(50,'w',[50,100, 1000,800],'Name',figtitlestr{1});
%Good,Bad,  Cue appear in peripheral
%subplot(2,2,3)
subplot(2,2,1)
title('Good Bad Cue in Peripheral','FontSize',15);
%axis off

%axes('Position',[0.1,0.7,0.4,0.15]);
%Good peripheral

shadedErrorBar(PSTH_Time_GoodPeripheral,PSTH_Mean_GoodPeripheral,PSTH_Sem_GoodPeripheral,'lineprops',{ColorForPeripheralMode{1}},'transparent',1,'patchSaturation',0.3);
hold on
%Bad peripheral
shadedErrorBar(PSTH_Time_BadPeripheral,PSTH_Mean_BadPeripheral,PSTH_Sem_BadPeripheral,'lineprops',{ColorForPeripheralMode{2}},'transparent',1,'patchSaturation',0.3);

shadedErrorBar(PSTH_Time_CuePeripheral,PSTH_Mean_CuePeripheral,PSTH_Sem_CuePeripheral,'lineprops',{ColorForPeripheralMode{3}},'transparent',1,'patchSaturation',0.3);

xlim([-PlotStartOffset(1), 600]);

set(findall(gca, 'Type', 'Line'),'LineWidth',3);
ylabel(sprintf('Firing Rate(Hz)'));
%ylabel('Firing Rate(Hz)')
%xticks([]);
%xticklabels([]);
legend({'Good','Bad','Cue'});
xlabel('Time from object onset;')
set(gca,'LineWidth',3,'FontSize',15);



subplot(2,2,2)
title('Good Bad in central','FontSize',15);
%axis off

%axes('Position',[0.1,0.7,0.4,0.15]);
%Good peripheral
shadedErrorBar(PSTH_Time_GoodCentral,PSTH_Mean_GoodCentral,PSTH_Sem_GoodCentral,'lineprops',{'#EA1A41'},'transparent',1,'patchSaturation',0.3);
hold on
%Bad peripheral
shadedErrorBar(PSTH_Time_BadCentral,PSTH_Mean_BadCentral,PSTH_Sem_BadCentral,'lineprops',{'#648FFF'},'transparent',1,'patchSaturation',0.3);



xlim([-PlotStartOffset(2), 600]);

set(findall(gca, 'Type', 'Line'),'LineWidth',3);
ylabel(sprintf('Firing Rate(Hz)'));
%ylabel('Firing Rate(Hz)')
%xticks([]);
%xticklabels([]);
legend({'Good','Bad'});
xlabel('Time from object onset;')
set(gca,'LineWidth',3,'FontSize',15);


subplot(2,2,3)
title('Good in Peripheral and Good in Central','FontSize',25);
%axis off

%axes('Position',[0.1,0.7,0.4,0.15]);
%Good peripheral
shadedErrorBar(PSTH_Time_GoodPeripheral,PSTH_Mean_GoodPeripheral,PSTH_Sem_GoodPeripheral,'lineprops',{'#EA1A41'},'transparent',1,'patchSaturation',0.3);
hold on
%Bad peripheral
shadedErrorBar(PSTH_Time_GoodCentral,PSTH_Mean_GoodCentral,PSTH_Sem_GoodCentral,'lineprops',{'#FFB000'},'transparent',1,'patchSaturation',0.3);


xlim([-PlotStartOffset(1), 600]);

% plot([min(TimeMarker(TimeMarker(:,2)>1,2)),min(TimeMarker(TimeMarker(:,2)>1,2))],yrange,'--k');
%xlim([-PlotStartOffset, median(TimeMarker(:,2))+PlotStopOffset]);
%xlim([-PlotStartOffset,max(min(TimeMarker(:,2))+PlotStopOffset,400)]);
%xlim([-PlotStartOffset,min(TimeMarker(:,2))+PlotStopOffset]);
%xlim([-PlotStartOffset,mean(TimeMarker(:,2))+PlotStopOffset]);
%xlim([-PlotStartOffset,max(mean(TimeMarker(:,2))+PlotStopOffset,800)]);
%{
%Add significance level
yrange=ylim;
h_AllvsSponEye(h_AllvsSponEye==0)=NaN;
sig=area(PSTH_Time_All(1):PSTH_Time_All(1)+length(h_AllvsSponEye)-1,h_AllvsSponEye*max(yrange(1),1)*0.6,0,'FaceColor','#FFC107','EdgeColor','#FFC107','FaceAlpha',0.5,'EdgeAlpha',0.5,'ShowBaseLine','off');
%}


set(findall(gca, 'Type', 'Line'),'LineWidth',3);
ylabel(sprintf('Firing Rate(Hz)'));
%ylabel('Firing Rate(Hz)')
%xticks([]);
%xticklabels([]);
legend({'Good P','Good C'});
xlabel('Time from object onset;')
set(gca,'LineWidth',3,'FontSize',15);


subplot(2,2,4)
title('Bad in Peripheral and Bad in Central','FontSize',15);
%axis off

%axes('Position',[0.1,0.7,0.4,0.15]);
%Good peripheral
shadedErrorBar(PSTH_Time_BadPeripheral,PSTH_Mean_BadPeripheral,PSTH_Sem_BadPeripheral,'lineprops',{'#648FFF'},'transparent',1,'patchSaturation',0.3);
hold on
%Bad peripheral
shadedErrorBar(PSTH_Time_BadCentral,PSTH_Mean_BadCentral,PSTH_Sem_BadCentral,'lineprops',{'#785EF0'},'transparent',1,'patchSaturation',0.3);



xlim([-PlotStartOffset(2), 600]);

% plot([min(TimeMarker(TimeMarker(:,2)>1,2)),min(TimeMarker(TimeMarker(:,2)>1,2))],yrange,'--k');
%xlim([-PlotStartOffset, median(TimeMarker(:,2))+PlotStopOffset]);
%xlim([-PlotStartOffset,max(min(TimeMarker(:,2))+PlotStopOffset,400)]);
%xlim([-PlotStartOffset,min(TimeMarker(:,2))+PlotStopOffset]);
%xlim([-PlotStartOffset,mean(TimeMarker(:,2))+PlotStopOffset]);
%xlim([-PlotStartOffset,max(mean(TimeMarker(:,2))+PlotStopOffset,800)]);
%{
%Add significance level
yrange=ylim;
h_AllvsSponEye(h_AllvsSponEye==0)=NaN;
sig=area(PSTH_Time_All(1):PSTH_Time_All(1)+length(h_AllvsSponEye)-1,h_AllvsSponEye*max(yrange(1),1)*0.6,0,'FaceColor','#FFC107','EdgeColor','#FFC107','FaceAlpha',0.5,'EdgeAlpha',0.5,'ShowBaseLine','off');
%}
legend({'Bad P','Bad C'});

set(findall(gca, 'Type', 'Line'),'LineWidth',3);
ylabel(sprintf('Firing Rate(Hz)'));
%ylabel('Firing Rate(Hz)')
%xticks([]);
%xticklabels([]);
xlabel('Time from object onset;')
set(gca,'LineWidth',3,'FontSize',15);





%%
%%%%%%%%%%%%%%%%%%%%%%%% Figure 2 PSTH for Each periperal fractal,Cue in each scene
%%%%%%%%%%%%%%%%%%%%%%%% bad%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figtitlestr{2}='PSTH_CueFractal_Peripheral';
%title('PSTH for each fractals','FontSize',25,'FontWeight','Bold');

figtitlestr{3}='PSTH_Fractal_Central';

fig{2}=PrepareFigure(51,'w',[70,100, 1000,800],'Name',figtitlestr{2});
fig{3}=PrepareFigure(52,'w',[90,100, 1000,800],'Name',figtitlestr{3});






for i=1:NumberOfScenes
   
    
set(groot,'CurrentFigure',fig{2});
subplot(ceil(NumberOfScenes/3),3,i);

title(SceneStr{i});

%Good peripheral
shadedErrorBar(PSTH_Time_GoodPeripheral_currScene{i},PSTH_Mean_GoodPeripheral_currScene{i},PSTH_Sem_GoodPeripheral_currScene{i},'lineprops',{'#EA1A41'},'transparent',1,'patchSaturation',0.3);
hold on
%Bad peripheral
shadedErrorBar(PSTH_Time_BadPeripheral_currScene{i},PSTH_Mean_BadPeripheral_currScene{i},PSTH_Sem_BadPeripheral_currScene{i},'lineprops',{'#648FFF'},'transparent',1,'patchSaturation',0.3);

shadedErrorBar(PSTH_Time_CuePeripheral_currScene{i},PSTH_Mean_CuePeripheral_currScene{i},PSTH_Sem_CuePeripheral_currScene{i},'lineprops',{'#36AF64'},'transparent',1,'patchSaturation',0.3);

xlim([-PlotStartOffset(1), 600]);

set(findall(gca, 'Type', 'Line'),'LineWidth',3);

%ylabel('Firing Rate(Hz)')
%xticks([]);
%xticklabels([]);
if i==1
legend({'Good','Bad','Cue'});
end
if mod(i,4)==1
    ylabel(sprintf('Firing Rate(Hz)'));
end
xlabel('Time from object onset;')
set(gca,'LineWidth',3,'FontSize',15);








set(groot,'CurrentFigure',fig{3});
%axes('Position',[0.1,0.7,0.4,0.15]);
subplot(ceil(NumberOfScenes/3),3,i);

%Good peripheral
shadedErrorBar(PSTH_Time_GoodCentral_currScene{i},PSTH_Mean_GoodCentral_currScene{i},PSTH_Sem_GoodCentral_currScene{i},'lineprops',{'#EA1A41'},'transparent',1,'patchSaturation',0.3);
hold on
%Bad peripheral
shadedErrorBar(PSTH_Time_BadCentral_currScene{i},PSTH_Mean_BadCentral_currScene{i},PSTH_Sem_BadCentral_currScene{i},'lineprops',{'#648FFF'},'transparent',1,'patchSaturation',0.3);



xlim([-PlotStartOffset(2), 600]);

set(findall(gca, 'Type', 'Line'),'LineWidth',3);
ylabel(sprintf('Firing Rate(Hz)'));
%ylabel('Firing Rate(Hz)')
%xticks([]);
%xticklabels([]);
legend({'Good','Bad'});
xlabel('Time from object onset;')
set(gca,'LineWidth',3,'FontSize',15);



end

%%
%%%%%%%%%%%%%%%%%%%%%%%% Figure 4 PSTH for Each periperal fractal,Cue for
%%%%%%%%%%%%%%%%%%%%%%%% each angle
%%%%%%%%%%%%%%%%%%%%%%%% bad%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figtitlestr{4}='PSTH_CueFractal_PeripheralCentral_EachAngle';
%title('PSTH for each fractals','FontSize',25,'FontWeight','Bold');



fig{4}=PrepareFigure(53,'w',[110,100, 1000,800],'Name',figtitlestr{4});


for i=1:NumUniqueAngle
   
    
subplot(2,3,1);
if i==1
title('GoodPeripheral');
end

%Good peripheral
shadedErrorBar(PSTH_Time_GoodPeripheral_currAngle{i},PSTH_Mean_GoodPeripheral_currAngle{i},PSTH_Sem_GoodPeripheral_currAngle{i},'lineprops',{Color(i,:)},'transparent',1,'patchSaturation',0.3);
hold on
xlim([-PlotStartOffset(1), 600]);

set(findall(gca, 'Type', 'Line'),'LineWidth',3);

ylabel('Firing Rate(Hz)');
if i==1
xlabel('Time from object onset');
ylabel('Firing Rate(Hz)');

 

 
end

if i==NumUniqueAngle
    %Add a shaded area to indicate analysis area
 
 yrange=ylim;
area(AverageInterval_Peripheral,[yrange(2) yrange(2)],...
    'FaceColor','k','FaceAlpha',.1,'EdgeAlpha',.1,'HandleVisibility','off');

end


set(gca,'LineWidth',3,'FontSize',15);

subplot(2,3,2);
if i==1
title('BadPeripheral');
end
%Bad peripheral
shadedErrorBar(PSTH_Time_BadPeripheral_currAngle{i},PSTH_Mean_BadPeripheral_currAngle{i},PSTH_Sem_BadPeripheral_currAngle{i},'lineprops',{Color(i,:)},'transparent',1,'patchSaturation',0.3);
hold on
xlim([-PlotStartOffset(1), 600]);
if i==1
xlabel('Time from object onset');
%Add a shaded area to indicate analysis area
 
 yrange=ylim;
area(AverageInterval_Peripheral,[yrange(2) yrange(2)],...
    'FaceColor','k','FaceAlpha',.1,'EdgeAlpha',.1,'HandleVisibility','off');

end
set(findall(gca, 'Type', 'Line'),'LineWidth',3);
set(gca,'LineWidth',3,'FontSize',15);

subplot(2,3,3);
%Cue Peripheral
if i==1
title('CuePeripheral');

end

shadedErrorBar(PSTH_Time_CuePeripheral_currAngle{i},PSTH_Mean_CuePeripheral_currAngle{i},PSTH_Sem_CuePeripheral_currAngle{i},'lineprops',{Color(i,:)},'transparent',1,'patchSaturation',0.3);
hold on
xlim([-PlotStartOffset(1), 600]);

set(findall(gca, 'Type', 'Line'),'LineWidth',3);

%ylabel('Firing Rate(Hz)')
%xticks([]);
%xticklabels([]);


if i==1
xlabel('Time from object onset');


end
set(gca,'LineWidth',3,'FontSize',15);
 %Add a shaded area to indicate analysis area
 if i==NumUniqueAngle
    %Add a shaded area to indicate analysis area
 
 yrange=ylim;
area(AverageInterval_Peripheral,[yrange(2) yrange(2)],...
    'FaceColor','k','FaceAlpha',.1,'EdgeAlpha',.1,'HandleVisibility','off');

end
 hold on










subplot(2,3,4);
if i==1
title('GoodCentral');

end


%Good Central
shadedErrorBar(PSTH_Time_GoodCentral_currAngle{i},PSTH_Mean_GoodCentral_currAngle{i},PSTH_Sem_GoodCentral_currAngle{i},'lineprops',{Color(i,:)},'transparent',1,'patchSaturation',0.3);
hold on
hold on
xlim([-PlotStartOffset(1), 600]);

set(findall(gca, 'Type', 'Line'),'LineWidth',3);
if i==1
    xlabel('Time from object onset;');
     ylabel('Firing Rate(Hz)');
end
set(gca,'LineWidth',3,'FontSize',15);
if i==NumUniqueAngle
    %Add a shaded area to indicate analysis area
 
 yrange=ylim;
area(AverageInterval_Peripheral,[yrange(2) yrange(2)],...
    'FaceColor','k','FaceAlpha',.1,'EdgeAlpha',.1,'HandleVisibility','off');

end


subplot(2,3,5);
title('BadCentral');
%Bad Central
shadedErrorBar(PSTH_Time_BadCentral_currAngle{i},PSTH_Mean_BadCentral_currAngle{i},PSTH_Sem_BadCentral_currAngle{i},'lineprops',{Color(i,:)},'transparent',1,'patchSaturation',0.3);
hold on

xlim([-PlotStartOffset(2), 600]);
if i==1
    xlabel('Time from object onset;');
   

end
set(findall(gca, 'Type', 'Line'),'LineWidth',3);
%ylabel('Firing Rate(Hz)')
%xticks([]);
%xticklabels([]);
%legend({'Good','Bad'});
if i==NumUniqueAngle
    %Add a shaded area to indicate analysis area
 
 yrange=ylim;
area(AverageInterval_Peripheral,[yrange(2) yrange(2)],...
    'FaceColor','k','FaceAlpha',.1,'EdgeAlpha',.1,'HandleVisibility','off');

end

set(gca,'LineWidth',3,'FontSize',15);



end
legend({num2str(round(UniqueTargetAngle))});





%%
%%%%%%%%%%%%%%%%%%%%%%%% Figure 5 Polar plot for each location tuning
figtitlestr{5}='LocationTuning_DifferentConditions';

fig{5}=PrepareFigure(54,'w',[110,200, 500,500],'Name',figtitlestr{5});

ConditionHeader_P={'GP','BP','CP'};
%Pre-saccade , During, Post saccade tuning
%subplot(1,3,j);
%{
p= polarplot((SaccadeAngle_First)*pi/180',SpikeRate_PreSac,'or');
p.MarkerFaceColor='r';
p.MarkerEdgeColor='r';
%}
%Peripheral condition
subplot(1,2,1);
for j=1:3
%{    
FR_Mean_Peripheral=[FR_Mean_GoodPeripheral,FR_Mean_BadPeripheral,FR_Mean_CuePeripheral];
FR_Sem_Peripheral=[FR_Sem_GoodPeripheral,FR_Sem_BadPeripheral,FR_Sem_CuePeripheral];
Peripheral_Tuning_p=[GoodPeripheral_Tuning_p,BadPeripheral_Tuning_p,CuePeripheral_Tuning_p];
LocationVector_Peripheral=[LocationVector_GoodPeripheral',LocationVector_BadPeripheral',LocationVector_CuePeripheral'];
PrefLoc_Peripheral=[PrefLoc_GoodPeripheral,PrefLoc_BadPeripheral,PrefLoc_CuePeripheral];
PrefLoc_Amp_Peripheral=[PrefLoc_Amp_GoodPeripheral,PrefLoc_Amp_BadPeripheral,PrefLoc_Amp_CuePeripheral]; 
   %}
%{  
Selection_temp=~isnan(SaccadeVector);
SaccadeVector=SaccadeVector(Selection_temp);
FiringRate_PreDurPost_mean=FiringRate_PreDurPost_mean(Selection_temp,:);
FiringRate_PreDurPost_sem=FiringRate_PreDurPost_sem(Selection_temp,:);
%}    
%hold on

style_str=strcat('-',ColorForPeripheralMode{j});

polarwitherrorbar([LocationVector_Peripheral(:,j)'/180*pi (LocationVector_Peripheral(1,j))/180*pi],[FR_Mean_Peripheral(:,j)' FR_Mean_Peripheral(1,j)],[FR_Sem_Peripheral(:,j)' FR_Sem_Peripheral(1,j)],...
    style_str,3);
%set(gca,'Color',ColorForPeripheralMode{j});
hold on
h_p = polarplot([PrefLoc_Peripheral(j) PrefLoc_Peripheral(j)]/180*pi,[0 PrefLoc_Amp_Peripheral(j)],'color',ColorForPeripheralMode{j});
    set(h_p,'LineWidth',2.5);

set(gca,'ThetaDir' , 'counterclockwise');
set(gca,'ThetaZeroLocation','top')
set(gca,'FontSize',25,'FontWeight','Bold','LineWidth',3);

thetaticks(sort(LocationVector_Peripheral(:,j)));
thetaticklabels(sort(LocationVector_Peripheral(:,j)));
if j==1
title('Peripheral','FontSize',20,'FontWeight','Bold');
end
pstr=sprintf('%s: p= %2.3f',ConditionHeader_P{j},Peripheral_Tuning_p(j));
a=annotation('textbox',[0.15,max(0.06-(j-1)/35,0.001),0.1,0.1] ,'String',pstr,'FitBoxToText','on','LineStyle','none');
a.FontSize=15;
a.FontWeight='Bold';
if Peripheral_Tuning_p(j)<0.05
    a.Color=ColorForPeripheralMode{j};
end

end
%Center condition
subplot(1,2,2);

ConditionHeader_C={'GC','BC'};
for j=1:2
%{    
FR_Mean_Peripheral=[FR_Mean_GoodPeripheral,FR_Mean_BadPeripheral,FR_Mean_CuePeripheral];
FR_Sem_Peripheral=[FR_Sem_GoodPeripheral,FR_Sem_BadPeripheral,FR_Sem_CuePeripheral];
Peripheral_Tuning_p=[GoodPeripheral_Tuning_p,BadPeripheral_Tuning_p,CuePeripheral_Tuning_p];
LocationVector_Peripheral=[LocationVector_GoodPeripheral',LocationVector_BadPeripheral',LocationVector_CuePeripheral'];
PrefLoc_Peripheral=[PrefLoc_GoodPeripheral,PrefLoc_BadPeripheral,PrefLoc_CuePeripheral];
PrefLoc_Amp_Peripheral=[PrefLoc_Amp_GoodPeripheral,PrefLoc_Amp_BadPeripheral,PrefLoc_Amp_CuePeripheral]; 
   %}
%{  
Selection_temp=~isnan(SaccadeVector);
SaccadeVector=SaccadeVector(Selection_temp);
FiringRate_PreDurPost_mean=FiringRate_PreDurPost_mean(Selection_temp,:);
FiringRate_PreDurPost_sem=FiringRate_PreDurPost_sem(Selection_temp,:);
%}    
%hold on
style_str=strcat('-',ColorForCentralMode{j});

polarwitherrorbar([LocationVector_Central(:,j)'/180*pi (LocationVector_Central(1,j))/180*pi],[FR_Mean_Central(:,j)' FR_Mean_Central(1,j)],[FR_Sem_Central(:,j)' FR_Sem_Central(1,j)],...
    style_str,3);
%set(gca,'Color',ColorForPeripheralMode{j});
hold on
h_p = polarplot([PrefLoc_Central(j) PrefLoc_Central(j)]/180*pi,[0 PrefLoc_Amp_Central(j)],'color',ColorForCentralMode{j});
    set(h_p,'LineWidth',2.5);

set(gca,'ThetaDir' , 'counterclockwise');
set(gca,'ThetaZeroLocation','top')
set(gca,'FontSize',25,'FontWeight','Bold','LineWidth',3);

thetaticks(sort(LocationVector_Central(:,j)));
thetaticklabels(sort(LocationVector_Central(:,j)));
if j==1
title('Central','FontSize',20,'FontWeight','Bold');
end
pstr=sprintf('%s: p= %2.3f',ConditionHeader_C{j},Central_Tuning_p(j));
a=annotation('textbox',[0.65,max(0.06-(j-1)/35,0.001),0.1,0.1] ,'String',pstr,'FitBoxToText','on','LineStyle','none');
a.FontSize=15;
a.FontWeight='Bold';
if Central_Tuning_p(j)<0.05
    a.Color=ColorForCentralMode{j};
end

end


%{
%%%%%%%%%%%%%%%%%%%%%%%% Figure 5 PSTH for each object
%%%%%%%%%%%%%%%%%%%%%%%% bad%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for j=1:ceil(NumberOfObject/16)
figtitlestr{5+j-1}=sprintf('PSTH_EachFractal_%d',j);

fig{5+j-1}=PrepareFigure(108+j,'w',[160+5*(j-1),100, 1000,800],'Name',figtitlestr{5+j-1});

NumberOfObjectInEachFigure=16;
Indicator=1;
for ii=(j-1)*NumberOfObjectInEachFigure+1:min(j*NumberOfObjectInEachFigure,NumberOfObject)
    subplot(ceil(NumberOfObjectInEachFigure/4),4,Indicator);
   
shadedErrorBar(PSTH_Time_All,Mean_PSTH_FR_EachObject(ii,:),SEM_PSTH_FR_EachObject(ii,:),'lineprops',{Color(ii,:)},'transparent',1,'patchSaturation',0.3);
set(findall(gca, 'Type', 'Line'),'LineWidth',3);
%ylim([-0.5*min(min(Mean_PSTH_FR_EachObject)),max(max(Mean_PSTH_FR_EachObject+SEM_PSTH_FR_EachObject))*1.02]);
%ylim([-0.5*min(min(Mean_PSTH_FR_EachObject)),max(max(Mean_PSTH_FR_EachObject(ceil(PlotStartOffset):ceil(PlotStartOffset)+1000)+SEM_PSTH_FR_EachObject(ceil(PlotStartOffset):ceil(PlotStartOffset)+1000)))*1.02])

%xlim([min(min(PSTH_Time_All))-10,max(max(PSTH_Time_All))+10]);
%xlim([-PlotStartOffset,min(TimeMarker(:,2))+PlotStopOffset]);
%xlim([-PlotStartOffset, median(TimeMarker(:,2))+PlotStopOffset]);
%xlim([-PlotStartOffset,max(min(TimeMarker(:,2))+PlotStopOffset,400)]);
xlim([-PlotStartOffset,min(TimeMarker(:,2))+PlotStopOffset]);
%ylabel('Spikes/s');

set(gca,'LineWidth',3,'FontSize',20,'FontWeight','Bold');
%ylim([0,60]);

hold on
%Add significance level
yrange=ylim;
h_EachScene(ii,h_EachScene(ii,:)==0)=NaN;
h_EachScene(ii,PSTH_Time_All<0)=NaN;

sig=area(PSTH_Time_All,h_EachScene(ii,:)*max(yrange(1)*0.6,3),0,'FaceColor','#FFC107','EdgeColor','#FFC107','FaceAlpha',0.5,'EdgeAlpha',0.5,'ShowBaseLine','off');



set(findall(gca, 'Type', 'Line'),'LineWidth',3);
%ylabel(sprintf('Spikes/s'));
%ylabel('Firing Rate(Hz)')
%xticks([]);
%xticklabels([]);
set(gca,'LineWidth',3,'FontSize',25);
title(FractalStr{ii},'FontSize',20,'FontWeight','Bold');
  
Indicator=Indicator+1;  
    
end
end



%{ 
%Temporally not need to show eye trace 
%%%%%%%%%%%%%%%%%%%%%%%% Figure 6 Eye Trace for each object
%%%%%%%%%%%%%%%%%%%%%%%% bad%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for j=1:ceil(NumberOfObject/16)
figtitlestr{8+j-1}=sprintf('PSTH_EachFractal_%d',j);

fig{8+j-1}=PrepareFigure(111+j,'w',[175+10*(j-1),100, 1000,800],'Name',figtitlestr{8+j-1});

NumberOfObjectInEachFigure=16;
Indicator=1;
for ii=(j-1)*NumberOfObjectInEachFigure+1:min(j*NumberOfObjectInEachFigure,NumberOfObject)
    subplot(ceil(NumberOfObjectInEachFigure/4),4,Indicator);
 
shadedErrorBar([-PlotStartOffset:-PlotStartOffset+length( Mean_Eye_Eccentricity_ForObject(ii,:))-1], Mean_Eye_Eccentricity_ForObject(ii,:), SEM_Eye_Eccentricity_ForObject(ii,:),'lineprops',{Color(ii,:)},'transparent',1,'patchSaturation',0.3);
 
  
%shadedErrorBar(PSTH_Time_All,Mean_PSTH_FR_EachObject(ii,:),SEM_PSTH_FR_EachObject(ii,:),'lineprops',{Color(ii,:)},'transparent',1,'patchSaturation',0.3);
set(findall(gca, 'Type', 'Line'),'LineWidth',3);
%ylim([-0.5*min(min(Mean_PSTH_FR_EachObject)),max(max(Mean_PSTH_FR_EachObject+SEM_PSTH_FR_EachObject))*1.02]);
%ylim([-0.5*min(min(Mean_Eye_Eccentricity_ForObject)),max(max(Mean_Eye_Eccentricity_ForObject(ceil(PlotStartOffset):ceil(PlotStopOffset))+SEM_Eye_Eccentricity_ForObject(ceil(PlotStartOffset):ceil(PlotStopOffset))))*1.02])

%xlim([min(min(PSTH_Time_All))-10,max(max(PSTH_Time_All))+10]);
%xlim([-PlotStartOffset,min(TimeMarker(:,2))+PlotStopOffset]);
%xlim([-PlotStartOffset, median(TimeMarker(:,2))+PlotStopOffset]);

%xlim([-PlotStartOffset,max(min(TimeMarker(:,2))+PlotStopOffset,400)]);
xlim([-PlotStartOffset,min(TimeMarker(:,2))+PlotStopOffset]);
%ylabel('Spikes/s');

set(gca,'LineWidth',3,'FontSize',20,'FontWeight','Bold');
%ylim([0,60]);

hold on
%{
%Add significance level
yrange=ylim;
h_EachScene(ii,h_EachScene(ii,:)==0)=NaN;
h_EachScene(ii,PSTH_Time_All<0)=NaN;

sig=area(PSTH_Time_All,h_EachScene(ii,:)*max(yrange(1)*0.6,3),0,'FaceColor','#FFC107','EdgeColor','#FFC107','FaceAlpha',0.5,'EdgeAlpha',0.5,'ShowBaseLine','off');

%}

set(findall(gca, 'Type', 'Line'),'LineWidth',3);
%ylabel(sprintf('Spikes/s'));

%xticks([]);
%xticklabels([]);
set(gca,'LineWidth',3,'FontSize',25);
ylabel('Eccentricity','FontSize',15)
title(FractalStr{ii},'FontSize',20,'FontWeight','Bold');
  
Indicator=Indicator+1;  
    
end
end
%}
%}
end %End of if showfigureflag




%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Output to files%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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





%MarkerFolder='DataAnalysis';

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
    if isfield(OutputFileName,'CPTargetResponse')
    ExistFlag=1;
    end
end


%%%%%%%%%%%%%%%%%%%%OutputData For current analysis
OutputData_New=[];

OutputData_New.CPTargetResponse.TrialType='GoodTrials';

%%%%%%%%%Save overall scene response characteristics%%%%%%%%%%%%%%%%%%%%%%%
%Organize the data
NameStr={'TrialNum',...
    'PSTH_Good_Peripheral_Average',...
    'PSTH_Bad_Peripheral_Average',...
    'PSTH_Cue_Peripheral_Average',...
    'PSTH_Good_Central_Average',...
    'PSTH_Bad_Central_Average',...
    'PSTH_Good_Peripheral_Average_EachScene',...
    'PSTH_Bad_Peripheral_Average_EachScene',...
    'PSTH_Cue_Peripheral_Average_EachScene',...
    'PSTH_Good_Central_Average_EachScene',...
    'PSTH_Bad_Central_Average_EachScene',...
    
    

};



DataLib={TrialNum,...
    [PSTH_Time_GoodPeripheral',PSTH_Mean_GoodPeripheral',PSTH_Sem_GoodPeripheral'],...
    [PSTH_Time_BadPeripheral',PSTH_Mean_BadPeripheral',PSTH_Sem_BadPeripheral'],...
    [PSTH_Time_CuePeripheral',PSTH_Mean_CuePeripheral',PSTH_Sem_CuePeripheral'],...
    [PSTH_Time_GoodCentral',PSTH_Mean_GoodCentral',PSTH_Sem_GoodCentral'],...
    [PSTH_Time_BadCentral',PSTH_Mean_BadCentral',PSTH_Sem_BadCentral'],...    
    {PSTH_Time_GoodPeripheral_currScene,PSTH_Mean_GoodPeripheral_currScene,PSTH_Sem_GoodPeripheral_currScene},...
    {PSTH_Time_BadPeripheral_currScene,PSTH_Mean_BadPeripheral_currScene,PSTH_Sem_BadPeripheral_currScene},...
    {PSTH_Time_CuePeripheral_currScene,PSTH_Mean_CuePeripheral_currScene,PSTH_Sem_CuePeripheral_currScene},...
    {PSTH_Time_GoodCentral_currScene,PSTH_Mean_GoodCentral_currScene,PSTH_Sem_GoodCentral_currScene},...
     {PSTH_Time_BadCentral_currScene,PSTH_Mean_BadCentral_currScene,PSTH_Sem_BadCentral_currScene}
    
};
DataStamp=containers.Map(NameStr,DataLib);


OutputData_New.CPTargetResponse.Task=TaskType;
OutputData_New.CPTargetResponse.TaskCode=TaskCode;
OutputData_New.CPTargetResponse.DataStamp=DataStamp;


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
        
        OutputFigureName_Final=strcat(OutputFigureName,figtitlestr{jj},'.jpg');
        saveas(fig{jj},OutputFigureName_Final);
        
            
        
    end

    disp('Figures have been exported to the neuron folder');

end

if ExistFlag
%Compare the old one with the new one
Task_Old=OutputData.CPTargetResponse.TaskCode;

     if Task_Old~=TaskCode
        %Not the same task, add the new dataset into the old one
        OutputData.CPTargetResponse(numel(OutputData.CPTargetResponse)+1)=OutputData_New.CPTargetResponse;
     else
         %If the same task,the same protocol,replace the old one with the new one
         
         OutputData.CPTargetResponse=OutputData_New.CPTargetResponse;
         
    
    
     end
else
    %Output the current file
   OutputData.CPTargetResponse=OutputData_New.CPTargetResponse; 
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