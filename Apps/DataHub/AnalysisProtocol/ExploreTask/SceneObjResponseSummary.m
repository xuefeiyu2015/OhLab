function SceneObjResponseSummary(Data, StartTrial, EndTrial,ShowFigureFlag,OutputFlag);
%Summary analysis for each neuron for an overview of scene and object
%responses;
%By Xuefei Yu 05012021
ECode;
FileName=Data.FileName;
TaskType=Data.TaskType;
TaskCode=Data.TaskCode;
%TrialNum=EndTrial-StartTrial+1;
Selection=StartTrial:EndTrial;%Select the trials according to the pannel

%RemoveTrial=129;
%Selection=Selection(~ismember(Selection,RemoveTrial));

%Abnormal check %Especially For Robin081120FORAGEXPLO4.rst
CheckCode=SCENE_ON;%This code should appear once in each trial
AbnormalTrials=AbnormalTrialCheck(Data,CheckCode);
Selection=Selection(~ismember(Selection,AbnormalTrials));


CheckCode=SCENE_OFF;%This code should appear once in each trial
AbnormalTrials=AbnormalTrialCheck(Data,CheckCode);
Selection=Selection(~ismember(Selection,AbnormalTrials));


CheckCodePair=[TGTCD,TGTOFF];%This code should be the same as target on for good trials
AbnormalTrials=AbnormalTrialCheck(Data,CheckCodePair);
Selection=Selection(~ismember(Selection,AbnormalTrials));

%{
if strcmp(FileName,'Robin092920FORAGEXPLO1.rst') ||strcmp(FileName,'Robin092820FORAGEXPLO1.rst')||strcmp(FileName,'Robin012621FORAGEXPLO1.rst')||strcmp(FileName,'Robin040221FORAGEXPLOITL1.rst')||strcmp(FileName,'Robin111520FORAGEXPLO3.rst')||strcmp(FileName,'Robin020721FORAGEXPLO1.rst')
    %Remove bad trial for this file

CheckCodePair=[FIXOFF,SHOWFIXCD];%This code should be the same as target on for good trials
AbnormalTrials=AbnormalTrialCheck(Data,CheckCodePair);
Selection=Selection(~ismember(Selection,AbnormalTrials));

if strcmp(FileName,'Robin092820FORAGEXPLO1.rst') ||strcmp(FileName,'Robin012621FORAGEXPLO1.rst')||strcmp(FileName,'Robin040221FORAGEXPLOITL1.rst')||strcmp(FileName,'Robin111520FORAGEXPLO3.rst')||strcmp(FileName,'Robin020721FORAGEXPLO1.rst')
    CheckCodePair=[TGTCD,TGTOFF];%This code should be the same as target on for good trials
    AbnormalTrials=AbnormalTrialCheck(Data,CheckCodePair);
    Selection=Selection(~ismember(Selection,AbnormalTrials));
    
end


%TGTCD
end
%}
TrialIndexRaw=Data.BlockType(:,1);

SceneOnAveTime=600;%Average Interval For SceneOnResponses

%Special request
remove_tmp=[];
SceneRemoval=[];
BlockRemoval=[];

%Remove the intermediate scene to focus on the pure scenes if needed
%SceneRemoval=[1054:1057,1060:1063];
if ~isempty(SceneRemoval)
    SceneName_tmp=Data.SceneIndex(Selection,3);
    remove_tmp=ismember(SceneName_tmp,SceneRemoval);
      
end

%Remove the not focus block type to focus on the random blocks if needed
%BlockRemoval=[0,10];
if ~isempty(BlockRemoval)
    Block_tmp=Data.BlockType(Selection,2);
    if ~isempty(remove_tmp)
        remove_tmp=remove_tmp | ismember(Block_tmp,BlockRemoval);
    else
        remove_tmp=ismember(Block_tmp,BlockRemoval);
    end
    
   
end
if ~isempty(remove_tmp)
    Selection=Selection(~remove_tmp);
end
TrialIndex=TrialIndexRaw(Selection);

TrialNum=length(Selection);
%User defined parameter for the neural/behavioral analysis
TargetOnOffset=50;%in ms, the time selection bin after target onset; 
StartDelay=0;%in ms, the time selection bin after fixation offset;

%For segmentation of  the neural responses 
%SegEvent=[EXTRAFIXSTARTCD,EXTRAFIXOFF,SCENE_ON,SCENE_OFF,SHOWFIXCD,FIXOFF,TGTCD];%Event Flow of interest
%SegEventEnd=[EXTRAFIXOFF,SHOWFIXCD,SHOWFIXCD,NaN,FIXOFF,TGTCD,TGTOFF];%Event end; %Automatic screening later
%{
SegEvent=[EXTRAFIXSTARTCD,EXTRAFIXOFF,SCENE_ON,SCENE_OFF,FIXOFF,CUEON,TGTCD, TGTOFF,RWDCD,   SCENE_OFF,ENABLECD];%Event Flow of interest
SegEventEnd=[EXTRAFIXOFF, PREFIX	,SHOWFIXCD,NaN,      CUEON, TGTCD,TGTOFF,NaN,   RWDOFFCD,NaN,      EXTRAFIXSTARTCD];%Event end; %Automatic screening later
%SegEventEnd=[EXTRAFIXOFF,FIXOFF,FIXOFF,TGTCD,TGTOFF];%Event end
SegEventStr={'ExtraFPOn','ExtraFPOff','SceneOn','SceneOff','FixOff','CueOn','TargetOnset','TargetOffset','Reward','PostTrialSpon','BeforeTrialSpon'};

PlotStartOffset=[200,200,600,500,1000,500,100,100,100,-1500,-100];%in ms
PlotStopOffset=[150,200,500,1500,100,50,100,500,100,2000,-100];%in ms

EventTypeFlag=[1,1,1,1,2,2,2,2,1,1,1];%1:Number same with the trial number, 2:Number same with the target number
%}
SegEvent=[SCENE_ON,    SCENE_OFF,TGTCD,  TGTOFF,RWDCD,  ];%Event Flow of interest
SegEventEnd=[SHOWFIXCD,NaN,      TGTOFF, NaN,   NaN];%Event end; %Automatic screening later
%SegEventEnd=[EXTRAFIXOFF,FIXOFF,FIXOFF,TGTCD,TGTOFF];%Event end
SegEventStr={'SceneOn','SceneOff','TargetOnset','TargetOffset','Reward'};

PlotStartOffset=[1000,1800,800,500,100];%in ms
PlotStopOffset=[500,1500,100,1000,300];%in ms

EventTypeFlag=[1,1,2,2,1];%1:Number same with the trial number, 2:Number same with the target number

 
RasterOn=0;%Show raster, but it will be very slow to do so

%Parameter to define the sliding window to measure firing rates
StepSize=10;%in ms
BinWidth=20;%in ms

%Parameter for raster plot
Barlength=1;
LineWidth=1;


%ResponseLatencyThreshold=3;

%Some user defined parameter for visualization
%Color Definition for each scenes
%ColorForScene={'r','g','b','m',[0.4,0.4,0],[0.28,0.78,0.89],'k',[0.6,0.1,1],[0.6,0.1,0.4],[0.3,0.5,1],[0.6,0.5,0.1],[0,0.8,0.4]};
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
ColorForOrder={'#FB5441','#4173fb','#FBBB41','#0EC335',[0.28,0.78,0.89],'k',[0.6,0.1,1],[0.6,0.1,0.4],[0.3,0.5,1],[0.6,0.5,0.1],[0,0.8,0.4],'r','g','b','m'};

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



%Automatically chanage the marker

StartCheck=FindOutTime(EventChannel,EventBin,SegEvent(SegEvent==SCENE_ON));
EndCheck=FindOutTime(EventChannel,EventBin,SegEventEnd(SegEvent==SCENE_ON));
EndCheck=EndCheck(:,1);
DiffCheck=EndCheck-StartCheck;
if sum(DiffCheck<200)>ceil(length(StartCheck)/2)
    SegEventEnd(SegEvent==SCENE_ON)=[FIXOFF];%Event end;change from showfix cd to fix off for contigency
    disp('Change to FixOff for FPOn First Mode');
    
end




if RecordDate<datestr('10,28,2020') && strcmp(MonkeyName,'Robin')
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
end %End of the if data & monkey
%}


%Remove extraFP marker if it is not existed
CheckExtraFP=FindOutTime(EventChannel,EventBin,EXTRAFIXSTARTCD);
if sum(isnan(CheckExtraFP))==length(CheckExtraFP)
    Remove=ismember(SegEvent,[EXTRAFIXSTARTCD,EXTRAFIXOFF]);
   % SetNaNStart=ismember(SegEvent,EXTRAFIXOFF);
    
    
    PlotStartOffset(Remove)=[];
    PlotStopOffset(Remove)=[];
    EventTypeFlag(Remove)=[];
    SegEvent(Remove)=[];
    SegEventEnd(Remove)=[];
    SegEventStr(Remove)=[];
    
   % PlotStartOffset(SetNaNStart)=1000;
    %PlotStopOffset(Remove)=[];
   % SegEvent(SetNaNStart)=NaN;
   
    SetNaNEnd=ismember(SegEventEnd,EXTRAFIXSTARTCD);
    PlotStopOffset(SetNaNEnd)=1000;
    SegEventEnd(SetNaNEnd)=NaN;
    disp('Remove extraFP Marker');
    
   

    %SegEvent=[EXTRAFIXSTARTCD,EXTRAFIXOFF,SCENE_ON,SCENE_OFF,FIXOFF,CUEON,TGTCD, TGTOFF,RWDCD,   SCENE_OFF,ENABLECD];%Event Flow of interest
%SegEventEnd=[EXTRAFIXOFF, PREFIX	,SHOWFIXCD,NaN,      CUEON, TGTCD,TGTOFF,NaN,   RWDOFFCD,NaN,      EXTRAFIXSTARTCD];%Event end; %Automatic screening later
    
end

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
 %screen out the searching trials
 FixOffTime=FixOffTime(AdvanceInTargetTrial,:);
 OutScreenTime=OutScreenTime(AdvanceInTargetTrial,:);
 
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
 
 ChoiceOutcome=ChoiceOutcome(AdvanceInTargetTrial,:);
 
 TrialIndex=TrialIndex(AdvanceInTargetTrial);
 
% GOODCHOICECD    
 
 TargetNum=TargetNum(AdvanceInTargetTrial);
 BlockType=BlockType(AdvanceInTargetTrial);
 
 [UniqueBlockType,orderindex,reprojectorderindex]=unique(BlockType,'row');
 UniqueBlockType=sort(UniqueBlockType,'ascend');%Flexible first, then stable;
NumberOfBlocks=max(reprojectorderindex);
 
 %If the target onset time and fixofftime dimension are not aligned 
 if size(TargetOnTime,2)~=size(FixOffTime,2)
     FixOffTime=FixOffTime(:,1:size(TargetOnTime,2));
     disp('Warning!FixOffTime Corrected!!')
 end
     
 
SearchingTime= (TargetOnTime-FixOffTime)/1000;%Change into s for presentation

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

%First calculate average firing rate for the neuron in this task
ISI=nanmean(diff(SpikeTime,1,2),'all');
AverageFR=1/ISI*1000;%Average FR in the whole tak\s;l]]

%Segregate spike data from Fixation onset to target onset
SegmentData=EventSegTool(SegEvent,SegEventEnd,PlotStartOffset,PlotStopOffset,EventChannel,EventTimeChannel,SpikeTime,StepSize,BinWidth);

%Retrive data from the field
%Total average for each event
MeanFR=SegmentData.MeanFR; %[Mean,sem]%Averall mean and sem
MeanFRRaw=SegmentData.MeanFRRaw;%=AverageFiringRate_Event;

%PSTH(Spike Density Function) for each event
MeanPSTH=SegmentData.MeanPSTH;%=[PSTH_Event_Time',PSTH_Event'];
MeanPSTH_Time=MeanPSTH(:,1);
MeanPSTH_FR=MeanPSTH(:,2);

%Raster and raw spike count PSTH
AlignedRaster=SegmentData.RasterAfterAlign;
PSTH_Event_SpikeCount=SegmentData.PSTH_Event_SpikeCount;


%Relative time bin for each event
RelativeTimeMarker=SegmentData.RelativeTimeMarker;%=EventTimeMarker;
%Absolute time marker for each event
AbosoluteTimeMarker=SegmentData.AbsoluteTimeMarker;


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

%Outcome:Accept/Reject
Outcome=TransformToOutcome(TargetIndex);
Outcome_Line=Outcome;
Outcome_Line=reshape(Outcome_Line,numel(Outcome_Line),1);
Outcome_Line=Outcome_Line(~isnan(Outcome_Line));%1:Accept;0:Reject


%Rearrange scene name into lines
SceneName_Line=repmat(SceneName,1,max(TargetNum));
SceneName_Line=reshape(SceneName_Line,numel(SceneName_Line),1);
SceneName_Line=SceneName_Line(NoNNaNSelector);

%Rearrange block type into lines
BlockType_Line=repmat(BlockType,1,max(TargetNum));
BlockType_Line=reshape(BlockType_Line,numel(BlockType_Line),1);
BlockType_Line=BlockType_Line(NoNNaNSelector);

%TargetOrderLine
TargetOrder=1:max(TargetNum);
TotalNum=length(TargetNum);
TargetOrder_Mat=repmat(TargetOrder,TotalNum,1);
TargetOrder_Line=reshape(TargetOrder_Mat,numel(TargetOrder_Mat),1);
TargetOrder_Line=TargetOrder_Line(NoNNaNSelector);

%Rearrange searching time into lines
SearchingTime_Line=reshape(SearchingTime,numel(SearchingTime),1);
SearchingTime_Line=SearchingTime_Line(NoNNaNSelector);

%Choice Outcome
ChoiceOutcome_Line=repmat(ChoiceOutcome,1,max(TargetNum));
ChoiceOutcome_Line=reshape(ChoiceOutcome_Line,numel(ChoiceOutcome_Line),1);
ChoiceOutcome_Line=ChoiceOutcome_Line(NoNNaNSelector);



%Separate scene types 
SceneWithStableBlock=SceneName(BlockType==11);
SceneWithFlexibleBlock=SceneName(BlockType==1);

UniqueScene_Stable=unique(SceneWithStableBlock);
UniqueScene_Flexible=unique(SceneWithFlexibleBlock);

StableScene=1010:1017;
FlexibleScene=1000:1007;
DuelScene=[1054:1057,1060:1063];

%First: Overall Scene and Object responses

%UniqueScenes=unique(SceneName);
SceneOnEventFlag=ismember(SegEventStr,'SceneOn');
SceneOffEventFlag=ismember(SegEventStr,'SceneOff');

ResponseOn=MeanPSTH_FR{SceneOnEventFlag};
ResponseOff=MeanPSTH_FR{SceneOffEventFlag};
Time_PSTH_SceneOn=nanmean(MeanPSTH_Time{SceneOnEventFlag},1);
Time_PSTH_SceneOff=nanmean(MeanPSTH_Time{SceneOffEventFlag},1);

PSTH_SceneOn_Mean=nanmean(ResponseOn,1);
PSTH_SceneOn_Sem=nanstd(ResponseOn,[],1)./sqrt(sum(~isnan(ResponseOn),1));
%Few number correction
ColNum=sum(~isnan(ResponseOn),1);
RemoveCol=ColNum<=3;
PSTH_SceneOn_Mean(RemoveCol)=NaN;
PSTH_SceneOn_Sem(RemoveCol)=NaN;


PSTH_SceneOff_Mean=nanmean(ResponseOff,1);
PSTH_SceneOff_Sem=nanstd(ResponseOff,[],1)./sqrt(sum(~isnan(ResponseOff),1));

%Few number correction
ColNum=sum(~isnan(ResponseOn),1);
RemoveCol=ColNum<=3;
PSTH_SceneOff_Mean(RemoveCol)=NaN;
PSTH_SceneOff_Sem(RemoveCol)=NaN;



%Find out other markers

SceneOnMarker=RelativeTimeMarker{SceneOnEventFlag};
SceneOffMarker=RelativeTimeMarker{SceneOffEventFlag};

SceneOnMarkerAbs=AbosoluteTimeMarker{SceneOnEventFlag};
SceneOffMarkerAbs=AbosoluteTimeMarker{SceneOffEventFlag};

PreFixTimeAbs=FindOutTime(EventChannel,EventTimeChannel,PREFIX);
PreFixTime=PreFixTimeAbs-SceneOnMarkerAbs(:,1);

SceneOnMarkerMean=nanmean(SceneOnMarker,1);
PreFixTimeMean=nanmean(PreFixTime,1);
if length(PreFixTimeMean)>1
PreFixTimeMean=PreFixTimeMean(PreFixTimeMean<0);
end

SceneOffMarkerMean=nanmean(SceneOffMarker,1);
TgtOfftmp=FindOutTime(EventChannel,EventTimeChannel,TGTOFF);
TgtOffTimeLastAbs=LastValidVal(TgtOfftmp);
TgtOffTimeLast=TgtOffTimeLastAbs-SceneOffMarkerAbs(:,1);
TgtOffTimeLastMean=nanmean(TgtOffTimeLast,1);

TgtOntmp=FindOutTime(EventChannel,EventTimeChannel,TGTCD);
TgtOnTimeLastAbs=LastValidVal(TgtOntmp);
TgtOnTimeLast=TgtOnTimeLastAbs-SceneOffMarkerAbs(:,1);
TgtOnTimeLastMean=nanmean(TgtOnTimeLast,1);



RwdOffAbs=FindOutTime(EventChannel,EventTimeChannel,RWDOFFCD);
RwdOffTime=RwdOffAbs-SceneOffMarkerAbs(:,1);

RwdOnAbs=FindOutTime(EventChannel,EventTimeChannel,RWDCD);
RwdOnTime=RwdOnAbs-SceneOffMarkerAbs(:,1);
RwdOnTimeMean=nanmean(RwdOnTime,1);

%%
%Object Responses
TgtOnEventFlag=ismember(SegEventStr,'TargetOnset');
TgtOffEventFlag=ismember(SegEventStr,'TargetOffset');

TgtOnMarker=RelativeTimeMarker{TgtOnEventFlag};
TgtOffMarker=RelativeTimeMarker{TgtOffEventFlag};

TgtOnMarkerAbs=AbosoluteTimeMarker{TgtOnEventFlag};
TgtOffMarkerAbs=AbosoluteTimeMarker{TgtOffEventFlag};

%%%All Tgt Responses
TgtResponseOn=MeanPSTH_FR{TgtOnEventFlag};
TgtResponseOff=MeanPSTH_FR{TgtOffEventFlag};
Time_PSTH_TgtOn=nanmean(MeanPSTH_Time{TgtOnEventFlag},1);
Time_PSTH_TgtOff=nanmean(MeanPSTH_Time{TgtOffEventFlag},1);


PSTH_TgtOn_Mean=nanmean(TgtResponseOn,1);
PSTH_TgtOn_Sem=nanstd(TgtResponseOn,[],1)./sqrt(sum(~isnan(TgtResponseOn),1));
%Few number correction
ColNum=sum(~isnan(TgtResponseOn),1);
RemoveCol=ColNum<=3;
PSTH_TgtOn_Mean(RemoveCol)=NaN;
PSTH_TgtOn_Sem(RemoveCol)=NaN;



%Separate into good and bad object
GoodFlag=ObjDefine(TargetName_Line);
BadFlag=~GoodFlag;

%Separate into Accept and Reject


GoodAcceptFlag=GoodFlag & Outcome_Line;
BadAcceptFlag=BadFlag & Outcome_Line;

GoodRejectFlag=GoodFlag & (~Outcome_Line);
BadRejectFlag=BadFlag & (~Outcome_Line);


%Only look at the reject tgtoffs
PSTH_TgtOff_Mean=nanmean(TgtResponseOff(~Outcome_Line,:),1);
PSTH_TgtOff_Sem=nanstd(TgtResponseOff(~Outcome_Line,:),[],1)./sqrt(sum(~isnan(TgtResponseOff(~Outcome_Line,:)),1));

%Few number correction
ColNum=sum(~isnan(TgtResponseOff(~Outcome_Line,:)),1);
RemoveCol=ColNum<=3;
PSTH_TgtOff_Mean(RemoveCol)=NaN;
PSTH_TgtOff_Sem(RemoveCol)=NaN;

%Tgt off marker
TgtOffMarkerMean=nanmean(TgtOffMarker,1);

%Only use reject marker
TgtOnMarkerReject=nanmean(TgtOnMarker(~Outcome_Line,:),1);

%CueOnTime
CueOnTimetmp=FindOutTime(EventChannel,EventTimeChannel,CUEON);

if sum(isnan(CueOnTimetmp))~=length(CueOnTimetmp)
CueOnTimetmp_Line=reshape(CueOnTimetmp,numel(CueOnTimetmp),1);

CueOnTimeAbs=CueOnTimetmp_Line(NoNNaNSelector);
else
    CueOnTimeAbs=NaN*ones(sum(NoNNaNSelector),1);
end



CueOnTime=CueOnTimeAbs-TgtOnMarkerAbs(:,1);
CueOnTimeMean=nanmean(CueOnTime);

%Tgton relative to tgt off
TgtOnMarkertmp=TgtOnMarkerAbs(:,1)-TgtOffMarkerAbs(:,1);
TgtOnMarkerRel=TgtOnMarkertmp(~Outcome_Line);
TgtOnMarkerRelMean=nanmean(TgtOnMarkerRel);

%Label the next FP onset 
tmp=FindOutTime(EventChannel,EventTimeChannel,SHOWFIXCD );
FPOnTimetmp=[tmp(:,2:end),NaN*ones(size(tmp,1),1)];
FPOnTimetmp_Line=reshape(FPOnTimetmp,numel(FPOnTimetmp),1);
FPOnTimetmp_Line=FPOnTimetmp_Line(NoNNaNSelector);
AfterObjFPOnTime=FPOnTimetmp_Line-TgtOffMarkerAbs(:,1);

AfterObjFPOnTimeMean=nanmean(AfterObjFPOnTime);

%Separate Object Responses into different conditions

GoodAccTgtOn=TgtResponseOn(GoodAcceptFlag,:);
GoodAccTgtOnMean=nanmean(GoodAccTgtOn,1);
GoodAccTgtOnSem=nanstd(GoodAccTgtOn,[],1)./sqrt(sum(~isnan(GoodAccTgtOn),1));

GoodRejTgtOn=TgtResponseOn(GoodRejectFlag,:);
GoodRejTgtOnMean=nanmean(GoodRejTgtOn,1);
GoodRejTgtOnSem=nanstd(GoodRejTgtOn,[],1)./sqrt(sum(~isnan(GoodRejTgtOn),1));

BadRejTgtOn=TgtResponseOn(BadRejectFlag,:);
BadRejTgtOnMean=nanmean(BadRejTgtOn,1);
BadRejTgtOnSem=nanstd(BadRejTgtOn,[],1)./sqrt(sum(~isnan(BadRejTgtOn),1));

BadAccTgtOn=TgtResponseOn(BadAcceptFlag,:);
BadAccTgtOnMean=nanmean(BadAccTgtOn,1);
BadAccTgtOnSem=nanstd(BadAccTgtOn,[],1)./sqrt(sum(~isnan(BadAccTgtOn),1));

%Fix Loc v.s. Rand Loc
StableSceneSel=ismember(SceneName,StableScene) & BlockType==11;
FlexibleSceneSel=ismember(SceneName,FlexibleScene)  & BlockType==1;
DuelSceneStableSel=ismember(SceneName,DuelScene) & BlockType==11;
DuelSceneFlexibleSel=ismember(SceneName,DuelScene) & BlockType==1;

StableSceneOnResponseMean=nanmean(ResponseOn(StableSceneSel,:),1);
StableSceneOnResponseSem=nanstd(ResponseOn(StableSceneSel,:),[],1)./sqrt(sum(~isnan(ResponseOn(StableSceneSel,:))));

FlexibleSceneOnResponseMean=nanmean(ResponseOn(FlexibleSceneSel,:),1);
FlexibleSceneOnResponseSem=nanstd(ResponseOn(FlexibleSceneSel,:),[],1)./sqrt(sum(~isnan(ResponseOn(FlexibleSceneSel,:))));

DuelStableSceneOnResponseMean=nanmean(ResponseOn(DuelSceneStableSel,:),1);
DuelStableSceneOnResponseSem=nanstd(ResponseOn(DuelSceneStableSel,:),[],1)./sqrt(sum(~isnan(ResponseOn(DuelSceneStableSel,:))));

DuelFlexibleSceneOnResponseMean=nanmean(ResponseOn(DuelSceneFlexibleSel,:),1);
DuelFlexibleSceneOnResponseSem=nanstd(ResponseOn(DuelSceneFlexibleSel,:),[],1)./sqrt(sum(~isnan(ResponseOn(DuelSceneFlexibleSel,:))));


%{
% Loop over the blocks
for i=1:NumberOfBlocks
    
  
    %Organize the raster, spikecount PSTH, spike density function PSTH
    for j=1:length(AlignedRaster)
        if EventTypeFlag(j)==1 %Trial event
            currBlocks_ori=BlockType==UniqueBlockType(i);
            
        elseif EventTypeFlag(j)==2 %Target event
            currBlocks_ori=BlockType_Line==UniqueBlockType(i);
        end
        %Special for reward
        if strcmp( SegEventStr(j),'Reward')
            RewardTrials=ChoiceOutcome==GOODCHOICECD;
            currBlocks=currBlocks_ori & RewardTrials;
        else
            
            currBlocks=currBlocks_ori;
        end
       
       
        %Raster
        Raster_currEvent_Block{i,j}=AlignedRaster{j}(currBlocks,:);
       
        
       
        %Spike count PSTH
        PSTH_Event_SpikeCount_currEvent{i,j}=PSTH_Event_SpikeCount{j}(currBlocks,:);
        
        %Spike density function PSTH
        MeanPSTH_FR_currEvent{i,j}=MeanPSTH_FR{j}(currBlocks,:);
        
        %Timing 
        PSTH_Time_currEvent{i,j}=nanmean(MeanPSTH_Time{j}(currBlocks,:),1);
        RelativeTimeMarker_currEvent{i,j}=RelativeTimeMarker{j}(currBlocks,:);
        
        %Total Average
        MeanFR_currEvent(i,j)=nanmean(MeanFRRaw{j}(currBlocks));
        RawFR_currEvent{i,j}=MeanFRRaw{j}(currBlocks);
        
        %{
        StartLoc=1;
        Barlength=1;
        LineWidth=1;
        SpikeRaterPlot(Raster_currEvent_Block,StartLoc,Barlength,[1,0,0],LineWidth);
        %}
    end
    %SegEventStr={'ExtraFPOn','ExtraFPOff','SceneOn','SceneOff','FixOff','CueOn','TargetOnset','TargetOffset','Reward','PostTrialSpon','BeforeTrialSpon'};
     
    
 %{   
    FixOffSel=ismember(SegEventStr,'FixOff');
    PSTH_fixoff= MeanPSTH_FR_currEvent{i,FixOffSel};
    PSTH_fixoff_time{i}=nanmean(PSTH_Time_currEvent{i,FixOffSel},1);
    %}
    Block_Sel=BlockType_Line==UniqueBlockType(i);
    SearchingTime_currBlock=SearchingTime_Line(Block_Sel);
    
   
end
%}

%%Plot range for sceneon
PlotRangeSceneOn=[-1000,1500];
PlotRangeSceneOff=[-1800,1000];
PlotRangeTgtOn=[-800,500];
PlotRangeTgtOff=[-500,1000];

%%
FigureStartNum=110;
FigureIndex=1;
if ShowFigureFlag
   
%%%%%%%%%%%%%%%%%%Next figure: Spike Density Function PSTH
 
 figtitlestr{FigureIndex}='SpikeDensityFunctionPSTH';
 fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[100,100, 1500,800],'Name',figtitlestr{FigureIndex});   
 
 subplot(3,2,1);%Scene Onset
 
 Sel=Time_PSTH_SceneOn>PlotRangeSceneOn(1) & Time_PSTH_SceneOn<PlotRangeSceneOn(2);
 
 shadedErrorBar(Time_PSTH_SceneOn(Sel),PSTH_SceneOn_Mean(Sel),PSTH_SceneOn_Sem(Sel),'lineprops',{'r'},'transparent',1,'patchSaturation',0.5);%PSTH for fix
 
 
  hold on
  %Add markers
  yrange=ylim;
 plot([SceneOnMarkerMean(1),SceneOnMarkerMean(1)],yrange,'-.r');
 plot([SceneOnMarkerMean(2),SceneOnMarkerMean(2)],yrange,'--k');
 plot([PreFixTimeMean,PreFixTimeMean],yrange,'--m');
 xlim(PlotRangeSceneOn);
 
 set(findall(gca, 'Type', 'Line'),'LineWidth',3);
 
 xlabel('Time from scene onset(ms)');

 ylabel('Firing Rates(Hz)');
            
set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
 

subplot(3,2,2); %SceneOffset
Sel=Time_PSTH_SceneOff>PlotRangeSceneOff(1) & Time_PSTH_SceneOff<PlotRangeSceneOff(2);
 
shadedErrorBar(Time_PSTH_SceneOff(Sel),PSTH_SceneOff_Mean(Sel),PSTH_SceneOff_Sem(Sel),'lineprops',{'r'},'transparent',1,'patchSaturation',0.5);%PSTH for fix
 
 
  hold on
  %Add markers
  yrange=ylim;
 plot([SceneOffMarkerMean(1),SceneOffMarkerMean(1)],yrange,'-.k');
 hold on
 plot([TgtOffTimeLastMean,TgtOffTimeLastMean],yrange,'--m');
 plot([RwdOnTimeMean,RwdOnTimeMean],yrange,'--y')
 plot([TgtOnTimeLastMean,TgtOnTimeLastMean],yrange,'--p')

set(findall(gca, 'Type', 'Line'),'LineWidth',3);
 
 xlim(PlotRangeSceneOff);
 
 xlabel('Time from scene offset(ms)');

 ylabel('Firing Rates(Hz)');
            
set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');


%PSTH_TgtOn_Mean
subplot(3,2,3);%Target Onset\
Sel=Time_PSTH_TgtOn>PlotRangeTgtOn(1) & Time_PSTH_TgtOn<PlotRangeTgtOn(2);

shadedErrorBar(Time_PSTH_TgtOn(Sel),PSTH_TgtOn_Mean(Sel),PSTH_TgtOn_Sem(Sel),'lineprops',{'b'},'transparent',1,'patchSaturation',0.5);%PSTH for fix
 
%Add markers
  yrange=ylim;
  hold on
 plot([TgtOnMarkerReject(1),TgtOnMarkerReject(1)],yrange,'-.k');
 plot([TgtOnMarkerReject(2),TgtOnMarkerReject(2)],yrange,'--m');
 plot([CueOnTimeMean,CueOnTimeMean],yrange,'--g');
 
 set(findall(gca, 'Type', 'Line'),'LineWidth',3);
 xlim(PlotRangeTgtOn);
 
 xlabel('Time from target onset(ms)');

 ylabel('Firing Rates(Hz)');
            
set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');

 
 %PSTH_TgtOff_Mean
subplot(3,2,4);%Target Offset\
Sel=Time_PSTH_TgtOff>PlotRangeTgtOff(1) & Time_PSTH_TgtOff<PlotRangeTgtOff(2);

shadedErrorBar(Time_PSTH_TgtOff(Sel),PSTH_TgtOff_Mean(Sel),PSTH_TgtOff_Sem(Sel),'lineprops',{'b'},'transparent',1,'patchSaturation',0.5);%PSTH for fix
 
%Add markers
  yrange=ylim;
  hold on
plot([SceneOffMarkerMean(1),SceneOffMarkerMean(1)],yrange,'-.k');
 hold on
 plot([TgtOffTimeLastMean,TgtOffTimeLastMean],yrange,'--m');
 plot([RwdOnTimeMean,RwdOnTimeMean],yrange,'--y')
 plot([TgtOnTimeLastMean,TgtOnTimeLastMean],yrange,'--p')
 set(findall(gca, 'Type', 'Line'),'LineWidth',3);
 
 xlim(PlotRangeTgtOff);
 
 xlabel('Time from target offset(ms)');

 ylabel('Firing Rates(Hz)');
            
set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');



%Object on responses separate for each conditions
subplot(3,2,5);%Target On for each condition

Sel=Time_PSTH_TgtOn>PlotRangeTgtOn(1) & Time_PSTH_TgtOn<PlotRangeTgtOn(2);





%shadedErrorBar(Time_PSTH_TgtOn(Sel),GoodRejTgtOnMean(Sel),GoodRejTgtOnSem(Sel),'lineprops',{'#EED4A4'},'transparent',1,'patchSaturation',0.5);%PSTH for fix
hold on
%shadedErrorBar(Time_PSTH_TgtOn(Sel),BadAccTgtOnMean(Sel),BadAccTgtOnSem(Sel),'lineprops',{'#BBC7F2'},'transparent',1,'patchSaturation',0.5);%PSTH for fix

shadedErrorBar(Time_PSTH_TgtOn(Sel),GoodAccTgtOnMean(Sel),GoodAccTgtOnSem(Sel),'lineprops',{'r'},'transparent',1,'patchSaturation',0.5);%PSTH for fix
 
shadedErrorBar(Time_PSTH_TgtOn(Sel),BadRejTgtOnMean(Sel),BadRejTgtOnSem(Sel),'lineprops',{'b'},'transparent',1,'patchSaturation',0.5);%PSTH for fix

%Add markers
  yrange=ylim;
  hold on
 plot([TgtOnMarkerReject(1),TgtOnMarkerReject(1)],yrange,'-.k');
 plot([TgtOnMarkerReject(2),TgtOnMarkerReject(2)],yrange,'--k');
 plot([CueOnTimeMean,CueOnTimeMean],yrange,'--k');
 
%legend({'GoodRej','BadAcc','GoodAcc','BadRej'});

legend({'GoodAcc','BadRej'});
 set(findall(gca, 'Type', 'Line'),'LineWidth',3);

xlim(PlotRangeTgtOn);
 
 xlabel('Time from target onset(ms)');

 ylabel('Firing Rates(Hz)');
            
set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');

subplot(3,2,6);%SceneOnResponse for different type of scenes
Sel=Time_PSTH_SceneOn>PlotRangeSceneOn(1) & Time_PSTH_SceneOn<PlotRangeSceneOn(2);
shadedErrorBar(Time_PSTH_SceneOn(Sel),StableSceneOnResponseMean(Sel),StableSceneOnResponseSem(Sel),'lineprops',{'r'},'transparent',1,'patchSaturation',0.5);%PSTH for fix
 
 
  hold on
  shadedErrorBar(Time_PSTH_SceneOn(Sel),FlexibleSceneOnResponseMean(Sel),FlexibleSceneOnResponseSem(Sel),'lineprops',{'b'},'transparent',1,'patchSaturation',0.5);%PSTH for fix
 
hold on
shadedErrorBar(Time_PSTH_SceneOn(Sel),DuelStableSceneOnResponseMean(Sel),DuelStableSceneOnResponseSem(Sel),'lineprops',{'m'},'transparent',1,'patchSaturation',0.5);%PSTH for fix
shadedErrorBar(Time_PSTH_SceneOn(Sel),DuelFlexibleSceneOnResponseMean(Sel),DuelFlexibleSceneOnResponseSem(Sel),'lineprops',{'g'},'transparent',1,'patchSaturation',0.5);%PSTH for fix

  
  
  %Add markers
  yrange=ylim;
 plot([SceneOnMarkerMean(1),SceneOnMarkerMean(1)],yrange,'-.k');
 plot([SceneOnMarkerMean(2),SceneOnMarkerMean(2)],yrange,'--k');
 plot([PreFixTimeMean,PreFixTimeMean],yrange,'--k');
 xlim(PlotRangeSceneOn);
 legend('Sta','Fle','StaDuel','FleDuel');
  
 set(findall(gca, 'Type', 'Line'),'LineWidth',3);
 
 xlabel('Time from scene onset(ms)');

 ylabel('Firing Rates(Hz)');
            
set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
 


    %Add data information
  %  data_msg=sprintf(' %s, SpikeChannel: %d, TrialNum: %d ',FileName,SpikeChannel,TrialNum);
    data_msg=strcat(FileName,' SpikeChannel:',num2str(SpikeChannel),';ChannelID:',num2str(ChannelID),';ClusterID:',num2str(ClusterID),';TrialNum: ',num2str(TrialNum));

    annotation('textbox', [0.17, 0.95, 0.7, 0.04], 'string', data_msg,'FontSize',20,'FontWeight','Bold','LineStyle','none')

    FigureIndex=FigureIndex+1;
    FigureStartNum=FigureStartNum+1; 


   %{ 
   for i=1:length(AlignedRaster)
        MeanPSTH_FR_currEvent=MeanPSTH_FR{i};
         
      %  for j=1:NumberOfBlocks
            
            %Spike-Count PSTH
           
            currBlock=UniqueBlockType(j);
            
            subplot(ceil((length(AlignedRaster)+1)/2),2,i);
            
            
            curr_Response=nanmean(MeanPSTH_FR_currEvent{j,i},1);
            curr_Response_err=nanstd(MeanPSTH_FR_currEvent{j,i},[],1)./sqrt(sum(~isnan(MeanPSTH_FR_currEvent{j,i}),1));
             %Num correction (set number of responses too few to be NaN)
            NumOfValid=sum(~isnan(MeanPSTH_FR_currEvent{j,i}),1);
            RemovalFlag=NumOfValid<3;
            curr_Response(RemovalFlag)=NaN;
            curr_Response_err(RemovalFlag)=NaN;
            
            plot([0,0],[0,max(curr_Response)],'--k','Color','#929498','LineWidth',1);
            hold on
            %plot(PSTH_Time_currEvent{j,i},curr_Response,'-k','Color',ColorForEachMode(AssociateBlock==currBlock,:),'LineWidth',2);
           % bar(PSTH_Time_currEvent{j,i},nanmean(PSTH_Event_SpikeCount_currEvent{j,i},1),'FaceColor',ColorForEachMode(AssociateBlock==currBlock,:));
            shadedErrorBar(PSTH_Time_currEvent{j,i},curr_Response,curr_Response_err,'lineprops',{ColorForEachMode(AssociateBlock==currBlock,:)},'transparent',1,'patchSaturation',0.5);%PSTH for fix
             set(findall(gca, 'Type', 'Line'),'LineWidth',3);
            hold on 
            %Add event-end marker
            plot([min(RelativeTimeMarker_currEvent{j,i}(:,2)), min(RelativeTimeMarker_currEvent{j,i}(:,2))],[0,max(curr_Response)],'-k','Color','#707276','LineWidth',1.5);
           % SegEventStr={'ExtraFPOn','ExtraFPOff','SceneOn','SceneOff','FixOff','CueOn','TargetOnset','TargetOffset','Reward','Spon'};
          %SegEventEnd=[EXTRAFIXOFF,SHOWFIXCD,SHOWFIXCD,NaN,CUEON,TGTCD,TGTOFF,NaN,RWDOFFCD,NaN ];%Event end; %Automatic screening later
          %Add significance
          if ~isnan(h_Event_PSTH{i})
              sig=area(PSTH_Time_currEvent{j,i},h_Event_PSTH{i}*0.5,'FaceColor','#FFC107','EdgeColor','#FFC107','FaceAlpha',0.5,'EdgeAlpha',0.5);
              set( get( get( sig, 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );
          end
  
         % h_Event_PSTH
            if isnan(SegEventEnd(i))
                PlotStopOffsetPlot(i)=PlotStopOffset(i)/5;
            else
                PlotStopOffsetPlot(i)=PlotStopOffset(i);
            end
            xlim([-PlotStartOffset(i),PlotStopOffsetPlot(i)+nanmean(RelativeTimeMarker_currEvent{j,i}(:,2))]);
            
            %Set ylim
            %Set ylim
            Interval=curr_Response(:,PSTH_Time_currEvent{j,i}<PlotStopOffsetPlot(i)+nanmean(RelativeTimeMarker_currEvent{j,i}(:,2)) & PlotStopOffsetPlot(i)+nanmean(RelativeTimeMarker_currEvent{j,i}(:,2))>-PlotStartOffset(i));
              Interval_max(j)=max( Interval)*1.12;
            
           
            if i==3
             ylabel('Spike Count per s(Hz)');
            end
            if i==length(AlignedRaster)
                xlabel('Time from event onset');
            end
             title(SegEventStr{i});
            set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
            
       % end
        
         ylim([0,max(max(  Interval_max),1)]);
       
   end
   %Add spon with block
   PlotIndex=1:length(SponResponse);
    subplot(ceil((length(AlignedRaster)+1)/2),2,length(AlignedRaster)+1);
    for i=1:NumberOfBlocks
      currBlock=UniqueBlockType(i);
      Sel=BlockType==currBlock;
      plot(PlotIndex(Sel),SponResponse(Sel),'-k','Color',ColorForEachMode(AssociateBlock==currBlock,:),'LineWidth',3);
      hold on
   
    end

set(gca,'XTick',PlotIndex(1:50:end));
set(gca,'XTickLabel',TrialIndex(1:50:end));
xlabel('Trial Index');
ylabel('Spon');
title('PostTrialSpon with trials');
set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');


  
    %Add data information
  %  data_msg=sprintf(' %s, SpikeChannel: %d, TrialNum: %d ',FileName,SpikeChannel,TrialNum);
    data_msg=strcat(FileName,' SpikeChannel:',num2str(SpikeChannel),';ChannelID:',num2str(ChannelID),';ClusterID:',num2str(ClusterID),';TrialNum: ',num2str(TrialNum),celltype);

    annotation('textbox', [0.17, 0.95, 0.7, 0.04], 'string', data_msg,'FontSize',20,'FontWeight','Bold','LineStyle','none')

    FigureIndex=FigureIndex+1;
    FigureStartNum=FigureStartNum+1; 
  %}  



end %End of the if ShowFigureFlag



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
OutputFigureName=sprintf('%s_%s_N%s_C%s_%CID_%NID',MonkeyName,string(RecordDateOriginal),string(NeuronNum),string(ChannelNum),string(ChannelID),string(ClusterID));



ExistFlag=0;
%Load the old file if exist
if exist(OutputFileName)

    load(OutputFileName);
    %Load OutputData into the memory
    if isfield(OutputFileName,'SceneObjSummary')
    ExistFlag=1;
    end
end

%Store figures in common folders
ProtocolMarker='SceneObjSummary';
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
%FileName=strcat(FileName,'_RearchingTime.mat');
if SaveFlag
save(FileName,'OutputData');
disp('Seaching Time Data Exported');
else
    disp('No data export');
end
%}

%{
%%%%%%%%%%%%%%%%%%%%OutputData For current analysis
OutputData_New=[];

OutputData_New.SceneObjSummary.TrialType='GoodTrials';

%%%%%%%%%Save overall scene response characteristics%%%%%%%%%%%%%%%%%%%%%%%
%Organize the data
NameStr={'NaN'};
%{
NameStr={'AverageFR','AverageSpon','AverageSceneDiffBlock','p_SceneDiffBetweenBlocks',...
    'PSTH_SceneOn_FleSta_TimeMeanSem',...
    'PSTH_SceneOff_FleSta_TimeMeanSem',...
    'SI_FleSta','p_AnovaSceneOnAllFleSta',...
  'SceneOn500MeanMinusSpon','p_SceneOn500withSpon','SceneOnMeanMinusSponWhole','p_SceneOnWithSponWhole',...
  'pSceneOnWhole','SceneOnDifference_WholeMean','pSceneOn500','SceneOnDifference_500Mean',...
  'SD_RelSceneOn','Mean_RelSceneOn','UniqueScene'...
  'NumOfStableSigScenes','NumOfFlexibleSigScenes','StableSigScenes','FlexibleSigScenes'...
  'PSTH_MostPreferedStableScene_MeanSem','PSTH_MostPreferedFlexibleScene_MeanSem',...
  'PSTH_OtherStableScene_MeanSem','PSTH_OtherFlexibleScene_MeanSem',...
  'R2WholeFleSta','p_SceneOnAnova_WholeStaFle','pSceneDiffWithTime','rsquare_psth',...
  'pSceneDiffWithTime_Stable','rsquare_psth_Stable','pSceneDiffWithTime_Flexible','rsquare_psth_Flexible',...
  'PropOfSig_Stable','PropOfSig_Flexible',...
  'p_SceneOnAnova_StableHalf','p_SceneOnAnova_FlexibleHalf',...
  'R2_StableHalf','R2_FlexibleHalf',...
    
 
};
%}

DataLib={NaN};
%{
DataLib={AverageFR,AverageSpon,SceneOnAverage',p_SceneDiffBetweenBlocks,...
    [Time_PSTH_SceneOn',MeanPSTH_FR_SceneOn_Flexible',MeanPSTH_FR_SceneOn_Stable',SemPSTH_FR_SceneOn_Flexible',SemPSTH_FR_SceneOn_Stable'],...
     [Time_PSTH_SceneOff',MeanPSTH_FR_SceneOff_Flexible',MeanPSTH_FR_SceneOff_Stable',SemPSTH_FR_SceneOff_Flexible',SemPSTH_FR_SceneOff_Stable'],...
   [SI_Stable,SI_Flexible],[p_SceneOnAnova,p_SceneOnAnova_Stable,p_SceneOnAnova_Flexible],...
SceneOn500MeanMinusSpon,p_SceneOn500WithSpon,SceneOnMeanMinusSponWhole,p_SceneOnWithSponWhole,...
pSceneOnWhole,SceneOnDifference_WholeMean,pSceneOn500,SceneOnDifference_500Mean,...
SD_RelSceneOn,Mean_RelSceneOn,UniqueScenes...
NumOfStableSigScenes,NumOfFlexibleSigScenes,StableSigScenes,FlexibleSigScenes,...
[PSTH_MostPreferedStableScene_Mean',PSTH_MostPreferedStableScene_Sem'],[PSTH_MostPreferedFlexibleScene_Mean',PSTH_MostPreferedFlexibleScene_Sem'],...
[PSTH_OtherStableScene_Mean',PSTH_OtherStableScene_Sem'],[PSTH_OtherFlexibleScene_Mean',PSTH_OtherFlexibleScene_Sem'],...
[R2_Whole,R2_Flexible,R2_Stable],[p_SceneOnAnova,p_SceneOnAnova_Stable,p_SceneOnAnova_Flexible],pSceneDiffWithTime,rsquare_psth,...
pSceneDiffWithTime_Stable,rsquare_psth_Stable,pSceneDiffWithTime_Flexible,rsquare_psth_Flexible,...
 PropOfSig_Stable,PropOfSig_Flexible,...
 [p_SceneOnAnova_StableHalf1,p_SceneOnAnova_StableHalf2],[p_SceneOnAnova_FlexibleHalf1,p_SceneOnAnova_FlexibleHalf2],...
 [R2_StableHalf1,R2_StableHalf2],[R2_FlexibleHalf1,R2_FlexibleHalf2]
};
%}


DataStamp=containers.Map(NameStr,DataLib);

OutputData_New.SceneObjSummary.Task=TaskType;
OutputData_New.SceneObjSummary.TaskCode=TaskCode;
OutputData_New.SceneObjSummary.DataStamp=DataStamp;
%}



%{
if ExistFlag
%Compare the old one with the new one
Task_Old=OutputData.SceneObjSummary.TaskCode;

     if Task_Old~=TaskCode
        %Not the same task, add the new dataset into the old one
        OutputData.SceneObjSummary(numel(OutputData.SceneObjSummary)+1)=OutputData_New.SceneObjSummary;
     else
         %If the same task,the same protocol,replace the old one with the new one
         
         OutputData.SceneObjSummary=OutputData_New.SceneObjSummary;
         
    
    
     end
else
    %Output the current file
   OutputData.SceneObjSummary=OutputData_New.SceneObjSummary; 
end
if OutputFlag
    cd(OutPath);
    save(OutputFileName,'OutputData');
    disp('Data Exported');
else
    disp('No data export');
end  

%}

end %End of the loop of the spike channel

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
function val=LastValidVal(data);
%Find the last nonNaN data
for i=1:size(data,1)
    tmp=data(i,:);
    try
    val(i)=tmp(find(~isnan(tmp)==1,1,'last'));
    catch
        val(i)=NaN;
    end
end
val=val';
end

function index=ObjDefine(Name);
%Function to define good and bad objects
Set=floor(Name/10);
num=mod(Name,Set);
index=num<4;


end