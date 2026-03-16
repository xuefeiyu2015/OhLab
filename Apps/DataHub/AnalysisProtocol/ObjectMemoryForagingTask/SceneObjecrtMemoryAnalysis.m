function SceneObjecrtMemoryAnalysis(Data, StartTrial, EndTrial,ShowFigureFlag,OutputFlag);
%Analysis for the object working memory task;
%Xuefei Yu 05142021
ECode;
%Parameter to define the gaussian window to measure firing rates
StepSize=5;%in ms
BinWidth=15;%in ms

%Currently focus on the scene and object responses
SegEvent=[SCENE_ON,TGTCD,SCENE_OFF,FIXOFF];%Event Flow of interest
%SegEventEnd=[SHOWFIXCD,FIXOFF,TGTOFF,RWDOFFCD ,NaN,EXTRAFP];%Event end
SegEventEnd=[SHOWFIXCD,TGTOFF,NaN,NaN];%Event end
SegEventStr={'SceneOn','TgtOn','SceneOff','FixOff'};

PlotStartOffset=[1000,1000,1800,2000];%in ms
PlotStopOffset=[1000,1000,1500,2000];%in ms

%Indicator for scene and objects 
PlotIndividualSceneOn=0;
PlotIndividualObjectOn=0;


%Color Definition for each events

Color=[lines(6);prism(6);pink(6)];
%Color(2,:)=copper(1)+0.5;
Color=[Color;hsv(10)];
Color(18,:)=[0.1,0.457,0.194];

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




EventChannel=Data.EventChannel(Selection,:);
EventTimeChannel=Data.EventTimeChannel(Selection,:);
EventBin=Data.EventBin(Selection,:);

%Automatically chanage the marker

StartCheck=FindOutTime(EventChannel,EventBin,SegEvent(SegEvent==SCENE_ON));
EndCheck=FindOutTime(EventChannel,EventBin,SegEventEnd(SegEvent==SCENE_ON));
EndCheck=EndCheck(:,1);
DiffCheck=EndCheck-StartCheck;
if sum(DiffCheck<100)>ceil(length(StartCheck)/2)
    SegEventEnd(SegEvent==SCENE_ON)=[FIXOFF];%Event end;change from showfix cd to fix off for contigency
    disp('Change to FixOff for FPOn First Mode');
    
end


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
             0.25 0.80 0.54 
             0.47 0.25 0.80         
             1.00 0.54 0.00  %1
             0.63 0.34 0.94
             0.24 0.74 0.88
             0.37 0.9 0.20 %100
             0.5 0.20 0.84
             ];
%AssociateBlock=[10,11,0,1];
%Blockstr={'FixLoc,FixSceneSeq','FixLoc,RandSceneSeq','RandLoc,FixSceneSeq','RandLoc,RandSceneSeq'};

AssociateBlock=[10,11,0,1,110,111,100,101];
Blockstr={'F-SceneLoc,R-Obj','F-Loc,R-SceneObj','F-Scene,R-LocObj,','R-LocSceneObj','F-ObjSceneLoc','F-LocObj,R-Scene','F-ObjScene,R-Loc','F-ObjLoc,R-Scene'};


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


BlockType=BlockType(Selection,2);

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


%Overview of the performance
NumCorrectTrials=sum(ChoiceOutcome==GOODCHOICECD);
NumErrorTrials=sum(ChoiceOutcome==BADCHOICECD);
TotalTrials=length(ChoiceOutcome);
disp(['CorrectTrials=',num2str(NumCorrectTrials),', ErrorTrials=',num2str(NumErrorTrials),', TotalTrials=',num2str(TotalTrials)]);



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

MeanPSTHEach=SegmentData.MeanPSTHEach;

MeanPSTHEach_Time=MeanPSTHEach(:,1);
MeanPSTHEach_FR=MeanPSTHEach(:,2);

%Raster and raw spike count PSTH
AlignedRaster=SegmentData.RasterAfterAlign;
PSTH_Event_SpikeCount=SegmentData.PSTH_Event_SpikeCount;


%Relative time bin for each event
RelativeTimeMarker=SegmentData.RelativeTimeMarker;%=EventTimeMarker;
%Absolute time marker for each event
AbosoluteTimeMarker=SegmentData.AbsoluteTimeMarker;

%First focus on the scene responses 

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


%%
%%Objects
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

%Matched ObjectIndex Line
MatchedIndex=TargetName==TargetName(:,1);
MatchedIndex_Line=reshape(MatchedIndex,numel(MatchedIndex),1);
MatchedIndex_Line=MatchedIndex_Line(NoNNaNSelector);




%Choice Outcome
ChoiceOutcome_Line=repmat(ChoiceOutcome,1,max(TargetNum));
ChoiceOutcome_Line=reshape(ChoiceOutcome_Line,numel(ChoiceOutcome_Line),1);
ChoiceOutcome_Line=ChoiceOutcome_Line(NoNNaNSelector);


%Incorrect trials 
%Either Reject match and accept non match happend in a trial
RejectMatchTrial=(Outcome==0) &MatchedIndex;
RejectMatchTrial=RejectMatchTrial(:,2:end);
RejectMatchTrialIndex=sum(RejectMatchTrial,2)>=1;
%Both Reject 
%Outcome
IncorrectIndex=RejectMatchTrialIndex | ChoiceOutcome==BADCHOICECD;

%Next is the objects
TgtOnEventFlag=ismember(SegEventStr,'TgtOn');


TgtOnMarker=RelativeTimeMarker{TgtOnEventFlag};


TgtOnMarkerAbs=AbosoluteTimeMarker{TgtOnEventFlag};

%Tgt Off marker
TgtOffTimetmp=FindOutTime(EventChannel,EventTimeChannel,TGTOFF);
TgtOnTimetmp=FindOutTime(EventChannel,EventTimeChannel,TGTCD);
TgtOffTimeRel=TgtOffTimetmp-TgtOnTimetmp;

SampleTgtOff=TgtOffTimeRel(:,1);

SampleTgtOff_Mean=nanmean(SampleTgtOff);

%%%All Tgt Responses
TgtResponseOn=MeanPSTH_FR{TgtOnEventFlag};

Time_PSTH_TgtOn=nanmean(MeanPSTH_Time{TgtOnEventFlag},1);

%%%Separate into sample objects and target objects
SampleObjectIndex=TargetOrder_Line==1;
TargetObjectIndex=TargetOrder_Line>1;

MatchedTargetIndex=MatchedIndex_Line & TargetObjectIndex;
NonMatchedTargetIndex=~MatchedIndex_Line & TargetObjectIndex;


%Separate test according to the choice
%Outcome_Line
RejectMatchTarget=MatchedTargetIndex &(Outcome_Line==0);
RejectNonMatchTarget=NonMatchedTargetIndex &(Outcome_Line==0);

AcceptMatchTarget=MatchedTargetIndex &(Outcome_Line==1);
AcceptNonMatchTarget=NonMatchedTargetIndex &(Outcome_Line==1);

disp(['Number of RejMatch',num2str(sum(RejectMatchTarget)),'; Number of AccNonmatch',num2str(sum(AcceptNonMatchTarget))])





SampleResponseOn=TgtResponseOn(SampleObjectIndex,:);
TestResponseOn=TgtResponseOn(TargetObjectIndex,:);



PSTH_SampleOn_Mean=nanmean(SampleResponseOn,1);
PSTH_SampleOn_Sem=nanstd(SampleResponseOn,[],1)./sqrt(sum(~isnan(SampleResponseOn),1));

%Few number correction
ColNum=sum(~isnan(SampleResponseOn),1);
RemoveCol=ColNum<=3;
PSTH_SampleOn_Mean(RemoveCol)=NaN;
PSTH_SampleOn_Sem(RemoveCol)=NaN;

%Sample for correct trials and incorrect trials 
CorrectTrials=ChoiceOutcome==GOODCHOICECD;
%InCorrectTrials=ChoiceOutcome==BADCHOICECD;


CorrectSampleResponseOn=SampleResponseOn(CorrectTrials,:);
PSTH_CorrectSampleOn_Mean=nanmean(CorrectSampleResponseOn,1);
PSTH_CorrectSampleOn_Sem=nanstd(CorrectSampleResponseOn,[],1)./sqrt(sum(~isnan(CorrectSampleResponseOn),1));
%Few number correction
ColNum=sum(~isnan(CorrectSampleResponseOn),1);
RemoveCol=ColNum<=3;
PSTH_CorrectSampleOn_Mean(RemoveCol)=NaN;
PSTH_CorrectSampleOn_Sem(RemoveCol)=NaN;


InCorrectSampleResponseOn=SampleResponseOn(IncorrectIndex,:);
PSTH_InCorrectSampleOn_Mean=nanmean(InCorrectSampleResponseOn,1);
PSTH_InCorrectSampleOn_Sem=nanstd(InCorrectSampleResponseOn,[],1)./sqrt(sum(~isnan(InCorrectSampleResponseOn),1));
%Few number correction
ColNum=sum(~isnan(InCorrectSampleResponseOn),1);
RemoveCol=ColNum<=3;
PSTH_InCorrectSampleOn_Mean(RemoveCol)=NaN;
PSTH_InCorrectSampleOn_Sem(RemoveCol)=NaN;





PSTH_TestOn_Mean=nanmean(TestResponseOn,1);
PSTH_TestOn_Sem=nanstd(TestResponseOn,[],1)./sqrt(sum(~isnan(TestResponseOn),1));

%Few number correction
ColNum=sum(~isnan(TestResponseOn),1);
RemoveCol=ColNum<=3;
PSTH_TestOn_Mean(RemoveCol)=NaN;
PSTH_TestOn_Sem(RemoveCol)=NaN;

%Compared matched test with nonmatched test
MatchedResponseOn=TgtResponseOn(MatchedTargetIndex,:);
NonMatchedResponseOn=TgtResponseOn(NonMatchedTargetIndex,:);

PSTH_MatchOn_Mean=nanmean(MatchedResponseOn,1);
PSTH_MatchOn_Sem=nanstd(MatchedResponseOn,[],1)./sqrt(sum(~isnan(MatchedResponseOn),1));

%Few number correction
ColNum=sum(~isnan(MatchedResponseOn),1);
RemoveCol=ColNum<=3;
PSTH_MatchOn_Mean(RemoveCol)=NaN;
PSTH_MatchOn_Sem(RemoveCol)=NaN;


PSTH_NonMatchOn_Mean=nanmean(NonMatchedResponseOn,1);
PSTH_NonMatchOn_Sem=nanstd(NonMatchedResponseOn,[],1)./sqrt(sum(~isnan(NonMatchedResponseOn),1));

%Few number correction
ColNum=sum(~isnan(NonMatchedResponseOn),1);
RemoveCol=ColNum<=3;
PSTH_NonMatchOn_Mean(RemoveCol)=NaN;
PSTH_NonMatchOn_Sem(RemoveCol)=NaN;

%Compare responses different with choices
%Reject Match
RejMatchResponseOn=TgtResponseOn(RejectMatchTarget,:);
RejNonMatchResponseOn=TgtResponseOn(RejectNonMatchTarget,:);

RejMatchResponseOn_Mean=nanmean(RejMatchResponseOn,1);
RejMatchResponseOn_Sem=nanstd(RejMatchResponseOn,[],1)./sqrt(sum(~isnan(RejMatchResponseOn)));
%Few number correction
ColNum=sum(~isnan(RejMatchResponseOn),1);
RemoveCol=ColNum<=3;
RejMatchResponseOn_Mean(RemoveCol)=NaN;
RejMatchResponseOn_Sem(RemoveCol)=NaN;



RejNonMatchResponseOn_Mean=nanmean(RejNonMatchResponseOn,1);
RejNonMatchResponseOn_Sem=nanstd(RejNonMatchResponseOn,[],1)./sqrt(sum(~isnan(RejNonMatchResponseOn)));
%Few number correction
ColNum=sum(~isnan(RejNonMatchResponseOn),1);
RemoveCol=ColNum<=3;
RejMatchResponseOn_Mean(RemoveCol)=NaN;
RejMatchResponseOn_Sem(RemoveCol)=NaN;



%AcceptNonmatch
AccMatchResponseOn=TgtResponseOn(AcceptMatchTarget,:);
AccNonMatchResponseOn=TgtResponseOn(AcceptNonMatchTarget,:);

AccMatchResponseOn_Mean=nanmean(AccMatchResponseOn,1);
AccMatchResponseOn_Sem=nanstd(AccMatchResponseOn,[],1)./sqrt(sum(~isnan(AccMatchResponseOn)));
%Few number correction
ColNum=sum(~isnan(AccMatchResponseOn),1);
RemoveCol=ColNum<=3;
AccMatchResponseOn_Mean(RemoveCol)=NaN;
AccMatchResponseOn_Sem(RemoveCol)=NaN;



AccNonMatchResponseOn_Mean=nanmean(AccNonMatchResponseOn,1);
AccNonMatchResponseOn_Sem=nanstd(AccNonMatchResponseOn,[],1)./sqrt(sum(~isnan(AccNonMatchResponseOn)));
%Few number correction
ColNum=sum(~isnan(AccNonMatchResponseOn),1);
RemoveCol=ColNum<=3;
AccNonMatchResponseOn_Mean(RemoveCol)=NaN;
AccNonMatchResponseOn_Sem(RemoveCol)=NaN;







%CueOnTime+FirstFPOnTimefor sample
CueOnTimetmp=FindOutTime(EventChannel,EventTimeChannel,CUEON);

if sum(isnan(CueOnTimetmp))~=length(CueOnTimetmp)
CueOnTimetmp_Line=reshape(CueOnTimetmp,numel(CueOnTimetmp),1);
RemoveNan=isnan(CueOnTimetmp_Line);
CueOnTimeAbs=CueOnTimetmp_Line(~RemoveNan);
else
    CueOnTimeAbs=NaN*ones(sum(NoNNaNSelector),1);
end



CueOnTime=CueOnTimeAbs-TgtOnMarkerAbs(TargetObjectIndex,1);
CueOnTimeMean=nanmean(CueOnTime);


%Loop over the block type
UniqueBlockType=unique(BlockType_Line);

for i=1:length(UniqueBlockType)
    currBlock=UniqueBlockType(i);
    BlockIndex(i)=find(AssociateBlock==currBlock);
    
    selscene=BlockType==currBlock;
    selobj=BlockType_Line==currBlock;
    
    %For Scene
    SceneResponsesCurrBlock=ResponseOn(selscene,:);
    
    tmpmean=nanmean(SceneResponsesCurrBlock,1);
    tmpsem=nanstd(SceneResponsesCurrBlock,[],1)./sqrt(sum(~isnan(SceneResponsesCurrBlock)));
    %Fewer number correction
    SceneResponsesCurrBlock_Mean{i}=FewerCorrection(SceneResponsesCurrBlock,tmpmean,3);
     SceneResponsesCurrBlock_Sem{i}=FewerCorrection(SceneResponsesCurrBlock,tmpsem,3);
     
     %For Object
     %Sample objects
     SampleResponseCurrBlock=TgtResponseOn(SampleObjectIndex&selobj,:);
    
    tmpmean=nanmean(SampleResponseCurrBlock,1);
    tmpsem=nanstd(SampleResponseCurrBlock,[],1)./sqrt(sum(~isnan(SampleResponseCurrBlock)));
    %Fewer number correction
    SampleResponsesCurrBlock_Mean{i}=FewerCorrection(SampleResponseCurrBlock,tmpmean,3);
     SampleResponsesCurrBlock_Sem{i}=FewerCorrection(SampleResponseCurrBlock,tmpsem,3);
    
     %Test objects
     %correct matched objects
    MatchedResponseCurrBlock=TgtResponseOn(AcceptMatchTarget&selobj,:);
    
    tmpmean=nanmean(MatchedResponseCurrBlock,1);
    tmpsem=nanstd(MatchedResponseCurrBlock,[],1)./sqrt(sum(~isnan(MatchedResponseCurrBlock)));
    %Fewer number correction
    MatchedResponsesCurrBlock_Mean{i}=FewerCorrection(MatchedResponseCurrBlock,tmpmean,3);
     MatchedResponsesCurrBlock_Sem{i}=FewerCorrection(MatchedResponseCurrBlock,tmpsem,3);
     
      %correct reject nonmatch objects
    NonMatchedResponseCurrBlock=TgtResponseOn(RejectNonMatchTarget&selobj,:);
    
    tmpmean=nanmean(NonMatchedResponseCurrBlock,1);
    tmpsem=nanstd(NonMatchedResponseCurrBlock,[],1)./sqrt(sum(~isnan(NonMatchedResponseCurrBlock)));
    %Fewer number correction
    NonMatchedResponsesCurrBlock_Mean{i}=FewerCorrection(NonMatchedResponseCurrBlock,tmpmean,3);
     NonMatchedResponsesCurrBlock_Sem{i}=FewerCorrection(NonMatchedResponseCurrBlock,tmpsem,3);
     
     
     
    
end

%Loop over scenes and objects
UniqueScene=unique(SceneName);
UniqueBlockType=unique(BlockType);

for i=1:length(UniqueScene)
    selscene=SceneName==UniqueScene(i);
    
   
    for j=1:length(UniqueBlockType)
        
        BlockType_currScene=UniqueBlockType(j);
        sel_block=BlockType==BlockType_currScene;
        
        sel=selscene & sel_block;
        sceneresponse_curr=ResponseOn(sel,:);
        
        tmpmean=nanmean(sceneresponse_curr,1);
        tmpsem=nanstd(sceneresponse_curr,[],1)./sqrt(sum(~isnan(sceneresponse_curr)));
        
         SceneOn_Mean_Each_Block{i,j}=FewerCorrection(sceneresponse_curr,tmpmean,3);
         SceneOn_Sem_Each_Block{i,j}=FewerCorrection(sceneresponse_curr,tmpsem,3);
         
         
   
        
    end
    
    
    
end

%Loop over the objects

%Loop over the sample objects
SampleObjectName=TargetName(:,1);

UniqueSample=unique(SampleObjectName);


%     SceneCurrBlock=SceneName(selscene);
%     uniqueSceneCurrBlock=unique()
for i=1:length(UniqueSample)
    
    for j=1:length(UniqueBlockType)
        
        sel=SampleObjectName==UniqueSample(i) & BlockType==UniqueBlockType(j);
        
        rsp=SampleResponseOn(sel,:);
        
        tmpmean=nanmean(rsp,1);
        tmpsem=nanstd(rsp,[],1)./sqrt(sum(~isnan(rsp)));
        
         SampleOn_Mean_Each_Block{i,j}=FewerCorrection(rsp,tmpmean,2);
         SampleOn_Sem_Each_Block{i,j}=FewerCorrection(rsp,tmpsem,2);
        
        
    end
    
    
    
    
    
end

%Loop over test objects

TestObjectName=TargetName_Line(TargetObjectIndex);

UniqueTest=unique(TestObjectName);


%     SceneCurrBlock=SceneName(selscene);
%     uniqueSceneCurrBlock=unique()
for i=1:length(UniqueTest)
    
    for j=1:length(UniqueBlockType)
        
        sel=TargetName_Line==UniqueTest(i) & BlockType_Line==UniqueBlockType(j);
        
        rsp=TgtResponseOn(sel,:);
        
        tmpmean=nanmean(rsp,1);
        tmpsem=nanstd(rsp,[],1)./sqrt(sum(~isnan(rsp)));
        
         TestOn_Mean_Each_Block{i,j}=FewerCorrection(rsp,tmpmean,2);
         TestOn_Sem_Each_Block{i,j}=FewerCorrection(rsp,tmpsem,2);
        
        
    end
    
    
    
    
    
end
%%
%%Other period of focus
%The period of fixation while cue showed up

FixOffEventFlag=ismember(SegEventStr,'FixOff');

 ResponseFixOff=MeanPSTH_FR{FixOffEventFlag};

Time_PSTH_FixOff=MeanPSTH_Time{FixOffEventFlag};

NumOfFixOff=max(TargetOrder_Line)-1;
FixOff_Linescreen=TargetOrder_Line-1;
FixOff_Line=FixOff_Linescreen(FixOff_Linescreen>0);


%Find out the markers for fix off


FixOffMarkerAbs=FindOutTime(EventChannel,EventTimeChannel,FIXOFF);

ObjectOntmp=FindOutTime(EventChannel,EventTimeChannel,TGTCD);
ObjectOfftmp=FindOutTime(EventChannel,EventTimeChannel,TGTOFF);

FixStartMarkerAbstmp=FindOutTime(EventChannel,EventTimeChannel,FIXSTARTCD);
FixStartMarkerAbs=FixStartMarkerAbstmp(:,2:end);

%FixStartMarkerRelFixOff=FixStartMarkerAbs-FixOffMarkerAbs;


SampleOnRelFixoff=ObjectOntmp(:,1)-FixOffMarkerAbs(:,1);
SampleOffRelFixoff=ObjectOfftmp(:,1)-FixOffMarkerAbs(:,1);

SampleOnRelFixoff_Mean=nanmean(SampleOnRelFixoff);
SampleOffRelFixoff_Mean=nanmean(SampleOffRelFixoff);

FirstTestOnRelFixoff=ObjectOntmp(:,2)-FixOffMarkerAbs(:,1);
FirstTestOffRelFixoff=ObjectOfftmp(:,2)-FixOffMarkerAbs(:,1);

FirstTestOnRelFixoff_Mean=median(FirstTestOnRelFixoff,'omitnan' );
FirstTestOffRelFixoff_Mean=median(FirstTestOffRelFixoff,'omitnan' );



NumOfFixOff_ori=NumOfFixOff;

%Loop over the time for fix-off
for j=1:NumOfFixOff_ori

  

Time_PSTH_FixOffEach{j}=nanmean(Time_PSTH_FixOff(FixOff_Line==j,:),1);

for i=1:length(UniqueBlockType)
    
     fixOffNumEachBlock(i)=NumOfFixOff;
     fixOffNumEachBlock_Match(i)=NumOfFixOff;
     fixOffNumEachBlock_NonMatch(i)=NumOfFixOff;
        
        sel_block= BlockType_Line==UniqueBlockType(i);
        sel=sel_block(FixOff_Linescreen>0) & (FixOff_Line==j);
        
    
        rsp=ResponseFixOff(sel,:);
        
        tmpmean=nanmean(rsp,1);
        tmpsem=nanstd(rsp,[],1)./sqrt(sum(~isnan(rsp)));
        
        if sum(sel)>=3
         FixOff_Mean_Each_EachRoundBlock{i,j}=FewerCorrection(rsp,tmpmean,3);
         FixOff_Sem_Each_EachRoundBlock{i,j}=FewerCorrection(rsp,tmpsem,3);
        
        
         
        else
            
          fixOffNumEachBlock(i)=fixOffNumEachBlock(i)-1;
           
        end
        
        
        
        %Separate into match and nonmatch
        sel_block_match=sel_block & MatchedTargetIndex;
        sel_block_nonmatch=sel_block & NonMatchedTargetIndex;
        
        sel_match_tmp=sel_block_match(FixOff_Linescreen>0)& (FixOff_Line==j);
        sel_nonmatch_tmp=sel_block_nonmatch(FixOff_Linescreen>0)& (FixOff_Line==j);
        
        rsp_match=ResponseFixOff(sel_match_tmp,:);
        rsp_nonmatch=ResponseFixOff(sel_nonmatch_tmp,:);
        
         tmpmean_match=nanmean(rsp_match,1);
         tmpsem_match=nanstd(rsp_match,[],1)./sqrt(sum(~isnan(rsp_match)));
         
          tmpmean_nonmatch=nanmean(rsp_nonmatch,1);
         tmpsem_nonmatch=nanstd(rsp_nonmatch,[],1)./sqrt(sum(~isnan(rsp_nonmatch)));
         
         
        if sum(sel_match_tmp)>=3
         FixOff_Mean_Match_EachRoundBlock{i,j}=FewerCorrection(rsp_match,tmpmean_match,3);
         FixOff_Sem_Match_EachRoundBlock{i,j}=FewerCorrection(rsp_match,tmpsem_match,3);
        else
             fixOffNumEachBlock_Match(i)=fixOffNumEachBlock_Match(i)-1;
    
        end
         
         if sum(sel_match_tmp)>=3
         FixOff_Mean_NonMatch_EachRoundBlock{i,j}=FewerCorrection(rsp_nonmatch,tmpmean_nonmatch,3);
         FixOff_Sem_NonMatch_EachRoundBlock{i,j}=FewerCorrection(rsp_nonmatch,tmpsem_nonmatch,3);
         else
              fixOffNumEachBlock_NonMatch(i)=fixOffNumEachBlock_NonMatch(i)-1;
    
        end
         
        
        
        
end


 %FixStartMarkerRelFixOff_Mean(j)=nanmean(FixStartMarkerRelFixOff(j,:));

   if j>1
       %For the test
       ObjOnRelTgt=ObjectOntmp(:,j+1)-FixOffMarkerAbs(:,j);
       ObjOffRelTgt=ObjectOfftmp(:,j+1)-FixOffMarkerAbs(:,j);
       
       ObjOnRelTgtMean(j)=median(ObjOnRelTgt,'omitnan' );
       
       ObjOffRelTgtMean(j)=median(ObjOffRelTgt,'omitnan' );
       
   end
end


MaxNumOfFix=max(fixOffNumEachBlock);

%Last analysis
%Fix off before matched object

%MatchedTargetIndex_FixOff=MatchedTargetIndex(FixOff_Linescreen>0);
%NonMatchedTargetIndex_FixOff=NonMatchedTargetIndex(FixOff_Linescreen>0);

%ResponseFixOff_MatchTest=ResponseFixOff(MatchedTargetIndex_FixOff,:);

%ResponseFixOff_NonMatchTest=ResponseFixOff(NonMatchedTargetIndex_FixOff,:);

for i=1:length(UniqueBlockType)
    
       sel_block_match=(BlockType_Line==UniqueBlockType(i))&(MatchedTargetIndex);
       sel_block_nonmatch=(BlockType_Line==UniqueBlockType(i))&(NonMatchedTargetIndex);
       
       
       sel_match=sel_block_match(FixOff_Linescreen>0);
       sel_nonmatch=sel_block_nonmatch(FixOff_Linescreen>0);
    
        rsp=ResponseFixOff(sel_match,:);
        
         rsp_nonmatch=ResponseFixOff(sel_nonmatch,:);
     
        
        tmpmean=nanmean(rsp,1);
        tmpsem=nanstd(rsp,[],1)./sqrt(sum(~isnan(rsp)));
        
        ResponseFixOff_MatchTest_Mean{i}=FewerCorrection(rsp,tmpmean,3);
        
        ResponseFixOff_MatchTest_Sem{i}=FewerCorrection(rsp,tmpsem,3);
        
      %  ResponseFixOff_NonMatchTest_Mean{i}=FewerCorrection(rsp,tmpsem,3);
      
      tmpmean=nanmean(rsp_nonmatch,1);
        tmpsem=nanstd(rsp_nonmatch,[],1)./sqrt(sum(~isnan(rsp_nonmatch)));
        
        ResponseFixOff_NonMatchTest_Mean{i}=FewerCorrection(rsp_nonmatch,tmpmean,3);
        
        ResponseFixOff_NonMatchTest_Sem{i}=FewerCorrection(rsp_nonmatch,tmpsem,3);
    
    
    
    
end




%%Plot range for sceneon
PlotRangeSceneOn=[-1000,2000];
PlotRangeSceneOff=[-1800,1000];
PlotRangeTgtOn=[-900,1000];
PlotRangeTgtOff=[-500,1000];
PlotRangeFixOff=[-1000,2500];




%Time_SamplePSTH_FixOff
%SamplFixOff_Mean_Each_Block



%%
FigureStartNum=110;
FigureIndex=1;
if ShowFigureFlag
   
%%%%%%%%%%%%%%%%%%Next figure: Spike Density Function PSTH
 
 figtitlestr{FigureIndex}='SpikeDensityFunctionPSTH';
 fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[50,100, 1500,800],'Name',figtitlestr{FigureIndex});   
 
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
% plot([TgtOffTimeLastMean,TgtOffTimeLastMean],yrange,'--m');
% plot([RwdOnTimeMean,RwdOnTimeMean],yrange,'--y')
% plot([TgtOnTimeLastMean,TgtOnTimeLastMean],yrange,'--p')

set(findall(gca, 'Type', 'Line'),'LineWidth',3);
 
 xlim(PlotRangeSceneOff);
 
 xlabel('Time from scene offset(ms)');

 ylabel('Firing Rates(Hz)');
            
set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');


%
subplot(3,2,3);%Sample and TestObjects
Sel=Time_PSTH_TgtOn>PlotRangeTgtOn(1) & Time_PSTH_TgtOn<PlotRangeTgtOn(2);
%Sample
shadedErrorBar(Time_PSTH_TgtOn(Sel),PSTH_SampleOn_Mean(Sel),PSTH_SampleOn_Sem(Sel),'lineprops',{'r'},'transparent',1,'patchSaturation',0.5);%PSTH for fix
 hold on
 %Test
 shadedErrorBar(Time_PSTH_TgtOn(Sel),PSTH_TestOn_Mean(Sel),PSTH_TestOn_Sem(Sel),'lineprops',{'b'},'transparent',1,'patchSaturation',0.5);%PSTH for fix

%Add markers
  yrange=ylim;
  hold on
  plot([0,0],yrange,'-.k');
 plot([SampleTgtOff_Mean,SampleTgtOff_Mean],yrange,'--m');
 
 plot([CueOnTimeMean,CueOnTimeMean],yrange,'--g');
 
 set(findall(gca, 'Type', 'Line'),'LineWidth',3);
 xlim(PlotRangeTgtOn);
 legend({'Sample','Test'});
 xlabel('Time from target onset(ms)');

 ylabel('Firing Rates(Hz)');
            
set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');

 %
 %TestObjects
subplot(3,2,4);
Sel=Time_PSTH_TgtOn>PlotRangeTgtOn(1) & Time_PSTH_TgtOn<PlotRangeTgtOn(2);

shadedErrorBar(Time_PSTH_TgtOn(Sel),PSTH_MatchOn_Mean(Sel),PSTH_MatchOn_Sem(Sel),'lineprops',{'r'},'transparent',1,'patchSaturation',0.5);%PSTH for fix
 hold on
shadedErrorBar(Time_PSTH_TgtOn(Sel),PSTH_NonMatchOn_Mean(Sel),PSTH_NonMatchOn_Sem(Sel),'lineprops',{'b'},'transparent',1,'patchSaturation',0.5);%PSTH for fix
 
 
%Add markers
  yrange=ylim;
  hold on
%plot([SceneOffMarkerMean(1),SceneOffMarkerMean(1)],yrange,'-.k');
 %hold on
% plot([TgtOffTimeLastMean,TgtOffTimeLastMean],yrange,'--m');
% plot([RwdOnTimeMean,RwdOnTimeMean],yrange,'--y')
% plot([TgtOnTimeLastMean,TgtOnTimeLastMean],yrange,'--p')
 set(findall(gca, 'Type', 'Line'),'LineWidth',3);
 legend({'Match','NonMatch'});
 
 xlim(PlotRangeTgtOn);
 
 xlabel('Time from target onset(ms)');

 ylabel('Firing Rates(Hz)');
            
set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');




%RejMatchResponseOn_Mean

%
%Object on responses separate for each conditions
subplot(3,2,5);%Target On for each condition

Sel=Time_PSTH_TgtOn>PlotRangeTgtOn(1) & Time_PSTH_TgtOn<PlotRangeTgtOn(2);





shadedErrorBar(Time_PSTH_TgtOn(Sel),RejMatchResponseOn_Mean(Sel),RejMatchResponseOn_Sem(Sel),'lineprops',{'#EED4A4'},'transparent',1,'patchSaturation',0.5);%PSTH for fix
hold on
shadedErrorBar(Time_PSTH_TgtOn(Sel),RejNonMatchResponseOn_Mean(Sel),RejNonMatchResponseOn_Sem(Sel),'lineprops',{'b'},'transparent',1,'patchSaturation',0.5);%PSTH for fix

shadedErrorBar(Time_PSTH_TgtOn(Sel),AccMatchResponseOn_Mean(Sel),AccMatchResponseOn_Sem(Sel),'lineprops',{'r'},'transparent',1,'patchSaturation',0.5);%PSTH for fix
 
shadedErrorBar(Time_PSTH_TgtOn(Sel),AccNonMatchResponseOn_Mean(Sel),AccNonMatchResponseOn_Sem(Sel),'lineprops',{'#BBC7F2'},'transparent',1,'patchSaturation',0.5);%PSTH for fix

%Add markers
  yrange=ylim;
  hold on
% plot([TgtOnMarkerReject(1),TgtOnMarkerReject(1)],yrange,'-.k');
% plot([TgtOnMarkerReject(2),TgtOnMarkerReject(2)],yrange,'--k');
% plot([CueOnTimeMean,CueOnTimeMean],yrange,'--k');
 
%legend({'GoodRej','BadAcc','GoodAcc','BadRej'});

legend({'RejMatch','RejNonMatch','AccMatch','AccNonMatch'});
 set(findall(gca, 'Type', 'Line'),'LineWidth',3);

xlim(PlotRangeTgtOn);
 
 xlabel('Time from target onset(ms)');

 ylabel('Firing Rates(Hz)');
            
set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
%
subplot(3,2,6);%SceneOnResponse for different type of scenes
Sel=Time_PSTH_SceneOn>PlotRangeTgtOn(1) & Time_PSTH_SceneOn<PlotRangeTgtOn(2);

 shadedErrorBar(Time_PSTH_TgtOn(Sel),PSTH_CorrectSampleOn_Mean(Sel),PSTH_CorrectSampleOn_Sem(Sel),'lineprops',{'r'},'transparent',1,'patchSaturation',0.5);%PSTH for fix
hold on

 shadedErrorBar(Time_PSTH_TgtOn(Sel),PSTH_InCorrectSampleOn_Mean(Sel),PSTH_InCorrectSampleOn_Sem(Sel),'lineprops',{'b'},'transparent',1,'patchSaturation',0.5);%PSTH for fix

xlim(PlotRangeTgtOn);
  
  %Add markers
  yrange=ylim;
  hold on
  h=plot([0,0],yrange,'-.k');
  h.Annotation.LegendInformation.IconDisplayStyle = 'off';
  
 h=plot([SampleTgtOff_Mean,SampleTgtOff_Mean],yrange,'--m');
 h.Annotation.LegendInformation.IconDisplayStyle = 'off';
 
 
  legend({'CorrectSample','InCorrectSample'});
 set(findall(gca, 'Type', 'Line'),'LineWidth',3);
 
 xlabel('Time from Sample onset(ms)');

 ylabel('Firing Rates(Hz)');
            
set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
 

%}
    %Add data information
  %  data_msg=sprintf(' %s, SpikeChannel: %d, TrialNum: %d ',FileName,SpikeChannel,TrialNum);
    data_msg=strcat(FileName,' SpikeChannel:',num2str(SpikeChannel),';ChannelID:',num2str(ChannelID),';ClusterID:',num2str(ClusterID),';TrialNum: ',num2str(TrialNum));

    annotation('textbox', [0.17, 0.95, 0.7, 0.04], 'string', data_msg,'FontSize',20,'FontWeight','Bold','LineStyle','none')

    FigureIndex=FigureIndex+1;
    FigureStartNum=FigureStartNum+1; 


 %Next figure; focus on the objects response compare in different blocks
 %SampleResponsesCurrBlock_Mean
  figtitlestr{FigureIndex}='SceneObjectResponsesInEachBlock';
 fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[100,100, 1500,800],'Name',figtitlestr{FigureIndex});   
 
 subplot(2,2,1);%Sample Objects
 
    
 Sel=Time_PSTH_TgtOn>PlotRangeTgtOn(1) & Time_PSTH_TgtOn<PlotRangeTgtOn(2);
%Sample
for i=1:length(UniqueBlockType)
    % SceneResponsesCurrBlock_Mean
    shadedErrorBar(Time_PSTH_TgtOn(Sel),SampleResponsesCurrBlock_Mean{i}(Sel),SampleResponsesCurrBlock_Sem{i}(Sel),'lineprops',{ColorForEachMode(BlockIndex(i),:)},'transparent',1,'patchSaturation',0.5);%PSTH for fix
hold on
end

xlim(PlotRangeTgtOn);
 legend(Blockstr{BlockIndex});
 
 
 %Add markers
  yrange=ylim;
  hold on
  h=plot([0,0],yrange,'-.k');
  h.Annotation.LegendInformation.IconDisplayStyle = 'off';
  
 h=plot([SampleTgtOff_Mean,SampleTgtOff_Mean],yrange,'--m');
 h.Annotation.LegendInformation.IconDisplayStyle = 'off';
 
 %h=plot([CueOnTimeMean,CueOnTimeMean],yrange,'--g');
 %h.Annotation.LegendInformation.IconDisplayStyle = 'off';
 
 set(findall(gca, 'Type', 'Line'),'LineWidth',3);
 
 xlabel('Time from sample object onset(ms)');

 ylabel('Firing Rates(Hz)');
            
set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');

subplot(2,2,2);%Matched Objects
 
    
 Sel=Time_PSTH_TgtOn>PlotRangeTgtOn(1) & Time_PSTH_TgtOn<PlotRangeTgtOn(2);
%Sample
for i=1:length(UniqueBlockType)
    % SceneResponsesCurrBlock_Mean
    shadedErrorBar(Time_PSTH_TgtOn(Sel),MatchedResponsesCurrBlock_Mean{i}(Sel),MatchedResponsesCurrBlock_Sem{i}(Sel),'lineprops',{ColorForEachMode(BlockIndex(i),:)},'transparent',1,'patchSaturation',0.5);%PSTH for fix
hold on
end

xlim(PlotRangeTgtOn);
 legend(Blockstr{BlockIndex});
 
 
 %Add markers
  yrange=ylim;
  hold on
  h=plot([0,0],yrange,'-.k');
  h.Annotation.LegendInformation.IconDisplayStyle = 'off';
  
 %h=plot([SampleTgtOff_Mean,SampleTgtOff_Mean],yrange,'--m');
 %h.Annotation.LegendInformation.IconDisplayStyle = 'off';
 
 h=plot([CueOnTimeMean,CueOnTimeMean],yrange,'--g');
 h.Annotation.LegendInformation.IconDisplayStyle = 'off';
 
 set(findall(gca, 'Type', 'Line'),'LineWidth',3);
 
 xlabel('Time from correct matched object onset(ms)');

 ylabel('Firing Rates(Hz)');
            
set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');

subplot(2,2,3);%NonMatched Objects
 
    
 Sel=Time_PSTH_TgtOn>PlotRangeTgtOn(1) & Time_PSTH_TgtOn<PlotRangeTgtOn(2);
%Sample
for i=1:length(UniqueBlockType)
    % SceneResponsesCurrBlock_Mean
    shadedErrorBar(Time_PSTH_TgtOn(Sel),NonMatchedResponsesCurrBlock_Mean{i}(Sel),NonMatchedResponsesCurrBlock_Sem{i}(Sel),'lineprops',{ColorForEachMode(BlockIndex(i),:)},'transparent',1,'patchSaturation',0.5);%PSTH for fix
hold on
end

xlim(PlotRangeTgtOn);
 legend(Blockstr{BlockIndex});
 
 
 %Add markers
  yrange=ylim;
  hold on
  h=plot([0,0],yrange,'-.k');
  h.Annotation.LegendInformation.IconDisplayStyle = 'off';
  
 %h=plot([SampleTgtOff_Mean,SampleTgtOff_Mean],yrange,'--m');
 %h.Annotation.LegendInformation.IconDisplayStyle = 'off';
 
 h=plot([CueOnTimeMean,CueOnTimeMean],yrange,'--g');
 h.Annotation.LegendInformation.IconDisplayStyle = 'off';
 
 set(findall(gca, 'Type', 'Line'),'LineWidth',3);
 
 xlabel('Time from reject nonmatched object onset(ms)');

 ylabel('Firing Rates(Hz)');
            
set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');

subplot(2,2,4);%SceneOnResponse for different type of scenes
Sel=Time_PSTH_SceneOn>PlotRangeSceneOn(1) & Time_PSTH_SceneOn<PlotRangeSceneOn(2);


for i=1:length(UniqueBlockType)
    % SceneResponsesCurrBlock_Mean
    shadedErrorBar(Time_PSTH_SceneOn(Sel),SceneResponsesCurrBlock_Mean{i}(Sel),SceneResponsesCurrBlock_Sem{i}(Sel),'lineprops',{ColorForEachMode(BlockIndex(i),:)},'transparent',1,'patchSaturation',0.5);%PSTH for fix
hold on
end


  
  %
  %Add markers
  yrange=ylim;
 h=plot([SceneOnMarkerMean(1),SceneOnMarkerMean(1)],yrange,'-.k');
 h.Annotation.LegendInformation.IconDisplayStyle = 'off';
 h=plot([SceneOnMarkerMean(2),SceneOnMarkerMean(2)],yrange,'--k');
 h.Annotation.LegendInformation.IconDisplayStyle = 'off';
 h=plot([PreFixTimeMean,PreFixTimeMean],yrange,'--k');
 h.Annotation.LegendInformation.IconDisplayStyle = 'off';
%}
 xlim(PlotRangeSceneOn);
 legend(Blockstr{BlockIndex});
  
 set(findall(gca, 'Type', 'Line'),'LineWidth',3);
 
 xlabel('Time from scene onset(ms)');

 ylabel('Firing Rates(Hz)');
            
set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
  %Add data information
  %  data_msg=sprintf(' %s, SpikeChannel: %d, TrialNum: %d ',FileName,SpikeChannel,TrialNum);
    data_msg=strcat(FileName,' SpikeChannel:',num2str(SpikeChannel),';ChannelID:',num2str(ChannelID),';ClusterID:',num2str(ClusterID),';TrialNum: ',num2str(TrialNum));

    annotation('textbox', [0.17, 0.95, 0.7, 0.04], 'string', data_msg,'FontSize',20,'FontWeight','Bold','LineStyle','none');
    
   FigureIndex=FigureIndex+1;
    FigureStartNum=FigureStartNum+1; 

 %%%%%%%%%%%%%%%Next figure; Scene responses in each block for each scene
 if PlotIndividualSceneOn
   %Determine the number of figure according to the number of scenes
   NumberOfScenes=length(UniqueScene);
   %Assign two figures if over 4 scene
   NumOfFig=ceil(NumberOfScenes/4);
   Sel=Time_PSTH_SceneOn>PlotRangeSceneOn(1) & Time_PSTH_SceneOn<PlotRangeSceneOn(2);
   for ii=1:NumOfFig
 %SampleResponsesCurrBlock_Mean
  figtitlestr{FigureIndex}=strcat('EachSceneResponsesInEachBlock',ii);
 fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[50+30*ii,100, 1000,800],'Name',figtitlestr{FigureIndex});   
 %Scene number range
 StartIndex=1+(ii-1)*4;
 EndIndex=min(ii*4,NumberOfScenes);
 
 
 UniqueSceneForPlot=UniqueScene(StartIndex:EndIndex);
 
 
 
 for i=1:length(UniqueSceneForPlot)
     currBlockindex=[];
     subplot(2,2,i)
     for j=1:length(UniqueBlockType)
     
         currmean=SceneOn_Mean_Each_Block{StartIndex+i-1,j};
         currsem=SceneOn_Sem_Each_Block{StartIndex+i-1,j};
         
         
         
         if ~isempty(currmean)
             currBlockindex(j)=find(AssociateBlock==UniqueBlockType(j));
             
          shadedErrorBar(Time_PSTH_SceneOn(Sel),currmean(Sel),currsem(Sel),'lineprops',{ColorForEachMode(currBlockindex(j),:)},'transparent',1,'patchSaturation',0.5);%PSTH for fix
             
          hold on
         end
         
         
     end
     title(strcat('Scene',num2str(UniqueSceneForPlot(i))));
     legend(Blockstr{currBlockindex});
  
 set(findall(gca, 'Type', 'Line'),'LineWidth',3);
 
 xlabel('Time from scene onset(ms)');

 ylabel('Firing Rates(Hz)');
            
set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
     
     
 end
 %Add data information
  %  data_msg=sprintf(' %s, SpikeChannel: %d, TrialNum: %d ',FileName,SpikeChannel,TrialNum);
    data_msg=strcat(FileName,' SpikeChannel:',num2str(SpikeChannel),';ChannelID:',num2str(ChannelID),';ClusterID:',num2str(ClusterID),';TrialNum: ',num2str(TrialNum));

    annotation('textbox', [0.17, 0.95, 0.7, 0.04], 'string', data_msg,'FontSize',20,'FontWeight','Bold','LineStyle','none')
 
    
    
FigureIndex=FigureIndex+1;
    FigureStartNum=FigureStartNum+1; 
    
   end%End of the number if fig
 end %End of the if plot each scene on
 
   %%%%%Next figure: Compare sample between blocks %%%%%%%%%%%%%%%%%
 
if PlotIndividualObjectOn
   %Determine the number of figure according to the number of scenes
   NumberOfSamples=length(UniqueSample);
   %Assign two figures if over 4 scene
   NumOfFig=ceil(NumberOfSamples/12);
   Sel=Time_PSTH_TgtOn>PlotRangeTgtOn(1) & Time_PSTH_TgtOn<PlotRangeTgtOn(2);
   for ii=1:NumOfFig
 %SampleResponsesCurrBlock_Mean
  figtitlestr{FigureIndex}=strcat('SampleObjectResponsesInEachBlock',ii);
 fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[50+30*ii,100, 1000,800],'Name',figtitlestr{FigureIndex});   
 %Scene number range
 StartIndex=1+(ii-1)*12;
 EndIndex=min(ii*12,NumberOfSamples);
 
 
 UniqueSamplesForPlot=UniqueSample(StartIndex:EndIndex);
 
 
 
 for i=1:length(UniqueSamplesForPlot)
     currBlockindex=[];
     NumOfSamples_curr=EndIndex-StartIndex+1;
     subplot(ceil(NumOfSamples_curr/3),3,i)
     for j=1:length(UniqueBlockType)
     
         currmean=SampleOn_Mean_Each_Block{StartIndex+i-1,j};
         currsem=SampleOn_Sem_Each_Block{StartIndex+i-1,j};
         
         
         
         if ~isempty(currmean)
             currBlockindex(j)=find(AssociateBlock==UniqueBlockType(j));
             
          shadedErrorBar(Time_PSTH_TgtOn(Sel),currmean(Sel),currsem(Sel),'lineprops',{ColorForEachMode(currBlockindex(j),:)},'transparent',1,'patchSaturation',0.5);%PSTH for fix
             
          hold on
         end
         
         
     end
     title(strcat('Sample',num2str(UniqueSamplesForPlot(i))));
     legend(Blockstr{currBlockindex});
  
 set(findall(gca, 'Type', 'Line'),'LineWidth',3);
 
 xlabel('Time from sample object onset(ms)');

 ylabel('Firing Rates(Hz)');
            
set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
     
     
 end
 
 %Add data information
  %  data_msg=sprintf(' %s, SpikeChannel: %d, TrialNum: %d ',FileName,SpikeChannel,TrialNum);
    data_msg=strcat(FileName,' SpikeChannel:',num2str(SpikeChannel),';ChannelID:',num2str(ChannelID),';ClusterID:',num2str(ClusterID),';TrialNum: ',num2str(TrialNum));

    annotation('textbox', [0.17, 0.95, 0.7, 0.04], 'string', data_msg,'FontSize',20,'FontWeight','Bold','LineStyle','none')   
    
FigureIndex=FigureIndex+1;
    FigureStartNum=FigureStartNum+1; 
    
   end%End of the number if fig 
   
   
    
   %%%%%Next figure: Compare test between blocks %%%%%%%%%%%%%%%%%
   
   %Determine the number of figure according to the number of scenes
   NumberOfTests=length(UniqueTest);
   %Assign two figures if over 4 scene
   NumOfFig=ceil(NumberOfTests/12);
   Sel=Time_PSTH_TgtOn>PlotRangeTgtOn(1) & Time_PSTH_TgtOn<PlotRangeTgtOn(2);
   for ii=1:NumOfFig
 %SampleResponsesCurrBlock_Mean
  figtitlestr{FigureIndex}=strcat('TestObjectResponsesInEachBlock',ii);
 fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[50+30*ii,100, 1000,800],'Name',figtitlestr{FigureIndex});   
 %Scene number range
 StartIndex=1+(ii-1)*12;
 EndIndex=min(ii*12,NumberOfTests);
 
 
 UniqueTestForPlot=UniqueTest(StartIndex:EndIndex);
 
 
 
 for i=1:length(UniqueTestForPlot)
     currBlockindex=[];
     numtestcurr=EndIndex-StartIndex+1;
     subplot(ceil(numtestcurr/3),3,i);
     for j=1:length(UniqueBlockType)
     
         currmean=TestOn_Mean_Each_Block{StartIndex+i-1,j};
         currsem=TestOn_Sem_Each_Block{StartIndex+i-1,j};
         
         
         
         if ~isempty(currmean)
             currBlockindex(j)=find(AssociateBlock==UniqueBlockType(j));
             
          shadedErrorBar(Time_PSTH_TgtOn(Sel),currmean(Sel),currsem(Sel),'lineprops',{ColorForEachMode(currBlockindex(j),:)},'transparent',1,'patchSaturation',0.5);%PSTH for fix
             
          hold on
         end
         
         
     end
     title(strcat('Test',num2str(UniqueTestForPlot(i))));
     legend(Blockstr{currBlockindex});
  
 set(findall(gca, 'Type', 'Line'),'LineWidth',3);
 
 xlabel('Time from test object onset(ms)');

 ylabel('Firing Rates(Hz)');
            
set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
     
     
 end
 
 %Add data information
  %  data_msg=sprintf(' %s, SpikeChannel: %d, TrialNum: %d ',FileName,SpikeChannel,TrialNum);
    data_msg=strcat(FileName,' SpikeChannel:',num2str(SpikeChannel),';ChannelID:',num2str(ChannelID),';ClusterID:',num2str(ClusterID),';TrialNum: ',num2str(TrialNum));

    annotation('textbox', [0.17, 0.95, 0.7, 0.04], 'string', data_msg,'FontSize',20,'FontWeight','Bold','LineStyle','none')   
    
FigureIndex=FigureIndex+1;
    FigureStartNum=FigureStartNum+1; 
    
   end%End of the number if fig   
   
end %end of the if PlotIndividualObjectOn 
   



%%%%%Next figure:Period from sample fixation off %%%%%%%%%%%%%%%%%
 
 figtitlestr{FigureIndex}='SampleFixOffInEachBlock';
 fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[120,100, 1000,800],'Name',figtitlestr{FigureIndex});   
 
 for j=1:MaxNumOfFix
     subplot(MaxNumOfFix,1,j)
     currBlockindex=[];
     for i=1:length(UniqueBlockType)
         Sel=Time_PSTH_FixOffEach{j}>PlotRangeFixOff(1) & Time_PSTH_FixOffEach{j}<PlotRangeFixOff(2);
         
          %if ~isempty(FixOff_Mean_Each_EachRoundBlock{i,j})
          if j<=size(FixOff_Mean_Each_EachRoundBlock,2) & ~isempty(FixOff_Mean_Each_EachRoundBlock{i,j})
             currBlockindex(i)=find(AssociateBlock==UniqueBlockType(i));
 
         shadedErrorBar(Time_PSTH_FixOffEach{j}(Sel),FixOff_Mean_Each_EachRoundBlock{i,j}(Sel),FixOff_Sem_Each_EachRoundBlock{i,j}(Sel),'lineprops',{ColorForEachMode(currBlockindex(i),:)},'transparent',1,'patchSaturation',0.5);%PSTH for fix
 
         hold on
         
          end
         
         
     end
     
     %Add marker
     if j==1
        % SampleOnRelFixoff_Mean
        %For sample
        yrange=ylim;
         h=plot([SampleOnRelFixoff_Mean,SampleOnRelFixoff_Mean],yrange,'--r');
         h.Annotation.LegendInformation.IconDisplayStyle = 'off';
         
         h=plot([SampleOffRelFixoff_Mean,SampleOffRelFixoff_Mean],yrange,'--b');
         h.Annotation.LegendInformation.IconDisplayStyle = 'off';
         
         h=plot([0,0],yrange,'-.k');
         
         h.Annotation.LegendInformation.IconDisplayStyle = 'off';
         
         h=plot([FirstTestOnRelFixoff_Mean,FirstTestOnRelFixoff_Mean],yrange,'--r');
         h.Annotation.LegendInformation.IconDisplayStyle = 'off';
         h=plot([FirstTestOffRelFixoff_Mean,FirstTestOffRelFixoff_Mean],yrange,'--b');
         h.Annotation.LegendInformation.IconDisplayStyle = 'off';
         
     else
         %Test marker
         
         yrange=ylim;
         h=plot([ObjOnRelTgtMean(j),ObjOnRelTgtMean(j)],yrange,'--r');
         h.Annotation.LegendInformation.IconDisplayStyle = 'off';
         
          h=plot([0,0],yrange,'-.k');
         h.Annotation.LegendInformation.IconDisplayStyle = 'off';
         
         
         h=plot([ObjOffRelTgtMean(j),ObjOffRelTgtMean(j)],yrange,'--b');
         h.Annotation.LegendInformation.IconDisplayStyle = 'off';
         
         
         
     end
     
     
     
     legend(Blockstr{currBlockindex});
  set(findall(gca, 'Type', 'Line'),'LineWidth',3);
 
 xlabel(strcat('Time from ',num2str(j),' fixation off(ms)'));

 ylabel('Firing Rates(Hz)');
            
set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
     
     
 end
 
 %Add data information
  %  data_msg=sprintf(' %s, SpikeChannel: %d, TrialNum: %d ',FileName,SpikeChannel,TrialNum);
    data_msg=strcat(FileName,' SpikeChannel:',num2str(SpikeChannel),';ChannelID:',num2str(ChannelID),';ClusterID:',num2str(ClusterID),';TrialNum: ',num2str(TrialNum));

    annotation('textbox', [0.17, 0.95, 0.7, 0.04], 'string', data_msg,'FontSize',20,'FontWeight','Bold','LineStyle','none')

FigureIndex=FigureIndex+1;
FigureStartNum=FigureStartNum+1; 
    
    
    %%%%%Next figure:Period from fixation off for match and nonmatch%%%%%%%%%%%%%%%%%
 
 figtitlestr{FigureIndex}='FixOffInNonMatchMatchEachBlock';
 fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[140,100, 1000,800],'Name',figtitlestr{FigureIndex});   
 
    for i=1:length(UniqueBlockType)
        subplot(length(UniqueBlockType),1,i);%Matched object
        
        Sel=Time_PSTH_FixOffEach{2}>PlotRangeFixOff(1) & Time_PSTH_FixOffEach{2}<PlotRangeFixOff(2);
        
        
          if ~isempty(ResponseFixOff_MatchTest_Mean{i})
             currBlockindex(i)=find(AssociateBlock==UniqueBlockType(i));
 
         shadedErrorBar(Time_PSTH_FixOffEach{2}(Sel),ResponseFixOff_MatchTest_Mean{i}(Sel),ResponseFixOff_MatchTest_Sem{i}(Sel),'lineprops',{'r'},'transparent',1,'patchSaturation',0.5);%PSTH for fix
      hold on
       shadedErrorBar(Time_PSTH_FixOffEach{2}(Sel),ResponseFixOff_NonMatchTest_Mean{i}(Sel),ResponseFixOff_NonMatchTest_Sem{i}(Sel),'lineprops',{'b'},'transparent',1,'patchSaturation',0.5);%PSTH for fix
     
         title(Blockstr{currBlockindex(i)});
          end
   
    
    legend({'Match','NonMatch'});
  set(findall(gca, 'Type', 'Line'),'LineWidth',3);
 
 xlabel(strcat('Time from fixation off before matched object(ms)'));

 ylabel('Firing Rates(Hz)');
            
set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');

    end  
    %Add data information
  %  data_msg=sprintf(' %s, SpikeChannel: %d, TrialNum: %d ',FileName,SpikeChannel,TrialNum);
    data_msg=strcat(FileName,' SpikeChannel:',num2str(SpikeChannel),';ChannelID:',num2str(ChannelID),';ClusterID:',num2str(ClusterID),';TrialNum: ',num2str(TrialNum));

    annotation('textbox', [0.17, 0.95, 0.7, 0.04], 'string', data_msg,'FontSize',20,'FontWeight','Bold','LineStyle','none')
   
FigureIndex=FigureIndex+1;
FigureStartNum=FigureStartNum+1;  

%%Next figure: Fix off responses for match and nonmatch object in each block
for i=1:length(UniqueBlockType)
    
    currBlockindex(i)=find(AssociateBlock==UniqueBlockType(i));
    
    MaxNumOfFixCurr=max(fixOffNumEachBlock_Match(i),fixOffNumEachBlock_NonMatch(i));
    
 figtitlestr{FigureIndex}=strcat('FixOffForBlock',Blockstr{currBlockindex(i)});
 fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[100+30*i,100, 1000,800],'Name',figtitlestr{FigureIndex});   
 
 
 for j=1:MaxNumOfFixCurr
     subplot(MaxNumOfFixCurr,1,j)
     currBlockindex=[];
     
         Sel=Time_PSTH_FixOffEach{j}>PlotRangeFixOff(1) & Time_PSTH_FixOffEach{j}<PlotRangeFixOff(2);
         %Match
         
          if  j<=size(FixOff_Mean_Match_EachRoundBlock,2) && i<=size(FixOff_Mean_Match_EachRoundBlock,1)
           %  currBlockindex(i)=find(AssociateBlock==UniqueBlockType(i));
           
           if  ~isempty(FixOff_Mean_Match_EachRoundBlock{i,j}) 
 
         shadedErrorBar(Time_PSTH_FixOffEach{j}(Sel),FixOff_Mean_Match_EachRoundBlock{i,j}(Sel),FixOff_Sem_Match_EachRoundBlock{i,j}(Sel),'lineprops',{'r'},'transparent',1,'patchSaturation',0.5);%PSTH for fix
      hold on
 
           end
           
          end
           
          
          %NonMatch
         if  j<=size(FixOff_Mean_NonMatch_EachRoundBlock,2)&& i<=size(FixOff_Mean_NonMatch_EachRoundBlock,1)
          if ~isempty(FixOff_Mean_NonMatch_EachRoundBlock{i,j})
           %  currBlockindex(i)=find(AssociateBlock==UniqueBlockType(i));
 
         shadedErrorBar(Time_PSTH_FixOffEach{j}(Sel),FixOff_Mean_NonMatch_EachRoundBlock{i,j}(Sel),FixOff_Sem_NonMatch_EachRoundBlock{i,j}(Sel),'lineprops',{'b'},'transparent',1,'patchSaturation',0.5);%PSTH for fix
      hold on
         
          end
         end
         
         
     
     
     %Add marker
     if j==1
        % SampleOnRelFixoff_Mean
        %For sample
        yrange=ylim;
         h=plot([SampleOnRelFixoff_Mean,SampleOnRelFixoff_Mean],yrange,'--r');
         h.Annotation.LegendInformation.IconDisplayStyle = 'off';
         
         h=plot([SampleOffRelFixoff_Mean,SampleOffRelFixoff_Mean],yrange,'--b');
         h.Annotation.LegendInformation.IconDisplayStyle = 'off';
         
         h=plot([0,0],yrange,'-.k');
         
         h.Annotation.LegendInformation.IconDisplayStyle = 'off';
         
         h=plot([FirstTestOnRelFixoff_Mean,FirstTestOnRelFixoff_Mean],yrange,'--r');
         h.Annotation.LegendInformation.IconDisplayStyle = 'off';
         h=plot([FirstTestOffRelFixoff_Mean,FirstTestOffRelFixoff_Mean],yrange,'--b');
         h.Annotation.LegendInformation.IconDisplayStyle = 'off';
         
     else
         %Test marker
         
         yrange=ylim;
         h=plot([ObjOnRelTgtMean(j),ObjOnRelTgtMean(j)],yrange,'--r');
         h.Annotation.LegendInformation.IconDisplayStyle = 'off';
         
          h=plot([0,0],yrange,'-.k');
         h.Annotation.LegendInformation.IconDisplayStyle = 'off';
         
         
         h=plot([ObjOffRelTgtMean(j),ObjOffRelTgtMean(j)],yrange,'--b');
         h.Annotation.LegendInformation.IconDisplayStyle = 'off';
         
         
         
     end
     
     
     
     
  set(findall(gca, 'Type', 'Line'),'LineWidth',3);
 
 xlabel(strcat('Time from ',num2str(j),' fixation off(ms)'));

 ylabel('Firing Rates(Hz)');
            
set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
     
     
 end
 
 %Add data information
  %  data_msg=sprintf(' %s, SpikeChannel: %d, TrialNum: %d ',FileName,SpikeChannel,TrialNum);
    data_msg=strcat(FileName,' SpikeChannel:',num2str(SpikeChannel),';ChannelID:',num2str(ChannelID),';ClusterID:',num2str(ClusterID),';TrialNum: ',num2str(TrialNum));

    annotation('textbox', [0.17, 0.95, 0.7, 0.04], 'string', data_msg,'FontSize',20,'FontWeight','Bold','LineStyle','none')
    
 FigureIndex=FigureIndex+1;
FigureStartNum=FigureStartNum+1; 
       

end %End of the block type
 


    
    
end %End of the if ShowFigureFlag





end %End of the if spike channel




end %End of the whole function

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
%Fewernumber correction
function val=FewerCorrection(data,tmp,the)
%Remove the tmo data in which the number is lesser than the
val=tmp;
     ColNum=sum(~isnan(data),1);
RemoveCol=ColNum<=the;
val(RemoveCol)=NaN;
end