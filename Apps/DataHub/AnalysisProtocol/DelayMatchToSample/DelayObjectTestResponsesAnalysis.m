function DelayObjectTestResponsesAnalysis(Data, StartTrial, EndTrial,ShowFigureFlag,OutputFlag);
%Analysis protocol for delay match to sample tasks
%By xuefei yu 05222021
ECode;
SegEvent=[TGTCD,RWDCD  ];%Event Flow of interest
%SegEventEnd=[SHOWFIXCD,FIXOFF,TGTOFF,RWDOFFCD ,NaN,EXTRAFP];%Event end
SegEventEnd=[NaN,NaN];%Event end
SegEventStr={'TargetOn','RewardOn'};
AmpMatchFlag=0;%1:Amplitude match on ;0:Amplitude Match off


PlotStartOffset=[2000,1000];%in ms
PlotStopOffset=[5000,2000];%in ms

PlotRangeSample=[-1000,3000];
PlotRangeTest=[-1000,3000];
PlotRangeRwd=[-1000,2000];


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
PlotSacPreBuffer=500;%in ms
PlotSacPostBuffer=4000;%in ms




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
      
Color_RGB2=[
          136,204,238;
          204,102,119;
          120,94,240;
          100,143,255;
           254,97,0;
           0,158,115;
          ]/255;
      
%Color for each object
Color=[lines(6);prism(6);pink(6)];
%Color(2,:)=copper(1)+0.5;
Color=[Color;hsv(10)];
Color(18,:)=[0.1,0.457,0.194];


ColorG=Color;
ColorB=Color;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Unpack all the data base and task settings
ECode;
FileName=Data.FileName;
TaskType=Data.TaskType;
TaskCode=Data.TaskCode;
TrialNum=EndTrial-StartTrial+1;
Selection=StartTrial:EndTrial;

%SpikeChannel=Data.SpikeChannel;

%{
if strcmp(FileName,'Robin040321FRACSAC1.rst') & SpikeChannel==22
    Selection=Selection(~ismember(Selection,[14,40]));
end
%}
%Trial outcome;
EventChannel=Data.EventChannel(Selection,:);
EventTimeChannel=Data.EventTimeChannel(Selection,:);
EventBin=Data.EventBin(Selection,:);



CorrectTrials=~(sum(ismember(EventChannel,UNMATCHOBJ),2)==1);

disp(['Number of error trials=',num2str(sum(~CorrectTrials))]);

%Target Information
TargetIndex=Data.TargetIndex(Selection,2:end,1);%For this task, there is only target name, but no target index

TargetName=Data.TargetIndex(Selection,2:end,2);%For this task, there is only target name, but no target index

%TargetLoc=Data.TargetLoc(Selection,:,:);

TargetNum=Data.TargetNum(Selection,2:end);

%Separate single target trials with two target trials

TaskType_ori=Data.ParaLib('TaskType');

TaskType=TaskType_ori(Selection,2:end);

TgtOnOffDur_ori=Data.ParaLib('TargetOnDuration,TargetOnDurationVar,TargetOffDuration');
TgtOffDur=TgtOnOffDur_ori(:,end);


SingleTgtTrial=TaskType==1;
ChoiceTrial=TaskType>1;

%SingleTgtTrial=TargetNum==1;
%ChoiceTrial=TargetNum>1;

SampleTargetIndex=TargetIndex(:,1);
SampleTargetName=TargetName(:,1);

%SampleGood=SampleTargetIndex<=4;
%SampleBad=SampleTargetIndex>4;

ObjectIndex=ReproduceFromEvent(EventChannel,[GOODOBJ  BADOBJ])';


SampleGood=(ObjectIndex==GOODOBJ);
SampleBad=(ObjectIndex==BADOBJ);

%For flexible value 
SampleOriGood=SampleTargetIndex<=4;
SampleOriBad=SampleTargetIndex>4;

FlipIndex=(SampleGood==SampleOriBad) | (SampleBad==SampleOriGood);


SingleTgtTrialGood=SingleTgtTrial & SampleGood;
SingleTgtTrialBad=SingleTgtTrial & SampleBad;

ChoiceTrialGood=ChoiceTrial & SampleGood & CorrectTrials;
ChoiceTrialBad=ChoiceTrial & SampleBad & CorrectTrials;

if size(TargetIndex,2)<3
    %All single targets
    TestTargetIndex=[TargetIndex(:,2:end),NaN*ones(size(TargetIndex,1),1)];
else
    TestTargetIndex=TargetIndex(:,2:end);
    
end



%Separate into 4 conditions

TestGood_Ori=TestTargetIndex(:,2)<=4;

TestGood=TestGood_Ori;
TestGood(FlipIndex)=~TestGood(FlipIndex);
TestBad=~TestGood;

GoodGoodTest=TestGood & ChoiceTrial & SampleGood;
SampleGoodTestBad=TestBad & ChoiceTrial & SampleGood;
SampleBadTestGood=TestGood & ChoiceTrial & SampleBad;
BadBadTest=TestBad & ChoiceTrial & SampleBad;

disp(['GoodGoodTest=',num2str(sum(GoodGoodTest)),';SampleGoodTestBad=',num2str(sum(SampleGoodTestBad)),';SampleBadTestGood=',num2str(sum(SampleBadTestGood)),';SampleGoodTestBad=',num2str(sum(BadBadTest))])


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

%PSTH
SegmentData=EventSegTool(SegEvent,SegEventEnd,PlotStartOffset,PlotStopOffset,EventChannel,EventTimeChannel,SpikeTime,StepSize,BinWidth);
%Relative time bin for each event
RelativeTimeMarkerEach=SegmentData.RelativeTimeMarkerEach;%=EventTimeMarker;
%Absolute time marker for each event
AbsoluteTimeMarkerEach=SegmentData.AbsoluteTimeMarkerEach;

ObjOnTimeMarker=RelativeTimeMarkerEach{1};
%ObjOffTimeMarker=nanmean(RelativeTimeMarker{2});

%PSTH for each event
MeanPSTHEach=SegmentData.MeanPSTHEach;%=[PSTH_Event_Time',PSTH_Event'];

%This protocol focus on target onset

%Align on target onset

%PSTH for target onset
PSTH_Test_Time=nanmean(MeanPSTHEach{1}{2},1);
PSTH_Test_Response=MeanPSTHEach{3}{2};


PSTH_GG_Response=PSTH_Test_Response(GoodGoodTest,:);
PSTH_GB_Response=PSTH_Test_Response(SampleGoodTestBad,:);

PSTH_BG_Response=PSTH_Test_Response(SampleBadTestGood,:);
PSTH_BB_Response=PSTH_Test_Response(BadBadTest,:);

[PSTH_GG_Response_Mean,PSTH_GG_Response_Sem]=StatisticsPSTH(PSTH_GG_Response,3);
[PSTH_GB_Response_Mean,PSTH_GB_Response_Sem]=StatisticsPSTH(PSTH_GB_Response,3);
[PSTH_BG_Response_Mean,PSTH_BG_Response_Sem]=StatisticsPSTH(PSTH_BG_Response,3);
[PSTH_BB_Response_Mean,PSTH_BB_Response_Sem]=StatisticsPSTH(PSTH_BB_Response,3);
%Test target on one target

PSTH_GSingleTest_Response=PSTH_Test_Response(SingleTgtTrialGood,:);
PSTH_BSingleTest_Response=PSTH_Test_Response(SingleTgtTrialBad,:);

[PSTH_GSingleTest_Response_Mean,PSTH_GSingleTest_Response_Sem]=StatisticsPSTH(PSTH_GSingleTest_Response,3);
[PSTH_BSingleTest_Response_Mean,PSTH_BSingleTest_Response_Sem]=StatisticsPSTH(PSTH_BSingleTest_Response,3);

RwdMarkerTest=RwdMarkerAbs-TgtOnMarkerAbs(:,2);
RwdMarker_Test_SingleMean=nanmean(RwdMarkerTest(SingleTgtTrial));


TgtOffMarker_Test=TgtOffMarkerAbs(:,2)-TgtOnMarkerAbs(:,2);
TgtOffMarker_Test_Single=TgtOffMarker_Test(SingleTgtTrial);
TgtOffMarker_Test_Single_Mean=nanmean(TgtOffMarker_Test_Single);

%PSTH for reward onset
PSTH_Reward_Time=nanmean(MeanPSTHEach{2}{1},1);
PSTH_Reward_Response=MeanPSTHEach{4}{1};

RewardAmount_ori=Data.ParaLib('RewardAmount');

RewardCondition_ori=Data.ParaLib('RewardAfterTgtOffFlag');

RewardAmount=RewardAmount_ori(:,2);
RewardCondition=RewardCondition_ori(:,2);

RewardAmountCorrect=RewardAmount(CorrectTrials);


RewardAmountUnqiue=uniquetol(RewardAmount,0.1);



%Loop over each sets

SampleSetName=floor(SampleTargetName/10);
UniqueSet=unique(SampleSetName);

for i=1:length(UniqueSet)
    
    currSet=UniqueSet(i);
    sel=SampleSetName==currSet;
    
    SampleGoodInSet=sel& SampleGood & CorrectTrials & ChoiceTrial;%Only looks at the correct trials
    SampleBadInSet=sel& SampleBad & CorrectTrials & ChoiceTrial;%Only looks at the correct trials
    
    SampleGoodInSetSingle=sel& SampleGood & CorrectTrials & SingleTgtTrial;%Only looks at the correct trials
    SampleBadInSetSingle=sel& SampleBad & CorrectTrials & SingleTgtTrial;%Only looks at the correct trials
    
    %Delay length
    DelayLength=TgtOffDur;
    %Only choose long delays;
    DelayChoose=abs(DelayLength-max(TgtOffDur))<0.1;
    
    %{
    
    PSTH_SampleInSet_Good=PSTH_Sample_Response(SampleGoodInSet & DelayChoose,:);
    PSTH_SampleInSet_Bad=PSTH_Sample_Response(SampleBadInSet & DelayChoose,:);
    
    
    PSTH_SampleInSetSingle_Good=PSTH_Sample_Response(SampleGoodInSetSingle & DelayChoose,:);
    PSTH_SampleInSetSingle_Bad=PSTH_Sample_Response(SampleBadInSetSingle & DelayChoose,:);
    
    
    [PSTH_SampleInSet_Good_Mean{i},PSTH_SampleInSet_Good_Sem{i}]=StatisticsPSTH(PSTH_SampleInSet_Good,3);
    [PSTH_SampleInSet_Bad_Mean{i},PSTH_SampleInSet_Bad_Sem{i}]=StatisticsPSTH(PSTH_SampleInSet_Bad,3);
    
    [PSTH_SampleInSetSingle_Good_Mean{i},PSTH_SampleInSetSingle_Good_Sem{i}]=StatisticsPSTH(PSTH_SampleInSetSingle_Good,3);
    [PSTH_SampleInSetSingle_Bad_Mean{i},PSTH_SampleInSetSingle_Bad_Sem{i}]=StatisticsPSTH(PSTH_SampleInSetSingle_Bad,3);
    %}
    
    
    %Loop over each object
    
    SampleObjInSet=SampleTargetName(sel);
    
    SampleObjInSet_Unqiue=unique(SampleObjInSet);
    
    
    for j=1:length(SampleObjInSet_Unqiue)
        
        sel_obj=SampleTargetName==SampleObjInSet_Unqiue(j);
        
        
        PSTH_SampleInSet_Each=PSTH_Sample_Response(sel_obj & DelayChoose & sel & CorrectTrials ,:);
        
        [PSTH_SampleInSet_Each_Mean{i,j},PSTH_SampleInSet_Each_Sem{i,j}]=StatisticsPSTH(PSTH_SampleInSet_Each,3);
        
        ObjName{i,j}=num2str(SampleObjInSet_Unqiue(j));
        
        
        %For flexible value, each object has a good and a bad period
        
        sel_obj_Good=sel_obj & SampleGoodInSet & CorrectTrials;
        sel_obj_Bad=sel_obj & SampleBadInSet & CorrectTrials;
        
        PSTH_SampleInSetGood_Each=PSTH_Sample_Response(sel_obj_Good & DelayChoose ,:);
        PSTH_SampleInSetBad_Each=PSTH_Sample_Response(sel_obj_Bad & DelayChoose ,:);
        
         [ PSTH_SampleInSetGood_Each_Mean{i,j}, PSTH_SampleInSetGood_Each_Sem{i,j}]=StatisticsPSTH(PSTH_SampleInSetGood_Each,3);
         [PSTH_SampleInSetBad_Each_Mean{i,j},PSTH_SampleInSetBad_Each_Sem{i,j}]=StatisticsPSTH(PSTH_SampleInSetBad_Each,3);
        
        
        
        
    end
    
    
    
    
end








end
end

function data=FewNumCorrection(mat,data_ori,the);
data=data_ori;
colnum=sum(~isnan(mat),1);
inavidcol=colnum<the;
data(inavidcol)=NaN;


end
function [mean,sem]=StatisticsPSTH(psth,the);

mean=nanmean(psth,1);
sem=nanstd(psth,[],1)./sqrt(sum(~isnan(psth),1));

mean=FewNumCorrection(psth,mean,the);
sem=FewNumCorrection(psth,sem,the);


end
function ObjectIndex=ReproduceFromEvent(event,code)
for i=1:size(event,1)
    
    ObjectIndex(i)=event(i,ismember(event(i,:),code));
    
  
    
end



end


