function DelayObjectResponsesAnalysis(Data, StartTrial, EndTrial,ShowFigureFlag,OutputFlag);
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
%Abnormal check
CheckCodePair=[TGTCD,TGTOFF];%This code should be the same as target on for good trials
AbnormalTrials=AbnormalTrialCheck(Data,CheckCodePair);
Selection=Selection(~ismember(Selection,AbnormalTrials));
    





if strcmp(FileName,'Robin040321FRACSAC1.rst') & SpikeChannel==22
    Selection=Selection(~ismember(Selection,[14,40]));
end

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
TgtOffDur=TgtOnOffDur_ori(Selection,end);


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


%PSTH for sample onset
PSTH_Sample_Time=nanmean(MeanPSTHEach{1}{1},1);
PSTH_Sample_Response=MeanPSTHEach{3}{1};

%Single target trials
PSTH_Sample_Response_SingleGood=PSTH_Sample_Response(SingleTgtTrialGood,:);
PSTH_Sample_Response_SingleBad=PSTH_Sample_Response(SingleTgtTrialBad,:);



[PSTH_Sample_Response_SingleGood_Mean,PSTH_Sample_Response_SingleGood_Sem]=StatisticsPSTH(PSTH_Sample_Response_SingleGood,3);
[PSTH_Sample_Response_SingleBad_Mean,PSTH_Sample_Response_SingleBad_Sem]=StatisticsPSTH(PSTH_Sample_Response_SingleBad,3);

%Two target trials

PSTH_Sample_Response_ChoiceTrialGood=PSTH_Sample_Response(ChoiceTrialGood,:);
PSTH_Sample_Response_ChoiceTrialBad=PSTH_Sample_Response(ChoiceTrialBad,:);



[PSTH_Sample_Response_ChoiceGood_Mean,PSTH_Sample_Response_ChoiceGood_Sem]=StatisticsPSTH(PSTH_Sample_Response_ChoiceTrialGood,3);
[PSTH_Sample_Response_ChoiceBad_Mean,PSTH_Sample_Response_ChoiceBad_Sem]=StatisticsPSTH(PSTH_Sample_Response_ChoiceTrialBad,3);


%Time marker
ObjectOffMarker=FindOutTime(EventChannel,EventTimeChannel,TGTOFF )-FindOutTime(EventChannel,EventTimeChannel,TGTCD);
SampleOffMarker=ObjectOffMarker(:,1);

SampleOffMarker_Mean=nanmean(SampleOffMarker);
TgtOnMarkerAbs=FindOutTime(EventChannel,EventTimeChannel,TGTCD);
TgtOffMarkerAbs=FindOutTime(EventChannel,EventTimeChannel,TGTOFF);

ChoiceTgtOnMarker=TgtOnMarkerAbs(:,2)-TgtOnMarkerAbs(:,1);
ChoiceTgtOffMarker=TgtOffMarkerAbs(:,2)-TgtOnMarkerAbs(:,1);

ChoiceTgtOnMarker_Mean=nanmean(ChoiceTgtOnMarker);
ChoiceTgtOffMarker_Mean=nanmean(ChoiceTgtOffMarker);

RwdMarkerAbs=FindOutTime(EventChannel,EventTimeChannel,RWDCD);
RwdMarker=RwdMarkerAbs-TgtOnMarkerAbs(:,1);
RwdMarker_Mean=nanmean(RwdMarker);

ChoiceTgtOffMarker_SingleMean=nanmean(ChoiceTgtOffMarker(SingleTgtTrial));
RwdMarker_Mean_SingleMean=nanmean(RwdMarker(SingleTgtTrial));


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

RewardAmount=RewardAmount_ori(Selection,2);
RewardCondition=RewardCondition_ori(Selection,2);

RewardAmountCorrect=RewardAmount(CorrectTrials);


RewardAmountUnqiue=uniquetol(RewardAmount,0.1);

AfterRewardIndex=1;
AfterReward=[];

for i=1:length(RewardAmountUnqiue)
   currReward= RewardAmountUnqiue(i);
   sel=RewardAmount-currReward<0.1;
   %{
   if currReward<0.1
       %This part of response is not included in the reward marker
       %Reproduced from the normal trials
       ErrorTrialResponse=PSTH_Test_Response(sel,:);
       
   end
   %}
   
   if  currReward>0.1
   currResponse=PSTH_Reward_Response(sel(CorrectTrials)&(RewardCondition(CorrectTrials)==0),:);
   [PSTH_Reward_Response_Mean{i},PSTH_Reward_Response_Sem{i}]=StatisticsPSTH(currResponse,3);
   
   %Reward After block
   sel_rwdafter=RewardCondition(CorrectTrials)==1 & sel(CorrectTrials);
   if sum(sel_rwdafter)>0
       
      AfterReward(AfterRewardIndex)=currReward;
     
      
      currResponse=PSTH_Reward_Response(sel_rwdafter,:);
       [PSTH_RewardAfter_Response_Mean{AfterRewardIndex},PSTH_RewardAfter_Response_Sem{AfterRewardIndex}]=StatisticsPSTH(currResponse,3);
      
       AfterRewardIndex=AfterRewardIndex+1;
       
       
   end
   
  
   end
    
    
end

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
    
    
    
    PSTH_SampleInSet_Good=PSTH_Sample_Response(SampleGoodInSet & DelayChoose,:);
    PSTH_SampleInSet_Bad=PSTH_Sample_Response(SampleBadInSet & DelayChoose,:);
    
    
    PSTH_SampleInSetSingle_Good=PSTH_Sample_Response(SampleGoodInSetSingle & DelayChoose,:);
    PSTH_SampleInSetSingle_Bad=PSTH_Sample_Response(SampleBadInSetSingle & DelayChoose,:);
    
    
    [PSTH_SampleInSet_Good_Mean{i},PSTH_SampleInSet_Good_Sem{i}]=StatisticsPSTH(PSTH_SampleInSet_Good,3);
    [PSTH_SampleInSet_Bad_Mean{i},PSTH_SampleInSet_Bad_Sem{i}]=StatisticsPSTH(PSTH_SampleInSet_Bad,3);
    
    [PSTH_SampleInSetSingle_Good_Mean{i},PSTH_SampleInSetSingle_Good_Sem{i}]=StatisticsPSTH(PSTH_SampleInSetSingle_Good,3);
    [PSTH_SampleInSetSingle_Bad_Mean{i},PSTH_SampleInSetSingle_Bad_Sem{i}]=StatisticsPSTH(PSTH_SampleInSetSingle_Bad,3);
    
    
    %Loop over each object
    
    SampleObjInSet=SampleTargetName(sel);
    
    SampleObjInSet_Unqiue=unique(SampleObjInSet);
    
    
    for j=1:length(SampleObjInSet_Unqiue)
        
        sel_obj=SampleTargetName==SampleObjInSet_Unqiue(j);
        
        response_sel=sel_obj & DelayChoose & sel & CorrectTrials & ChoiceTrial;
        
        PSTH_SampleInSet_Each=PSTH_Sample_Response(response_sel,:);
        
        [PSTH_SampleInSet_Each_Mean{i,j},PSTH_SampleInSet_Each_Sem{i,j}]=StatisticsPSTH(PSTH_SampleInSet_Each,3);
        
        ObjName{i,j}=num2str(SampleObjInSet_Unqiue(j));
        
        
        %For flexible value, each object has a good and a bad period
        
        sel_obj_Good=SampleGoodInSet & response_sel ;
        sel_obj_Bad=SampleBadInSet & response_sel ;
        
        PSTH_SampleInSetGood_Each=PSTH_Sample_Response(sel_obj_Good,:);
        PSTH_SampleInSetBad_Each=PSTH_Sample_Response(sel_obj_Bad ,:);
        
         [ PSTH_SampleInSetGood_Each_Mean{i,j}, PSTH_SampleInSetGood_Each_Sem{i,j}]=StatisticsPSTH(PSTH_SampleInSetGood_Each,3);
         [PSTH_SampleInSetBad_Each_Mean{i,j},PSTH_SampleInSetBad_Each_Sem{i,j}]=StatisticsPSTH(PSTH_SampleInSetBad_Each,3);
        
        
        
        
    end
    
    
    
    
end





 StartIndex=1;
 StartFig=100;

if ShowFigureFlag
   
    
    
   figtitlestr{StartIndex}='PSTHFracOnAndFracOff';

     fig{StartIndex}=PrepareFigure(StartFig,'w',[50,100, 1200,800],'Name',figtitlestr{StartIndex});
    
     sel=PSTH_Sample_Time>PlotRangeSample(1) & PSTH_Sample_Time<PlotRangeSample(2);
    subplot(2,3,1);%Sample On two target
     shadedErrorBar(PSTH_Sample_Time( sel),PSTH_Sample_Response_ChoiceGood_Mean( sel),PSTH_Sample_Response_ChoiceGood_Sem( sel),'lineprops',{[1    0    0]},'transparent',1,'patchSaturation',0.3);

    hold on
     shadedErrorBar(PSTH_Sample_Time( sel),PSTH_Sample_Response_ChoiceBad_Mean( sel),PSTH_Sample_Response_ChoiceBad_Sem( sel),'lineprops',{[0    0    1]},'transparent',1,'patchSaturation',0.3);
     
     %Add marker line
     yrange=ylim;
     xlim(PlotRangeSample);
     plot([0,0],yrange,'--k');  
     plot([SampleOffMarker_Mean,SampleOffMarker_Mean],yrange,'--k');
     plot([ChoiceTgtOnMarker_Mean,  ChoiceTgtOnMarker_Mean],yrange,'--r');
     plot([ChoiceTgtOffMarker_SingleMean,  ChoiceTgtOffMarker_SingleMean],yrange,'--r');
     plot([ RwdMarker_Mean,   RwdMarker_Mean],yrange,'--y');
     
     
     
     
     set(findall(gca, 'Type', 'Line'),'LineWidth',3);
set(gca,'FontSize',15,'FontWeight','Bold');
set(gca,'LineWidth',3);

xlabel('Time from Sample Onset(ms)','FontSize',25,'FontWeight','Bold');
ylabel('Spikes/s','FontSize',25,'FontWeight','Bold');
title('Choice Trials');


   subplot(2,3,2);%Sample On one target
     shadedErrorBar(PSTH_Sample_Time( sel),PSTH_Sample_Response_SingleGood_Mean( sel),PSTH_Sample_Response_SingleGood_Sem( sel),'lineprops',{[1    0    0]},'transparent',1,'patchSaturation',0.3);

    hold on
     shadedErrorBar(PSTH_Sample_Time( sel),PSTH_Sample_Response_SingleBad_Mean( sel),PSTH_Sample_Response_SingleBad_Sem( sel),'lineprops',{[0    0    1]},'transparent',1,'patchSaturation',0.3);
     
     %Add marker line
     yrange=ylim;
     xlim(PlotRangeSample);
     plot([0,0],yrange,'--k');  
     plot([SampleOffMarker_Mean,SampleOffMarker_Mean],yrange,'--k');
     plot([ChoiceTgtOnMarker_Mean,  ChoiceTgtOnMarker_Mean],yrange,'--r');
     plot([ChoiceTgtOffMarker_Mean,  ChoiceTgtOffMarker_Mean],yrange,'--r');
     plot([RwdMarker_Mean_SingleMean,   RwdMarker_Mean_SingleMean],yrange,'--y');
     
     
     
     
     set(findall(gca, 'Type', 'Line'),'LineWidth',3);
set(gca,'FontSize',15,'FontWeight','Bold');
set(gca,'LineWidth',3);

xlabel('Time from Sample Onset(ms)','FontSize',25,'FontWeight','Bold');
ylabel('Spikes/s','FontSize',25,'FontWeight','Bold');
title('Single Target Trials');  
  


%
    subplot(2,3,3);%Test On one target
    sel=PSTH_Test_Time>PlotRangeTest(1) & PSTH_Test_Time<PlotRangeTest(2);
    shadedErrorBar(PSTH_Test_Time( sel),PSTH_GSingleTest_Response_Mean( sel),PSTH_GSingleTest_Response_Sem( sel),'lineprops',{[1    0    0]},'transparent',1,'patchSaturation',0.3);

    hold on
     shadedErrorBar(PSTH_Test_Time( sel),PSTH_BSingleTest_Response_Mean( sel),PSTH_BSingleTest_Response_Sem( sel),'lineprops',{[0    0    1]},'transparent',1,'patchSaturation',0.3);
     
    
     %Add marker line
     yrange=ylim;
     xlim(PlotRangeTest);
     plot([0,0],yrange,'--k');  
    %{
    plot([SampleOffMarker_Mean,SampleOffMarker_Mean],yrange,'--k');
     
      plot([ChoiceTgtOnMarker_Mean,  ChoiceTgtOnMarker_Mean],yrange,'--r');
     plot([ChoiceTgtOffMarker_SingleMean,  ChoiceTgtOffMarker_SingleMean],yrange,'--r');
     
     %}
     plot([ RwdMarker_Test_SingleMean,  RwdMarker_Test_SingleMean],yrange,'--y');
      plot([TgtOffMarker_Test_Single_Mean,  TgtOffMarker_Test_Single_Mean],yrange,'--r');
     
     
      set(findall(gca, 'Type', 'Line'),'LineWidth',3);
set(gca,'FontSize',15,'FontWeight','Bold');
set(gca,'LineWidth',3);

xlabel('Time from Test Onset(ms)','FontSize',25,'FontWeight','Bold');
ylabel('Spikes/s','FontSize',25,'FontWeight','Bold');
title('Sinlge Target Trials');
 
    %}
    
    
  subplot(2,3,4);%Test On two target
  sel=PSTH_Test_Time>PlotRangeTest(1) & PSTH_Test_Time<PlotRangeTest(2);
     shadedErrorBar(PSTH_Test_Time( sel),PSTH_GG_Response_Mean( sel),PSTH_GG_Response_Sem( sel),'lineprops',{[1    0    0]},'transparent',1,'patchSaturation',0.3);

    hold on
     shadedErrorBar(PSTH_Test_Time( sel),PSTH_GB_Response_Mean( sel),PSTH_GB_Response_Sem( sel),'lineprops',{[0.5    0    0]},'transparent',1,'patchSaturation',0.3);
     
     %{
     %Add marker line
     yrange=ylim;
     plot([0,0],yrange,'--k');  
     plot([SampleOffMarker_Mean,SampleOffMarker_Mean],yrange,'--k');
     plot([ChoiceTgtOnMarker_Mean,  ChoiceTgtOnMarker_Mean],yrange,'--r');
     plot([ChoiceTgtOffMarker_Mean,  ChoiceTgtOffMarker_Mean],yrange,'--r');
     plot([ RwdMarker_Mean,   RwdMarker_Mean],yrange,'--y');  
    %}
    set(findall(gca, 'Type', 'Line'),'LineWidth',3);
     legend({'GoodGood','GoodBad'});
     %
set(gca,'FontSize',15,'FontWeight','Bold');
set(gca,'LineWidth',3);

xlabel('Time from Test Onset(ms)','FontSize',25,'FontWeight','Bold');
ylabel('Spikes/s','FontSize',25,'FontWeight','Bold');
title('Two Target Trials');
     %}

 subplot(2,3,5);%Test On two target
     shadedErrorBar(PSTH_Test_Time( sel),PSTH_BG_Response_Mean( sel),PSTH_BG_Response_Sem( sel),'lineprops',{[0    0    0.5]},'transparent',1,'patchSaturation',0.3);

    hold on
     shadedErrorBar(PSTH_Test_Time( sel),PSTH_BB_Response_Mean( sel),PSTH_BB_Response_Sem( sel),'lineprops',{[0    0    1]},'transparent',1,'patchSaturation',0.3);
     
     %{
     %Add marker line
     yrange=ylim;
     plot([0,0],yrange,'--k');  
     plot([SampleOffMarker_Mean,SampleOffMarker_Mean],yrange,'--k');
     plot([ChoiceTgtOnMarker_Mean,  ChoiceTgtOnMarker_Mean],yrange,'--r');
     plot([ChoiceTgtOffMarker_Mean,  ChoiceTgtOffMarker_Mean],yrange,'--r');
     plot([ RwdMarker_Mean,   RwdMarker_Mean],yrange,'--y');  
    %}
    set(findall(gca, 'Type', 'Line'),'LineWidth',3);
    legend({'BadGood','BadBad'});
set(gca,'FontSize',15,'FontWeight','Bold');
set(gca,'LineWidth',3);

xlabel('Time from Test Onset(ms)','FontSize',25,'FontWeight','Bold');
ylabel('Spikes/s','FontSize',25,'FontWeight','Bold');
title('Two Target Trials');
 
    
%{
[PSTH_GG_Response_Mean,PSTH_GG_Response_Sem]=StatisticsPSTH(PSTH_GG_Response,3);
[PSTH_GB_Response_Mean,PSTH_GB_Response_Sem]=StatisticsPSTH(PSTH_GB_Response,3);
[PSTH_BG_Response_Mean,PSTH_BG_Response_Sem]=StatisticsPSTH(PSTH_BG_Response,3);
[PSTH_BB_Response_Mean,PSTH_BB_Response_Sem]=StatisticsPSTH(PSTH_BB_Response,3);
%}
    
 subplot(2,3,6);%Reward Onset
 
 rewardafterindex=1;
 order=1:length(AfterReward);
 
 rwd_str={};
 index=1;
 
 for i=1:length(RewardAmountUnqiue)
   currReward= RewardAmountUnqiue(i);
   sel=RewardAmount-currReward<0.1;
   %{
   if currReward<0.1
       %This part of response is not included in the reward marker
       %Reproduced from the normal trials
       ErrorTrialResponse=PSTH_Test_Response(sel,:);
       
   end
   %}
   %
   if  currReward>0.1
  % currResponse=PSTH_Reward_Response(sel(CorrectTrials),:);
  % [PSTH_Reward_Response_Mean{i},PSTH_Reward_Response_Sem{i}]=StatisticsPSTH(currResponse,3);
   sel=PSTH_Reward_Time>PlotRangeRwd(1) & PSTH_Reward_Time<PlotRangeRwd(2);
   shadedErrorBar(PSTH_Reward_Time( sel),PSTH_Reward_Response_Mean{i}(sel),PSTH_Reward_Response_Sem{i}(sel),'lineprops',{Color_RGB(i,:)},'transparent',1,'patchSaturation',0.3);

    hold on
 set(findall(gca, 'Type', 'Line'),'LineWidth',3);
 
 rwd_str{index}=sprintf('BeforeTgt:%1.0d',RewardAmountUnqiue(i));
index=index+1;
   end
   %
   %
   if ismember(currReward,AfterReward)
       shadedErrorBar(PSTH_Reward_Time( sel),PSTH_RewardAfter_Response_Mean{ismember(AfterReward,currReward)}(sel),PSTH_RewardAfter_Response_Sem{ismember(AfterReward,currReward)}(sel),'lineprops',{Color_RGB2(i,:)},'transparent',1,'patchSaturation',0.3);
 
       order(rewardafterindex)=find(ismember(AfterReward,currReward)==1);
       
       rewardafterindex=rewardafterindex+1;
       hold on
       set(findall(gca, 'Type', 'Line'),'LineWidth',3);
      
        rwd_str{index}=sprintf('AfterTgt:%1.0d',RewardAmountUnqiue(i));
    index=index+1;
       
       
   end
   %}
  
    
end

   
 legend(rwd_str) ;
 
  set(findall(gca, 'Type', 'Line'),'LineWidth',3);
 set(gca,'FontSize',15,'FontWeight','Bold');
set(gca,'LineWidth',3);

xlabel('Time from Reward Onset(ms)','FontSize',25,'FontWeight','Bold');
ylabel('Spikes/s','FontSize',25,'FontWeight','Bold');

 
 

  %Add title

title_msg=sprintf(' %s, SpikeChannel: %d, TrialNum: %d ',FileName,SpikeChannel,TrialNum);
a=annotation('textbox',[0.25,0.68,0.3,0.3],'String',title_msg,'FitBoxToText','on');
a.FontSize=25;
a.FontWeight='Bold';
a.LineStyle='none';  


    
end

 StartIndex=StartIndex+1;
 StartFig=StartFig+1;



 figtitlestr{StartIndex}='PSTHSamplePeriodInEachSet';

 fig{StartIndex}=PrepareFigure(StartFig,'w',[60,100, 1200,800],'Name',figtitlestr{StartIndex});
 for i=1:length(UniqueSet)
     subplot(ceil(sqrt(length(UniqueSet))),ceil(sqrt(length(UniqueSet))),i);
     
     sel=PSTH_Sample_Time>PlotRangeSample(1) & PSTH_Sample_Time<PlotRangeSample(2);
    
     shadedErrorBar(PSTH_Sample_Time( sel),PSTH_SampleInSet_Good_Mean{i}( sel),PSTH_SampleInSet_Good_Sem{i}(sel),'lineprops',{[1    0    0]},'transparent',1,'patchSaturation',0.3);

    hold on
     shadedErrorBar(PSTH_Sample_Time( sel),PSTH_SampleInSet_Bad_Mean{i}( sel),PSTH_SampleInSet_Bad_Sem{i}( sel),'lineprops',{[0    0    1]},'transparent',1,'patchSaturation',0.3);
     
     %Add marker line
     yrange=ylim;
     xlim(PlotRangeSample);
     plot([0,0],yrange,'--k');  
     plot([SampleOffMarker_Mean,SampleOffMarker_Mean],yrange,'--k');
     plot([ChoiceTgtOnMarker_Mean,  ChoiceTgtOnMarker_Mean],yrange,'--r');
     plot([ChoiceTgtOffMarker_SingleMean,  ChoiceTgtOffMarker_SingleMean],yrange,'--r');
     plot([ RwdMarker_Mean,   RwdMarker_Mean],yrange,'--y');
     
     
     
     
     set(findall(gca, 'Type', 'Line'),'LineWidth',3);
set(gca,'FontSize',15,'FontWeight','Bold');
set(gca,'LineWidth',3);

xlabel('Time from Sample Onset(ms)','FontSize',25,'FontWeight','Bold');
ylabel('Spikes/s','FontSize',25,'FontWeight','Bold');
title(sprintf('Set %d',UniqueSet(i)));

     
     
 end
 StartIndex=StartIndex+1;
 StartFig=StartFig+1;

 
 %%%%  Next figure: PSTH of sample period for each object
 figtitlestr{StartIndex}='PSTHSamplePeriodInEachSet';

 fig{StartIndex}=PrepareFigure(StartFig,'w',[80,100, 1200,800],'Name',figtitlestr{StartIndex});
 for i=1:length(UniqueSet)
     subplot(ceil(sqrt(length(UniqueSet))),ceil(sqrt(length(UniqueSet))),i);
     
      currSet=UniqueSet(i);
       sel=SampleSetName==currSet;
    
  %   for j=1:
    SampleObjInSet=SampleTargetName(sel);
    
    SampleObjInSet_Unqiue=unique(SampleObjInSet);
    
     for j=1:length(SampleObjInSet_Unqiue)
          sel=PSTH_Sample_Time>PlotRangeSample(1) & PSTH_Sample_Time<PlotRangeSample(2);  
          shadedErrorBar(PSTH_Sample_Time( sel),PSTH_SampleInSet_Each_Mean{i,j}( sel),PSTH_SampleInSet_Each_Sem{i,j}(sel),'lineprops',{Color(j,:)},'transparent',1,'patchSaturation',0.3);

          hold on
        
  
        
    end
    
   
     %Add marker line
     yrange=ylim;
     xlim(PlotRangeSample);
     plot([0,0],yrange,'--k');  
     plot([SampleOffMarker_Mean,SampleOffMarker_Mean],yrange,'--k');
     plot([ChoiceTgtOnMarker_Mean,  ChoiceTgtOnMarker_Mean],yrange,'--r');
     plot([ChoiceTgtOffMarker_SingleMean,  ChoiceTgtOffMarker_SingleMean],yrange,'--r');
     plot([ RwdMarker_Mean,   RwdMarker_Mean],yrange,'--y');
     
     currName=ObjName(i,:);
     currName=currName(~cellfun(@isempty, currName));
    
     
     legend(currName);
     
     
     
     set(findall(gca, 'Type', 'Line'),'LineWidth',3);
set(gca,'FontSize',15,'FontWeight','Bold');
set(gca,'LineWidth',3);

xlabel('Time from Sample Onset(ms)','FontSize',25,'FontWeight','Bold');
ylabel('Spikes/s','FontSize',25,'FontWeight','Bold');
title(sprintf('Set %d',UniqueSet(i)));

     
     
 end
 
 StartIndex=StartIndex+1;
 StartFig=StartFig+1;
 
 
 
 %%%%  Next figure: PSTH of sample period for each object for flexible
 %%%%  value
 for i=1:length(UniqueSet)
 figtitlestr{StartIndex}=strcat('PSTHSamplePeriodInEachSetValueChange',num2str(UniqueSet(i)));
 
 fig{StartIndex}=PrepareFigure(StartFig,'w',[100+(i-1)*20,100, 1200,800],'Name',figtitlestr{StartIndex});
 
    % subplot(ceil(sqrt(length(UniqueSet))),ceil(sqrt(length(UniqueSet))),i);
     
      currSet=UniqueSet(i);
       sel=SampleSetName==currSet;
    
  %   for j=1:
    SampleObjInSet=SampleTargetName(sel);
    
    SampleObjInSet_Unqiue=unique(SampleObjInSet);
    
    
     for j=1:length(SampleObjInSet_Unqiue)
         subplot(ceil(sqrt(length(SampleObjInSet_Unqiue))),ceil(sqrt(length(SampleObjInSet_Unqiue))),j);
         
          sel=PSTH_Sample_Time>PlotRangeSample(1) & PSTH_Sample_Time<PlotRangeSample(2);  
          shadedErrorBar(PSTH_Sample_Time( sel),PSTH_SampleInSetGood_Each_Mean{i,j}( sel),PSTH_SampleInSetGood_Each_Sem{i,j}(sel),'lineprops',{[1,0,0]},'transparent',1,'patchSaturation',0.3);

          hold on
             shadedErrorBar(PSTH_Sample_Time( sel),PSTH_SampleInSetBad_Each_Mean{i,j}( sel),PSTH_SampleInSetBad_Each_Sem{i,j}(sel),'lineprops',{[0,0,1]},'transparent',1,'patchSaturation',0.3);

        
    %Add marker line
     yrange=ylim;
     xlim(PlotRangeSample);
     plot([0,0],yrange,'--k');  
     plot([SampleOffMarker_Mean,SampleOffMarker_Mean],yrange,'--k');
     plot([ChoiceTgtOnMarker_Mean,  ChoiceTgtOnMarker_Mean],yrange,'--r');
     plot([ChoiceTgtOffMarker_SingleMean,  ChoiceTgtOffMarker_SingleMean],yrange,'--r');
     plot([ RwdMarker_Mean,   RwdMarker_Mean],yrange,'--y');
     
     currName=ObjName(i,j);
     
    
     
     legend({'Good','Bad'});
     
     
     
     set(findall(gca, 'Type', 'Line'),'LineWidth',3);

set(gca,'LineWidth',3);

xlabel('Time from Sample Onset(ms)','FontSize',25,'FontWeight','Bold');
ylabel('Spikes/s','FontSize',25,'FontWeight','Bold');
title(strcat('Obj',currName));
set(gca,'FontSize',15,'FontWeight','Bold');

        
    end
    
   
    
StartIndex=StartIndex+1;
 StartFig=StartFig+1;
     
 end
 
 
 
 
 %%%%  Next figure: PSTH of sample period for each object
 figtitlestr{StartIndex}='PSTHSamplePeriodInEachSet_Single';

 fig{StartIndex}=PrepareFigure(StartFig,'w',[80,100, 1200,800],'Name',figtitlestr{StartIndex});
 for i=1:length(UniqueSet)
     subplot(ceil(sqrt(length(UniqueSet))),ceil(sqrt(length(UniqueSet))),i);
     
      currSet=UniqueSet(i);
       sel=SampleSetName==currSet;
    
  %   for j=1:
    SampleObjInSet=SampleTargetName(sel);
    
   % SampleObjInSet_Unqiue=unique(SampleObjInSet);
    
    % for j=1:length(SampleObjInSet_Unqiue)
          sel=PSTH_Sample_Time>PlotRangeSample(1) & PSTH_Sample_Time<PlotRangeSample(2);  
          shadedErrorBar(PSTH_Sample_Time( sel),PSTH_SampleInSetSingle_Good_Mean{i}( sel),PSTH_SampleInSetSingle_Good_Sem{i}(sel),'lineprops',{[1    0    0]},'transparent',1,'patchSaturation',0.3);

    hold on
     shadedErrorBar(PSTH_Sample_Time( sel),PSTH_SampleInSetSingle_Bad_Mean{i}( sel),PSTH_SampleInSetSingle_Bad_Sem{i}( sel),'lineprops',{[0    0    1]},'transparent',1,'patchSaturation',0.3);
     
        
  
        
  %  end
    
   
     %Add marker line
     yrange=ylim;
     xlim(PlotRangeSample);
     plot([0,0],yrange,'--k');  
     plot([SampleOffMarker_Mean,SampleOffMarker_Mean],yrange,'--k');
     plot([ChoiceTgtOnMarker_Mean,  ChoiceTgtOnMarker_Mean],yrange,'--r');
     plot([ChoiceTgtOffMarker_SingleMean,  ChoiceTgtOffMarker_SingleMean],yrange,'--r');
     plot([ RwdMarker_Mean,   RwdMarker_Mean],yrange,'--y');
     
   %  currName=ObjName(i,:);
    % currName=currName(~cellfun(@isempty, currName));
    
     
    % legend(currName);
     
     
     
     set(findall(gca, 'Type', 'Line'),'LineWidth',3);
set(gca,'FontSize',15,'FontWeight','Bold');
set(gca,'LineWidth',3);

xlabel('Time from Sample Onset(ms)','FontSize',25,'FontWeight','Bold');
ylabel('Spikes/s','FontSize',25,'FontWeight','Bold');
title(sprintf('Set %d',UniqueSet(i)));

     
     
 end
 
 StartIndex=StartIndex+1;
 StartFig=StartFig+1;
 
 
 
 
 







end %End of the spike loop



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
