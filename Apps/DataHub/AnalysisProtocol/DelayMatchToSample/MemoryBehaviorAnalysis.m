function  MemoryBehaviorAnalysis(Data, StartTrial, EndTrial,ShowFigureFlag,OutputFlag);
%Analysis protocol for the behavior of delay match to sample tasks
%By xuefei yu 05222021

ECode;


%PlotStartOffset=[2000];%in ms
%PlotStopOffset=[5000];%in ms
%For plots
PlotRangeSample=[-500,2000];
PlotRangeTest=[-1000,3000];

%Time bin pre saccade and time bin offset after sample onset
SampleOnPreBuffer=500;%in ms
SampleOnPostBuffer=5000;%in ms


ExamplePlot=0;%PLots for presentations, turn off for regular use

ResponseLatencyThreshold=3;

%Parameter to define the sliding window to measure firing rates
StepSize=5;%in ms
BinWidth=20;%in ms


%Criterias for sexwlecting saccades
V_Threshold=40;
ContinueBin=5;
ISI_Threshold=20;






%Color Definition for each fractals
%{
Color=[lines(6);prism(6) ];
Color(2,:)=copper(1)+0.5;
Color=[Color;Color*0.5];
%}
Color={'#648FFF','#FE6100','#009E73','#88CCEE','#CC6677','#785EF0'};
Color_RGB=[
           254,97,0;
           100,143,255;
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

%SpikeChannel=Data.SpikeChannel;


if strcmp(FileName,'Robin040321FRACSAC1.rst') & SpikeChannel==22
    Selection=Selection(~ismember(Selection,[14,40]));
end

%Trial outcome;
EventChannel=Data.EventChannel(Selection,:);
EventTimeChannel=Data.EventTimeChannel(Selection,:);
EventBin=Data.EventBin(Selection,:);

%GetEyeData
EyeData=Data.EyeDataRaw;
EyeChannel_X_L=EyeData(Selection,1);
EyeChannel_Y_L=EyeData(Selection,2);
%EyeChannel_X_R=EyeData(Selection,3);
%EyeChannel_Y_R=EyeData(Selection,4);

PuipleThreshold=-50;
PuipleSize=ValidPuiple(EyeData(Selection,3),PuipleThreshold);%Only valid after some date....


EyeBinWidth=Data.EyeBinWidth;



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

SingleTgtTrial=TaskType==1;
ChoiceTrial=TaskType>1;

%SingleTgtTrial=TargetNum==1;
%ChoiceTrial=TargetNum>1;

SampleTargetIndex=TargetIndex(:,1);
SampleTargetName=TargetName(:,1);

ObjectIndex=ReproduceFromEvent(EventChannel,[GOODOBJ  BADOBJ])';
%GoodObjectIndex=(ObjectIndex==GOODOBJ)';

%SampleGood=SampleTargetIndex<=4;
%SampleBad=SampleTargetIndex>4;

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

%TestTargetIndex=TargetIndex(:,2:end);
%Separate into 4 conditions
if size(TargetIndex,2)<3
    %All single targets
    TestTargetIndex=[TargetIndex(:,2:end),NaN*ones(size(TargetIndex,1),1)];
else
    TestTargetIndex=TargetIndex(:,2:end);
    
end

TestGood_Ori=TestTargetIndex(:,2)<=4;

TestGood=TestGood_Ori;
TestGood(FlipIndex)=~TestGood(FlipIndex);
TestBad=~TestGood;

GoodGoodTestAll=TestGood & ChoiceTrial & SampleGood  ;
SampleGoodTestBadAll=TestBad & ChoiceTrial & SampleGood ;
SampleBadTestGoodAll=TestGood & ChoiceTrial & SampleBad ;
BadBadTestAll=TestBad & ChoiceTrial & SampleBad & CorrectTrials;

GoodGoodTest=TestGood & ChoiceTrial & SampleGood  & CorrectTrials;
SampleGoodTestBad=TestBad & ChoiceTrial & SampleGood  & CorrectTrials;
SampleBadTestGood=TestGood & ChoiceTrial & SampleBad  & CorrectTrials;
BadBadTest=TestBad & ChoiceTrial & SampleBad & CorrectTrials;


disp(['GoodGoodTest=',num2str(sum(GoodGoodTest)),';SampleGoodTestBad=',num2str(sum(SampleGoodTestBad)),';SampleBadTestGood=',num2str(sum(SampleBadTestGood)),';SampleGoodTestBad=',num2str(sum(BadBadTest))])

%First check the puiple size difference between good and bad objects

%PuipleSize_SampleGood=PuipleSize(SampleGood,:);


FractalOnsetTime = FindOutTime(EventChannel,EventTimeChannel,TGTCD );
SampleOnsetTime=FractalOnsetTime(:,1);

%PuipleSize_SampleSelection_tmp= SelectEyeInterval(PuipleSize,SampleOnsetBinRange ,SampleOffsetBinRange);
%PuipleSize_SampleSelection=ReorganizeEye(PuipleSize_SampleSelection_tmp);
PuipleSize_SampleSelection=ReorganizeEye({PuipleSize});

%Align puiple data
EyeTime=0:EyeBinWidth:size(PuipleSize_SampleSelection,2)-1;

PuipleSize_Aligned=AlignEyeData(PuipleSize_SampleSelection,EyeTime,SampleOnsetTime,[-SampleOnPreBuffer,SampleOnPostBuffer]);

EyeTime_Sample=[-SampleOnPreBuffer:-SampleOnPreBuffer+size(PuipleSize_Aligned,2)-1];

SponIntervalSample=[-400,-100];
SponPreSample=EyeTime_Sample<SponIntervalSample(2) & EyeTime_Sample>SponIntervalSample(1);

PuipleSize_Aligned_SponEachTrial=nanmean(PuipleSize_Aligned(:,SponPreSample),2);
PuipleSize_AlignedOverspon=PuipleSize_Aligned-PuipleSize_Aligned_SponEachTrial;


PuipleSize_SampleAlign_Good=PuipleSize_AlignedOverspon(SampleGood,:);
PuipleSize_SampleAlign_Bad=PuipleSize_AlignedOverspon(SampleBad,:);

[PuipleSize_SampleAlign_Good_Mean,PuipleSize_SampleAlign_Good_Sem]=StatisticsPSTH(PuipleSize_SampleAlign_Good,3);
[PuipleSize_SampleAlign_Bad_Mean,PuipleSize_SampleAlign_Bad_Sem]=StatisticsPSTH(PuipleSize_SampleAlign_Bad,3);

%%
%Analysis the pre-choice wandering period
%First, the choice confirming duration;
%From fix off to final choice

FixOffTime = FindOutTime(EventChannel,EventTimeChannel,FIXOFF );
GoodChoiceTime= FindOutTime(EventChannel,EventTimeChannel,MATCHGOOD);
BadChoiceTime= FindOutTime(EventChannel,EventTimeChannel,MATCHBAD);

GoodChoiceDur=GoodChoiceTime-FixOffTime;
BadChoiceDur=BadChoiceTime-FixOffTime;


%Time for signle trials
GoodChoiceDurSingle=GoodChoiceDur(SingleTgtTrial & SampleGood);
BadChoiceDurSingle=BadChoiceDur(SingleTgtTrial & SampleBad);

%Time for four conditions in choice trials

ChoiceDur_GG=GoodChoiceDur(GoodGoodTest);
ChoiceDur_GB=GoodChoiceDur(SampleGoodTestBad);

ChoiceDur_BG=BadChoiceDur(SampleBadTestGood);
ChoiceDur_BB=BadChoiceDur(BadBadTest);


ChoiceDur_GG_Mean=nanmean(ChoiceDur_GG);
ChoiceDur_GB_Mean=nanmean(ChoiceDur_GB);
ChoiceDur_BG_Mean=nanmean(ChoiceDur_BG);
ChoiceDur_BB_Mean=nanmean(ChoiceDur_BB);

GoodChoiceDurSingle_Mean=nanmean(GoodChoiceDurSingle);
BadChoiceDurSingle_Mean=nanmean(BadChoiceDurSingle);

ChoiceDur_GG_Sem=nanstd(ChoiceDur_GG)/sqrt(sum(~isnan(ChoiceDur_GG)));
ChoiceDur_GB_Sem=nanstd(ChoiceDur_GB)/sqrt(sum(~isnan(ChoiceDur_GB)));
ChoiceDur_BG_Sem=nanstd(ChoiceDur_BG)/sqrt(sum(~isnan(ChoiceDur_BG)));
ChoiceDur_BB_Sem=nanstd(ChoiceDur_BB)/sqrt(sum(~isnan(ChoiceDur_BB)));

GoodChoiceDurSingle_Sem=nanstd(GoodChoiceDurSingle)/sqrt(sum(~isnan(GoodChoiceDurSingle)));
BadChoiceDurSingle_Sem=nanstd(BadChoiceDurSingle)/sqrt(sum(~isnan(BadChoiceDurSingle)));




ChoiceDur_MeanArray=[ChoiceDur_GG_Mean,ChoiceDur_GB_Mean,ChoiceDur_BG_Mean,ChoiceDur_BB_Mean,GoodChoiceDurSingle_Mean,BadChoiceDurSingle_Mean];
ChoiceDur_SemArray=[ChoiceDur_GG_Sem,ChoiceDur_GB_Sem,ChoiceDur_BG_Sem,ChoiceDur_BB_Sem,GoodChoiceDurSingle_Sem,BadChoiceDurSingle_Sem];

ChoiceDurArray_Str={'GG','GB','BG','BB','G-Only','B-Only'};

%First Saccade Into
FirstChoiceIndex=ReproduceFromEvent(EventChannel,[INTOMATCH  INTONONMATCH])';
NumOfWandering=cellfun(@numel,FirstChoiceIndex);

%Separate into different condition

NumOfWandering_GG=NumOfWandering(GoodGoodTest);
NumOfWandering_GB=NumOfWandering(SampleGoodTestBad);

NumOfWandering_BG=NumOfWandering(SampleBadTestGood);
NumOfWandering_BB=NumOfWandering(BadBadTest);


NumOfWandering_GG_Mean=nanmean(NumOfWandering_GG);
NumOfWandering_GB_Mean=nanmean(NumOfWandering_GB);

NumOfWandering_BG_Mean=nanmean(NumOfWandering_BG);
NumOfWandering_BB_Mean=nanmean(NumOfWandering_BB);



NumOfWandering_GG_Sem=nanstd(NumOfWandering_GG)/sqrt(length(NumOfWandering_GG));
NumOfWandering_GB_Sem=nanstd(NumOfWandering_GB)/sqrt(length(NumOfWandering_GB));

NumOfWandering_BG_Sem=nanstd(NumOfWandering_BG)/sqrt(length(NumOfWandering_BG));
NumOfWandering_BB_Sem=nanstd(NumOfWandering_BB)/sqrt(length(NumOfWandering_BG));

NumOfWandering_MeanArray=[NumOfWandering_GG_Mean,NumOfWandering_GB_Mean,NumOfWandering_BG_Mean,NumOfWandering_BB_Mean];
NumOfWandering_SemArray=[NumOfWandering_GG_Sem,NumOfWandering_GB_Sem,NumOfWandering_BG_Sem,NumOfWandering_BB_Sem];

NumOfWandering_Str={'GG','GB','BG','BB'};

%First saccade towards 
FirstChoice=cellfun(@(x) x(1),FirstChoiceIndex );

FirstChoiceObjIndex=NaN*ones(size(FirstChoice,1),1);
FirstChoiceObjIndex(FirstChoice==INTOMATCH)=TestTargetIndex(FirstChoice==INTOMATCH,1);
FirstChoiceObjIndex(FirstChoice==INTONONMATCH)=TestTargetIndex(FirstChoice==INTONONMATCH,2);

FirstChoiceObjIndexGood=FirstChoiceObjIndex<=4;
FirstChoiceObjIndexGood(FlipIndex)=~FirstChoiceObjIndexGood(FlipIndex);


FirstChoiceObjIndexBad=FirstChoiceObjIndex>4;
FirstChoiceObjIndexBad(FlipIndex)=~FirstChoiceObjIndexBad(FlipIndex);



FirstChoiceObjIndexGood_GG=FirstChoiceObjIndexGood(GoodGoodTest);
FirstChoiceObjIndexGood_GB=FirstChoiceObjIndexGood(SampleGoodTestBad);

FirstChoiceObjIndexGood_BG=FirstChoiceObjIndexGood(SampleBadTestGood);
FirstChoiceObjIndexGood_BB=FirstChoiceObjIndexGood(BadBadTest);



FirstChoiceObjIndexGood_GG_Prop=sum(FirstChoiceObjIndexGood_GG)/length(FirstChoiceObjIndexGood_GG);
FirstChoiceObjIndexGood_GB_Prop=sum(FirstChoiceObjIndexGood_GB)/length(FirstChoiceObjIndexGood_GB);
FirstChoiceObjIndexGood_BG_Prop=sum(FirstChoiceObjIndexGood_BG)/length(FirstChoiceObjIndexGood_BG);

FirstChoiceObjIndexGood_BB_Prop=sum(FirstChoiceObjIndexGood_BB)/length(FirstChoiceObjIndexGood_BB);

FirstChoiceObjIndexGood_PropArray=[FirstChoiceObjIndexGood_GG_Prop,FirstChoiceObjIndexGood_GB_Prop,FirstChoiceObjIndexGood_BG_Prop,FirstChoiceObjIndexGood_BB_Prop];


%Finally saccade latency

%Select the period between target onset and target off set
FractalOnsetBin = FindOutTime(EventChannel,EventBin,TGTCD )-100;%-100: offset
FractalOffsetBin = FindOutTime(EventChannel,EventBin,TGTOFF);
FractalOnsetTime = FindOutTime(EventChannel,EventTimeChannel,TGTCD );
FractalOffsetTime = FindOutTime(EventChannel,EventTimeChannel,TGTOFF );

EyeChannel_X_Selection_All= SelectEyeInterval(EyeChannel_X_L,FractalOnsetBin ,FractalOffsetBin )';
EyeChannel_Y_Selection_All= SelectEyeInterval(EyeChannel_Y_L,FractalOnsetBin ,FractalOffsetBin )';

EyeChannel_X_Selection=ReorganizeEye(EyeChannel_X_Selection_All{2});%Only need the test period
EyeChannel_Y_Selection=ReorganizeEye(EyeChannel_Y_Selection_All{2});%Only need the test period



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

TimeThreshold=-100;

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


RT_First_Target=RT_First_Target-100;

RT_First_Target_SingleGood=RT_First_Target(SingleTgtTrial & SampleGood);
RT_First_Target_SingleBad=RT_First_Target(SingleTgtTrial & SampleBad);

RT_First_Target_GG=RT_First_Target(GoodGoodTest);
RT_First_Target_GB=RT_First_Target(SampleGoodTestBad);
RT_First_Target_BG=RT_First_Target(SampleBadTestGood);
RT_First_Target_BB=RT_First_Target(BadBadTest);


RT_First_Target_SingleGood_Mean=nanmean(RT_First_Target_SingleGood);
RT_First_Target_SingleBad_Mean=nanmean(RT_First_Target_SingleBad);

RT_First_Target_GG_Mean=nanmean(RT_First_Target_GG);
RT_First_Target_GB_Mean=nanmean(RT_First_Target_GB);
RT_First_Target_BG_Mean=nanmean(RT_First_Target_BG);
RT_First_Target_BB_Mean=nanmean(RT_First_Target_BB);


RT_First_Target_SingleGood_Sem=nanstd(RT_First_Target_SingleGood)./sqrt(length(RT_First_Target_SingleGood));
RT_First_Target_SingleBad_Sem=nanstd(RT_First_Target_SingleBad)./sqrt(length(RT_First_Target_SingleBad));

RT_First_Target_GG_Sem=nanstd(RT_First_Target_GG)./sqrt(length(RT_First_Target_GG));
RT_First_Target_GB_Sem=nanstd(RT_First_Target_GB)./sqrt(length(RT_First_Target_GB));
RT_First_Target_BG_Sem=nanstd(RT_First_Target_BG)./sqrt(length(RT_First_Target_BG));
RT_First_Target_BB_Sem=nanstd(RT_First_Target_BB)./sqrt(length(RT_First_Target_BB));



RT_Array_Mean=[RT_First_Target_GG_Mean,RT_First_Target_GB_Mean,RT_First_Target_BG_Mean,RT_First_Target_BB_Mean,RT_First_Target_SingleGood_Mean,RT_First_Target_SingleBad_Mean];
RT_Array_Sem=[RT_First_Target_GG_Sem,RT_First_Target_GB_Sem,RT_First_Target_BG_Sem,RT_First_Target_BB_Sem,RT_First_Target_SingleGood_Sem,RT_First_Target_SingleBad_Sem];

RT_Str={'GG','GB','BG','BB','Single_G','Single_B'};


%Finally correct rate
CorrectRate_GG=sum(CorrectTrials(GoodGoodTestAll))/sum(GoodGoodTestAll);
CorrectRate_GB=sum(CorrectTrials(SampleGoodTestBadAll))/sum(SampleGoodTestBadAll);
CorrectRate_BG=sum(CorrectTrials(SampleBadTestGoodAll))/sum(SampleBadTestGoodAll);
CorrectRate_BB=sum(CorrectTrials(BadBadTestAll))/sum(BadBadTestAll);

CorrectRate_Array=[CorrectRate_GG,CorrectRate_GB,CorrectRate_BG,CorrectRate_BB];

CorrectRate_Str={'GG','GB','BG','BB'};

%Separate for each set

SampleTargetSet=floor(SampleTargetName/10);
UniqueSet=unique(SampleTargetSet);

if length(UniqueSet)>1
    %Loop over the set if more than one set
    for i=1:length(UniqueSet)
        currSet=UniqueSet(i);
        setsel=SampleTargetSet==currSet;
        
        %Choice Duration
        GoodChoiceDurSingle=GoodChoiceDur(SingleTgtTrial & SampleGood & setsel);
        BadChoiceDurSingle=BadChoiceDur(SingleTgtTrial & SampleBad &setsel);

        %Time for four conditions in choice trials

        ChoiceDur_GG=GoodChoiceDur(GoodGoodTest& setsel);
        ChoiceDur_GB=GoodChoiceDur(SampleGoodTestBad& setsel);
        ChoiceDur_BG=BadChoiceDur(SampleBadTestGood& setsel);
        ChoiceDur_BB=BadChoiceDur(BadBadTest& setsel);


        ChoiceDur_GG_Mean=nanmean(ChoiceDur_GG);
        ChoiceDur_GB_Mean=nanmean(ChoiceDur_GB);
        ChoiceDur_BG_Mean=nanmean(ChoiceDur_BG);
        ChoiceDur_BB_Mean=nanmean(ChoiceDur_BB);

        GoodChoiceDurSingle_Mean=nanmean(GoodChoiceDurSingle);
        BadChoiceDurSingle_Mean=nanmean(BadChoiceDurSingle);

        ChoiceDur_GG_Sem=nanstd(ChoiceDur_GG)/sqrt(sum(~isnan(ChoiceDur_GG)));
        ChoiceDur_GB_Sem=nanstd(ChoiceDur_GB)/sqrt(sum(~isnan(ChoiceDur_GB)));
        ChoiceDur_BG_Sem=nanstd(ChoiceDur_BG)/sqrt(sum(~isnan(ChoiceDur_BG)));
        ChoiceDur_BB_Sem=nanstd(ChoiceDur_BB)/sqrt(sum(~isnan(ChoiceDur_BB)));

         GoodChoiceDurSingle_Sem=nanstd(GoodChoiceDurSingle)/sqrt(sum(~isnan(GoodChoiceDurSingle)));
         BadChoiceDurSingle_Sem=nanstd(BadChoiceDurSingle)/sqrt(sum(~isnan(BadChoiceDurSingle)));




        ChoiceDur_MeanArray_set(i,:)=[ChoiceDur_GG_Mean,ChoiceDur_GB_Mean,ChoiceDur_BG_Mean,ChoiceDur_BB_Mean,GoodChoiceDurSingle_Mean,BadChoiceDurSingle_Mean];
        ChoiceDur_SemArray_Set(i,:)=[ChoiceDur_GG_Sem,ChoiceDur_GB_Sem,ChoiceDur_BG_Sem,ChoiceDur_BB_Sem,GoodChoiceDurSingle_Sem,BadChoiceDurSingle_Sem];

       %Puiple Change
       PuipleSize_SampleAlign_Good=PuipleSize_AlignedOverspon(SampleGood& setsel,:);
       PuipleSize_SampleAlign_Bad=PuipleSize_AlignedOverspon(SampleBad& setsel,:);

       [PuipleSize_SampleAlign_Good_Mean_set(i,:),PuipleSize_SampleAlign_Good_Sem_set(i,:)]=StatisticsPSTH(PuipleSize_SampleAlign_Good,3);
       [PuipleSize_SampleAlign_Bad_Mean_set(i,:),PuipleSize_SampleAlign_Bad_Sem_set(i,:)]=StatisticsPSTH(PuipleSize_SampleAlign_Bad,3);
       
       %Reaction time

    RT_First_Target_SingleGood=RT_First_Target(SingleTgtTrial & SampleGood& setsel);
    RT_First_Target_SingleBad=RT_First_Target(SingleTgtTrial & SampleBad& setsel);

    RT_First_Target_GG=RT_First_Target(GoodGoodTest& setsel);
    RT_First_Target_GB=RT_First_Target(SampleGoodTestBad& setsel);
    RT_First_Target_BG=RT_First_Target(SampleBadTestGood& setsel);
    RT_First_Target_BB=RT_First_Target(BadBadTest& setsel);


    RT_First_Target_SingleGood_Mean=nanmean(RT_First_Target_SingleGood);
    RT_First_Target_SingleBad_Mean=nanmean(RT_First_Target_SingleBad);

    RT_First_Target_GG_Mean=nanmean(RT_First_Target_GG);
    RT_First_Target_GB_Mean=nanmean(RT_First_Target_GB);
    RT_First_Target_BG_Mean=nanmean(RT_First_Target_BG);
    RT_First_Target_BB_Mean=nanmean(RT_First_Target_BB);


    RT_First_Target_SingleGood_Sem=nanstd(RT_First_Target_SingleGood)./sqrt(length(RT_First_Target_SingleGood));
    RT_First_Target_SingleBad_Sem=nanstd(RT_First_Target_SingleBad)./sqrt(length(RT_First_Target_SingleBad));

    RT_First_Target_GG_Sem=nanstd(RT_First_Target_GG)./sqrt(length(RT_First_Target_GG));
    RT_First_Target_GB_Sem=nanstd(RT_First_Target_GB)./sqrt(length(RT_First_Target_GB));
    RT_First_Target_BG_Sem=nanstd(RT_First_Target_BG)./sqrt(length(RT_First_Target_BG));
    RT_First_Target_BB_Sem=nanstd(RT_First_Target_BB)./sqrt(length(RT_First_Target_BB));



RT_Array_Mean_set(i,:)=[RT_First_Target_GG_Mean,RT_First_Target_GB_Mean,RT_First_Target_BG_Mean,RT_First_Target_BB_Mean,RT_First_Target_SingleGood_Mean,RT_First_Target_SingleBad_Mean];
RT_Array_Sem_set(i,:)=[RT_First_Target_GG_Sem,RT_First_Target_GB_Sem,RT_First_Target_BG_Sem,RT_First_Target_BB_Sem,RT_First_Target_SingleGood_Sem,RT_First_Target_SingleBad_Sem];

       
   %Correct Rate:
   
    CorrectRate_GG=sum(CorrectTrials(GoodGoodTestAll& setsel))/sum(GoodGoodTestAll& setsel);
CorrectRate_GB=sum(CorrectTrials(SampleGoodTestBadAll& setsel))/sum(SampleGoodTestBadAll& setsel);
CorrectRate_BG=sum(CorrectTrials(SampleBadTestGoodAll& setsel))/sum(SampleBadTestGoodAll& setsel);
CorrectRate_BB=sum(CorrectTrials(BadBadTestAll& setsel))/sum(BadBadTestAll& setsel);

CorrectRate_Array_set(i,:)=[CorrectRate_GG,CorrectRate_GB,CorrectRate_BG,CorrectRate_BB];

%CorrectRate_Str={'GG','GB','BG','BB'};
    
        
        
        
    end
    
    
    
end



StartIndex=1;
 StartFig=100;

if ShowFigureFlag
     figtitlestr{StartIndex}='BehaviorMeasurement';

     fig{StartIndex}=PrepareFigure(StartFig,'w',[50,100, 1200,800],'Name',figtitlestr{StartIndex});
    
     sel=EyeTime_Sample>PlotRangeSample(1) & EyeTime_Sample<PlotRangeSample(2);
    subplot(3,2,1);%Sample On two target
     shadedErrorBar(EyeTime_Sample( sel),PuipleSize_SampleAlign_Good_Mean( sel),PuipleSize_SampleAlign_Good_Sem( sel),'lineprops',{[1    0    0]},'transparent',1,'patchSaturation',0.3);

    hold on
     shadedErrorBar(EyeTime_Sample( sel),PuipleSize_SampleAlign_Bad_Mean( sel),PuipleSize_SampleAlign_Bad_Sem( sel),'lineprops',{[0    0    1]},'transparent',1,'patchSaturation',0.3);
     %{
     %Add marker line
     yrange=ylim;
     xlim(PlotRangeSample);
     plot([0,0],yrange,'--k');  
     plot([SampleOffMarker_Mean,SampleOffMarker_Mean],yrange,'--k');
     plot([ChoiceTgtOnMarker_Mean,  ChoiceTgtOnMarker_Mean],yrange,'--r');
     plot([ChoiceTgtOffMarker_Mean,  ChoiceTgtOffMarker_Mean],yrange,'--r');
     plot([ RwdMarker_Mean,   RwdMarker_Mean],yrange,'--y');
     %}
    
     set(findall(gca, 'Type', 'Line'),'LineWidth',3);
     set(gca,'FontSize',15,'FontWeight','Bold');
set(gca,'LineWidth',3);

xlabel('Time from Sample Onset(ms)','FontSize',25,'FontWeight','Bold');
ylabel('PuipleSize(a.u.)','FontSize',25,'FontWeight','Bold');
title('Sample Period');
    
    subplot(3,2,2);%Prechoice wandering duration
    g=bar([ChoiceDur_MeanArray],'LineWidth',3);
%g=bar([FR_Scene_mean_Sorted,Mean_Spon(1)],'LineWidth',3);
g.FaceColor = 'flat';
g.EdgeColor = 'flat';

hold on
%errorbar([FR_Scene_mean_Sorted,Mean_Spon(1)],[FR_Scene_sem(Index_Sort),Mean_Spon(2)],'.k','LineWidth',3);
errorbar([ChoiceDur_MeanArray],[ChoiceDur_SemArray],'.k','LineWidth',3);

ColorChosen= Color_RGB(1:6,:);

g.CData = ColorChosen;

xticklabels(ChoiceDurArray_Str);
ylabel('Searching Duration(ms)');
set(gca,'FontSize',15,'FontWeight','Bold');
set(gca,'LineWidth',3);
box off;


subplot(3,2,3);%Number of choice wandering 
 g=bar([NumOfWandering_MeanArray],'LineWidth',3);
%g=bar([FR_Scene_mean_Sorted,Mean_Spon(1)],'LineWidth',3);
g.FaceColor = 'flat';
g.EdgeColor = 'flat';

hold on
%errorbar([FR_Scene_mean_Sorted,Mean_Spon(1)],[FR_Scene_sem(Index_Sort),Mean_Spon(2)],'.k','LineWidth',3);
errorbar([NumOfWandering_MeanArray],[NumOfWandering_SemArray],'.k','LineWidth',3);

ColorChosen= Color_RGB(1:4,:);

g.CData = ColorChosen;

xticklabels(NumOfWandering_Str);
ylabel('NumberOfWandering');
set(gca,'FontSize',15,'FontWeight','Bold');
set(gca,'LineWidth',3);
box off;


subplot(3,2,4);%Proportion of first choice direction
 g=bar([FirstChoiceObjIndexGood_PropArray],'LineWidth',3);
%g=bar([FR_Scene_mean_Sorted,Mean_Spon(1)],'LineWidth',3);
g.FaceColor = 'flat';
g.EdgeColor = 'flat';


ColorChosen= Color_RGB(1:4,:);

g.CData = ColorChosen;

xticklabels(NumOfWandering_Str);
ylabel('ProportionOfFirstChoiceTowardsGood');
set(gca,'FontSize',15,'FontWeight','Bold');
set(gca,'LineWidth',3);
box off;


subplot(3,2,5);%Saccade Latency
 g=bar([RT_Array_Mean],'LineWidth',3);
%g=bar([FR_Scene_mean_Sorted,Mean_Spon(1)],'LineWidth',3);
g.FaceColor = 'flat';
g.EdgeColor = 'flat';

hold on
errorbar([RT_Array_Mean],[RT_Array_Sem],'.k','LineWidth',3);

ColorChosen= Color_RGB(1:length(RT_Array_Mean),:);

g.CData = ColorChosen;

xticklabels(RT_Str);
ylabel('Reaction Time(ms)');
set(gca,'FontSize',15,'FontWeight','Bold');
set(gca,'LineWidth',3);
box off;


subplot(3,2,6);%Correct Rate
 g=bar([CorrectRate_Array],'LineWidth',3);
%g=bar([FR_Scene_mean_Sorted,Mean_Spon(1)],'LineWidth',3);
g.FaceColor = 'flat';
g.EdgeColor = 'flat';

hold on

ColorChosen= Color_RGB(1:length(CorrectRate_Array),:);

g.CData = ColorChosen;

xticklabels(CorrectRate_Str);
ylabel('Correct Rate');
set(gca,'FontSize',15,'FontWeight','Bold');
set(gca,'LineWidth',3);
box off;




%NumOfWandering_MeanArray
%CorrectRate_Array=
 
%Add title

title_msg=sprintf(' %s, TrialNum: %d ',FileName,TrialNum);
a=annotation('textbox',[0.25,0.68,0.3,0.3],'String',title_msg,'FitBoxToText','on');
a.FontSize=25;
a.FontWeight='Bold';
a.LineStyle='none';  

StartIndex=StartIndex+1;
StartFig=StartFig+1;
 

%%%Next figure for each set
if length(UniqueSet)>1
figtitlestr{StartIndex}='BehaviorQuantityForEachSet';
fig{StartIndex}=PrepareFigure(StartFig,'w',[70,100, 1200,800],'Name',figtitlestr{StartIndex});

for i=1:length(UniqueSet)
    %second column
    %Prechoice wandering period
axes('Position',[0.3,0.8-(i-1)*(0.6/max(length(UniqueSet)-1,1)),0.23,0.15]);
 g=bar([ChoiceDur_MeanArray_set(i,:)],'LineWidth',3);
%g=bar([FR_Scene_mean_Sorted,Mean_Spon(1)],'LineWidth',3);
g.FaceColor = 'flat';
g.EdgeColor = 'flat';

hold on
%errorbar([FR_Scene_mean_Sorted,Mean_Spon(1)],[FR_Scene_sem(Index_Sort),Mean_Spon(2)],'.k','LineWidth',3);
errorbar(ChoiceDur_MeanArray_set(i,:),ChoiceDur_SemArray_Set(i,:),'.k','LineWidth',3);

ColorChosen= Color_RGB(1:6,:);

g.CData = ColorChosen;

xticklabels(ChoiceDurArray_Str);
ylabel('Searching Duration(ms)');
title(strcat('Set:',num2str(UniqueSet(i))));
set(gca,'FontSize',15,'FontWeight','Bold');
set(gca,'LineWidth',3);
box off;

%First column
%Puiple Size Change
    
axes('Position',[0.05,0.8-(i-1)*(0.6/max(length(UniqueSet)-1,1)),0.2,0.15]);
%PuipleSize
 sel=EyeTime_Sample>PlotRangeSample(1) & EyeTime_Sample<PlotRangeSample(2);
   
     shadedErrorBar(EyeTime_Sample( sel),PuipleSize_SampleAlign_Good_Mean_set(i, sel),PuipleSize_SampleAlign_Good_Sem_set(i, sel),'lineprops',{[1    0    0]},'transparent',1,'patchSaturation',0.3);

    hold on
     shadedErrorBar(EyeTime_Sample( sel),PuipleSize_SampleAlign_Bad_Mean_set(i, sel),PuipleSize_SampleAlign_Bad_Sem_set(i, sel),'lineprops',{[0    0    1]},'transparent',1,'patchSaturation',0.3);
     %{
     %Add marker line
     yrange=ylim;
     xlim(PlotRangeSample);
     plot([0,0],yrange,'--k');  
     plot([SampleOffMarker_Mean,SampleOffMarker_Mean],yrange,'--k');
     plot([ChoiceTgtOnMarker_Mean,  ChoiceTgtOnMarker_Mean],yrange,'--r');
     plot([ChoiceTgtOffMarker_Mean,  ChoiceTgtOffMarker_Mean],yrange,'--r');
     plot([ RwdMarker_Mean,   RwdMarker_Mean],yrange,'--y');
     %}
    
     set(findall(gca, 'Type', 'Line'),'LineWidth',3);
     set(gca,'FontSize',15,'FontWeight','Bold');
set(gca,'LineWidth',3);
if i==length(UniqueSet)
xlabel('Time from Sample Onset(ms)','FontSize',15,'FontWeight','Bold');
end
ylabel('PuipleSize(a.u.)','FontSize',15,'FontWeight','Bold');
title(strcat('Set:',num2str(UniqueSet(i))));

%Third column: Reaction Time(ms)

axes('Position',[0.55,0.8-(i-1)*(0.6/max(length(UniqueSet)-1,1)),0.23,0.15]);

 g=bar([RT_Array_Mean_set(i,:)],'LineWidth',3);
%g=bar([FR_Scene_mean_Sorted,Mean_Spon(1)],'LineWidth',3);
g.FaceColor = 'flat';
g.EdgeColor = 'flat';

hold on
errorbar([RT_Array_Mean_set(i,:)],[RT_Array_Sem_set(i,:)],'.k','LineWidth',3);

ColorChosen= Color_RGB(1:length(RT_Array_Mean_set(i,:)),:);

g.CData = ColorChosen;

xticklabels(RT_Str);
ylabel('Reaction Time(ms)');
title(strcat('Set:',num2str(UniqueSet(i))));
set(gca,'FontSize',15,'FontWeight','Bold');
set(gca,'LineWidth',3);
box off;

%fourth column: Correct Rate
axes('Position',[0.8,0.8-(i-1)*(0.6/max(length(UniqueSet)-1,1)),0.21,0.15]);

 g=bar([CorrectRate_Array_set(i,:)],'LineWidth',3);
%g=bar([FR_Scene_mean_Sorted,Mean_Spon(1)],'LineWidth',3);
g.FaceColor = 'flat';
g.EdgeColor = 'flat';

hold on

ColorChosen= Color_RGB(1:length(CorrectRate_Array_set(i,:)),:);

g.CData = ColorChosen;

xticklabels(CorrectRate_Str);
ylabel('Correct Rate');
title(strcat('Set:',num2str(UniqueSet(i))));
set(gca,'FontSize',15,'FontWeight','Bold');
set(gca,'LineWidth',3);
box off;




end


end %End of if length(UniqueSet)>1
    
StartIndex=StartIndex+1;
StartFig=StartFig+1;


end



end

function organized=ReorganizeEye(data);
organized=[];

for i=1:length(data)
    data_curr=data{i};
    try
    maxnumel=max(cellfun(@numel,data_curr));
    
    organized{i}=cell2mat(cellfun(@(x) [x,NaN*ones(1,maxnumel-length(x))],data_curr,'uniform',0));
    catch
    organized{i}=data_curr;
    end
    
end


 maxnumel=max(cellfun(@(x) size(x,2),organized));
  organized=cell2mat(cellfun(@(x) [x,NaN*ones(size(x,1),maxnumel-size(x,2))],organized,'uniform',0)');

end

function output=ValidPuiple(data,the);

output=cellfun(@(x) addNaN(x,the),data,'uniform',0);


end

function output=addNaN(x,the);
output=x;
output(output<the)=NaN;

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
    
    ObjectIndex{i}=event(i,ismember(event(i,:),code));
    
  
    
end
numeach=cellfun(@numel,ObjectIndex);
nummax=max(numeach);
if nummax==1
    ObjectIndex=cell2mat(ObjectIndex);

    
end


end
