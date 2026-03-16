function SequenceMemoryBehavior(Data, StartTrial, EndTrial,ShowFigureFlag,OutputFlag);
%Analysis protocol for sequence memory task; by YXF 07212021

ECode;
FileName=Data.FileName;
TaskType=Data.TaskType;
TaskCode=Data.TaskCode;
%TrialNum=EndTrial-StartTrial+1;
Selection=StartTrial:EndTrial;%Select the trials according to the pannel

%GetEyeData
EyeData=Data.EyeDataRaw;
EyeChannel_X_L=EyeData(Selection,1);
EyeChannel_Y_L=EyeData(Selection,2);
%EyeChannel_X_R=EyeData(Selection,3);
%EyeChannel_Y_R=EyeData(Selection,4);

EyeBinWidth=Data.EyeBinWidth;
EyeCalInfo=Data.EyeCalInfo;%Each row: X_Gain,X_Bias; Y_Gain; Y_bias; X=(raw+X_bias)*X_Gain
EyeCalInfo_L=EyeCalInfo(1,:);

TrialNum=length(Selection);


%Load event data
EventChannel=Data.EventChannel(Selection,:);
EventTimeChannel=Data.EventTimeChannel(Selection,:);
EventBin=Data.EventBin(Selection,:);



%Get the choice task information
ParaLib=Data.ParaLib;
ChoiceTaskTypeRaw=ParaLib('TaskType');

ChoiceTaskType=ChoiceTaskTypeRaw(Selection,2);

BehaviorOutcomeRaw=ParaLib('BehaviorOutCome');
BehaviorOutcome=BehaviorOutcomeRaw(Selection,2);

%Select out choice task
ChoiceTrialSel=ChoiceTaskType==2;


BehaviorOutcome_ChoiceTrials=BehaviorOutcome(ChoiceTrialSel);

CorrectChoice=BehaviorOutcome_ChoiceTrials==CORRECTCD;



TotalChoiceTrials=sum(ChoiceTrialSel);
TotalCorrectChoice=sum(CorrectChoice);
ProportionOfCorrectChoice=sum(CorrectChoice)/sum(ChoiceTrialSel);



%Find out the sequence length

ObjectIndex=ReproduceFromEvent(EventChannel,[MATCHOBJ  NONMATCHOBJ]);

SeqLengthPerformed=cellfun(@numel,ObjectIndex);


SeqLengthPerformed(~CorrectChoice)=SeqLengthPerformed(~CorrectChoice)-1;

SeqLengthChoice=SeqLengthPerformed(ChoiceTrialSel);

SqeAll=1:5;
SeqCounts_Seq=histcounts(SeqLengthChoice,SqeAll);
CorrectCounts_Seq(1)=sum(SeqCounts_Seq(2:end));
CorrectCounts_Seq(2)=sum(SeqCounts_Seq(3:end));
CorrectCounts_Seq(3)=sum(SeqCounts_Seq(4:end));
CorrectRate_Seq=CorrectCounts_Seq/length(SeqLengthChoice);


p_Correct(1)=myBinomTest(CorrectCounts_Seq(1),length(SeqLengthChoice),0.5);
p_Correct(2)=myBinomTest(CorrectCounts_Seq(2),length(SeqLengthChoice),0.25);
p_Correct(3)=myBinomTest(CorrectCounts_Seq(3),length(SeqLengthChoice),0.125);




CorrectRate_Condition(1)=CorrectCounts_Seq(1)/length(SeqLengthChoice);
CorrectRate_Condition(2)=CorrectCounts_Seq(2)/CorrectCounts_Seq(1);
CorrectRate_Condition(3)=CorrectCounts_Seq(3)/CorrectCounts_Seq(2);


p_ConditionedCorrect(1)=myBinomTest(CorrectCounts_Seq(1),length(SeqLengthChoice),0.5);
p_ConditionedCorrect(2)=myBinomTest(CorrectCounts_Seq(2),CorrectCounts_Seq(1),0.5);
p_ConditionedCorrect(3)=myBinomTest(CorrectCounts_Seq(3),CorrectCounts_Seq(2),0.5);

%Event based correct rate: Correct Rate for each fractal in each location
FractalIndex=Data.TargetIndex(Selection,2:end,1);
FractalIndexChoice=FractalIndex(ChoiceTrialSel,:);

FractalIndexChoice_Target=FractalIndexChoice(:,1:2:8);
FractalIndexChoice_Distractor=FractalIndexChoice(:,2:2:8);
FractalIndexChoice_Distance=FractalIndexChoice_Distractor-FractalIndexChoice_Target;
FractalIndexChoice_Distance(FractalIndexChoice_Distance<0)=4+FractalIndexChoice_Distance(FractalIndexChoice_Distance<0);%Symmetry

ObjectIndexMat=cell2mat(cellfun(@(x) [x,NaN*ones(1,4-numel(x))],ObjectIndex,'uniform',0)');
ObjectIndexMat_Choice=ObjectIndexMat(ChoiceTrialSel,:);

CorrectMat_Choice=ObjectIndexMat_Choice;
CorrectMat_Choice(ObjectIndexMat_Choice==MATCHOBJ)=1;
CorrectMat_Choice(ObjectIndexMat_Choice== NONMATCHOBJ)=0;

FractalTempoLocMat_Choice=repmat([1:4],size(CorrectMat_Choice,1),1);

%CorrectRateForEachObject in each location
for i=1:4
    for j=2:4
        
        CorrectCount_Obj_TmpLoc(i,j-1)=sum(sum(FractalIndexChoice_Target==i & FractalTempoLocMat_Choice==j & CorrectMat_Choice==1));
        ObjConditionCount_Obj_TmpLoc(i,j-1)=sum(sum(FractalIndexChoice_Target==i & FractalTempoLocMat_Choice==j));
        p_Obj_TmpLoc(i,j-1)=myBinomTest(CorrectCount_Obj_TmpLoc(i,j-1),ObjConditionCount_Obj_TmpLoc(i,j-1),0.5);
    end
end

CorrectRate_Obj_TmpLoc=CorrectCount_Obj_TmpLoc./ObjConditionCount_Obj_TmpLoc;

CorrectRate_Obj=sum(CorrectCount_Obj_TmpLoc,2)./sum(ObjConditionCount_Obj_TmpLoc,2);
p_CorrectRate_Obj=myBinomTest(sum(CorrectCount_Obj_TmpLoc,2),sum(ObjConditionCount_Obj_TmpLoc,2),0.5);

CorrectRate_TmpLoc=sum(CorrectCount_Obj_TmpLoc,1)./sum(ObjConditionCount_Obj_TmpLoc,1);
p_CorrectRate_TmpLoc=myBinomTest(sum(CorrectCount_Obj_TmpLoc,1),sum(ObjConditionCount_Obj_TmpLoc,1),0.5);

%Correct Rate for distractor distance

for i=1:3
    
        
    seldist=FractalIndexChoice_Distance(:,2:end)==i;
    CorrectCount_Obj_Dist(i)=sum(sum(seldist & CorrectMat_Choice(:,2:end)==1));
    ConditionCount_Obj_Dist(i)=sum(sum(seldist));
    
    p_CorrectRate_Obj_Dist(i)=myBinomTest(CorrectCount_Obj_Dist(i),ConditionCount_Obj_Dist(i),0.5);
    
  
    
    
    
end
distance_unique=[1:3];
CorrectRate_Obj_Dist=CorrectCount_Obj_Dist./ConditionCount_Obj_Dist;


%For each fractal
for j=1:4
    for i=1:3
    
    
        
    seldist=FractalIndexChoice_Distance(:,2:end)==i & FractalIndexChoice_Target(:,2:end)==j;
    CorrectCount_Obj_Dist_Obj(j,i)=sum(sum(seldist & CorrectMat_Choice(:,2:end)==1));
    ConditionCount_Obj_Dist_Obj(j,i)=sum(sum(seldist));
    
    p_CorrectRate_Obj_Dist_Obj(j,i)=myBinomTest(CorrectCount_Obj_Dist_Obj(j,i),ConditionCount_Obj_Dist_Obj(j,i),0.5);
    
  
    
    
    
    end
end
CorrectRate_Obj_Dist_Obj=CorrectCount_Obj_Dist_Obj./ConditionCount_Obj_Dist_Obj;

%SeqLengthCorrect()

%Error analysis

SequenceRelatedRaw=ParaLib('Entry,SequenceLength,SequenceType');
SequenceRelated=SequenceRelatedRaw(Selection,2:end);

Entry=SequenceRelated(:,1)+1;
EntryWhole=1:4+1;



if ismember('DistractroArray',keys(ParaLib))
     DistractorTypeRaw=ParaLib('DistractroArray');
    DistractorType=DistractorTypeRaw(Selection,2)+1;

    DistractortTypeWhole=1:12;
    
    
    
else
    
    DistractorTypeRaw=ParaLib('DistractorType');
    DistractorType=DistractorTypeRaw(Selection,2);

    DistractortTypeWhole=0:12;


    
    
    
    
end


    DistractorTypeError=DistractorType(~CorrectChoice);

    Dist_Error_DistractorType=histcounts(DistractorTypeError,DistractortTypeWhole);

    DistractorTypeChoice=DistractorType(ChoiceTrialSel);
 
    Dist_Choice_DistractorType=histcounts(DistractorTypeChoice,DistractortTypeWhole);
    
    Dist_PropError=Dist_Error_DistractorType./Dist_Choice_DistractorType;
    
    
    %Error among different entry
    Entry_Error=Entry(~CorrectChoice);
    Dist_Entry_Error=histcounts(Entry_Error,EntryWhole);
    
    Entry_Choice=Entry(ChoiceTrialSel);
    Dist_Entry_Choice=histcounts(Entry_Choice,EntryWhole);
    
    Dist_PropEntryError=Dist_Entry_Error./Dist_Entry_Choice;
    
    
%Decision time

FixationOffsetTime = FindOutTime(EventChannel,EventTimeChannel,FIXOFF);
DecisionTimePoint_Match = FindOutTime(EventChannel,EventTimeChannel,MATCHOBJ);    
DecisionDur_Match=DecisionTimePoint_Match- FixationOffsetTime(:,1:end-1);

DecisionTimePoint_NonMatch = FindOutTime(EventChannel,EventTimeChannel,NONMATCHOBJ);    
%DecisionDur_NonMatch=DecisionTimePoint_NonMatch- FixationOffsetTime(:,1:end-1);

DecisionTimePoint_Both=fillinNonMatch(DecisionTimePoint_Match,DecisionTimePoint_NonMatch);
DecisionDur_Both=DecisionTimePoint_Both-FixationOffsetTime(:,1:end-1);


DecisionDur_Match_Choice=DecisionDur_Match(ChoiceTrialSel,:);

DecisionDur_Match_Correct=DecisionDur_Match(CorrectChoice,:);

Ave_DecisionDur_Match_Correct=mean(DecisionDur_Match_Correct(:,2:end),2);

CorrectTrialIndex=Selection(CorrectChoice&ChoiceTrialSel);

DecisionDur_Match_Correct_Seq_Mean=mean(DecisionDur_Match_Correct(:,2:end),1);
DecisionDur_Match_Correct_Seq_Sem=std(DecisionDur_Match_Correct(:,2:end),1)./sqrt(3);

%FP Aquisition time
FixationOnsetTime = FindOutTime(EventChannel,EventTimeChannel,SHOWFIXCD );
StartFixationTime = FindOutTime(EventChannel,EventTimeChannel,FIXSTARTCD  );
FPAquisitionTime=StartFixationTime-FixationOnsetTime;
%{
FPAquisitionTime_Correct=FPAquisitionTime(CorrectChoice&ChoiceTrialSel,:);
FPAquisitionTime_Correct_Mean=mean(FPAquisitionTime_Correct,1);
FPAquisitionTime_Correct_Sem=std(FPAquisitionTime_Correct,[],1)./sqrt(size(FPAquisitionTime_Correct,1));
%}
FPAquisitionTime_Choice=FPAquisitionTime(ChoiceTrialSel,:);

%Loop over the FP sequence
FPSeqLength_Choice=sum(~isnan(StartFixationTime(ChoiceTrialSel,:)),2);


%DecisionDur_Both

for i=1:5
    sel=FPSeqLength_Choice==i;
    
    FPAquisitionTime_Choice_curr=FPAquisitionTime_Choice(sel,:);
    
    FPAquisitionTime_Choice_Order_Mean(i,:)=nanmean(FPAquisitionTime_Choice_curr,1);
    FPAquisitionTime_Choice_Order_Sem(i,:)=nanstd(FPAquisitionTime_Choice_curr,1)./sqrt(sum(~isnan(FPAquisitionTime_Choice_curr),1));
    
    
    
    %Decision time 
    if i<=4 & i>=1
        if i==4
             sel=SeqLengthChoice'==i;      
        else
            sel=SeqLengthChoice'==i;
        end
        
       
        tmp=DecisionDur_Both(sel,:);
        
        
        
        Mean_DecisionDur_Both_Seq(i,:)=nanmean(tmp,1);
        Sem_DecisionDur_Both_Seq(i,:)=nanstd(tmp,[],1)./sqrt(sum(~isnan(tmp),1));
        
    end
    
    
    
end



%Plot correct rate and sequence length
ChoiceTrialIndex=Selection(ChoiceTrialSel);
PlotEveryNumTrials=20;%Plot every 5 trials

index=1;

for i=1:PlotEveryNumTrials:length(ChoiceTrialIndex)
    
    currCorrectIndex=CorrectChoice(i:min(i+PlotEveryNumTrials-1,length(ChoiceTrialIndex)));
    
   
    PropCorrectEvery(index)=sum(currCorrectIndex)/(min(i+PlotEveryNumTrials-1,length(ChoiceTrialIndex))-i+1)*100;
    
    TrialEveryIndex(index)=min(i+PlotEveryNumTrials-1,length(ChoiceTrialIndex));
    
    CumulativeCorrectRate(index)=sum(CorrectChoice(1:min(i+PlotEveryNumTrials-1,length(ChoiceTrialIndex))))/min(i+PlotEveryNumTrials-1,length(ChoiceTrialIndex));
    
    p_CumulativeCR(index)=myBinomTest(sum(CorrectChoice(1:min(i+PlotEveryNumTrials-1,length(ChoiceTrialIndex)))),min(i+PlotEveryNumTrials-1,length(ChoiceTrialIndex)),0.125);
    %Ave sequence length
    
     currSequenceIndex=SeqLengthChoice(i:min(i+PlotEveryNumTrials-1,length(ChoiceTrialIndex)));
    
    MeanSeqLengthEvery(index)=nanmean(currSequenceIndex);
    SemSeqLengthEvery(index)=nanstd(currSequenceIndex)/sqrt(TrialEveryIndex(index));
    
    %CorrectRate Each Loc
    tmp=[];
    SeqCounts_Seqcurr=histcounts(currSequenceIndex,SqeAll);
    tmp(1)=sum(SeqCounts_Seqcurr(2:end));
    tmp(2)=sum(SeqCounts_Seqcurr(3:end));
    tmp(3)=sum(SeqCounts_Seqcurr(4:end));
    CorrectRate_SeqEvery(index,:)=tmp/length(currSequenceIndex);
    
    %Conditioned correct rate
    tmp2(1)=tmp(1)/length(currSequenceIndex);
    tmp2(2)=tmp(2)/tmp(1);
    tmp2(3)=tmp(3)/tmp(2);
    ConditionedCorrectRate_SeqEvery(index,:)=tmp2;

    
    
    index=index+1;
end



if ShowFigureFlag
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Figure 1: BehaviorAnalysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figtitlestr{1}='CorrectRateAlongWithTrial';

fig{1}=PrepareFigure(20,'w',[50,100, 1200,800],'Name',figtitlestr{1});

%Overall performance
axes('Position',[0.1,0.7,0.12,0.2]);
bar([1,2],[TotalChoiceTrials,TotalCorrectChoice],'LineWidth',3);
xlim([0.4,2.6]);
box off

ylabel('Number of trials');
set(gca, 'XTickLabel', {'Choices' 'Correct'},'LineWidth',3,'FontSize',15,'FontWeight','Bold');

axes('Position',[0.28,0.7,0.05,0.2]);
bar(1,[ProportionOfCorrectChoice],'LineWidth',3,'Facecolor','r');
xlim([0.4,1.6]);
box off
ylabel('Correct Rate');
set(gca, 'XTickLabel', {'CorrectRate'},'LineWidth',3,'FontSize',15,'FontWeight','Bold');


%Dist_PropEntryError
axes('Position',[0.39,0.7,0.15,0.2]);
bar(EntryWhole(1:end-1),Dist_PropEntryError,'LineWidth',3,'Facecolor','r');
xlim([0.4,4.6]);

box off
ylabel('PropError');
xlabel('Entry Type');
set(gca, 'XTickLabel', string(DistractortTypeWhole(1:end-1)),'LineWidth',3,'FontSize',15,'FontWeight','Bold');




axes('Position',[0.6,0.7,0.3,0.2]);
bar(DistractortTypeWhole(1:end-1),Dist_PropError,'LineWidth',3,'Facecolor','r');

box off
ylabel('PropError');
xlabel('Distractor Type');
set(gca, 'XTickLabel', string(DistractortTypeWhole(1:end-1)),'LineWidth',3,'FontSize',15,'FontWeight','Bold');


%Dist_PropError

axes('Position',[0.1,0.5,0.45,0.13]);
plot(TrialEveryIndex,PropCorrectEvery,'-b','LineWidth',3);
hold on
%Add chance level
ChanceLevel=1/8*100;
plot(TrialEveryIndex,ChanceLevel*ones(length(TrialEveryIndex),1),'--k','LineWidth',2);
yrange=ylim;
yrange(1)=0;
ylim(yrange);
xlabel('Cumulative Trials');
ylabel('Correct Rate(%)');
legend(sprintf('Every %d Trials',PlotEveryNumTrials));
set(gca,'LineWidth',2,'FontSize',15,'FontWeight','Bold');
hold on


axes('Position',[0.1,0.28,0.45,0.15]);

errorbar(TrialEveryIndex,MeanSeqLengthEvery,SemSeqLengthEvery,'-r','LineWidth',2);
xlabel('Cumulative Trials');
ylabel('Ave Seq Length');
legend(sprintf('Every %d Trials',PlotEveryNumTrials));
set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');


axes('Position',[0.1,0.06,0.45,0.15]);

plot(ChoiceTrialIndex,SeqLengthChoice,'-g','LineWidth',2);
xlabel('Trial Index');
ylabel('Seq Length');

set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');


axes('Position',[0.62,0.5,0.3,0.12]);

plot(CorrectTrialIndex,log(Ave_DecisionDur_Match_Correct),'-b','LineWidth',2);
xlabel('Correct Trial Index');
ylabel('Log(RT(Correct))');

set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');

axes('Position',[0.62,0.28,0.3,0.12]);

plot(repmat(CorrectTrialIndex',1,size(DecisionDur_Match_Correct,2)),DecisionDur_Match_Correct,'LineWidth',2);
xlabel('Correct Trial Index');
ylabel('Log(RT(Correct))');

legend(string(1:4));
set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');


axes('Position',[0.62,0.06,0.3,0.12]);


plot(DecisionDur_Match_Correct(:,2:end)','-k','Color',[100 100 100]/255);
hold on
errorbar([1:3],DecisionDur_Match_Correct_Seq_Mean,DecisionDur_Match_Correct_Seq_Sem,'-b','LineWidth',3);

xlabel('Choice Order');
ylabel('RT');
xlim([0.8,3.2]);
set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Figure 2:Sequence Analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figtitlestr{2}='SequenceAnalysis';

fig{2}=PrepareFigure(21,'w',[70,100, 1200,800],'Name',figtitlestr{2});
subplot(3,2,1)
%plot(FPAquisitionTime_Correct','-k','Color',[100 100 100]/255);
%hold on
ColorEach={'m','k','g','b','r'};
for i=1:5
    
    errorbar([1:5],FPAquisitionTime_Choice_Order_Mean(i,:), FPAquisitionTime_Choice_Order_Sem(i,:),'-r','Color',ColorEach{i},'LineWidth',3);
    hold on
end

legend(string([1:5]));

box off
xlim([0.5,5.5]);
xlabel('FP Order');
ylabel('FP Aquisition Time (ms)')
set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
title('FP Aquisition Time (Each Seq Length)');

subplot(3,2,2)
%plot(FPAquisitionTime_Correct','-k','Color',[100 100 100]/255);
%hold on
ColorEach={'k','g','b','r'};
for i=1:4
    
    errorbar([1:4],Mean_DecisionDur_Both_Seq(i,:), Sem_DecisionDur_Both_Seq(i,:),'-r','Color',ColorEach{i},'LineWidth',3);
    hold on
end

legend(string([1:4]));

box off
xlim([0.5,4.5]);
xlabel('Choice Order');
ylabel('DecisionDuration (ms)')
set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
title('Decision Duration (Each Seq Length)');


subplot(3,2,3)
%plot(FPAquisitionTime_Correct','-k','Color',[100 100 100]/255);
%hold on
ColorEach={'k','g','b','r'};
bar([1:3],CorrectRate_Seq);
hold on
plot([0.5,1.5],[0.5,0.5],'--k');
plot([1.5,2.5],[0.25,0.25],'--k');
plot([2.5,3.5],[0.125,0.125],'--k');

%p_Correct
%Add significant level
for i=1:3
    if p_Correct(i)<0.05 & p_Correct(i)>0.01
        plot(i,1,'*','color','r','MarkerSize',15,'LineWidth',3);
        
    elseif p_Correct(i)<0.01 &  p_Correct(i)>0.001
        plot(i,1,'*','color','r','MarkerSize',15,'LineWidth',3);
         plot(i,0.9,'*','color','r','MarkerSize',15,'LineWidth',3);
        
    elseif p_Correct(i)<0.001
         plot(i,1,'*','color','r','MarkerSize',15,'LineWidth',3);
         plot(i,0.9,'*','color','r','MarkerSize',15,'LineWidth',3);
         plot(i,0.8,'*','color','r','MarkerSize',15,'LineWidth',3);
               
    end
end


box off
xlim([0.5,3.5]);
xlabel('Choice Order');
ylabel('Correct Rate');
set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
title('Correct Rate (EachChoice)');

subplot(3,2,4)
%plot(FPAquisitionTime_Correct','-k','Color',[100 100 100]/255);
%hold on
ColorEach={'g','b','r'};
for i=1:3
    plot(TrialEveryIndex,CorrectRate_SeqEvery(:,i),'-k','Color',ColorEach{i},'LineWidth',3);
    hold on
    
end


plot(TrialEveryIndex,0.5*ones(1,length(TrialEveryIndex)),'--k','LineWidth',3);
plot(TrialEveryIndex,0.25*ones(1,length(TrialEveryIndex)),'--k','LineWidth',3);
plot(TrialEveryIndex,0.125*ones(1,length(TrialEveryIndex)),'--k','LineWidth',3);

legend([string([1:3]),'0.5Chance','0.25Chance','0.125Chance']);

box off
%xlim([0.5,3.5]);
xlabel('Cumulative Trial');
ylabel('Correct Rate');
set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
title('Correct Rate (EachChoice)');

subplot(3,2,5)
ColorEach={'k','g','b','r'};
bar([1:3],CorrectRate_Condition);
hold on
plot([0.5,3.5],[0.5,0.5],'--k');
box off

%Add significant level
for i=1:3
    if p_ConditionedCorrect(i)<0.05 & p_ConditionedCorrect(i)>0.01
        plot(i,1,'*','color','r','MarkerSize',15,'LineWidth',3);
        
    elseif p_ConditionedCorrect(i)<0.01 &  p_ConditionedCorrect(i)>0.001
        plot(i,1,'*','color','r','MarkerSize',15,'LineWidth',3);
         plot(i,0.9,'*','color','r','MarkerSize',15,'LineWidth',3);
        
    elseif p_ConditionedCorrect(i)<0.001
         plot(i,1,'*','color','r','MarkerSize',15,'LineWidth',3);
         plot(i,0.9,'*','color','r','MarkerSize',15,'LineWidth',3);
         plot(i,0.8,'*','color','r','MarkerSize',15,'LineWidth',3);
               
    end
end




xlim([0.5,3.5]);
xlabel('Choice Order');
ylabel('Correct Rate');
set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
title('Conditioned Correct Rate');

subplot(3,2,6)
%plot(FPAquisitionTime_Correct','-k','Color',[100 100 100]/255);
%hold on
ColorEach={'g','b','r'};
for i=1:3
    plot(TrialEveryIndex,ConditionedCorrectRate_SeqEvery(:,i),'-k','Color',ColorEach{i},'LineWidth',3);
    hold on
    
end


plot(TrialEveryIndex,0.5*ones(1,length(TrialEveryIndex)),'--k','LineWidth',3);

legend([string([1:3]),'0.5Chance','0.25Chance','0.125Chance']);

box off
%xlim([0.5,3.5]);
xlabel('Cumulative Trial');
ylabel('Conditioned Correct Rate');
set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
title('Conditioned Correct Rate (EachChoice)');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Figure 3:Fractal Based CorrectRate Analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figtitlestr{3}='FractalBasedCorrectRate';

fig{3}=PrepareFigure(22,'w',[90,100, 1400,800],'Name',figtitlestr{3});
subplot(2,4,1)
%Heat Map of all fractals in all temporal locations
h=heatmap(CorrectRate_Obj_TmpLoc);
colormap default;
h.ColorLimits=[0.5,1];
h.Title='CorrectRate perFractal perTemporal Location';
h.YDisplayData=[1,2,3,4];
h.YDisplayLabels={'4to1','1to2','2to3','3to4'};

h.XDisplayData=[1,2,3];
h.XDisplayLabels={'FirstChoice','SecondChoice','ThirdChoice'};
h.FontSize=15;


subplot(2,4,2)
%Heat Map of p_Value of all fractals in all temporal locations

p_ValMat=p_Obj_TmpLoc;
p_ValMat(p_Obj_TmpLoc<=0.001)=1;
p_ValMat(p_Obj_TmpLoc<0.01 & p_Obj_TmpLoc>=0.001)=0.9;
p_ValMat(p_Obj_TmpLoc<0.05 & p_Obj_TmpLoc>=0.01)=0.8;
p_ValMat(p_Obj_TmpLoc>=0.05)=0;


h=heatmap(p_ValMat);
colormap default;
h.ColorLimits=[0.5,1];
h.Title='SigLevel';
h.YDisplayData=[1,2,3,4];
h.YDisplayLabels={'4to1','1to2','2to3','3to4'};

h.XDisplayData=[1,2,3];
h.XDisplayLabels={'FirstChoice','SecondChoice','ThirdChoice'};
h.FontSize=15;


subplot(2,4,3) %Marginal fractal

ColorEach={'k','g','b','r'};
bar([1:4],CorrectRate_Obj);
hold on
plot([0.5,4.5],[0.5,0.5],'--k');
box off

%Add significant level
for i=1:4
    if p_CorrectRate_Obj(i)<0.05 & p_CorrectRate_Obj(i)>0.01
        plot(i,1,'*','color','r','MarkerSize',15,'LineWidth',3);
        
    elseif p_CorrectRate_Obj(i)<0.01 &  p_CorrectRate_Obj(i)>0.001
        plot(i,1,'*','color','r','MarkerSize',15,'LineWidth',3);
         plot(i,0.9,'*','color','r','MarkerSize',15,'LineWidth',3);
        
    elseif p_CorrectRate_Obj(i)<0.001
         plot(i,1,'*','color','r','MarkerSize',15,'LineWidth',3);
         plot(i,0.9,'*','color','r','MarkerSize',15,'LineWidth',3);
         plot(i,0.8,'*','color','r','MarkerSize',15,'LineWidth',3);
               
    end
end



xlim([0.5,4.5]);
xlabel('Fractal Index');
ylabel('Correct Rate');
set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
title('Marginal Fractal Correct Rate');

subplot(2,4,4) %Marginal temporal location

ColorEach={'k','g','b','r'};
bar([1:3],CorrectRate_TmpLoc);
hold on
plot([0.5,3.5],[0.5,0.5],'--k');
box off

%Add significant level
for i=1:3
    if p_CorrectRate_TmpLoc(i)<0.05 & p_CorrectRate_TmpLoc(i)>0.01
        plot(i,1,'*','color','r','MarkerSize',15,'LineWidth',3);
        
    elseif p_CorrectRate_TmpLoc(i)<0.01 &  p_CorrectRate_TmpLoc(i)>0.001
        plot(i,1,'*','color','r','MarkerSize',15,'LineWidth',3);
         plot(i,0.9,'*','color','r','MarkerSize',15,'LineWidth',3);
        
    elseif p_CorrectRate_TmpLoc(i)<0.001
         plot(i,1,'*','color','r','MarkerSize',15,'LineWidth',3);
         plot(i,0.9,'*','color','r','MarkerSize',15,'LineWidth',3);
         plot(i,0.8,'*','color','r','MarkerSize',15,'LineWidth',3);
               
    end
end



xlim([0.5,4.5]);
xlabel('TemporalLocation');
ylabel('Correct Rate');
set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
title('Marginal Temporal Location Correct Rate');


subplot(2,4,5) %CorrectRate relation with distance
%Heat Map of all fractals with all distance distractors
h=heatmap(CorrectRate_Obj_Dist_Obj);
colormap default;
h.ColorLimits=[0.5,1];
h.Title='CorrectRate perFractal perDistance';
h.YDisplayData=[1,2,3,4];
h.YDisplayLabels={'4to1','1to2','2to3','3to4'};

h.XDisplayData=[1,2,3];
h.XDisplayLabels={'Distance1','Distance2','Distance3'};
h.FontSize=15;


%Heat Map of p_Value of all fractals with distractors in all distance
subplot(2,4,6) %CorrectRate relation with distance
p_ValMat=p_CorrectRate_Obj_Dist_Obj;
p_ValMat(p_CorrectRate_Obj_Dist_Obj<=0.001)=1;
p_ValMat(p_CorrectRate_Obj_Dist_Obj<0.01 & p_CorrectRate_Obj_Dist_Obj>=0.001)=0.9;
p_ValMat(p_CorrectRate_Obj_Dist_Obj<0.05 & p_CorrectRate_Obj_Dist_Obj>=0.01)=0.8;
p_ValMat(p_CorrectRate_Obj_Dist_Obj>=0.05)=0;


h=heatmap(p_ValMat);
colormap default;
h.ColorLimits=[0.5,1];
h.Title='SigLevel';
h.YDisplayData=[1,2,3,4];
h.YDisplayLabels={'4to1','1to2','2to3','3to4'};

h.XDisplayData=[1,2,3];
h.XDisplayLabels={'Distance1','Distance2','Distance3'};
h.FontSize=15;



    

subplot(2,4,7) %CorrectRate relation with distance(Marginal)

ColorEach={'k','g','b','r'};
bar(distance_unique,CorrectRate_Obj_Dist);
hold on
plot([0.5,3.5],[0.5,0.5],'--k');
box off

%Add significant level
for i=1:3
    if p_CorrectRate_Obj_Dist(i)<0.05 & p_CorrectRate_Obj_Dist(i)>0.01
        plot(i,1,'*','color','r','MarkerSize',15,'LineWidth',3);
        
    elseif p_CorrectRate_Obj_Dist(i)<0.01 &  p_CorrectRate_Obj_Dist(i)>0.001
        plot(i,1,'*','color','r','MarkerSize',15,'LineWidth',3);
         plot(i,0.9,'*','color','r','MarkerSize',15,'LineWidth',3);
        
    elseif p_CorrectRate_Obj_Dist(i)<0.001
         plot(i,1,'*','color','r','MarkerSize',15,'LineWidth',3);
         plot(i,0.9,'*','color','r','MarkerSize',15,'LineWidth',3);
         plot(i,0.8,'*','color','r','MarkerSize',15,'LineWidth',3);
               
    end
end

xlim([0.5,3.5]);
xticklabels(string(distance_unique));
xlabel('Distractor Distance');
ylabel('Correct Rate');
set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
title('CorrectRate Each Dist');

subplot(2,4,8) %Cumulative Correct Rate for the day
plot(TrialEveryIndex,CumulativeCorrectRate,'-r','LineWidth',3);
hold on
plot(TrialEveryIndex,0.125*ones(1,length(TrialEveryIndex)),'--k','LineWidth',3);

box off

%Add significant level
hplot=p_CumulativeCR;
hplot(p_CumulativeCR<0.05)=1;
hplot(p_CumulativeCR>=0.05)=NaN;

sig=area(TrialEveryIndex,hplot*0.02,'FaceColor','#FFC107','EdgeColor','#FFC107','FaceAlpha',0.5,'EdgeAlpha',0.5);



%xlim([0.5,3.5]);
xlabel('Cumulative Trial');
ylabel('Cumulative CorrectRate');
set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
title('Cumulative Correct Rate (Final)');

end%End of the if show figure flag






end

function ObjectIndex=ReproduceFromEvent(event,code)

for i=1:size(event,1)
    
   tmp=event(i,ismember(event(i,:),code));
    
 %  number(i)=length(tmp);
   ObjectIndex{i}=tmp;
    
end

end

function output=fillinNonMatch(match,nonmatch);
for i=1:size(match,1)
  tmp=match(i,:);
  output(i,:)=tmp;
  nancol=find(isnan(tmp)==1,1,'first');
  if ~isempty(nancol)
      output(i,nancol)=nonmatch(i);
  
  end
 
end

end



