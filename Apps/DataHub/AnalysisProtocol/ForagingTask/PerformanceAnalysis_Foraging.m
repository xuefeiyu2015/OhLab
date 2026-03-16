function PerformanceAnalysis_Foraging(Data, StartTrial, EndTrial,ShowFigureFlag,OutputFlag);

%This function is to summarise the behavior performance during all kinds of
%foraging task.
%Xuefei Yu, 01,07,2020
%Note: this program can not discriminate the block with exactly the same
%scene set and order condition.
%ie. This program can only analysis the unique blocks 
%This program can be used for both good data and the whole data

%Analyze rejection rate and acception rate for every 5 repetitions

%Unpack all the data base and task settings
ECode;
TaskType=Data.TaskType;
TrialNum=EndTrial-StartTrial+1;
Selection=StartTrial:EndTrial;%Select the trials according to the pannel

%Real name of all the objects and scenes.
GoodObjectLib=Data.GoodObjectLib;
BadObjectLib=Data.BadObjectLib;
SceneLib=Data.SceneLib;

%Real name of each scene and object with their index
SceneIndex=Data.SceneIndex(:,2);
TargetIndex=Data.TargetIndex(:,2:end,1);

PlotEveryNumRep=2;%Plot performance every 5 repetitions

if OutputFlag==1
h=questdlg('Output data?',...
    'Yes','No');
waitfor(h);
if strcmp(h,'No')
    OutputFlag=0;
    
end
end




%Calculate number of blocks and the block separation index

%BlockSepTrial=unique(SceneLib(:,1));
BlockType=Data.BlockType(:,:);


BlockType=BlockType(Selection,:);

BlockType(BlockType(:,3)>1,3)=999;

[UniqueBlockType,orderindex,reprojectorderindex]=unique(BlockType(:,2:end),'row');

%Definition of ploting parameters
ColorForEachMode=[0.83 0.14 0.14
             1.00 0.54 0.00
             0.47 0.25 0.80
             0.25 0.80 0.54
             0.1373    0.9333    0.0745
             ];
AssociateBlock=[10,11,0,1,-1];
Blockstr={'FLocFScene','FLocRScene','RLocFScene','RLocRScene','ReverseOrderFScene'};


NumberOfBlocks=max(reprojectorderindex);
%For above: The first column: 0:SceneSetA;1:SceneSetB;
%Second column:0:Fixed,1:Random


%Revisit behavior outcome
ParameterLib=Data.ParaLib;
BehaviorOutCome=ParameterLib('BehaviorOutCome');
BehaviorOutCome=BehaviorOutCome(Selection,:);
ChoiceOutcome=Data.ChoiceOutcome;%This is only for good trials





%Collect information for each unique block
NumberOfTrialsInCurrentBlock=nan*ones(1,NumberOfBlocks);
NumberOfGoodTrialsInCurrentBlock=nan*ones(1,NumberOfBlocks);
NumberOfBreakTrialsInCurrentBlock=nan*ones(1,NumberOfBlocks);
NumberOfNoFixTrialsInCurrentBlock=nan*ones(1,NumberOfBlocks);
NumberOfFailSacTrialsInCurrentBlock=nan*ones(1,NumberOfBlocks);
NumberOfOtherFailTrialsInCurrentBlock=nan*ones(1,NumberOfBlocks);
NumberOfGoodChoiceInCurrentBlock=nan*ones(1,NumberOfBlocks);
NumberOfBadChoiceInCurrentBlock=nan*ones(1,NumberOfBlocks);

BlockTypeInOrder=nan*UniqueBlockType;
SceneNameInCurrentBlock=cell(1,NumberOfBlocks);
GoodObjectNameInCurrentBlock=cell(1,NumberOfBlocks);
BadObjectNameInCurrentBlock=cell(1,NumberOfBlocks);

PropOfGoodChoiceInCurrentBlock=nan*ones(1,NumberOfBlocks);
PropOfBadRejectInCurrentBlock=nan*ones(1,NumberOfBlocks);
PropOfGoodAcceptInCurrentBlock=nan*ones(1,NumberOfBlocks);



for i=1:NumberOfBlocks
    CurrentBlockType=UniqueBlockType(i,:);
    BlockLogicSelection=sum(BlockType(:,2:end)==CurrentBlockType,2)==2;
    TrialIndex=BlockType(BlockLogicSelection,1);
    BlockOrder=find(sort(orderindex)==orderindex(i));
    BlockOrderSave(i)=find(CurrentBlockType(:,1)==AssociateBlock);
    
    PlotOrder(i)=find(AssociateBlock==CurrentBlockType(1));
    
    NumberOfTrialsInCurrentBlock(BlockOrder)=length(TrialIndex);%Number of total trials in this block
    BlockTypeInOrder(BlockOrder,:)=CurrentBlockType;%Store the current block type in the original order
    
    BehaviorOutComeCurrentBlock=BehaviorOutCome(ismember(BehaviorOutCome(:,1),TrialIndex),2);%Load the behavioroutcome
    
    NumberOfGoodTrials=sum(BehaviorOutComeCurrentBlock==CORRECTCD);
    NumberOfBreakTrials=sum(BehaviorOutComeCurrentBlock==BROKEFIX);
    NumberOfNoFixTrials=sum(BehaviorOutComeCurrentBlock==NOFIX);
    NumberOfFailSacTrials=sum(BehaviorOutComeCurrentBlock==FAILED_SAC);
    NumberOfOtherFailTrials=NumberOfTrialsInCurrentBlock(BlockOrder)-NumberOfGoodTrials-NumberOfBreakTrials-NumberOfNoFixTrials-NumberOfFailSacTrials;
 
    NumberOfGoodTrialsInCurrentBlock(BlockOrder)=NumberOfGoodTrials;
    NumberOfBreakTrialsInCurrentBlock(BlockOrder)=NumberOfBreakTrials;
    NumberOfNoFixTrialsInCurrentBlock(BlockOrder)=NumberOfNoFixTrials;
    NumberOfFailSacTrialsInCurrentBlock(BlockOrder)=NumberOfFailSacTrials;
    NumberOfOtherFailTrialsInCurrentBlock(BlockOrder)=NumberOfOtherFailTrials;
    
    ChoiceOutcomeBlock=ChoiceOutcome(ismember(ChoiceOutcome(:,1),TrialIndex),2);
    NumberOfBadChoiceInCurrentBlock(BlockOrder)=sum(ChoiceOutcomeBlock==BADCHOICECD);
    NumberOfGoodChoiceInCurrentBlock(BlockOrder)=sum(ChoiceOutcomeBlock==GOODCHOICECD);
    
    PropOfGoodChoiceInCurrentBlock(BlockOrder)=sum(ChoiceOutcomeBlock==GOODCHOICECD)/length(ChoiceOutcomeBlock)*100;
    
    %The scene and target name for the current block
    SceneNameInCurrentBlock{BlockOrder}=unique(SceneLib(ismember(SceneLib(:,1),TrialIndex),3),'stable');
    GoodObjectNameInCurrentBlock{BlockOrder}=unique(GoodObjectLib(ismember(GoodObjectLib(:,1),TrialIndex),3),'stable');
    BadObjectNameInCurrentBlock{BlockOrder}=unique(BadObjectLib(ismember(BadObjectLib(:,1),TrialIndex),3),'stable');
  
    UniqueSceneIndexInCurrentBlock=unique(SceneIndex( BlockLogicSelection,:));
     UniqueTargetIndexInCurrentBlock=unique(TargetIndex( BlockLogicSelection,:));
     UniqueTargetIndexInCurrentBlock=  UniqueTargetIndexInCurrentBlock(~isnan(  UniqueTargetIndexInCurrentBlock));
     
     SceneIndexCurrBlock=SceneIndex( BlockLogicSelection,:);
     TargetIndexCurrBlock=TargetIndex( BlockLogicSelection,:);
     
     %Overall performance of good acception and bad rejection
     
     GoodObjectIndexCurrBlock=mod(TargetIndexCurrBlock,2);
     BadObjectIndexCurrBlock=mod(TargetIndexCurrBlock,2)==0;
     
     NumOfGoodObjCurrBlock=sum(GoodObjectIndexCurrBlock,'all','omitnan' );
     NumOfBadObjCurrBlock=sum(BadObjectIndexCurrBlock,'all','omitnan');
     
     RejectAcceptOutcomeCurrBlock=TransformToOutcome(TargetIndexCurrBlock);%1:Accept;0:Reject
     
     BadRejectCurrBlock=RejectAcceptOutcomeCurrBlock==0 & BadObjectIndexCurrBlock;
     GoodAcceptCurrBlock=RejectAcceptOutcomeCurrBlock==1 & GoodObjectIndexCurrBlock==1;
     
     PropOfBadRejectInCurrentBlock(BlockOrder)=sum(BadRejectCurrBlock,'all','omitnan')/NumOfBadObjCurrBlock*100;
     
     PropOfGoodAcceptInCurrentBlock(BlockOrder)=sum(GoodAcceptCurrBlock,'all','omitnan')/NumOfGoodObjCurrBlock*100;
     
     
     %Number of repetition
     RepNum=ceil(sum(BlockLogicSelection)/length(UniqueSceneIndexInCurrentBlock));
     
     %Number of condition in each repetition
    [temp_loc,NumInRep,RepinScene]=FindNumOfRepetition(SceneIndexCurrBlock);
    
    RepStart=1;
    
    NumGoodChoiceRep=[];
    NumBadChoiceRep=[];
    NumGoodObjectRep=[];
    NumBadObjectRep=[];
    NumGoodAcceptRep=[];
    NumBadAcceptRep=[];
    NumGoodRejectRep=[];
    NumBadRejectRep=[];
    
    
    
    for j=1:RepNum %Loop over the repetition
        RepSelection=RepStart:RepStart+NumInRep(j)-1;
        ChoiceOutcomeRep=ChoiceOutcomeBlock(RepSelection);
        NumGoodChoiceRep(j)=sum(ChoiceOutcomeRep==GOODCHOICECD);
        NumBadChoiceRep(j)=sum(ChoiceOutcomeRep==BADCHOICECD);
        
        RejectAcceptOutcome=TransformToOutcome(TargetIndexCurrBlock(RepSelection,:));%1:Accept;0:Reject
        
        RejectAcceptOutcome_Reverse=RejectAcceptOutcome;
        RejectAcceptOutcome_Reverse( ~isnan(RejectAcceptOutcome_Reverse))=~RejectAcceptOutcome_Reverse( ~isnan(RejectAcceptOutcome_Reverse));
        
        GoodObjectIndex=mod(TargetIndexCurrBlock(RepSelection,:),2);
        BadObjectIndex=mod(TargetIndexCurrBlock(RepSelection,:),2)==0;
        
        NumGoodObjectRep(j)=sum(sum(mod(TargetIndexCurrBlock(RepSelection,:),2),2,'omitnan' ));
        
        NumBadObjectRep(j)=sum(sum(mod(TargetIndexCurrBlock(RepSelection,:),2)==0,2,'omitnan' ));
        
        GoodAcceptIndex=RejectAcceptOutcome.*GoodObjectIndex;
        
        GoodRejectIndex=RejectAcceptOutcome_Reverse.*GoodObjectIndex;
        
        BadAcceptIndex=RejectAcceptOutcome.*BadObjectIndex;
        
       BadRejectIndex=RejectAcceptOutcome_Reverse.*BadObjectIndex;
       
       NumGoodAcceptRep(j)=sum(sum(GoodAcceptIndex,2,'omitnan'));
       NumBadAcceptRep(j)=sum(sum(BadAcceptIndex,2,'omitnan'));
       NumGoodRejectRep(j)=sum(sum(GoodRejectIndex,2,'omitnan'));
       NumBadRejectRep(j)=sum(sum(BadRejectIndex,2,'omitnan'));
        
        
        RepStart=RepStart+NumInRep(j);

    end
    
    PerformancePerRepBlock{i}=[NumInRep',NumGoodChoiceRep',NumBadChoiceRep',NumGoodObjectRep',NumBadObjectRep',NumGoodAcceptRep',NumBadAcceptRep',NumGoodRejectRep',NumBadRejectRep'];
    
    %For plot, every num repetitions
    
    NumInRep_Fill=[NumInRep,NaN*ones(1,PlotEveryNumRep*ceil(length(NumInRep)/PlotEveryNumRep)-length(NumInRep))];
    NumGoodChoiceRep_Fill=[NumGoodChoiceRep,NaN*ones(1,PlotEveryNumRep*ceil(length(NumInRep)/PlotEveryNumRep)-length(NumInRep))];
    
    NumBadChoiceRep_Fill=[NumBadChoiceRep,NaN*ones(1,PlotEveryNumRep*ceil(length(NumInRep)/PlotEveryNumRep)-length(NumInRep))];
    
    NumGoodObjectRep_Fill=[NumGoodObjectRep,NaN*ones(1,PlotEveryNumRep*ceil(length(NumInRep)/PlotEveryNumRep)-length(NumInRep))];
    NumBadObjectRep_Fill=[NumBadObjectRep,NaN*ones(1,PlotEveryNumRep*ceil(length(NumInRep)/PlotEveryNumRep)-length(NumInRep))];
    
    NumGoodAcceptRep_Fill=[NumGoodAcceptRep,NaN*ones(1,PlotEveryNumRep*ceil(length(NumInRep)/PlotEveryNumRep)-length(NumInRep))];
    NumBadAcceptRep_Fill=[NumBadAcceptRep,NaN*ones(1,PlotEveryNumRep*ceil(length(NumInRep)/PlotEveryNumRep)-length(NumInRep))];
    NumGoodRejectRep_Fill=[NumGoodRejectRep,NaN*ones(1,PlotEveryNumRep*ceil(length(NumInRep)/PlotEveryNumRep)-length(NumInRep))];
    NumBadRejectRep_Fill=[NumBadRejectRep,NaN*ones(1,PlotEveryNumRep*ceil(length(NumInRep)/PlotEveryNumRep)-length(NumInRep))];
    
    
    
    NumInNRep=sum(reshape(NumInRep_Fill,PlotEveryNumRep,ceil(length(NumInRep)/PlotEveryNumRep)),1,'omitnan');
    
    NumGoodChoiceInNRep=sum(reshape(NumGoodChoiceRep_Fill,PlotEveryNumRep,ceil(length(NumInRep)/PlotEveryNumRep)),1,'omitnan');
    NumBadChoiceInNRep=sum(reshape(NumBadChoiceRep_Fill,PlotEveryNumRep,ceil(length(NumInRep)/PlotEveryNumRep)),1,'omitnan');
    
    NumGoodObjectInNRep=sum(reshape(NumGoodObjectRep_Fill,PlotEveryNumRep,ceil(length(NumInRep)/PlotEveryNumRep)),1,'omitnan');
    NumBadObjectInNRep=sum(reshape(NumBadObjectRep_Fill,PlotEveryNumRep,ceil(length(NumInRep)/PlotEveryNumRep)),1,'omitnan');
    
    NumGoodAcceptInNRep=sum(reshape(NumGoodAcceptRep_Fill,PlotEveryNumRep,ceil(length(NumInRep)/PlotEveryNumRep)),1,'omitnan');
    NumBadAcceptInNRep=sum(reshape(NumBadAcceptRep_Fill,PlotEveryNumRep,ceil(length(NumInRep)/PlotEveryNumRep)),1,'omitnan');
    
     NumGoodRejectInNRep=sum(reshape(NumGoodRejectRep_Fill,PlotEveryNumRep,ceil(length(NumInRep)/PlotEveryNumRep)),1,'omitnan');
     NumBadRejectInNRep=sum(reshape(NumBadRejectRep_Fill,PlotEveryNumRep,ceil(length(NumInRep)/PlotEveryNumRep)),1,'omitnan');
    
    
    
    PropGoodChoiceInNRep{i}=NumGoodChoiceInNRep./NumInNRep;
    PropBadChoiceInNRep{i}=NumBadChoiceInNRep./NumInNRep;
    
    PropGoodAcceptInNRep{i}=NumGoodAcceptInNRep./NumGoodObjectInNRep;
    PropBadAcceptInNRep{i}=NumBadAcceptInNRep./NumBadObjectInNRep;
    
    PropGoodRejectInNRep{i}=NumGoodRejectInNRep./NumGoodObjectInNRep;
    PropBadRejectInNRep{i}=NumBadRejectInNRep./NumBadObjectInNRep;
    
    
end

if ShowFigureFlag
%Bar plot for each block

if ishandle(103)
    close 103
end

titlestr{1}='Behavior Summary For Foraging Task';

fig{1}=figure(103);

set(103,'color','w','Position', [100,50 1300,800],'Name',titlestr{1})
str=strcat('Behavior analysis of ',TaskType);
title(str)
axis off

axes('position',[0.2,0.7,0.3,0.2])
bar(NumberOfTrialsInCurrentBlock)
hold on
if NumberOfBlocks==1
bar(NumberOfGoodTrialsInCurrentBlock+NumberOfBreakTrialsInCurrentBlock+NumberOfNoFixTrialsInCurrentBlock,'r')
hold on
bar(NumberOfGoodTrialsInCurrentBlock+NumberOfBreakTrialsInCurrentBlock,'g')
bar(NumberOfGoodTrialsInCurrentBlock,'b')
legend('Total Trials','NoFixTrials','BreakTrials','ValidTrials')
else
bar([NumberOfGoodTrialsInCurrentBlock',NumberOfBreakTrialsInCurrentBlock',NumberOfNoFixTrialsInCurrentBlock'],'stack')
legend('Total Trials','ValidTrials','BreakTrials','NoFixTrials')
end

xlabel('Block Number')
ylabel('Number of trials')
box off

axes('position',[0.6,0.7,0.3,0.2])
if NumberOfBlocks==1
bar(NumberOfGoodChoiceInCurrentBlock+NumberOfBadChoiceInCurrentBlock,'b')
hold on
bar(NumberOfGoodChoiceInCurrentBlock,'r')
legend('Bad Choice','Good Choice')
else
bar([NumberOfGoodChoiceInCurrentBlock',NumberOfBadChoiceInCurrentBlock'],'stack')
legend('Good Choice','Bad Choice')
end
xlabel('Block Number')
ylabel('Number of trials in valid trials')
box off


data=[NumberOfTrialsInCurrentBlock;NumberOfGoodTrialsInCurrentBlock;NumberOfBreakTrialsInCurrentBlock;NumberOfNoFixTrialsInCurrentBlock;NumberOfFailSacTrialsInCurrentBlock;NumberOfOtherFailTrialsInCurrentBlock;NumberOfGoodChoiceInCurrentBlock;NumberOfBadChoiceInCurrentBlock;BlockTypeInOrder'];
row_name={'Total Trials','ValidTrials','BreakTrials','NoFixTrials','FailSacTrials','OtherFailTrials','Good Choice','Bad Choice','SequenceType','SceneSet'};
column_name=[];
for i=1:NumberOfBlocks
    column_name=[column_name,{strcat('Block',num2str(i))}];
end
%column_name={'Block1','Block2','Block3','Block4'};
uitable(gcf,'Data',data,'Position',[200 200 500 200],'Columnname',column_name,'Rowname',row_name);

annotation('textbox',[0.58,0.45,0.3,0.2],'String','Scene Name','FitBoxToText','on')
annotation('textbox',[0.6,0.4,0.3,0.2],'String',SceneNameInCurrentBlock,'FitBoxToText','on')
%;GoodObjectNameInCurrentBlock;BadObjectNameInCurrentBlock]
annotation('textbox',[0.67,0.45,0.3,0.2],'String','Good Object Name','FitBoxToText','on')
annotation('textbox',[0.7,0.4,0.3,0.2],'String',GoodObjectNameInCurrentBlock,'FitBoxToText','on')
annotation('textbox',[0.78,0.45,0.3,0.2],'String','Bad Object Name','FitBoxToText','on')
annotation('textbox',[0.8,0.4,0.3,0.2],'String',BadObjectNameInCurrentBlock,'FitBoxToText','on')


%%%%%%%%%%%%%%%%%%%%%%%%%%%Plot performance every N repetition
%%%%%%%%%%%%%%%%%%%%%%%%%%%
titlestr{2}='LearningCurveOfChoiceOutcome';

fig{2}=PrepareFigure(104,'w',[200,100, 1000,600],'Name',titlestr{2}); 

%First Choice outcome
%subplot(1,2,1);

hold on
p=cellfun(@(x) plot([PlotEveryNumRep:PlotEveryNumRep:length(x)*PlotEveryNumRep],x),PropGoodChoiceInNRep,'uniform',0);
p_bad=cellfun(@(x) plot([PlotEveryNumRep:PlotEveryNumRep:length(x)*PlotEveryNumRep],x),PropBadChoiceInNRep,'uniform',0);
for i=1:length(p)
    p{i}.Color=ColorForEachMode(PlotOrder(i),:);
    p{i}.LineWidth=3;
    p{i}.Marker='s';
    p{i}.MarkerSize=5;
    
    p_bad{i}.Color=ColorForEachMode(PlotOrder(i),:);
    p_bad{i}.LineWidth=2;
    p_bad{i}.Marker='.';
    p_bad{i}.MarkerSize=5;
    p_bad{i}.LineStyle='--';
    
end
legend(Blockstr(PlotOrder));
ylim([0,1]);
xlabel(sprintf('Every %d repetition',PlotEveryNumRep));
ylabel('Proportion of Good/Bad Choices');
set(gca,'LineWidth',3,'FontSize',20,'FontWeight','Bold');



%%%%%%%%%%%%%%%%%%%%%%%%%%For Good and Bad Object Separately
titlestr{3}='LearningCurveOfGoodAcceptAndBadReject';

fig{3}=PrepareFigure(105,'w',[200,100, 1000,600],'Name',titlestr{3}); 

for i=1:NumberOfBlocks

subplot(2,ceil(NumberOfBlocks/2),i)

plot([PlotEveryNumRep:PlotEveryNumRep:length(PropGoodAcceptInNRep{i})*PlotEveryNumRep], PropGoodAcceptInNRep{i},'-sr','LineWidth',3);
hold on
plot([PlotEveryNumRep:PlotEveryNumRep:length(PropGoodAcceptInNRep{i})*PlotEveryNumRep], PropBadRejectInNRep{i},'-sb','LineWidth',3);
plot([PlotEveryNumRep:PlotEveryNumRep:length(PropGoodAcceptInNRep{i})*PlotEveryNumRep], PropGoodRejectInNRep{i},'--.r','LineWidth',3);
plot([PlotEveryNumRep:PlotEveryNumRep:length(PropGoodAcceptInNRep{i})*PlotEveryNumRep], PropBadAcceptInNRep{i},'--.b','LineWidth',3);

%xticklabels([PlotEveryNumRep:PlotEveryNumRep:length(PropGoodAcceptInNRep{i})*PlotEveryNumRep]);
xticks([PlotEveryNumRep:PlotEveryNumRep:length(PropGoodAcceptInNRep{i})*PlotEveryNumRep]);
ylim([0,1]);

ylabel('Proportion of choices')
if i==1
g=legend({'GoodAccept','Bad Reject','GoodReject','BadAccept'});
end
xlabel(sprintf('Every %d repetition',PlotEveryNumRep));
set(gca,'LineWidth',3,'FontSize',20,'FontWeight','Bold');

end

%Bar Plot of the good choice rate
%%%%%%%%%%%%%%%%%%%%%%%%%%For Good and Bad Object Separately
titlestr{4}='Proportion of Good Choices';

fig{4}=PrepareFigure(106,'w',[200,100, 1000,400],'Name',titlestr{4});
subplot(1,3,1)
b=bar(PropOfGoodChoiceInCurrentBlock);
b.FaceColor='flat';
b.EdgeColor=b.FaceColor;

b.CData=ColorForEachMode(BlockOrderSave,:);

b.BarWidth=0.8;
xticks([1,2]);

x_labels={Blockstr{BlockOrderSave}};
xticklabels(x_labels);
ylabel('Proportion of Good Choices');
set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
box off

subplot(1,3,2)
b=bar(PropOfBadRejectInCurrentBlock);
b.FaceColor='flat';
b.EdgeColor=b.FaceColor;

b.CData=ColorForEachMode(BlockOrderSave,:);

b.BarWidth=0.8;
xticks([1,2]);

x_labels={Blockstr{BlockOrderSave}};
xticklabels(x_labels);
ylabel('Proportion of Bad Rejection');
set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
box off

subplot(1,3,3)
b=bar(PropOfGoodAcceptInCurrentBlock);
b.FaceColor='flat';
b.EdgeColor=b.FaceColor;

b.CData=ColorForEachMode(BlockOrderSave,:);

b.BarWidth=0.8;
xticks([1,2]);

x_labels={Blockstr{BlockOrderSave}};
xticklabels(x_labels);
ylabel('Proportion of Good Acception');
set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
box off


%Bar Plot of the Bad Reject Rate and Good Accept Rate
%%%%%%%%%%%%%%%%%%%%%%%%%%For Good and Bad Object Separately




%PropOfGoodChoiceInCurrentBlock(BlockOrder)
%PropOfBadRejectInCurrentBlock(BlockOrder)=sum(BadRejectCurrBlock,'all','omitnan')/NumOfBadObjCurrBlock;
     
 %    PropOfGoodAcceptInCurrentBlock(BlockOrder)=sum(GoodAcceptCurrBlock,'all','omitnan')/NumOfGoodObjCurrBlock;

end %End of the if plotfigure

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Export Files
%Output to files
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
TaskShortFlag=find(TaskType=='_',1,'last');
TaskShort=TaskType(TaskShortFlag+1:end);
OutputFolerName=strcat('TrainingData_',TaskShort);
OutPath=strcat(BasicDirectory,'Results',BasicDirectory(end),OutputFolerName);
if ~exist(OutPath)
    mkdir(OutPath)
end
cd(OutPath);
FileName=Data.FileName;
FileName=FileName(1:end-4);


%Get File Information
num = regexp(FileName, '\d+', 'match');
RecordDateOriginal=cell2mat(num(1));


%To get the monkey name
RemainFiles=erase(FileName,RecordDateOriginal);
Delimiter=find(isstrprop(RemainFiles,'upper')==1);
MonkeyName=RemainFiles(Delimiter(1):Delimiter(2)-1);


%Output file name
OutputFileName=sprintf('%s_%s_%s.mat',MonkeyName,string(RecordDateOriginal),TaskShort);
%OutputFigureName=sprintf('%s_%s_N%s_C%s',MonkeyName,string(RecordDateOriginal));


%%%%%%%%%Save overall results %%%%%%%%%%%%%%%%%%%%%%%
%Organize the data
NameStr={'PlotEveryNumRep','PropGoodChoiceInNRep','PropBadChoiceInNRep','PropGoodAcceptInNRep','PropBadAcceptInNRep','PropGoodRejectInNRep','PropBadRejectInNRep'

};

DataLib={PlotEveryNumRep,PropGoodChoiceInNRep,PropBadChoiceInNRep,PropGoodAcceptInNRep,PropBadAcceptInNRep,PropGoodRejectInNRep,PropBadRejectInNRep
};


DataStamp=containers.Map(NameStr,DataLib);

OutputData.Performance.Task=TaskType;
OutputData.Performance.DataStamp=DataStamp;

if OutputFlag
   cd(OutPath);
    save(OutputFileName,'OutputData');
    disp('Data Exported');  
else
    disp('No data export');
end





end