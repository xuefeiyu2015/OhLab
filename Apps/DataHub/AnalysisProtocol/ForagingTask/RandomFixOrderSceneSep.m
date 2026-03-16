function  RandomFixOrderSceneSep(Data, StartTrial, EndTrial,ShowFigureFlag,OutputFlag);
%Function to separate scene tuning with temporal location tuning
%Xuefei Yu 04192020

%Xuefei Yu, 02252020
%Parameter to define the sliding window to measure firing rates
StepSize=10;%in ms
BinWidth=50;%in ms

%Parameter seting for interval selecting plot period
StartOffset=500;%in ms
StopOffset=50;%in ms
%For scene offset
StartOffset2=500;%in ms
StopOffset2=50;%in ms
IntervalAfterMarkerStart=1000;%in ms;%This for align marker
IntervalAfterMarkerStartOff=1000;%in ms;%For scene offset

%Color for fix and rand
Color_Fix=[	251,97,0]/255;
Color_Rand=[100,143,255]/255;

%Color Definition for each scenes

Color=[lines(6);prism(6);pink(6)];
%Color(2,:)=copper(1)+0.5;
Color=[Color;hsv(10)];
Color(18,:)=[0.1,0.457,0.194];



ColorFix=repmat([1,0,0],26,1);
ColorRand=repmat([0,0,1],26,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Unpack all the data base and task settings
ECode;
FileName=Data.FileName;
TaskType=Data.TaskType;
TrialNum=EndTrial-StartTrial+1;
Selection=StartTrial:EndTrial;

%FixationPointIndex=Data.FixPointIndex(Selection,:);
%FixLoc=Data.FixLoc(Selection,:,:);


TargetIndex=Data.TargetIndex(Selection,:,:);
BlockType=Data.BlockType(Selection,:,:);%0:fixed sequence;1;random sequence; 0:ScenesetA;1:SceneSetB
%BlockType(BlockType(:,3)>2,:)=999;

%%
SceneIndex=Data.SceneIndex(Selection,:);

%EyeChannel_X=Data.EyeChannel_X_L(Selection,:); %Only import the first two channels
%EyeChannel_Y=Data.EyeChannel_Y_L(Selection,:);
%EyeBinWidth=Data.EyeBinWidth;
%EyeTime=Data.EyeTime(Selection,:);

EventChannel=Data.EventChannel(Selection,:);
EventTimeChannel=Data.EventTimeChannel(Selection,:);
EventBin=Data.EventBin(Selection,:);

SpikeChannel=Data.SpikeChannel;
SpikeTime=squeeze(Data.SpikeTimeData(:,:,Selection))';
SpikeBin=squeeze(Data.SpikeBinData(:,:,Selection))';




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Get the PSTH for all neurons
SegEvent=[SCENE_ON];%Event Flow of interest
%SegEventEnd=[SHOWFIXCD,FIXOFF,TGTOFF,RWDOFFCD ,NaN,EXTRAFP];%Event end
SegEventEnd=[SHOWFIXCD];%Event end
SegEventStr={'SceneOn'};



PlotStartOffset=[200];%in ms
PlotStopOffset=[100];%in ms

ResponseLatencyThreshold=3;

%Parameter to define the sliding window to measure firing rates
StepSize=10;%in ms
BinWidth=50;%in ms

%Color Definition for each temporal order/ scene

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
TrialNum=EndTrial-StartTrial+1;
Selection=StartTrial:EndTrial;



SceneIndex=Data.SceneIndex(Selection,:);



EventChannel=Data.EventChannel(Selection,:);
EventTimeChannel=Data.EventTimeChannel(Selection,:);
%EventBin=Data.EventBin(Selection,:);

SpikeChannel=Data.SpikeChannel;
SpikeTime=squeeze(Data.SpikeTimeData(:,:,Selection))';
%SpikeBin=squeeze(Data.SpikeBinData(:,:,Selection))';

%Segregate spike data from scene onset to fixation onset

SegmentData=EventSegTool(SegEvent,SegEventEnd,PlotStartOffset,PlotStopOffset,EventChannel,EventTimeChannel,SpikeTime,StepSize,BinWidth);

%Retrive data from the field
%Total average for each event
%MeanFR=SegmentData.MeanFR; %[Mean,sem]%Averall mean and sem
MeanFRRaw=SegmentData.MeanFRRaw{1};%=AverageFiringRate_Event;

%PSTH for each event
MeanPSTH=SegmentData.MeanPSTH;%=[PSTH_Event_Time',PSTH_Event'];
MeanPSTH_Time=MeanPSTH{1};
MeanPSTH_FR=MeanPSTH{2};


%Relative time bin for each event
RelativeTimeMarker=SegmentData.RelativeTimeMarker{1};%=EventTimeMarker;
%Absolute time marker for each event
AbosoluteTimeMarker=SegmentData.AbsoluteTimeMarker{1};




%#################Separate trials according to scene name and block type

%identify block type

UniqueBlockType=unique(BlockType(:,2:end),'row','stable');

RandFixFlag=unique(UniqueBlockType(:,1),'stable');

if length(RandFixFlag)==1
    strFile=sprintf('No random/fix comparison data in the file %s\n',FileName);
    disp(strFile);
    disp('You may need to change to Protocol Scene Tuning');
    return  
end


FixSeqTrials=BlockType(:,2)==0;
RandSeqTrials=BlockType(:,2)==1;


FixSeqSceneName=SceneIndex(FixSeqTrials,3);
RandSeqSceneName=SceneIndex(RandSeqTrials,3);




FixSeqCommonPart=ismember(FixSeqSceneName,unique(RandSeqSceneName));
ConmmonSceneName_FixRand=FixSeqSceneName(FixSeqCommonPart);

RandomSeqCommonPart=ismember(RandSeqSceneName,unique(FixSeqSceneName));
ConmmonSceneName_Rand=RandSeqSceneName(RandomSeqCommonPart);

if isempty(ConmmonSceneName_FixRand)
     strFile=sprintf('No overlap scenes between Fix Seq and Rand Seq %s\n',FileName);
      disp(strFile);
    disp('You may need to change to Protocol Scene Tuning');
    return
end








%%%%%%%%%%%%%%%%%%%%%%%First, overall responses%%%%%%%%%%%%%%%%%%%%%%%%%%

MeanFR_Fix_trial_temp=MeanFRRaw(FixSeqTrials);
MeanFR_Rand_trial_temp=MeanFRRaw(RandSeqTrials);


MeanFR_Fix_trial=MeanFR_Fix_trial_temp(FixSeqCommonPart);%FR for fix sequence
MeanFR_Rand_trial=MeanFR_Rand_trial_temp(RandomSeqCommonPart);%FR for random sequence



%%%%%%%%%%%%%%%%%%%%%%%Then, PSTH difference %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

PSTH_Fix_trial_temp=MeanPSTH_FR(FixSeqTrials,:);
PSTH_Rand_trial_temp=MeanPSTH_FR(RandSeqTrials,:);


PSTH_Fix_trial=PSTH_Fix_trial_temp(FixSeqCommonPart,:);
PSTH_Rand_trial=PSTH_Rand_trial_temp(RandomSeqCommonPart,:);



PSTH_Time=nanmean(MeanPSTH_Time,1);




%%%%%%%%%%%%%%%%%%%%%%%%% Next scene tuning and temporal location
UniqueScene=unique(ConmmonSceneName_FixRand,'stable');
NumberOfUniqueScene=length(UniqueScene);

%Number of condition in each repetition
[temp_locFix,NumInRepFix,RepinSceneFix]=FindNumOfRepetition(ConmmonSceneName_FixRand);
[temp_locRand,NumInRepRand,RepinSceneRand]=FindNumOfRepetition(ConmmonSceneName_Rand);

RandTrialOrder=NaN*size(RepinSceneRand);
FixTrialOrder=NaN*size(RepinSceneFix);

SceneTransform_Rand=NaN*ones(size(ConmmonSceneName_Rand));


%Group trials into scene and temporal order
for i=1:NumberOfUniqueScene %Loop over scene
    RandTrialSelection=ConmmonSceneName_Rand==UniqueScene(i);
    FixTrialSelection=ConmmonSceneName_FixRand==UniqueScene(i);
    %Average
    FR_Scene_Rand_Whole=MeanFR_Rand_trial(RandTrialSelection);
    FR_Scene_Fix_Whole=MeanFR_Fix_trial(FixTrialSelection);
    %PSTH
    PSTH_Fix=PSTH_Fix_trial(FixTrialSelection,:);
    PSTH_Rand=PSTH_Rand_trial(RandTrialSelection,:);
    %Order
    [RandTrialOrder_temp,~,~]=find(RepinSceneRand==UniqueScene(i)); 
    [FixTrialOrder_temp,~,~]=find(RepinSceneFix==UniqueScene(i)); 
    
    RandTrialOrder(1:length(RandTrialOrder_temp),i)=RandTrialOrder_temp;
    FixTrialOrder(1:length(FixTrialOrder_temp),i)=FixTrialOrder_temp;
    
    %Transfrom
    SceneTransform_Rand(RandTrialSelection)=i;
        
    for j=1:NumberOfUniqueScene %Loop over the repetition
         
         
          RandTrialOrder_Selection=RandTrialOrder_temp==j;
          
          
          FixTrialOrder_Selection=FixTrialOrder_temp==j;
          
          
          %Average 
          FR_Scene_Rand_Whole_Order{j,i}=FR_Scene_Rand_Whole(RandTrialOrder_Selection);
          FR_Scene_Fix_Whole_Order{j,i}=FR_Scene_Fix_Whole(FixTrialOrder_Selection);
          %PSTH
        PSTH_Fix_Order{j,i}=PSTH_Fix(FixTrialOrder_Selection,:);
        PSTH_Rand_Order{j,i}=PSTH_Rand(RandTrialOrder_Selection,:);
        
        
    end
    
    
    RelativeTimeMarkerForeachScene(i)=min(RelativeTimeMarker(SceneIndex(:,3)==UniqueScene(i),2));
 
end


%Fit the average firing rate
%using original form;


RepinSceneRand_fit=ConmmonSceneName_Rand;
RepinSceneFix_fit=ConmmonSceneName_FixRand;

OrderFix_fit=temp_locFix(1:length(RepinSceneFix_fit))';
OrderRand_fit=temp_locRand(1:length(RepinSceneRand_fit))';


MeanFR_Rand_trial_fit=MeanFR_Rand_trial;
MeanFR_Fix_trial_fit=MeanFR_Fix_trial;


p = anovan(MeanFR_Rand_trial_fit,[OrderRand_fit,RepinSceneRand_fit],'model','full','display','off');
p_Order=p(1)
p_Scene=p(2)
p_Interact=p(3)



%Control for the order to measure the scene tuning 
[rho_rand_scene,pval_rand_scene] = partialcorr(MeanFR_Rand_trial_fit,RepinSceneRand_fit,OrderRand_fit)
%Control the scene to measure the order tuning
[rho_rand_order,pval_rand_order] = partialcorr(MeanFR_Rand_trial_fit,OrderRand_fit,RepinSceneRand_fit)






%Normalize 

MeanFR_Rand_Order_Norm=MeanFR_Rand_trial;
MeanFR_Rand_Scene_Norm=MeanFR_Rand_trial;
PSTH_Rand_Scene_Norm=PSTH_Rand_trial;
PSTH_Rand_Order_Norm=PSTH_Rand_trial;



 for i=1:NumberOfUniqueScene %Loop over the scene
     SceneFR_mean(i)=nanmean(MeanFR_Rand_trial_fit(ConmmonSceneName_Rand==UniqueScene(i)));
     
     SceneFR_std(i)=nanstd(MeanFR_Rand_trial_fit(ConmmonSceneName_Rand==UniqueScene(i)));
     
     ScenePSTH_mean=nanmean(PSTH_Rand_trial(ConmmonSceneName_Rand==UniqueScene(i),:));
     ScenePSTH_std=nanstd(PSTH_Rand_trial(ConmmonSceneName_Rand==UniqueScene(i),:));
     
     
     %Normed fr
     MeanFR_Rand_Order_Norm(ConmmonSceneName_Rand==UniqueScene(i))= (MeanFR_Rand_Order_Norm(ConmmonSceneName_Rand==UniqueScene(i))-SceneFR_mean(i))/SceneFR_std(i);
     %Normed PSTH
     PSTH_Rand_Order_Norm(ConmmonSceneName_Rand==UniqueScene(i),:)=(PSTH_Rand_Order_Norm(ConmmonSceneName_Rand==UniqueScene(i),:)-ScenePSTH_mean)./ScenePSTH_std
     
    
 end

 for j=1:NumberOfUniqueScene %Loop over the order
     OrderFR_mean(j)=nanmean(MeanFR_Rand_trial_fit(OrderRand_fit==j));
  
    OrderFR_std(j)=nanstd(MeanFR_Rand_trial_fit(OrderRand_fit==j));
    
     OrderPSTH_mean=nanmean(PSTH_Rand_trial(OrderRand_fit==j,:));
     OrderPSTH_std=nanstd(PSTH_Rand_trial(OrderRand_fit==j,:));
     
    
    MeanFR_Rand_Scene_Norm(OrderRand_fit==j)= (MeanFR_Rand_Scene_Norm(OrderRand_fit==j)-OrderFR_mean(j))/OrderFR_std(j);
  
      %Normed PSTH
     PSTH_Rand_Scene_Norm(OrderRand_fit==j,:)=(PSTH_Rand_Order_Norm(OrderRand_fit==j,:)-OrderPSTH_mean)./OrderPSTH_std;
     
    
 end

%Normed fr

for i=1:NumberOfUniqueScene %Loop over the scene
     MeanFR_Rand_Scene_Norm_mean(i)=nanmean(MeanFR_Rand_Scene_Norm(ConmmonSceneName_Rand==UniqueScene(i)));
     
     MeanFR_Rand_Scene_Norm_sem(i)=nanstd(MeanFR_Rand_Scene_Norm(ConmmonSceneName_Rand==UniqueScene(i)))/sqrt(sum(ConmmonSceneName_Rand==UniqueScene(i)));
     
     PSTH_Rand_Scene_Norm_mean(i,:)=nanmean(PSTH_Rand_Scene_Norm(ConmmonSceneName_Rand==UniqueScene(i),:));
  PSTH_Rand_Scene_Norm_sem(i,:)=nanstd(PSTH_Rand_Scene_Norm(ConmmonSceneName_Rand==UniqueScene(i),:))/sqrt(sum(ConmmonSceneName_Rand==UniqueScene(i)));
    
 end

 for j=1:NumberOfUniqueScene %Loop over the order
     MeanFR_Rand_Order_Norm_mean(j)=nanmean(MeanFR_Rand_Order_Norm(OrderRand_fit==j));
  
    MeanFR_Rand_Order_Norm_sem(j)=nanstd(MeanFR_Rand_Order_Norm(OrderRand_fit==j))/sqrt(sum(OrderRand_fit==j));
    
  
  PSTH_Rand_Order_Norm_mean(j,:)=nanmean(PSTH_Rand_Order_Norm(OrderRand_fit==j,:));
  PSTH_Rand_Order_Norm_sem(j,:)=nanstd(PSTH_Rand_Order_Norm(OrderRand_fit==j,:))/sqrt(sum(OrderRand_fit==j));
  
     
 end
 
 
 
 

%p_Scene_single=anova1(MeanFR_Rand_Scene_Norm,ConmmonSceneName_Rand,'display','off');
%p_Order_single=anova1(MeanFR_Rand_Scene_Norm,OrderRand_fit,'display','off');
%}
%{
%Average
FR_Scene_Rand_Whole_Order_mean=cell2mat(cellfun(@nanmean,FR_Scene_Rand_Whole_Order,'uniform',0));
FR_Scene_Fix_Whole_Order_mean=cell2mat(cellfun(@nanmean,FR_Scene_Fix_Whole_Order,'uniform',0));

FR_Scene_Rand_Whole_Order_sem=cell2mat(cellfun(@(x) nanstd(x)/sqrt(sum(~isnan(x))),FR_Scene_Rand_Whole_Order,'uniform',0));
FR_Scene_Fix_Whole_Order_sem=cell2mat(cellfun(@(x) nanstd(x)/sqrt(sum(~isnan(x))),FR_Scene_Fix_Whole_Order,'uniform',0));
%}






%{
%Repetition changes

NumMaxRand=max(cell2mat(cellfun(@numel,FR_Scene_Rand_Whole_Rep_Cell,'Uniform',0)));
%Change cell into a matrix
FR_Scene_Rand_Whole_Rep=cell2mat(cellfun(@(x) [x,ones(1,NumMaxRand-numel(x))]',FR_Scene_Rand_Whole_Rep_Cell,'Uniform',0));

NumMaxFix=max(cell2mat(cellfun(@numel,FR_Scene_Fix_Whole_Rep_Cell,'Uniform',0)));
%Change cell into a matrix
FR_Scene_Fix_Whole_Rep=cell2mat(cellfun(@(x) [x,ones(1,NumMaxFix-numel(x))]', FR_Scene_Fix_Whole_Rep_Cell,'Uniform',0));

%%Calculate mean and sem over scenes


FR_Scene_Rand_Whole_Rep_Mean=nanmean(FR_Scene_Rand_Whole_Rep,2);
FR_Scene_Fix_Whole_Rep_Rep_Mean=nanmean(FR_Scene_Fix_Whole_Rep,2);

FR_Scene_Rand_Whole_Rep_Sem=nanstd(FR_Scene_Rand_Whole_Rep,[],2)/sqrt(NumberOfUniqueScene);
FR_Scene_Fix_Whole_Rep_Sem=nanstd(FR_Scene_Fix_Whole_Rep,[],2)/sqrt(NumberOfUniqueScene);

%}




%%
%If show figures 
if ShowFigureFlag
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Bar Plot for the averagre firing rate between fixed and random sequence condition%%%%%%%%%%%%%%%%%%%%%%%%

titlestr{1}='Normed_FR_SceneOrder';
fig{1}=PrepareFigure(113,'w',[200,100, 800,600],'Name',titlestr{1});
subplot(2,1,1);
errorbar(MeanFR_Rand_Scene_Norm_mean,MeanFR_Rand_Scene_Norm_sem,'-sb','LineWidth',3);
xticks([1:NumberOfUniqueScene]);
xticklabels(string(UniqueScene));
box off
title('Normed scene tuning');
set(gca,'LineWidth',3,'FontSize',20,'FontWeight','Bold');
subplot(2,1,2);
errorbar(MeanFR_Rand_Order_Norm_mean,MeanFR_Rand_Order_Norm_sem,'-sb','LineWidth',3);
box off
title('Normed temporal order tuning');
set(gca,'LineWidth',3,'FontSize',20,'FontWeight','Bold');


%{
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%PSTH for scene and order
titlestr{2}='Normed_PSTH_Order';

fig{2}=PrepareFigure(114,'w',[200,100, 800,800],'Name',titlestr{2});

titlestr{3}='Normed_PSTH_Scene';

fig{3}=PrepareFigure(115,'w',[200,100, 800,800],'Name',titlestr{3});

for i=1:NumberOfUniqueScene
     set(groot,'CurrentFigure',fig{2});
     subplot(ceil(NumberOfUniqueScene/2),2,i);
    shadedErrorBar(PSTH_Time,PSTH_Rand_Order_Norm_mean(i,:),PSTH_Rand_Order_Norm_sem(i,:),'lineprops',{Color(i,:)},'transparent',1);   
   
   hold on 
  set(findall(gca, 'Type', 'Line'),'LineWidth',3);
  
   set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
   
    set(groot,'CurrentFigure',fig{3});
   subplot(ceil(NumberOfUniqueScene/2),2,i);
    shadedErrorBar(PSTH_Time,PSTH_Rand_Scene_Norm_mean(i,:),PSTH_Rand_Scene_Norm_sem(i,:),'lineprops',{Color(i,:)},'transparent',1);   
   
   hold on 
  set(findall(gca, 'Type', 'Line'),'LineWidth',3);
  
   set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
    
end

%}





%%
%}
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
    if isfield(OutputFileName,'SceneTempOrderSep')
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

OutputData_New.SceneTempOrderSep.TrialType='GoodTrials';

%%%%%%%%%Save overall scene response characteristics%%%%%%%%%%%%%%%%%%%%%%%
%Organize the data
%{
NameStr={'TrialNum','AverageFR_Fix_Rand','p_Rand_Fix_Whole','SceneSelectiveIndex_Fix','SceneSelectiveIndex_Rand',...
    'PSTH_SSI_Fix','PSTH_SSI_Rand','pSSI_fix','pSSI_rand','PSTH_Time','PSTH_Fix_All_Mean','PSTH_Fix_All_Sem',...
'PSTH_Rand_All_Mean','PSTH_Rand_All_Sem'

};





DataLib={TrialNum,[AverageFR_Fix_Rand_mean',AverageFR_Fix_Rand_sem'],[h_Rand_Fix_Whole',p_Rand_Fix_Whole'],SceneSelectiveIndex_Fix,SceneSelectiveIndex_Rand,...  
SSI_time_Fix,SSI_time_Rand,[h_SSI_fix',pSSI_fix'],[h_SSI_rand',pSSI_rand'],PSTH_Time,PSTH_Fix_All_Mean,PSTH_Fix_All_Sem,...
PSTH_Rand_All_Mean,PSTH_Rand_All_Sem
};


DataStamp=containers.Map(NameStr,DataLib);


OutputData_New.SceneTempOrderSep.Task=TaskType;
OutputData_New.SceneTempOrderSep.TaskCode=TaskCode;
OutputData_New.SceneTempOrderSep.DataStamp=DataStamp;
%}
%Store figures
if ShowFigureFlag
   
    %OutputData_New.SceneTuning.Figures={fig99,fig100,fig101,fig102};
    %Output figures in a folder;
    FolderName=OutputFigureName;
    if ~exist(FolderName,'dir')
        mkdir(FolderName);
    end
    cd(FolderName);
    
    for j=1:length(titlestr)
        OutputFigureName_curr=strcat(OutputFigureName,strcat(titlestr{j},'.jpg'));
        saveas(fig{i},OutputFigureName_curr);
    end
    
    disp('Figures have been exported to the neuron folder');

end







if ExistFlag
%Compare the old one with the new one
Task_Old=OutputData.SceneTempOrderSep.TaskCode;

     if Task_Old~=TaskCode
        %Not the same task, add the new dataset into the old one
        OutputData.SceneTempOrderSep(numel(OutputData.SceneTempOrderSep)+1)=OutputData_New.SceneTempOrderSep;
     else
         %If the same task,the same protocol,replace the old one with the new one
         
         OutputData.SceneTempOrderSep=OutputData_New.SceneTempOrderSep;
         
    
    
     end
else
    %Output the current file
   OutputData.SceneTempOrderSep=OutputData_New.SceneTempOrderSep; 
end



if OutputFlag
    cd(OutPath);
    save(OutputFileName,'OutputData');
    disp('Data Exported');
else
    disp('No data export');
end  




end