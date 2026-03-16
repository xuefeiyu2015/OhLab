function  Batch_EventRT_Exploration(Data,  StartFile, EndFile,ShowFigureFlag,OutputFlag);
%Summary for the event segregation analysis 
%Xuefei Yu 01192021
BatchFileName=Data.BatchFileName;
%Batch data file path
FilesPath=Data.ResultFilePath(StartFile: EndFile);
FilesName=Data.FileName;
RecordDate=Data.RecordDate(StartFile: EndFile);
FileNum=EndFile-StartFile+1;
TotalNeuronNum=sum(cellfun(@numel,FilesPath));
ChannelID=Data.ChannelNumber;

PlotExample=1;
%Color Definition for each events
Color=[lines(6);prism(6);pink(6)];
%Color(2,:)=copper(1)+0.5;
Color=[Color;hsv(10)];
Color(18,:)=[0.1,0.457,0.194];

Blockstr={'FixLoc,FixSceneSeq','FixLoc,RandSceneSeq','RandLoc,FixSceneSeq','RandLoc,RandSceneSeq'};
ColorForEachMode=[0.83 0.14 0.14
             1.00 0.54 0.00
             0.47 0.25 0.80
             0.25 0.80 0.54];
         
StableSceneColor=[ 1.00 0.54 0.00];
FlexibleSceneColor=[ 0.25 0.80 0.54];
DifferenceColor=[0.4471    0.3020    0.9176];
PassiveIndex=NaN;
ChannelFull=[];
%Load file
FileIndex=1;
for i=1:length(FilesPath)
  clear OutputData;
  ChannelFull=[ChannelFull,ChannelID(i,1):ChannelID(i,2)];
  
  for j=1:length(FilesPath{i})
      currFileName=FilesPath{i}{j};
      load(currFileName);
     OutputDataTemp=OutputData; 
     
     Task(FileIndex)=OutputData.EventSegRT(1).TaskCode;
     TaskName{FileIndex}=OutputData.EventSegRT(1).Task;  
     FilesNameEach{FileIndex}=FilesName{i};
     
    
     %Separate Keyvalues in each container;
    % Stamp_tmp=OutputData.EventSegRT(1).DataStamp;
    % KeyWhole=keys(Stamp_tmp);
     
     EventSegRT{FileIndex}=OutputData.EventSegRT(1).DataStamp;
     
     %For active passive comparison
     if isfield(OutputData,'SceneTuningPassive')
         PassiveIndex(FileIndex)=1;
          SceneTuningPassive{FileIndex}=OutputData.SceneTuningPassive(1).DataStamp;
     else
         PassiveIndex(FileIndex)=0;
     
     end
     
          
      FileIndex=FileIndex+1;
      
  end
  
    
     
     
end
%Assume that keys are all the same for each file, use the first file to
%abstract keys for all files;
KeysLib=keys(EventSegRT{1});
AllValues=cellfun(@(x) values(x),EventSegRT,'uniform',0);
%Reorganize the value according to each key
for i=1:length(KeysLib)
  %  currKey=KeysLib{i};
   KeyValues{i}=cellfun(@(x) x{i},AllValues,'uniform',0);
 
end

%Load scene passive files
KeysLibPassive=keys(SceneTuningPassive{find(PassiveIndex==1,1,'first')});
AllValuesPassive=cellfun(@(x) values(x),SceneTuningPassive(PassiveIndex==1),'uniform',0);
%Reorganize the value according to each key
for i=1:length(KeysLibPassive)
  %  currKey=KeysLib{i};
   KeyValuesPassive{i}=cellfun(@(x) x{i},AllValuesPassive,'uniform',0);
 
end


%Form the table for the neuron basic
FocusKey={'AverageFR','AverageSceneDiffBlock','AverageSpon','p_SceneDiffBetweenBlocks'};

 Sel=ismember(KeysLib,FocusKey);
 FocusKeyData=cell2mat(cellfun(@(x) cell2mat(x)',KeyValues(Sel),'uniform',0));
 KeySel=KeysLib(Sel);
 TableTitle={'Average','AverageSpon','FlexibleSceneAverage','StableSceneAverage','p_SceneDiffBetweenBlocks'};
 TableData=FocusKeyData(:,[1,4,2,3,5]);
 TableOutput=array2table(TableData,'VariableNames',TableTitle);
 
 AverageFR=TableData(:,1);
 
 %Separate neuron into three different groups according to the difference
 %between two blocks;
 p_SceneBlockDiff=TableData(:,end);
 SceneResponseDiff_Overall=(TableData(:,3)-TableData(:,4))./TableData(:,2);%Flexible-Stable
 FlexDominace=(SceneResponseDiff_Overall>0) & (p_SceneBlockDiff<0.05);
 StableDominance=(SceneResponseDiff_Overall<0) & (p_SceneBlockDiff<0.05);
 NoDominance=~(FlexDominace |StableDominance);
 
 
 
 
 
 
 %About the PSTH for the scene responses of different block
 PSTHSceneKey={'PSTH_SceneOff_FleSta_TimeMeanSem','PSTH_SceneOn_FleSta_TimeMeanSem'};
 Sel_tmp=ismember(KeysLib,PSTHSceneKey);
 PSTHSceneKey_Sel=KeysLib(Sel_tmp);
 PSTHSceneData_Mean=cellfun(@(x) cellfun(@(y) y(:,[2:3]),x,'uniform',0),KeyValues(Sel_tmp),'uniform',0);%Normalize
 PSTHSceneData_Time=cellfun(@(x) cellfun(@(y) y(:,1),x,'uniform',0),KeyValues(Sel_tmp),'uniform',0);
 
 %Denomiter for normalization
 SDKey='SD_RelSceneOn';
 Sel_tmp=ismember(KeysLib,SDKey);
 SDKey_SD=cell2mat(KeyValues{Sel_tmp})';
 
 
 Norm_Denormiter=SDKey_SD;
 
MeanKey='Mean_RelSceneOn';
Sel_tmp=ismember(KeysLib,MeanKey);
MeanKey_Mean=cell2mat(KeyValues{Sel_tmp})';
 
 
 
 %%%%%%%%%%First for scene on responses
 SelOn=ismember(PSTHSceneKey_Sel,'PSTH_SceneOn_FleSta_TimeMeanSem');
 PSTHSceneOnFlexible_Mean=cellfun(@(x) x(:,1),PSTHSceneData_Mean{SelOn},'uniform',0);
 PSTHSceneOnStable_Mean=cellfun(@(x) x(:,2),PSTHSceneData_Mean{SelOn},'uniform',0);
 %Put time on the x-axis for more ituition;
 MaxRow=max(cellfun(@(x) size(x,1),PSTHSceneOnFlexible_Mean));
 PSTHSceneOnFlexible_MeanMat=(cell2mat(cellfun(@(x) [x; NaN*ones(MaxRow-size(x,1),size(x,2))],PSTHSceneOnFlexible_Mean,'uniform',0))'-MeanKey_Mean)./repmat(Norm_Denormiter,1,MaxRow);%Normalized by SD
 PSTHSceneOnStable_MeanMat=(cell2mat(cellfun(@(x) [x; NaN*ones(MaxRow-size(x,1),size(x,2))],PSTHSceneOnStable_Mean,'uniform',0))'-MeanKey_Mean)./repmat(Norm_Denormiter,1,MaxRow);%Normalized by SD

 %Time for psth
 PSTHSceneOn_Time=PSTHSceneData_Time{SelOn};
 %Put time on the x-axis for more ituition;
 MaxRow=max(cellfun(@(x) size(x,1),PSTHSceneOn_Time));
 PSTHSceneOn_TimeMat=cell2mat(cellfun(@(x) [x; NaN*ones(MaxRow-size(x,1),size(x,2))],PSTHSceneOn_Time,'uniform',0))';
 PSTHSceneOn_TimeMean=nanmean(PSTHSceneOn_TimeMat,1);
 
 %%%%%%%%%%Then for scene off responses;
SelOff=ismember(PSTHSceneKey_Sel,'PSTH_SceneOff_FleSta_TimeMeanSem');
 PSTHSceneOffFlexible_Mean=cellfun(@(x) x(:,1),PSTHSceneData_Mean{SelOff},'uniform',0);
 PSTHSceneOffStable_Mean=cellfun(@(x) x(:,2),PSTHSceneData_Mean{SelOff},'uniform',0);
 %Put time on the x-axis for more ituition;
 MaxRow=max(cellfun(@(x) size(x,1),PSTHSceneOffFlexible_Mean));
 PSTHSceneOffFlexible_MeanMat=(cell2mat(cellfun(@(x) [x; NaN*ones(MaxRow-size(x,1),size(x,2))],PSTHSceneOffFlexible_Mean,'uniform',0))'-MeanKey_Mean)./repmat(Norm_Denormiter,1,MaxRow);%Normalized by SD
 PSTHSceneOffStable_MeanMat=(cell2mat(cellfun(@(x) [x; NaN*ones(MaxRow-size(x,1),size(x,2))],PSTHSceneOffStable_Mean,'uniform',0))'-MeanKey_Mean)./repmat(Norm_Denormiter,1,MaxRow);%Normalized by SD

 %Time for psth
 PSTHSceneOff_Time=PSTHSceneData_Time{SelOff};
 %Put time on the x-axis for more ituition;
 MaxRow=max(cellfun(@(x) size(x,1),PSTHSceneOff_Time));
 PSTHSceneOff_TimeMat=cell2mat(cellfun(@(x) [x; NaN*ones(MaxRow-size(x,1),size(x,2))],PSTHSceneOff_Time,'uniform',0))';
 PSTHSceneOff_TimeMean=nanmean(PSTHSceneOff_TimeMat,1);
 
 
 %Select out the scene responsive neurons
 CriteriumKey={'SceneOn500Mean','p_SceneOn500withSpon','SceneOnMeanMinusSponWhole','p_SceneOnWithSponWhole',...
  'pSceneOnWhole','SceneOnDifference_WholeMean','pSceneOn500','SceneOnDifference_500Mean'};

%SceneOnDiff500Mean=cell2mat(KeyValues{ismember(KeysLib,'SceneOnDifference_500Mean')});
 Sel=ismember(KeysLib,CriteriumKey);
 FocusKeyData=KeyValues(Sel);
 FocusKey_Sel=KeysLib(Sel);
 %First use the 500 criterium
 SceneOnDifference_500Mean=cell2mat(FocusKeyData{ismember(FocusKey_Sel,'SceneOnDifference_500Mean')});
 p_SceneOnDifference_500=cell2mat(FocusKeyData{ismember(FocusKey_Sel,'pSceneOn500')});
 
 %Select out the scene respinsive neuron
 SceneResponsiveFlag=(p_SceneOnDifference_500<0.05)';
 NumOfSceneResponsiveCell=sum(SceneResponsiveFlag);
 
 %Another criterium: any one of the scene has responses
 p_SceneOn500withSpon=cell2mat(cellfun(@(x) x',FocusKeyData{ismember(FocusKey_Sel,'p_SceneOn500withSpon')},'uniform',0));
 
 NumOfSignificantScene500=sum(p_SceneOn500withSpon<0.05);
 
 SceneResponsiveFlag500=(NumOfSignificantScene500>0)';
 
  NumOfSceneResponsiveCell500=sum(SceneResponsiveFlag500);
  
  %Try the whole scene response version
  p_SceneOnwithSpon=cell2mat(cellfun(@(x) x',FocusKeyData{ismember(FocusKey_Sel,'p_SceneOnWithSponWhole')},'uniform',0));
 
 NumOfSignificantScene=sum(p_SceneOnwithSpon<0.05);
 
 SceneResponsiveFlag=(NumOfSignificantScene>0)';
 
  NumOfSceneResponsiveCellWhole=sum(SceneResponsiveFlag);
  
%% 
%%eparate scene-on/off responses for the most responsive scene for each group;
   CriteriumKey2={ 'NumOfStableSigScenes','NumOfFlexibleSigScenes','StableSigScenes','FlexibleSigScenes'...
  'PSTH_MostPreferedStableScene_MeanSem','PSTH_MostPreferedFlexibleScene_MeanSem',...
  'PSTH_OtherStableScene_MeanSem','PSTH_OtherFlexibleScene_MeanSem',};

%SceneOnDiff500Mean=cell2mat(KeyValues{ismember(KeysLib,'SceneOnDifference_500Mean')});
 Sel=ismember(KeysLib,CriteriumKey2);
 FocusKeyData=KeyValues(Sel);
 FocusKey_Sel=KeysLib(Sel);
  
 NumOfStableSigScenes=cell2mat(cellfun(@(x) x',FocusKeyData{ismember(FocusKey_Sel,'NumOfStableSigScenes')},'uniform',0));
 NumOfFlexibleSigScenes=cell2mat(cellfun(@(x) x',FocusKeyData{ismember(FocusKey_Sel,'NumOfFlexibleSigScenes')},'uniform',0));
 
PropOfNeuronOnlyStable=sum(NumOfStableSigScenes>0 & NumOfFlexibleSigScenes==0)/length(NumOfStableSigScenes);
PropOfNeuronOnlyFlexible=sum(NumOfFlexibleSigScenes>0 & NumOfStableSigScenes==0)/length(NumOfStableSigScenes);
PropOfNeuronBoth=sum(NumOfFlexibleSigScenes>0 & NumOfStableSigScenes>0)/length(NumOfStableSigScenes);
PropOfNeuronNone=sum(NumOfFlexibleSigScenes==0 & NumOfStableSigScenes==0)/length(NumOfStableSigScenes);
 
 
PSTH_MostPreferedStableScene_MeanSem=cellfun(@(x) x,FocusKeyData{ismember(FocusKey_Sel,'PSTH_MostPreferedStableScene_MeanSem')},'uniform',0);
PSTH_MostPreferedStableScene_Mean=cellfun(@(x) x(:,1),PSTH_MostPreferedStableScene_MeanSem,'uniform',0);

PSTH_MostPreferedFlexibleScene_MeanSem=cellfun(@(x) x,FocusKeyData{ismember(FocusKey_Sel,'PSTH_MostPreferedFlexibleScene_MeanSem')},'uniform',0);
PSTH_MostPreferedFlexibleScene_Mean=cellfun(@(x) x(:,1),PSTH_MostPreferedFlexibleScene_MeanSem,'uniform',0);

PSTH_OtherStableScene_MeanSem=cellfun(@(x) x,FocusKeyData{ismember(FocusKey_Sel,'PSTH_OtherStableScene_MeanSem')},'uniform',0);
PSTH_OtherStableScene_Mean=cellfun(@(x) x(:,1),PSTH_OtherStableScene_MeanSem,'uniform',0);

PSTH_OtherFlexibleScene_MeanSem=cellfun(@(x) x,FocusKeyData{ismember(FocusKey_Sel,'PSTH_OtherFlexibleScene_MeanSem')},'uniform',0);
PSTH_OtherFlexibleScene_Mean=cellfun(@(x) x(:,1),PSTH_OtherFlexibleScene_MeanSem,'uniform',0);
  
%AverageTimeBin
TimeBin=diff(PSTHSceneOn_TimeMean);
TimeBin=TimeBin(1);
TimeAverageDur=500;
AverageTime=PSTHSceneOn_TimeMean>0 & PSTHSceneOn_TimeMean<TimeAverageDur;

AveMostPreferResponseFromSpon_Stable_abs=abs(cellfun(@(x) nanmean(x(AverageTime),1),PSTH_MostPreferedStableScene_Mean));
AveMostPreferResponseFromSpon_Flexible_abs=abs(cellfun(@(x) nanmean(x(AverageTime),1),PSTH_MostPreferedFlexibleScene_Mean));


AveOtherResponseFromSpon_Stable_abs=abs(cellfun(@(x) nanmean(x(AverageTime),1),PSTH_OtherStableScene_Mean));
AveOtherResponseFromSpon_Flexible_abs=abs(cellfun(@(x) nanmean(x(AverageTime),1),PSTH_OtherFlexibleScene_Mean));

AveMostPreferResponseFromSpon_Stable_raw=cellfun(@(x) nanmean(x(AverageTime),1),PSTH_MostPreferedStableScene_Mean);
AveMostPreferResponseFromSpon_Flexible_raw=cellfun(@(x) nanmean(x(AverageTime),1),PSTH_MostPreferedFlexibleScene_Mean);

%Zscore
AveMostPreferResponseFromSpon_Stable_ZScore=(AveMostPreferResponseFromSpon_Stable_raw'-MeanKey_Mean)./Norm_Denormiter;

AveMostPreferResponseFromSpon_Flexible_ZScore=(AveMostPreferResponseFromSpon_Flexible_raw'-MeanKey_Mean)./Norm_Denormiter;


AveOtherResponseFromSpon_Stable_raw=cellfun(@(x) nanmean(x(AverageTime),1),PSTH_OtherStableScene_Mean);
AveOtherResponseFromSpon_Flexible_raw=cellfun(@(x) nanmean(x(AverageTime),1),PSTH_OtherFlexibleScene_Mean);

AveOtherResponseFromSpon_Stable_ZScore=(AveOtherResponseFromSpon_Stable_raw'-MeanKey_Mean)./Norm_Denormiter;

AveOtherResponseFromSpon_Flexible_ZScore=(AveOtherResponseFromSpon_Flexible_raw'-MeanKey_Mean)./Norm_Denormiter;


%-MeanKey_Mean)./repmat(Norm_Denormiter,1,MaxRow);
 
 %Separate scene-on/off responses as a whole for each group;
 
 FlexDominace=(SceneResponseDiff_Overall>0) & (p_SceneBlockDiff<0.05) & SceneResponsiveFlag500;
 StableDominance=(SceneResponseDiff_Overall<0) & (p_SceneBlockDiff<0.05)& SceneResponsiveFlag500;
 NoDominance=~(FlexDominace |StableDominance)&SceneResponsiveFlag500;
 
 %Flexible Dominance
 %Scene on
 PSTHSceneOnFlexible_FlexDo_Mean=nanmean(PSTHSceneOnFlexible_MeanMat(FlexDominace,:),1);
 PSTHSceneOnFlexible_FlexDo_Sem=nanstd(PSTHSceneOnFlexible_MeanMat(FlexDominace,:),[],1)./sqrt(sum(~isnan(PSTHSceneOnFlexible_MeanMat(FlexDominace,:)),1));
 
 PSTHSceneOnStable_FlexDo_Mean=nanmean(PSTHSceneOnStable_MeanMat(FlexDominace,:),1);
 PSTHSceneOnStable_FlexDo_Sem=nanstd(PSTHSceneOnStable_MeanMat(FlexDominace,:),[],1)./sqrt(sum(~isnan(PSTHSceneOnFlexible_MeanMat(FlexDominace,:)),1));
 
 %Scene off
PSTHSceneOffFlexible_FlexDo_Mean=nanmean(PSTHSceneOffFlexible_MeanMat(FlexDominace,:),1);
 PSTHSceneOffFlexible_FlexDo_Sem=nanstd(PSTHSceneOffFlexible_MeanMat(FlexDominace,:),[],1)./sqrt(sum(~isnan(PSTHSceneOffFlexible_MeanMat(FlexDominace,:)),1));
 
 PSTHSceneOffStable_FlexDo_Mean=nanmean(PSTHSceneOffStable_MeanMat(FlexDominace,:),1);
 PSTHSceneOffStable_FlexDo_Sem=nanstd(PSTHSceneOffStable_MeanMat(FlexDominace,:),[],1)./sqrt(sum(~isnan(PSTHSceneOffFlexible_MeanMat(FlexDominace,:)),1));

 %Stable Dominance
 %Scene on
 PSTHSceneOnFlexible_StaDo_Mean=nanmean(PSTHSceneOnFlexible_MeanMat(StableDominance,:),1);
 PSTHSceneOnFlexible_StaDo_Sem=nanstd(PSTHSceneOnFlexible_MeanMat(StableDominance,:),[],1)./sqrt(sum(~isnan(PSTHSceneOnFlexible_MeanMat(StableDominance,:)),1));
 PSTHSceneOnStable_StaDo_Mean=nanmean(PSTHSceneOnStable_MeanMat(StableDominance,:),1);
 PSTHSceneOnStable_StaDo_Sem=nanstd(PSTHSceneOnStable_MeanMat(StableDominance,:),[],1)./sqrt(sum(~isnan(PSTHSceneOnFlexible_MeanMat(StableDominance,:)),1));
 
 %Scene off
 PSTHSceneOffFlexible_StaDo_Mean=nanmean(PSTHSceneOffFlexible_MeanMat(StableDominance,:),1);
 PSTHSceneOffFlexible_StaDo_Sem=nanstd(PSTHSceneOffFlexible_MeanMat(StableDominance,:),[],1)./sqrt(sum(~isnan(PSTHSceneOffFlexible_MeanMat(StableDominance,:)),1));
 PSTHSceneOffStable_StaDo_Mean=nanmean(PSTHSceneOffStable_MeanMat(StableDominance,:),1);
 PSTHSceneOffStable_StaDo_Sem=nanstd(PSTHSceneOffStable_MeanMat(StableDominance,:),[],1)./sqrt(sum(~isnan(PSTHSceneOffFlexible_MeanMat(StableDominance,:)),1));
 
 
 %No Dominance
 
 PSTHSceneOnFlexible_NoDo_Mean=nanmean(PSTHSceneOnFlexible_MeanMat(NoDominance,:),1);
 PSTHSceneOnFlexible_NoDo_Sem=nanstd(PSTHSceneOnFlexible_MeanMat(NoDominance,:),[],1)./sqrt(sum(~isnan(PSTHSceneOnFlexible_MeanMat(NoDominance,:)),1));
 PSTHSceneOnStable_NoDo_Mean=nanmean(PSTHSceneOnStable_MeanMat(NoDominance,:),1);
 PSTHSceneOnStable_NoDo_Sem=nanstd(PSTHSceneOnStable_MeanMat(NoDominance,:),[],1)./sqrt(sum(~isnan(PSTHSceneOnFlexible_MeanMat(NoDominance,:)),1));
 
 PSTHSceneOffFlexible_NoDo_Mean=nanmean(PSTHSceneOffFlexible_MeanMat(NoDominance,:),1);
 PSTHSceneOffFlexible_NoDo_Sem=nanstd(PSTHSceneOffFlexible_MeanMat(NoDominance,:),[],1)./sqrt(sum(~isnan(PSTHSceneOffFlexible_MeanMat(NoDominance,:)),1));
 PSTHSceneOffStable_NoDo_Mean=nanmean(PSTHSceneOffStable_MeanMat(NoDominance,:),1);
 PSTHSceneOffStable_NoDo_Sem=nanstd(PSTHSceneOffStable_MeanMat(NoDominance,:),[],1)./sqrt(sum(~isnan(PSTHSceneOffFlexible_MeanMat(NoDominance,:)),1));
 
 
 %Scene Selectivity Index Distribution
 
 SIKey={'SI_FleSta'};
 Sel_tmp=ismember(KeysLib,SIKey);
 SISceneKey_Sel=KeysLib(Sel_tmp);
 SISceneOnData_Mean=cell2mat(cellfun(@(x) x',KeyValues{Sel_tmp},'uniform',0));
 SISceneOn_Fle=SISceneOnData_Mean(1,:);
 SISceneOn_Sta=SISceneOnData_Mean(2,:);
 
 %Significance quantified by anova
 sigKey={'p_AnovaSceneOnAllFleSta'};
 Sel_tmp=ismember(KeysLib,sigKey);
 sig_Anova_SceneKey_Sel=KeysLib(Sel_tmp);
 sig_Anova_SceneOnData=cell2mat(cellfun(@(x) x',KeyValues{Sel_tmp},'uniform',0));
 
 sigSceneOn_Who=sig_Anova_SceneOnData(1,:);
 sigSceneOn_Fle=sig_Anova_SceneOnData(2,:);
 sigSceneOn_Sta=sig_Anova_SceneOnData(3,:);
 
 NumOfSig_Fle=sum(sigSceneOn_Fle<0.05);
 NumOfSig_Sta=sum(sigSceneOn_Sta<0.05);
 
 %Separate with different scene response type
 SI_SceneOn_Fle_FleDo=SISceneOn_Fle(FlexDominace);
 SI_SceneOn_Sta_FleDo=SISceneOn_Sta(FlexDominace);
  
 
 sigSceneOn_Fle_FleDo=sigSceneOn_Fle(FlexDominace);
 sigSceneOn_Sta_FleDo=sigSceneOn_Fle(FlexDominace);
 
 
 SI_SceneOn_Fle_StaDo=SISceneOn_Fle(StableDominance);
 SI_SceneOn_Sta_StaDo=SISceneOn_Sta(StableDominance);
 
 sigSceneOn_Fle_StaDo=sigSceneOn_Fle(StableDominance);
 sigSceneOn_Sta_StaDo=sigSceneOn_Fle(StableDominance);
 
 
 SI_SceneOn_Fle_NoDo=SISceneOn_Fle(NoDominance);
 SI_SceneOn_Sta_NoDo=SISceneOn_Sta(NoDominance);
 
 sigSceneOn_Fle_NoDo=sigSceneOn_Fle(NoDominance);
 sigSceneOn_Sta_NoDo=sigSceneOn_Fle(NoDominance);
 
 
 %R square
  R2CriteriumKey={ 'R2WholeFleSta','p_SceneOnAnova_WholeStaFle','pSceneDiffWithTime','rsquare_psth',...
  'pSceneDiffWithTime_Stable','rsquare_psth_Stable','pSceneDiffWithTime_Flexible','rsquare_psth_Flexible'
    };

%SceneOnDiff500Mean=cell2mat(KeyValues{ismember(KeysLib,'SceneOnDifference_500Mean')});
 Sel=ismember(KeysLib,R2CriteriumKey);
 FocusKeyDataR2=KeyValues(Sel);
 FocusKey_Sel=KeysLib(Sel);
 %R2 significance
 
 p_SSI_WholeStaFle=FocusKeyDataR2{ismember(FocusKey_Sel,'p_SceneOnAnova_WholeStaFle')};
 p_SSI_Whole=cellfun(@(x) x(1),p_SSI_WholeStaFle);
 p_SSI_Sta=cellfun(@(x) x(2),p_SSI_WholeStaFle);
 p_SSI_Fle=cellfun(@(x) x(3),p_SSI_WholeStaFle);
 
 %R2
 R2SSI_WholeFleSta=FocusKeyDataR2{ismember(FocusKey_Sel,'R2WholeFleSta')};
 R2SSI_Whole=cellfun(@(x) x(1),R2SSI_WholeFleSta);
 R2SSI_Fle=cellfun(@(x) x(2),R2SSI_WholeFleSta);
 R2SSI_Sta=cellfun(@(x) x(3),R2SSI_WholeFleSta);
 
 %PSTH of R2
 %Whole
 rsquare_psth=FocusKeyDataR2{ismember(FocusKey_Sel,'rsquare_psth')};

 MaxCol=max(cellfun(@(x) size(x,2),rsquare_psth));
 rsquare_psthMat=cell2mat(cellfun(@(x) [x, NaN*ones(size(x,1),MaxCol-size(x,2))]',rsquare_psth,'uniform',0))';
 
 prsquare_psth=FocusKeyDataR2{ismember(FocusKey_Sel,'pSceneDiffWithTime')};

 MaxCol=max(cellfun(@(x) size(x,2), prsquare_psth));
 prsquare_psthMat=cell2mat(cellfun(@(x) [x, NaN*ones(size(x,1),MaxCol-size(x,2))]', prsquare_psth,'uniform',0))';
 
 
 %Stable

 rsquare_psth_stable=FocusKeyDataR2{ismember(FocusKey_Sel,'rsquare_psth_Stable')};
 MaxCol=max(cellfun(@(x) size(x,2),rsquare_psth_stable));
 rsquare_psth_stable_Mat=cell2mat(cellfun(@(x) [x, NaN*ones(size(x,1),MaxCol-size(x,2))]',rsquare_psth_stable,'uniform',0))';
 
 prsquare_psth_stable=FocusKeyDataR2{ismember(FocusKey_Sel,'pSceneDiffWithTime_Stable')};

 MaxCol=max(cellfun(@(x) size(x,2), prsquare_psth_stable));
 prsquare_psth_stableMat=cell2mat(cellfun(@(x) [x, NaN*ones(size(x,1),MaxCol-size(x,2))]', prsquare_psth_stable,'uniform',0))';
 
 sigrsquare_psth_stableMat=prsquare_psth_stableMat<0.05;
 
 PropOfSig_psth_stableMat=sum(sigrsquare_psth_stableMat,1)/size(sigrsquare_psth_stableMat,1)*100;
 
 
 %Flexible
 rsquare_psth_flexible=FocusKeyDataR2{ismember(FocusKey_Sel,'rsquare_psth_Flexible')};
 MaxCol=max(cellfun(@(x) size(x,2),rsquare_psth_flexible));
 rsquare_psth_flexible_Mat=cell2mat(cellfun(@(x) [x, NaN*ones(size(x,1),MaxCol-size(x,2))]',rsquare_psth_flexible,'uniform',0))';
 
 prsquare_psth_flexible=FocusKeyDataR2{ismember(FocusKey_Sel,'pSceneDiffWithTime_Flexible')};

 MaxCol=max(cellfun(@(x) size(x,2), prsquare_psth_flexible));
 prsquare_psth_flexibleMat=cell2mat(cellfun(@(x) [x, NaN*ones(size(x,1),MaxCol-size(x,2))]', prsquare_psth_flexible,'uniform',0))';
 
 sigrsquare_psth_flexibleMat=prsquare_psth_flexibleMat<0.05;
 PropOfSig_psth_flexibleMat=sum(sigrsquare_psth_flexibleMat,1)/size(sigrsquare_psth_flexibleMat,1)*100;
 
 %{
 
 figure
 plot(PropOfSig_psth_stableMat,'-r');
 hold on
 plot(PropOfSig_psth_flexibleMat,'-b');
 %}
 
 
 %Separate into different types
 SSI_Sig_Whole=p_SSI_Whole<0.05;
 
 rsquare_psthMat_WholeSig=rsquare_psthMat(SSI_Sig_Whole,:);
 rsquare_psthMat_WholeNonSig=rsquare_psthMat(~SSI_Sig_Whole,:);
 
 rsquare_psthMat_StableSig=rsquare_psth_stable_Mat(SSI_Sig_Whole,:);
 
 rsquare_psthMat_FlexibleSig=rsquare_psth_flexible_Mat(SSI_Sig_Whole,:);
 
 rsquare_psthMat_WholeSig_Mean=nanmean(rsquare_psthMat_WholeSig,1);
 rsquare_psthMat_WholeSig_Sem=nanstd(rsquare_psthMat_WholeSig,[],1)./sqrt(sum(~isnan(rsquare_psthMat_WholeSig),1));
 
 rsquare_psthMat_WholeNonSig_Mean=nanmean(rsquare_psthMat_WholeNonSig,1);
 rsquare_psthMat_WholeNonSig_Sem=nanstd(rsquare_psthMat_WholeNonSig,[],1)./sqrt(sum(~isnan(rsquare_psthMat_WholeNonSig),1));
 
 
 
 rsquare_psthMat_WholeSig_StableMean=nanmean(rsquare_psthMat_StableSig,1);
 rsquare_psthMat_WholeSig_StableSem=nanstd(rsquare_psthMat_StableSig,[],1)./sqrt(sum(~isnan(rsquare_psthMat_StableSig),1));
 
 rsquare_psthMat_WholeSig_FlexibleMean=nanmean(rsquare_psthMat_FlexibleSig,1);
 rsquare_psthMat_WholeSig_FlexibleSem=nanstd(rsquare_psthMat_FlexibleSig,[],1)./sqrt(sum(~isnan(rsquare_psthMat_FlexibleSig),1));
 
 %Stable Selective Neuron
 SSI_Sig_StableOnly=p_SSI_Sta<0.05 & p_SSI_Fle>=0.05 ;
 SSI_Sig_FlexibleOnly=p_SSI_Sta>=0.05 & p_SSI_Fle<0.05;
 SSI_Sig_Both=p_SSI_Sta<0.05 & p_SSI_Fle<0.05;
 SSI_Sig_None=p_SSI_Sta>=0.05 & p_SSI_Fle>=0.05;
 
 rsquare_psthMat_StableSig_Stable=rsquare_psth_stable_Mat(SSI_Sig_StableOnly,:);
 
 rsquare_psthMat_StableSig_Flexible=rsquare_psth_flexible_Mat(SSI_Sig_StableOnly,:);
 
 rsquare_psthMat_FlexibleSig_Stable=rsquare_psth_stable_Mat(SSI_Sig_FlexibleOnly,:);
 
 rsquare_psthMat_FlexibleSig_Flexible=rsquare_psth_flexible_Mat(SSI_Sig_FlexibleOnly,:);
 
 rsquare_psthMat_BothSig_Stable=rsquare_psth_stable_Mat(SSI_Sig_Both,:);
 
 rsquare_psthMat_BothSig_Flexible=rsquare_psth_flexible_Mat(SSI_Sig_Both,:);
 
 rsquare_psthMat_StableSig_Stable_Mean=nanmean(rsquare_psthMat_StableSig_Stable,1);
 rsquare_psthMat_StableSig_Stable_Sem=nanstd(rsquare_psthMat_StableSig_Stable,[],1)./sqrt(sum(~isnan(rsquare_psthMat_StableSig_Stable),1));
 
 
 rsquare_psthMat_StableSig_Flexible_Mean=nanmean(rsquare_psthMat_StableSig_Flexible,1);
 rsquare_psthMat_StableSig_Flexible_Sem=nanstd(rsquare_psthMat_StableSig_Flexible,[],1)./sqrt(sum(~isnan(rsquare_psthMat_StableSig_Flexible),1));
 
 diff_rsquare_psthMat_StableSig=rsquare_psthMat_StableSig_Stable-rsquare_psthMat_StableSig_Flexible;
 diff_rsquare_psthMat_StableSig_Mean=nanmean(diff_rsquare_psthMat_StableSig,1);
 diff_rsquare_psthMat_StableSig_Sem=nanstd(diff_rsquare_psthMat_StableSig,[],1)./sqrt(sum(~isnan(diff_rsquare_psthMat_StableSig),1));  
   
   
 
 rsquare_psthMat_FlexibleSig_Stable_Mean=nanmean(rsquare_psthMat_FlexibleSig_Stable,1);
 rsquare_psthMat_FlexibleSig_Stable_Sem=nanstd(rsquare_psthMat_FlexibleSig_Stable,[],1)./sqrt(sum(~isnan(rsquare_psthMat_FlexibleSig_Stable),1));
 
 rsquare_psthMat_FlexibleSig_Flexible_Mean=nanmean(rsquare_psthMat_FlexibleSig_Flexible,1);
 rsquare_psthMat_FlexibleSig_Flexible_Sem=nanstd(rsquare_psthMat_FlexibleSig_Flexible,[],1)./sqrt(sum(~isnan(rsquare_psthMat_FlexibleSig_Flexible),1));
 
 
 diff_rsquare_psthMat_FlexibleSig=rsquare_psthMat_FlexibleSig_Stable-rsquare_psthMat_FlexibleSig_Flexible;
 diff_rsquare_psthMat_FlexibleSig_Mean=nanmean(diff_rsquare_psthMat_FlexibleSig,1);
 diff_rsquare_psthMat_FlexibleSig_Sem=nanstd(diff_rsquare_psthMat_FlexibleSig,[],1)./sqrt(sum(~isnan(diff_rsquare_psthMat_FlexibleSig),1));  
 
 
 rsquare_psthMat_BothSig_Stable_Mean=nanmean(rsquare_psthMat_BothSig_Stable,1);
 rsquare_psthMat_BothSig_Stable_Sem=nanstd(rsquare_psthMat_BothSig_Stable,[],1)./sqrt(sum(~isnan(rsquare_psthMat_BothSig_Stable),1));
 
 rsquare_psthMat_BothSig_Flexible_Mean=nanmean(rsquare_psthMat_BothSig_Flexible,1);
 rsquare_psthMat_BothSig_Flexible_Sem=nanstd(rsquare_psthMat_BothSig_Flexible,[],1)./sqrt(sum(~isnan(rsquare_psthMat_BothSig_Flexible),1));
 
 diff_rsquare_psthMat_BothSig=rsquare_psthMat_BothSig_Stable-rsquare_psthMat_BothSig_Flexible;
 diff_rsquare_psthMat_BothSig_Mean=nanmean(diff_rsquare_psthMat_BothSig,1);
 diff_rsquare_psthMat_BothSig_Sem=nanstd(diff_rsquare_psthMat_BothSig,[],1)./sqrt(sum(~isnan(diff_rsquare_psthMat_BothSig),1));  
 
 %Check their raw responses
 %Stable Only
 psth_SceneOn_StableSig_Stable=(PSTHSceneOnStable_MeanMat(SSI_Sig_StableOnly,:));

 psth_SceneOn_StableSig_Flexible=(PSTHSceneOnFlexible_MeanMat(SSI_Sig_StableOnly,:));
 
 psth_SceneOn_StableSig_Stable_Mean=nanmean(psth_SceneOn_StableSig_Stable,1);
 psth_SceneOn_StableSig_Stable_Sem=nanstd(psth_SceneOn_StableSig_Stable,[],1)./sqrt(sum(~isnan(psth_SceneOn_StableSig_Stable),1));  
 
 psth_SceneOn_StableSig_Flexible_Mean=nanmean(psth_SceneOn_StableSig_Flexible,1);
 psth_SceneOn_StableSig_Flexible_Sem=nanstd(psth_SceneOn_StableSig_Flexible,[],1)./sqrt(sum(~isnan(psth_SceneOn_StableSig_Flexible),1));  
 
 diff_psth_SceneOn_StableSig=psth_SceneOn_StableSig_Stable-psth_SceneOn_StableSig_Flexible;
 diff_psth_SceneOn_StableSig_Mean=nanmean(diff_psth_SceneOn_StableSig,1);
 diff_psth_SceneOn_StableSig_Sem=nanstd(diff_psth_SceneOn_StableSig,[],1)./sqrt(sum(~isnan(diff_psth_SceneOn_StableSig),1));
 
 %{
 AveMostPrefer_Stable_ZScore=AveMostPreferResponseFromSpon_Stable_ZScore(SSI_Sig_StableOnly);
 
 AveMostPrefer_Flexible_ZScore=AveMostPreferResponseFromSpon_Flexible_ZScore(SSI_Sig_StableOnly);
 
 AveOtherResponse_Stable_ZScore=AveOtherResponseFromSpon_Stable_ZScore(SSI_Sig_StableOnly);
 
 AveOtherResponse_Flexible_ZScore=AveOtherResponseFromSpon_Flexible_ZScore(SSI_Sig_StableOnly);
 
 %}
 
 
 
 %Flexible Only
 psth_SceneOn_FlexibleSig_Stable=(PSTHSceneOnStable_MeanMat(SSI_Sig_FlexibleOnly,:));

 psth_SceneOn_FlexibleSig_Flexible=(PSTHSceneOnFlexible_MeanMat(SSI_Sig_FlexibleOnly,:));
 
 psth_SceneOn_FlexibleSig_Stable_Mean=nanmean(psth_SceneOn_FlexibleSig_Stable,1);
 psth_SceneOn_FlexibleSig_Stable_Sem=nanstd(psth_SceneOn_FlexibleSig_Stable,[],1)./sqrt(sum(~isnan(psth_SceneOn_FlexibleSig_Stable),1));  
 
 psth_SceneOn_FlexibleSig_Flexible_Mean=nanmean(psth_SceneOn_FlexibleSig_Flexible,1);
 psth_SceneOn_FlexibleSig_Flexible_Sem=nanstd(psth_SceneOn_FlexibleSig_Flexible,[],1)./sqrt(sum(~isnan(psth_SceneOn_FlexibleSig_Flexible),1)); 
 
 diff_psth_SceneOn_FlexibleSig=psth_SceneOn_FlexibleSig_Stable-psth_SceneOn_FlexibleSig_Flexible;
 diff_psth_SceneOn_FlexibleSig_Mean=nanmean(diff_psth_SceneOn_FlexibleSig,1);
 diff_psth_SceneOn_FlexibleSig_Sem=nanstd(diff_psth_SceneOn_FlexibleSig,[],1)./sqrt(sum(~isnan(diff_psth_SceneOn_FlexibleSig),1));
 
 %Both
 
 psth_SceneOn_BothSig_Stable=PSTHSceneOnStable_MeanMat(SSI_Sig_Both,:);

 psth_SceneOn_BothSig_Flexible=PSTHSceneOnFlexible_MeanMat(SSI_Sig_Both,:);
 
 psth_SceneOn_BothSig_Stable_Mean=nanmean(psth_SceneOn_BothSig_Stable,1);
 psth_SceneOn_BothSig_Stable_Sem=nanstd(psth_SceneOn_BothSig_Stable,[],1)./sqrt(sum(~isnan(psth_SceneOn_BothSig_Stable),1));  
 
 psth_SceneOn_BothSig_Flexible_Mean=nanmean(psth_SceneOn_BothSig_Flexible,1);
 psth_SceneOn_BothSig_Flexible_Sem=nanstd(psth_SceneOn_BothSig_Flexible,[],1)./sqrt(sum(~isnan(psth_SceneOn_BothSig_Flexible),1)); 
 
 diff_psth_SceneOn_BothSig=psth_SceneOn_BothSig_Stable-psth_SceneOn_BothSig_Flexible;
 diff_psth_SceneOn_BothSig_Mean=nanmean(diff_psth_SceneOn_BothSig,1);
 diff_psth_SceneOn_BothSig_Sem=nanstd(diff_psth_SceneOn_BothSig,[],1)./sqrt(sum(~isnan(diff_psth_SceneOn_BothSig),1));
 
 %% Active Passive Comparison
 %First compare the response; the pearson correlation between active and
 %passive responses
 
  ActiveRawResponseKey={ 'SceneOnMeanMinusSponWhole','UniqueScene'};

%SceneOnDiff500Mean=cell2mat(KeyValues{ismember(KeysLib,'SceneOnDifference_500Mean')});
 Sel=ismember(KeysLib,ActiveRawResponseKey);
 FocusKeyResponse=KeyValues(Sel);
 FocusKey_Sel=KeysLib(Sel);
 
 %Responses
 SceneOnMeanMinusSponWhole=FocusKeyResponse{ismember(FocusKey_Sel,'SceneOnMeanMinusSponWhole')};
  SceneOnMeanMinusSponWhole_ForPassive=SceneOnMeanMinusSponWhole(PassiveIndex==1);
  %Scenes
  
  UniqueSceneActive=FocusKeyResponse{ismember(FocusKey_Sel,'UniqueScene')};
  UniqueSceneActive_ForPassive=UniqueSceneActive(PassiveIndex==1);
  
  %Load the passive raw scene responses
   PassiveRawResponseKey={ 'StabelScene','FlexibleScene','FRSceneMean_Stable','FRSceneMean_Flexible','FRSceneMean_StableScramble','FRSceneMean_FlexibleScramble','p_ContextVSScramble','r_ContextVSScramble'};
 Sel=ismember(KeysLibPassive,PassiveRawResponseKey);
 FocusKeyPassive=KeyValuesPassive(Sel);
 FocusKey_Sel=KeysLibPassive(Sel);
 
 SceneOnMeanMinusSponPassive_Stable=FocusKeyPassive{ismember(FocusKey_Sel,'FRSceneMean_Stable')};
 SceneOnMeanMinusSponPassive_Flexible=FocusKeyPassive{ismember(FocusKey_Sel,'FRSceneMean_Flexible')};
 
 ScenePassive_Stable=FocusKeyPassive{ismember(FocusKey_Sel,'StabelScene')};
 ScenePassive_Flexible=FocusKeyPassive{ismember(FocusKey_Sel,'FlexibleScene')};
 
  
 SceneOnMeanMinusSponPassive=cellfun(@(x,y) [x,y],SceneOnMeanMinusSponPassive_Flexible,SceneOnMeanMinusSponPassive_Stable,'uniform',0);
 
 ScenePassiveForComparison=cellfun(@(x,y) [x;y],ScenePassive_Flexible,ScenePassive_Stable,'uniform',0);
 
 %Select out the overlap between active and passive scenes
OverlapActivePassive=cellfun(@(x,y) ismember(x,y,'legacy'),UniqueSceneActive_ForPassive, ScenePassiveForComparison,'uniform',0);
ResponseActiveCompare=cellfun(@(x,y) x(y),SceneOnMeanMinusSponWhole_ForPassive,OverlapActivePassive,'uniform',0);

tmp=cellfun(@(x,y) corrcoef_out(x,y), SceneOnMeanMinusSponPassive,ResponseActiveCompare,'uniform',0);

r_PassiveActive=cellfun(@(x) x(1), tmp);
p_PassiveActive=cellfun(@(x) x(2), tmp);

%Distribution of scene and scramble
p_ContextVSScramble=cell2mat(FocusKeyPassive{ismember(FocusKey_Sel,'p_ContextVSScramble')});
r_ContextVSScramble=cell2mat(FocusKeyPassive{ismember(FocusKey_Sel,'r_ContextVSScramble')});

%Prop of significant ssi
 PassiveR2Key={'R2_Whole','R2_Scramble','R2_Stable','R2_Flexible','p_SceneOnAnova','p_SceneOnAnovaScramble','p_SceneOnAnovaStable','p_SceneOnAnovaFlexible',...
'PropOfSig_Stable','PropOfSig_Flexible'};

Sel=ismember(KeysLibPassive,PassiveR2Key);
 FocusKeyR2Passive=KeyValuesPassive(Sel);
 FocusKey_Sel=KeysLibPassive(Sel);
 
 PropOfSig_Stable_Passive=cell2mat(FocusKeyR2Passive{ismember(FocusKey_Sel,'PropOfSig_Stable')});
 
 PropOfSig_Flexible_Passive=cell2mat(FocusKeyR2Passive{ismember(FocusKey_Sel,'PropOfSig_Flexible')});
 
  ActiveR2ResponseKey={ 'PropOfSig_Stable','PropOfSig_Flexible'};

%SceneOnDiff500Mean=cell2mat(KeyValues{ismember(KeysLib,'SceneOnDifference_500Mean')});
 Sel=ismember(KeysLib, ActiveR2ResponseKey);
 FocusKeyResponse=KeyValues(Sel);
 FocusKey_Sel=KeysLib(Sel);
 
 %Responses
 PropOfSig_Stable_Active=cell2mat(FocusKeyResponse{ismember(FocusKey_Sel,'PropOfSig_Stable')});
  PropOfSig_Flexible_Active=cell2mat(FocusKeyResponse{ismember(FocusKey_Sel,'PropOfSig_Flexible')});
 
 PropOfSig_Stable_ActiveCompare=PropOfSig_Stable_Active(PassiveIndex==1);
 PropOfSig_Flexible_ActiveCompare=PropOfSig_Flexible_Active(PassiveIndex==1);
 
 %Separate for three type of neurons
 DifferencePropOfSig_Stable=PropOfSig_Stable_ActiveCompare-PropOfSig_Stable_Passive;
 DifferencePropOfSig_Flexibe=PropOfSig_Flexible_ActiveCompare-PropOfSig_Flexible_Passive;
 
 %Stable only
 SSI_Sig_StableOnly_Compare=SSI_Sig_StableOnly(PassiveIndex==1);
 
 PropOfSig_Stable_ActiveCompare_StableOnly=PropOfSig_Stable_ActiveCompare(SSI_Sig_StableOnly_Compare);
 PropOfSig_Flexible_ActiveCompare_StableOnly=PropOfSig_Flexible_ActiveCompare(SSI_Sig_StableOnly_Compare);
 
 PropOfSig_Stable_PassiveCompare_StableOnly=PropOfSig_Stable_Passive(SSI_Sig_StableOnly_Compare);
 PropOfSig_Flexible_PassiveCompare_StableOnly=PropOfSig_Flexible_Passive(SSI_Sig_StableOnly_Compare);
 
 DifferencePropOfSig_StableOnly_Stable=PropOfSig_Stable_ActiveCompare_StableOnly-PropOfSig_Stable_PassiveCompare_StableOnly;
 DifferencePropOfSig_StableOnly_Flexible=PropOfSig_Flexible_ActiveCompare_StableOnly-PropOfSig_Flexible_PassiveCompare_StableOnly;
 
 NumStaCompare=sum(SSI_Sig_StableOnly_Compare)
 
 %Flexible only
 SSI_Sig_FlexibleOnly_Compare=SSI_Sig_FlexibleOnly(PassiveIndex==1);
 
PropOfSig_Stable_ActiveCompare_FlexibleOnly=PropOfSig_Stable_ActiveCompare(SSI_Sig_FlexibleOnly_Compare);
 PropOfSig_Flexible_ActiveCompare_FlexibleOnly=PropOfSig_Flexible_ActiveCompare(SSI_Sig_FlexibleOnly_Compare);
 
 PropOfSig_Stable_PassiveCompare_FlexibleOnly=PropOfSig_Stable_Passive(SSI_Sig_FlexibleOnly_Compare);
 PropOfSig_Flexible_PassiveCompare_FlexibleOnly=PropOfSig_Flexible_Passive(SSI_Sig_FlexibleOnly_Compare);
 
 DifferencePropOfSig_FlexibleOnly_Stable=PropOfSig_Stable_ActiveCompare_FlexibleOnly-PropOfSig_Stable_PassiveCompare_FlexibleOnly;
 DifferencePropOfSig_FlexibleOnly_Flexible=PropOfSig_Flexible_ActiveCompare_FlexibleOnly-PropOfSig_Flexible_PassiveCompare_FlexibleOnly;
 
 NumFleCompare=sum(SSI_Sig_FlexibleOnly_Compare)
 
 %Both Selective cell
 
 SSI_Sig_Both_Compare=SSI_Sig_Both(PassiveIndex==1);
 
PropOfSig_Stable_ActiveCompare_Both=PropOfSig_Stable_ActiveCompare(SSI_Sig_Both_Compare);
 PropOfSig_Flexible_ActiveCompare_Both=PropOfSig_Flexible_ActiveCompare(SSI_Sig_Both_Compare);
 
 PropOfSig_Stable_PassiveCompare_Both=PropOfSig_Stable_Passive(SSI_Sig_Both_Compare);
 PropOfSig_Flexible_PassiveCompare_Both=PropOfSig_Flexible_Passive(SSI_Sig_Both_Compare);
 
 DifferencePropOfSig_Both_Stable=PropOfSig_Stable_ActiveCompare_Both-PropOfSig_Stable_PassiveCompare_Both;
 DifferencePropOfSig_Both_Flexible=PropOfSig_Flexible_ActiveCompare_Both-PropOfSig_Flexible_PassiveCompare_Both;
 
 NumBothCompare=sum(SSI_Sig_Both_Compare)
 
 %None Selective Cell
 
 SSI_Sig_None_Compare=SSI_Sig_None(PassiveIndex==1);
 
PropOfSig_Stable_ActiveCompare_None=PropOfSig_Stable_ActiveCompare(SSI_Sig_None_Compare);
 PropOfSig_Flexible_ActiveCompare_None=PropOfSig_Flexible_ActiveCompare(SSI_Sig_None_Compare);
 
 PropOfSig_Stable_PassiveCompare_None=PropOfSig_Stable_Passive(SSI_Sig_None_Compare);
 PropOfSig_Flexible_PassiveCompare_None=PropOfSig_Flexible_Passive(SSI_Sig_None_Compare);
 
 DifferencePropOfSig_None_Stable=PropOfSig_Stable_ActiveCompare_None-PropOfSig_Stable_PassiveCompare_None;
 DifferencePropOfSig_None_Flexible=PropOfSig_Flexible_ActiveCompare_None-PropOfSig_Flexible_PassiveCompare_None;
 
NumNoneCompare=sum(SSI_Sig_None_Compare)
 

%Plot One example

if PlotExample
ExampleName='Robin111620FORAGEXPLO1';
ExampleChannel=10;

ExampleSel_tmp=strcmp(FilesNameEach,ExampleName) & ChannelFull==ExampleChannel ;
ExampleSel=ExampleSel_tmp(PassiveIndex==1);

ScenePassiveExample=ScenePassiveForComparison{ExampleSel};

ResponseActive_Example=ResponseActiveCompare{ExampleSel};
ResponsePassive_Example=SceneOnMeanMinusSponPassive{ExampleSel};

r_example=r_PassiveActive(ExampleSel)
p_example=p_PassiveActive(ExampleSel)

PropOfSig_Stable_Passive_Example=PropOfSig_Stable_Passive(ExampleSel)
PropOfSig_Flexible_Passive_Example=PropOfSig_Flexible_Passive(ExampleSel)

PropOfSig_Stable_Active_Example=PropOfSig_Stable_ActiveCompare(ExampleSel)
PropOfSig_Flexible_Active_Example=PropOfSig_Flexible_ActiveCompare(ExampleSel)


end
 %%
 %Plots
 %%%%%%% Basic response properties%%%%%%%%
FigureStartNum=500;
FigureIndex=1;
if ShowFigureFlag
%%First Figure: Basic response comparisons
figtitlestr{FigureIndex}='BasicResponseProperties';
fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[50,100, 1000,800],'Name',figtitlestr{FigureIndex});

subplot(2,2,1);
histogram(TableData(:,1));
xlabel('AverageFRInTask');
ylabel('Number of neurons');
title('Distribution of the Average FR in task');
set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');

subplot(2,2,2);
histogram(TableData(:,2));
xlabel('AverageSponITI');
ylabel('Number of neurons');
title('Distribution of the Average Spon FR in task');
set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');

subplot(2,2,3);
plot(TableData(:,3),TableData(:,4),'or');
hold on 
x=0:max(TableData(:,[3:4]),[],'all');
y=x;
plot(x,y,'--k');
%Fill in the dots that are significant;
SigPairs=TableData(TableData(:,5)<0.05,[3:4]);
plot(SigPairs(:,1),SigPairs(:,2),'.r','MarkerSize',15);

xlabel('Flexible');
ylabel('Stable');
axis equal
set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');

subplot(2,2,4);
%histogram((TableData(:,3)-TableData(:,4))./TableData(:,2));
%ValueForHistogram=(TableData(:,3)-TableData(:,4))./(TableData(:,3)+TableData(:,4));
ValueForHistogram=(TableData(:,3)-TableData(:,4))./TableData(:,2);
pVal=TableData(:,5);
SigVal=ValueForHistogram(pVal<0.05);
NonSigVal=ValueForHistogram(~(pVal<0.05));
%SepBin=[linspace(min(ValueForHistogram)-1,0,10),linspace(0.1,max(ValueForHistogram)+1,10)];
SepBin=[linspace(min(ValueForHistogram)-1,max(ValueForHistogram)+1,20)];

Val_Sig_Sep=histcounts(SigVal,SepBin);
NonVal_Sig_Sep=histcounts(NonSigVal,SepBin);
BinCenter=SepBin(1:end-1)+diff(SepBin);

b=bar([BinCenter',BinCenter'],[Val_Sig_Sep',NonVal_Sig_Sep'],'stacked');
b(1).FaceColor=[0,0,0];
b(2).FaceColor=[1,1,1];

%histogram((TableData(:,3)-TableData(:,4))./(TableData(:,3)+TableData(:,4)));
xlabel('Flexible-Stable');
set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
FigureIndex=FigureIndex+1;
FigureStartNum=FigureStartNum+1; 
    

%%Figure2: Separate scene responses for difference type of neurons 

figtitlestr{FigureIndex}='PSTH_SceneOnEachType';
fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[70,100, 1000,800],'Name',figtitlestr{FigureIndex});

subplot(3,2,1)
shadedErrorBar(PSTHSceneOn_TimeMean,PSTHSceneOnFlexible_FlexDo_Mean, PSTHSceneOnFlexible_FlexDo_Sem,'lineprops',{FlexibleSceneColor},'transparent',1,'patchSaturation',0.5);%PSTH for Flexible
hold on
shadedErrorBar(PSTHSceneOn_TimeMean,PSTHSceneOnStable_FlexDo_Mean, PSTHSceneOnStable_FlexDo_Sem,'lineprops',{StableSceneColor},'transparent',1,'patchSaturation',0.5);%PSTH for Stable
  set(findall(gca, 'Type', 'Line'),'LineWidth',3);
  
  ylabel('Z-Score');
title('Flexible Dominance Type');
set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');

subplot(3,2,3)
shadedErrorBar(PSTHSceneOn_TimeMean,PSTHSceneOnFlexible_StaDo_Mean, PSTHSceneOnFlexible_StaDo_Sem,'lineprops',{FlexibleSceneColor},'transparent',1,'patchSaturation',0.5);%PSTH for Flexible
hold on
shadedErrorBar(PSTHSceneOn_TimeMean,PSTHSceneOnStable_StaDo_Mean, PSTHSceneOnStable_StaDo_Sem,'lineprops',{StableSceneColor},'transparent',1,'patchSaturation',0.5);%PSTH for Stable
  set(findall(gca, 'Type', 'Line'),'LineWidth',3);
title('Stable Dominance Type');
ylabel('Z-Score');
set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');

subplot(3,2,5)
shadedErrorBar(PSTHSceneOn_TimeMean,PSTHSceneOnFlexible_NoDo_Mean, PSTHSceneOnFlexible_NoDo_Sem,'lineprops',{FlexibleSceneColor},'transparent',1,'patchSaturation',0.5);%PSTH for Flexible
hold on
shadedErrorBar(PSTHSceneOn_TimeMean,PSTHSceneOnStable_NoDo_Mean, PSTHSceneOnStable_NoDo_Sem,'lineprops',{StableSceneColor},'transparent',1,'patchSaturation',0.5);%PSTH for Stable
set(findall(gca, 'Type', 'Line'),'LineWidth',3);
title('No Dominance Type');
xlabel('Time from scene onset');
ylabel('Z-Score');
set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');

%Scene offset
subplot(3,2,2)
shadedErrorBar(PSTHSceneOff_TimeMean,PSTHSceneOffFlexible_FlexDo_Mean, PSTHSceneOffFlexible_FlexDo_Sem,'lineprops',{FlexibleSceneColor},'transparent',1,'patchSaturation',0.5);%PSTH for Flexible
hold on
shadedErrorBar(PSTHSceneOff_TimeMean,PSTHSceneOffStable_FlexDo_Mean, PSTHSceneOffStable_FlexDo_Sem,'lineprops',{StableSceneColor},'transparent',1,'patchSaturation',0.5);%PSTH for Stable
  set(findall(gca, 'Type', 'Line'),'LineWidth',3);
  ylabel('Z-Score');
title('Flexible Dominance Type');
set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');


subplot(3,2,4)
shadedErrorBar(PSTHSceneOff_TimeMean,PSTHSceneOffFlexible_StaDo_Mean, PSTHSceneOffFlexible_StaDo_Sem,'lineprops',{FlexibleSceneColor},'transparent',1,'patchSaturation',0.5);%PSTH for Flexible
hold on
shadedErrorBar(PSTHSceneOff_TimeMean,PSTHSceneOffStable_StaDo_Mean, PSTHSceneOffStable_StaDo_Sem,'lineprops',{StableSceneColor},'transparent',1,'patchSaturation',0.5);%PSTH for Stable
  set(findall(gca, 'Type', 'Line'),'LineWidth',3);
 ylabel('Z-Score');
title('Stable Dominance Type');
set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');


subplot(3,2,6)
shadedErrorBar(PSTHSceneOff_TimeMean,PSTHSceneOffFlexible_NoDo_Mean, PSTHSceneOffFlexible_NoDo_Sem,'lineprops',{FlexibleSceneColor},'transparent',1,'patchSaturation',0.5);%PSTH for Flexible
hold on
shadedErrorBar(PSTHSceneOff_TimeMean,PSTHSceneOffStable_NoDo_Mean, PSTHSceneOffStable_NoDo_Sem,'lineprops',{StableSceneColor},'transparent',1,'patchSaturation',0.5);%PSTH for Stable
set(findall(gca, 'Type', 'Line'),'LineWidth',3);
 xlabel('Time from scene onset');
 ylabel('Z-Score');
  title('No Dominance Type');
set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');

FigureIndex=FigureIndex+1;
FigureStartNum=FigureStartNum+1; 


%{
%%Figure3: Scene Selective Index 


figtitlestr{FigureIndex}='SceneSelectiveIndex';
fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[90,100, 500,500],'Name',figtitlestr{FigureIndex});


%Flexible Dominance
subplot(2,2,1);
ValueForHistogram=SISceneOn_Fle;
pVal=sigSceneOn_Fle;
SigVal=ValueForHistogram(pVal<0.05);
NonSigVal=ValueForHistogram(~(pVal<0.05));
%SepBin=[linspace(min(ValueForHistogram)-1,0,10),linspace(0.1,max(ValueForHistogram)+1,10)];
SepBin=linspace(0,1,20);

Val_Sig_Sep=histcounts(SigVal,SepBin);
NonVal_Sig_Sep=histcounts(NonSigVal,SepBin);
BinCenter=SepBin(1:end-1)+diff(SepBin);

b=bar([BinCenter',BinCenter'],[Val_Sig_Sep',NonVal_Sig_Sep'],'stacked');
b(1).FaceColor=[0,0,0];
b(2).FaceColor=[1,1,1];
title('SI For Flexible Scene');
xlabel('SI');
ylabel('Firing Rate(Hz)');
set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');

%Distribition of SI on stable scene
subplot(2,2,2);

ValueForHistogram=SISceneOn_Sta;
pVal=sigSceneOn_Sta;
SigVal=ValueForHistogram(pVal<0.05);
NonSigVal=ValueForHistogram(~(pVal<0.05));
%SepBin=[linspace(min(ValueForHistogram)-1,0,10),linspace(0.1,max(ValueForHistogram)+1,10)];
SepBin=linspace(0,1,20);

Val_Sig_Sep=histcounts(SigVal,SepBin);
NonVal_Sig_Sep=histcounts(NonSigVal,SepBin);
BinCenter=SepBin(1:end-1)+diff(SepBin);

b=bar([BinCenter',BinCenter'],[Val_Sig_Sep',NonVal_Sig_Sep'],'stacked');
b(1).FaceColor=[0,0,0];
b(2).FaceColor=[1,1,1];
title('SI For Stable Scene');
xlabel('SI');
ylabel('Firing Rate(Hz)');
set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');


%sigSceneOn_Fle

subplot(2,2,3);
EitherOneSig=(sigSceneOn_Sta<0.05)|(sigSceneOn_Fle<0.05);

plot(SISceneOn_Fle(EitherOneSig),SISceneOn_Sta(EitherOneSig),'.r','MarkerSize',25);
hold on
plot(SISceneOn_Fle(~EitherOneSig),SISceneOn_Sta(~EitherOneSig),'or','MarkerSize',10);
axis equal;
hold on 
x=0:0.1:1;
y=x;
plot(x,y,'--k');

xlabel('SI SceneOn Flexible');
ylabel('SI SceneOn Stable');
title('Scene Selective Index');

set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');

FigureIndex=FigureIndex+1;
FigureStartNum=FigureStartNum+1; 
%}
%%
% Figure 4: Histogram of number of significant scenes 
figtitlestr{FigureIndex}='HistogramOfSigScenes';
fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[110,100, 800,500],'Name',figtitlestr{FigureIndex});
%Overall
subplot(2,2,1)
[counts,edge]=histcounts(NumOfSignificantScene500,[0:17]);

b=bar(edge(1:end-1),counts,'FaceColor','#94969A','EdgeColor','k','LineWidth',3);

%histogram(NumOfSignificantScene500,'BinWidth',1,'FaceColor','k','LineWidth',3);
xlabel('Number of significant scenes');
ylabel('Number of neurons');
set(gca,'LineWidth',3,'FontSize',20,'FontWeight','Bold');
box off
xlim([-1,17]);
%Stable Scenes
subplot(2,2,2)
[counts,edge]=histcounts(NumOfStableSigScenes,[0:9]);

b=bar(edge(1:end-1),counts,'FaceColor',StableSceneColor,'EdgeColor','k','LineWidth',3);

%histogram(NumOfStableSigScenes,'BinWidth',1,'FaceColor',StableSceneColor,'LineWidth',3);
xlabel('Number of significant stable scenes');
ylabel('Number of neurons');
set(gca,'LineWidth',3,'FontSize',20,'FontWeight','Bold');
box off
xlim([-1,9]);

%Flexible Scenes
subplot(2,2,3)
[counts,edge]=histcounts(NumOfFlexibleSigScenes,[0:9]);

b=bar(edge(1:end-1),counts,'FaceColor',FlexibleSceneColor,'EdgeColor','k','LineWidth',3);

%histogram(NumOfFlexibleSigScenes,'BinWidth',1,'FaceColor',FlexibleSceneColor,'LineWidth',3);
xlabel('Number of significant flexible scenes');
ylabel('Number of neurons');
set(gca,'LineWidth',3,'FontSize',20,'FontWeight','Bold');
box off
xlim([-1,9]);

%Proportion of neuron types
subplot(2,2,4)
b=bar([1,2,3,4],[PropOfNeuronOnlyStable,PropOfNeuronOnlyFlexible,PropOfNeuronBoth,PropOfNeuronNone],'LineWidth',3);

b.FaceColor = 'flat';
%b.EdgeColor = 'flat';
ColorChosen=[StableSceneColor;FlexibleSceneColor;[0.9843    0.2549    0.3647];[ 0.5216    0.5137    0.5137]];

%g.CData = ColorChosen([Index_Sort,NumOfCondition+1],:);
b.CData = ColorChosen;
xlim([0.2,4+0.8]);
box off
xticks(1:4);
%xticklabels(SceneStr(Index_Sort));
xticklabels({'StaOnly','FleOnly','Both','None'});
ylabel('Prop of neurons');

set(gca,'LineWidth',3,'FontSize',20,'FontWeight','Bold');


FigureIndex=FigureIndex+1;
FigureStartNum=FigureStartNum+1; 


% Figure 5: ResponseAverageBetweenStableAndFlexible_MostPrefered
figtitlestr{FigureIndex}='MostPreferSceneStableVSFlexible';
fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[120,100, 800,400],'Name',figtitlestr{FigureIndex});

subplot(1,2,1)
title('Most prefer stable v.s. flexible')
plot(AveMostPreferResponseFromSpon_Stable_abs,AveMostPreferResponseFromSpon_Flexible_abs,'or','LineWidth',3,'MarkerSize',15);
%Add a diagnol axis
x=1:max([AveMostPreferResponseFromSpon_Stable_abs,AveMostPreferResponseFromSpon_Flexible_abs])+1;
y=x;
hold on
plot(x,y,'--k','LineWidth',3);
axis equal
xlabel('Most prefer stable response');
ylabel('Most prefer flexible response');
set(gca,'LineWidth',3,'FontSize',20,'FontWeight','Bold');

subplot(1,2,2)
title('Other stable v.s. flexible')
plot(AveOtherResponseFromSpon_Stable_abs,AveOtherResponseFromSpon_Flexible_abs,'or','LineWidth',3,'MarkerSize',15);
%Add a diagnol axis
x=1:max([AveOtherResponseFromSpon_Stable_abs,AveOtherResponseFromSpon_Flexible_abs])+1;
y=x;
hold on
plot(x,y,'--k','LineWidth',3);
axis equal
xlabel('Other stable response');
ylabel('Other flexible response');
set(gca,'LineWidth',3,'FontSize',20,'FontWeight','Bold');



FigureIndex=FigureIndex+1;
FigureStartNum=FigureStartNum+1; 

%% Figure 6 Distribution of the R2 under all conditions
figtitlestr{FigureIndex}='PSTH_SceneOnEachType';
fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[70,100, 1000,800],'Name',figtitlestr{FigureIndex});

subplot(2,2,1);
ValueForHistogram=R2SSI_Whole;
pVal=p_SSI_Whole;
SigVal=ValueForHistogram(pVal<0.05);
NonSigVal=ValueForHistogram(~(pVal<0.05));
%SepBin=[linspace(min(ValueForHistogram)-1,0,10),linspace(0.1,max(ValueForHistogram)+1,10)];
SepBin=[linspace(0,1,20)];

Val_Sig_Sep=histcounts(SigVal,SepBin);
NonVal_Sig_Sep=histcounts(NonSigVal,SepBin);
BinCenter=SepBin(1:end-1)+diff(SepBin);

b=bar([BinCenter',BinCenter'],[Val_Sig_Sep',NonVal_Sig_Sep'],'stacked');
b(1).FaceColor=[0,0,0];
b(2).FaceColor=[1,1,1];

%histogram((TableData(:,3)-TableData(:,4))./(TableData(:,3)+TableData(:,4)));
xlabel('Distribution of SSI(Whole)');
ylabel('Number of neurons');
set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
box off;

subplot(2,2,2);
ValueForHistogram=R2SSI_Sta;

pVal=p_SSI_Sta;
SigVal=ValueForHistogram(pVal<0.05);
NonSigVal=ValueForHistogram(~(pVal<0.05));
%SepBin=[linspace(min(ValueForHistogram)-1,0,10),linspace(0.1,max(ValueForHistogram)+1,10)];
SepBin=[linspace(0,1,20)];

Val_Sig_Sep=histcounts(SigVal,SepBin);
NonVal_Sig_Sep=histcounts(NonSigVal,SepBin);
BinCenter=SepBin(1:end-1)+diff(SepBin);

b=bar([BinCenter',BinCenter'],[Val_Sig_Sep',NonVal_Sig_Sep'],'stacked');
b(1).FaceColor=StableSceneColor;
b(2).FaceColor=[1,1,1];

%histogram((TableData(:,3)-TableData(:,4))./(TableData(:,3)+TableData(:,4)));
xlabel('Distribution of SSI(Sta)');
ylabel('Number of neurons');
set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
box off;

subplot(2,2,3);
ValueForHistogram=R2SSI_Fle;

pVal=p_SSI_Fle;
SigVal=ValueForHistogram(pVal<0.05);
NonSigVal=ValueForHistogram(~(pVal<0.05));
%SepBin=[linspace(min(ValueForHistogram)-1,0,10),linspace(0.1,max(ValueForHistogram)+1,10)];
SepBin=[linspace(0,1,20)];

Val_Sig_Sep=histcounts(SigVal,SepBin);
NonVal_Sig_Sep=histcounts(NonSigVal,SepBin);
BinCenter=SepBin(1:end-1)+diff(SepBin);

b=bar([BinCenter',BinCenter'],[Val_Sig_Sep',NonVal_Sig_Sep'],'stacked');
b(1).FaceColor=FlexibleSceneColor;
b(2).FaceColor=[1,1,1];

%histogram((TableData(:,3)-TableData(:,4))./(TableData(:,3)+TableData(:,4)));
xlabel('Distribution of SSI(Fle)');
ylabel('Number of neurons');
set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
box off;

subplot(2,2,4);
%Comparison between Sta and Fle
%{
SSI_Sig_StableOnly=p_SSI_Sta<0.05 & p_SSI_Fle>=0.05 ;
 SSI_Sig_FlexibleOnly=p_SSI_Sta>=0.05 & p_SSI_Fle<0.05;
 SSI_Sig_Both=p_SSI_Sta<0.05 & p_SSI_Fle<0.05;
 SSI_Sig_None=p_SSI_Sta>=0.05 & p_SSI_Fle>=0.05;
%}
%First Stable Only
MarkerSize=150;
NoNeColor='#C8C7C9';

scatter(R2SSI_Fle(SSI_Sig_None),R2SSI_Sta(SSI_Sig_None),MarkerSize*ones(1,sum(SSI_Sig_None)),'o','MarkerFaceColor', NoNeColor,'MarkerEdgeColor','w');

hold on
scatter(R2SSI_Fle(SSI_Sig_StableOnly),R2SSI_Sta(SSI_Sig_StableOnly),MarkerSize*ones(sum(SSI_Sig_StableOnly),1),'MarkerFaceColor', StableSceneColor,'MarkerEdgeColor','w');
hold on
scatter(R2SSI_Fle(SSI_Sig_FlexibleOnly),R2SSI_Sta(SSI_Sig_FlexibleOnly),MarkerSize*ones(sum(SSI_Sig_FlexibleOnly),1),'MarkerFaceColor', FlexibleSceneColor,'MarkerEdgeColor','w');
scatter(R2SSI_Fle(SSI_Sig_Both),R2SSI_Sta(SSI_Sig_Both),MarkerSize*ones(sum(SSI_Sig_Both),1),'^','MarkerFaceColor', '#8F6BD6','MarkerEdgeColor','w');



plot([0:0.1:1],[0:0.1:1],'--k','LineWidth',3);
xlabel('SSI Fle');
ylabel('SSI Sta');
set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
box off;

FigureIndex=FigureIndex+1;
FigureStartNum=FigureStartNum+1; 


%% Figure 7 PSTH of R2 
figtitlestr{FigureIndex}='PSTH_SceneOnRSquare';
fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[90,100, 1000,800],'Name',figtitlestr{FigureIndex});
subplot(2,1,1);
shadedErrorBar(PSTHSceneOn_TimeMean,rsquare_psthMat_WholeSig_Mean, rsquare_psthMat_WholeSig_Sem,'lineprops',{'k'},'transparent',1,'patchSaturation',0.5);%PSTH for Flexible
set(findall(gca, 'Type', 'Line'),'LineWidth',3);
hold on
shadedErrorBar(PSTHSceneOn_TimeMean,rsquare_psthMat_WholeNonSig_Mean, rsquare_psthMat_WholeNonSig_Sem,'lineprops',{[0.6980    0.6863    0.6863]},'transparent',1,'patchSaturation',0.5);%PSTH for Flexible
set(findall(gca, 'Type', 'Line'),'LineWidth',3);
xlim([-500,1500]);


set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
xlabel('Time from scene onset');
ylabel('SSI');
title('All scene selective neurons');

subplot(2,1,2);
shadedErrorBar(PSTHSceneOn_TimeMean,rsquare_psthMat_WholeSig_StableMean, rsquare_psthMat_WholeSig_StableSem,'lineprops',{StableSceneColor},'transparent',1,'patchSaturation',0.5);%PSTH for Flexible
hold on
shadedErrorBar(PSTHSceneOn_TimeMean,rsquare_psthMat_WholeSig_FlexibleMean, rsquare_psthMat_WholeSig_FlexibleSem,'lineprops',{FlexibleSceneColor},'transparent',1,'patchSaturation',0.5);%PSTH for Flexible
set(findall(gca, 'Type', 'Line'),'LineWidth',3);
set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
xlim([-500,1500]);
xlabel('Time from scene onset');
ylabel('SSI');
title('All scene selective neurons');
set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
box off;


FigureIndex=FigureIndex+1;
FigureStartNum=FigureStartNum+1; 
    

%% Figure 8 PSTH of R2 under each condition

figtitlestr{FigureIndex}='PSTH_SceneOnRSquareForEachType';
fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[100,100, 1000,800],'Name',figtitlestr{FigureIndex});

%StableOnly
subplot(3,2,1);
shadedErrorBar(PSTHSceneOn_TimeMean,rsquare_psthMat_StableSig_Stable_Mean, rsquare_psthMat_StableSig_Stable_Sem,'lineprops',{StableSceneColor},'transparent',1,'patchSaturation',0.5);%PSTH for Flexible
set(findall(gca, 'Type', 'Line'),'LineWidth',3);
hold on
shadedErrorBar(PSTHSceneOn_TimeMean,rsquare_psthMat_StableSig_Flexible_Mean, rsquare_psthMat_StableSig_Flexible_Sem,'lineprops',{FlexibleSceneColor},'transparent',1,'patchSaturation',0.5);%PSTH for Flexible
set(findall(gca, 'Type', 'Line'),'LineWidth',3);
xlim([-500,1500]);

xlabel('Time from scene onset');
ylabel('SSI');
title(sprintf('Stable Selective Only (N= %d)',sum(SSI_Sig_StableOnly)));
set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');


%Flexible Only
subplot(3,2,3);
shadedErrorBar(PSTHSceneOn_TimeMean,rsquare_psthMat_FlexibleSig_Stable_Mean, rsquare_psthMat_FlexibleSig_Stable_Sem,'lineprops',{StableSceneColor},'transparent',1,'patchSaturation',0.5);%PSTH for Flexible
set(findall(gca, 'Type', 'Line'),'LineWidth',3);
hold on
shadedErrorBar(PSTHSceneOn_TimeMean,rsquare_psthMat_FlexibleSig_Flexible_Mean, rsquare_psthMat_FlexibleSig_Flexible_Sem,'lineprops',{FlexibleSceneColor},'transparent',1,'patchSaturation',0.5);%PSTH for Flexible
set(findall(gca, 'Type', 'Line'),'LineWidth',3);

xlim([-500,1500]);
xlabel('Time from scene onset');
ylabel('SSI');
title(sprintf('Flexible Selective Only (N= %d)',sum(SSI_Sig_FlexibleOnly)));
set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');


%Both
subplot(3,2,5);
shadedErrorBar(PSTHSceneOn_TimeMean,rsquare_psthMat_BothSig_Stable_Mean, rsquare_psthMat_BothSig_Stable_Sem,'lineprops',{StableSceneColor},'transparent',1,'patchSaturation',0.5);%PSTH for Flexible
set(findall(gca, 'Type', 'Line'),'LineWidth',3);
hold on
shadedErrorBar(PSTHSceneOn_TimeMean,rsquare_psthMat_BothSig_Flexible_Mean, rsquare_psthMat_BothSig_Flexible_Sem,'lineprops',{FlexibleSceneColor},'transparent',1,'patchSaturation',0.5);%PSTH for Flexible
set(findall(gca, 'Type', 'Line'),'LineWidth',3);

xlim([-500,1500]);
xlabel('Time from scene onset');
ylabel('SSI');
title(sprintf('Both Selective (N= %d)',sum(SSI_Sig_Both)));
set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');

%Alos the difference
%Stable only difference
subplot(3,2,2);

shadedErrorBar(PSTHSceneOn_TimeMean,diff_rsquare_psthMat_StableSig_Mean, diff_rsquare_psthMat_StableSig_Sem,'lineprops',{DifferenceColor},'transparent',1,'patchSaturation',0.5);%PSTH for Flexible
set(findall(gca, 'Type', 'Line'),'LineWidth',3);


xlim([-500,1500]);
xlabel('Time from scene onset');
ylabel('Differene of SSI');
title(sprintf('Stable Selective (N= %d)',sum(SSI_Sig_StableOnly)));
set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');



subplot(3,2,4);

shadedErrorBar(PSTHSceneOn_TimeMean,diff_rsquare_psthMat_FlexibleSig_Mean, diff_rsquare_psthMat_FlexibleSig_Sem,'lineprops',{DifferenceColor},'transparent',1,'patchSaturation',0.5);%PSTH for Flexible
set(findall(gca, 'Type', 'Line'),'LineWidth',3);

xlim([-500,1500]);
xlabel('Time from scene onset');
ylabel('Differene of SSI');
title(sprintf('Flexible Selective (N= %d)',sum(SSI_Sig_FlexibleOnly)));
set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');


subplot(3,2,6);

shadedErrorBar(PSTHSceneOn_TimeMean,diff_rsquare_psthMat_BothSig_Mean, diff_rsquare_psthMat_BothSig_Sem,'lineprops',{DifferenceColor},'transparent',1,'patchSaturation',0.5);%PSTH for Flexible
set(findall(gca, 'Type', 'Line'),'LineWidth',3);

xlim([-500,1500]);
xlabel('Time from scene onset');
ylabel('Differene of SSI');
title(sprintf('Both Selective (N= %d)',sum(SSI_Sig_Both)));
set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
FigureIndex=FigureIndex+1;
FigureStartNum=FigureStartNum+1;

%% Figure 9 PSTH sparated by SSI under each condition

figtitlestr{FigureIndex}='PSTH_SceneOnResponseSepSSForEachType';
fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[100,100, 1000,800],'Name',figtitlestr{FigureIndex});

%StableOnly
subplot(3,2,1);
shadedErrorBar(PSTHSceneOn_TimeMean,psth_SceneOn_StableSig_Stable_Mean, psth_SceneOn_StableSig_Stable_Sem,'lineprops',{StableSceneColor},'transparent',1,'patchSaturation',0.5);%PSTH for Flexible
set(findall(gca, 'Type', 'Line'),'LineWidth',3);
hold on
shadedErrorBar(PSTHSceneOn_TimeMean,psth_SceneOn_StableSig_Flexible_Mean, psth_SceneOn_StableSig_Flexible_Sem,'lineprops',{FlexibleSceneColor},'transparent',1,'patchSaturation',0.5);%PSTH for Flexible
set(findall(gca, 'Type', 'Line'),'LineWidth',3);
xlim([-500,1500]);

xlabel('Time from scene onset');
ylabel('Z-Score');
title(sprintf('Stable Selective Only (N= %d)',sum(SSI_Sig_StableOnly)));
set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');


%Flexible Only
subplot(3,2,3);
shadedErrorBar(PSTHSceneOn_TimeMean,psth_SceneOn_FlexibleSig_Stable_Mean, psth_SceneOn_FlexibleSig_Stable_Sem,'lineprops',{StableSceneColor},'transparent',1,'patchSaturation',0.5);%PSTH for Flexible
set(findall(gca, 'Type', 'Line'),'LineWidth',3);
hold on
shadedErrorBar(PSTHSceneOn_TimeMean,psth_SceneOn_FlexibleSig_Flexible_Mean, psth_SceneOn_FlexibleSig_Flexible_Sem,'lineprops',{FlexibleSceneColor},'transparent',1,'patchSaturation',0.5);%PSTH for Flexible
set(findall(gca, 'Type', 'Line'),'LineWidth',3);

xlim([-500,1500]);
xlabel('Time from scene onset');
ylabel('Z-Score');
title(sprintf('Flexible Selective Only (N= %d)',sum(SSI_Sig_FlexibleOnly)));
set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');


%Both
subplot(3,2,5);
shadedErrorBar(PSTHSceneOn_TimeMean,psth_SceneOn_BothSig_Stable_Mean, psth_SceneOn_BothSig_Stable_Sem,'lineprops',{StableSceneColor},'transparent',1,'patchSaturation',0.5);%PSTH for Flexible
set(findall(gca, 'Type', 'Line'),'LineWidth',3);
hold on
shadedErrorBar(PSTHSceneOn_TimeMean,psth_SceneOn_BothSig_Flexible_Mean, psth_SceneOn_BothSig_Flexible_Sem,'lineprops',{FlexibleSceneColor},'transparent',1,'patchSaturation',0.5);%PSTH for Flexible
set(findall(gca, 'Type', 'Line'),'LineWidth',3);

xlim([-500,1500]);
xlabel('Time from scene onset');
ylabel('Z-Score');
title(sprintf('Both Selective (N= %d)',sum(SSI_Sig_Both)));
set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');

%Alos the difference
%Stable only difference
subplot(3,2,2);

shadedErrorBar(PSTHSceneOn_TimeMean,diff_psth_SceneOn_StableSig_Mean, diff_psth_SceneOn_StableSig_Sem,'lineprops',{DifferenceColor},'transparent',1,'patchSaturation',0.5);%PSTH for Flexible
set(findall(gca, 'Type', 'Line'),'LineWidth',3);


xlim([-500,1500]);
xlabel('Time from scene onset');
ylabel('Differene of responses');
title(sprintf('Stable Selective (N= %d)',sum(SSI_Sig_StableOnly)));
set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');



subplot(3,2,4);

shadedErrorBar(PSTHSceneOn_TimeMean,diff_psth_SceneOn_FlexibleSig_Mean, diff_psth_SceneOn_FlexibleSig_Sem,'lineprops',{DifferenceColor},'transparent',1,'patchSaturation',0.5);%PSTH for Flexible
set(findall(gca, 'Type', 'Line'),'LineWidth',3);

xlim([-500,1500]);
xlabel('Time from scene onset');
ylabel('Differene of responses');
title(sprintf('Flexible Selective (N= %d)',sum(SSI_Sig_FlexibleOnly)));
set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');


subplot(3,2,6);

shadedErrorBar(PSTHSceneOn_TimeMean,diff_psth_SceneOn_BothSig_Mean, diff_psth_SceneOn_BothSig_Sem,'lineprops',{DifferenceColor},'transparent',1,'patchSaturation',0.5);%PSTH for Flexible
set(findall(gca, 'Type', 'Line'),'LineWidth',3);

xlim([-500,1500]);
xlabel('Time from scene onset');
ylabel('Differene of responses');
title(sprintf('Both Selective (N= %d)',sum(SSI_Sig_Both)));
set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');


FigureIndex=FigureIndex+1;
FigureStartNum=FigureStartNum+1;


%% Passive and Active comparison

figtitlestr{FigureIndex}='PearsonCorrPassiveActiveSvsScramble';
fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[100,100, 1000,400],'Name',figtitlestr{FigureIndex});
subplot(1,2,1);
ValueForHistogram=r_PassiveActive;
pVal=p_PassiveActive;
SigVal=ValueForHistogram(pVal<0.05);
NonSigVal=ValueForHistogram(~(pVal<0.05));
%SepBin=[linspace(min(ValueForHistogram)-1,0,10),linspace(0.1,max(ValueForHistogram)+1,10)];
SepBin=linspace(-1,1,30);

MeanVal=nanmean(ValueForHistogram)
%Compare with zero
[pRPassiveActiveResponse,h]=ttest(ValueForHistogram,0)

Val_Sig_Sep=histcounts(SigVal,SepBin);
NonVal_Sig_Sep=histcounts(NonSigVal,SepBin);
BinCenter=SepBin(1:end-1)+diff(SepBin);

b=bar([BinCenter',BinCenter'],[Val_Sig_Sep',NonVal_Sig_Sep'],'stacked');
b(1).FaceColor=[0,0,0];
b(2).FaceColor=[1,1,1];
yrange=ylim;
hold on
plot(MeanVal,yrange(2),'vk','MarkerFaceColor','k','MarkerSize',20);
hold on
plot([MeanVal,MeanVal],yrange,'--k','LineWidth',3);


%histogram((TableData(:,3)-TableData(:,4))./(TableData(:,3)+TableData(:,4)));
xlabel('Distribution of Pearson Correlation(PV)');
ylabel('Number of neurons');
set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
box off;

subplot(1,2,2);
ValueForHistogram=r_ContextVSScramble;
pVal=p_ContextVSScramble;
SigVal=ValueForHistogram(pVal<0.05);
NonSigVal=ValueForHistogram(~(pVal<0.05));
%SepBin=[linspace(min(ValueForHistogram)-1,0,10),linspace(0.1,max(ValueForHistogram)+1,10)];
SepBin=linspace(-1,1,30);

MeanVal=nanmean(ValueForHistogram)
%Compare with zero
[pContextvsScramble,h]=ttest(ValueForHistogram,0)

Val_Sig_Sep=histcounts(SigVal,SepBin);
NonVal_Sig_Sep=histcounts(NonSigVal,SepBin);
BinCenter=SepBin(1:end-1)+diff(SepBin);

b=bar([BinCenter',BinCenter'],[Val_Sig_Sep',NonVal_Sig_Sep'],'stacked');
b(1).FaceColor=[0,0,0];
b(2).FaceColor=[1,1,1];
yrange=ylim;
hold on
plot(MeanVal,yrange(2),'vk','MarkerFaceColor','k','MarkerSize',20);
hold on
plot([MeanVal,MeanVal],yrange,'--k','LineWidth',3);


%histogram((TableData(:,3)-TableData(:,4))./(TableData(:,3)+TableData(:,4)));
xlabel('Distribution of Pearson Correlation(CS)');
ylabel('Number of neurons');
set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
box off;


FigureIndex=FigureIndex+1;
FigureStartNum=FigureStartNum+1;

%Example plot
if PlotExample
   figtitlestr{FigureIndex}='PearsonCorrPassiveActiveExample';
fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[110,100, 1000,800],'Name',figtitlestr{FigureIndex});

 MarkerSize=300;
 scatter(ResponsePassive_Example(5:8),ResponseActive_Example(5:8),MarkerSize*ones(1,4),'o','MarkerFaceColor',StableSceneColor,'MarkerEdgeColor','w');
 hold on
  scatter(ResponsePassive_Example(1:4),ResponseActive_Example(1:4),MarkerSize*ones(1,4),MarkerSize*ones(1,4),'o','MarkerFaceColor',FlexibleSceneColor,'MarkerEdgeColor','w');
 axis equal;
 %Add diagnal line
 x=min([ResponseActive_Example,ResponsePassive_Example])-0.1:max([ResponseActive_Example,ResponsePassive_Example]+1);
 
 plot(x,x,'--k','LineWidth',2);
 xlim([min(x),max(x)]);
 ylim([min(x),max(x)])
 
 set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
box off;
xlabel('Over spon response of passive scene');
ylabel('Over spon response of actice scene');
    
 FigureIndex=FigureIndex+1;
FigureStartNum=FigureStartNum+1;   
end




%Scatter plot between active and passive 
figtitlestr{FigureIndex}='PropOfSig';
fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[100,100, 600,600],'Name',figtitlestr{FigureIndex});
%First, compare with itself% Active stable v.s. Flexible
%Overall 
subplot(2,2,1); %Overall comparison between flexible and stable in active condition
MarkerSize=200;
%None
scatter(PropOfSig_Stable_Active(SSI_Sig_None),PropOfSig_Stable_Active(SSI_Sig_None),MarkerSize*ones(1,sum(SSI_Sig_None)),'o','MarkerFaceColor', NoNeColor,'MarkerEdgeColor','w');
hold on
%Stable only
scatter(PropOfSig_Stable_Active(SSI_Sig_StableOnly),PropOfSig_Flexible_Active(SSI_Sig_StableOnly),MarkerSize*size(PropOfSig_Stable_Active(SSI_Sig_StableOnly),1),'o','MarkerFaceColor',StableSceneColor,'MarkerEdgeColor','w');

%Flexible only
scatter(PropOfSig_Stable_Active(SSI_Sig_FlexibleOnly),PropOfSig_Flexible_Active(SSI_Sig_FlexibleOnly),MarkerSize*size(PropOfSig_Stable_Active(SSI_Sig_FlexibleOnly),1),'o','MarkerFaceColor',FlexibleSceneColor,'MarkerEdgeColor','w');

%Both
scatter(PropOfSig_Stable_Active(SSI_Sig_Both),PropOfSig_Flexible_Active(SSI_Sig_Both),MarkerSize*size(PropOfSig_Flexible_Active(SSI_Sig_Both),1),'^','MarkerFaceColor', '#8F6BD6','MarkerEdgeColor','w');


%Add diagnal line
 x=linspace(0,1,20);
 
 plot(x,x,'--k','LineWidth',2);
 xlim([min(x),max(x)]);
 ylim([min(x),max(x)]);
 
 set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
box off;
xlabel('PropOfSig_Stable');
ylabel('PropOfSig_Flexible');
title('Active');

%Overall Passive
subplot(2,2,2); %Overall comparison between flexible and stable in active condition

%None
scatter(PropOfSig_Stable_Passive(SSI_Sig_None_Compare),PropOfSig_Stable_Passive(SSI_Sig_None_Compare),MarkerSize*ones(1,sum(SSI_Sig_None_Compare)),'o','MarkerFaceColor', NoNeColor,'MarkerEdgeColor','w');
hold on
%Stable only
scatter(PropOfSig_Stable_Passive(SSI_Sig_StableOnly_Compare),PropOfSig_Flexible_Passive(SSI_Sig_StableOnly_Compare),MarkerSize*size(PropOfSig_Stable_Passive(SSI_Sig_StableOnly_Compare),1),'o','MarkerFaceColor',StableSceneColor,'MarkerEdgeColor','w');

%Flexible only
scatter(PropOfSig_Stable_Passive(SSI_Sig_FlexibleOnly_Compare),PropOfSig_Flexible_Passive(SSI_Sig_FlexibleOnly_Compare),MarkerSize*size(PropOfSig_Stable_Passive(SSI_Sig_FlexibleOnly_Compare),1),'o','MarkerFaceColor',FlexibleSceneColor,'MarkerEdgeColor','w');

%Both
scatter(PropOfSig_Stable_Passive(SSI_Sig_Both_Compare),PropOfSig_Flexible_Passive(SSI_Sig_Both_Compare),MarkerSize*size(PropOfSig_Flexible_Passive(SSI_Sig_Both_Compare),1),'^','MarkerFaceColor', '#8F6BD6','MarkerEdgeColor','w');


%Add diagnal line
 x=linspace(0,1,20);
 
 plot(x,x,'--k','LineWidth',2);
 xlim([min(x),max(x)]);
 ylim([min(x),max(x)]);
 
 set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
box off;
xlabel('PropOfSig_Stable');
ylabel('PropOfSig_Flexible');
title('Passive');



    
%Stable Scenes
subplot(2,2,3); %Overall comparison between active and passive


%None
scatter(PropOfSig_Stable_ActiveCompare_None,PropOfSig_Stable_PassiveCompare_None,MarkerSize*size(PropOfSig_Stable_ActiveCompare_None,1),'o','MarkerFaceColor', NoNeColor,'MarkerEdgeColor','w');
hold on
%Stable only
scatter(PropOfSig_Stable_ActiveCompare_StableOnly,PropOfSig_Stable_PassiveCompare_StableOnly,MarkerSize*size(PropOfSig_Stable_PassiveCompare_StableOnly,1),'o','MarkerFaceColor',StableSceneColor,'MarkerEdgeColor','w');

%Flexible only
scatter(PropOfSig_Stable_ActiveCompare_FlexibleOnly,PropOfSig_Stable_PassiveCompare_FlexibleOnly,MarkerSize*size(PropOfSig_Stable_ActiveCompare_FlexibleOnly,1),'o','MarkerFaceColor',FlexibleSceneColor,'MarkerEdgeColor','w');

%Both
scatter(PropOfSig_Stable_ActiveCompare_Both,PropOfSig_Stable_PassiveCompare_Both,MarkerSize*size(PropOfSig_Stable_ActiveCompare_Both,1),'^','MarkerFaceColor', '#8F6BD6','MarkerEdgeColor','w');


%Add diagnal line
 x=linspace(0,1,20);
 
 plot(x,x,'--k','LineWidth',2);
 xlim([min(x),max(x)]);
 ylim([min(x),max(x)]);
 
  set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
box off;
xlabel('PropOfSig Active');
ylabel('PropOfSig Passive');


 title('Stable Scenes');
 
 
 %Flexible Scenes
subplot(2,2,4); %Overall comparison between active and passive


%None
scatter(PropOfSig_Flexible_ActiveCompare_None,PropOfSig_Flexible_PassiveCompare_None,MarkerSize*size(PropOfSig_Flexible_ActiveCompare_None,1),'o','MarkerFaceColor', NoNeColor,'MarkerEdgeColor','w');
hold on
%Stable only
scatter(PropOfSig_Flexible_ActiveCompare_StableOnly,PropOfSig_Flexible_PassiveCompare_StableOnly,MarkerSize*size(PropOfSig_Flexible_PassiveCompare_StableOnly,1),'o','MarkerFaceColor',StableSceneColor,'MarkerEdgeColor','w');

%Flexible only
scatter(PropOfSig_Flexible_ActiveCompare_FlexibleOnly,PropOfSig_Flexible_PassiveCompare_FlexibleOnly,MarkerSize*size(PropOfSig_Flexible_ActiveCompare_FlexibleOnly,1),'o','MarkerFaceColor',FlexibleSceneColor,'MarkerEdgeColor','w');

%Both
scatter(PropOfSig_Flexible_ActiveCompare_Both,PropOfSig_Flexible_PassiveCompare_Both,MarkerSize*size(PropOfSig_Flexible_ActiveCompare_Both,1),'^','MarkerFaceColor', '#8F6BD6','MarkerEdgeColor','w');


%Add diagnal line
 x=linspace(0,1,20);
 
 plot(x,x,'--k','LineWidth',2);
 xlim([min(x),max(x)]);
 ylim([min(x),max(x)]);
 
  set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
box off;
xlabel('PropOfSig Active');
ylabel('PropOfSig Passive');


 title('Flexible Scenes');
 
FigureIndex=FigureIndex+1;
FigureStartNum=FigureStartNum+1;   

%% Distribution of the difference between prop of sig between active and passive condition
figtitlestr{FigureIndex}='DistrDiffPropOfSig';
fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[120,100, 600,1000],'Name',figtitlestr{FigureIndex}); 
 
edge=linspace(-1,1,19);

%First for stable scenes
subplot(4,2,1);
%Stable Only-Stable

[counts,edge]=histcounts(DifferencePropOfSig_StableOnly_Stable,edge);
b=bar(edge(1:end-1),counts,'FaceColor',StableSceneColor,'EdgeColor','k','LineWidth',3);
hold on
MeanVal=nanmean(DifferencePropOfSig_StableOnly_Stable);
[p,h]=ttest(DifferencePropOfSig_StableOnly_Stable,0)
yrange=ylim;

NumberOfStableOnlyHere=length(DifferencePropOfSig_StableOnly_Stable)

plot(MeanVal,yrange(2),'v','MarkerFaceColor','k','MarkerSize',10);
%Add a line

plot([MeanVal,MeanVal],yrange,'--k','LineWidth',3);

xlabel('Difference(Stable Only) ');
ylabel('Num of neurons');
title('Stable Scenes');
set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
box off

subplot(4,2,2);
%Stable Only-Flexible

[counts,edge]=histcounts(DifferencePropOfSig_StableOnly_Flexible,edge);
b=bar(edge(1:end-1),counts,'FaceColor',StableSceneColor,'EdgeColor','k','LineWidth',3);
hold on
MeanVal=nanmean(DifferencePropOfSig_StableOnly_Flexible);
[p,h]=ttest(DifferencePropOfSig_StableOnly_Flexible,0)
yrange=ylim;

plot(MeanVal,yrange(2),'v','MarkerFaceColor','k','MarkerSize',10);
%Add a line

plot([MeanVal,MeanVal],yrange,'--k','LineWidth',3);

xlabel('Difference(Stable Only) ');
ylabel('Num of neurons');
title('Flexible Scenes');
set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
box off

subplot(4,2,3);
%Flexible Only-Stable

[counts,edge]=histcounts(DifferencePropOfSig_FlexibleOnly_Stable,edge);
b=bar(edge(1:end-1),counts,'FaceColor',FlexibleSceneColor,'EdgeColor','k','LineWidth',3);
hold on
MeanVal=nanmean(DifferencePropOfSig_FlexibleOnly_Stable);
[p,h]=ttest(DifferencePropOfSig_FlexibleOnly_Stable,0)

NumberOfFlexibleOnlyHere=length(DifferencePropOfSig_FlexibleOnly_Stable)

yrange=ylim;

plot(MeanVal,yrange(2),'v','MarkerFaceColor','k','MarkerSize',10);
%Add a line

plot([MeanVal,MeanVal],yrange,'--k','LineWidth',3);

xlabel('Difference(Flexible Only) ');
ylabel('Num of neurons');
set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
box off

subplot(4,2,4);
%Flexible Only-Flexible

[counts,edge]=histcounts(DifferencePropOfSig_FlexibleOnly_Flexible,edge);
b=bar(edge(1:end-1),counts,'FaceColor',FlexibleSceneColor,'EdgeColor','k','LineWidth',3);
hold on
MeanVal=nanmean(DifferencePropOfSig_FlexibleOnly_Flexible);
yrange=ylim;
[p,h]=ttest(DifferencePropOfSig_FlexibleOnly_Flexible,0)




plot(MeanVal,yrange(2),'v','MarkerFaceColor','k','MarkerSize',10);
%Add a line

plot([MeanVal,MeanVal],yrange,'--k','LineWidth',3);

xlabel('Difference(Flexible Only) ');
ylabel('Num of neurons');
set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
box off

subplot(4,2,5);
%Both-Stable

[counts,edge]=histcounts(DifferencePropOfSig_Both_Stable,edge);
b=bar(edge(1:end-1),counts,'FaceColor','#8F6BD6','EdgeColor','k','LineWidth',3);
hold on
MeanVal=nanmean(DifferencePropOfSig_Both_Stable);
yrange=ylim;
[p,h]=ttest(DifferencePropOfSig_Both_Stable,0)

NumberOfBothHere=length(DifferencePropOfSig_Both_Stable)

plot(MeanVal,yrange(2),'v','MarkerFaceColor','k','MarkerSize',10);
%Add a line

plot([MeanVal,MeanVal],yrange,'--k','LineWidth',3);

xlabel('Difference(Both) ');
ylabel('Num of neurons');
set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
box off

subplot(4,2,6);
%Both-Flexible

[counts,edge]=histcounts(DifferencePropOfSig_Both_Flexible,edge);
b=bar(edge(1:end-1),counts,'FaceColor','#8F6BD6','EdgeColor','k','LineWidth',3);
hold on
MeanVal=nanmean(DifferencePropOfSig_Both_Flexible);
yrange=ylim;
[p,h]=ttest(DifferencePropOfSig_Both_Flexible,0)

plot(MeanVal,yrange(2),'v','MarkerFaceColor','k','MarkerSize',10);
%Add a line

plot([MeanVal,MeanVal],yrange,'--k','LineWidth',3);

xlabel('Difference(Both) ');
ylabel('Num of neurons');
set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
box off
subplot(4,2,7);
%None

[counts,edge]=histcounts(DifferencePropOfSig_None_Stable,edge);
b=bar(edge(1:end-1),counts,'FaceColor',NoNeColor,'EdgeColor','k','LineWidth',3);
hold on
MeanVal=nanmean(DifferencePropOfSig_None_Stable);
yrange=ylim;

[p,h]=ttest(DifferencePropOfSig_None_Stable,0)

NumberOfNoneHere=length(DifferencePropOfSig_None_Stable)

plot(MeanVal,yrange(2),'v','MarkerFaceColor','k','MarkerSize',10);
%Add a line

plot([MeanVal,MeanVal],yrange,'--k','LineWidth',3);

xlabel('Difference(None) ');
ylabel('Num of neurons');
set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
box off

subplot(4,2,8);
%None

[counts,edge]=histcounts(DifferencePropOfSig_None_Flexible,edge);
b=bar(edge(1:end-1),counts,'FaceColor',NoNeColor,'EdgeColor','k','LineWidth',3);
hold on
MeanVal=nanmean(DifferencePropOfSig_None_Flexible);
yrange=ylim;
[p,h]=ttest(DifferencePropOfSig_None_Flexible,0)

plot(MeanVal,yrange(2),'v','MarkerFaceColor','k','MarkerSize',10);
%Add a line

plot([MeanVal,MeanVal],yrange,'--k','LineWidth',3);

xlabel('Difference(None) ');
ylabel('Num of neurons');
set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');

box off


 
 
FigureIndex=FigureIndex+1;
FigureStartNum=FigureStartNum+1;   

%scatter(PropOfSig_Stable_Active,PropOfSig_Flexible_Active,'or');
%Stable only-Stable


%DifferencePropOfSig_StableOnly_Stable=PropOfSig_Stable_ActiveCompare_StableOnly-PropOfSig_Stable_PassiveCompare_StableOnly;
% DifferencePropOfSig_StableOnly_Flexibe=PropOfSig_Flexible_ActiveCompare_StableOnly-PropOfSig_Flexible_PassiveCompare_StableOnly;

%{
DifferencePropOfSig_Stable=PropOfSig_Stable_ActiveCompare-PropOfSig_Stable_Passive;
 DifferencePropOfSig_Flexibe=PropOfSig_Flexible_ActiveCompare-PropOfSig_Flexible_Passive;
%}
%{
%%Figure5: SSI Scene Responses for each type

figtitlestr{FigureIndex}='RawSceneSSI_StableOnly';
fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[90,100, 500,500],'Name',figtitlestr{FigureIndex});

%rawResponse_Flexible=PSTHSceneOnFlexible_MeanMat(FlexDominace,:);
%rawResponse_Stable=PSTHSceneOnStable_MeanMat(FlexDominace,:);


for i=1:sum(SSI_Sig_StableOnly)
    subplot(ceil(sqrt(sum(SSI_Sig_StableOnly))),ceil(sqrt(sum(SSI_Sig_StableOnly))),i);
    plot(rsquare_psthMat_StableSig_Stable(i,:),'-r');
    hold on
    plot(rsquare_psthMat_StableSig_Flexible(i,:),'-b');
end


FigureIndex=FigureIndex+1;
FigureStartNum=FigureStartNum+1; 

%%Figure6: Raw Scene Responses for each type

figtitlestr{FigureIndex}='RawSceneSSI_FlexibleOnly';
fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[90,100, 500,500],'Name',figtitlestr{FigureIndex});

for i=1:sum(SSI_Sig_FlexibleOnly)
    subplot(ceil(sqrt(sum(SSI_Sig_FlexibleOnly))),ceil(sqrt(sum(SSI_Sig_FlexibleOnly))),i);
    plot(rsquare_psthMat_FlexibleSig_Stable(i,:),'-r');
    hold on
    plot(rsquare_psthMat_FlexibleSig_Flexible(i,:),'-b');
end


FigureIndex=FigureIndex+1;
FigureStartNum=FigureStartNum+1; 

%%Figure7: Raw Scene Responses for each type

figtitlestr{FigureIndex}='RawSceneRepsonses_Both';
fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[90,100, 500,500],'Name',figtitlestr{FigureIndex});

%rawResponse_Flexible=PSTHSceneOnFlexible_MeanMat(NoDominance,:);
%rawResponse_Stable=PSTHSceneOnStable_MeanMat(NoDominance,:);

for i=1:sum(SSI_Sig_Both)
    subplot(ceil(sqrt(sum(SSI_Sig_Both))),ceil(sqrt(sum(SSI_Sig_Both))),i);
    plot(rsquare_psthMat_BothSig_Stable(i,:),'-r');
    hold on
    plot(rsquare_psthMat_BothSig_Flexible(i,:),'-b');
end


FigureIndex=FigureIndex+1;
FigureStartNum=FigureStartNum+1; 

%}

%NumOfSignificantScene500
%%
%{
%%Figure4: Scene Selective Index for each type
figtitlestr{FigureIndex}='SceneSelectiveIndex_EachTypeOfNeuron';
fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[110,100, 500,500],'Name',figtitlestr{FigureIndex});

%Flexible Dominance
subplot(2,2,1);

ValueForHistogram=SI_SceneOn_Fle_FleDo;
pVal=sigSceneOn_Fle_FleDo;

Pro_Fle_FleDo=sum(pVal<0.05)/length(ValueForHistogram);
Mean_Fle_FleDo_Whole=nanmean(ValueForHistogram);
Mean_Fle_FleDo_Sig=nanmean(ValueForHistogram(pVal<0.05));
Mean_Fle_FleDo_NonSig=nanmean(ValueForHistogram(~(pVal<0.05)));


SigVal=ValueForHistogram(pVal<0.05);
NonSigVal=ValueForHistogram(~(pVal<0.05));
%SepBin=[linspace(min(ValueForHistogram)-1,0,10),linspace(0.1,max(ValueForHistogram)+1,10)];
SepBin=linspace(0,1,20);

Val_Sig_Sep=histcounts(SigVal,SepBin);
NonVal_Sig_Sep=histcounts(NonSigVal,SepBin);
BinCenter=SepBin(1:end-1)+diff(SepBin);

b=bar([BinCenter',BinCenter'],[Val_Sig_Sep',NonVal_Sig_Sep'],'stacked');
b(1).FaceColor=[0,0,0];
b(2).FaceColor=[1,1,1];
title('Flexible Dominance');
xlabel('SI for Flexible Scene');
ylabel('Firing Rate(Hz)');
set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');

%Distribition of SI on stable scene
subplot(2,2,2);
ValueForHistogram=SI_SceneOn_Sta_StaDo;
pVal=sigSceneOn_Sta_StaDo;
Pro_Sta_StaDo=sum(pVal<0.05)/length(ValueForHistogram);
Mean_Sta_StaDo_Whole=nanmean(ValueForHistogram);
Mean_Sta_StaDo_Sig=nanmean(ValueForHistogram(pVal<0.05));
Mean_Sta_StaDo_NonSig=nanmean(ValueForHistogram(~(pVal<0.05)));

SigVal=ValueForHistogram(pVal<0.05);
NonSigVal=ValueForHistogram(~(pVal<0.05));
%SepBin=[linspace(min(ValueForHistogram)-1,0,10),linspace(0.1,max(ValueForHistogram)+1,10)];
SepBin=linspace(0,1,20);

Val_Sig_Sep=histcounts(SigVal,SepBin);
NonVal_Sig_Sep=histcounts(NonSigVal,SepBin);
BinCenter=SepBin(1:end-1)+diff(SepBin);

b=bar([BinCenter',BinCenter'],[Val_Sig_Sep',NonVal_Sig_Sep'],'stacked');
b(1).FaceColor=[0,0,0];
b(2).FaceColor=[1,1,1];
title('Stable Dominance');
xlabel('SI For Stable Scene');
ylabel('Firing Rate(Hz)');
set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');

%sigSceneOn_Fle

subplot(2,2,3);
plot(SI_SceneOn_Fle_NoDo,SI_SceneOn_Sta_NoDo,'or','MarkerSize',10);
hold on
%plot(SISceneOn_Fle(~EitherOneSig),SISceneOn_Sta(~EitherOneSig),'or','MarkerSize',10);
axis equal;
hold on 
x=0:0.1:1;
y=x;
plot(x,y,'--k');

xlabel('SI SceneOn Flexible');
ylabel('SI SceneOn Stable');
title('No Dominance');

set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');



FigureIndex=FigureIndex+1;
FigureStartNum=FigureStartNum+1; 

%}


%%
%{
%%Figure5: Raw Scene Responses for each type

figtitlestr{FigureIndex}='RawSceneRepsonses_FlexibleDominance';
fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[90,100, 500,500],'Name',figtitlestr{FigureIndex});

rawResponse_Flexible=PSTHSceneOnFlexible_MeanMat(FlexDominace,:);
rawResponse_Stable=PSTHSceneOnStable_MeanMat(FlexDominace,:);

for i=1:sum(FlexDominace)
    subplot(ceil(sqrt(sum(FlexDominace))),ceil(sqrt(sum(FlexDominace))),i);
    plot(rawResponse_Flexible(i,:),'-r');
    hold on
    plot(rawResponse_Stable(i,:),'-b');
end


FigureIndex=FigureIndex+1;
FigureStartNum=FigureStartNum+1; 

%%Figure6: Raw Scene Responses for each type

figtitlestr{FigureIndex}='RawSceneRepsonses_StableDominance';
fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[90,100, 500,500],'Name',figtitlestr{FigureIndex});

rawResponse_Flexible=PSTHSceneOnFlexible_MeanMat(StableDominance,:);
rawResponse_Stable=PSTHSceneOnStable_MeanMat(StableDominance,:);

for i=1:sum(StableDominance)
    subplot(ceil(sqrt(sum(StableDominance))),ceil(sqrt(sum(StableDominance))),i);
    plot(rawResponse_Flexible(i,:),'-r');
    hold on
    plot(rawResponse_Stable(i,:),'-b');
end


FigureIndex=FigureIndex+1;
FigureStartNum=FigureStartNum+1; 

%%Figure7: Raw Scene Responses for each type

figtitlestr{FigureIndex}='RawSceneRepsonses_NoDominance';
fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[90,100, 500,500],'Name',figtitlestr{FigureIndex});

rawResponse_Flexible=PSTHSceneOnFlexible_MeanMat(NoDominance,:);
rawResponse_Stable=PSTHSceneOnStable_MeanMat(NoDominance,:);

for i=1:sum(NoDominance)
    subplot(ceil(sqrt(sum(NoDominance))),ceil(sqrt(sum(NoDominance))),i);
    plot(rawResponse_Flexible(i,:),'-r');
    hold on
    plot(rawResponse_Stable(i,:),'-b');
end


FigureIndex=FigureIndex+1;
FigureStartNum=FigureStartNum+1; 

%}
end

%{

SI_SceneOn_Fle_FleDo=SISceneOn_Fle(FlexDominace);
 SI_SceneOn_Sta_FleDo=SISceneOn_Sta(FlexDominace);
 
 SI_SceneOn_Fle_StaDo=SISceneOn_Fle(StableDominance);
 SI_SceneOn_Sta_StaDo=SISceneOn_Sta(StableDominance);
 
 SI_SceneOn_Fle_NoDo=SISceneOn_Fle(NoDominance);
 SI_SceneOn_Sta_nODo=SISceneOn_Sta(NoDominance);
%}

end
function stats=corrcoef_out(x,y)

[mat_ContextVSScramble,pmat]=corrcoef(x,y);
stats(1)=mat_ContextVSScramble(1,2);
stats(2)=pmat(1,2);

end