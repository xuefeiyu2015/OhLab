function  SceneOnResponseAcrossModes(Data, StartTrial, EndTrial,ShowFigureFlag,OutputFlag);
%Function to analyze neural responses during exploratory behavior for both
%fix loc series mode and rand loc series mode
%by Xuefei Yu 08222020
ECode;
TaskType=Data.TaskType;
TaskCode=Data.TaskCode;
TrialNum=EndTrial-StartTrial+1;
Selection=StartTrial:EndTrial;%Select the trials according to the pannel
%Abnormal check %Especially For Robin081120FORAGEXPLO4.rst
CheckCode=SCENE_ON;%This code should appear once in each trial
AbnormalTrials=AbnormalTrialCheck(Data,CheckCode);
Selection=Selection(~ismember(Selection,AbnormalTrials));

CheckCode=SCENE_OFF;%This code should appear once in each trial
AbnormalTrials=AbnormalTrialCheck(Data,CheckCode);
Selection=Selection(~ismember(Selection,AbnormalTrials));


%User defined parameter for the neural/behavioral analysis
TargetOnOffset=50;%in ms, the time selection bin after target onset; 
StartDelay=0;%in ms, the time selection bin after fixation offset;

%For segmentation of  the neural responses 
%{
SegEvent=[SCENE_ON,FIXOFF];%Event Flow of interest
SegEventEnd=[SHOWFIXCD,NaN];%Event end
%SegEventEnd=[FIXOFF];%Event end %Will change automatically
SegEventStr={'SceneOn','FixOff'};

PlotStartOffset=[500,1200];%in ms
PlotStopOffset=[100,200];%in ms

EventTypeFlag=[1,1];%1:Number same with the trial number, 2:Number same with the target number
%}
SegEvent=[SCENE_ON];%Event Flow of interest
SegEventEnd=[SHOWFIXCD];%Event end
%SegEventEnd=[FIXOFF];%Event end %Will change automatically
SegEventStr={'SceneOn'};

PlotStartOffset=[600];%in ms
PlotStopOffset=[600];%in ms

EventTypeFlag=[1];%1:Number same with the trial number, 2:Number same with the target number

RasterOn=0;%Show raster, but it will be very slow to do so

%Scene Type
SceneTypeSep=1;
StableLocSeqScene=1010:1017;
FlexibleLocSeqScene=1000:1007;
IntermediateStableScene=[1054:1057,1060:1063];

%Parameter to define the sliding window to measure firing rates
StepSize=5;%in ms
BinWidth=15;%in ms

%Parameter for raster plot
Barlength=1;
LineWidth=1;
%Spon Event
%
SponEvent=[SCENE_OFF];
SponPreTimeBuffer=-500;
SponEventEnd=NaN;
SponStartOffset=0;
SponStopOffset=1000;
%}
%{
SponEvent=[SCENE_ON];
SponPreTimeBuffer=200;
SponEventEnd=NaN;
SponStartOffset=0;
SponStopOffset=200;
%}
%ResponseLatencyThreshold=3;

%Some user defined parameter for visualization
%Color Definition for each scenes
Color={'r','g','b','m',[0.4,0.4,0],[0.28,0.78,0.89],'k',[0.6,0.1,1],[0.6,0.1,0.4],[0.3,0.5,1],[0.6,0.5,0.1],[0,0.8,0.4]};
ColorForEachMode=[0.83 0.14 0.14
             1.00 0.54 0.00
             0.47 0.25 0.80
             0.25 0.80 0.54];
AssociateBlock=[10,11,0,1];
Blockstr={'FixLoc,FixSceneSeq','FixLoc,RandSceneSeq','RandLoc,FixSceneSeq','RandLoc,RandSceneSeq'};
InterMediate_FlexibleColor=[0.2627    0.5137    0.3451];
InterMediate_StableColor=[ 0.8980    0.7373    0.4941];

%BlockSepTrial=unique(SceneLib(:,1));
BlockType=Data.BlockType(:,:);

BlockType=BlockType(Selection,2);

%Get File Information
num = regexp(Data.FileName, '\d+', 'match');
RecordDateOriginal=cell2mat(num(1));
RecordDate=datetime(RecordDateOriginal,'InputFormat','MMddyy');

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
%
% This analysis is constrained on the random sequence mode
BlockType(BlockType==10)=NaN;
BlockType(BlockType==0)=NaN;
%}
%BlockType(BlockType==0)=NaN;
%BlockType(BlockType==1)=NaN;


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
%Get the information about the extra FP
%
%Recalibrate eyes according to the extra FP
ParaLib=Data.ParaLib;

%Change period automatically

StartCheck=FindOutTime(EventChannel,EventBin,SegEvent(1));
EndCheck=FindOutTime(EventChannel,EventBin,SegEventEnd(1));
EndCheck=EndCheck(:,1);
DiffCheck=EndCheck-StartCheck;
if sum(DiffCheck<100)>ceil(length(StartCheck)/2)
    SegEventEnd(1)=[FIXOFF];%Event end;change from showfix cd to fix off for contigency
    disp('Change to FixOff for FPOn First Mode');
    
end


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

%}

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
 
 TargetNum=TargetNum(AdvanceInTargetTrial);
 BlockType=BlockType(AdvanceInTargetTrial);
 
 [UniqueBlockType,orderindex,reprojectorderindex]=unique(BlockType,'row');
 RemoveUnique=isnan(UniqueBlockType);
 UniqueBlockType=UniqueBlockType(~RemoveUnique);
 reprojectorderindex=reprojectorderindex(~RemoveUnique);
%NumberOfBlocks=max(reprojectorderindex);
NumberOfBlocks=length(UniqueBlockType);
 
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
    
    

%Segregate spike data from Fixation onset to target onset
SegmentData=EventSegTool(SegEvent,SegEventEnd,PlotStartOffset,PlotStopOffset,EventChannel,EventTimeChannel,SpikeTime,StepSize,BinWidth);

%Retrive data from the field
%Total average for each event
MeanFR=SegmentData.MeanFR; %[Mean,sem]%Averall mean and sem
MeanFRRaw=SegmentData.MeanFRRaw;%=AverageFiringRate_Event;

%PSTH(Spike Density Function) for each event
%MeanPSTH=SegmentData.MeanPSTH;%=[PSTH_Event_Time',PSTH_Event'];
%MeanPSTH_Time=MeanPSTH(:,1);
%MeanPSTH_FR=MeanPSTH(:,2);

MeanPSTH=SegmentData.MeanPSTHEach;
MeanPSTH_Time=cellfun(@(x) x{1},MeanPSTH(:,1),'uniform',0);%Only need responeses from the first FP off
MeanPSTH_FR=cellfun(@(x) x{1},MeanPSTH(:,2),'uniform',0);%Only need responeses from the first FP off


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

%Rearrange block type into lines
BlockType_Line=repmat(BlockType,1,max(TargetNum));
BlockType_Line=reshape(BlockType_Line,numel(BlockType_Line),1);
BlockType_Line=BlockType_Line(NoNNaNSelector);

%FR for spon
EventTimeChannelForSpon=EventTimeChannel-SponPreTimeBuffer;

SegmentDataSpon=EventSegTool(SponEvent,SponEventEnd,SponStartOffset,SponStopOffset,EventChannel,EventTimeChannelForSpon,SpikeTime,StepSize,BinWidth);

Mean_Spon=SegmentDataSpon.MeanFR;
Mean_Spon_ForEachTrial=SegmentDataSpon.MeanFRRaw{1};

Mean_Spon_ForEachTrial_Line=repmat(Mean_Spon_ForEachTrial,1,max(TargetNum));
Mean_Spon_ForEachTrial_Line=reshape(Mean_Spon_ForEachTrial_Line,numel(Mean_Spon_ForEachTrial_Line),1);
Mean_Spon_ForEachTrial_Line=Mean_Spon_ForEachTrial_Line(NoNNaNSelector);




% Loop over the blocks
for i=1:NumberOfBlocks
    
  
    %Organize the raster, spikecount PSTH, spike density function PSTH
    for j=1:length(AlignedRaster)
        if EventTypeFlag(j)==1 %Trial event
            currBlocks=BlockType==UniqueBlockType(i);
            spon_tmp{i}=Mean_Spon_ForEachTrial(currBlocks);
            
        elseif EventTypeFlag(j)==2 %Target event
            currBlocks=BlockType_Line==UniqueBlockType(i);
             spon_tmp{i}=Mean_Spon_ForEachTrial_Line(currBlocks);
        end
        
        %Raster
        Raster_currEvent_Block{i,j}=AlignedRaster{j}(currBlocks,:);
        
        
        %Spike count PSTH
        PSTH_Event_SpikeCount_currEvent{i,j}=PSTH_Event_SpikeCount{j}(currBlocks,:);
        %Spike density function PSTH
        MeanPSTH_FR_currEvent{i,j}=MeanPSTH_FR{j}(currBlocks,:);%-spon_tmp{i}; %Minus spon;
        
        %Timing 
        PSTH_Time_currEvent{i,j}=nanmean(MeanPSTH_Time{j}(currBlocks,:),1);
        RelativeTimeMarker_currEvent{i,j}=RelativeTimeMarker{j}(currBlocks,:);
        
        %{
        StartLoc=1;
        Barlength=1;
        LineWidth=1;
        SpikeRaterPlot(Raster_currEvent_Block,StartLoc,Barlength,[1,0,0],LineWidth);
        %}
    end
   
 
end
%Get scene related responses and time being
 Response=MeanPSTH_FR{1};%-Mean_Spon_ForEachTrial;
 TimeForAll=nanmean(MeanPSTH_Time{1},1);
 AverageInterval=nanmean(RelativeTimeMarker{1},1);
 IntervalSel=find((TimeForAll>=AverageInterval(1))&(TimeForAll<=AverageInterval(2)));
%Loop over all scenes

Intermediate_Flexible_Scene_Mean_Whole=[];
    Intermediate_Flexible_Scene_Sem_Whole=[];
    Intermediate_Stable_Scene_Mean_Whole=[];
    Intermediate_Stable_Scene_Sem_Whole=[];
     Intermediate_Stable_Resp_Mean_Each=[];
    Intermediate_Flexible_Resp_Mean_Each=[];
    
UniqueScenes=unique(SceneName);

 DuralSceneIndexS=1;
 DuralSceneIndexF=1;
for i=1:length(UniqueScenes)
    Sel=(SceneName==UniqueScenes(i))&ismember(BlockType,[1,11]);
    Resp_tmp=nanmean(Response(Sel,IntervalSel),2);
    Scene_Resp_Mean_Whole(i)=nanmean(Resp_tmp);
    Scene_Resp_Sem_Whole(i)=nanstd(Resp_tmp)/sqrt(sum(~isnan(Resp_tmp)));
    
    
    Resp_tmp_dynamics=Response(Sel,:);
    Scene_Resp_Mean_Each{i}=nanmean( Resp_tmp_dynamics,1);
    Scene_Resp_Sem_Each{i}=nanstd( Resp_tmp_dynamics,[],1)./sqrt(sum(~isnan( Resp_tmp_dynamics),1));
    
    %For Dual scene
    
   
   
    
    if ismember(UniqueScenes(i),IntermediateStableScene)
         Sel_Intermediate_Stable=(SceneName==UniqueScenes(i))& BlockType==11;
         Resp_tmp=nanmean(Response(Sel_Intermediate_Stable,IntervalSel),2);
         Intermediate_Stable_Scene_Mean_Whole(DuralSceneIndexS)=nanmean(Resp_tmp);
         Intermediate_Stable_Scene_Sem_Whole(DuralSceneIndexS)=nanstd(Resp_tmp)/sqrt(sum(~isnan(Resp_tmp)));
         
         
         Resp_Dynamics_tmpS=Response(Sel_Intermediate_Stable,:);
         Intermediate_Stable_Resp_Mean_Each{DuralSceneIndexS}=nanmean( Resp_Dynamics_tmpS,1)';
         
         Intermediate_Stable_Resp_Sem_Each{DuralSceneIndexS}=nanstd( Resp_Dynamics_tmpS,[],1)./sqrt(sum(~isnan( Resp_Dynamics_tmpS),1));
         
         DuralSceneIndexS=DuralSceneIndexS+1;
  
         Sel_Intermediate_Flexible=(SceneName==UniqueScenes(i))& BlockType==1;
         Resp_tmp=nanmean(Response(Sel_Intermediate_Flexible,IntervalSel),2);
         Intermediate_Flexible_Scene_Mean_Whole(DuralSceneIndexF)=nanmean(Resp_tmp);
         Intermediate_Flexible_Scene_Sem_Whole(DuralSceneIndexF)=nanstd(Resp_tmp)/sqrt(sum(~isnan(Resp_tmp)));
         
         Resp_Dynamics_tmpF=Response(Sel_Intermediate_Flexible,:);
         Intermediate_Flexible_Resp_Mean_Each{DuralSceneIndexF}=nanmean( Resp_Dynamics_tmpF,1)';
         
         Intermediate_Flexible_Resp_Sem_Each{DuralSceneIndexF}=nanstd( Resp_Dynamics_tmpF,[],1)./sqrt(sum(~isnan( Resp_Dynamics_tmpF),1));
        DuralSceneIndexF=DuralSceneIndexF+1;
        
    end
    
    
    
    
    
end



if SceneTypeSep
    %Plot scene onset responses for each type of scene
   
    
    StableSceneSel=ismember(SceneName,StableLocSeqScene)&ismember(BlockType,[1,11]);
    FlexibleSceneSel=ismember(SceneName,FlexibleLocSeqScene)&ismember(BlockType,[1,11]);
    IntermediateSceneSel=ismember(SceneName,IntermediateStableScene)&ismember(BlockType,[1,11]);%Here only consider for random scene sequence condition
    
    StableResp=Response(StableSceneSel,:);
    FlexibleResp=Response(FlexibleSceneSel,:);
    IntermedaiteResp=Response(IntermediateSceneSel,:);
    
    IntermedaiteResp_Flexible=Response(IntermediateSceneSel&BlockType==1,:);
    IntermedaiteResp_Stable=Response(IntermediateSceneSel&BlockType==11,:);
    
    StableResp_Mean=nanmean(StableResp,1);
    FlexibleResp_Mean=nanmean(FlexibleResp,1);
    IntermedaiteResp_Mean=nanmean(IntermedaiteResp,1);
    
    IntermedaiteResp_Flexible_Mean=nanmean(IntermedaiteResp_Flexible,1);
    IntermedaiteResp_Stable_Mean=nanmean(IntermedaiteResp_Stable,1);
    
    StableResp_Sem=nanstd(StableResp,[],1)./sqrt(sum(~isnan(StableResp),1));
    FlexibleResp_Sem=nanstd(FlexibleResp,1)./sqrt(sum(~isnan(FlexibleResp),1));
    IntermedaiteResp_Sem=nanstd(IntermedaiteResp,1)./sqrt(sum(~isnan(IntermedaiteResp),1));
    
    IntermedaiteResp_Flexible_Sem=nanstd(IntermedaiteResp_Flexible,[],1)./sqrt(sum(~isnan(IntermedaiteResp_Flexible),1));
    IntermedaiteResp_Stable_Sem=nanstd(IntermedaiteResp_Stable,[],1)./sqrt(sum(~isnan(IntermedaiteResp_Stable),1));
    
    
    
    StableResp_Mean_Whole=nanmean(nanmean(StableResp(:,IntervalSel),2));
    FlexibleResp_Mean_Whole=nanmean(nanmean(FlexibleResp(:,IntervalSel),2));
    IntermedaiteResp_Mean_Whole=nanmean(nanmean(IntermedaiteResp(:,IntervalSel),2));
    
    IntermedaiteResp_Flexible_Mean_Whole=nanmean(nanmean(IntermedaiteResp_Flexible(:,IntervalSel),2));
    IntermedaiteResp_Stable_Mean_Whole=nanmean(nanmean(IntermedaiteResp_Stable(:,IntervalSel),2));
    
    
    StableResp_Sem_Whole=nanstd(nanmean(StableResp(:,IntervalSel),2))/sqrt(length(nanmean(StableResp(:,IntervalSel),2)));
    FlexibleResp_Sem_Whole=nanstd(nanmean(FlexibleResp(:,IntervalSel),2))/sqrt(length(nanmean(FlexibleResp(:,IntervalSel),2)));
    IntermedaiteResp_Sem_Whole=nanstd(nanmean(IntermedaiteResp(:,IntervalSel),2))/sqrt(length(nanmean(IntermedaiteResp(:,IntervalSel),2)));
    
    IntermedaiteResp_Flexible_Sem_Whole=nanstd(nanmean(IntermedaiteResp_Flexible(:,IntervalSel),2))/sqrt(length(nanmean(IntermedaiteResp_Flexible(:,IntervalSel),2)));
    IntermedaiteResp_Stable_Sem_Whole=nanstd(nanmean(IntermedaiteResp_Stable(:,IntervalSel),2))/sqrt(length(nanmean(IntermedaiteResp_Flexible(:,IntervalSel),2)));
    
    [h_Stable_Flexible,p_Stable_Flexible]=ttest2(nanmean(StableResp(:,IntervalSel),2),nanmean(FlexibleResp(:,IntervalSel),2));
    [h_Stable_Flexible_Dual,p_Stable_Flexible_Dual]=ttest2(nanmean(IntermedaiteResp_Flexible(:,IntervalSel),2),nanmean(IntermedaiteResp_Stable(:,IntervalSel),2));
    
    %Scene Selective Index, from ?Wirth, Sylvia et.al 2009 Neuron
    % SI =(n-sum(lamda(i)-lamda(max)))/(n-1);
    %Stable Scene
    Stable_Unique_Sel=ismember(UniqueScenes,StableLocSeqScene);
    N_Stable_Unique=sum(Stable_Unique_Sel);
    Resp_Stable_Unique_Mean=Scene_Resp_Mean_Whole(Stable_Unique_Sel);
    SI_Stable=(N_Stable_Unique-sum(Resp_Stable_Unique_Mean/max(Resp_Stable_Unique_Mean)))/(N_Stable_Unique-1);
    %Flexible Scene
    Flexible_Unique_Sel=ismember(UniqueScenes,FlexibleLocSeqScene);
    N_Flexible_Unique=sum(Flexible_Unique_Sel);
    Resp_Flexible_Unique_Mean=Scene_Resp_Mean_Whole(Flexible_Unique_Sel);
    SI_Flexible=(N_Flexible_Unique-sum(Resp_Flexible_Unique_Mean/max(Resp_Flexible_Unique_Mean)))/(N_Flexible_Unique-1);
    %Dural Scene
    %Whole
    Intermediate_Unique_Sel=ismember(UniqueScenes,IntermediateStableScene);
    N_Intermediate_Unique=sum(Intermediate_Unique_Sel);
    Resp_Intermediate_Unique_Mean=Scene_Resp_Mean_Whole(Intermediate_Unique_Sel);
    SI_Intermediate=(N_Intermediate_Unique-sum(Resp_Intermediate_Unique_Mean/max(Resp_Intermediate_Unique_Mean)))/(N_Intermediate_Unique-1);
    
    
    %Separate
    N_IntermediateStable_Unique=NaN;
    N_IntermediateFlexible_Unique=NaN;
    SI_IntermediateStable=NaN;
    SI_IntermediateFlexible=NaN;
    
    if ~isempty(Intermediate_Flexible_Scene_Mean_Whole) && ~isempty(Intermediate_Stable_Scene_Mean_Whole) 
        %Intermediate Stable Scene
        N_IntermediateStable_Unique=length(Intermediate_Stable_Scene_Mean_Whole);
        SI_IntermediateStable=(N_IntermediateStable_Unique-sum(Intermediate_Stable_Scene_Mean_Whole/max(Intermediate_Stable_Scene_Mean_Whole)))/(N_IntermediateStable_Unique-1);
    %IntermediateFlexible Scene
        N_IntermediateFlexible_Unique=length(Intermediate_Flexible_Scene_Mean_Whole);
        SI_IntermediateFlexible=(N_IntermediateFlexible_Unique-sum(Intermediate_Flexible_Scene_Mean_Whole/max(Intermediate_Flexible_Scene_Mean_Whole)))/(N_IntermediateFlexible_Unique-1);
 
    end
    
    %Temporal dynamics of SI
    %Stable
    Scene_Resp_Mean_Each_Stable=cell2mat(Scene_Resp_Mean_Each(Stable_Unique_Sel)');
    N_Stable_Each=sum(~isnan(Scene_Resp_Mean_Each_Stable),1);
    Scene_Resp_Mean_Each_Stable(:,N_Stable_Each<3)=NaN;
    if ~isempty(Scene_Resp_Mean_Each_Stable)
        SI_Stable_SceneEach=(N_Stable_Each-sum(Scene_Resp_Mean_Each_Stable./max(Scene_Resp_Mean_Each_Stable,[],1),1))./(N_Stable_Each-1);
    
    else
         SI_Stable_SceneEach=NaN*ones(1,length(TimeForAll));
    end
        
    
    %Flexible
    Scene_Resp_Mean_Each_Flexible=cell2mat(Scene_Resp_Mean_Each(Flexible_Unique_Sel)');
    N_Flexible_Each=sum(~isnan(Scene_Resp_Mean_Each_Flexible),1);
    Scene_Resp_Mean_Each_Flexible(:,N_Flexible_Each<3)=NaN;
    if ~isempty(Scene_Resp_Mean_Each_Flexible)
        SI_Flexible_SceneEach=(N_Flexible_Each-sum(Scene_Resp_Mean_Each_Flexible./max(Scene_Resp_Mean_Each_Flexible,[],1),1))./(N_Flexible_Each-1);
    else
        SI_Flexible_SceneEach=NaN*ones(1,length(TimeForAll));
    end
        
    %Intermediate
    Scene_Resp_Mean_Each_Intermediate=cell2mat(Scene_Resp_Mean_Each(Intermediate_Unique_Sel)');
    N_Intermediate_Each=sum(~isnan(Scene_Resp_Mean_Each_Intermediate),1);
    Scene_Resp_Mean_Each_Intermediate(:,N_Intermediate_Each<3)=NaN;
    if ~isempty(Scene_Resp_Mean_Each_Intermediate)
     SI_Intermediate_SceneEach=(N_Intermediate_Each-sum(Scene_Resp_Mean_Each_Intermediate./max(Scene_Resp_Mean_Each_Intermediate,[],1),1))./(N_Intermediate_Each-1);
    else
       SI_Intermediate_SceneEach=NaN*ones(1,length(TimeForAll)); 
    end
    %Separate
    %Stable
    Scene_Resp_Mean_Each_IntermediateS=cell2mat(Intermediate_Stable_Resp_Mean_Each)';
    N_IntermediateS_Each=sum(~isnan(Scene_Resp_Mean_Each_IntermediateS),1);
    Scene_Resp_Mean_Each_IntermediateS(:,N_IntermediateS_Each<3)=NaN;
    if ~isempty(Scene_Resp_Mean_Each_IntermediateS)
        SI_IntermediateS_SceneEach=(N_IntermediateS_Each-sum(Scene_Resp_Mean_Each_IntermediateS./max(Scene_Resp_Mean_Each_IntermediateS,[],1),1))./(N_IntermediateS_Each-1);
    else
         SI_IntermediateS_SceneEach=NaN*ones(1,length(TimeForAll)); 
    end
   %Flexible
    Scene_Resp_Mean_Each_IntermediateF=cell2mat(Intermediate_Flexible_Resp_Mean_Each)';
    N_IntermediateF_Each=sum(~isnan(Scene_Resp_Mean_Each_IntermediateF),1);
    Scene_Resp_Mean_Each_IntermediateF(:,N_IntermediateF_Each<3)=NaN;
    if ~isempty(Scene_Resp_Mean_Each_IntermediateF)
        SI_IntermediateF_SceneEach=(N_IntermediateF_Each-sum(Scene_Resp_Mean_Each_IntermediateF./max(Scene_Resp_Mean_Each_IntermediateF,[],1),1))./(N_IntermediateF_Each-1);
    else
        SI_IntermediateF_SceneEach=NaN*ones(1,length(TimeForAll)); 
    end
    
   % Intermediate_Flexible_Resp_Mean_Each
    
 
    NumStable=sum(StableSceneSel);
    NumFlexible=sum(FlexibleSceneSel);
    NumIntermediate=sum(IntermediateSceneSel);
    NumIntermediate_Flexible=sum(IntermediateSceneSel&BlockType==1);
    NumIntermediate_Stable=sum(IntermediateSceneSel&BlockType==11);
    
end


%%
FigureStartNum=110;
FigureIndex=1;
if ShowFigureFlag
   if RasterOn
    figtitlestr{FigureIndex}='Event_Raster_Check';
    fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[50,100, 2000,800],'Name',figtitlestr{FigureIndex});
   
    for i=1:length(AlignedRaster)
        
         StartLoc=1;
        for j=1:NumberOfBlocks
            %First Figure: raster
            
            subplot(1,length(AlignedRaster),i);
            NumTrials=size(Raster_currEvent_Block{j,i},1);
            currBlock=UniqueBlockType(j);
            
            %Sort Event according to the event persist time
            [Val_Sort,SortIndex]=sort(RelativeTimeMarker_currEvent{j,i}(:,2));
            
            currResponse=Raster_currEvent_Block{j,i}(SortIndex,:);
            %Event-on line
            plot([0,0],[StartLoc-1,StartLoc+NumTrials],'--k','Color','#929498','LineWidth',1);
            hold on
            %Spike raster
            SpikeRaterPlot(currResponse,StartLoc,Barlength,ColorForEachMode(AssociateBlock==currBlock,:),LineWidth);
            
            %Add event-end marker
            y_tmp=[StartLoc:StartLoc+NumTrials-1]';
            x_tmp=Val_Sort;
            plot(x_tmp,y_tmp,'>','MarkerEdgeColor','#707276','MarkerFaceColor','#707276','MarkerSize',2);
            
            xlabel('Time from event onset');
            ylabel('Trial #');
            title(SegEventStr{i});
            set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
          %  legend(Blockstr{AssociateBlock==currBlock});
        
            StartLoc=StartLoc+NumTrials;
            %{
            %Second row: Spike-Count PSTH
            subplot(3,length(AlignedRaster),i+length(AlignedRaster));
            %Event-on line
            plot([0,0],[0,max(nanmean(PSTH_Event_SpikeCount_currEvent{i,j})/StepSize*1000)],'--k','Color','#929498','LineWidth',1);
            bar(PSTH_Time_currEvent{i,j},nanmean(PSTH_Event_SpikeCount_currEvent{i,j},1)/StepSize*1000,'FaceColor',ColorForEachMode(AssociateBlock==currBlock,:));
            hold on 
            %Add event-end marker
            plot([nanmean(RelativeTimeMarker_currEvent{j,i}(:,2)), nanmean(RelativeTimeMarker_currEvent{j,i}(:,2))],[0,max(nanmean(PSTH_Event_SpikeCount_currEvent{i,j})/StepSize*1000)],'-k','Color','#707276','LineWidth',3);
             ylabel('Spike Count per s(Hz)');
            set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
            %}
        end
    end
    FigureIndex=FigureIndex+1;
    FigureStartNum=FigureStartNum+1;
  end %End of the if raster on
  
  
 
%%%%%%%%%%%%%%%%%%Next figure: Spike Density Function PSTH
 if ~MultiChannel || (MultiChannel && spk==1)
 figtitlestr{1}='SpikeDensityFunctionPSTH';
 fig{1}=PrepareFigure(110,'w',[100,100, 800,600+MultiChannel*1000],'Name',figtitlestr{1});   
 else
     
     fig{1}=figure(110);
       
 end
 if MultiChannel 
     %Plot each spike on each line
     subplot(length(SpikeChannelWhole),1,spk);
 end
 
   for i=1:length(AlignedRaster)
        
         
        for j=1:NumberOfBlocks
            
            %Spike-Count PSTH
           
            currBlock=UniqueBlockType(j);
            BlockOrder(j)=find(AssociateBlock==currBlock);
            
            subplot(length(AlignedRaster),1,i);
            %Event-on line
            curr_Response=nanmean(MeanPSTH_FR_currEvent{j,i},1);
            curr_Response_err=nanstd(MeanPSTH_FR_currEvent{j,i},[],1)./sqrt(sum(~isnan(MeanPSTH_FR_currEvent{j,i}),1));
            
            h=plot([0,0],[0,max(curr_Response)],'--k','Color','#929498','LineWidth',1);
            hold on
            %plot(PSTH_Time_currEvent{j,i},curr_Response,'-k','Color',ColorForEachMode(AssociateBlock==currBlock,:),'LineWidth',2);
           % bar(PSTH_Time_currEvent{j,i},nanmean(PSTH_Event_SpikeCount_currEvent{j,i},1),'FaceColor',ColorForEachMode(AssociateBlock==currBlock,:));
            shadedErrorBar(PSTH_Time_currEvent{j,i},curr_Response,curr_Response_err,'lineprops',{ColorForEachMode(AssociateBlock==currBlock,:)},'transparent',1,'patchSaturation',0.5);%PSTH for fix
             set(findall(gca, 'Type', 'Line'),'LineWidth',3);
            hold on 
            %Add event-end marker
            h2=plot([min(RelativeTimeMarker_currEvent{j,i}(:,2)), min(RelativeTimeMarker_currEvent{j,i}(:,2))],[0,max(curr_Response)],'--k','Color','#707276','LineWidth',1.5);
            xlim([-PlotStartOffset(i),PlotStopOffset(i)+nanmean(RelativeTimeMarker_currEvent{j,i}(:,2))]);
            
            %Set ylim
            %Set ylim
            Interval=curr_Response(:,PSTH_Time_currEvent{j,i}<PlotStopOffset(i)+nanmean(RelativeTimeMarker_currEvent{j,i}(:,2)) & PlotStopOffset(i)+nanmean(RelativeTimeMarker_currEvent{j,i}(:,2))>-PlotStartOffset(i));
              Interval_max(j)=max( Interval)*1.12;
            
           %Some legend off
            set( get( get( h, 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );
            set( get( get( h2, 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );
           if ~MultiChannel 
            if i==3
             ylabel('Spike Count per s(Hz)');
            end
            if i==length(AlignedRaster)
                
                xlabel('Time from event onset');
            end
             title(SegEventStr{i});
            
            set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
           end
            %}
        end
       %  ylim([0,max(  Interval_max)]);
        if ~MultiChannel 
         ylabel('Firing Rate(Hz)');% relative to spon');
          legend(Blockstr{BlockOrder});
        end
         
          box off;
    end
    FigureIndex=FigureIndex+1;
    FigureStartNum=FigureStartNum+1; 
     
     
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Next figure: PSTH for each type of scenes
 if SceneTypeSep
     if ~MultiChannel || (MultiChannel && spk==1)
     
     figtitlestr{2}='SpikeDensityFunctionFor_Different_TypeOf_Scenes';
     fig{2}=PrepareFigure(111,'w',[120,100, 800,600+MultiChannel*1000],'Name',figtitlestr{2});   
     else
          fig{2}=figure(111);
     end
     if MultiChannel
       %Plot each spike on each line
       subplot(length(SpikeChannelWhole),1,spk);
     else
   
        subplot(1,2,1)
     end
     shadedErrorBar(TimeForAll,StableResp_Mean,StableResp_Sem,'lineprops',{ColorForEachMode(AssociateBlock==11,:)},'transparent',1,'patchSaturation',0.5);
     
     hold on 
     shadedErrorBar(TimeForAll,FlexibleResp_Mean,FlexibleResp_Sem,'lineprops',{ColorForEachMode(AssociateBlock==1,:)},'transparent',1,'patchSaturation',0.5);
     
     shadedErrorBar(TimeForAll,IntermedaiteResp_Mean,IntermedaiteResp_Sem,'lineprops',{(ColorForEachMode(AssociateBlock==1,:)+ColorForEachMode(AssociateBlock==11,:))/2},'transparent',1,'patchSaturation',0.5);
     set(findall(gca, 'Type', 'Line'),'LineWidth',3);
     
     xlim([-PlotStartOffset(1),1000]);
     if ~MultiChannel
     xlabel('Time from Scene Onset(ms)');
     ylabel('Firing Rate(Hz)');
     set(gca,'FontSize',15,'FontWeight','Bold','LineWidth',3);
     legend('StableFixLoc','FlexibleLoc','Intermediate');
     end
     
     if ~MultiChannel
         %Dont't show this figure if multichannel, too crowded
     subplot(1,2,2)
      shadedErrorBar(TimeForAll,StableResp_Mean,StableResp_Sem,'lineprops',{ColorForEachMode(AssociateBlock==11,:)},'transparent',1,'patchSaturation',0.5);
     
     hold on 
     shadedErrorBar(TimeForAll,FlexibleResp_Mean,FlexibleResp_Sem,'lineprops',{ColorForEachMode(AssociateBlock==1,:)},'transparent',1,'patchSaturation',0.5);
     
    
      shadedErrorBar(TimeForAll,IntermedaiteResp_Stable_Mean,IntermedaiteResp_Stable_Sem,'lineprops',{InterMediate_StableColor},'transparent',1,'patchSaturation',0.5);
       shadedErrorBar(TimeForAll,IntermedaiteResp_Flexible_Mean,IntermedaiteResp_Flexible_Sem,'lineprops',{InterMediate_FlexibleColor},'transparent',1,'patchSaturation',0.5);
     
     set(findall(gca, 'Type', 'Line'),'LineWidth',3);
     
     xlim([-PlotStartOffset(1),1000]);
     xlabel('Time from Scene Onset(ms)');
    % ylabel('Firing Rate(Hz)');
     set(gca,'FontSize',15,'FontWeight','Bold','LineWidth',3);
       legend('StableFixLoc','FlexibleLoc','IntermediateStable','IntermediateFlexible');
     end
     
     FigureIndex=FigureIndex+1;
    FigureStartNum=FigureStartNum+1; 
       
    if ~MultiChannel
         %Dont't show this figure if multichannel, too crowded
         
   %%%%%%%%%%%%%%%%%%Next figure:Scene Selective Index Comparison betweeen
   %%%%%%%%%%%%%%%%%%stable and flexible scenes
    figtitlestr{3}='SIForStableAndFlexibleScenes';
     fig{3}=PrepareFigure(112,'w',[140,100, 800,400],'Name',figtitlestr{3});   
      subplot(1,2,1)
      g=bar([StableResp_Mean_Whole,FlexibleResp_Mean_Whole,IntermedaiteResp_Mean_Whole,IntermedaiteResp_Stable_Mean_Whole,IntermedaiteResp_Flexible_Mean_Whole],'LineWidth',3);
      hold on
      errorbar([StableResp_Mean_Whole,FlexibleResp_Mean_Whole,IntermedaiteResp_Mean_Whole,IntermedaiteResp_Stable_Mean_Whole,IntermedaiteResp_Flexible_Mean_Whole],[StableResp_Sem_Whole,FlexibleResp_Sem_Whole,IntermedaiteResp_Sem_Whole,IntermedaiteResp_Stable_Sem_Whole,IntermedaiteResp_Flexible_Sem_Whole],'sk','LineWidth',3);
       
      g.FaceColor = 'flat';
      g.EdgeColor = 'flat';
     

%g.CData = ColorChosen([Index_Sort,NumOfCondition+1],:);
ColorShow(1,:)=ColorForEachMode(AssociateBlock==11,:);
ColorShow(2,:)=ColorForEachMode(AssociateBlock==1,:);
ColorShow(3,:)=(ColorForEachMode(AssociateBlock==1,:)+ColorForEachMode(AssociateBlock==11,:))/2;
ColorShow(4,:)=InterMediate_StableColor;
ColorShow(5,:)=InterMediate_FlexibleColor;

g.CData = ColorShow;
set(gca,'XTickLabel',{'Stable','Flexible','Interm','IntermS','IntermF'});
ylabel('Reponses(Hz)');
xlim([0.2,5.8]);
title('Average Response');
box off
set(gca,'FontWeight','Bold','FontSize',15,'LineWidth',3);
  
     subplot(1,2,2)
       b=bar([1,2,3,4,5],[SI_Stable,SI_Flexible,SI_Intermediate,SI_IntermediateStable,SI_IntermediateFlexible]);
       b.FaceColor = 'flat';
      b.EdgeColor = 'flat';

      b.CData = ColorShow;
      set(gca,'XTickLabel',{'Stable','Flexible','Interm','IntermS','IntermF'});
ylabel('SI');

      xlim([0.2,5.8]);
      title('Scene Selective Index');
      box off
      set(gca,'FontWeight','Bold','FontSize',15,'LineWidth',3);
       
     FigureIndex=FigureIndex+1;
    FigureStartNum=FigureStartNum+1;   
    
    
    
    
    %SI_Stable_SceneEach
    figtitlestr{4}='DynamicsOfSI';
    fig{4}=PrepareFigure(113,'w',[140,100, 800,400],'Name',figtitlestr{4});
    subplot(1,2,1);
    plot(TimeForAll,SI_Stable_SceneEach,'-k','Color',ColorForEachMode(AssociateBlock==11,:),'LineWidth',3);
    hold on 
    plot(TimeForAll,SI_Flexible_SceneEach,'-k','Color',ColorForEachMode(AssociateBlock==1,:),'LineWidth',3);
    plot(TimeForAll,SI_Intermediate_SceneEach,'-k','Color',(ColorForEachMode(AssociateBlock==1,:)+ColorForEachMode(AssociateBlock==11,:))/2,'LineWidth',3);
      box off
      xlabel('Time from scene onset(ms)');
      ylabel('Scene Selective Index');
      legend('Stable','Flexible','Inter');
      title('Dynamics of SI (S.vs.F.vs.I)');
      set(gca,'FontWeight','Bold','FontSize',15,'LineWidth',3);
      subplot(1,2,2);
       plot(TimeForAll,SI_IntermediateS_SceneEach,'-k','Color',InterMediate_StableColor,'LineWidth',3);
         hold on 
       plot(TimeForAll,SI_IntermediateF_SceneEach,'-k','Color',InterMediate_FlexibleColor,'LineWidth',3);
   box off
      xlabel('Time from scene onset(ms)');
      ylabel('Scene Selective Index');
      legend('InterS','InterF');
      title('Dynamics of SI (IS.vs.IF)')
      set(gca,'FontWeight','Bold','FontSize',15,'LineWidth',3);
      
    
    
    FigureIndex=FigureIndex+1;
    FigureStartNum=FigureStartNum+1;   
    end%End of the if multichannle
 
 end%End of the if scenetype sep
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
OutputFigureName=sprintf('%s_%s_N%s_C%s',MonkeyName,string(RecordDateOriginal),string(NeuronNum),string(ChannelNum));



ExistFlag=0;
%Load the old file if exist
if exist(OutputFileName)

    load(OutputFileName);
    %Load OutputData into the memory
    if isfield(OutputFileName,'SceneTuningModes')
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
%%%%%%%%%%%%%%%%%%%%OutputData For current analysis
OutputData_New=[];

OutputData_New.SceneTuningModes.TrialType='GoodTrials';

%%%%%%%%%Save overall scene response characteristics%%%%%%%%%%%%%%%%%%%%%%%
%Organize the data
NameStr={'StableScenePSTH','FlexibleScenePSTH','IntermediateScenePSTH',...
    'IntermediateStablePSTH','IntermediateStablePSTH',...
    'StableFlexibleIntermediateISIF_MeanSEM',...
    'SI_SFIISIF',...
    'SI_Dynamics_SFIISIF',...
    'p_SFD',...
    'ResponseEachSceneWhole',...
    'ResponseInterSceneSFEachModeWhole'
    
    
    
    
    

};

DataLib={[TimeForAll',StableResp_Mean',StableResp_Sem'],[TimeForAll',FlexibleResp_Mean',FlexibleResp_Sem'],[TimeForAll',IntermedaiteResp_Mean',IntermedaiteResp_Sem'],...
[TimeForAll',IntermedaiteResp_Stable_Mean',IntermedaiteResp_Stable_Sem'],[TimeForAll',IntermedaiteResp_Flexible_Mean',IntermedaiteResp_Flexible_Sem'],...
[[StableResp_Mean_Whole;FlexibleResp_Mean_Whole;IntermedaiteResp_Mean_Whole;IntermedaiteResp_Stable_Mean_Whole;IntermedaiteResp_Flexible_Mean_Whole],[StableResp_Sem_Whole;FlexibleResp_Sem_Whole;IntermedaiteResp_Sem_Whole;IntermedaiteResp_Stable_Sem_Whole;IntermedaiteResp_Flexible_Sem_Whole]],...
[SI_Stable,SI_Flexible,SI_Intermediate,SI_IntermediateStable,SI_IntermediateFlexible],...
[TimeForAll',SI_Stable_SceneEach',SI_Flexible_SceneEach',SI_Intermediate_SceneEach',SI_IntermediateS_SceneEach',SI_IntermediateF_SceneEach'],...
[p_Stable_Flexible,p_Stable_Flexible_Dual],...
[UniqueScenes,Scene_Resp_Mean_Whole',Scene_Resp_Sem_Whole'],...
[Intermediate_Stable_Scene_Mean_Whole',Intermediate_Stable_Scene_Sem_Whole',Intermediate_Flexible_Scene_Mean_Whole',Intermediate_Flexible_Scene_Sem_Whole']

};





DataStamp=containers.Map(NameStr,DataLib);

OutputData_New.SceneTuningModes.Task=TaskType;
OutputData_New.SceneTuningModes.TaskCode=TaskCode;
OutputData_New.SceneTuningModes.DataStamp=DataStamp;
%{
%Store figures
if ShowFigureFlag
  if ~MultiChannel
        
    %OutputData_New.SceneTuning.Figures={fig99,fig100,fig101,fig102};
    %Output figures in a folder;
    FolderName=OutputFigureName;
    if ~exist(FolderName,'dir')
        mkdir(FolderName);
    end
    cd(FolderName);
    for i=1:length(fig)
        OutputFigureNameFig=[];
        titlestr=figtitlestr{i};
        OutputFigureNameFig=strcat(OutputFigureName,titlestr,'.jpg');
        saveas(fig{i},OutputFigureNameFig);
        
        
    end
        
    
    disp('Figures have been exported to the neuron folder');
  else
      %Only save for the last one
      if spk==length(SpikeChannelWhole)
          FolderName=OutputFigureName;
         if ~exist(FolderName,'dir')
              mkdir(FolderName);
         end
            cd(FolderName);
         for i=1:length(fig)
              OutputFigureNameFig=[];
               titlestr=figtitlestr{i};
              OutputFigureNameFig=strcat(OutputFigureName,titlestr,'.jpg');
              saveas(fig{i},OutputFigureNameFig);
        
        
         end
        
    
      disp('Figures have been exported to the neuron folder');
          
      end
  end
  

end %End of the if figure show flag
%}






if ExistFlag
%Compare the old one with the new one
Task_Old=OutputData.SceneTuningModes.TaskCode;

     if Task_Old~=TaskCode
        %Not the same task, add the new dataset into the old one
        OutputData.SceneTuningModes(numel(OutputData.SceneTuningModes)+1)=OutputData_New.SceneTuningModes;
     else
         %If the same task,the same protocol,replace the old one with the new one
         
         OutputData.SceneTuningModes=OutputData_New.SceneTuningModes;
         
    
    
     end
else
    %Output the current file
   OutputData.SceneTuningModes=OutputData_New.SceneTuningModes; 
end



if OutputFlag
    cd(OutPath);
    save(OutputFileName,'OutputData');
    disp('Data Exported');
else
    disp('No data export');
end  

end %End of the spike channel loop

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