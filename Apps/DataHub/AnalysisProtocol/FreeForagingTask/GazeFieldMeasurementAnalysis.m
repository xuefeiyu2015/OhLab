function GazeFieldMeasurementAnalysis(Data, StartTrial, EndTrial,ShowFigureFlag,OutputFlag);
%Gaze field neuronal responses analysis 
%by YXF 06302021
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




 %Parameter for 2D histogram for spikes
%2D Density map
TwoDBinRange=[-25,25;-25,25];%Width;Height
BinDim=[2;2];%Width;Height
Param.Range=TwoDBinRange;
Param.Bin=BinDim;
%Param.FR_Threshold=100;%ms/deg2
Param.FR_Threshold=100;%ms/degree
Param.GaussianSmooth=1;
Param.GaussianSigma=4;
Param.GaussianWindow=[13;13];
Param.OriBinWidth=EyeBinWidth;

%Load the spike channel
SpikeChannelWhole=Data.SpikeChannel;
SpikeTimeWhole=Data.SpikeTimeData(:,:,Selection);
SpikeBinWhole=Data.SpikeBinData(:,:,Selection);

ChannelIDWhole=Data.ChannelID;
ClusterIDWhole=Data.ClusterID;

MultiChannel=0;
if length(SpikeChannelWhole)>1
    MultiChannel=1;
    
end
for spk=1:length(SpikeChannelWhole)
    SpikeChannel=SpikeChannelWhole(spk);
    SpikeTime=squeeze(SpikeTimeWhole(spk,:,:))';
    SpikeBin=squeeze(SpikeBinWhole(spk,:,:))';
    ChannelID=ChannelIDWhole(spk);
    ClusterID=ClusterIDWhole(spk);
    
    
 %Find out the eyetrace from search onset to scene offset
%EyeInBin=FindOutTime(EventChannel,EventBin,SCENE_ON )+1000;
EyeInBin_whole=FindOutTime(EventChannel,EventBin,SEARCHSTARTCD);
EyeInBin=EyeInBin_whole(:,1);


EyeOffBin=FindOutTime(EventChannel,EventBin,SCENE_OFF);

Eye_Explore_X=SelectEyeInterval(EyeChannel_X_L,EyeInBin,EyeOffBin);
Eye_Explore_Y=SelectEyeInterval(EyeChannel_Y_L,EyeInBin,EyeOffBin);

Eye_Explore_X_Line=ReorganizeEye(Eye_Explore_X);
Eye_Explore_Y_Line=ReorganizeEye(Eye_Explore_Y);

CommonTimeSeq=1:length(Eye_Explore_X_Line);
TimeSeqExploreRaw=[EyeInBin,EyeOffBin];

%Convert spike into a digital sequence
SpikeEventWhole=ConvertSpikeEvent(SpikeBin);


%Select out the interval between scene on and scene off
SpikeEventExplore=SelectSpikeInterval(SpikeEventWhole,EyeInBin,EyeOffBin);

SpikeEventExplore_ForPlot=SpikeEventExplore;
SpikeEventExplore_ForPlot(SpikeEventExplore_ForPlot==0)=NaN;

%Separate trials into Good Object Representation, Bad Object
%Representation, Before RewardAgain, During RewardGain and After Reward
%Gain

%Basic Information about the scene and objects
SceneName=Data.SceneIndex(Selection,3);
SceneIndex=Data.SceneIndex(Selection,2);
TargetName=Data.TargetIndex(Selection,2:end,2);
TargetIndex=Data.TargetIndex(Selection,2:end,1);

GoodObjectTrial=TargetIndex(:,1)<=4;
BadObjectTrial=TargetIndex(:,1)>4;

%Get the reward gain information
ParaLib=Data.ParaLib;
RewardGainRaw=ParaLib('RewardGain');
RearchingOutcomeRaw=ParaLib('RearchingOutcome');

if ismember('NeutralObjectMode',keys(ParaLib))
    NeutralObjModeRaw=ParaLib('NeutralObjectMode');
    NeutralObjMode=NeutralObjModeRaw(Selection,2);
else
     NeutralObjMode=0*ones(length(Selection),1);
end

PositionIndexRaw=ParaLib('PositionIndex');


NeutralObjModeUnique=unique(NeutralObjMode);

%First, plot three raw maps about the gaze field before reward change with
%either good or bad object representation

RewardGain=RewardGainRaw(Selection,2);
RewardGainDiff=diff(RewardGain);
RewardGainChangeTrial=find(RewardGainDiff~=0)+1;


RawTrialIndex=Data.SceneIndex(Selection,1);
RearchingOutcome=RearchingOutcomeRaw(ismember(RearchingOutcomeRaw(:,1),RawTrialIndex),:);

RearchingOutcomeMat=ReOrganizeData(RearchingOutcome);
RearchingOutcomeMat=RearchingOutcomeMat(:,2:end);

PositionIndex=PositionIndexRaw(ismember(PositionIndexRaw(:,1),RawTrialIndex),2);

PositionIndexMatSuc=reshape(PositionIndex,9,TrialNum)';



DefaultChangeTrial=15;
MaxTrialSeg=ceil(EndTrial/DefaultChangeTrial)-floor(StartTrial/DefaultChangeTrial);

%Find out the reward gain position 

RewardGainLoc=unique(SceneIndex)-FOURTH_SCENE;

%Now Focus on the gaze period


GazeOnBin=FindOutTime(EventChannel,EventBin,TGTCD)-800;%Forgot the name of the event marker
%GazeOnBin=FindOutTime(EventChannel,EventBin,GAZESTARTCD);

GazeOffBin=FindOutTime(EventChannel,EventBin,TGTCD)-1;

AveFR_GazePeriod=FRForEachEvent(SpikeEventWhole,GazeOnBin,GazeOffBin);


binsize=50;%in ms
PSTH_GazePeriod=PSTHForEachEvent(SpikeEventWhole,GazeOnBin,GazeOffBin,binsize);


%Reorder the FR according to the location index

[AveFR_GazePeriod_GazeLoc,PositionIndex_Reorder]=ReorderFR(PositionIndexMatSuc,AveFR_GazePeriod);

%Reorder the PSTH according to the location index
[PSTH_GazePeriod_GazeLoc,PositionIndex_Reorder]=ReorderPSTH(PositionIndexMatSuc,PSTH_GazePeriod);

%Now Focus on the object 

ObjOnOffset=800;%in ms
ObjOffOffset=200;%in ms


ObjOnBin=FindOutTime(EventChannel,EventBin,TGTCD);%Forgot the name of the event marker
ObjOffBin=FindOutTime(EventChannel,EventBin,TGTOFF);

ObjOnBinPSTH=ObjOnBin-ObjOnOffset;
ObjOffBinPSTH=ObjOffBin+ObjOffOffset;

ObjPresentDur=max(max(ObjOffBin-ObjOnBin));

AveFR_ObjectPeriod=FRForEachEvent(SpikeEventWhole,ObjOnBin,ObjOffBin);
[AveFR_ObjectPeriod_GazeLoc,PositionIndex_Reorder]=ReorderFR(PositionIndexMatSuc,AveFR_ObjectPeriod);


PSTH_ObjectPeriod=PSTHForEachEvent(SpikeEventWhole,ObjOnBinPSTH,ObjOffBinPSTH,binsize);
%Reorder the PSTH according to the location index
[PSTH_ObjectPeriod_GazeLoc,PositionIndex_Reorder]=ReorderPSTH(PositionIndexMatSuc,PSTH_ObjectPeriod);

%Average every num trials
AveNum=15;
Index=1;
for i=1:AveNum:size(AveFR_GazePeriod_GazeLoc,1)
    
   AveFR_GazePeriod_GazeLoc_Ave(Index,:)=nanmean(AveFR_GazePeriod_GazeLoc(i:min(i+AveNum-1,size(AveFR_GazePeriod_GazeLoc,1)),:),1);
   AveFR_GazePeriod_Ave(Index,:)=nanmean(AveFR_GazePeriod(i:min(i+AveNum-1,size(AveFR_GazePeriod_GazeLoc,1)),:),1);
   
   AveFR_GazePeriod_GazeLoc_Sem(Index,:)=nanstd(AveFR_GazePeriod_GazeLoc(i:min(i+AveNum-1,size(AveFR_GazePeriod_GazeLoc,1)),:),[],1)./sqrt(sum(~isnan(AveFR_GazePeriod_GazeLoc(i:min(i+AveNum-1,size(AveFR_GazePeriod_GazeLoc,1)),:)) ,1));
   AveFR_GazePeriod_Sem(Index,:)=nanstd(AveFR_GazePeriod(i:min(i+AveNum-1,size(AveFR_GazePeriod_GazeLoc,1)),:),1)./sqrt(sum(~isnan(AveFR_GazePeriod(i:min(i+AveNum-1,size(AveFR_GazePeriod_GazeLoc,1)),:))));
    
    Index=Index+1;
end



%{
%Reproduce averageFR

%selInterval=[ceil(ObjOnOffset/binsize)+1,ceil((ObjOnOffset+ObjPresentDur)/binsize)];
selInterval=[ceil(ObjOnOffset/binsize)+1,ceil((ObjOnOffset+500)/binsize)];

tmp=cell2mat(cellfun(@(x) cellfun(@(y) nanmean(y(selInterval(1):selInterval(2))),x)',PSTH_ObjectPeriod_GazeLoc,'uniform',0))';
AveFR_ObjectPeriod_GazeLoc=tmp;
%}

% Remove the Time Period of cue on--cue off;tgt on-tgt off from the raw eye
%  trace and spike event
CueOnBin=FindOutTime(EventChannel,EventBin,CUEON);
CueOffBin=FindOutTime(EventChannel,EventBin,CUEOFF);%Forgot the name of the event marker
ObjOffBin=FindOutTime(EventChannel,EventBin,TGTOFF);%Forgot the name of the event marker
RemoveOffBin=ObjOffBin+500;%100ms offset
RemoveInterval=CombineInterval(CueOnBin,CueOffBin,RemoveOffBin,RearchingOutcomeMat);
%1,:,: StartBin;2::EndBin

PureExploreEyeTrace_X=RemoveFromInterval(Eye_Explore_X_Line,TimeSeqExploreRaw,RemoveInterval);
PureExploreEyeTrace_Y=RemoveFromInterval(Eye_Explore_Y_Line,TimeSeqExploreRaw,RemoveInterval);
SpikeEvent_PureExplore=RemoveFromInterval(SpikeEventExplore_ForPlot,TimeSeqExploreRaw,RemoveInterval);



MaxValObj=0;
MinDiffVal=1000;
%Loop over the trial seg
for i=1:MaxTrialSeg
    if i==1
        TrialSegCurr=StartTrial:15;
    else
        TrialSegCurr=(i-1)*15+1:min(i*15,EndTrial);
    end
    try
    RewardGainSeg(i)=unique(RewardGain(TrialSegCurr));
    catch
        RewardGainSeg(i)=1;%Will solve this tomorrow
        disp('Auto reward gain');
    end
        
    
    GoodSel=ismember(Selection,TrialSegCurr)' & GoodObjectTrial;
    BadSel=ismember(Selection,TrialSegCurr)' & BadObjectTrial;
    
    EyeExplore_X_Good{i}=Eye_Explore_X_Line(GoodSel,:);
    EyeExplore_Y_Good{i}=Eye_Explore_Y_Line(GoodSel,:);
    
    EyeExplore_X_Bad{i}=Eye_Explore_X_Line(BadSel,:);
    EyeExplore_Y_Bad{i}=Eye_Explore_Y_Line(BadSel,:);
    
    SpikeEventExplore_Good{i}=SpikeEventExplore(GoodSel,:);
    SpikeEventExplore_Bad{i}=SpikeEventExplore(BadSel,:);
    
    %For raw plotting
    
    
     EyeExplore_X_GoodLine_tmp=[Eye_Explore_X_Line(GoodSel,:),NaN*ones(size(Eye_Explore_X_Line(GoodSel,:),1),1)];
     EyeExplore_Y_GoodLine_tmp=[Eye_Explore_Y_Line(GoodSel,:),NaN*ones(size(Eye_Explore_Y_Line(GoodSel,:),1),1)];
   
     
     EyeExplore_X_Good_Line{i}=reshape(EyeExplore_X_GoodLine_tmp',1,numel(EyeExplore_X_GoodLine_tmp));
      EyeExplore_Y_Good_Line{i}=reshape(EyeExplore_Y_GoodLine_tmp',1,numel(EyeExplore_Y_GoodLine_tmp));
      
      
     
     EyeExplore_X_BadLine_tmp=[Eye_Explore_X_Line(BadSel,:),NaN*ones(size(Eye_Explore_X_Line(BadSel,:),1),1)];
     EyeExplore_Y_BadLine_tmp=[Eye_Explore_Y_Line(BadSel,:),NaN*ones(size(Eye_Explore_Y_Line(BadSel,:),1),1)];
   
     
     EyeExplore_X_Bad_Line{i}=reshape(EyeExplore_X_BadLine_tmp',1,numel(EyeExplore_X_BadLine_tmp));
     EyeExplore_Y_Bad_Line{i}=reshape(EyeExplore_Y_BadLine_tmp',1,numel(EyeExplore_Y_BadLine_tmp));  
     
     
    SpikeEventExplore_ForPlot_Good_tmp=[SpikeEventExplore_ForPlot(GoodSel,:),NaN*ones(size(SpikeEventExplore_ForPlot(GoodSel,:),1),1)];
     SpikeEventExplore_ForPlot_Bad_tmp=[SpikeEventExplore_ForPlot(BadSel,:),NaN*ones(size(SpikeEventExplore_ForPlot(BadSel,:),1),1)];
   
     SpikeEventExplore_ForPlot_Good{i}=reshape(SpikeEventExplore_ForPlot_Good_tmp',1,numel(SpikeEventExplore_ForPlot_Good_tmp));
    SpikeEventExplore_ForPlot_Bad{i}=reshape(SpikeEventExplore_ForPlot_Bad_tmp',1,numel(SpikeEventExplore_ForPlot_Bad_tmp));
    
    
    %Select out the eye trace for pure exploration
    
    
    
     PureExploreEyeTrace_X_GoodLine=[PureExploreEyeTrace_X(GoodSel,:),NaN*ones(size(PureExploreEyeTrace_X(GoodSel,:),1),1)];
     PureExploreEyeTrace_Y_GoodLine=[PureExploreEyeTrace_Y(GoodSel,:),NaN*ones(size(PureExploreEyeTrace_Y(GoodSel,:),1),1)];
   
     
     EyePureExplore_X_Good_Line{i}=reshape(PureExploreEyeTrace_X_GoodLine',1,numel(PureExploreEyeTrace_X_GoodLine));
     EyePureExplore_Y_Good_Line{i}=reshape(PureExploreEyeTrace_Y_GoodLine',1,numel(PureExploreEyeTrace_Y_GoodLine));
     
     SpikeEventPureExplore_ForPlot_Good_tmp=[SpikeEvent_PureExplore(GoodSel,:),NaN*ones(size(SpikeEvent_PureExplore(GoodSel,:),1),1)];
     SpikeEventPureExplore_ForPlot_Bad_tmp=[SpikeEvent_PureExplore(BadSel,:),NaN*ones(size(SpikeEvent_PureExplore(BadSel,:),1),1)];
     
     SpikeEventPureExplore_ForPlot_Good{i}=reshape(SpikeEventPureExplore_ForPlot_Good_tmp',1,numel(SpikeEventPureExplore_ForPlot_Good_tmp));
    SpikeEventPureExplore_ForPlot_Bad{i}=reshape(SpikeEventPureExplore_ForPlot_Bad_tmp',1,numel(SpikeEventPureExplore_ForPlot_Bad_tmp));
    
  

    %2D Density map raw
     data_tmp=TwoDimHist(EyePureExplore_X_Good_Line{i},EyePureExplore_Y_Good_Line{i}, SpikeEventPureExplore_ForPlot_Good{i},Param);
      
      FR_Mat_Block_Raw{i}=data_tmp.HistData;
 
      tmp=data_tmp.HistX;
      tmp=tmp(1:end-1)+diff(tmp)/2;
      ScreenLoc_X_Block_Raw{i}=tmp;
      
      tmp=data_tmp.HistY;
      tmp=tmp(1:end-1)+diff(tmp)/2;
      ScreenLoc_Y_Block_Raw{i}=tmp;
      %{
      figure
     % heatmap(FR_Mat_Block_Raw{i});
      contourf(FR_Mat_Block_Raw{i});
      colormap default;
      %}
    
    %Now focus only on the gaze period 
    
   % PositionIndexMatSuc_Good=PositionIndexMatSuc(GoodSel,:);
   % PositionIndexMatSuc_Bad=PositionIndexMatSuc(BadSel,:);
    
    
    AveFR_GazePeriod_GazeLoc_Good=AveFR_GazePeriod_GazeLoc(GoodSel,:);
    AveFR_GazePeriod_GazeLoc_Bad=AveFR_GazePeriod_GazeLoc(BadSel,:);
    
    AveFR_GazePeriod_GazeLoc_Good_Mean=nanmean(AveFR_GazePeriod_GazeLoc_Good,1);
    AveFR_GazePeriod_GazeLoc_Bad_Mean=nanmean(AveFR_GazePeriod_GazeLoc_Bad,1);
    
    AveFR_GazePeriod_GazeLoc_Good_Sem=nanstd(AveFR_GazePeriod_GazeLoc_Good,[],1)./sqrt(sum(~isnan(AveFR_GazePeriod_GazeLoc_Good),1));
    AveFR_GazePeriod_GazeLoc_Bad_Sem=nanstd(AveFR_GazePeriod_GazeLoc_Bad,[],1)./sqrt(sum(~isnan(AveFR_GazePeriod_GazeLoc_Bad),1));
    
    
    tmp=[AveFR_GazePeriod_GazeLoc_Good_Mean(3),AveFR_GazePeriod_GazeLoc_Good_Mean(2),AveFR_GazePeriod_GazeLoc_Good_Mean(9);...
        AveFR_GazePeriod_GazeLoc_Good_Mean(4),AveFR_GazePeriod_GazeLoc_Good_Mean(1),AveFR_GazePeriod_GazeLoc_Good_Mean(8);...
        AveFR_GazePeriod_GazeLoc_Good_Mean(5),AveFR_GazePeriod_GazeLoc_Good_Mean(6),AveFR_GazePeriod_GazeLoc_Good_Mean(7);];
   % tmp=tmp/max(max(tmp));
    AveFR_GazePeriod_GazeLoc_Good_Mean_Reshape{i}=tmp;
    
     AveFR_GazePeriod_GazeLoc_Bad_Mean_Reshape{i}=[AveFR_GazePeriod_GazeLoc_Bad_Mean(3),AveFR_GazePeriod_GazeLoc_Bad_Mean(2),AveFR_GazePeriod_GazeLoc_Bad_Mean(9);...
        AveFR_GazePeriod_GazeLoc_Bad_Mean(4),AveFR_GazePeriod_GazeLoc_Bad_Mean(1),AveFR_GazePeriod_GazeLoc_Bad_Mean(8);...
        AveFR_GazePeriod_GazeLoc_Bad_Mean(5),AveFR_GazePeriod_GazeLoc_Bad_Mean(6),AveFR_GazePeriod_GazeLoc_Bad_Mean(7);];
    
     %Statistics: one way anova
     %Gaze for good and bad
     p_GazePeriod_Good(i)=anova1(AveFR_GazePeriod_GazeLoc_Good,[],'off');
     p_GazePeriod_Bad(i)=anova1(AveFR_GazePeriod_GazeLoc_Bad,[],'off');
     
 
    
    %Now focus on the whole object period
     AveFR_ObjectPeriod_GazeLoc_Good=AveFR_ObjectPeriod(GoodSel,:);
    AveFR_ObjectPeriod_GazeLoc_Bad=AveFR_ObjectPeriod(BadSel,:);
    
    AveFR_ObjectPeriod_GazeLoc_Good_Mean=nanmean(AveFR_ObjectPeriod_GazeLoc_Good,1);
    AveFR_ObjectPeriod_GazeLoc_Bad_Mean=nanmean(AveFR_ObjectPeriod_GazeLoc_Bad,1);
    
    AveFR_ObjectPeriod_GazeLoc_Good_Sem=nanstd(AveFR_ObjectPeriod_GazeLoc_Good,[],1)./sqrt(sum(~isnan(AveFR_ObjectPeriod_GazeLoc_Good),1));
    AveFR_ObjectPeriod_GazeLoc_Bad_Sem=nanstd(AveFR_ObjectPeriod_GazeLoc_Bad,[],1)./sqrt(sum(~isnan(AveFR_ObjectPeriod_GazeLoc_Bad),1));
    
    
    AveFR_ObjectPeriod_GazeLoc_Good_Mean_Reshape{i}=[AveFR_ObjectPeriod_GazeLoc_Good_Mean(3),AveFR_ObjectPeriod_GazeLoc_Good_Mean(2),AveFR_ObjectPeriod_GazeLoc_Good_Mean(9);...
        AveFR_ObjectPeriod_GazeLoc_Good_Mean(4),AveFR_ObjectPeriod_GazeLoc_Good_Mean(1),AveFR_ObjectPeriod_GazeLoc_Good_Mean(8);...
        AveFR_ObjectPeriod_GazeLoc_Good_Mean(5),AveFR_ObjectPeriod_GazeLoc_Good_Mean(6),AveFR_ObjectPeriod_GazeLoc_Good_Mean(7);];
    
     AveFR_ObjectPeriod_GazeLoc_Bad_Mean_Reshape{i}=[AveFR_ObjectPeriod_GazeLoc_Bad_Mean(3),AveFR_ObjectPeriod_GazeLoc_Bad_Mean(2),AveFR_ObjectPeriod_GazeLoc_Bad_Mean(9);...
        AveFR_ObjectPeriod_GazeLoc_Bad_Mean(4),AveFR_ObjectPeriod_GazeLoc_Bad_Mean(1),AveFR_ObjectPeriod_GazeLoc_Bad_Mean(8);...
        AveFR_ObjectPeriod_GazeLoc_Bad_Mean(5),AveFR_ObjectPeriod_GazeLoc_Bad_Mean(6),AveFR_ObjectPeriod_GazeLoc_Bad_Mean(7);];
    
     %Statistics: one way anova
     %Object for good and bad
     p_ObjectPeriod_Good(i)=anova1(AveFR_ObjectPeriod_GazeLoc_Good,[],'off');
     p_ObjectPeriod_Bad(i)=anova1(AveFR_ObjectPeriod_GazeLoc_Bad,[],'off');
     
     
     %Object Response Gain
     
     AveFR_ObjectGain_GazeLoc_Good=AveFR_ObjectPeriod_GazeLoc_Good-AveFR_GazePeriod_GazeLoc_Good;
     AveFR_ObjectGain_GazeLoc_Bad=AveFR_ObjectPeriod_GazeLoc_Bad-AveFR_GazePeriod_GazeLoc_Bad;
     
     AveFR_ObjectGain_GazeLoc_Good_Mean=AveFR_ObjectPeriod_GazeLoc_Good_Mean-AveFR_GazePeriod_GazeLoc_Good_Mean;
     AveFR_ObjectGain_GazeLoc_Bad_Mean=AveFR_ObjectPeriod_GazeLoc_Bad_Mean-AveFR_GazePeriod_GazeLoc_Bad_Mean;
    
      AveFR_ObjectGain_GazeLoc_Good_Mean_Reshape{i}=[AveFR_ObjectGain_GazeLoc_Good_Mean(3),AveFR_ObjectGain_GazeLoc_Good_Mean(2),AveFR_ObjectGain_GazeLoc_Good_Mean(9);...
        AveFR_ObjectGain_GazeLoc_Good_Mean(4),AveFR_ObjectGain_GazeLoc_Good_Mean(1),AveFR_ObjectGain_GazeLoc_Good_Mean(8);...
        AveFR_ObjectGain_GazeLoc_Good_Mean(5),AveFR_ObjectGain_GazeLoc_Good_Mean(6),AveFR_ObjectGain_GazeLoc_Good_Mean(7);];
    
     AveFR_ObjectGain_GazeLoc_Bad_Mean_Reshape{i}=[AveFR_ObjectGain_GazeLoc_Bad_Mean(3),AveFR_ObjectGain_GazeLoc_Bad_Mean(2),AveFR_ObjectGain_GazeLoc_Bad_Mean(9);...
        AveFR_ObjectGain_GazeLoc_Bad_Mean(4),AveFR_ObjectGain_GazeLoc_Bad_Mean(1),AveFR_ObjectGain_GazeLoc_Bad_Mean(8);...
        AveFR_ObjectGain_GazeLoc_Bad_Mean(5),AveFR_ObjectGain_GazeLoc_Bad_Mean(6),AveFR_ObjectGain_GazeLoc_Bad_Mean(7);];
    
    
     %Statistics: one way anova
     %ObjectGain for good and bad
     p_ObjectGain_Good(i)=anova1(AveFR_ObjectGain_GazeLoc_Good,[],'off');
     p_ObjectGain_Bad(i)=anova1(AveFR_ObjectGain_GazeLoc_Bad,[],'off');
     
     
    
    %psth
    PSTH_GazePeriod_GazeLoc_Good=PSTH_GazePeriod_GazeLoc(GoodSel');
    PSTH_GazePeriod_GazeLoc_Bad=PSTH_GazePeriod_GazeLoc(BadSel');
  
    
    [PSTH_GazePeriod_GazeLoc_Good_Mean,PSTH_GazePeriod_GazeLoc_Good_Sem]=PSTH_LocAve(PSTH_GazePeriod_GazeLoc_Good);
    [PSTH_GazePeriod_GazeLoc_Bad_Mean,PSTH_GazePeriod_GazeLoc_Bad_Sem]=PSTH_LocAve(PSTH_GazePeriod_GazeLoc_Bad);
    
     PSTH_ObjectPeriod_GazeLoc_Good=PSTH_ObjectPeriod_GazeLoc(GoodSel');
    PSTH_ObjectPeriod_GazeLoc_Bad=PSTH_ObjectPeriod_GazeLoc(BadSel');
  
    
    [PSTH_ObjectPeriod_GazeLoc_Good_Mean,PSTH_ObjectPeriod_GazeLoc_Good_Sem]=PSTH_LocAve(PSTH_ObjectPeriod_GazeLoc_Good);
    [PSTH_ObjectPeriod_GazeLoc_Bad_Mean,PSTH_ObjectPeriod_GazeLoc_Bad_Sem]=PSTH_LocAve(PSTH_ObjectPeriod_GazeLoc_Bad);
    
    
    PlotOrder=[3,2,9,4,1,8,5,6,7];
    
    PSTH_ObjectPeriod_GazeLoc_Good_Mean_Plot{i}=PSTH_ObjectPeriod_GazeLoc_Good_Mean(PlotOrder);
    PSTH_ObjectPeriod_GazeLoc_Bad_Mean_Plot{i}=PSTH_ObjectPeriod_GazeLoc_Bad_Mean(PlotOrder);
    
    PSTH_ObjectPeriod_GazeLoc_Good_Sem_Plot{i}=PSTH_ObjectPeriod_GazeLoc_Good_Mean(PlotOrder);
    PSTH_ObjectPeriod_GazeLoc_Bad_Sem_Plot{i}=PSTH_ObjectPeriod_GazeLoc_Bad_Mean(PlotOrder);
    
    MaxValObjSum1=cellfun(@(x,y) x+y,PSTH_ObjectPeriod_GazeLoc_Good_Mean,PSTH_ObjectPeriod_GazeLoc_Good_Sem,'uniform',0);
    MaxValObjSum2=cellfun(@(x,y) x+y,PSTH_ObjectPeriod_GazeLoc_Bad_Mean,PSTH_ObjectPeriod_GazeLoc_Bad_Sem,'uniform',0);
    
    MaxValObjtmp1=cellfun(@max,MaxValObjSum1);
    MaxValObjtmp2=cell2mat(cellfun(@max,MaxValObjSum2,'uniform',0));
    
    MaxValObjtmp=max([max(MaxValObjtmp1),max(MaxValObjtmp2)]);
    if MaxValObjtmp>MaxValObj
        MaxValObj=MaxValObjtmp;
        
    end
    %{
     MinValObjDiff1=cellfun(@(x,y) x-y,PSTH_ObjectPeriod_GazeLoc_Good_Mean,PSTH_ObjectPeriod_GazeLoc_Bad_Mean,'uniform',0);
    MinDiffValtmp=min(cellfun(@min,MinValObjDiff1));
    if MinDiffValtmp<MinDiffVal
        MinDiffVal=MinDiffValtmp;
        
    end
    %}
    
 
end

MaxValObj=ceil(MaxValObj);

%
if ShowFigureFlag
    
FigureStartNum=150;
FigureIndex=1;  


figtitlestr{FigureIndex}='RawSpikeEventMap';
fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[80,100, 1200,800],'Name',figtitlestr{FigureIndex});
%Size of the spike marker
SpikeSizePlot=15;

for i=1:MaxTrialSeg
    %{
    if NeutralObjModeUnique==1
        rownum=1;
    else
        rownum=2;
    end
    %}
subplot(3,MaxTrialSeg,i);
%Search for good
plot(EyeExplore_X_Good_Line{i}',EyeExplore_Y_Good_Line{i}','-w');
hold on
%Overlay spike
s=scatter(EyeExplore_X_Good_Line{i}',EyeExplore_Y_Good_Line{i}',SpikeSizePlot*SpikeEventExplore_ForPlot_Good{i}',[238,60,11]/255,'filled');

xlim([-40,40]);
ylim([-40,40]);
set(gca,'Color','k','xtick',[],'ytick',[]);
box off

title(sprintf('Search For Good: %d',i),'FontSize',15);

subplot(3,MaxTrialSeg,i+MaxTrialSeg);
%Search for bad
plot(EyePureExplore_X_Good_Line{i}',EyePureExplore_Y_Good_Line{i}','-w');
hold on
%Overlay spike
s=scatter(EyePureExplore_X_Good_Line{i}',EyePureExplore_Y_Good_Line{i}',SpikeSizePlot*SpikeEventPureExplore_ForPlot_Good{i}',[238,60,11]/255,'filled');

xlim([-40,40]);
ylim([-40,40]);
set(gca,'Color','k','xtick',[],'ytick',[]);
box off

title(sprintf('Pure Search For Good: %d',i),'FontSize',15);

subplot(3,MaxTrialSeg,i+MaxTrialSeg*2);
%Search for bad
[x,y]=meshgrid(ScreenLoc_X_Block_Raw{i},ScreenLoc_Y_Block_Raw{i});
%[m,c]=contourf(x,y,FR_Mat_Block_Raw{i});
h=heatmap(ScreenLoc_X_Block_Raw{i},-ScreenLoc_Y_Block_Raw{i},flipud(FR_Mat_Block_Raw{i}),'GridVisible','off');
colormap default;
h.Title = sprintf('Pure Search For Good: %d',i);
colorbar;
%h.XDisplayLabels={};
%h.YDisplayLabels={};
%xlim([-25,25]);
%ylim([-25,25]);
%set(gca,'Color','k','xtick',[],'ytick',[]);
%box off

%title(sprintf('Pure Search For Good: %d',i),'FontSize',15);


%{
      FR_Mat_Block_Raw{i}=data_tmp.HistData;
 
      ScreenLoc_X_Block_Raw{i}=data_tmp.HistX;
      ScreenLoc_Y_Block_Raw{i}=data_tmp.HistY;
 %}

if ~NeutralObjModeUnique

subplot(2,MaxTrialSeg,i+MaxTrialSeg);
%Search for bad
plot(EyeExplore_X_Bad_Line{i}',EyeExplore_Y_Bad_Line{i}','-w');
hold on
%Overlay spike
s=scatter(EyeExplore_X_Bad_Line{i}',EyeExplore_Y_Bad_Line{i}',SpikeSizePlot*SpikeEventExplore_ForPlot_Bad{i}',[11,60,238]/255,'filled');

xlim([-40,40]);
ylim([-40,40]);
set(gca,'Color','k','xtick',[],'ytick',[]);
box off

title(sprintf('Search For Bad: %d',i),'FontSize',15);

end
 
end %End of the max trial segment 
FigureStartNum=FigureStartNum+1;
FigureIndex=FigureIndex+1;  

%Figure 2
figtitlestr{FigureIndex}='PSTHOfObjectPeriodOnEachLocation';
fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[80,100, 1200,500],'Name',figtitlestr{FigureIndex});

ColorForThreePeriod={'#CED8F3','#FFC5D7','#E4D9C0','#CED8F3','#FFC5D7','#E4D9C0','#CED8F3','#FFC5D7','#E4D9C0','#CED8F3','#FFC5D7','#E4D9C0','#CED8F3','#FFC5D7','#E4D9C0'};
ColorForGoodSpot={'#EFAF88'};

for j=1:MaxTrialSeg
    for i=1:9
   
        indexorder=[1,2,3,1+3*MaxTrialSeg,2+3*MaxTrialSeg,3+3*MaxTrialSeg,1+3*MaxTrialSeg*2,2+3*MaxTrialSeg*2,3+3*MaxTrialSeg*2]+(j-1)*3;
        LocOrder=[2,1,8,3,0,7,4,5,6];
    subplot(3,3*MaxTrialSeg,indexorder(i));

    %Good Object
    val1=PSTH_ObjectPeriod_GazeLoc_Good_Mean_Plot{j}{i};
    errorplot1=PSTH_ObjectPeriod_GazeLoc_Good_Sem_Plot{j}{i};
    timeplot=[1:length(val1)]*binsize-ObjOnOffset;
   % errorbar(timeplot,val1,errorplot1,'-r','LineWidth',3);
   % hold on
    plot(timeplot,val1,'-r','LineWidth',3);
    hold on
    
    
   
    
    
    if ~NeutralObjModeUnique
    val2=PSTH_ObjectPeriod_GazeLoc_Bad_Mean_Plot{j}{i};
    errorplot2=PSTH_ObjectPeriod_GazeLoc_Bad_Sem_Plot{j}{i};
    timeplot2=[1:length(val2)]*binsize-ObjOnOffset;
  % errorbar(timeplot2,val2,errorplot2,'-b','LineWidth',3);
   hold on
     plot(timeplot2,val2,'-b','LineWidth',3);
     
     %Add the difference
   %  plot(timeplot2,val1-val2,'-y','LineWidth',3);
     
    end
    
    
    
    % xlim([-800,-400]);
    xlim([min(timeplot),max(timeplot)]);
   % ylim([min(MinDiffVal,0),MaxValObj]);
   ylim([0,MaxValObj]);
    
    
    if ismember(indexorder(i),[2:3:(2+3*(MaxTrialSeg-1))]) 
         title(sprintf('Period:%d',j));
    end
    
    set(gca,'LineWidth',3,'FontSize',10,'FontWeight','Bold');
    
    if  RewardGainSeg(j)>1
        order=find(RewardGainSeg>1);
        if ismember(LocOrder(i),RewardGainLoc(find(order==j))) 
           set(gca,'color',ColorForGoodSpot{1});
        else
           set(gca,'color', ColorForThreePeriod{j});
           
        end
    else
        set(gca,'color', ColorForThreePeriod{j});
       
        
    end
    
    if ismember(indexorder(i),indexorder(end-1)) 
        xlabel('Time from object onset(ms)');
        set(gca,'FontSize',10,'FontWeight','Bold');
    end
    
     if ismember(indexorder(i),[4+3*(MaxTrialSeg-1)]) 
        ylabel('Firing Rate(Hz)');
        set(gca,'FontSize',10,'FontWeight','Bold');
     end
     
     if ~ismember(indexorder(i),[1,1+3*(MaxTrialSeg),1+3*(MaxTrialSeg)*2])
         set(gca,'YTick',[]);
     
     end
    
    
    end
   
  
end

FigureStartNum=FigureStartNum+1;
FigureIndex=FigureIndex+1;  

 
%Figure 3
figtitlestr{FigureIndex}='HeatMapForGazePeriod_ForGood';
fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[100,100, 900,900],'Name',figtitlestr{FigureIndex});
for j=1:MaxTrialSeg
    %Gaze Period
    subplot(3,MaxTrialSeg,j);
    h=heatmap(AveFR_GazePeriod_GazeLoc_Good_Mean_Reshape{j});
    colormap default;
   % caxis([0, 1]);
   % p_GazePeriod_Good
   h.Title=sprintf('GazeForGood,Period %d (p=%1.4d)',j,p_GazePeriod_Good(j));
   h.XData=[-15,0,15];
   h.YData=[15,0,-15];
   
   if p_GazePeriod_Good(j)<0.05
      h.FontColor='r';
   else
      h.FontColor='k';
   end
   
   %Object Period
    subplot(3,MaxTrialSeg,MaxTrialSeg+j);
    h=heatmap(AveFR_ObjectPeriod_GazeLoc_Good_Mean_Reshape{j});
    colormap default;
   % caxis([0, 1]);
   % p_GazePeriod_Good
   h.Title=sprintf('ObjectPeriod,Period %d (p=%1.4d)',j,p_ObjectPeriod_Good(j));
   h.XData=[-15,0,15];
   h.YData=[15,0,-15];
   
   if p_ObjectPeriod_Good(j)<0.05
      h.FontColor='r';
   else
      h.FontColor='k';
   end
   
    %Object Gain Period
    subplot(3,MaxTrialSeg,MaxTrialSeg*2+j);
    h=heatmap(AveFR_ObjectGain_GazeLoc_Good_Mean_Reshape{j});
    colormap default;
   % p_GazePeriod_Good
   h.Title=sprintf('ObjectGain,Period %d (p=%1.4d)',j,p_ObjectGain_Good(j));
   h.XData=[-15,0,15];
   h.YData=[15,0,-15];
   
   if p_ObjectGain_Good(j)<0.05
      h.FontColor='r';
   else
      h.FontColor='k';
   end
    
 
end

FigureStartNum=FigureStartNum+1;
FigureIndex=FigureIndex+1;  


%Add Another Plot for bad if there is such a condition
if ~NeutralObjModeUnique
  
%Figure 4
figtitlestr{FigureIndex}='HeatMapForGazePeriod_ForBad';
fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[110,100, 900,900],'Name',figtitlestr{FigureIndex});
for j=1:MaxTrialSeg
    %Gaze Period
    subplot(3,MaxTrialSeg,j);
    h=heatmap(AveFR_GazePeriod_GazeLoc_Bad_Mean_Reshape{j});
    colormap default;
   
   h.Title=sprintf('GazeForBad,Period %d (p=%1.4d)',j,p_GazePeriod_Bad(j));
   h.XData=[-15,0,15];
   h.YData=[15,0,-15];
   
   if p_GazePeriod_Bad(j)<0.05
      h.FontColor='r';
   else
      h.FontColor='k';
   end
   
   %Object Period
    subplot(3,MaxTrialSeg,MaxTrialSeg+j);
    h=heatmap(AveFR_ObjectPeriod_GazeLoc_Bad_Mean_Reshape{j});
    colormap default;
   
   h.Title=sprintf('ObjectPeriod,Period %d (p=%1.4d)',j,p_ObjectPeriod_Bad(j));
   h.XData=[-15,0,15];
   h.YData=[15,0,-15];
   
   if p_ObjectPeriod_Bad(j)<0.05
      h.FontColor='r';
   else
      h.FontColor='k';
   end
   
    %Object Gain Period
    subplot(3,MaxTrialSeg,MaxTrialSeg*2+j);
    h=heatmap(AveFR_ObjectGain_GazeLoc_Bad_Mean_Reshape{j});
    colormap default;
   
   h.Title=sprintf('ObjectGain,Period %d (p=%1.4d)',j,p_ObjectGain_Bad(j));
   h.XData=[-15,0,15];
   h.YData=[15,0,-15];
   
   if p_ObjectGain_Bad(j)<0.05
      h.FontColor='r';
   else
      h.FontColor='k';
   end
    
 
end

FigureStartNum=FigureStartNum+1;
FigureIndex=FigureIndex+1;  
end


%Figure 5
figtitlestr{FigureIndex}='AverageForEachLocationAndEachTimeOrder';
fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[130,100, 900,900],'Name',figtitlestr{FigureIndex});    
 

Color=[lines(6);prism(6);pink(6)];
%Color(2,:)=copper(1)+0.5;
Color(9,:)=[0.28,0.78,0.89];
Color=[Color;hsv(10)];
Color(18,:)=[0.1,0.457,0.194];

EveryForPlot=3;
NumRep=ceil(size(AveFR_GazePeriod_GazeLoc_Ave,1)/EveryForPlot);
index=1;
for i=1:NumRep
    %Spatial Location
    subplot(NumRep,2,index);
    for j=(i-1)*EveryForPlot+1:i*EveryForPlot
     
     
    errorbar(AveFR_GazePeriod_GazeLoc_Ave(j,:),AveFR_GazePeriod_GazeLoc_Sem(j,:),'-k','Color',Color(j,:),'LineWidth',3);
    hold on
    
    end
    legend({'Before','Dur','After'});
    title('Spacial Location');
    xlabel('Gaze Location');
    ylabel('Firing Rate(Hz)');
    set(gca,'FontSize',15,'FontWeight','Bold','LineWidth',3);
    %Temporal Location
    subplot(NumRep,2,index+1);
    for j=(i-1)*EveryForPlot+1:i*EveryForPlot
    errorbar(AveFR_GazePeriod_Ave(j,:),AveFR_GazePeriod_Sem(j,:),'-k','Color',Color(j,:),'LineWidth',3);
    hold on
    end
     title('Temporal Order');
      xlabel('Order');
      ylabel('Firing Rate(Hz)');
      set(gca,'FontSize',15,'FontWeight','Bold','LineWidth',3);
  index=index+2;
end


end %End of the if show figure flag   
    
    
end % End of the spike channel

end %End of the function


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


function SpikeEvent=ConvertSpikeEvent(SpikeBin);

SpikeEvent=0*ones(size(SpikeBin,1),max(SpikeBin(:,end)));

for i=1:size(SpikeBin,1)
    currSpikeBin=SpikeBin(i,1:end-2);
    currSpikeBin=currSpikeBin(~isnan(currSpikeBin));
    
    %Solve this problem later
   % currSpikeBin(currSpikeBin<1)=1;
   SpikeEvent(i,currSpikeBin)=1;
    
    
    SpikeEvent(i,SpikeBin(i,end)+1:end)=NaN;
 
end



end

function SpikeEvent=SelectSpikeInterval(Spike,Start,End);

for i=1:length(Start)
   SpikeEventCell{i}= Spike(i,Start(i):End(i));

end

%Convert to matrix
MaxCol=max(cellfun(@numel, SpikeEventCell));
SpikeEvent=cell2mat(cellfun(@(x) [x,NaN*ones(1,MaxCol-numel(x))]',SpikeEventCell,'uniform',0))';

end

function fr=FRForEachEvent(spk,startbin,endbin)

for i=1:size(startbin,1)
    for j=1:length(startbin(i,:))
        
        tmp=spk(i,startbin(i,j):endbin(i,j));
        fr(i,j)=sum(tmp)/length(tmp(~isnan(tmp)))*1000;
        
    end
    
    
end

end

function psth=PSTHForEachEvent(spk,startbin,endbin,binsize);

for i=1:size(startbin,1)
    for j=1:length(startbin(i,:))
        
        tmp=spk(i,startbin(i,j):endbin(i,j));
        numofbin=ceil(length(tmp)/binsize);
        for ii=1:numofbin
            tmp2(ii)=sum(tmp((ii-1)*binsize+1:min(ii*binsize,length(tmp))));
            timebin(ii)=length(tmp((ii-1)*binsize+1:min(ii*binsize,length(tmp))));
 
            timebinstore(i,j,ii)=length(tmp((ii-1)*binsize+1:min(ii*binsize,length(tmp))));
            
        end
        
       
        psth{i,j}=tmp2./timebin*1000;
    end
  
    
end


end

function [fr_reorder,mat_reorder]=ReorderFR(mat,fr);
for i=1:size(mat,1)
    [ordermat,index]=sort(mat(i,:),'ascend');
    
    fr_reorder(i,:)=fr(i,index);
    mat_reorder(i,:)=ordermat;
    
end




end

function [psth_reorder,mat_reorder]=ReorderPSTH(mat,psth);
for i=1:size(mat,1)
    [ordermat,index]=sort(mat(i,:),'ascend');
    
    psth_reorder{i}=psth(i,index);
    mat_reorder(i,:)=ordermat;
    
   
    
end




end

function [psth_mean,psth_sem]=PSTH_LocAve(psth);

for i=1:9
    psth_tmp_tmp=cellfun(@(x) x(i),psth)';
    
    if ~isempty(psth_tmp_tmp)
    
    maxcol=max(cellfun(@numel,psth_tmp_tmp));
    
    
    psth_tmp_tmp_tmp=cellfun(@(x) [x,NaN*ones(1,maxcol-length(x))],psth_tmp_tmp,'uniform',0);
    
    
    psth_tmp=cell2mat(psth_tmp_tmp_tmp);
    
    colnum=sum(~isnan(psth_tmp),1);
    
    tmp=nanmean(psth_tmp,1);
    tmp(colnum<3)=NaN;
    psth_mean{i}=tmp;
    tmp=nanstd(psth_tmp,[],1)./sqrt(sum(~isnan(psth_tmp),1));
    tmp(colnum<3)=NaN;
    psth_sem{i}=tmp;
    else
        psth_mean{i}=[];
        psth_sem{i}=[];
    end
    
    
end



end

function interval=CombineInterval(CueOnBin,CueOffBin,ObjOffBin,RearchingOutcomeMat);

StartBin=CueOnBin;
EndBin=NaN*ones(size(CueOnBin));

for i=1:size(StartBin,1)
    currOutcome=RearchingOutcomeMat(i,:);
    EndBin(i,currOutcome==1)=ObjOffBin(i,:);
    EndBin(i,currOutcome==0)=CueOffBin(i,~isnan(CueOffBin(i,:)));
    
    
end


interval(1,:,:)=StartBin;
interval(2,:,:)=EndBin;

end

function RemovedTrace=RemoveFromInterval(RawTrace,timeseq,RmvInterval);

for i=1:size(RawTrace,1)
    invStart=squeeze(RmvInterval(1,i,:));
    invEnd=squeeze(RmvInterval(2,i,:));
    
    trace_tmp=RawTrace(i,:);
    timeseq_curr=timeseq(i,1):timeseq(i,2);
    for j=1:sum(~isnan(invStart))
        currinv=[invStart(j),invEnd(j)];
        
        trace_tmp(timeseq_curr>=currinv(1) & timeseq_curr<=currinv(2))=NaN;
   
    end
    RemovedTrace(i,:)=trace_tmp;
    
end

end


