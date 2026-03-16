function RepetitionChangeAnalysis(Data, StartTrial, EndTrial,ShowFigureFlag,OutputFlag);
%Protocol to analysis the neural response change with repetitions
%Xuefei Yu  04092021
ECode;
TaskType=Data.TaskType;
TaskCode=Data.TaskCode;
DataPath= Data.Path;
FileName=Data.FileName;

Selection=StartTrial:EndTrial;%Select the trials according to the pannel
%Abnormal check %Especially For Robin081120FORAGEXPLO4.rst
CheckCode=SCENE_ON;%This code should appear once in each trial
AbnormalTrials=AbnormalTrialCheck(Data,CheckCode);
Selection=Selection(~ismember(Selection,AbnormalTrials));

CheckCode=SCENE_OFF;%This code should appear once in each trial
AbnormalTrials=AbnormalTrialCheck(Data,CheckCode);
Selection=Selection(~ismember(Selection,AbnormalTrials));

LegendOn=0;%Flag to turn legend on
EachSceneObjOn=0;%Flag to turn each scene object figure on

if strcmp(FileName,'Robin092920FORAGEXPLO1.rst') ||strcmp(FileName,'Robin092820FORAGEXPLO1.rst')||strcmp(FileName,'Robin012621FORAGEXPLO1.rst')||strcmp(FileName,'Robin032021SCENEOBJNEW1.rst')
    %Remove bad trial for this file

CheckCodePair=[FIXOFF,SHOWFIXCD];%This code should be the same as target on for good trials
AbnormalTrials=AbnormalTrialCheck(Data,CheckCodePair);
Selection=Selection(~ismember(Selection,AbnormalTrials));

if strcmp(FileName,'Robin092820FORAGEXPLO1.rst') ||strcmp(FileName,'Robin012621FORAGEXPLO1.rst')||strcmp(FileName,'Robin032021SCENEOBJNEW1.rst')
    %Remove bad trial for this file
    CheckCodePair=[TGTCD,TGTOFF];%This code should be the same as target on for good trials
    AbnormalTrials=AbnormalTrialCheck(Data,CheckCodePair);
    Selection=Selection(~ismember(Selection,AbnormalTrials));
    
end

if strcmp(FileName,'Robin032021SCENEOBJNEW1.rst')
    %Remove bad trial for this file
    CheckCodePair=[SHOWFIXCD,FIXOFF];%This code should be the same as target on for good trials
    AbnormalTrials=AbnormalTrialCheck(Data,CheckCodePair);
    Selection=Selection(~ismember(Selection,AbnormalTrials));
    
end


%TGTCD
end





TrialNum=length(Selection);
%User defined parameter for the neural/behavioral analysis
%TargetOnOffset=100;%in ms, the time selection bin after target onset; 
%StartDelay=-50;%in ms, the time selection bin after fixation offset;




%For segmentation of  the neural responses 
SegEvent=[SCENE_ON,SCENE_OFF,TGTCD,TGTOFF];%Event Flow of interest
SegEventEnd=[SHOWFIXCD,NaN,TGTOFF,NaN ];%Event end; %Automatic screening later
%SegEventEnd=[EXTRAFIXOFF,FIXOFF,FIXOFF,TGTCD,TGTOFF];%Event end
SegEventStr={'SceneOn','SceneOff','FractalOn','FractalOff'};

PlotStartOffset=[1000,500,700,500];%in ms
PlotStopOffset=[500,1000,200,600];%in ms


EventTypeFlag=[1,1,2,2];%1:Number same with the trial number, 2:Number same with the target number

%Spon Event
SponEvent=[SCENE_OFF];
SponPreTimeBuffer=-500;
SponEventEnd=NaN;
SponStartOffset=0;
SponStopOffset=1000;

%Parameter to define the sliding window to measure firing rates
StepSize=5;%in ms
BinWidth=15;%in ms

%Plot Raster Flag
PlotRasterFlag=0;




%ResponseLatencyThreshold=3;

%Some user defined parameter for visualization
%Color Definition for each saccade nums
ColorRep={'#FB5441','#4173fb','#FBBB41','#0EC335',[0.28,0.78,0.89],'k',[0.6,0.1,1],[0.6,0.1,0.4],[0.3,0.5,1],[0.6,0.5,0.1],[0,0.8,0.4],'r','g','b','m'};
ColorForEachMode=[0.83 0.14 0.14
             1.00 0.54 0.00
             0.47 0.25 0.80
             0.25 0.80 0.54];
AssociateBlock=[10,11,0,1];
Blockstr={'FixLoc,FixSceneSeq','FixLoc,RandSceneSeq','RandLoc,FixSceneSeq','RandLoc,RandSceneSeq'};
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
if RecordDate<datestr('04,03,2021') && TaskCode==FORAGING_EXPLORE_TASK_ITL
    fix_tmp=BlockType(:,3)>0;
    BlockType(fix_tmp,2)=BlockType(fix_tmp,2)-10;
    
  
end


BlockType=BlockType(Selection,2);

if strcmp(FileName,'Robin040321FORAGEXPLOITL1.rst')
    
    %Correct for trial random 4
    sel=(Selection>320)';
    correct_b=BlockType==11|BlockType==1;
    
    BlockType(sel&correct_b)=BlockType(sel&correct_b)-1;
  %  fix_tmp=BlockType(:,3)>0;
  %  BlockType(fix_tmp,2)=BlockType(fix_tmp,2)-10;
end
if strcmp(FileName,'Robin040721FORAGEXPLOITL1.rst')
    
    %Correct for trial random 4
    sel=(Selection>560)';
    correct_b=BlockType==11|BlockType==1;
    
    BlockType(sel&correct_b)=BlockType(sel&correct_b)-1;
  %  fix_tmp=BlockType(:,3)>0;
  %  BlockType(fix_tmp,2)=BlockType(fix_tmp,2)-10;
end


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

%Automatically chanage the marker

StartCheck=FindOutTime(EventChannel,EventBin,SegEvent(SegEvent==SCENE_ON));
EndCheck=FindOutTime(EventChannel,EventBin,SegEventEnd(SegEvent==SCENE_ON));
EndCheck=EndCheck(:,1);
DiffCheck=EndCheck-StartCheck;
if sum(DiffCheck<200)>ceil(length(StartCheck)/2)
    SegEventEnd(SegEvent==SCENE_ON)=[FIXOFF];%Event end;change from showfix cd to fix off for contigency
    disp('Change to FixOff for FPOn First Mode');
    
end




%Get the information about the extra FP
%
if RecordDate<datestr('10,28,2020') && strcmp(MonkeyName,'Robin')
%Recalibrate eyes according to the extra FP
ParaLib=Data.ParaLib;

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
end



%Get scene and fractals information
SceneName=Data.SceneIndex(Selection,3);
TargetName=Data.TargetIndex(Selection,2:end,2);
TargetIndex=Data.TargetIndex(Selection,2:end,1);
TargetLoc=Data.TargetLoc(Selection,2:end,:);

%Find the repetition
RepIndex=DefineRep(SceneName)';

%TargetNum=Data.TargetNum(Selection,2:end);
TargetNum=sum(~isnan(TargetName),2);
TargetSize=Data.TargetSize(Selection,2:end);
TargetPolarAngle=Data.TargetPolarLoc(Selection,2:end,1);
TargetPolarEccentricity=Data.TargetPolarLoc(Selection,2:end,2);

ChoiceOutcome=Data.ChoiceOutcome(Selection,2:end);

%Select the trials at lease has one case of advance eye in
if ~(TaskCode==SCENEOBJ_NEW_LEARNING_TASK)
AdvanceInTargetNum=sum(EventChannel==ADVEYE_IN,2);
%AdvanceInTargetTrial=AdvanceInTargetNum>0;
AdvanceInTargetTrial=AdvanceInTargetNum==TargetNum;

strAdvanceIn=fprintf('There are %d Advance-Eye-In trials\n',sum(AdvanceInTargetTrial));
else
    AdvanceInTargetTrial=1:TrialNum;
    disp('Not a exploration task');
    
end

%Select Searching time
 FixOffTime=FindOutTime(EventChannel,EventTimeChannel,FIXOFF);
 %TargetOnTime=FixOffTime+300;
 AdvanceEyeInTime=FindOutTime(EventChannel,EventTimeChannel,ADVEYE_IN );
 
 TargetOnTime=FindOutTime(EventChannel,EventTimeChannel,TGTCD);
 OutScreenTime=FindOutTime(EventChannel,EventTimeChannel,EYEOUTSCREEN);
 
 FPOnTime=FindOutTime(EventChannel,EventTimeChannel,SHOWFIXCD);
 %screen out the searching trials
 FixOffTime=FixOffTime(AdvanceInTargetTrial,:);
 OutScreenTime=OutScreenTime(AdvanceInTargetTrial,:);
 FPOnTime=FPOnTime(AdvanceInTargetTrial,:);
 
 FixOffTime=ScreenOutOverScreen(FixOffTime,OutScreenTime);
 
 
 TargetOnTime=TargetOnTime(AdvanceInTargetTrial,:);
 %Remove the col
 NaNCol_Check=sum(isnan(TargetOnTime),1);

 TargetOnTime(:,NaNCol_Check==size(TargetOnTime,1))=[];
 
 SceneName=SceneName(AdvanceInTargetTrial,:);
 
 RepIndex=RepIndex(AdvanceInTargetTrial,:);
 
 EyeChannel_X_L=EyeChannel_X_L(AdvanceInTargetTrial);
 EyeChannel_Y_L=EyeChannel_Y_L(AdvanceInTargetTrial);
 
 TargetLoc=TargetLoc(AdvanceInTargetTrial,:,:);
 TargetPolarAngle=TargetPolarAngle(AdvanceInTargetTrial,:);
 
 TargetNum=TargetNum(AdvanceInTargetTrial);
 BlockType=BlockType(AdvanceInTargetTrial);
 
 [UniqueBlockType,orderindex,reprojectorderindex]=unique(BlockType,'row');
NumberOfBlocks=max(reprojectorderindex);
 
 %If the target onset time and fixofftime dimension are not aligned 
 if size(TargetOnTime,2)~=size(FixOffTime,2)
     FixOffTime=FixOffTime(:,1:size(TargetOnTime,2));
     disp('Warning!FixOffTime Corrected!!')
 end
     
 
SearchingTime= (TargetOnTime-FixOffTime);%/1000;%Change into s for presentation

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

%TargetNumLine
TargetNum_Line=repmat(TargetNum,1,max(TargetNum));
TargetNum_Line=reshape(TargetNum_Line,numel(TargetNum_Line),1);
TargetNum_Line=TargetNum_Line(NoNNaNSelector);

GoodObj_Line=mod(TargetIndex_Line,2)==1;
BadObj_Line=mod(TargetIndex_Line,2)==0;

TargetPolarAngle_Line=TargetPolarAngle;
TargetPolarAngle_Line=reshape(TargetPolarAngle,numel(TargetPolarAngle),1);
TargetPolarAngle_Line=TargetPolarAngle_Line(NoNNaNSelector);


%TargetOrderLine
TargetOrder=1:max(TargetNum);
TotalNum=length(TargetNum);
TargetOrder_Mat=repmat(TargetOrder,TotalNum,1);
TargetOrder_Line=reshape(TargetOrder_Mat,numel(TargetOrder_Mat),1);
TargetOrder_Line=TargetOrder_Line(NoNNaNSelector);


%Rearrange block type into lines
BlockType_Line=repmat(BlockType,1,max(TargetNum));
BlockType_Line=reshape(BlockType_Line,numel(BlockType_Line),1);
BlockType_Line=BlockType_Line(NoNNaNSelector);

%RepIndex line
RepIndex_Line=repmat(RepIndex,1,max(TargetNum));
RepIndex_Line=reshape(RepIndex_Line,numel(RepIndex_Line),1);
RepIndex_Line=RepIndex_Line(NoNNaNSelector);

%Rearrange searching time into lines
SearchingTime_Line=reshape(SearchingTime,numel(SearchingTime),1);
SearchingTime_Line=SearchingTime_Line(NoNNaNSelector);

%Rearrange_trialIndex into lines
TrialIndex=[1:size(SearchingTime,1)]';
TrialIndex_Line=repmat(TrialIndex,1,max(TargetNum));
TrialIndex_Line=reshape(TrialIndex_Line,numel(TrialIndex_Line),1);
TrialIndex_Line=TrialIndex_Line(NoNNaNSelector);



UniqueSceneName=unique(SceneName);



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
MeanPSTHEach=SegmentData.MeanPSTHEach;%=[PSTH_Event_Time',PSTH_Event'];

%First measure the scene responses change

MeanPSTH_SceneOnTime=MeanPSTH{1,1};
MeanPSTH_SceneOnFR=MeanPSTH{1,2};

PSTH_SceneOnTime_Mean=nanmean(MeanPSTH_SceneOnTime,1);


MeanPSTH_SceneOffTime=MeanPSTH{2,1};
MeanPSTH_SceneOffFR=MeanPSTH{2,2};

PSTH_SceneOffTime_Mean=nanmean(MeanPSTH_SceneOffTime,1);

%Relative time bin for each event
RelativeTimeMarkerSceneOn=SegmentData.RelativeTimeMarker{1};%=EventTimeMarker;
%Absolute time marker for each event
AbosoluteTimeMarkerSceneOn=SegmentData.AbsoluteTimeMarker{1};
SceneOnMarker=min(RelativeTimeMarkerSceneOn,[],1);



RelativeTimeMarkerSceneOff=SegmentData.RelativeTimeMarker{2};%=EventTimeMarker;
%Absolute time marker for each event
AbosoluteTimeMarkerSceneOff=SegmentData.AbsoluteTimeMarker{2};
SceneOffMarker=min(RelativeTimeMarkerSceneOff,[],1);

%Object responses
MeanPSTH_ObjOnTime=MeanPSTH{3,1};
MeanPSTH_ObjOnFR=MeanPSTH{3,2};

PSTH_ObjOnTime_Mean=nanmean(MeanPSTH_ObjOnTime,1);

%Relative time bin for each event
RelativeTimeMarkerObjOn=SegmentData.RelativeTimeMarker{3};%=EventTimeMarker;
%Absolute time marker for each event
AbosoluteTimeMarkerObjOn=SegmentData.AbsoluteTimeMarker{3};
ObjOnMarker=min(RelativeTimeMarkerObjOn,[],1);


MeanPSTH_ObjOffTime=MeanPSTH{4,1};
MeanPSTH_ObjOffFR=MeanPSTH{4,2};

PSTH_ObjOffTime_Mean=nanmean(MeanPSTH_ObjOffTime,1);

%Relative time bin for each event
RelativeTimeMarkerObjOff=SegmentData.RelativeTimeMarker{4};%=EventTimeMarker;
%Absolute time marker for each event
AbosoluteTimeMarkerObjOff=SegmentData.AbsoluteTimeMarker{4};
ObjOffMarker=min(RelativeTimeMarkerObjOff,[],1);


SearchingTimeFirst=SearchingTime(:,1)/1000;


%{
%Raster and raw spike count PSTH
AlignedRaster=SegmentData.RasterAfterAlign{1};
PSTH_Event_SpikeCount=SegmentData.PSTH_Event_SpikeCount{1};
%}


%Average the firing rate for each interval
%First scene on
%MeanPSTH_SceneOnFR
SceneOnInterval=[0,1000];% 1000 ms after scene onset
Sel_SceneOnInterval=PSTH_SceneOnTime_Mean>SceneOnInterval(1) & PSTH_SceneOnTime_Mean<SceneOnInterval(2);

FR_SceneOn_Ave=nanmean(MeanPSTH_SceneOnFR(:,Sel_SceneOnInterval),2);

PreSceneOnInterval=[-500,0];% 500 ms before scene onset
Sel_PreSceneOnInterval=PSTH_SceneOnTime_Mean>PreSceneOnInterval(1) & PSTH_SceneOnTime_Mean<PreSceneOnInterval(2);

FR_PreSceneOn_Ave=nanmean(MeanPSTH_SceneOnFR(:,Sel_PreSceneOnInterval),2);

%Loop over blocks 
EveryRep=5;%Average every EveryRep repetitions
UniqueBlock=unique(BlockType);



for i=1:length(UniqueBlock)
    currBlock=UniqueBlock(i);
    RepIndex_CurrBlock=RepIndex(BlockType==currBlock);
    MinRep=min(RepIndex_CurrBlock);
    MaxRep=max(RepIndex_CurrBlock);
    
    SceneOnPSTH_currBlock=MeanPSTH_SceneOnFR(BlockType==currBlock,:);
    SceneOffPSTH_currBlock=MeanPSTH_SceneOffFR(BlockType==currBlock,:);
    
    ObjOnPSTH_currBlock=MeanPSTH_ObjOnFR(BlockType_Line==currBlock,:);
   
    ObjOffPSTH_currBlock=MeanPSTH_ObjOffFR(BlockType_Line==currBlock,:);
   
    
    RepIndex_CurrBlock_Line=RepIndex_Line(BlockType_Line==currBlock);
    
    SceneCurrBlock=SceneName(BlockType==currBlock);
    UnqiueSceneCurrBlock=unique(SceneCurrBlock);
    
    RepIndex_Block{i}=DefineRep(SceneCurrBlock)';
    RepIndex_Block_tmp=DefineRep(SceneCurrBlock)';
    
    
    ObjectCurrBlock=TargetName_Line(BlockType_Line==currBlock);
    UnqiueObjectCurrBlock=unique(ObjectCurrBlock);
    
    FR_SceneOn_Ave_currBlock=FR_SceneOn_Ave(BlockType==currBlock);
    
    AveFR_Rep_SceneOnBlock_Mean_tmp=[];
    AveFR_Rep_SceneOnBlock_Sem_tmp=[];
    
    %Loop over each rep
     for m=1:length(unique(RepIndex_Block_tmp))
         sel_tmp=RepIndex_Block_tmp==RepIndex_Block_tmp(m);
         %Correction for too few number
         %{
         if sum(sel_tmp)<3
             AveFR_Rep_SceneOnBlock_Mean_tmp(m)=NaN;
             AveFR_Rep_SceneOnBlock_Sem_tmp(m)=NaN;
         else
         %}
         AveFR_Rep_SceneOnBlock_Mean_tmp(m)=nanmean(FR_SceneOn_Ave_currBlock(sel_tmp));
         AveFR_Rep_SceneOnBlock_Sem_tmp(m)=nanstd(FR_SceneOn_Ave_currBlock(sel_tmp))/sqrt(sum(~isnan(FR_SceneOn_Ave_currBlock(sel_tmp))));
        % end
         
     end
    
    AveFR_Rep_SceneOnBlock_Mean{i}=AveFR_Rep_SceneOnBlock_Mean_tmp;
    AveFR_Rep_SceneOnBlock_Sem{i}=AveFR_Rep_SceneOnBlock_Sem_tmp;
    
    
    %First average over all scenes 
    indexflag=0;
    for j=MinRep:EveryRep:MaxRep
        RepSel=j:min(j+EveryRep-1,MaxRep);
        Sel_currBlockRep=ismember(RepIndex_CurrBlock,RepSel);
        indexflag=indexflag+1;
        %Scene responses 
        AllSceneOnRepBlock=SceneOnPSTH_currBlock(Sel_currBlockRep,:);
        AllSceneOnRepBlock_Mean{i,indexflag}=nanmean(AllSceneOnRepBlock,1);
        AllSceneOnRepBlock_Sem{i,indexflag}=nanstd(AllSceneOnRepBlock,[],1)./sqrt(sum(~isnan(AllSceneOnRepBlock),1));

        AllSceneOffRepBlock=SceneOffPSTH_currBlock(Sel_currBlockRep,:);
        AllSceneOffRepBlock_Mean{i,indexflag}=nanmean(AllSceneOffRepBlock,1);
        AllSceneOffRepBlock_Sem{i,indexflag}=nanstd(AllSceneOffRepBlock,[],1)./sqrt(sum(~isnan(AllSceneOffRepBlock),1));
        
        
        RepStore{i,indexflag}=num2str(j);
        %Object responses
        Sel_currBlockRepLine=ismember(RepIndex_CurrBlock_Line,RepSel);
        AllObjOnRepBlock=ObjOnPSTH_currBlock(Sel_currBlockRepLine,:);
        AllObjOnRepBlock_Mean{i,indexflag}=nanmean(AllObjOnRepBlock,1);
        AllObjOnRepBlock_Sem{i,indexflag}=nanstd(AllObjOnRepBlock,[],1)./sqrt(sum(~isnan(AllObjOnRepBlock),1));
        
        AllObjOffRepBlock=ObjOffPSTH_currBlock(Sel_currBlockRepLine,:);
        AllObjOffRepBlock_Mean{i,indexflag}=nanmean(AllObjOffRepBlock,1);
        AllObjOffRepBlock_Sem{i,indexflag}=nanstd(AllObjOffRepBlock,[],1)./sqrt(sum(~isnan(AllObjOffRepBlock),1));
        
        
        %Select some interval as an average 
        
        %Loop over each scene 
        for m=1:length(UnqiueSceneCurrBlock)
            scenesel=SceneCurrBlock==UnqiueSceneCurrBlock(m);
            
            sel_withrep= Sel_currBlockRep &scenesel;
            
            EachSceneOnRepBlock=SceneOnPSTH_currBlock(sel_withrep,:);
            EachSceneOnRepBlock_Mean{i,m,indexflag}=nanmean(EachSceneOnRepBlock,1);
            EachSceneOnRepBlock_Sem{i,m,indexflag}=nanstd(EachSceneOnRepBlock,[],1)./sqrt(sum(~isnan(EachSceneOnRepBlock),1));
            
            EachSceneOffRepBlock=SceneOffPSTH_currBlock(sel_withrep,:);
            EachSceneOffRepBlock_Mean{i,m,indexflag}=nanmean(EachSceneOffRepBlock,1);
            EachSceneOffRepBlock_Sem{i,m,indexflag}=nanstd(EachSceneOffRepBlock,[],1)./sqrt(sum(~isnan(EachSceneOffRepBlock),1));
            
            SceneNameRepBlock(i,m)=UnqiueSceneCurrBlock(m);
        end
        
        %Loop over each object
        for m=1:length(UnqiueObjectCurrBlock)
            objsel=ObjectCurrBlock==UnqiueObjectCurrBlock(m);
            
            sel_withrep= Sel_currBlockRepLine & objsel;
            
            EachObjOnRepBlock=ObjOnPSTH_currBlock(sel_withrep,:);
            EachObjOnRepBlock_Mean{i,m,indexflag}=nanmean(EachObjOnRepBlock,1);
            EachObjOnRepBlock_Sem{i,m,indexflag}=nanstd(EachObjOnRepBlock,[],1)./sqrt(sum(~isnan(EachObjOnRepBlock),1));
            
            EachObjOffRepBlock=ObjOffPSTH_currBlock(sel_withrep,:);
            EachObjOffRepBlock_Mean{i,m,indexflag}=nanmean(EachObjOffRepBlock,1);
            EachObjOffRepBlock_Sem{i,m,indexflag}=nanstd(EachObjOffRepBlock,[],1)./sqrt(sum(~isnan(EachObjOffRepBlock),1));
            
            ObjNameRepBlock(i,m)=UnqiueObjectCurrBlock(m);
        end
        
        
        
        
    end
    
    
end


%PSTH_SceneOnTime_Mean

%%
FigureStartNum=57;
FigureIndex=1;
if ShowFigureFlag
    %%FR_SceneOn_Ave
 
 
 %%%%%%%%%%%%%%%%%%Next figure: Raw Average Response Change with time
 
 figtitlestr{FigureIndex}='AverageSceneResponses';
 fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[30,100, 1500,600],'Name',figtitlestr{FigureIndex});  
 
 TrialIndex=1:length(FR_SceneOn_Ave);
 %Scene onset
 subplot(3,2,1);
 for i=1:length(UniqueSceneName)
     sel=SceneName==UniqueSceneName(i);
      choosecolor=mod(i,length(ColorRep));
         if choosecolor==0
             choosecolor=1;
         end
         
     plot([TrialIndex(sel)],[FR_SceneOn_Ave(sel)'],'-or','Color',ColorRep{choosecolor},'LineWidth',3);
     hold on
     
     
     

     
     
 end
 %Add block type
 UniqueBlockType=unique(BlockType);
 StartPlotIndex=-1;
for i=1:length(UniqueBlockType)
    sel=BlockType==UniqueBlockType(i);
    order(i)=find(AssociateBlock==UniqueBlockType(i));
    colorhere=find(AssociateBlock==UniqueBlockType(i));
 plot([TrialIndex(sel)],[StartPlotIndex*ones(1,sum(sel))],'s','Color',ColorForEachMode(colorhere,:));
 
 StartPlotIndex=StartPlotIndex-1;
 hold on
end
if LegendOn
legend([cellstr(string(num2str(UniqueSceneName)));Blockstr(order)']);
end
xlabel('Trial Number');
ylabel('Average FR of scene onset');
title(' Scene onset')
 set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');

 %PreScene onset
 subplot(3,2,2);
 
 for i=1:length(UniqueSceneName)
     sel=SceneName==UniqueSceneName(i);
      choosecolor=mod(i,length(ColorRep));
         if choosecolor==0
             choosecolor=1;
         end
         
     plot([TrialIndex(sel)],[FR_PreSceneOn_Ave(sel)'],'-or','Color',ColorRep{choosecolor},'LineWidth',3);
     hold on
     
     
     

     
     
 end
 %Add block type
 UniqueBlockType=unique(BlockType);
 StartPlotIndex=-1;
for i=1:length(UniqueBlockType)
    sel=BlockType==UniqueBlockType(i);
    order(i)=find(AssociateBlock==UniqueBlockType(i));
    colorhere=find(AssociateBlock==UniqueBlockType(i));
 plot([TrialIndex(sel)],[StartPlotIndex*ones(1,sum(sel))],'s','Color',ColorForEachMode(colorhere,:));
 
 StartPlotIndex=StartPlotIndex-1;
 hold on
end
if LegendOn
legend([cellstr(string(num2str(UniqueSceneName)));Blockstr(order)']);
end
xlabel('Trial Number');
ylabel('Average FR before scene onset');
title('500ms Before Scene onset');
 set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
 
 
 %RepIndex
 subplot(3,2,3);%Align each scene together
 
 for i=1:length(UniqueSceneName)
     sel=SceneName==UniqueSceneName(i);
      choosecolor=mod(i,length(ColorRep));
         if choosecolor==0
             choosecolor=1;
         end
         
     plot([RepIndex(sel)],[FR_SceneOn_Ave(sel)'],'-or','Color',ColorRep{choosecolor},'LineWidth',3);
     hold on
     
     
     

     
     
 end
 %Add block type
 UniqueBlockType=unique(BlockType);
 StartPlotIndex=-1;
for i=1:length(UniqueBlockType)
    sel=BlockType==UniqueBlockType(i);
    order(i)=find(AssociateBlock==UniqueBlockType(i));
    colorhere=find(AssociateBlock==UniqueBlockType(i));
 plot([RepIndex(sel)],[StartPlotIndex*ones(1,sum(sel))],'s','Color',ColorForEachMode(colorhere,:));
 
 StartPlotIndex=StartPlotIndex-1;
 hold on
end
if LegendOn
legend([cellstr(string(num2str(UniqueSceneName)));Blockstr(order)']);
end

xlabel('Repetition Number');
ylabel('Average FR of scene onset');
title(' Scene onset')
 set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
 
 
  %RepIndex
 subplot(3,2,4);%Align each pre scene together
 
 for i=1:length(UniqueSceneName)
     sel=SceneName==UniqueSceneName(i);
      choosecolor=mod(i,length(ColorRep));
         if choosecolor==0
             choosecolor=1;
         end
         
     plot([RepIndex(sel)],[FR_PreSceneOn_Ave(sel)'],'-or','Color',ColorRep{choosecolor},'LineWidth',3);
     hold on
     
     
     

     
     
 end
 %Add block type
 UniqueBlockType=unique(BlockType);
 StartPlotIndex=-1;
for i=1:length(UniqueBlockType)
    sel=BlockType==UniqueBlockType(i);
    order(i)=find(AssociateBlock==UniqueBlockType(i));
    colorhere=find(AssociateBlock==UniqueBlockType(i));
 plot([RepIndex(sel)],[StartPlotIndex*ones(1,sum(sel))],'s','Color',ColorForEachMode(colorhere,:));
 
 StartPlotIndex=StartPlotIndex-1;
 hold on
end
if LegendOn
legend([cellstr(string(num2str(UniqueSceneName)));Blockstr(order)']);
end

xlabel('Repetition Number');
ylabel('Average FR of scene onset');
title(' Scene onset');
 set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
 
 
 %Behavior: Searching time change
 subplot(3,2,5);%Raw trial number
 for i=1:length(UniqueSceneName)
     sel=SceneName==UniqueSceneName(i);
      choosecolor=mod(i,length(ColorRep));
         if choosecolor==0
             choosecolor=1;
         end
         
    % semilogy([TrialIndex(sel)],[SearchingTimeFirst(sel)'],'-or','Color',ColorRep{choosecolor},'LineWidth',3);
    plot([TrialIndex(sel)],[SearchingTimeFirst(sel)'],'-or','Color',ColorRep{choosecolor},'LineWidth',3);
     hold on
     
     
     

     
     
 end
 %Add block type
 UniqueBlockType=unique(BlockType);
 StartPlotIndex=-1;
for i=1:length(UniqueBlockType)
    sel=BlockType==UniqueBlockType(i);
    order(i)=find(AssociateBlock==UniqueBlockType(i));
    colorhere=find(AssociateBlock==UniqueBlockType(i));
 plot([TrialIndex(sel)],[StartPlotIndex*ones(1,sum(sel))],'s','Color',ColorForEachMode(colorhere,:));
 
 StartPlotIndex=StartPlotIndex-1;
 hold on
end
if LegendOn
legend([cellstr(string(num2str(UniqueSceneName)));Blockstr(order)']);
end

xlabel('Trial Number');
ylabel('First Rearching Time(s)');
title(' Behavior');
 set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
 
 
 %Behavior: Searching time change
 subplot(3,2,6);%Repetition
 for i=1:length(UniqueSceneName)
     sel=SceneName==UniqueSceneName(i);
      choosecolor=mod(i,length(ColorRep));
         if choosecolor==0
             choosecolor=1;
         end
         
   %  semilogy([RepIndex(sel)],[SearchingTimeFirst(sel)'],'-or','Color',ColorRep{choosecolor},'LineWidth',3);
      plot([RepIndex(sel)],[SearchingTimeFirst(sel)'],'-or','Color',ColorRep{choosecolor},'LineWidth',3);
     hold on
     
     
     

     
     
 end
 %Add block type
 UniqueBlockType=unique(BlockType);
 StartPlotIndex=-1;
for i=1:length(UniqueBlockType)
    sel=BlockType==UniqueBlockType(i);
    order(i)=find(AssociateBlock==UniqueBlockType(i));
    colorhere=find(AssociateBlock==UniqueBlockType(i));
plot([RepIndex(sel)],[StartPlotIndex*ones(1,sum(sel))],'s','Color',ColorForEachMode(colorhere,:));
 
 StartPlotIndex=StartPlotIndex-1;
 hold on
end
if LegendOn
legend([cellstr(string(num2str(UniqueSceneName)));Blockstr(order)']);
end 
xlabel('Repetition Number');
ylabel('First Rearching Time(s)');
title(' Behavior');
 set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
 
  data_msg=strcat(FileName,' SpikeChannel:',num2str(SpikeChannel),';ChannelID:',num2str(ChannelID),';ClusterID:',num2str(ClusterID),';TrialNum: ',num2str(TrialNum));

    annotation('textbox', [0.17, 0.95, 0.7, 0.04], 'string', data_msg,'FontSize',20,'FontWeight','Bold','LineStyle','none');
    
 
 FigureIndex=FigureIndex+1;
FigureStartNum=FigureStartNum+1;  

%%%%%%%%%%%%%%%%%%Next figure: Raw Average Response Change with time for
%%%%%%%%%%%%%%%%%%each block

figtitlestr{FigureIndex}='AverageSceneResponsesChangeForEachBlock';
 fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[50,100, 1500,600],'Name',figtitlestr{FigureIndex});  
 
 TrialIndex=1:length(FR_SceneOn_Ave);
 %Scene onset
 %subplot(2,2,1);
 
 
 for i=1:length(UniqueBlockType)
     subplot(3,length(UniqueBlockType),i);
     for j=1:length(UniqueSceneName)
         sel=BlockType==UniqueBlockType(i) & SceneName==UniqueSceneName(j);
         sel_scene=SceneName(BlockType==UniqueBlockType(i))==UniqueSceneName(j);
        choosecolor=mod(j,length(ColorRep));
         if choosecolor==0
             choosecolor=1;
         end
         
         
     plot([RepIndex_Block{i}(sel_scene)],[FR_SceneOn_Ave(sel)'],'-or','Color',ColorRep{choosecolor},'LineWidth',3);
     hold on
     end
     
     xlabel('Repetition Number');
     ylabel('Average FR of scene onset');
     title(Blockstr(find(AssociateBlock==UniqueBlockType(i))));
     set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
     
 end
 
 %PreScene On
 for i=1:length(UniqueBlockType)
     subplot(3,length(UniqueBlockType),length(UniqueBlockType)+i);
     for j=1:length(UniqueSceneName)
         sel=BlockType==UniqueBlockType(i) & SceneName==UniqueSceneName(j);
         sel_scene=SceneName(BlockType==UniqueBlockType(i))==UniqueSceneName(j);
        choosecolor=mod(j,length(ColorRep));
         if choosecolor==0
             choosecolor=1;
         end
         
         
     plot([RepIndex_Block{i}(sel_scene)],[FR_PreSceneOn_Ave(sel)'],'-or','Color',ColorRep{choosecolor},'LineWidth',3);
     hold on
     end
     
     xlabel('Repetition Number');
     ylabel('Average FR before scene onset');
     title(Blockstr(find(AssociateBlock==UniqueBlockType(i))));
     set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
     
 end
 
 %Behavior
 for i=1:length(UniqueBlockType)
     subplot(3,length(UniqueBlockType),length(UniqueBlockType)*2+i);
     for j=1:length(UniqueSceneName)
         sel=BlockType==UniqueBlockType(i) & SceneName==UniqueSceneName(j);
         sel_scene=SceneName(BlockType==UniqueBlockType(i))==UniqueSceneName(j);
        choosecolor=mod(j,length(ColorRep));
         if choosecolor==0
             choosecolor=1;
         end
         
         
     plot([RepIndex_Block{i}(sel_scene)],[SearchingTimeFirst(sel)'],'-or','Color',ColorRep{choosecolor},'LineWidth',3);
     hold on
     end
     
     xlabel('Repetition Number');
     ylabel('SearchingTime(s)');
     title(Blockstr(find(AssociateBlock==UniqueBlockType(i))));
     set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
     
 end
 
 
  data_msg=strcat(FileName,' SpikeChannel:',num2str(SpikeChannel),';ChannelID:',num2str(ChannelID),';ClusterID:',num2str(ClusterID),';TrialNum: ',num2str(TrialNum));

    annotation('textbox', [0.17, 0.95, 0.7, 0.04], 'string', data_msg,'FontSize',20,'FontWeight','Bold','LineStyle','none');
    
 FigureIndex=FigureIndex+1;
FigureStartNum=FigureStartNum+1;  
 
%%%%%%%%%%%%%%%%%%Next figure: Mean of Average Response Change with time for
%%%%%%%%%%%%%%%%%%each block, compared to behavior

figtitlestr{FigureIndex}='AverageSceneResponsesChangeForEachBlock';
 fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[50,100, 1500,600],'Name',figtitlestr{FigureIndex});  
 for i=1:length(UniqueBlock)
     errorbar([1:length(AveFR_Rep_SceneOnBlock_Mean{i})],AveFR_Rep_SceneOnBlock_Mean{i},AveFR_Rep_SceneOnBlock_Sem{i},'-or','Color',ColorForEachMode(find(AssociateBlock==UniqueBlock(i)),:),'LineWidth',3);
     
     hold on
 end
 xlabel('Repetition Number');
     ylabel('Average Scene On responses');
     title('Scene On Responses');
    if LegendOn
        legend(Blockstr(find(AssociateBlock==UniqueBlock(i))));
    end
     
 set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
 
 
 FigureIndex=FigureIndex+1;
FigureStartNum=FigureStartNum+1;  
%AveFR_Rep_SceneOnBlock_Mean


    %%%%%%%%%%%%%%%%%%Next figure: Scene on PSTH for all the divied
    %%%%%%%%%%%%%%%%%%repetitions

 figtitlestr{FigureIndex}='SceneOnOffForRepEachBlock';
 fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[70,100, 1500,800],'Name',figtitlestr{FigureIndex});  
 for i=1:length(UniqueBlock)
     subplot(ceil(length(UniqueBlock)),3,(i-1)*3+1);
     maxval=0;
     for j=1:size(AllSceneOnRepBlock_Mean,2)
         currMean=AllSceneOnRepBlock_Mean{i,j};
         currSem=AllSceneOnRepBlock_Sem{i,j};
         
         
         choosecolor=mod(j,length(ColorRep));
         if choosecolor==0
             choosecolor=1;
         end
         if ~isempty(currMean)
          shadedErrorBar(PSTH_SceneOnTime_Mean,currMean,currSem,'lineprops',{ColorRep{choosecolor}},'transparent',1,'patchSaturation',0.5);
         end
           set(findall(gca, 'Type', 'Line'),'LineWidth',3);
           hold on
     
         max_tmp=max(currMean+currSem);
         if max_tmp>maxval
             maxval=max_tmp;
         end
         
     end
     
      h1=plot([SceneOnMarker(1), SceneOnMarker(1)],[0, maxval],'--k','LineWidth',1.5);
      h2=plot([SceneOnMarker(2), SceneOnMarker(2)],[0, maxval],'--k','LineWidth',1.5);
  
     set( get( get( h1, 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );
     set( get( get( h2, 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );
     
     currLegend=RepStore(i,:);
     numin=cellfun(@numel,currLegend);
     currLegend=currLegend(numin>0);
if LegendOn
     legend(currLegend);
end
     title(Blockstr{find(AssociateBlock==UniqueBlock(i))});
     xlabel('Time from scene onset(ms)');
      ylabel('Firing Rate(Hz)');
  
     set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
      xlim([-PlotStartOffset(1),PlotStopOffset(1)+SceneOnMarker(2)]);
     
      %Scene Off
     subplot(ceil(length(UniqueBlock)),3,(i-1)*3+2);
     maxval=0;
     for j=1:size(AllSceneOffRepBlock_Mean,2)
         currMean=AllSceneOffRepBlock_Mean{i,j};
         currSem=AllSceneOffRepBlock_Sem{i,j};
         
         
         choosecolor=mod(j,length(ColorRep));
         if choosecolor==0
             choosecolor=1;
         end
         if ~isempty(currMean)
          shadedErrorBar(PSTH_SceneOffTime_Mean,currMean,currSem,'lineprops',{ColorRep{choosecolor}},'transparent',1,'patchSaturation',0.5);
         end
           set(findall(gca, 'Type', 'Line'),'LineWidth',3);
           hold on
    
          max_tmp=max(currMean+currSem);
         if max_tmp>maxval
             maxval=max_tmp;
         end
         
     end
     
      h1=plot([SceneOffMarker(1), SceneOffMarker(1)],[0,maxval],'--k','LineWidth',1.5);
      h2=plot([SceneOffMarker(2), SceneOffMarker(2)],[0,maxval],'--k','LineWidth',1.5);
  
     set( get( get( h1, 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );
     set( get( get( h2, 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );

     currLegend=RepStore(i,:);
     numin=cellfun(@numel,currLegend);
     currLegend=currLegend(numin>0);
if LegendOn
     legend(currLegend);
end
     title(Blockstr{find(AssociateBlock==UniqueBlock(i))});
     xlabel('Time from scene offset(ms)');
      ylabel('Firing Rate(Hz)');
  
     set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
      xlim([-PlotStartOffset(2),PlotStopOffset(2)*1.5]);
      
  %Object On   
subplot(ceil(length(UniqueBlock)),3,(i-1)*3+3);

maxval=0;
     for j=1:size(AllObjOnRepBlock_Mean,2)
         currMean=AllObjOnRepBlock_Mean{i,j};
         currSem=AllObjOnRepBlock_Sem{i,j};
         
         
         choosecolor=mod(j,length(ColorRep));
         if choosecolor==0
             choosecolor=1;
         end
         if ~isempty(currMean)
          shadedErrorBar(PSTH_ObjOnTime_Mean,currMean,currSem,'lineprops',{ColorRep{choosecolor}},'transparent',1,'patchSaturation',0.5);
         end
           set(findall(gca, 'Type', 'Line'),'LineWidth',3);
           hold on
     
         max_tmp=max(currMean+currSem);
         if max_tmp>maxval
             maxval=max_tmp;
         end
         
     end
     
      h1=plot([ObjOnMarker(1), ObjOnMarker(1)],[0, maxval],'--k','LineWidth',1.5);
      h2=plot([ObjOnMarker(2), ObjOnMarker(2)],[0, maxval],'--k','LineWidth',1.5);
  
     set( get( get( h1, 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );
     set( get( get( h2, 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );
     
     currLegend=RepStore(i,:);
     numin=cellfun(@numel,currLegend);
     currLegend=currLegend(numin>0);
     if LegendOn
     legend(currLegend);
     end
     title(Blockstr{find(AssociateBlock==UniqueBlock(i))});
     xlabel('Time from object onset(ms)');
      ylabel('Firing Rate(Hz)');
  
     set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
      xlim([-PlotStartOffset(3),PlotStopOffset(3)+ObjOnMarker(2)]);

%{
      %Remove off, not usefull
%Object Off   
subplot(ceil(length(UniqueBlock)),4,(i-1)*4+4);

maxval=0;
     for j=1:size(AllObjOffRepBlock_Mean,2)
         currMean=AllObjOffRepBlock_Mean{i,j};
         currSem=AllObjOffRepBlock_Sem{i,j};
         
         
         choosecolor=mod(j,length(ColorRep));
         if choosecolor==0
             choosecolor=1;
         end
         if ~isempty(currMean)
          shadedErrorBar(PSTH_ObjOffTime_Mean,currMean,currSem,'lineprops',{ColorRep{choosecolor}},'transparent',1,'patchSaturation',0.5);
         end
           set(findall(gca, 'Type', 'Line'),'LineWidth',3);
           hold on
     
         max_tmp=max(currMean+currSem);
         if max_tmp>maxval
             maxval=max_tmp;
         end
         
     end
     
      h1=plot([ObjOffMarker(1), ObjOffMarker(1)],[0, maxval],'--k','LineWidth',1.5);
      h2=plot([ObjOffMarker(2), ObjOffMarker(2)],[0, maxval],'--k','LineWidth',1.5);
  
     set( get( get( h1, 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );
     set( get( get( h2, 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );
     
     currLegend=RepStore(i,:);
     numin=cellfun(@numel,currLegend);
     currLegend=currLegend(numin>0);

     legend(currLegend);
     title(Blockstr{find(AssociateBlock==UniqueBlock(i))});
     xlabel('Time from object offset(ms)');
      ylabel('Firing Rate(Hz)');
  
     set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
      xlim([-PlotStartOffset(4),PlotStopOffset(4)+ObjOffMarker(2)]);

%}
     
 end
 
   data_msg=strcat(FileName,' SpikeChannel:',num2str(SpikeChannel),';ChannelID:',num2str(ChannelID),';ClusterID:',num2str(ClusterID),';TrialNum: ',num2str(TrialNum));

    annotation('textbox', [0.17, 0.95, 0.7, 0.04], 'string', data_msg,'FontSize',20,'FontWeight','Bold','LineStyle','none');
 

 
FigureIndex=FigureIndex+1;
FigureStartNum=FigureStartNum+1;  

if EachSceneObjOn

%%%%%%%%%%%%%%%%%%Next figure: Loop over each scene for each block
 for i=1:length(UniqueBlock)


 figtitlestr{FigureIndex}=strcat('SceneOnOffForRepBlock((',Blockstr{find(AssociateBlock==UniqueBlock(i))},'))');
 
 fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[50,100, 1500,800],'Name',figtitlestr{FigureIndex});  
 
 dataonmean=squeeze(EachSceneOnRepBlock_Mean(i,:,:));
 dataonsem=squeeze(EachSceneOnRepBlock_Sem(i,:,:));
 
 dataoffmean=squeeze(EachSceneOffRepBlock_Mean(i,:,:));
 dataoffsem=squeeze(EachSceneOffRepBlock_Sem(i,:,:));
 
 
 
 
 if size(dataonmean,1)==1
     dataonmean=dataonmean';
     dataonsem=dataonsem';
     
 end
 if size(dataoffmean,1)==1
     dataoffmean=dataoffmean';
     dataoffsem=dataoffsem';
     
 end
 
 
 
 Scenecurr=SceneNameRepBlock(i,:);
%Loop over the scene
for m=1:length(Scenecurr)
    %Put each scene on each figure
    subplot(ceil(length(Scenecurr)/2),4,(m-1)*2+1);
    maxval=0;
     for j=1:size(dataonmean,2)
         
         currMean=dataonmean{m,j};
         currSem=dataonsem{m,j};
         
         
         
         choosecolor=mod(j,length(ColorRep));
         if choosecolor==0
             choosecolor=1;
         end
         if ~isempty(currMean)
          shadedErrorBar(PSTH_SceneOnTime_Mean,currMean,currSem,'lineprops',{ColorRep{choosecolor}},'transparent',1,'patchSaturation',0.5);
         end
           set(findall(gca, 'Type', 'Line'),'LineWidth',3);
           hold on
     
         max_tmp=max(currMean+currSem);
         if max_tmp>maxval
             maxval=max_tmp;
         end
         
     end
     
      h1=plot([SceneOnMarker(1), SceneOnMarker(1)],[0, maxval],'--k','LineWidth',1.5);
      h2=plot([SceneOnMarker(2), SceneOnMarker(2)],[0, maxval],'--k','LineWidth',1.5);
  
     set( get( get( h1, 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );
     set( get( get( h2, 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );
     
     currLegend=RepStore(i,:);
     numin=cellfun(@numel,currLegend);
     currLegend=currLegend(numin>0);
if LegendOn
     legend(currLegend);
end
     title(num2str(Scenecurr(m)));
     xlabel('Time from scene onset(ms)');
      ylabel('Firing Rate(Hz)');
  
     set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
      xlim([-PlotStartOffset(1),PlotStopOffset(1)+SceneOnMarker(2)]);
     %Scene off 
     subplot(ceil(length(Scenecurr)/2),4,(m-1)*2+2);
     for j=1:size(dataoffmean,2)
         
         currMean=dataoffmean{m,j};
         currSem=dataoffsem{m,j};
         
         
         
         choosecolor=mod(j,length(ColorRep));
         if choosecolor==0
             choosecolor=1;
         end
         if ~isempty(currMean)
          shadedErrorBar(PSTH_SceneOffTime_Mean,currMean,currSem,'lineprops',{ColorRep{choosecolor}},'transparent',1,'patchSaturation',0.5);
         end
           set(findall(gca, 'Type', 'Line'),'LineWidth',3);
           hold on
     
         max_tmp=max(currMean+currSem);
         if max_tmp>maxval
             maxval=max_tmp;
         end
         
     end
     
      h1=plot([SceneOffMarker(1), SceneOffMarker(1)],[0, maxval],'--k','LineWidth',1.5);
      h2=plot([SceneOffMarker(2), SceneOffMarker(2)],[0, maxval],'--k','LineWidth',1.5);
  
     set( get( get( h1, 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );
     set( get( get( h2, 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );
     
     currLegend=RepStore(i,:);
     numin=cellfun(@numel,currLegend);
     currLegend=currLegend(numin>0);
if LegendOn
     legend(currLegend);
end
     title(num2str(Scenecurr(m)));
     xlabel('Time from scene offset(ms)');
      ylabel('Firing Rate(Hz)');
  
     set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
      xlim([-PlotStartOffset(2),PlotStopOffset(2)+SceneOffMarker(2)]);
    
    
    
    
end


%EachSceneOnRepBlock_Mean{i,m,indexflag}=nanmean(EachSceneOnRepBlock,1);
%            EachSceneOnRepBlock_Sem{i,m,indexflag}=nanstd(EachSceneOnRepBlock,[],1)./sqrt(sum(~isnan(EachSceneOnRepBlock),1));



FigureIndex=FigureIndex+1;
FigureStartNum=FigureStartNum+1;  

 end% Loop over the block

 
 
 %%%%%%%%%%%%%%%%%%Next figure: Loop over each object for each block
 for i=1:length(UniqueBlock)


 figtitlestr{FigureIndex}=strcat('ObjectOnOffForRepBlock((',Blockstr{find(AssociateBlock==UniqueBlock(i))},'))');
 
 fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[70,100, 1500,800],'Name',figtitlestr{FigureIndex});  
 
 dataonmean=squeeze(EachObjOnRepBlock_Mean(i,:,:));
 dataonsem=squeeze(EachObjOnRepBlock_Sem(i,:,:));
 
 
 
 
 
 
 if size(dataonmean,1)==1
     dataonmean=dataonmean';
     dataonsem=dataonsem';
     
 end
 
 
 
 
 Objcurr=ObjNameRepBlock(i,:);
%Loop over the scene
for m=1:length(Objcurr)
    %Put each scene on each figure
    subplot(ceil(length(Objcurr)/4),4,m);
    maxval=0;
     for j=1:size(dataonmean,2)
         
         currMean=dataonmean{m,j};
         currSem=dataonsem{m,j};
         
         
         
         choosecolor=mod(j,length(ColorRep));
         if choosecolor==0
             choosecolor=1;
         end
         if ~isempty(currMean)
          shadedErrorBar(PSTH_ObjOnTime_Mean,currMean,currSem,'lineprops',{ColorRep{choosecolor}},'transparent',1,'patchSaturation',0.5);
         end
           set(findall(gca, 'Type', 'Line'),'LineWidth',3);
           hold on
     
         max_tmp=max(currMean+currSem);
         if max_tmp>maxval
             maxval=max_tmp;
         end
         
     end
     
      h1=plot([ObjOnMarker(1), ObjOnMarker(1)],[0, maxval],'--k','LineWidth',1.5);
      h2=plot([ObjOnMarker(2), ObjOnMarker(2)],[0, maxval],'--k','LineWidth',1.5);
  
     set( get( get( h1, 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );
     set( get( get( h2, 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );
     
     currLegend=RepStore(i,:);
     numin=cellfun(@numel,currLegend);
     currLegend=currLegend(numin>0);
if LegendOn
     legend(currLegend);
end
     title(num2str(Objcurr(m)));
     xlabel('Time from object onset(ms)');
      ylabel('Firing Rate(Hz)');
  
     set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
      xlim([-PlotStartOffset(3),PlotStopOffset(3)+ObjOnMarker(2)]);
     
    
    
    
    
end


%EachSceneOnRepBlock_Mean{i,m,indexflag}=nanmean(EachSceneOnRepBlock,1);
%            EachSceneOnRepBlock_Sem{i,m,indexflag}=nanstd(EachSceneOnRepBlock,[],1)./sqrt(sum(~isnan(EachSceneOnRepBlock),1));



FigureIndex=FigureIndex+1;
FigureStartNum=FigureStartNum+1;  

 end% Loop over the block
 
 
end %End of the each scene object on
   
end


end %End of the loop of spikechannel


end
function Rep=DefineRep(Index);
uindex=unique(Index,'Stable');
RepRaw=NaN*ones(1,length(Index));
for i=1:length(uindex)
    Sel=Index==uindex(i);
    
    RepRaw(Sel)=1:sum(Sel);
    
end

Rep=RepRaw;

end