function FourBlockTypeSceneResponseCompare(Data, StartTrial, EndTrial,ShowFigureFlag,OutputFlag);

ECode;
FileName=Data.FileName;
TaskType=Data.TaskType;
TaskCode=Data.TaskCode;
%TrialNum=EndTrial-StartTrial+1;
Selection=StartTrial:EndTrial;%Select the trials according to the pannel

%Abnormal check %Especially For Robin081120FORAGEXPLO4.rst
CheckCode=SCENE_ON;%This code should appear once in each trial
AbnormalTrials=AbnormalTrialCheck(Data,CheckCode);
Selection=Selection(~ismember(Selection,AbnormalTrials));


CheckCode=SCENE_OFF;%This code should appear once in each trial
AbnormalTrials=AbnormalTrialCheck(Data,CheckCode);
Selection=Selection(~ismember(Selection,AbnormalTrials));

if strcmp(FileName,'Robin092920FORAGEXPLO1.rst') ||strcmp(FileName,'Robin092820FORAGEXPLO1.rst')||strcmp(FileName,'Robin012621FORAGEXPLO1.rst')
    %Remove bad trial for this file

CheckCodePair=[FIXOFF,SHOWFIXCD];%This code should be the same as target on for good trials
AbnormalTrials=AbnormalTrialCheck(Data,CheckCodePair);
Selection=Selection(~ismember(Selection,AbnormalTrials));

if strcmp(FileName,'Robin092820FORAGEXPLO1.rst') ||strcmp(FileName,'Robin012621FORAGEXPLO1.rst')
    CheckCodePair=[TGTCD,TGTOFF];%This code should be the same as target on for good trials
    AbnormalTrials=AbnormalTrialCheck(Data,CheckCodePair);
    Selection=Selection(~ismember(Selection,AbnormalTrials));
    
end


%TGTCD
end

TrialIndexRaw=Data.BlockType(:,1);


%Special request
remove_tmp=[];
SceneRemoval=[];
BlockRemoval=[];

%Remove the intermediate scene to focus on the pure scenes if needed
%SceneRemoval=[1054:1057,1060:1063];
if ~isempty(SceneRemoval)
    SceneName_tmp=Data.SceneIndex(Selection,3);
    remove_tmp=ismember(SceneName_tmp,SceneRemoval);
      
end

%Remove the not focus block type to focus on the random blocks if needed
%BlockRemoval=[0,10];
if ~isempty(BlockRemoval)
    Block_tmp=Data.BlockType(Selection,2);
    if ~isempty(remove_tmp)
        remove_tmp=remove_tmp | ismember(Block_tmp,BlockRemoval);
    else
        remove_tmp=ismember(Block_tmp,BlockRemoval);
    end
    
   
end
if ~isempty(remove_tmp)
    Selection=Selection(~remove_tmp);
end
TrialIndex=TrialIndexRaw(Selection);

TrialNum=length(Selection);
%User defined parameter for the neural/behavioral analysis
TargetOnOffset=50;%in ms, the time selection bin after target onset; 
StartDelay=0;%in ms, the time selection bin after fixation offset;

%For segmentation of  the neural responses 
%SegEvent=[EXTRAFIXSTARTCD,EXTRAFIXOFF,SCENE_ON,SCENE_OFF,SHOWFIXCD,FIXOFF,TGTCD];%Event Flow of interest
%SegEventEnd=[EXTRAFIXOFF,SHOWFIXCD,SHOWFIXCD,NaN,FIXOFF,TGTCD,TGTOFF];%Event end; %Automatic screening later
SegEvent=[SCENE_ON,SCENE_OFF];%Event Flow of interest
SegEventEnd=[SHOWFIXCD,NaN];%Event end; %Automatic screening later
%SegEventEnd=[EXTRAFIXOFF,FIXOFF,FIXOFF,TGTCD,TGTOFF];%Event end
SegEventStr={'SceneOn','SceneOff'};

PlotStartOffset=[500,500];%in ms
PlotStopOffset=[500,1000];%in ms

EventTypeFlag=[1,1];%1:Number same with the trial number, 2:Number same with the target number


RasterOn=0;%Show raster, but it will be very slow to do so

%Parameter to define the sliding window to measure firing rates
StepSize=5;%in ms
BinWidth=20;%in ms

%Parameter for raster plot
Barlength=1;
LineWidth=1;


%ResponseLatencyThreshold=3;

%Some user defined parameter for visualization
%Color Definition for each scenes
%ColorForScene={'r','g','b','m',[0.4,0.4,0],[0.28,0.78,0.89],'k',[0.6,0.1,1],[0.6,0.1,0.4],[0.3,0.5,1],[0.6,0.5,0.1],[0,0.8,0.4]};
Color=[lines(6);prism(6);pink(6)];
%Color(2,:)=copper(1)+0.5;
Color(9,:)=[0.28,0.78,0.89];
Color=[Color;hsv(10)];
Color(18,:)=[0.1,0.457,0.194];
ColorForScene=Color;

ColorForEachMode=[0.83 0.14 0.14
             1.00 0.54 0.00
             0.47 0.25 0.80
             0.25 0.80 0.54];
AssociateBlock=[10,11,0,1];
Blockstr={'FixLoc,FixSceneSeq','FixLoc,RandSceneSeq','RandLoc,FixSceneSeq','RandLoc,RandSceneSeq'};
ColorForOrder={'#FB5441','#4173fb','#FBBB41','#0EC335',[0.28,0.78,0.89],'k',[0.6,0.1,1],[0.6,0.1,0.4],[0.3,0.5,1],[0.6,0.5,0.1],[0,0.8,0.4],'r','g','b','m'};

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





if RecordDate<datestr('10,28,2020') && strcmp(MonkeyName,'Robin')
%Get the information about the extra FP
%
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
end %End of the if data & monkey
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
 
 ChoiceOutcome=ChoiceOutcome(AdvanceInTargetTrial,:);
 
 TrialIndex=TrialIndex(AdvanceInTargetTrial);
 
% GOODCHOICECD    
 
 TargetNum=TargetNum(AdvanceInTargetTrial);
 BlockType=BlockType(AdvanceInTargetTrial);
 
 [UniqueBlockType,orderindex,reprojectorderindex]=unique(BlockType,'row');
 UniqueBlockType=sort(UniqueBlockType(~isnan(UniqueBlockType)),'ascend');%Flexible first, then stable;
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

%First calculate average firing rate for the neuron in this task
ISI=nanmean(diff(SpikeTime,1,2),'all');
AverageFR=1/ISI*1000;%Average FR in the whole tak\s;l]]

%Segregate spike data from Fixation onset to target onset
SegmentData=EventSegTool(SegEvent,SegEventEnd,PlotStartOffset,PlotStopOffset,EventChannel,EventTimeChannel,SpikeTime,StepSize,BinWidth);

%Retrive data from the field
%Total average for each event
MeanFR=SegmentData.MeanFR; %[Mean,sem]%Averall mean and sem
MeanFRRaw=SegmentData.MeanFRRaw;%=AverageFiringRate_Event;

%PSTH(Spike Density Function) for each event
MeanPSTH=SegmentData.MeanPSTH;%=[PSTH_Event_Time',PSTH_Event'];
MeanPSTH_Time=MeanPSTH(:,1);
MeanPSTH_FR=MeanPSTH(:,2);

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

%TargetOrderLine
TargetOrder=1:max(TargetNum);
TotalNum=length(TargetNum);
TargetOrder_Mat=repmat(TargetOrder,TotalNum,1);
TargetOrder_Line=reshape(TargetOrder_Mat,numel(TargetOrder_Mat),1);
TargetOrder_Line=TargetOrder_Line(NoNNaNSelector);

%Rearrange searching time into lines
SearchingTime_Line=reshape(SearchingTime,numel(SearchingTime),1);
SearchingTime_Line=SearchingTime_Line(NoNNaNSelector);



%Separate with different fractal order;
ObjOrderSep=[1,2];

% Loop over the blocks
for i=1:NumberOfBlocks
    
  
    %Organize the raster, spikecount PSTH, spike density function PSTH
    for j=1:length(AlignedRaster)
        if EventTypeFlag(j)==1 %Trial event
            currBlocks_ori=BlockType==UniqueBlockType(i);
            
        elseif EventTypeFlag(j)==2 %Target event
            currBlocks_ori=BlockType_Line==UniqueBlockType(i);
        end
        
       currBlocks=currBlocks_ori;
        
        %Raster
        Raster_currEvent_Block{i,j}=AlignedRaster{j}(currBlocks,:);
        
        
       
        %Spike count PSTH
        PSTH_Event_SpikeCount_currEvent{i,j}=PSTH_Event_SpikeCount{j}(currBlocks,:);
        
        %Spike density function PSTH
        MeanPSTH_FR_currEvent{i,j}=MeanPSTH_FR{j}(currBlocks,:);
        
        %Timing 
        PSTH_Time_currEvent{i,j}=nanmean(MeanPSTH_Time{j}(currBlocks,:),1);
        RelativeTimeMarker_currEvent{i,j}=RelativeTimeMarker{j}(currBlocks,:);
        
        %Total Average
        MeanFR_currEvent(i,j)=nanmean(MeanFRRaw{j}(currBlocks));
        RawFR_currEvent{i,j}=MeanFRRaw{j}(currBlocks);
        
        %{
        StartLoc=1;
        Barlength=1;
        LineWidth=1;
        SpikeRaterPlot(Raster_currEvent_Block,StartLoc,Barlength,[1,0,0],LineWidth);
        %}
    end
    %SegEventStr={'ExtraFPOn','ExtraFPOff','SceneOn','SceneOff','FixOff','CueOn','TargetOnset','TargetOffset','Reward','PostTrialSpon','BeforeTrialSpon'};
     
    
 
  
end

%[p_Stable_PSTH,h_Stable_PSTH]=ttest2_group(ResponseOn_HighRT_StableScene,ResponseOn_LowRT_StableScene,'both');
%[p_Flexible_PSTH,h_Flexible_PSTH]=ttest2_group(ResponseOn_HighRT_FlexibleScene,ResponseOn_LowRT_FlexibleScene,'both');



%Loop overt the scenes and blocks
%SceneName
%SegEventStr={'ExtraFPOn','ExtraFPOff','SceneOn','SceneOff','FixOff','CueOn','TargetOnset','TargetOffset','Reward','Spon'};
UniqueScenes=unique(SceneName);
SceneOnEventFlag=ismember(SegEventStr,'SceneOn');
SceneOffEventFlag=ismember(SegEventStr,'SceneOff');



ResponseOn=MeanPSTH_FR{SceneOnEventFlag};
ResponseOff=MeanPSTH_FR{SceneOffEventFlag};
for j=1:NumberOfBlocks
    for i=1:length(UniqueScenes)
    Sel=(SceneName==UniqueScenes(i))& BlockType==UniqueBlockType(j);
%    Resp_tmp=nanmean(ResponseOn(Sel,IntervalSel),2);
%    Scene_Resp_Mean_Whole(i)=nanmean(Resp_tmp);
%    Scene_Resp_Sem_Whole(i)=nanstd(Resp_tmp)/sqrt(sum(~isnan(Resp_tmp)));
    RespOn_tmp_dynamics=ResponseOn(Sel,:);
    SceneOn_Resp_Mean_Each{i,j}=nanmean( RespOn_tmp_dynamics,1);
    SceneOn_Resp_Sem_Each{i,j}=nanstd( RespOn_tmp_dynamics,[],1)./sqrt(sum(~isnan( RespOn_tmp_dynamics),1));
    
    RespOff_tmp_dynamics=ResponseOff(Sel,:);
    SceneOff_Resp_Mean_Each{i,j}=nanmean( RespOff_tmp_dynamics,1);
    SceneOff_Resp_Sem_Each{i,j}=nanstd( RespOff_tmp_dynamics,[],1)./sqrt(sum(~isnan( RespOff_tmp_dynamics),1));
    
    %Num correction (set number of responses too few to be NaN)
      NumOfValidOn=sum(~isnan(RespOn_tmp_dynamics),1);
      RemovalFlagOn=NumOfValidOn<3;
      SceneOn_Resp_Mean_Each{i,j}(RemovalFlagOn)=NaN;
      SceneOn_Resp_Sem_Each{i,j}(RemovalFlagOn)=NaN;
      
      NumOfValidOff=sum(~isnan(RespOff_tmp_dynamics),1);
      RemovalFlagOff=NumOfValidOff<3;
      
      SceneOff_Resp_Mean_Each{i,j}(RemovalFlagOff)=NaN;
      SceneOff_Resp_Sem_Each{i,j}(RemovalFlagOff)=NaN;
      
      %Total mean for each sceneon/off responses;
      SceneOn_Resp_tmp=MeanFRRaw{SceneOnEventFlag}(Sel);
      SceneOff_Resp_tmp=MeanFRRaw{SceneOffEventFlag}(Sel);
      
      SceneOn_Resp_Mean_Whole(i,j)=nanmean(SceneOn_Resp_tmp);
      SceneOn_Resp_Sem_Whole(i,j)=nanstd(SceneOn_Resp_tmp)/sqrt(sum(~isnan(SceneOn_Resp_tmp)));
      
      SceneOff_Resp_Mean_Whole(i,j)=nanmean(SceneOff_Resp_tmp);
      SceneOff_Resp_Sem_Whole(i,j)=nanstd(SceneOff_Resp_tmp)/sqrt(sum(~isnan(SceneOff_Resp_tmp)));
      
     
      
  
    
    end
end
Time_PSTH_SceneOn=nanmean(MeanPSTH_Time{SceneOnEventFlag},1);
Time_PSTH_SceneOff=nanmean(MeanPSTH_Time{SceneOffEventFlag},1);
SceneOnMarker=nanmean(RelativeTimeMarker{SceneOnEventFlag},1);
SceneOffMarker=nanmean(RelativeTimeMarker{SceneOffEventFlag},1);
PlotSceneOnMarkerStart=PlotStartOffset(SceneOnEventFlag);
PlotSceneOnMarkerStop=PlotStopOffset(SceneOnEventFlag);

PlotSceneOffMarkerStart=PlotStartOffset(SceneOffEventFlag);
PlotSceneOffMarkerStop=PlotStopOffset(SceneOffEventFlag);

%%SSI for different scene groups 
SceneGroup={[1000:1007],[1010:1017],[1054:1057,1060:1063]};
SceneGroupStr={'FlexibleScenes','StableScenes','DuelScenes'};
SSI_Group={};
SSI_BlockIndex=[];
SceneTypeIndex=NaN*ones(3,1);
p_SSI_Group={};

SSI_GroupOff={};
p_SSI_GroupOff={};

PSTH_Mean4Block={};
PSTH_Sem4Block={};

PSTHOff_Mean4Block={};
PSTHOff_Sem4Block={};


for i=1:3
    currSceneGroup=SceneGroup{i};
    SelSceneGroup=ismember(SceneName,currSceneGroup);
    SceneGroupBlock=BlockType(SelSceneGroup);
    UniqueBlockTypeGroup=unique(SceneGroupBlock);
    if isempty(UniqueBlockTypeGroup)
        SSI_Group{i,1}=[];
        SSI_BlockIndex{i,1}=[];
        SceneTypeIndex(i)=0;
        SSI_GroupOff{i,1}=[];
        
        PSTH_Mean4Block{i,1}=[];
       PSTH_Sem4Block{i,1}=[];
       PSTHOff_Mean4Block{i,1}=[];
       PSTHOff_Sem4Block{i,1}=[];

   
    else
        SceneTypeIndex(i)=1;
        for j=1:length(UniqueBlockTypeGroup)
            currBlockSel=BlockType==UniqueBlockTypeGroup(j);
            %Scene On
            PSTH_tmp=ResponseOn(SelSceneGroup&currBlockSel,:);
            currScene=SceneName(SelSceneGroup&currBlockSel);
            
            %SSI
            [pSceneOnGroup,rsquare_SceneOnGroup,~]=anova1_group_stats(PSTH_tmp,currScene,'off');
            %Fewer number correction
            NumOfData=sum(~isnan(PSTH_tmp),1);
            NumOfSceneInGroup=length(unique(currScene));
            Criterium=3;
            Remove=NumOfData<NumOfSceneInGroup*Criterium;
            rsquare_SceneOnGroup(Remove)=NaN;
            pSceneOnGroup(Remove)=NaN;
            SSI_BlockIndex{i,j}=UniqueBlockTypeGroup(j);
            SSI_Group{i,j}=rsquare_SceneOnGroup;
            p_SSI_Group{i,j}=pSceneOnGroup;
            
            %Mean
            PSTH_Mean4Block{i,j}=nanmean(PSTH_tmp,1);
            PSTH_Sem4Block{i,j}=nanstd(PSTH_tmp,[],1)./sqrt(sum(~isnan(PSTH_tmp)));
            
            
            %SceneOff
            PSTH_tmpoff=ResponseOff(SelSceneGroup&currBlockSel,:);
            [pSceneOffGroup,rsquare_SceneOffGroup,~]=anova1_group_stats(PSTH_tmpoff,currScene,'off');
            %Fewer number correction
            NumOfDataoff=sum(~isnan(PSTH_tmpoff),1);
            NumOfSceneInGroup=length(unique(currScene));
            Criterium=3;
            Removeoff=NumOfDataoff<NumOfSceneInGroup*Criterium;
            rsquare_SceneOffGroup(Removeoff)=NaN;
            pSceneOffGroup(Removeoff)=NaN;
            
            
            SSI_GroupOff{i,j}=rsquare_SceneOffGroup;
            p_SSI_GroupOff{i,j}=pSceneOffGroup;
            
            %Mean
            PSTHOff_Mean4Block{i,j}=nanmean(PSTH_tmpoff,1);
            PSTHOff_Sem4Block{i,j}=nanstd(PSTH_tmpoff,[],1)./sqrt(sum(~isnan(PSTH_tmpoff)));
            
            
        end
        
    end
   % PSTH_tmp=ResponseOn(SelSceneGroup,:);
  
    
end




%%
FigureStartNum=110;
FigureIndex=1;
if ShowFigureFlag
   
  

  
%%%%%%%%%%%%%%%%%%Next figure: Scene Onset Offset PSTH for each scene

 figtitlestr{FigureIndex}='SceneOnEachBlock';
 fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[50,100, 1600,800],'Name',figtitlestr{FigureIndex});   
 
 for i=1:length(UniqueScenes)
     subplot(ceil(sqrt(length(UniqueScenes))),ceil(sqrt(length(UniqueScenes))),i);
     for j=1:NumberOfBlocks
          currBlock=UniqueBlockType(j);
           BlockOrder(j)=find(AssociateBlock==currBlock);
     shadedErrorBar(Time_PSTH_SceneOn,SceneOn_Resp_Mean_Each{i,j},SceneOn_Resp_Sem_Each{i,j},'lineprops',{ColorForEachMode(AssociateBlock==currBlock,:)},'transparent',1,'patchSaturation',0.5);
     hold on
     
     end
     
      h1=plot([SceneOnMarker(1), SceneOnMarker(1)],[0,max(max(cellfun(@(x,y) max(max(x+y)),SceneOn_Resp_Mean_Each,SceneOn_Resp_Sem_Each)))],'--k','LineWidth',1.5);
  h2=plot([SceneOnMarker(2), SceneOnMarker(2)],[0,max(max(cellfun(@(x,y) max(max(x+y)),SceneOn_Resp_Mean_Each,SceneOn_Resp_Sem_Each)))],'--k','LineWidth',1.5);
  
  set( get( get( h1, 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );
 set( get( get( h2, 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );

     set(findall(gca, 'Type', 'Line'),'LineWidth',3);
  xlabel('Time from scene onset(ms)');
  ylabel('Firing Rate(Hz)');
    title(strcat('Scene',string(UniqueScenes(i))));

 % title('SceneOnset');
  set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
   xlim([-PlotSceneOnMarkerStart,PlotSceneOnMarkerStop+SceneOnMarker(2)]);
   
   if i==length(UniqueScenes)
       ld=legend(Blockstr{BlockOrder});
   end     
   
 end
 

%ld=legend([string(UniqueScene);'SSI']);
 %Add data information
  %  data_msg=sprintf(' %s, SpikeChannel: %d, TrialNum: %d ',FileName,SpikeChannel,TrialNum);
    data_msg=strcat(FileName,' SpikeChannel:',num2str(SpikeChannel),';ChannelID:',num2str(ChannelID),';ClusterID:',num2str(ClusterID),';TrialNum: ',num2str(TrialNum));

    annotation('textbox', [0.17, 0.95, 0.7, 0.04], 'string', data_msg,'FontSize',20,'FontWeight','Bold','LineStyle','none')
 
    
    
    
FigureIndex=FigureIndex+1;
FigureStartNum=FigureStartNum+1;  


%%%%%%%%%%%%%%%%%%Next figure: Scene Offset PSTH for each scene

 figtitlestr{FigureIndex}='SceneOffEachBlock';
 fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[90,100, 1600,800],'Name',figtitlestr{FigureIndex});   
 
 for i=1:length(UniqueScenes)
     subplot(ceil(sqrt(length(UniqueScenes))),ceil(sqrt(length(UniqueScenes))),i);
     for j=1:NumberOfBlocks
          currBlock=UniqueBlockType(j);
           BlockOrder(j)=find(AssociateBlock==currBlock);
     shadedErrorBar(Time_PSTH_SceneOff,SceneOff_Resp_Mean_Each{i,j},SceneOff_Resp_Sem_Each{i,j},'lineprops',{ColorForEachMode(AssociateBlock==currBlock,:)},'transparent',1,'patchSaturation',0.5);
     hold on
     
     end
     
      h1=plot([SceneOffMarker(1), SceneOffMarker(1)],[0,max(max(cellfun(@(x,y) max(max(x+y)),SceneOff_Resp_Mean_Each,SceneOff_Resp_Sem_Each)))],'--k','LineWidth',1.5);
  h2=plot([SceneOffMarker(2), SceneOffMarker(2)],[0,max(max(cellfun(@(x,y) max(max(x+y)),SceneOff_Resp_Mean_Each,SceneOff_Resp_Sem_Each)))],'--k','LineWidth',1.5);
  
  set( get( get( h1, 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );
 set( get( get( h2, 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );

     set(findall(gca, 'Type', 'Line'),'LineWidth',3);
  xlabel('Time from scene offset(ms)');
  ylabel('Firing Rate(Hz)');
    title(strcat('Scene',string(UniqueScenes(i))));

 % title('SceneOnset');
  set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
   xlim([-PlotSceneOffMarkerStart,PlotSceneOffMarkerStop*1.25]);
   
   if i==length(UniqueScenes)
       ld=legend(Blockstr{BlockOrder});
   end     
   
 end
 

%ld=legend([string(UniqueScene);'SSI']);
 %Add data information
  %  data_msg=sprintf(' %s, SpikeChannel: %d, TrialNum: %d ',FileName,SpikeChannel,TrialNum);
    data_msg=strcat(FileName,' SpikeChannel:',num2str(SpikeChannel),';ChannelID:',num2str(ChannelID),';ClusterID:',num2str(ClusterID),';TrialNum: ',num2str(TrialNum));

    annotation('textbox', [0.17, 0.95, 0.7, 0.04], 'string', data_msg,'FontSize',20,'FontWeight','Bold','LineStyle','none')
 
    
    
    
FigureIndex=FigureIndex+1;
FigureStartNum=FigureStartNum+1;  



 %%%%%%%%%%%%%%%%%%Next figure:SSI for each block

figtitlestr{FigureIndex}='SSI_SceneOnEachBlock';
fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[110,100, 1200,800],'Name',figtitlestr{FigureIndex});   

%Organize according to the SceneTypeIndex
maxVal=0;
if sum(SceneTypeIndex)==1
    NumOfPlots=1;
else
    NumOfPlots=sum(SceneTypeIndex);
end
plotindex=1;
for i=1:3
    tmp=cell2mat(SSI_BlockIndex(i,:));
    if ~isempty(tmp)
        if NumOfPlots>1
            subplot(NumOfPlots,1,plotindex);
            plotindex=plotindex+1;
        end
        for j=1:length(tmp)
            
            currBlock=tmp(j);
            BlockOrder(j)=find(AssociateBlock==currBlock);
            plot(Time_PSTH_SceneOn,SSI_Group{i,j},'-r','Color',ColorForEachMode(AssociateBlock==currBlock,:))
            hold on
            
            maxVal_tmp=max(SSI_Group{i,j});
            if maxVal_tmp>maxVal
                maxVal=maxVal_tmp;
            end
            
            
        end
         xlabel('Time from scene onset(ms)');
         ylabel('R Square');
         title(SceneGroupStr{i});
          ld=legend(Blockstr{BlockOrder});
        h1=plot([SceneOnMarker(1), SceneOnMarker(1)],[0,maxVal+0.1],'--k','LineWidth',1.5);
  h2=plot([SceneOnMarker(2), SceneOnMarker(2)],[0,maxVal+0.1],'--k','LineWidth',1.5);
  
  set( get( get( h1, 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );
 set( get( get( h2, 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );

     set(findall(gca, 'Type', 'Line'),'LineWidth',3);
 
   

 % title('SceneOnset');
  set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
   xlim([-PlotSceneOnMarkerStart,PlotSceneOnMarkerStop+SceneOnMarker(2)]);
   
   
        
    end
    
    
    
    
    
    
end
 

   data_msg=strcat(FileName,' SpikeChannel:',num2str(SpikeChannel),';ChannelID:',num2str(ChannelID),';ClusterID:',num2str(ClusterID),';TrialNum: ',num2str(TrialNum));

    annotation('textbox', [0.17, 0.95, 0.7, 0.04], 'string', data_msg,'FontSize',20,'FontWeight','Bold','LineStyle','none')
 
      
  
   
 
FigureIndex=FigureIndex+1;
FigureStartNum=FigureStartNum+1;  


%%%%%%%%%%%%%%%%%%Next figure:SSI of scene off for each block

figtitlestr{FigureIndex}='SSI_SceneOffEachBlock';
fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[700,100, 1200,800],'Name',figtitlestr{FigureIndex});   

%Organize according to the SceneTypeIndex
maxVal=0;
if sum(SceneTypeIndex)==1
    NumOfPlots=1;
else
    NumOfPlots=sum(SceneTypeIndex);
end
plotindex=1;
for i=1:length(SceneGroup)
    tmp=cell2mat(SSI_BlockIndex(i,:));
    if ~isempty(tmp)
        if NumOfPlots>1
            subplot(NumOfPlots,1,plotindex);
            plotindex=plotindex+1;
        end
        for j=1:length(tmp)
            
            currBlock=tmp(j);
            BlockOrder(j)=find(AssociateBlock==currBlock);
            plot(Time_PSTH_SceneOff,SSI_GroupOff{i,j},'-r','Color',ColorForEachMode(AssociateBlock==currBlock,:))
            hold on
            
            maxVal_tmp=max(SSI_GroupOff{i,j});
            if maxVal_tmp>maxVal
                maxVal=maxVal_tmp;
            end
            
            
        end
         xlabel('Time from scene offset(ms)');
         ylabel('R Square');
         title(SceneGroupStr{i});
          ld=legend(Blockstr{BlockOrder});
        
          h1=plot([SceneOffMarker(1), SceneOffMarker(1)],[0,maxVal+0.1],'--k','LineWidth',1.5);
  h2=plot([SceneOffMarker(2), SceneOffMarker(2)],[0,maxVal+0.1],'--k','LineWidth',1.5);
  
  set( get( get( h1, 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );
 set( get( get( h2, 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );

     set(findall(gca, 'Type', 'Line'),'LineWidth',3);
 
   

 % title('SceneOnset');
  set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
   xlim([-PlotSceneOffMarkerStart,PlotSceneOffMarkerStop*1.25]);
   
        
    end
    
    
    
    
    
    
end
 

 
   
    data_msg=strcat(FileName,' SpikeChannel:',num2str(SpikeChannel),';ChannelID:',num2str(ChannelID),';ClusterID:',num2str(ClusterID),';TrialNum: ',num2str(TrialNum));

    annotation('textbox', [0.17, 0.95, 0.7, 0.04], 'string', data_msg,'FontSize',20,'FontWeight','Bold','LineStyle','none')
     
  
   
 
FigureIndex=FigureIndex+1;
FigureStartNum=FigureStartNum+1; 


%%%%%%%%%%%%%%%%%%Next figure:Mean onset responses of scene for each block

figtitlestr{FigureIndex}='Mean_SceneOnEachBlock';
fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[100,100, 1200,800],'Name',figtitlestr{FigureIndex});   
%Organize according to the SceneTypeIndex
maxVal=0;
if sum(SceneTypeIndex)==1
    NumOfPlots=1;
else
    NumOfPlots=sum(SceneTypeIndex);
end
plotindex=1;
for i=1:length(SceneGroup)
    tmp=cell2mat(SSI_BlockIndex(i,:));
    if ~isempty(tmp)
        if NumOfPlots>1
            subplot(NumOfPlots,1,plotindex);
            plotindex=plotindex+1;
        end
        for j=1:length(tmp)
            
            currBlock=tmp(j);
            BlockOrder(j)=find(AssociateBlock==currBlock);
            
            
            shadedErrorBar(Time_PSTH_SceneOn,PSTH_Mean4Block{i,j}, PSTH_Sem4Block{i,j},'lineprops',{ColorForEachMode(AssociateBlock==currBlock,:)},'transparent',1,'patchSaturation',0.5);%PSTH for Flexible
            set(findall(gca, 'Type', 'Line'),'LineWidth',3);
            
            
            %PSTH_Mean4Block{i,j}
            hold on
            
            maxVal_tmp=max(PSTH_Mean4Block{i,j});
            if maxVal_tmp>maxVal
                maxVal=maxVal_tmp;
            end
            
            
        end
         xlabel('Time from scene onset(ms)');
         ylabel('Mean reponseses(Hz)');
         title(SceneGroupStr{i});
          ld=legend(Blockstr{BlockOrder});
        
          h1=plot([SceneOnMarker(1), SceneOnMarker(1)],[0,maxVal*1.1],'--k','LineWidth',1.5);
  h2=plot([SceneOnMarker(2), SceneOnMarker(2)],[0,maxVal*1.1],'--k','LineWidth',1.5);
  
  set( get( get( h1, 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );
 set( get( get( h2, 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );

     set(findall(gca, 'Type', 'Line'),'LineWidth',3);
 
   

 % title('SceneOnset');
  set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
    xlim([-PlotSceneOnMarkerStart,PlotSceneOnMarkerStop+SceneOnMarker(2)]);
   
        
    end
    
    
    
    
    
    
end
 

 
   
    data_msg=strcat(FileName,' SpikeChannel:',num2str(SpikeChannel),';ChannelID:',num2str(ChannelID),';ClusterID:',num2str(ClusterID),';TrialNum: ',num2str(TrialNum));

    annotation('textbox', [0.17, 0.95, 0.7, 0.04], 'string', data_msg,'FontSize',20,'FontWeight','Bold','LineStyle','none')
     
  
   
 
FigureIndex=FigureIndex+1;
FigureStartNum=FigureStartNum+1; 


%%%%%%%%%%%%%%%%%%Next figure:Mean offset responses of scene for each block

figtitlestr{FigureIndex}='Mean_SceneOffEachBlock';
fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[100,100, 1200,800],'Name',figtitlestr{FigureIndex});   
%Organize according to the SceneTypeIndex
maxVal=0;
if sum(SceneTypeIndex)==1
    NumOfPlots=1;
else
    NumOfPlots=sum(SceneTypeIndex);
end
plotindex=1;
for i=1:length(SceneGroup)
    tmp=cell2mat(SSI_BlockIndex(i,:));
    if ~isempty(tmp)
        if NumOfPlots>1
            subplot(NumOfPlots,1,plotindex);
            plotindex=plotindex+1;
        end
        for j=1:length(tmp)
            
            currBlock=tmp(j);
            BlockOrder(j)=find(AssociateBlock==currBlock);
            
            
            shadedErrorBar(Time_PSTH_SceneOff,PSTHOff_Mean4Block{i,j}, PSTHOff_Sem4Block{i,j},'lineprops',{ColorForEachMode(AssociateBlock==currBlock,:)},'transparent',1,'patchSaturation',0.5);%PSTH for Flexible
            set(findall(gca, 'Type', 'Line'),'LineWidth',3);
            
            
            %PSTH_Mean4Block{i,j}
            hold on
            
            maxVal_tmp=max(PSTHOff_Mean4Block{i,j});
            if maxVal_tmp>maxVal
                maxVal=maxVal_tmp;
            end
            
            
        end
         xlabel('Time from scene offset(ms)');
         ylabel('Mean reponseses(Hz)');
         title(SceneGroupStr{i});
          ld=legend(Blockstr{BlockOrder});
        
          h1=plot([SceneOffMarker(1), SceneOffMarker(1)],[0,maxVal*1.1],'--k','LineWidth',1.5);
  h2=plot([SceneOffMarker(2), SceneOffMarker(2)],[0,maxVal*1.1],'--k','LineWidth',1.5);
  
  set( get( get( h1, 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );
 set( get( get( h2, 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );

     set(findall(gca, 'Type', 'Line'),'LineWidth',3);
 
   

 % title('SceneOnset');
  set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
    xlim([-PlotSceneOffMarkerStart,PlotSceneOffMarkerStop*1.25]);
   
        
    end
    
    
    
    
    
    
end
 

 
   
    data_msg=strcat(FileName,' SpikeChannel:',num2str(SpikeChannel),';ChannelID:',num2str(ChannelID),';ClusterID:',num2str(ClusterID),';TrialNum: ',num2str(TrialNum));

    annotation('textbox', [0.17, 0.95, 0.7, 0.04], 'string', data_msg,'FontSize',20,'FontWeight','Bold','LineStyle','none')
     
  
   
 
FigureIndex=FigureIndex+1;
FigureStartNum=FigureStartNum+1; 


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
OutputFigureName=sprintf('%s_%s_N%s_C%s_%CID_%NID',MonkeyName,string(RecordDateOriginal),string(NeuronNum),string(ChannelNum),string(ChannelID),string(ClusterID));



ExistFlag=0;
%Load the old file if exist
if exist(OutputFileName)

    load(OutputFileName);
    %Load OutputData into the memory
    if isfield(OutputFileName,'EventSegRT')
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

OutputData_New.FourBlockSceneCompare.TrialType='GoodTrials';

%%%%%%%%%Save overall scene response characteristics%%%%%%%%%%%%%%%%%%%%%%%
%Organize the data
NameStr={'Time_PSTH_SceneOn','SSI_Group','SSI_BlockIndex','SceneTypeIndex','SSI_GroupOff',' p_SSI_GroupOff',' p_SSI_GroupOff',... 
    'Time_PSTH_SceneOff'
 
};



DataLib={
   Time_PSTH_SceneOn,SSI_Group,SSI_BlockIndex,SceneTypeIndex,SSI_GroupOff, p_SSI_GroupOff, p_SSI_Group,...
   Time_PSTH_SceneOff

};

 
        
%{
%Use anova to statistically test scene selectivity
SceneOnResponseAll=MeanFRRaw{SceneOnEventFlag};
[p_SceneOnAnova,stats]=anova1(SceneOnResponseAll,SceneName,'off');
%Separate for each block
Stable_Raw_Sel=ismember(SceneName,UniqueScene_Stable);
SceneOnResponseStable_Raw=MeanFRRaw{SceneOnEventFlag}(Stable_Raw_Sel);
SceneName_Stable_Raw=SceneName(Stable_Raw_Sel);
[p_SceneOnAnova_Stable,stats]=anova1(SceneOnResponseStable_Raw,SceneName_Stable_Raw,'off');


Flexible_Raw_Sel=ismember(SceneName,UniqueScene_Flexible);
SceneOnResponseFlexible_Raw=MeanFRRaw{SceneOnEventFlag}(Flexible_Raw_Sel);
SceneName_Flexible_Raw=SceneName(Flexible_Raw_Sel);
[p_SceneOnAnova_Flexible,stats]=anova1(SceneOnResponseFlexible_Raw,SceneName_Flexible_Raw,'off');
%}


%SceneOnAverage=MeanFR_currEvent(:,SceneOnEventFlag);
%[h,p_SceneDiffBetweenBlocks]=ttest2(RawSceneOnAverage{1},RawSceneOnAverage{2});

%'MeanPSTH_Event','EventName','EventTimeMarker'
%MeanPSTH_FR_currEvent,SegEventStr,RelativeTimeMarker_currEvent


DataStamp=containers.Map(NameStr,DataLib);

OutputData_New.FourBlockSceneCompare.Task=TaskType;
OutputData_New.FourBlockSceneCompare.TaskCode=TaskCode;
OutputData_New.FourBlockSceneCompare.DataStamp=DataStamp;
%Store figures in common folders
ProtocolMarker='FourBlockSceneCompare';
FigurePath=strcat(BasicDirectory,'Figures',BasicDirectory(end),ProtocolMarker,BasicDirectory(end));%,OutputFolerName);
if ~exist(FigurePath)
    mkdir(FigurePath)
end


if ShowFigureFlag & OutputFlag
    for i=1:length(fig)
        cd(FigurePath);
        OutputFigureNameFig=[];
        titlestr=figtitlestr{i};
        FolderName=titlestr;
         if ~exist(FolderName,'dir')
          mkdir(FolderName);
         end
          cd(FolderName);
    
        OutputFigureNameFig=strcat(OutputFigureName,'.jpg');
        saveas(fig{i},OutputFigureNameFig);
        
        
    end
        
    
    disp('Figures have been exported to the common folder');
    
end
if ExistFlag
%Compare the old one with the new one
Task_Old=OutputData.FourBlockSceneCompare.TaskCode;

     if Task_Old~=TaskCode
        %Not the same task, add the new dataset into the old one
        OutputData.FourBlockSceneCompare(numel(OutputData.FourBlockSceneCompare)+1)=OutputData_New.FourBlockSceneCompare;
     else
         %If the same task,the same protocol,replace the old one with the new one
         
         OutputData.FourBlockSceneCompare=OutputData_New.FourBlockSceneCompare;
         
    
    
     end
else
    %Output the current file
   OutputData.FourBlockSceneCompare=OutputData_New.FourBlockSceneCompare; 
end
if OutputFlag
    cd(OutPath);
    save(OutputFileName,'OutputData');
    disp('Data Exported');
else
    disp('No data export');
end  
end %End of the loop of the spike channel


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