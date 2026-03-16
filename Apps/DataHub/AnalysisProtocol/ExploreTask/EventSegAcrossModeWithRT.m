function EventSegAcrossModeWithRT(Data, StartTrial, EndTrial,ShowFigureFlag,OutputFlag);
%summary of all marked event of particular neuron,
%And try to related them with RT;
%by Xuefei Yu 01102020
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

if strcmp(FileName,'Robin092920FORAGEXPLO1.rst') ||strcmp(FileName,'Robin092820FORAGEXPLO1.rst')||strcmp(FileName,'Robin012621FORAGEXPLO1.rst')||strcmp(FileName,'Robin040221FORAGEXPLOITL1.rst')
    %Remove bad trial for this file

CheckCodePair=[FIXOFF,SHOWFIXCD];%This code should be the same as target on for good trials
AbnormalTrials=AbnormalTrialCheck(Data,CheckCodePair);
Selection=Selection(~ismember(Selection,AbnormalTrials));

if strcmp(FileName,'Robin092820FORAGEXPLO1.rst') ||strcmp(FileName,'Robin012621FORAGEXPLO1.rst')||strcmp(FileName,'Robin040221FORAGEXPLOITL1.rst')
    CheckCodePair=[TGTCD,TGTOFF];%This code should be the same as target on for good trials
    AbnormalTrials=AbnormalTrialCheck(Data,CheckCodePair);
    Selection=Selection(~ismember(Selection,AbnormalTrials));
    
end


%TGTCD
end

TrialIndexRaw=Data.BlockType(:,1);

SceneOnAveTime=600;%Average Interval For SceneOnResponses

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
SegEvent=[EXTRAFIXSTARTCD,EXTRAFIXOFF,SCENE_ON,SCENE_OFF,FIXOFF,CUEON,TGTCD, TGTOFF,RWDCD,   SCENE_OFF,ENABLECD];%Event Flow of interest
SegEventEnd=[EXTRAFIXOFF, PREFIX	,SHOWFIXCD,NaN,      CUEON, TGTCD,TGTOFF,NaN,   RWDOFFCD,NaN,      EXTRAFIXSTARTCD];%Event end; %Automatic screening later
%SegEventEnd=[EXTRAFIXOFF,FIXOFF,FIXOFF,TGTCD,TGTOFF];%Event end
SegEventStr={'ExtraFPOn','ExtraFPOff','SceneOn','SceneOff','FixOff','CueOn','TargetOnset','TargetOffset','Reward','PostTrialSpon','BeforeTrialSpon'};

PlotStartOffset=[200,200,600,500,1000,500,100,100,100,-1500,-100];%in ms
PlotStopOffset=[150,200,500,1500,100,50,100,500,100,2000,-100];%in ms

EventTypeFlag=[1,1,1,1,2,2,2,2,1,1,1];%1:Number same with the trial number, 2:Number same with the target number


RasterOn=0;%Show raster, but it will be very slow to do so

%Parameter to define the sliding window to measure firing rates
StepSize=10;%in ms
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

CueOnCheck=FindOutTime(EventChannel,EventBin,SegEvent(SegEvent==CUEON));
if sum(isnan(CueOnCheck))==numel(CueOnCheck)
    SegEvent(SegEvent==CUEON)=[TGTCD];
    SegEventStr{ismember(SegEventStr,'CueOn')}='TgtOn';
    SegEventEnd(SegEventEnd==CUEON)=[TGTCD];
    disp('Change to TgtOn for lacking of cue');
    
end

PreFixCheck=FindOutTime(EventChannel,EventBin,SegEventEnd(SegEventEnd==PREFIX));
if sum(isnan(PreFixCheck))==numel(PreFixCheck)
    SegEventEnd(SegEventEnd==PREFIX)=[SCENE_ON];
    disp('Change to PreFix to SceneOnset');
    
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


%Remove extraFP marker if it is not existed
CheckExtraFP=FindOutTime(EventChannel,EventBin,EXTRAFIXSTARTCD);
if sum(isnan(CheckExtraFP))==length(CheckExtraFP)
    Remove=ismember(SegEvent,[EXTRAFIXSTARTCD,EXTRAFIXOFF]);
   % SetNaNStart=ismember(SegEvent,EXTRAFIXOFF);
    
    
    PlotStartOffset(Remove)=[];
    PlotStopOffset(Remove)=[];
    EventTypeFlag(Remove)=[];
    SegEvent(Remove)=[];
    SegEventEnd(Remove)=[];
    SegEventStr(Remove)=[];
    
   % PlotStartOffset(SetNaNStart)=1000;
    %PlotStopOffset(Remove)=[];
   % SegEvent(SetNaNStart)=NaN;
   
    SetNaNEnd=ismember(SegEventEnd,EXTRAFIXSTARTCD);
    PlotStopOffset(SetNaNEnd)=1000;
    SegEventEnd(SetNaNEnd)=NaN;
    disp('Remove extraFP Marker');
    
   

    %SegEvent=[EXTRAFIXSTARTCD,EXTRAFIXOFF,SCENE_ON,SCENE_OFF,FIXOFF,CUEON,TGTCD, TGTOFF,RWDCD,   SCENE_OFF,ENABLECD];%Event Flow of interest
%SegEventEnd=[EXTRAFIXOFF, PREFIX	,SHOWFIXCD,NaN,      CUEON, TGTCD,TGTOFF,NaN,   RWDOFFCD,NaN,      EXTRAFIXSTARTCD];%Event end; %Automatic screening later
    
end

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
 UniqueBlockType=sort(UniqueBlockType,'ascend');%Flexible first, then stable;
NumberOfBlocks=max(reprojectorderindex);
 
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

%Rearrange scene name into lines
SceneName_Line=repmat(SceneName,1,max(TargetNum));
SceneName_Line=reshape(SceneName_Line,numel(SceneName_Line),1);
SceneName_Line=SceneName_Line(NoNNaNSelector);

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

%Average of the spon
SponSel=ismember(SegEventStr,'PostTrialSpon');
SponResponse_PSTH_tmp=MeanPSTH_FR{SponSel};
SponTime_PSTH_tmp=nanmean(MeanPSTH_Time{SponSel},1);
TimeInterval=[abs(PlotStartOffset(SponSel)),abs(PlotStopOffset(SponSel))];
SelInterval=(SponTime_PSTH_tmp>=TimeInterval(1)) & (SponTime_PSTH_tmp<=TimeInterval(2));
SponResponse=nanmean(SponResponse_PSTH_tmp(:,SelInterval),2);

%Total Average of the spon
AverageSpon=nanmean(SponResponse);

%Average of the before spon
BeforeSponSel=ismember(SegEventStr,'BeforeTrialSpon');
BeforeSponResponse_PSTH_tmp=MeanPSTH_FR{BeforeSponSel};
BeforeSponTime_PSTH_tmp=nanmean(MeanPSTH_Time{BeforeSponSel},1);
MarkerTimeBeforeSpon=nanmean(RelativeTimeMarker{BeforeSponSel},1);
if ~isnan(MarkerTimeBeforeSpon)
BeforeTimeInterval=[abs(PlotStartOffset(BeforeSponSel)),MarkerTimeBeforeSpon(2)+abs(PlotStopOffset(BeforeSponSel))];
BeforeSelInterval=(BeforeSponTime_PSTH_tmp>=BeforeTimeInterval(1)) & (BeforeSponTime_PSTH_tmp<=BeforeTimeInterval(2));
BeforeSponResponse=nanmean(BeforeSponResponse_PSTH_tmp(:,BeforeSelInterval),2);
else
    BeforeSponResponse=NaN*ones(length(MarkerTimeBeforeSpon),1);
end




%Separate scene types 
SceneWithStableBlock=SceneName(BlockType==11);
SceneWithFlexibleBlock=SceneName(BlockType==1);

UniqueScene_Stable=unique(SceneWithStableBlock);
UniqueScene_Flexible=unique(SceneWithFlexibleBlock);

%Separate with different fractal order;
ObjOrderSep=[1,2];
FixOffSel=ismember(SegEventStr,'FixOff');
PSTH_fixoffAll= MeanPSTH_FR{FixOffSel};

% Loop over the blocks
for i=1:NumberOfBlocks
    
  
    %Organize the raster, spikecount PSTH, spike density function PSTH
    for j=1:length(AlignedRaster)
        if EventTypeFlag(j)==1 %Trial event
            currBlocks_ori=BlockType==UniqueBlockType(i);
            
        elseif EventTypeFlag(j)==2 %Target event
            currBlocks_ori=BlockType_Line==UniqueBlockType(i);
        end
        %Special for reward
        if strcmp( SegEventStr(j),'Reward')
            RewardTrials=ChoiceOutcome==GOODCHOICECD;
            currBlocks=currBlocks_ori & RewardTrials;
        else
            
            currBlocks=currBlocks_ori;
        end
       
       
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
     
    
    
    FixOffSel=ismember(SegEventStr,'FixOff');
    PSTH_fixoff= MeanPSTH_FR_currEvent{i,FixOffSel};
    PSTH_fixoff_time{i}=nanmean(PSTH_Time_currEvent{i,FixOffSel},1);
    Block_Sel=BlockType_Line==UniqueBlockType(i);
    SearchingTime_currBlock=SearchingTime_Line(Block_Sel);
    %Separate into different fractal orders for the fix off response
     SceneCurrBlock=SceneName_Line(Block_Sel);
  for j=1:length(ObjOrderSep)+1
          if j==length(ObjOrderSep)+1
             sel=TargetOrder_Line(Block_Sel)>ObjOrderSep(end);
             order_string{i,j}=sprintf('Fractal Order> %d (N=%d)',ObjOrderSep(end),sum(sel));
             
          else
            sel=TargetOrder_Line(Block_Sel)==ObjOrderSep(j);
            order_string{i,j}=sprintf('Fractal Order %d (N=%d)',ObjOrderSep(j),sum(sel));
           
          end
          SceneCurrOrder=SceneCurrBlock(sel);
         
     PSTH_tmp=PSTH_fixoff(sel,:);
 
     FixOffPSTH_currBlock_OrderSep_mean{i,j}=nanmean(PSTH_tmp,1);
     FixOffPSTH_currBlock_OrderSep_sem{i,j}=nanstd(PSTH_tmp,[],1)./sqrt(sum(~isnan(PSTH_tmp),1));
    
        
     RearchingTime_OrderSep(i,j)=nanmean(SearchingTime_currBlock( sel));
     
     %% SSI for fix off
     
      [pFixOffWithTime,rsquare_psth_fixoff,~]=anova1_group_stats(PSTH_tmp,SceneCurrOrder,'off');
     
      UniqueScenesHere=unique(SceneCurrOrder);
      %
%Few number correction
NumOfValidOn=sum(~isnan(PSTH_tmp),1);
RemovalFlagOn=NumOfValidOn<3*length(UniqueScenesHere);
pFixOffWithTime(RemovalFlagOn)=NaN;
rsquare_psth_fixoff(RemovalFlagOn)=NaN;

%}

h_sig_ssi_fixoffwithtime=floor(pFixOffWithTime<0.05);
h_sig_ssi_fixoffwithtime(h_sig_ssi_fixoffwithtime==0)=NaN;
      
rsquare_psth_fixoff_order{i,j}=rsquare_psth_fixoff;
pFixOffWithTime_order{i,j}=pFixOffWithTime;
h_sig_ssi_fixoffwithtime_order{i,j}=h_sig_ssi_fixoffwithtime;
   
     
  end
  
  
  
  %Next. Sort neural responses according to the Rearching time

%GroupStr={'Highest','Higher Medium','Mid Medium','Lowest'};


%PSTH_currBlock=PSTH_fixoff;%-spon_currBlocks;

%RT_First_Target=SacEnd_First_Target;
SepGroupAmp=4;
SepGroupAmpNum= SepGroupAmp;
%SepGroupAmp={[2300,1800,1100]};
%{
%SepGroupAmp={[500]};
GroupStr={'Highest','Higher Medium','Mid Medium','Lowest'};
%GroupStr={'Long','Medium','Short'};
SepGroupAmpNum=length(SepGroupAmp{1})+1;
%}
GroupStr={'Highest','Higher Medium','Mid Medium','Lowest'};
%Before Exploration
PSTH_group_temp=SortDataAccordingToGroup(PSTH_fixoff,SearchingTime_currBlock,SepGroupAmp,'descend');


FixOffPSTH_RTGroup_mean{i}=PSTH_group_temp.mean(1:SepGroupAmpNum,:);
FixOffPSTH_RTGroup_sem{i}=PSTH_group_temp.sem(1:SepGroupAmpNum,:);

RT_Group_mean{i}=PSTH_group_temp.MarkerMean;
RT_Group_sem{i}=PSTH_group_temp.MarkerSem;

NumInGroupRT=PSTH_group_temp.Num;
RT_SortIndex=PSTH_group_temp.SortIndex;




  
end
%Fix off marker time
FixOffStartOffset=PlotStartOffset(FixOffSel);
FixOffStopOffset=PlotStopOffset(FixOffSel);
%FixOffTimeMarker=nanmean(cell2mat(RelativeTimeMarker_currEvent(:,FixOffSel)),1);
%FixOffStopOffset=FixOffTimeMarker(2)+FixOffStopOffset;


%PlotStartOffset=[200,200,600,500,500,500,100,100,100,-1500,-100];%in ms
%PlotStopOffset=[150,200,500,1500,100,50,100,500,100,2000,-100];%in ms


%Statistics
for i=1:length(AlignedRaster)
    BlockSel=ismember(UniqueBlockType,[1,11]);
    currPSTH=MeanPSTH_FR_currEvent(BlockSel,i);
    if length(currPSTH)>=2
    [p_Event_PSTH{i},h_Event_PSTH{i}]=ttest2_group(currPSTH{1},currPSTH{2},'both');
    else
       p_Event_PSTH{i}=NaN;h_Event_PSTH{i} =NaN;
    end
    
    
    
end
%[p_Stable_PSTH,h_Stable_PSTH]=ttest2_group(ResponseOn_HighRT_StableScene,ResponseOn_LowRT_StableScene,'both');
%[p_Flexible_PSTH,h_Flexible_PSTH]=ttest2_group(ResponseOn_HighRT_FlexibleScene,ResponseOn_LowRT_FlexibleScene,'both');



%Loop overt the scenes
%SceneName
%SegEventStr={'ExtraFPOn','ExtraFPOff','SceneOn','SceneOff','FixOff','CueOn','TargetOnset','TargetOffset','Reward','Spon'};
UniqueScenes=unique(SceneName);
SceneOnEventFlag=ismember(SegEventStr,'SceneOn');
SceneOffEventFlag=ismember(SegEventStr,'SceneOff');



 

ResponseOn=MeanPSTH_FR{SceneOnEventFlag};
ResponseOff=MeanPSTH_FR{SceneOffEventFlag};
Time_PSTH_SceneOn=nanmean(MeanPSTH_Time{SceneOnEventFlag},1);
Time_PSTH_SceneOff=nanmean(MeanPSTH_Time{SceneOffEventFlag},1);


%Pre Scene On Spon
%PreSceneOnSpon:Response 500ms before scene onset
PreSceneTimeInterval=find(Time_PSTH_SceneOn<0 &Time_PSTH_SceneOn>-500);

AlignedRaster_SceneOn=AlignedRaster{SceneOnEventFlag};
SelectRaster_PreSceneOn=AlignedRaster_SceneOn<0 &AlignedRaster_SceneOn>-500;
PreSceneOnSpon=sum(SelectRaster_PreSceneOn,2)/500*1000;
SceneOn500Resp=sum(AlignedRaster_SceneOn>0 &AlignedRaster_SceneOn<SceneOnAveTime,2)/SceneOnAveTime*1000;
%{
%PreSceneOnSpon by direct mean
PreSceneResponseSel=ResponseOn(:,[PreSceneTimeInterval(1):PreSceneTimeInterval(end)]);
PreSceneOnSpon_DirectMean=nanmean(PreSceneResponseSel,2);
PreSceneOnSpon=PreSceneOnSpon_DirectMean;
%}





%Examine whether it's a scene responsive neuron
%Compare the scene responses with PreSceneOn Spon
EachSceneOn=MeanFRRaw{SceneOnEventFlag};

%Use the relative responses
%Not remove spon anymore
%ResponseOn=ResponseOn-PreSceneOnSpon;%Response changed to the responses relative to the spon
%ResponseOff=ResponseOff-PreSceneOnSpon;

%Average PSTH
MeanPSTH_FR_SceneOn_Flexible_PSTH=nanmean(ResponseOn(BlockType==1,:),1);
MeanPSTH_FR_SceneOn_Flexible_PSTH_SEM=nanstd(ResponseOn(BlockType==1,:),[],1)./sqrt(sum(~isnan(ResponseOn(BlockType==1,:)),1));

MeanPSTH_FR_SceneOn_Stable_PSTH=nanmean(ResponseOn(BlockType==11,:),1);
MeanPSTH_FR_SceneOn_Stable_PSTH_SEM=nanstd(ResponseOn(BlockType==11,:),[],1)./sqrt(sum(~isnan(ResponseOn(BlockType==11,:)),1));


%For DataStorage and export
%MeanPSTH_FR_SceneOn_Flexible=nanmean(MeanPSTH_FR_currEvent{UniqueBlockType==1,SceneOnEventFlag}-PreSceneOnSpon(BlockType==1),1);
%MeanPSTH_FR_SceneOn_Stable=nanmean(MeanPSTH_FR_currEvent{UniqueBlockType==11,SceneOnEventFlag}-PreSceneOnSpon(BlockType==11),1);

MeanPSTH_FR_SceneOn_Flexible=nanmean(MeanPSTH_FR_currEvent{UniqueBlockType==1,SceneOnEventFlag},1);
MeanPSTH_FR_SceneOn_Stable=nanmean(MeanPSTH_FR_currEvent{UniqueBlockType==11,SceneOnEventFlag},1);

%SemPSTH_FR_SceneOn_Flexible=nanstd(MeanPSTH_FR_currEvent{UniqueBlockType==1,SceneOnEventFlag}-PreSceneOnSpon(BlockType==1),[],1)./sqrt(sum(~isnan(MeanPSTH_FR_currEvent{UniqueBlockType==1,SceneOnEventFlag}-PreSceneOnSpon(BlockType==1)),1));
%SemPSTH_FR_SceneOn_Stable=nanstd(MeanPSTH_FR_currEvent{UniqueBlockType==11,SceneOnEventFlag}-PreSceneOnSpon(BlockType==11),[],1)./sqrt(sum(~isnan(MeanPSTH_FR_currEvent{UniqueBlockType==11,SceneOnEventFlag}-PreSceneOnSpon(BlockType==11)),1));
SemPSTH_FR_SceneOn_Flexible=nanstd(MeanPSTH_FR_currEvent{UniqueBlockType==1,SceneOnEventFlag},[],1)./sqrt(sum(~isnan(MeanPSTH_FR_currEvent{UniqueBlockType==1,SceneOnEventFlag}),1));
SemPSTH_FR_SceneOn_Stable=nanstd(MeanPSTH_FR_currEvent{UniqueBlockType==11,SceneOnEventFlag},[],1)./sqrt(sum(~isnan(MeanPSTH_FR_currEvent{UniqueBlockType==11,SceneOnEventFlag}),1));

%Few Num Correction
NumEachCol_Flexible=sum(~isnan(MeanPSTH_FR_currEvent{UniqueBlockType==1,SceneOnEventFlag}),1);
MeanPSTH_FR_SceneOn_Flexible(NumEachCol_Flexible<3)=NaN;
SemPSTH_FR_SceneOn_Flexible(NumEachCol_Flexible<3)=NaN;

NumEachCol_Stable=sum(~isnan(MeanPSTH_FR_currEvent{UniqueBlockType==11,SceneOnEventFlag}),1);
MeanPSTH_FR_SceneOn_Stable(NumEachCol_Stable<3)=NaN;
SemPSTH_FR_SceneOn_Stable(NumEachCol_Stable<3)=NaN;


%MeanPSTH_FR_SceneOff_Flexible=nanmean(MeanPSTH_FR_currEvent{UniqueBlockType==1,SceneOffEventFlag}-PreSceneOnSpon(BlockType==1),1);
%MeanPSTH_FR_SceneOff_Stable=nanmean(MeanPSTH_FR_currEvent{UniqueBlockType==11,SceneOffEventFlag}-PreSceneOnSpon(BlockType==11),1);

MeanPSTH_FR_SceneOff_Flexible=nanmean(MeanPSTH_FR_currEvent{UniqueBlockType==1,SceneOffEventFlag},1);
MeanPSTH_FR_SceneOff_Stable=nanmean(MeanPSTH_FR_currEvent{UniqueBlockType==11,SceneOffEventFlag},1);

%SemPSTH_FR_SceneOff_Flexible=nanstd(MeanPSTH_FR_currEvent{UniqueBlockType==1,SceneOffEventFlag}-PreSceneOnSpon(BlockType==1),[],1)./sqrt(sum(~isnan(MeanPSTH_FR_currEvent{UniqueBlockType==1,SceneOffEventFlag}-PreSceneOnSpon(BlockType==1)),1));
%SemPSTH_FR_SceneOff_Stable=nanstd(MeanPSTH_FR_currEvent{UniqueBlockType==11,SceneOffEventFlag}-PreSceneOnSpon(BlockType==11),[],1)./sqrt(sum(~isnan(MeanPSTH_FR_currEvent{UniqueBlockType==11,SceneOffEventFlag}-PreSceneOnSpon(BlockType==11)),1));
SemPSTH_FR_SceneOff_Flexible=nanstd(MeanPSTH_FR_currEvent{UniqueBlockType==1,SceneOffEventFlag},[],1)./sqrt(sum(~isnan(MeanPSTH_FR_currEvent{UniqueBlockType==1,SceneOffEventFlag}),1));
SemPSTH_FR_SceneOff_Stable=nanstd(MeanPSTH_FR_currEvent{UniqueBlockType==11,SceneOffEventFlag},[],1)./sqrt(sum(~isnan(MeanPSTH_FR_currEvent{UniqueBlockType==11,SceneOffEventFlag}),1));
     

%Few Num Correction
NumEachCol_Flexible=sum(~isnan(MeanPSTH_FR_currEvent{UniqueBlockType==1,SceneOffEventFlag}),1);
MeanPSTH_FR_SceneOff_Flexible(NumEachCol_Flexible<3)=NaN;
SemPSTH_FR_SceneOff_Flexible(NumEachCol_Flexible<3)=NaN;

NumEachCol_Stable=sum(~isnan(MeanPSTH_FR_currEvent{UniqueBlockType==11,SceneOffEventFlag}),1);
MeanPSTH_FR_SceneOff_Stable(NumEachCol_Stable<3)=NaN;
SemPSTH_FR_SceneOff_Stable(NumEachCol_Stable<3)=NaN;


%Average responses of scene on difference in each blocks
%RawSceneOnAverage=RawFR_currEvent(:,SceneOnEventFlag);
%SceneOnAverage=MeanFR_currEvent(:,SceneOnEventFlag);
%[h,p_SceneDiffBetweenBlocks]=ttest2(RawSceneOnAverage{1},RawSceneOnAverage{2});

SceneOnAverageFlexible=EachSceneOn(BlockType==1);%-PreSceneOnSpon(BlockType==1);
SceneOnAverageStable=EachSceneOn(BlockType==11);%-PreSceneOnSpon(BlockType==11);

SceneOnAverage=[nanmean(SceneOnAverageFlexible),nanmean(SceneOnAverageStable)];
[h,p_SceneDiffBetweenBlocks]=ttest2(SceneOnAverageFlexible,SceneOnAverageStable);

%Standard devidation for normalization
%SD_RelSceneOn=std(EachSceneOn-PreSceneOnSpon); 
%Mean_RelSceneOn=mean(EachSceneOn-PreSceneOnSpon);
SD_RelSceneOn=std(EachSceneOn); 
Mean_RelSceneOn=mean(EachSceneOn);

for i=1:length(UniqueScenes)
    Sel=(SceneName==UniqueScenes(i));
    %Sel=(SceneName==1);
%    Resp_tmp=nanmean(ResponseOn(Sel,IntervalSel),2);
%    Scene_Resp_Mean_Whole(i)=nanmean(Resp_tmp);
%    Scene_Resp_Sem_Whole(i)=nanstd(Resp_tmp)/sqrt(sum(~isnan(Resp_tmp)));
    RespOn_tmp_dynamics=ResponseOn(Sel,:);
    SceneOn_Resp_Mean_Each{i}=nanmean( RespOn_tmp_dynamics,1);
    SceneOn_Resp_Sem_Each{i}=nanstd( RespOn_tmp_dynamics,[],1)./sqrt(sum(~isnan( RespOn_tmp_dynamics),1));
    
    RespOff_tmp_dynamics=ResponseOff(Sel,:);
    SceneOff_Resp_Mean_Each{i}=nanmean( RespOff_tmp_dynamics,1);
    SceneOff_Resp_Sem_Each{i}=nanstd( RespOff_tmp_dynamics,[],1)./sqrt(sum(~isnan( RespOff_tmp_dynamics),1));
    
    %Num correction (set number of responses too few to be NaN)
      NumOfValidOn=sum(~isnan(RespOn_tmp_dynamics),1);
      RemovalFlagOn=NumOfValidOn<3;
      SceneOn_Resp_Mean_Each{i}(RemovalFlagOn)=NaN;
      SceneOn_Resp_Sem_Each{i}(RemovalFlagOn)=NaN;
      
      NumOfValidOff=sum(~isnan(RespOff_tmp_dynamics),1);
      RemovalFlagOff=NumOfValidOff<3;
      
      SceneOff_Resp_Mean_Each{i}(RemovalFlagOff)=NaN;
      SceneOff_Resp_Sem_Each{i}(RemovalFlagOff)=NaN;
      
      
      
      
      
      
      
      
      %Total mean for each sceneon/off responses;
      SceneOn_Resp_tmp=MeanFRRaw{SceneOnEventFlag}(Sel);%-PreSceneOnSpon(Sel);
      SceneOff_Resp_tmp=MeanFRRaw{SceneOffEventFlag}(Sel);%-PreSceneOnSpon(Sel);
      
      SceneOn_Resp_Mean_Whole(i)=nanmean(SceneOn_Resp_tmp);
      SceneOn_Resp_Sem_Whole(i)=nanstd(SceneOn_Resp_tmp)/sqrt(sum(~isnan(SceneOn_Resp_tmp)));
      
      SceneOff_Resp_Mean_Whole(i)=nanmean(SceneOff_Resp_tmp);
      SceneOff_Resp_Sem_Whole(i)=nanstd(SceneOff_Resp_tmp)/sqrt(sum(~isnan(SceneOff_Resp_tmp)));
      
       %Average Response on this scene compared to spon
      [h,p_SceneOnWithSponWhole(i)]=ttest(PreSceneOnSpon(Sel),EachSceneOn(Sel));
     % [h,p_SceneOnWithSpon(i)]=ttest2(PreSceneOnSpon,EachSceneOn(Sel));
      SceneOnMeanMinusSponWhole(i)=nanmean(EachSceneOn(Sel));%-nanmean(PreSceneOnSpon(Sel));
      
      
      %Choose the average fr for the first 500ms time window to represent the mean 
      
      SceneOnRaster_CurrScene=AlignedRaster_SceneOn(Sel,:);
      ChooseSceneOnSpikeFlag=SceneOnRaster_CurrScene>0 &SceneOnRaster_CurrScene<SceneOnAveTime;
      SceneOn_500Mean=sum(ChooseSceneOnSpikeFlag,2)/SceneOnAveTime*1000;
       [h,p_SceneOn500WithSpon(i)]=ttest(PreSceneOnSpon(Sel),SceneOn_500Mean);
       
      SceneOn500MeanMinusSpon(i)=nanmean( SceneOn_500Mean);%-nanmean(PreSceneOnSpon(Sel));
     
     
      
      
      %Rearching time for each scene
      %First look at the first searching time
      FirstSearchingTime_Each=SearchingTime(Sel,1);
      %Sorting according to the first searching time
      [sortedRT,sorting_index]=sort(FirstSearchingTime_Each,'descend');
       UpSel=sorting_index(1:floor(length(sorting_index)/2));
       DownSel=sorting_index(floor(length(sorting_index)/2)+1:end);
       
       ResponseOn_UpSel=RespOn_tmp_dynamics(UpSel,:);
       ResponseOn_DownSel=RespOn_tmp_dynamics(DownSel,:);
       
       ResponseOnMean_HighRT_currScene(i,:)=nanmean(ResponseOn_UpSel,1);
       ResponseOnSem_HighRT_currScene(i,:)=nanstd(ResponseOn_UpSel,[],1)./sqrt(sum(~isnan(ResponseOn_UpSel),1));
      
       ResponseOnMean_LowRT_currScene(i,:)=nanmean(ResponseOn_DownSel,1);
       ResponseOnSem_LowRT_currScene(i,:)=nanstd(ResponseOn_DownSel,[],1)./sqrt(sum(~isnan(ResponseOn_DownSel),1));
       
       
       ResponseOn_HighRT_currScene{i}=ResponseOn_UpSel;
       ResponseOn_LowRT_currScene{i}=ResponseOn_DownSel;
       
       HighRT_currScene{i}=sortedRT(UpSel);
       LowRT_currScene{i}=sortedRT(DownSel);
       
       
       %Rearching time for the second location
       SecondSearchingTime_Each=SearchingTime(Sel,2);
      %Sorting according to the first searching time
      [sortedRTSecond,sorting_indexSecond]=sort(SecondSearchingTime_Each,'descend');
      
      %Fix off responses separated for each scene and For each fractal order
      SelLine=(SceneName_Line==UniqueScenes(i));
      for j=1:length(ObjOrderSep)+1
          if j==length(ObjOrderSep)+1
             selOrder=TargetOrder_Line>ObjOrderSep(end);
            % order_string{i,j}=sprintf('Fractal Order> %d (N=%d)',ObjOrderSep(end),sum(selOrder));
             
          else
            selOrder=TargetOrder_Line==ObjOrderSep(j);
          %  order_string{i,j}=sprintf('Fractal Order %d (N=%d)',ObjOrderSep(j),sum(selOrder));
            
          end
          sel_here=selOrder&SelLine;
          tmp=PSTH_fixoffAll(sel_here,:);
 
          FixOffPSTH_Scene_OrderSep_mean{i,j}=nanmean(tmp,1);
          FixOffPSTH_Scene_OrderSep_sem{i,j}=nanstd(tmp,[],1)./sqrt(sum(~isnan(tmp),1));
          
          
          
          
          
          
      end
      
      
      
      
      
       
       
  
    
end %End of if for unique scenes

SceneOnMarker=nanmean(RelativeTimeMarker{SceneOnEventFlag},1);
SceneOffMarker=nanmean(RelativeTimeMarker{SceneOffEventFlag},1);
PlotSceneOnMarkerStart=PlotStartOffset(SceneOnEventFlag);
PlotSceneOnMarkerStop=PlotStopOffset(SceneOnEventFlag);

PlotSceneOffMarkerStart=PlotStartOffset(SceneOffEventFlag);
PlotSceneOffMarkerStop=PlotStopOffset(SceneOnEventFlag);

%Interval to calculate sig prop
PSTH_SelInterval_tmp=SceneOnMarker;
PSTH_SelInterval=Time_PSTH_SceneOn>=PSTH_SelInterval_tmp(1) & Time_PSTH_SceneOn<=PSTH_SelInterval_tmp(2);


%% PSTH R^2 with ANOVA for all scenes and each condition
%Overall
[pSceneDiffWithTime,rsquare_psth,~]=anova1_group_stats(ResponseOn,SceneName,'off');
%Few number correction
NumOfValidOn=sum(~isnan(ResponseOn),1);
RemovalFlagOn=NumOfValidOn<3*length(UniqueScenes);
pSceneDiffWithTime(RemovalFlagOn)=NaN;
rsquare_psth(RemovalFlagOn)=NaN;



h_sig_ssi_withtime=floor(pSceneDiffWithTime<0.05);
h_sig_ssi_withtime(h_sig_ssi_withtime==0)=NaN;


%Stable
[pSceneDiffWithTime_Stable,rsquare_psth_Stable,~]=anova1_group_stats(ResponseOn(ismember(SceneName,UniqueScene_Stable),:),SceneName(ismember(SceneName,UniqueScene_Stable)),'off');
%Few number correction
NumOfValidOn=sum(~isnan(ResponseOn(ismember(SceneName,UniqueScene_Stable),:)),1);
RemovalFlagOn=NumOfValidOn<3*length(UniqueScene_Stable);
pSceneDiffWithTime_Stable(RemovalFlagOn)=NaN;
rsquare_psth_Stable(RemovalFlagOn)=NaN;

h_sig_ssi_withtime_Stable=floor(pSceneDiffWithTime_Stable<0.05);
h_sig_ssi_withtime_Stable(h_sig_ssi_withtime_Stable==0)=NaN;

PropOfSig_Stable=sum(pSceneDiffWithTime_Stable(PSTH_SelInterval)<0.05)/sum(PSTH_SelInterval);


%Flexible
[pSceneDiffWithTime_Flexible,rsquare_psth_Flexible,~]=anova1_group_stats(ResponseOn(ismember(SceneName,UniqueScene_Flexible),:),SceneName(ismember(SceneName,UniqueScene_Flexible)),'off');

%Few number correction
NumOfValidOn=sum(~isnan(ResponseOn(ismember(SceneName,UniqueScene_Flexible),:)),1);
RemovalFlagOn=NumOfValidOn<3*length(UniqueScene_Flexible);
pSceneDiffWithTime_Flexible(RemovalFlagOn)=NaN;
rsquare_psth_Flexible(RemovalFlagOn)=NaN;


h_sig_ssi_withtime_Flexible=floor(pSceneDiffWithTime_Flexible<0.05);
h_sig_ssi_withtime_Flexible(h_sig_ssi_withtime_Flexible==0)=NaN;

PropOfSig_Flexible=sum(pSceneDiffWithTime_Flexible(PSTH_SelInterval)<0.05)/sum(PSTH_SelInterval);




%%
%Criterium to select out scene responsive neurons Part1
%Use pair t test to examine the significance between preSceenResponse and the Overall scene on responses 
[h,pSceneOnWhole]=ttest(PreSceneOnSpon,EachSceneOn);
SceneOnDifference_WholeMean=nanmean(EachSceneOn)-nanmean(PreSceneOnSpon);

[h,pSceneOn500]=ttest(PreSceneOnSpon,SceneOn500Resp);
SceneOnDifference_500Mean=nanmean(SceneOn500Resp)-nanmean(PreSceneOnSpon);
%Use former 500ms as a representative responses
pSceneOn=pSceneOn500;
SceneOnDifference_OverallMean=SceneOnDifference_500Mean;
p_SceneOnWithSpon=p_SceneOn500WithSpon;
SceneOnMeanMinusSpon=SceneOn500MeanMinusSpon;%For each scene

%Add cell type string to show on the figure
  %SceneOnAverage,p_SceneDiffBetweenBlocks
  if pSceneOn<0.05
      if SceneOnDifference_OverallMean>0
          Responsetype='Excitatory';
      else
          Responsetype='Inhibitory';
      end
      
 
  else
    Responsetype='(No Scene Response)';  
  end
  
   if p_SceneDiffBetweenBlocks<0.05
      if diff(SceneOnAverage)>0
          celltype=strcat(Responsetype,'(Sta Dominance)');
      else
          celltype=strcat(Responsetype,'(Fle Dominance)');
         
          
      end
  else
      celltype=strcat(Responsetype,'(No Dominance)');
  end
  
  
  

%Scene Selective Index(SI) for each type of scenes 
%Stable Scene
    Stable_Unique_Sel=ismember(UniqueScenes,UniqueScene_Stable);
    N_Stable_Unique=sum(Stable_Unique_Sel);
    Resp_Stable_Unique_Mean= SceneOn_Resp_Mean_Whole(Stable_Unique_Sel);
    SI_Stable=(N_Stable_Unique-sum(Resp_Stable_Unique_Mean/max(Resp_Stable_Unique_Mean)))/(N_Stable_Unique-1);
    %Flexible Scene
    Flexible_Unique_Sel=ismember(UniqueScenes,UniqueScene_Flexible);
    N_Flexible_Unique=sum(Flexible_Unique_Sel);
    Resp_Flexible_Unique_Mean= SceneOn_Resp_Mean_Whole(Flexible_Unique_Sel);
    SI_Flexible=(N_Flexible_Unique-sum(Resp_Flexible_Unique_Mean/max(Resp_Flexible_Unique_Mean)))/(N_Flexible_Unique-1);
%}

%Use anova to statistically test scene selectivity
SceneOnResponseAll=MeanFRRaw{SceneOnEventFlag};
[p_SceneOnAnova,t,stats]=anova1(SceneOnResponseAll,SceneName,'off');
%[c,m,h,gnames] = multcompare(stats);
%R square
ss=t{2,2};
ts=t{4,2};   
R2_Whole=ss/ts;

%Stable scenes
%Separate for each block
Stable_Raw_Sel=ismember(SceneName,UniqueScene_Stable);
SceneOnResponseStable_Raw=MeanFRRaw{SceneOnEventFlag}(Stable_Raw_Sel);
SceneName_Stable_Raw=SceneName(Stable_Raw_Sel);
[p_SceneOnAnova_Stable,t,stats]=anova1(SceneOnResponseStable_Raw,SceneName_Stable_Raw,'off');
%R square
ss=t{2,2};
ts=t{4,2};   
R2_Stable=ss/ts;

%Separate by half for comparison with passive
UniqueScene_StableHalf1=1010:1013;
UniqueScene_StableHalf2=1014:1017;

Stable_Raw_SelHalf=ismember(SceneName,UniqueScene_StableHalf1);
SceneOnResponseStable_Raw=MeanFRRaw{SceneOnEventFlag}(Stable_Raw_SelHalf);
SceneName_Stable_Raw=SceneName(Stable_Raw_SelHalf);
[p_SceneOnAnova_StableHalf1,t,stats]=anova1(SceneOnResponseStable_Raw,SceneName_Stable_Raw,'off');
%R square
ss=t{2,2};
ts=t{4,2};   
R2_StableHalf1=ss/ts;

Stable_Raw_SelHalf=ismember(SceneName,UniqueScene_StableHalf2);
SceneOnResponseStable_Raw=MeanFRRaw{SceneOnEventFlag}(Stable_Raw_SelHalf);
SceneName_Stable_Raw=SceneName(Stable_Raw_SelHalf);
[p_SceneOnAnova_StableHalf2,t,stats]=anova1(SceneOnResponseStable_Raw,SceneName_Stable_Raw,'off');
%R square
ss=t{2,2};
ts=t{4,2};   
R2_StableHalf2=ss/ts;


%Flexible scenes
Flexible_Raw_Sel=ismember(SceneName,UniqueScene_Flexible);
SceneOnResponseFlexible_Raw=MeanFRRaw{SceneOnEventFlag}(Flexible_Raw_Sel);
SceneName_Flexible_Raw=SceneName(Flexible_Raw_Sel);
[p_SceneOnAnova_Flexible,t,stats]=anova1(SceneOnResponseFlexible_Raw,SceneName_Flexible_Raw,'off');
%R square
ss=t{2,2};
ts=t{4,2};   
R2_Flexible=ss/ts;

UniqueScene_FlexibleHalf1=1000:1003;
UniqueScene_FlexibleHalf2=1004:1007;

Flexible_Raw_SelHalf=ismember(SceneName,UniqueScene_FlexibleHalf1);
SceneOnResponseFlexible_Raw=MeanFRRaw{SceneOnEventFlag}(Flexible_Raw_SelHalf);
SceneName_Flexible_Raw=SceneName(Flexible_Raw_SelHalf);
[p_SceneOnAnova_FlexibleHalf1,t,stats]=anova1(SceneOnResponseFlexible_Raw,SceneName_Flexible_Raw,'off');
%R square
ss=t{2,2};
ts=t{4,2};   
R2_FlexibleHalf1=ss/ts;

Flexible_Raw_SelHalf=ismember(SceneName,UniqueScene_FlexibleHalf2);
SceneOnResponseFlexible_Raw=MeanFRRaw{SceneOnEventFlag}(Flexible_Raw_SelHalf);
SceneName_Flexible_Raw=SceneName(Flexible_Raw_SelHalf);
[p_SceneOnAnova_FlexibleHalf2,t,stats]=anova1(SceneOnResponseFlexible_Raw,SceneName_Flexible_Raw,'off');
%R square
ss=t{2,2};
ts=t{4,2};   
R2_FlexibleHalf2=ss/ts;




%%
%Criterium to select out scene responsive neurons Part2

%Selecte Out the significantly above spon responsive scenes
%p_SceneOnWithSpon=p_SceneOn500WithSpon;
%SceneOnMeanMinusSpon=SceneOn500MeanMinusSpon;
SceneOnSigScenes=UniqueScenes(p_SceneOnWithSpon<0.05);
SceneOnSigSceneRelResponses=SceneOnMeanMinusSpon(p_SceneOnWithSpon<0.05);
NumOfSigScenes=sum(p_SceneOnWithSpon<0.05);
%Find out the most responsive scenes and least responsive scenes

if ~isempty(SceneOnSigScenes)
    [maxval,MostResponsiveSceneIndex]=max(abs(SceneOnSigSceneRelResponses));
    MostResponsiveScene=SceneOnSigScenes(MostResponsiveSceneIndex);
    MostResponsiveSceenResponses=SceneOnSigSceneRelResponses(MostResponsiveSceneIndex);
    
     [minval,LeastResponsiveSceneIndex]=min(abs(SceneOnSigSceneRelResponses));
     LeastResponsiveScene=SceneOnSigScenes(LeastResponsiveSceneIndex);
     LeastResponsiveSceenResponses=SceneOnSigSceneRelResponses(LeastResponsiveSceneIndex);
     
     %Average the PSTH for non responsive scene responses
NonResponsiveSceneSelector=~(p_SceneOnWithSpon<0.05);
NonResponsiveScenes=UniqueScenes(NonResponsiveSceneSelector);

PSTH_SceneOn_NonResponsiveScenes=ResponseOn(ismember(SceneName,NonResponsiveScenes),:);
PSTH_SceneOn_NonResponsiveScenes_Mean=nanmean(PSTH_SceneOn_NonResponsiveScenes,1);
PSTH_SceneOn_NonResponsiveScenes_Sem=nanstd(PSTH_SceneOn_NonResponsiveScenes,[],1)./sqrt(sum(~isnan(PSTH_SceneOn_NonResponsiveScenes),1));

%PSTH for Most responsive and Least Responsive Scene responses, and
%remaining responsive neurons
PSTH_SceneOn_MostResponsiveScenes_Mean=SceneOn_Resp_Mean_Each{ismember(UniqueScenes, MostResponsiveScene)};
PSTH_SceneOn_MostResponsiveScenes_Sem=SceneOn_Resp_Sem_Each{ismember(UniqueScenes, MostResponsiveScene)};

PSTH_SceneOn_LeastResponsiveScenes_Mean=SceneOn_Resp_Mean_Each{ismember(UniqueScenes, LeastResponsiveScene)};
PSTH_SceneOn_LeastResponsiveScenes_Sem=SceneOn_Resp_Sem_Each{ismember(UniqueScenes, LeastResponsiveScene)};

%Scene responses other than the most responsive scene
RemainingResponsiveScenesIndex=ismember(SceneName, SceneOnSigScenes)&(~ismember(SceneName, MostResponsiveScene));
PSTH_SceneOn_RemainingResponsiveScenes=ResponseOn(RemainingResponsiveScenesIndex,:);

PSTH_SceneOn_RemainingResponsiveScenes_Mean=nanmean(PSTH_SceneOn_RemainingResponsiveScenes,1);
PSTH_SceneOn_RemainingResponsiveScenes_Sem=nanstd(PSTH_SceneOn_RemainingResponsiveScenes,[],1)./sqrt(sum(~isnan(PSTH_SceneOn_RemainingResponsiveScenes),1));

%PSTH separate for excitotary responsive scenes and inhibitory responsive
%scenes
ExcitatoryResponsiveSceneIndex=SceneOnSigScenes(SceneOnSigSceneRelResponses>0);
PSTH_SceneOn_ExcitatoryResponsiveScene=ResponseOn(ismember(SceneName,ExcitatoryResponsiveSceneIndex),:);

PSTH_SceneOn_ExcitatoryResponsiveScene_Mean=nanmean(PSTH_SceneOn_ExcitatoryResponsiveScene,1);
PSTH_SceneOn_ExcitatoryResponsiveScene_Sem=nanstd(PSTH_SceneOn_ExcitatoryResponsiveScene,[],1)./sqrt(sum(~isnan(PSTH_SceneOn_ExcitatoryResponsiveScene),1));

InhibitoryResponsiveSceneIndex=SceneOnSigScenes(SceneOnSigSceneRelResponses<0);
PSTH_SceneOn_InhibitoryResponsiveScene=ResponseOn(ismember(SceneName,InhibitoryResponsiveSceneIndex),:);

PSTH_SceneOn_InhibitoryResponsiveScene_Mean=nanmean(PSTH_SceneOn_InhibitoryResponsiveScene,1);
PSTH_SceneOn_InhibitoryResponsiveScene_Sem=nanstd(PSTH_SceneOn_InhibitoryResponsiveScene,[],1)./sqrt(sum(~isnan(PSTH_SceneOn_InhibitoryResponsiveScene),1));

     
    
else
    MostResponsiveScene=NaN;
    MostResponsiveSceenResponses=NaN;
    LeastResponsiveScene=NaN;
    LeastResponsiveSceenResponses=NaN;
    
    
    PSTH_SceneOn_NonResponsiveScenes_Mean=NaN*ones(1,size(Time_PSTH_SceneOn,2));
    PSTH_SceneOn_NonResponsiveScenes_Sem=NaN*ones(1,size(Time_PSTH_SceneOn,2));
    
    PSTH_SceneOn_MostResponsiveScenes_Mean=NaN*ones(1,size(Time_PSTH_SceneOn,2));
    PSTH_SceneOn_MostResponsiveScenes_Sem=NaN*ones(1,size(Time_PSTH_SceneOn,2));
    
    PSTH_SceneOn_RemainingResponsiveScenes_Mean=NaN*ones(1,size(Time_PSTH_SceneOn,2));
    PSTH_SceneOn_RemainingResponsiveScenes_Sem=NaN*ones(1,size(Time_PSTH_SceneOn,2));
    
    PSTH_SceneOn_ExcitatoryResponsiveScene_Mean=NaN*ones(1,size(Time_PSTH_SceneOn,2));
    PSTH_SceneOn_ExcitatoryResponsiveScene_Sem=NaN*ones(1,size(Time_PSTH_SceneOn,2));
    
    PSTH_SceneOn_InhibitoryResponsiveScene_Mean=NaN*ones(1,size(Time_PSTH_SceneOn,2));
    PSTH_SceneOn_InhibitoryResponsiveScene_Sem=NaN*ones(1,size(Time_PSTH_SceneOn,2));
    
    
    
    
end

%Average of all scene responses
PSTH_SceneOn_Whole_Mean=nanmean(ResponseOn,1);
PSTH_SceneOn_Whole_Sem=nanstd(ResponseOn,[],1)./sqrt(sum(~isnan(ResponseOn),1));


%Separate significant scenes into stable scenes and flexible scenes
StableSigScenes=SceneOnSigScenes(ismember(SceneOnSigScenes,UniqueScene_Stable));
FlexibleSigScenes=SceneOnSigScenes(ismember(SceneOnSigScenes,UniqueScene_Flexible));

StableSigSceneResponseEachScene=SceneOnMeanMinusSpon(ismember(SceneOnSigScenes,UniqueScene_Stable));
FlexibleSigSceneResponseEachScene=SceneOnMeanMinusSpon(ismember(SceneOnSigScenes,UniqueScene_Flexible));


NumOfStableSigScenes=length(StableSigScenes);
NumOfFlexibleSigScenes=length(FlexibleSigScenes);

if ~isempty(StableSigSceneResponseEachScene)
    [val,MostPreferedStableSceneIndex]=max(abs(StableSigSceneResponseEachScene));
    MostPreferedStableScene=StableSigScenes(MostPreferedStableSceneIndex);
    PSTH_MostPreferedStableScene_Mean=SceneOn_Resp_Mean_Each{ismember(UniqueScenes, MostPreferedStableScene)};
    PSTH_MostPreferedStableScene_Sem=SceneOn_Resp_Sem_Each{ismember(UniqueScenes, MostPreferedStableScene)};
    
    %Scene responses other than the most responsive scene
    OtherResponsiveScenesIndex=ismember(SceneName, UniqueScene_Stable)&(~ismember(SceneName, MostPreferedStableScene));
     PSTH_OtherStableScene=ResponseOn(OtherResponsiveScenesIndex,:);

     PSTH_OtherStableScene_Mean=nanmean( PSTH_OtherStableScene,1);
     PSTH_OtherStableScene_Sem=nanstd(PSTH_OtherStableScene,[],1)./sqrt(sum(~isnan(PSTH_OtherStableScene),1));
     %Few Number Correction
     NumValid=sum(~isnan(PSTH_OtherStableScene),1);
     ColValid=NumValid>3;
     
     PSTH_OtherStableScene_Mean(~ColValid)=NaN;
     PSTH_OtherStableScene_Sem(~ColValid)=NaN;
     
   
    
else
    MostPreferedStableScene=NaN;
    PSTH_MostPreferedStableScene_Mean=NaN*ones(1,size(Time_PSTH_SceneOn,2));
    PSTH_MostPreferedStableScene_Sem=NaN*ones(1,size(Time_PSTH_SceneOn,2));
    
     PSTH_OtherStableScene_Mean=MeanPSTH_FR_SceneOn_Stable;
     PSTH_OtherStableScene_Sem=SemPSTH_FR_SceneOn_Stable;
    
    
end

if ~isempty(FlexibleSigSceneResponseEachScene)
    [val,MostPreferedFlexibleSceneIndex]=max(abs(FlexibleSigSceneResponseEachScene));
    MostPreferedFlexibleScene=FlexibleSigScenes(MostPreferedFlexibleSceneIndex);
    PSTH_MostPreferedFlexibleScene_Mean=SceneOn_Resp_Mean_Each{ismember(UniqueScenes, MostPreferedFlexibleScene)};
    PSTH_MostPreferedFlexibleScene_Sem=SceneOn_Resp_Sem_Each{ismember(UniqueScenes, MostPreferedFlexibleScene)};
    
     
    %Scene responses other than the most responsive scene
    OtherResponsiveScenesIndex=ismember(SceneName, UniqueScene_Flexible)&(~ismember(SceneName, MostPreferedFlexibleScene));
     PSTH_OtherFlexibleScene=ResponseOn(OtherResponsiveScenesIndex,:);

     PSTH_OtherFlexibleScene_Mean=nanmean( PSTH_OtherFlexibleScene,1);
     PSTH_OtherFlexibleScene_Sem=nanstd(PSTH_OtherFlexibleScene,[],1)./sqrt(sum(~isnan(PSTH_OtherFlexibleScene),1));
    
     %Few Number Correction
     NumValid=sum(~isnan(PSTH_OtherFlexibleScene),1);
     ColValid=NumValid>3;
     
     PSTH_OtherFlexibleScene_Mean(~ColValid)=NaN;
     PSTH_OtherFlexibleScene_Sem(~ColValid)=NaN;
    
    
    
else
    MostPreferedFlexibleScene=NaN;
    PSTH_MostPreferedFlexibleScene_Mean=NaN*ones(1,size(Time_PSTH_SceneOn,2));
    PSTH_MostPreferedFlexibleScene_Sem=NaN*ones(1,size(Time_PSTH_SceneOn,2));
    
     PSTH_OtherFlexibleScene_Mean=MeanPSTH_FR_SceneOn_Flexible;
     PSTH_OtherFlexibleScene_Sem=SemPSTH_FR_SceneOn_Flexible;
     
end







%% Next Section: Relate neural resonses with behaviors--RT
%SearchingTime
%FirstSearchingTime_EachScene=cell2mat(FirstSearchingTime_Each);
%Scene response sorting according to the first searching time

ResponseOn_HighRT_StableScene=cell2mat(ResponseOn_HighRT_currScene(ismember(UniqueScenes,UniqueScene_Stable))');
ResponseOn_LowRT_StableScene=cell2mat(ResponseOn_LowRT_currScene(ismember(UniqueScenes,UniqueScene_Stable))');

ResponseOn_HighRT_StableScene_mean=nanmean(ResponseOn_HighRT_StableScene,1);
ResponseOn_LowRT_StableScene_mean=nanmean(ResponseOn_LowRT_StableScene,1);

ResponseOn_HighRT_StableScene_sem=nanstd(ResponseOn_HighRT_StableScene,[],1)./sqrt(sum(~isnan(ResponseOn_HighRT_StableScene),1));
ResponseOn_LowRT_StableScene_sem=nanstd(ResponseOn_LowRT_StableScene,[],1)./sqrt(sum(~isnan(ResponseOn_LowRT_StableScene),1));

ResponseOn_HighRT_FlexibleScene=cell2mat(ResponseOn_HighRT_currScene(ismember(UniqueScenes,UniqueScene_Flexible))');
ResponseOn_LowRT_FlexibleScene=cell2mat(ResponseOn_LowRT_currScene(ismember(UniqueScenes,UniqueScene_Flexible))');

ResponseOn_HighRT_FlexibleScene_mean=nanmean(ResponseOn_HighRT_FlexibleScene,1);
ResponseOn_LowRT_FlexibleScene_mean=nanmean(ResponseOn_LowRT_FlexibleScene,1);

ResponseOn_HighRT_FlexibleScene_sem=nanstd(ResponseOn_HighRT_FlexibleScene,[],1)./sqrt(sum(~isnan(ResponseOn_HighRT_FlexibleScene),1));
ResponseOn_LowRT_FlexibleScene_sem=nanstd(ResponseOn_LowRT_FlexibleScene,[],1)./sqrt(sum(~isnan(ResponseOn_LowRT_FlexibleScene),1));

%Statistics
[p_Stable_PSTH,h_Stable_PSTH]=ttest2_group(ResponseOn_HighRT_StableScene,ResponseOn_LowRT_StableScene,'both');
[p_Flexible_PSTH,h_Flexible_PSTH]=ttest2_group(ResponseOn_HighRT_FlexibleScene,ResponseOn_LowRT_FlexibleScene,'both');



%RT statistics
emptyexclude=cellfun(@isempty,HighRT_currScene)';

HighRT_StableScene=cell2mat(HighRT_currScene(ismember(UniqueScenes,UniqueScene_Stable) & ~emptyexclude)');
LowRT_StableScene=cell2mat(LowRT_currScene(ismember(UniqueScenes,UniqueScene_Stable) & ~emptyexclude)');

HighRT_FlexibleScene=cell2mat(HighRT_currScene(ismember(UniqueScenes,UniqueScene_Flexible) & ~emptyexclude)');
LowRT_FlexibleScene=cell2mat(LowRT_currScene(ismember(UniqueScenes,UniqueScene_Flexible) & ~emptyexclude)');

HighRT_Stable_Mean=nanmean(HighRT_StableScene);
LowRT_Stable_Mean=nanmean(LowRT_StableScene);

HighRT_Stable_Sem=nanstd(HighRT_StableScene)/sqrt(length(HighRT_StableScene));
LowRT_Stable_Sem=nanstd(LowRT_StableScene)/sqrt(length(LowRT_StableScene));

HighRT_Flexible_Mean=nanmean(HighRT_FlexibleScene);
LowRT_Flexible_Mean=nanmean(LowRT_FlexibleScene);

HighRT_Flexible_Sem=nanstd(HighRT_FlexibleScene)/sqrt(length(HighRT_FlexibleScene));
LowRT_Flexible_Sem=nanstd(LowRT_FlexibleScene)/sqrt(length(LowRT_FlexibleScene));

[h_Stable_RT,p_Stable_RT]=ttest2(HighRT_StableScene,LowRT_StableScene);

[h_Flexible_RT,p_Flexible_RT]=ttest2(HighRT_FlexibleScene,LowRT_FlexibleScene);



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
 
 figtitlestr{FigureIndex}='SpikeDensityFunctionPSTH';
 fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[100,100, 1500,800],'Name',figtitlestr{FigureIndex});   
    
   for i=1:length(AlignedRaster)
        
         
        for j=1:NumberOfBlocks
            
            %Spike-Count PSTH
           
            currBlock=UniqueBlockType(j);
            
            subplot(ceil((length(AlignedRaster)+1)/2),2,i);
            
            
            curr_Response=nanmean(MeanPSTH_FR_currEvent{j,i},1);
            curr_Response_err=nanstd(MeanPSTH_FR_currEvent{j,i},[],1)./sqrt(sum(~isnan(MeanPSTH_FR_currEvent{j,i}),1));
             %Num correction (set number of responses too few to be NaN)
            NumOfValid=sum(~isnan(MeanPSTH_FR_currEvent{j,i}),1);
            RemovalFlag=NumOfValid<3;
            curr_Response(RemovalFlag)=NaN;
            curr_Response_err(RemovalFlag)=NaN;
            
            plot([0,0],[0,max(curr_Response)],'--k','Color','#929498','LineWidth',1);
            hold on
            %plot(PSTH_Time_currEvent{j,i},curr_Response,'-k','Color',ColorForEachMode(AssociateBlock==currBlock,:),'LineWidth',2);
           % bar(PSTH_Time_currEvent{j,i},nanmean(PSTH_Event_SpikeCount_currEvent{j,i},1),'FaceColor',ColorForEachMode(AssociateBlock==currBlock,:));
            shadedErrorBar(PSTH_Time_currEvent{j,i},curr_Response,curr_Response_err,'lineprops',{ColorForEachMode(AssociateBlock==currBlock,:)},'transparent',1,'patchSaturation',0.5);%PSTH for fix
             set(findall(gca, 'Type', 'Line'),'LineWidth',3);
            hold on 
            %Add event-end marker
            plot([min(RelativeTimeMarker_currEvent{j,i}(:,2)), min(RelativeTimeMarker_currEvent{j,i}(:,2))],[0,max(curr_Response)],'-k','Color','#707276','LineWidth',1.5);
           % SegEventStr={'ExtraFPOn','ExtraFPOff','SceneOn','SceneOff','FixOff','CueOn','TargetOnset','TargetOffset','Reward','Spon'};
          %SegEventEnd=[EXTRAFIXOFF,SHOWFIXCD,SHOWFIXCD,NaN,CUEON,TGTCD,TGTOFF,NaN,RWDOFFCD,NaN ];%Event end; %Automatic screening later
          %Add significance
          if ~isnan(h_Event_PSTH{i})
              sig=area(PSTH_Time_currEvent{j,i},h_Event_PSTH{i}*0.5,'FaceColor','#FFC107','EdgeColor','#FFC107','FaceAlpha',0.5,'EdgeAlpha',0.5);
              set( get( get( sig, 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );
          end
  
         % h_Event_PSTH
            if isnan(SegEventEnd(i))
                PlotStopOffsetPlot(i)=PlotStopOffset(i)/5;
            else
                PlotStopOffsetPlot(i)=PlotStopOffset(i);
            end
            xlim([-PlotStartOffset(i),PlotStopOffsetPlot(i)+nanmean(RelativeTimeMarker_currEvent{j,i}(:,2))]);
            
            %Set ylim
            %Set ylim
            Interval=curr_Response(:,PSTH_Time_currEvent{j,i}<PlotStopOffsetPlot(i)+nanmean(RelativeTimeMarker_currEvent{j,i}(:,2)) & PlotStopOffsetPlot(i)+nanmean(RelativeTimeMarker_currEvent{j,i}(:,2))>-PlotStartOffset(i));
              Interval_max(j)=max( Interval)*1.12;
            
           
            if i==3
             ylabel('Spike Count per s(Hz)');
            end
            if i==length(AlignedRaster)
                xlabel('Time from event onset');
            end
             title(SegEventStr{i});
            set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
            %}
        end
        
         ylim([0,max(max(  Interval_max),1)]);
       
   end
   %Add spon with block
   PlotIndex=1:length(SponResponse);
    subplot(ceil((length(AlignedRaster)+1)/2),2,length(AlignedRaster)+1);
    for i=1:NumberOfBlocks
      currBlock=UniqueBlockType(i);
      Sel=BlockType==currBlock;
      plot(PlotIndex(Sel),SponResponse(Sel),'-k','Color',ColorForEachMode(AssociateBlock==currBlock,:),'LineWidth',3);
      hold on
   
    end

set(gca,'XTick',PlotIndex(1:50:end));
set(gca,'XTickLabel',TrialIndex(1:50:end));
xlabel('Trial Index');
ylabel('Spon');
title('PostTrialSpon with trials');
set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');


  
    %Add data information
  %  data_msg=sprintf(' %s, SpikeChannel: %d, TrialNum: %d ',FileName,SpikeChannel,TrialNum);
    data_msg=strcat(FileName,' SpikeChannel:',num2str(SpikeChannel),';ChannelID:',num2str(ChannelID),';ClusterID:',num2str(ClusterID),';TrialNum: ',num2str(TrialNum),celltype);

    annotation('textbox', [0.17, 0.95, 0.7, 0.04], 'string', data_msg,'FontSize',20,'FontWeight','Bold','LineStyle','none')

    FigureIndex=FigureIndex+1;
    FigureStartNum=FigureStartNum+1; 
    
%%%%%%%%%%%%%%%%%%Next figure: Scene Onset Offset PSTH for each scene

 figtitlestr{FigureIndex}='SceneOnOffPSTH';
 fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[120,100, 1600,800],'Name',figtitlestr{FigureIndex});   
 %First half
for i=1:ceil(length(UniqueScenes)/2)
   
  subplot(2,2,1);
  shadedErrorBar(Time_PSTH_SceneOn,SceneOn_Resp_Mean_Each{i},SceneOn_Resp_Sem_Each{i},'lineprops',{ColorForScene(i,:)},'transparent',1,'patchSaturation',0.5);
  hold on
  %Add time marker
  if i==ceil(length(UniqueScenes)/2)
      hold on
 % plot([SceneOnMarker(1), SceneOnMarker(1)],[0,max(cellfun(@(x,y) max(x+y),SceneOn_Resp_Mean_Each,SceneOn_Resp_Sem_Each))],'-k','Color','#707276','LineWidth',1.5);
 % plot([SceneOnMarker(2), SceneOnMarker(2)],[0,max(cellfun(@(x,y) max(x+y),SceneOn_Resp_Mean_Each,SceneOn_Resp_Sem_Each))],'-k','Color','#707276','LineWidth',1.5);
   h1=plot([SceneOnMarker(1), SceneOnMarker(1)],[0,max(cellfun(@(x,y) max(x+y),SceneOn_Resp_Mean_Each,SceneOn_Resp_Sem_Each))],'--k','LineWidth',1.5);
  h2=plot([SceneOnMarker(2), SceneOnMarker(2)],[0,max(cellfun(@(x,y) max(x+y),SceneOn_Resp_Mean_Each,SceneOn_Resp_Sem_Each))],'--k','LineWidth',1.5);
  
  set( get( get( h1, 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );
 set( get( get( h2, 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );
 
  end
  
  

  set(findall(gca, 'Type', 'Line'),'LineWidth',3);
  xlabel('Time from scene onset(ms)');
  ylabel('Firing Rate(Hz)');
  

  title('SceneOnset');
  set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
   xlim([-PlotSceneOnMarkerStart,PlotSceneOnMarkerStop+SceneOnMarker(2)]);
   
  subplot(2,2,2);
   shadedErrorBar(Time_PSTH_SceneOff,SceneOff_Resp_Mean_Each{i},SceneOff_Resp_Sem_Each{i},'lineprops',{ColorForScene(i,:)},'transparent',1,'patchSaturation',0.5);
  %Add time marker
  if i==ceil(length(UniqueScenes)/2)
    h3=plot([SceneOffMarker(1), SceneOffMarker(1)],[0,max(cellfun(@(x,y) max(x+y),SceneOff_Resp_Mean_Each,SceneOff_Resp_Sem_Each))],'--k','LineWidth',1.5);
  h4=plot([SceneOffMarker(2), SceneOffMarker(2)],[0,max(cellfun(@(x,y) max(x+y),SceneOff_Resp_Mean_Each,SceneOff_Resp_Sem_Each))],'--k','LineWidth',1.5);
  
 set( get( get( h3, 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );
 set( get( get( h4, 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );
  end
  
  
    xlim([-PlotSceneOffMarkerStart,PlotSceneOffMarkerStop+SceneOffMarker(2)]);
    
   set(findall(gca, 'Type', 'Line'),'LineWidth',3);
    xlabel('Time from scene offset(ms)');
  ylabel('Firing Rate(Hz)');

   title('SceneOffset');
   set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
     hold on
    
end
ld=legend([string(UniqueScenes(1:ceil(length(UniqueScenes)/2)))]);
%ld=legend([string(UniqueScene);'SSI']);
ld.Box='off';
ld.Position=[0.9   0.7    0.0854    0.0880];

%Second half
for i=ceil(length(UniqueScenes)/2)+1:length(UniqueScenes)
   
  subplot(2,2,3);
  shadedErrorBar(Time_PSTH_SceneOn,SceneOn_Resp_Mean_Each{i},SceneOn_Resp_Sem_Each{i},'lineprops',{ColorForScene(i,:)},'transparent',1,'patchSaturation',0.5);
  hold on
  %Add time marker
  if i==length(UniqueScenes)
      hold on
 % plot([SceneOnMarker(1), SceneOnMarker(1)],[0,max(cellfun(@(x,y) max(x+y),SceneOn_Resp_Mean_Each,SceneOn_Resp_Sem_Each))],'-k','Color','#707276','LineWidth',1.5);
 % plot([SceneOnMarker(2), SceneOnMarker(2)],[0,max(cellfun(@(x,y) max(x+y),SceneOn_Resp_Mean_Each,SceneOn_Resp_Sem_Each))],'-k','Color','#707276','LineWidth',1.5);
   plot([SceneOnMarker(1), SceneOnMarker(1)],[0,max(cellfun(@(x,y) max(x+y),SceneOn_Resp_Mean_Each,SceneOn_Resp_Sem_Each))],'--k','LineWidth',1.5);
  plot([SceneOnMarker(2), SceneOnMarker(2)],[0,max(cellfun(@(x,y) max(x+y),SceneOn_Resp_Mean_Each,SceneOn_Resp_Sem_Each))],'--k','LineWidth',1.5);
  end
  set(findall(gca, 'Type', 'Line'),'LineWidth',3);
  xlabel('Time from scene onset(ms)');
  ylabel('Firing Rate(Hz)');
  xlim([-PlotSceneOnMarkerStart,PlotSceneOnMarkerStop+SceneOnMarker(2)]);

  title('SceneOnset');
  set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
     
     
  subplot(2,2,4);
   shadedErrorBar(Time_PSTH_SceneOff,SceneOff_Resp_Mean_Each{i},SceneOff_Resp_Sem_Each{i},'lineprops',{ColorForScene(i,:)},'transparent',1,'patchSaturation',0.5);
  %Add time marker
  if i==length(UniqueScenes)
    plot([SceneOffMarker(1), SceneOffMarker(1)],[0,max(cellfun(@(x,y) max(x+y),SceneOff_Resp_Mean_Each,SceneOff_Resp_Sem_Each))],'--k','LineWidth',1.5);
  plot([SceneOffMarker(2), SceneOffMarker(2)],[0,max(cellfun(@(x,y) max(x+y),SceneOff_Resp_Mean_Each,SceneOff_Resp_Sem_Each))],'--k','LineWidth',1.5);
  end
   set(findall(gca, 'Type', 'Line'),'LineWidth',3);
    xlabel('Time from scene offset(ms)');
  ylabel('Firing Rate(Hz)');
  xlim([-PlotSceneOffMarkerStart,PlotSceneOffMarkerStop+SceneOffMarker(2)]);
   title('SceneOffset');
   set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
     hold on
  

end
  ld=legend([string(UniqueScenes(ceil(length(UniqueScenes)/2)+1:end))]);
%ld=legend([string(UniqueScene);'SSI']);
ld.Box='off';
ld.Position=[0.9   0.5    0.0854    0.0880];
%ld=legend([string(UniqueScene);'SSI']);
 %Add data information
  %  data_msg=sprintf(' %s, SpikeChannel: %d, TrialNum: %d ',FileName,SpikeChannel,TrialNum);
    data_msg=strcat(FileName,' SpikeChannel:',num2str(SpikeChannel),';ChannelID:',num2str(ChannelID),';ClusterID:',num2str(ClusterID),';TrialNum: ',num2str(TrialNum),celltype);

    annotation('textbox', [0.17, 0.95, 0.7, 0.04], 'string', data_msg,'FontSize',20,'FontWeight','Bold','LineStyle','none')
 
FigureIndex=FigureIndex+1;
FigureStartNum=FigureStartNum+1;  


%%%%%%%%%%%%%%%%%%%%%Next figure: Scene response average according to the
%%%%%%%%%%%%%%%%%%%%%RT of the first object
figtitlestr{FigureIndex}='ScenePSTH_RT';
fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[150,100, 1000,400],'Name',figtitlestr{FigureIndex});   
subplot(2,2,1);
 shadedErrorBar(Time_PSTH_SceneOn,ResponseOn_HighRT_StableScene_mean,ResponseOn_HighRT_StableScene_sem,'lineprops',{'b'},'transparent',1,'patchSaturation',0.5);
 hold on
 shadedErrorBar(Time_PSTH_SceneOn,ResponseOn_LowRT_StableScene_mean,ResponseOn_LowRT_StableScene_sem,'lineprops',{'r'},'transparent',1,'patchSaturation',0.5);
  set(findall(gca, 'Type', 'Line'),'LineWidth',3);
  
  
    h3=plot([SceneOnMarker(1), SceneOnMarker(1)],[0,max(ResponseOn_HighRT_StableScene_mean+ResponseOn_HighRT_StableScene_sem)],'--k','LineWidth',1.5);
  h4=plot([SceneOnMarker(2), SceneOnMarker(2)],[0,max(ResponseOn_HighRT_StableScene_mean+ResponseOn_HighRT_StableScene_sem)],'--k','LineWidth',1.5);
  
  set( get( get( h3, 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );
 set( get( get( h4, 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );
 
  %Add significance
  sig=area(Time_PSTH_SceneOn,h_Stable_PSTH*0.5,'FaceColor','#FFC107','EdgeColor','#FFC107','FaceAlpha',0.5,'EdgeAlpha',0.5);
  set( get( get( sig, 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );
  
 
    xlim([-PlotSceneOnMarkerStart,PlotSceneOnMarkerStop+SceneOnMarker(2)]);
    
    
  
  
   ylabel('FR(Hz)');
   title('Stable Scene sep by RT');
   set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
  
  
 subplot(2,2,3);
 if ~isempty(ResponseOn_HighRT_FlexibleScene_mean)
 shadedErrorBar(Time_PSTH_SceneOn,ResponseOn_HighRT_FlexibleScene_mean,ResponseOn_HighRT_FlexibleScene_sem,'lineprops',{'b'},'transparent',1,'patchSaturation',0.5);
 hold on
 shadedErrorBar(Time_PSTH_SceneOn,ResponseOn_LowRT_FlexibleScene_mean,ResponseOn_LowRT_FlexibleScene_sem,'lineprops',{'r'},'transparent',1,'patchSaturation',0.5);
 
 set(findall(gca, 'Type', 'Line'),'LineWidth',3);
  
     h3=plot([SceneOnMarker(1), SceneOnMarker(1)],[0,max(ResponseOn_HighRT_StableScene_mean+ResponseOn_HighRT_StableScene_sem)],'--k','LineWidth',1.5);
  h4=plot([SceneOnMarker(2), SceneOnMarker(2)],[0,max(ResponseOn_HighRT_StableScene_mean+ResponseOn_HighRT_StableScene_sem)],'--k','LineWidth',1.5);
  
  set( get( get( h3, 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );
 set( get( get( h4, 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );
 
 sig=area(Time_PSTH_SceneOn,h_Flexible_PSTH*0.5,'FaceColor','#FFC107','EdgeColor','#FFC107','FaceAlpha',0.5,'EdgeAlpha',0.5);
  set( get( get( sig, 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );
 
    xlim([-PlotSceneOnMarkerStart,PlotSceneOnMarkerStop+SceneOnMarker(2)]);
  
    legend('HighRT','Low RT');
 
 
  ylabel('FR(Hz)');
   title('Flexible Scene sep by RT');
   
  set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
 end
     data_msg=strcat(FileName,' SpikeChannel:',num2str(SpikeChannel),';ChannelID:',num2str(ChannelID),';ClusterID:',num2str(ClusterID),';TrialNum: ',num2str(TrialNum),celltype);

    annotation('textbox', [0.1, 0.96, 0.8, 0.04], 'string', data_msg,'FontSize',15,'FontWeight','Bold','LineStyle','none');
    
 %Raw RT bar plot
 subplot(2,2,2);
 b=bar([HighRT_Stable_Mean,LowRT_Stable_Mean]);
 b.FaceColor = 'flat';
b.EdgeColor = 'flat';
ColorChosen= [0,0,1;1,0,0];
b.CData = ColorChosen;

 hold on
 errorbar([1,2],[HighRT_Stable_Mean,LowRT_Stable_Mean],[HighRT_Stable_Sem,LowRT_Stable_Sem],'sk');
 set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
 set(gca,'XTickLabel',{'HighRT','LowRT'});
 box off;
 
 %Add significant level between high RT and low RT
 yrange=ylim;
if h_Stable_RT==1
ctr2 = bsxfun(@plus, b.XData, [b.XOffset]');
hold on
plot(ctr2(1:2), [1 1]*yrange(2)*1.35, '-k', 'LineWidth',3,'Color','r')
if  p_Stable_RT<0.05 & p_Stable_RT>=0.01
   plot(mean(ctr2(1:2)), yrange(2)*1.4, '*k','MarkerSize',15,'Color','r','LineWidth',3);
elseif p_Stable_RT<0.01 & p_Stable_RT>=0.001
   plot(mean(ctr2(1:2))*0.9, yrange(2)*1.4, '*k','MarkerSize',15,'Color','r','LineWidth',3);
   plot(mean(ctr2(1:2))*1.1, yrange(2)*1.4, '*k','MarkerSize',15,'Color','r','LineWidth',3);
elseif p_Stable_RT<0.001 
   plot(mean(ctr2(1:2))*0.9, yrange(2)*1.4, '*k','MarkerSize',15,'Color','r','LineWidth',3);
   plot(mean(ctr2(1:2))*1.1, yrange(2)*1.4, '*k','MarkerSize',15,'Color','r','LineWidth',3);
   plot(mean(ctr2(1:2))*1.3, yrange(2)*1.4, '*k','MarkerSize',15,'Color','r','LineWidth',3);
end

end
 
 
 
 
 subplot(2,2,4);
 b=bar([HighRT_Flexible_Mean,LowRT_Flexible_Mean]);
 b.FaceColor = 'flat';
b.EdgeColor = 'flat';
ColorChosen= [0,0,1;1,0,0];
b.CData = ColorChosen;
box off



 hold on
 errorbar([1,2],[HighRT_Flexible_Mean,LowRT_Flexible_Mean],[HighRT_Flexible_Sem,LowRT_Flexible_Sem],'sk')
 set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
 set(gca,'XTickLabel',{'HighRT','LowRT'});
 
 %Add significant level between high RT and low RT
 yrange=ylim;
if h_Flexible_RT==1
ctr2 = bsxfun(@plus, b.XData, [b.XOffset]');
hold on
plot(ctr2(1:2), [1 1]*yrange(2)*1.35, '-k', 'LineWidth',3,'Color','r')
if  p_Flexible_RT<0.05 &  p_Flexible_RT>=0.01
   plot(mean(ctr2(1:2)), yrange(2)*1.4, '*k','MarkerSize',15,'Color','r','LineWidth',3);
elseif  p_Flexible_RT<0.01 &  p_Flexible_RT>=0.001
   plot(mean(ctr2(1:2))*0.9, yrange(2)*1.4, '*k','MarkerSize',15,'Color','r','LineWidth',3);
   plot(mean(ctr2(1:2))*1.1, yrange(2)*1.4, '*k','MarkerSize',15,'Color','r','LineWidth',3);
elseif  p_Flexible_RT<0.001 
   plot(mean(ctr2(1:2))*0.9, yrange(2)*1.4, '*k','MarkerSize',15,'Color','r','LineWidth',3);
   plot(mean(ctr2(1:2))*1.1, yrange(2)*1.4, '*k','MarkerSize',15,'Color','r','LineWidth',3);
   plot(mean(ctr2(1:2))*1.3, yrange(2)*1.4, '*k','MarkerSize',15,'Color','r','LineWidth',3);
end

end
  
FigureIndex=FigureIndex+1;
FigureStartNum=FigureStartNum+1;  

% FixOffPSTH_currBlock_OrderSep_mean{i,j}=nanmean(PSTH_tmp,1);
%     FixOffPSTH_currBlock_OrderSep_sem{i,j}=nanstd(PSTH_tmp,[],1)./sqrt(sum(~isnan(PSTH_tmp),1));
    
 %%%%%%%Next Figure, Response before fixation offset for different subsequent
    %object orders for each block
    figtitlestr{FigureIndex}='PSTH_DifferentFractalOrder_Blocks';
    fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[130,100, 1000,800],'Name',figtitlestr{FigureIndex});
    plotpointer=1;
     for j=1:NumberOfBlocks
    % subplot(ceil(NumberOfBlocks/4),4,j);
     subplot(NumberOfBlocks,2,(j-1)*2+1)
     currBlock=UniqueBlockType(j);
    BlockOrder(j)=find(AssociateBlock==currBlock);
     %Separated by object order       
        for i=1:(length(ObjOrderSep)+1)
             
           % Eye_V_Explore_currBlock_SacSep{i,j}
            
             shadedErrorBar(PSTH_fixoff_time{j}, FixOffPSTH_currBlock_OrderSep_mean{j,i}, FixOffPSTH_currBlock_OrderSep_sem{j,i},'lineprops',{ColorForOrder{i}},'transparent',1,'patchSaturation',0.5);%PSTH for fix
             set(findall(gca, 'Type', 'Line'),'LineWidth',3);
      
        end
        xlim([-FixOffStartOffset,FixOffStopOffset]);
        
         legend(order_string(j,:));
       
        
       %  ylim([min_Val*0.8,max_Val*1.2+2]);   
            
       % xlim([-PlotStartOffset,RearchingTime_SacSep{i,j}]);
     %  xlim([-PlotStartOffset,RearchingTime_SacSep{1,1}]);
        xlabel('Time from FixOff(ms)');
        ylabel('Firing Rate(Hz)');
       title({Blockstr{BlockOrder(j)}});
         set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
         
            
         subplot(NumberOfBlocks,2,(j)*2)
      %Separated by searching time
       for i=1:SepGroupAmpNum
             
           % Eye_V_Explore_currBlock_SacSep{i,j}
            
             shadedErrorBar(PSTH_fixoff_time{j}, FixOffPSTH_RTGroup_mean{j}(i,:), FixOffPSTH_RTGroup_sem{j}(i,:),'lineprops',{ColorForOrder{i}},'transparent',1,'patchSaturation',0.5);%PSTH for fix
             set(findall(gca, 'Type', 'Line'),'LineWidth',3);
      
       end
       xlim([-FixOffStartOffset,FixOffStopOffset]);
        legend(GroupStr);
         xlabel('Time from FixOff(ms)');
        ylabel('Firing Rate(Hz)');
        title({Blockstr{BlockOrder(j)}});
      set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
            
     end
    data_msg=strcat(FileName,' SpikeChannel:',num2str(SpikeChannel),';ChannelID:',num2str(ChannelID),';ClusterID:',num2str(ClusterID),';TrialNum: ',num2str(TrialNum),celltype);

    annotation('textbox', [0.1, 0.96, 0.8, 0.04], 'string', data_msg,'FontSize',15,'FontWeight','Bold','LineStyle','none');
        
   
FigureStartNum=FigureStartNum+1;
FigureIndex=FigureIndex+1;

%%%%%%%Next Figure, Response before fixation offset for different
%%%%%%%subsequent
figtitlestr{FigureIndex}='PSTH_MostLeastNonResponsiveScene';
fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[150,100, 1000,800],'Name',figtitlestr{FigureIndex});
%Compare the most responsive,remaining responsive, and the nonresponsive
%scene responses
subplot(2,1,1);
%Most responsive scene
shadedErrorBar(Time_PSTH_SceneOn, PSTH_SceneOn_MostResponsiveScenes_Mean, PSTH_SceneOn_MostResponsiveScenes_Sem,'lineprops',{'r'},'transparent',1,'patchSaturation',0.5);%PSTH for fix
set(findall(gca, 'Type', 'Line'),'LineWidth',3);

%Least responsive scene
shadedErrorBar(Time_PSTH_SceneOn, PSTH_SceneOn_RemainingResponsiveScenes_Mean, PSTH_SceneOn_RemainingResponsiveScenes_Sem,'lineprops',{'b'},'transparent',1,'patchSaturation',0.5);%PSTH for fix
set(findall(gca, 'Type', 'Line'),'LineWidth',3);

%Non Responsive Scene
shadedErrorBar(Time_PSTH_SceneOn, PSTH_SceneOn_NonResponsiveScenes_Mean, PSTH_SceneOn_NonResponsiveScenes_Sem,'lineprops',{'k'},'transparent',1,'patchSaturation',0.5);%PSTH for fix
set(findall(gca, 'Type', 'Line'),'LineWidth',3);

legend({'Most Responsive','RemainingResponsive','NonResponsive'});

xlabel('Time from SceneOn(ms)');
ylabel('Firing Rate(Hz)');
title('MostPrefer v.s. Remaining v.s. NonResponsive');

 set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
 
 subplot(2,1,2);
 %Compare the most excitatory responsive, and the inhibitory responsive
%scene responses

%Non Responsive Scene
shadedErrorBar(Time_PSTH_SceneOn, PSTH_SceneOn_NonResponsiveScenes_Mean, PSTH_SceneOn_NonResponsiveScenes_Sem,'lineprops',{'k'},'transparent',1,'patchSaturation',0.5);%PSTH for fix
set(findall(gca, 'Type', 'Line'),'LineWidth',3);

 
 %Excitatory responsive scene
shadedErrorBar(Time_PSTH_SceneOn, PSTH_SceneOn_ExcitatoryResponsiveScene_Mean, PSTH_SceneOn_ExcitatoryResponsiveScene_Sem,'lineprops',{'r'},'transparent',1,'patchSaturation',0.5);%PSTH for fix
set(findall(gca, 'Type', 'Line'),'LineWidth',3);

%Inhibitory responsive scene
shadedErrorBar(Time_PSTH_SceneOn, PSTH_SceneOn_InhibitoryResponsiveScene_Mean, PSTH_SceneOn_InhibitoryResponsiveScene_Sem,'lineprops',{'b'},'transparent',1,'patchSaturation',0.5);%PSTH for fix
set(findall(gca, 'Type', 'Line'),'LineWidth',3);
 
 legend({'NonResponsive','Excitatory Response','Inhibitory Response'});
 xlabel('Time from SceneOn(ms)');
ylabel('Firing Rate(Hz)');
 
title('Excitatory v.s. Inhibitory v.s. NonResponsive');
 
 set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
 
 
data_msg=strcat(FileName,' SpikeChannel:',num2str(SpikeChannel),';ChannelID:',num2str(ChannelID),';ClusterID:',num2str(ClusterID),';TrialNum: ',num2str(TrialNum),celltype);

 annotation('textbox', [0.1, 0.96, 0.8, 0.04], 'string', data_msg,'FontSize',15,'FontWeight','Bold','LineStyle','none');
   
 
       

FigureStartNum=FigureStartNum+1;
FigureIndex=FigureIndex+1;

%% Prefered scene v.s. Other scene responses in each context

figtitlestr{FigureIndex}='PSTHPreferScene in each context';
fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[150,100, 1000,800],'Name',figtitlestr{FigureIndex});

%Stable Scene Context
subplot(2,1,1);
shadedErrorBar(Time_PSTH_SceneOn, PSTH_MostPreferedStableScene_Mean, PSTH_MostPreferedStableScene_Sem,'lineprops',{ColorForEachMode(AssociateBlock==11,:)},'transparent',1,'patchSaturation',0.5);%PSTH for fix
set(findall(gca, 'Type', 'Line'),'LineWidth',3);
hold on

shadedErrorBar(Time_PSTH_SceneOn, PSTH_OtherStableScene_Mean, PSTH_OtherStableScene_Sem,'lineprops',{'k'},'transparent',1,'patchSaturation',0.5);%PSTH for fix
set(findall(gca, 'Type', 'Line'),'LineWidth',3);

legend({'MostPreferedScenes','OtherScenes'});
% xlabel('Time from SceneOn(ms)');
ylabel('Firing Rate(Hz)');

 
title('Stable Scenes');
 set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');


%Flexible Scene Context
subplot(2,1,2);
shadedErrorBar(Time_PSTH_SceneOn, PSTH_MostPreferedFlexibleScene_Mean, PSTH_MostPreferedFlexibleScene_Sem,'lineprops',{ColorForEachMode(AssociateBlock==1,:)},'transparent',1,'patchSaturation',0.5);%PSTH for fix
set(findall(gca, 'Type', 'Line'),'LineWidth',3);
hold on

shadedErrorBar(Time_PSTH_SceneOn, PSTH_OtherFlexibleScene_Mean, PSTH_OtherFlexibleScene_Sem,'lineprops',{'k'},'transparent',1,'patchSaturation',0.5);%PSTH for fix
set(findall(gca, 'Type', 'Line'),'LineWidth',3);

legend({'MostPreferedScenes','OtherScenes'});
 xlabel('Time from SceneOn(ms)');
ylabel('Firing Rate(Hz)');
 
title('Flexible Scenes');
 
 set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
 
 
data_msg=strcat(FileName,' SpikeChannel:',num2str(SpikeChannel),';ChannelID:',num2str(ChannelID),';ClusterID:',num2str(ClusterID),';TrialNum: ',num2str(TrialNum),celltype);

 annotation('textbox', [0.1, 0.96, 0.8, 0.04], 'string', data_msg,'FontSize',15,'FontWeight','Bold','LineStyle','none');
 

FigureStartNum=FigureStartNum+1;
FigureIndex=FigureIndex+1;

%PSTH_MostPreferedFlexibleScene_Mean

%% R square for each condition;
figtitlestr{FigureIndex}='R Square in each context';
fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[170,100, 1000,800],'Name',figtitlestr{FigureIndex});
%Overall
subplot(2,1,1)
plot(Time_PSTH_SceneOn,rsquare_psth,'-k','LineWidth',3);
%Add significance
hold on
sig=area(Time_PSTH_SceneOn,h_sig_ssi_withtime*0.01,'FaceColor','k','EdgeColor','none','FaceAlpha',0.5,'EdgeAlpha',0.5);
 xlim([-PlotSceneOnMarkerStart,PlotSceneOnMarkerStop+SceneOnMarker(2)]);

title('Scene Selectivity Dynamics for all the scenes');
xlabel('Time from scene onset');
ylabel('SSI');
set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
box off

subplot(2,1,2)
plot(Time_PSTH_SceneOn,rsquare_psth_Stable,'-k','Color',ColorForEachMode(AssociateBlock==11,:),'LineWidth',3);
hold on
plot(Time_PSTH_SceneOn,rsquare_psth_Flexible,'-k','Color',ColorForEachMode(AssociateBlock==1,:),'LineWidth',3);
 xlim([-PlotSceneOnMarkerStart,PlotSceneOnMarkerStop+SceneOnMarker(2)]);

title('Scene Selectivity Dynamics for stable scenes and flexible scenes');
xlabel('Time from scene onset');
ylabel('SSI');
set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');

%Add significance
sig=area(Time_PSTH_SceneOn,h_sig_ssi_withtime_Stable*0.01,'FaceColor',ColorForEachMode(AssociateBlock==11,:),'EdgeColor','none','FaceAlpha',0.5,'EdgeAlpha',0.5);
%Add significance
sig=area(Time_PSTH_SceneOn,-h_sig_ssi_withtime_Flexible*0.01,'FaceColor',ColorForEachMode(AssociateBlock==1,:),'EdgeColor','none','FaceAlpha',0.5,'EdgeAlpha',0.5);
box off

data_msg=strcat(FileName,' SpikeChannel:',num2str(SpikeChannel),';ChannelID:',num2str(ChannelID),';ClusterID:',num2str(ClusterID),';TrialNum: ',num2str(TrialNum),celltype);

 annotation('textbox', [0.1, 0.96, 0.8, 0.04], 'string', data_msg,'FontSize',15,'FontWeight','Bold','LineStyle','none');
 FigureStartNum=FigureStartNum+1;
FigureIndex=FigureIndex+1;

%{
%%Average scene responses for two types of scenes

figtitlestr{FigureIndex}='Average Scene Responses';
fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[190,100, 1000,800],'Name',figtitlestr{FigureIndex});

shadedErrorBar(Time_PSTH_SceneOn,MeanPSTH_FR_SceneOn_Stable_PSTH,MeanPSTH_FR_SceneOn_Stable_PSTH_SEM,'lineprops',{ColorForEachMode(AssociateBlock==11,:)},'transparent',1,'patchSaturation',0.5);
 hold on
 shadedErrorBar(Time_PSTH_SceneOn,MeanPSTH_FR_SceneOn_Flexible_PSTH,MeanPSTH_FR_SceneOn_Flexible_PSTH_SEM,'lineprops',{ColorForEachMode(AssociateBlock==1,:)},'transparent',1,'patchSaturation',0.5);
  set(findall(gca, 'Type', 'Line'),'LineWidth',3);
  
  
    h3=plot([SceneOnMarker(1), SceneOnMarker(1)],[0,max(max(MeanPSTH_FR_SceneOn_Stable_PSTH+MeanPSTH_FR_SceneOn_Stable_PSTH_SEM))],'--k','LineWidth',1.5);
  h4=plot([SceneOnMarker(2), SceneOnMarker(2)],[0,max(max(MeanPSTH_FR_SceneOn_Stable_PSTH+MeanPSTH_FR_SceneOn_Stable_PSTH_SEM))],'--k','LineWidth',1.5);
  
  set( get( get( h3, 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );
 set( get( get( h4, 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );
 
  %Add significance
 % sig=area(Time_PSTH_SceneOn,h_Stable_PSTH*0.5,'FaceColor','#FFC107','EdgeColor','#FFC107','FaceAlpha',0.5,'EdgeAlpha',0.5);
 % set( get( get( sig, 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );
  
 
    xlim([-PlotSceneOnMarkerStart,PlotSceneOnMarkerStop+SceneOnMarker(2)]);
   title('Scene response for stable scenes and flexible scenes');
xlabel('Time from scene onset');
ylabel('Response ');
set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');

   
     
FigureStartNum=FigureStartNum+1;
FigureIndex=FigureIndex+1;
%}
%% R^2 PSTH FOR before fix off period, for all orders
%rsquare_psth_fixoff_order{i,j}=rsquare_psth_fixoff;
%pFixOffWithTime_order{i,j}=pFixOffWithTime;
%h_sig_ssi_fixoffwithtime_order{i,j}=h_sig_ssi_fixoffwithtime;

 figtitlestr{FigureIndex}='SSI_DifferentFractalOrder_FixOffBlocks';
 fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[190,100, 1000,800],'Name',figtitlestr{FigureIndex});
    
     for j=1:NumberOfBlocks
    % subplot(ceil(NumberOfBlocks/4),4,j);
     subplot(ceil(NumberOfBlocks/2),2,j)
     currBlock=UniqueBlockType(j);
    BlockOrder(j)=find(AssociateBlock==currBlock);
     %Separated by object order       
        for i=1:(length(ObjOrderSep)+1)
             
           % Eye_V_Explore_currBlock_SacSep{i,j}
            plot(PSTH_fixoff_time{j},rsquare_psth_fixoff_order{j,i},'-k','Color',ColorForOrder{i});
            % shadedErrorBar(PSTH_fixoff_time{j}, FixOffPSTH_currBlock_OrderSep_mean{j,i}, FixOffPSTH_currBlock_OrderSep_sem{j,i},'lineprops',{ColorForOrder{i}},'transparent',1,'patchSaturation',0.5);%PSTH for fix
             set(findall(gca, 'Type', 'Line'),'LineWidth',3);
             hold on
      
        end
        xlim([-FixOffStartOffset,FixOffStopOffset]);
        
         legend(order_string(j,:));
       
        
       %  ylim([min_Val*0.8,max_Val*1.2+2]);   
            
       % xlim([-PlotStartOffset,RearchingTime_SacSep{i,j}]);
     %  xlim([-PlotStartOffset,RearchingTime_SacSep{1,1}]);
        xlabel('Time from FixOff(ms)');
        ylabel('SSI');
       title({Blockstr{BlockOrder(j)}});
         set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
         
           ylim([0,1]); 
 
            
     end
    data_msg=strcat(FileName,' SpikeChannel:',num2str(SpikeChannel),';ChannelID:',num2str(ChannelID),';ClusterID:',num2str(ClusterID),';TrialNum: ',num2str(TrialNum),celltype);

    annotation('textbox', [0.1, 0.96, 0.8, 0.04], 'string', data_msg,'FontSize',15,'FontWeight','Bold','LineStyle','none');
        
   
FigureStartNum=FigureStartNum+1;
FigureIndex=FigureIndex+1;

%%
figtitlestr{FigureIndex}='ResponsesFractalOrder_FixOffEachScene';
fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[200,100, 1000,800],'Name',figtitlestr{FigureIndex});
%First half
for i=1:ceil(length(UniqueScenes)/2)
    
  for j=1:3 %Fractal order loop 
  subplot(2,3,j);
  shadedErrorBar(PSTH_fixoff_time{1},FixOffPSTH_Scene_OrderSep_mean{i,j},FixOffPSTH_Scene_OrderSep_sem{i,j},'lineprops',{ColorForScene(i,:)},'transparent',1,'patchSaturation',0.5);
  hold on
  %Add time marker
 % if i==ceil(length(UniqueScenes)/2)
  %    hold on
 % plot([SceneOnMarker(1), SceneOnMarker(1)],[0,max(cellfun(@(x,y) max(x+y),SceneOn_Resp_Mean_Each,SceneOn_Resp_Sem_Each))],'-k','Color','#707276','LineWidth',1.5);
 % plot([SceneOnMarker(2), SceneOnMarker(2)],[0,max(cellfun(@(x,y) max(x+y),SceneOn_Resp_Mean_Each,SceneOn_Resp_Sem_Each))],'-k','Color','#707276','LineWidth',1.5);
 %  h1=plot([SceneOnMarker(1), SceneOnMarker(1)],[0,max(cellfun(@(x,y) max(x+y),SceneOn_Resp_Mean_Each,SceneOn_Resp_Sem_Each))],'--k','LineWidth',1.5);
 % h2=plot([SceneOnMarker(2), SceneOnMarker(2)],[0,max(cellfun(@(x,y) max(x+y),SceneOn_Resp_Mean_Each,SceneOn_Resp_Sem_Each))],'--k','LineWidth',1.5);
  
%  set( get( get( h1, 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );
% set( get( get( h2, 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );
 
  %end
  
  

  set(findall(gca, 'Type', 'Line'),'LineWidth',3);
  xlabel('Time from fixation off(ms)');
  ylabel('Firing Rate(Hz)');
  

  title(sprintf('Fix off Order %d',j));
  set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
   %xlim([-PlotSceneOnMarkerStart,PlotSceneOnMarkerStop+SceneOnMarker(2)]);
   xlim([-FixOffStartOffset,FixOffStopOffset]);
   %{
  subplot(2,2,2);
   shadedErrorBar(Time_PSTH_SceneOff,SceneOff_Resp_Mean_Each{i},SceneOff_Resp_Sem_Each{i},'lineprops',{ColorForScene(i,:)},'transparent',1,'patchSaturation',0.5);
  %Add time marker
  if i==ceil(length(UniqueScenes)/2)
    h3=plot([SceneOffMarker(1), SceneOffMarker(1)],[0,max(cellfun(@(x,y) max(x+y),SceneOff_Resp_Mean_Each,SceneOff_Resp_Sem_Each))],'--k','LineWidth',1.5);
  h4=plot([SceneOffMarker(2), SceneOffMarker(2)],[0,max(cellfun(@(x,y) max(x+y),SceneOff_Resp_Mean_Each,SceneOff_Resp_Sem_Each))],'--k','LineWidth',1.5);
  
 set( get( get( h3, 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );
 set( get( get( h4, 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );
  end
  
  
    xlim([-PlotSceneOffMarkerStart,PlotSceneOffMarkerStop+SceneOffMarker(2)]);
    
   set(findall(gca, 'Type', 'Line'),'LineWidth',3);
    xlabel('Time from scene offset(ms)');
  ylabel('Firing Rate(Hz)');

   title('SceneOffset');
   set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
     hold on
   %}
  end %End of the fractal order loop

    
end
ld=legend([string(UniqueScenes(1:ceil(length(UniqueScenes)/2)))]);
%ld=legend([string(UniqueScene);'SSI']);
ld.Box='off';
ld.Position=[0.9   0.7    0.0854    0.0880];

%Second half
for i=ceil(length(UniqueScenes)/2)+1:length(UniqueScenes)
  for j=1:3 %Fractal order loop  
  subplot(2,3,j+3);
%  shadedErrorBar(Time_PSTH_SceneOn,SceneOn_Resp_Mean_Each{i},SceneOn_Resp_Sem_Each{i},'lineprops',{ColorForScene(i,:)},'transparent',1,'patchSaturation',0.5);
shadedErrorBar(PSTH_fixoff_time{1},FixOffPSTH_Scene_OrderSep_mean{i,j},FixOffPSTH_Scene_OrderSep_sem{i,j},'lineprops',{ColorForScene(i,:)},'transparent',1,'patchSaturation',0.5);
  
hold on
%{
  %Add time marker
  if i==length(UniqueScenes)
      hold on
 % plot([SceneOnMarker(1), SceneOnMarker(1)],[0,max(cellfun(@(x,y) max(x+y),SceneOn_Resp_Mean_Each,SceneOn_Resp_Sem_Each))],'-k','Color','#707276','LineWidth',1.5);
 % plot([SceneOnMarker(2), SceneOnMarker(2)],[0,max(cellfun(@(x,y) max(x+y),SceneOn_Resp_Mean_Each,SceneOn_Resp_Sem_Each))],'-k','Color','#707276','LineWidth',1.5);
   plot([SceneOnMarker(1), SceneOnMarker(1)],[0,max(cellfun(@(x,y) max(x+y),SceneOn_Resp_Mean_Each,SceneOn_Resp_Sem_Each))],'--k','LineWidth',1.5);
  plot([SceneOnMarker(2), SceneOnMarker(2)],[0,max(cellfun(@(x,y) max(x+y),SceneOn_Resp_Mean_Each,SceneOn_Resp_Sem_Each))],'--k','LineWidth',1.5);
  end
%}
  set(findall(gca, 'Type', 'Line'),'LineWidth',3);
 % xlabel('Time from scene onset(ms)');
  xlabel('Time from fixation off(ms)');
  ylabel('Firing Rate(Hz)');
 % xlim([-PlotSceneOnMarkerStart,PlotSceneOnMarkerStop+SceneOnMarker(2)]);
 xlim([-FixOffStartOffset,FixOffStopOffset]);
 % title('SceneOnset');
  title(sprintf('Fix off Order %d',j));
  set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
     
 %{    
  subplot(2,2,4);
   shadedErrorBar(Time_PSTH_SceneOff,SceneOff_Resp_Mean_Each{i},SceneOff_Resp_Sem_Each{i},'lineprops',{ColorForScene(i,:)},'transparent',1,'patchSaturation',0.5);
  %Add time marker
  if i==length(UniqueScenes)
    plot([SceneOffMarker(1), SceneOffMarker(1)],[0,max(cellfun(@(x,y) max(x+y),SceneOff_Resp_Mean_Each,SceneOff_Resp_Sem_Each))],'--k','LineWidth',1.5);
  plot([SceneOffMarker(2), SceneOffMarker(2)],[0,max(cellfun(@(x,y) max(x+y),SceneOff_Resp_Mean_Each,SceneOff_Resp_Sem_Each))],'--k','LineWidth',1.5);
  end
   set(findall(gca, 'Type', 'Line'),'LineWidth',3);
    xlabel('Time from scene offset(ms)');
  ylabel('Firing Rate(Hz)');
  xlim([-PlotSceneOffMarkerStart,PlotSceneOffMarkerStop+SceneOffMarker(2)]);
   title('SceneOffset');
   set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
     hold on
  %}
  end

end
  ld=legend([string(UniqueScenes(ceil(length(UniqueScenes)/2)+1:end))]);
%ld=legend([string(UniqueScene);'SSI']);
ld.Box='off';
ld.Position=[0.9   0.5    0.0854    0.0880];
%ld=legend([string(UniqueScene);'SSI']);
 %Add data information
  %  data_msg=sprintf(' %s, SpikeChannel: %d, TrialNum: %d ',FileName,SpikeChannel,TrialNum);
data_msg=strcat(FileName,' SpikeChannel:',num2str(SpikeChannel),';ChannelID:',num2str(ChannelID),';ClusterID:',num2str(ClusterID),';TrialNum: ',num2str(TrialNum),celltype);

annotation('textbox', [0.17, 0.95, 0.7, 0.04], 'string', data_msg,'FontSize',20,'FontWeight','Bold','LineStyle','none')


 
FigureStartNum=FigureStartNum+1;
FigureIndex=FigureIndex+1;

%FixOffPSTH_Scene_OrderSep_mean{i,j}=nanmean(tmp,1);
%          FixOffPSTH_Scene_OrderSep_sem{i,j}=nanstd(tmp,[],1)./sqrt(sum(~isnan(tmp),1));


%%
%{
%Average PSTH
MeanPSTH_FR_SceneOn_Flexible_PSTH=nanmean(ResponseOn(BlockType==1,:),1);
MeanPSTH_FR_SceneOn_Flexible_PSTH_SEM=nanstd(ResponseOn(BlockType==1,:),[],1)./sqrt(sum(~isnan(ResponseOn(BlockType==1,:)),1));

MeanPSTH_FR_SceneOn_Stable_PSTH=nanmean(ResponseOn(BlockType==11,:),1);
MeanPSTH_FR_SceneOn_Stable_PSTH_SEM=nanstd(ResponseOn(BlockType==11,:),[],1)./sqrt(sum(~isnan(ResponseOn(BlockType==11,:)),1));
%}

%Below is combined to the first figures
%{
%%%%%%%%%%%%%%%%%%Next figure: Average Spon Change with block type 

 figtitlestr{FigureIndex}='SponChangeWithTrials';
 fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[140,100, 800,400],'Name',figtitlestr{FigureIndex});   
%plot(SponResponse)
PlotIndex=1:length(BeforeSponResponse);
for i=1:NumberOfBlocks
    currBlock=UniqueBlockType(i);
    Sel=BlockType==currBlock;
    plot(PlotIndex(Sel),BeforeSponResponse(Sel),'-k','Color',ColorForEachMode(AssociateBlock==currBlock,:),'LineWidth',3);
  
    
    hold on
   
end
box off;
set(gca,'XTick',PlotIndex(1:50:end));
set(gca,'XTickLabel',TrialIndex(1:50:end));
xlabel('Trial Index');
ylabel('Spontaneous Firing Rate(Hz)');
set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');

%}


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

OutputData_New.EventSegRT.TrialType='GoodTrials';

%%%%%%%%%Save overall scene response characteristics%%%%%%%%%%%%%%%%%%%%%%%
%Organize the data
NameStr={'AverageFR','AverageSpon','AverageSceneDiffBlock','p_SceneDiffBetweenBlocks',...
    'PSTH_SceneOn_FleSta_TimeMeanSem',...
    'PSTH_SceneOff_FleSta_TimeMeanSem',...
    'SI_FleSta','p_AnovaSceneOnAllFleSta',...
  'SceneOn500MeanMinusSpon','p_SceneOn500withSpon','SceneOnMeanMinusSponWhole','p_SceneOnWithSponWhole',...
  'pSceneOnWhole','SceneOnDifference_WholeMean','pSceneOn500','SceneOnDifference_500Mean',...
  'SD_RelSceneOn','Mean_RelSceneOn','UniqueScene'...
  'NumOfStableSigScenes','NumOfFlexibleSigScenes','StableSigScenes','FlexibleSigScenes'...
  'PSTH_MostPreferedStableScene_MeanSem','PSTH_MostPreferedFlexibleScene_MeanSem',...
  'PSTH_OtherStableScene_MeanSem','PSTH_OtherFlexibleScene_MeanSem',...
  'R2WholeFleSta','p_SceneOnAnova_WholeStaFle','pSceneDiffWithTime','rsquare_psth',...
  'pSceneDiffWithTime_Stable','rsquare_psth_Stable','pSceneDiffWithTime_Flexible','rsquare_psth_Flexible',...
  'PropOfSig_Stable','PropOfSig_Flexible',...
  'p_SceneOnAnova_StableHalf','p_SceneOnAnova_FlexibleHalf',...
  'R2_StableHalf','R2_FlexibleHalf',...
    
 
};




DataLib={AverageFR,AverageSpon,SceneOnAverage',p_SceneDiffBetweenBlocks,...
    [Time_PSTH_SceneOn',MeanPSTH_FR_SceneOn_Flexible',MeanPSTH_FR_SceneOn_Stable',SemPSTH_FR_SceneOn_Flexible',SemPSTH_FR_SceneOn_Stable'],...
     [Time_PSTH_SceneOff',MeanPSTH_FR_SceneOff_Flexible',MeanPSTH_FR_SceneOff_Stable',SemPSTH_FR_SceneOff_Flexible',SemPSTH_FR_SceneOff_Stable'],...
   [SI_Stable,SI_Flexible],[p_SceneOnAnova,p_SceneOnAnova_Stable,p_SceneOnAnova_Flexible],...
SceneOn500MeanMinusSpon,p_SceneOn500WithSpon,SceneOnMeanMinusSponWhole,p_SceneOnWithSponWhole,...
pSceneOnWhole,SceneOnDifference_WholeMean,pSceneOn500,SceneOnDifference_500Mean,...
SD_RelSceneOn,Mean_RelSceneOn,UniqueScenes...
NumOfStableSigScenes,NumOfFlexibleSigScenes,StableSigScenes,FlexibleSigScenes,...
[PSTH_MostPreferedStableScene_Mean',PSTH_MostPreferedStableScene_Sem'],[PSTH_MostPreferedFlexibleScene_Mean',PSTH_MostPreferedFlexibleScene_Sem'],...
[PSTH_OtherStableScene_Mean',PSTH_OtherStableScene_Sem'],[PSTH_OtherFlexibleScene_Mean',PSTH_OtherFlexibleScene_Sem'],...
[R2_Whole,R2_Flexible,R2_Stable],[p_SceneOnAnova,p_SceneOnAnova_Stable,p_SceneOnAnova_Flexible],pSceneDiffWithTime,rsquare_psth,...
pSceneDiffWithTime_Stable,rsquare_psth_Stable,pSceneDiffWithTime_Flexible,rsquare_psth_Flexible,...
 PropOfSig_Stable,PropOfSig_Flexible,...
 [p_SceneOnAnova_StableHalf1,p_SceneOnAnova_StableHalf2],[p_SceneOnAnova_FlexibleHalf1,p_SceneOnAnova_FlexibleHalf2],...
 [R2_StableHalf1,R2_StableHalf2],[R2_FlexibleHalf1,R2_FlexibleHalf2]
};

%{
[p_SceneOnAnova_StableHalf1,t,stats]=anova1(SceneOnResponseStable_Raw,SceneName_Stable_Raw,'off');
%R square
ss=t{2,2};
ts=t{4,2};   
R2_StableHalf1=ss/ts;
%}


DataStamp=containers.Map(NameStr,DataLib);

OutputData_New.EventSegRT.Task=TaskType;
OutputData_New.EventSegRT.TaskCode=TaskCode;
OutputData_New.EventSegRT.DataStamp=DataStamp;
%Store figures in common folders
ProtocolMarker='EventSegRT';
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
Task_Old=OutputData.EventSegRT.TaskCode;

     if Task_Old~=TaskCode
        %Not the same task, add the new dataset into the old one
        OutputData.EventSegRT(numel(OutputData.EventSegRT)+1)=OutputData_New.EventSegRT;
     else
         %If the same task,the same protocol,replace the old one with the new one
         
         OutputData.EventSegRT=OutputData_New.EventSegRT;
         
    
    
     end
else
    %Output the current file
   OutputData.EventSegRT=OutputData_New.EventSegRT; 
end
if OutputFlag
    cd(OutPath);
    save(OutputFileName,'OutputData');
    disp('Data Exported');
else
    disp('No data export');
end  



end %End of the loop of the spike channel

end %End of the whole 
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
