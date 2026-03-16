function ITI_CorrAnalysis(Data, StartTrial, EndTrial,ShowFigureFlag,OutputFlag);
%This protocol is to extract the whole ITI period and analysis whether it's
%correlated with the previous scene or predictive of the next scene
%added by Xuefei Yu 03202021
%This procotol is aim for analying of the whole data (including the
%breaktrials);
ECode;
FileName=Data.FileName;
TaskType=Data.TaskType;
TaskCode=Data.TaskCode;
%TrialNum=EndTrial-StartTrial+1;
Selection=StartTrial:EndTrial;%Select the trials according to the pannel

%User Control
PreStartPlotTime=1000;



%Abnormal check %Especially For Robin081120FORAGEXPLO4.rst
CheckCode=SCENE_ON;%This code should appear once in each trial
AbnormalTrials=AbnormalTrialCheck(Data,CheckCode);
Selection=Selection(~ismember(Selection,AbnormalTrials));


CheckCode=SCENE_OFF;%This code should appear once in each trial
AbnormalTrials=AbnormalTrialCheck(Data,CheckCode);
Selection=Selection(~ismember(Selection,AbnormalTrials));

if strcmp(FileName,'Robin092920FORAGEXPLO1.rst') ||strcmp(FileName,'Robin092820FORAGEXPLO1.rst')
    %Remove bad trial for this file

CheckCodePair=[FIXOFF,SHOWFIXCD];%This code should be the same as target on for good trials
AbnormalTrials=AbnormalTrialCheck(Data,CheckCodePair);
Selection=Selection(~ismember(Selection,AbnormalTrials));

if strcmp(FileName,'Robin092820FORAGEXPLO1.rst')
    CheckCodePair=[TGTCD,TGTOFF];%This code should be the same as target on for good trials
    AbnormalTrials=AbnormalTrialCheck(Data,CheckCodePair);
    Selection=Selection(~ismember(Selection,AbnormalTrials));
    
end


%TGTCD
end

TrialIndexRaw=Data.BlockType(:,1);



TrialIndex=TrialIndexRaw(Selection);

TrialNum=length(Selection);
%User defined parameter for the neural/behavioral analysis
TargetOnOffset=50;%in ms, the time selection bin after target onset; 
StartDelay=0;%in ms, the time selection bin after fixation offset;

%For segmentation of  the neural responses 
%SegEvent=[EXTRAFIXSTARTCD,EXTRAFIXOFF,SCENE_ON,SCENE_OFF,SHOWFIXCD,FIXOFF,TGTCD];%Event Flow of interest
%SegEventEnd=[EXTRAFIXOFF,SHOWFIXCD,SHOWFIXCD,NaN,FIXOFF,TGTCD,TGTOFF];%Event end; %Automatic screening later
SegEvent=[SCENE_ON,SCENE_OFF,];%Event Flow of interest
SegEventEnd=[SHOWFIXCD,NaN];%Event end; %Automatic screening later
%SegEventEnd=[EXTRAFIXOFF,FIXOFF,FIXOFF,TGTCD,TGTOFF];%Event end
SegEventStr={'SceneOn','SceneOff'};

PlotStartOffset=[NaN,500];%in ms
PlotStopOffset=[500,NaN];%in ms

EventTypeFlag=[1,1];%1:Number same with the trial number, 2:Number same with the target number


RasterOn=0;%Show raster, but it will be very slow to do so

%Parameter to define the sliding window to measure firing rates
StepSize=10;%in ms
BinWidth=15;%in ms

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
 UniqueBlockType=sort(UniqueBlockType,'ascend');%Flexible first, then stable;
NumberOfBlocks=max(reprojectorderindex);
 
 %If the target onset time and fixofftime dimension are not aligned 
 if size(TargetOnTime,2)~=size(FixOffTime,2)
     FixOffTime=FixOffTime(:,1:size(TargetOnTime,2));
     disp('Warning!FixOffTime Corrected!!')
 end
     
 
SearchingTime= (TargetOnTime-FixOffTime)/1000;%Change into s for presentation

RandFixFlag=UniqueBlockType;

if length(RandFixFlag)==1
    strFile=sprintf('No random/fix comparison data in the file %s\n',FileName);
    disp(strFile);
    disp('You may need to change to Protocol Scene Tuning');
    return  
end
FixSeqTrials=(BlockType==0)|(BlockType==10);
RandSeqTrials=(BlockType==1)|(BlockType==11);


FixSeqSceneName=SceneName(FixSeqTrials);
RandSeqSceneName=SceneName(RandSeqTrials);

%Separate scene types 
FixSceneFlag=ismember(FixSeqSceneName,RandSeqSceneName);
RandSceneFlag=ismember(RandSeqSceneName,FixSeqSceneName);
%CommonSceneFlag=FixSceneFlag & RandSceneFlag;



CommonScenes=FixSeqSceneName(FixSceneFlag);
UniqueSceneCommon=unique(CommonScenes);
%{
if isempty(UniqueSceneCommon)
    disp('No scene for comparison!');
    return
end
%}
%Load the spike channel
SpikeChannelWhole=Data.SpikeChannel;
SpikeTimeWhole=Data.SpikeTimeData(:,:,Selection);
ChannelIDWhole=Data.ChannelID;
ClusterIDWhole=Data.ClusterID;
PreSpikeTimeWhole=Data.SpikeTimeDataPre(Data.SpikeChannel,:,Selection);



MultiChannel=0;
if length(SpikeChannelWhole)>1
    MultiChannel=1;
    
end


for spk=1:length(SpikeChannelWhole)
    SpikeChannel=SpikeChannelWhole(spk);
    SpikeTime_Ori=squeeze(SpikeTimeWhole(spk,:,:))';
    PreSpikeTime=squeeze(PreSpikeTimeWhole(spk,:,:))';
    ChannelID=ChannelIDWhole(spk);
    ClusterID=ClusterIDWhole(spk);

%Combine PreSpikeTime with SpikeTime

  SpikeTimeBody=[PreSpikeTime(:,1:end-2),SpikeTime_Ori(:,1:end-2)];  
  NewStartTime=PreSpikeTime(:,end-1);
  %Except for the first trial
  NewStartTime(1)=SpikeTime_Ori(1,end-1);
  NewEndTime=SpikeTime_Ori(:,end);
  SpikeTime=[SpikeTimeBody,NewStartTime,NewEndTime];
    
%Segregate spike data
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

%Pause here temporally 
%Relative time bin for each event
RelativeTimeMarker=SegmentData.RelativeTimeMarker;%=EventTimeMarker;
%Absolute time marker for each event
AbosoluteTimeMarker=SegmentData.AbsoluteTimeMarker;



%Loop over the scenes for each kind of block:random/fix
ResponseOn=MeanPSTH_FR{ismember(SegEventStr,'SceneOn')};
ResponseOnTime=MeanPSTH_Time{ismember(SegEventStr,'SceneOn')};

SceneOnMarker=nanmean(RelativeTimeMarker{ismember(SegEventStr,'SceneOn')},1);
SceneOnMarkerAbsAll=AbosoluteTimeMarker{ismember(SegEventStr,'SceneOn')};


%Align the PSTH relative to scene onset(0)'
[AlignedPSTH_SceneOn,AlignedTime_SceneOn]= AlignPSTH(ResponseOn,ResponseOnTime,StepSize,0);

%Sceneoff
ResponseOff=MeanPSTH_FR{ismember(SegEventStr,'SceneOff')};
ResponseOffTime=MeanPSTH_Time{ismember(SegEventStr,'SceneOff')};
%Remove unnessary cols
Col_Valid=sum(~isnan(ResponseOff),1);
RemoveNaNCol=Col_Valid<floor(sqrt(size(ResponseOff,1)));
ResponseOff(:,RemoveNaNCol)=[];
ResponseOffTime(:,RemoveNaNCol)=[];


SceneOffMarker=nanmean(RelativeTimeMarker{ismember(SegEventStr,'SceneOff')},1);


%First check the intertrial interval before scene onset for all scenes 
MeanPSTHSceneOn=nanmean(AlignedPSTH_SceneOn,1);
SemPSTHSceneOn=nanstd(AlignedPSTH_SceneOn,[],1)./sqrt(sum(~isnan(AlignedPSTH_SceneOn),1));

%Few number correction
NumColtmp=sum(~isnan(AlignedPSTH_SceneOn),1);
RemoveCol=NumColtmp<3;
MeanPSTHSceneOn(RemoveCol)=NaN;
SemPSTHSceneOn(RemoveCol)=NaN;


MeanPSTHTime=nanmean(AlignedTime_SceneOn,1);

%Some markers
ExtraFPOnTimeAbs=FindOutTime(EventChannel,EventBin,EXTRAFP);
ExtraFPOffTimeAbs=FindOutTime(EventChannel,EventBin,EXTRAFIXOFF);
PreFPTimeAbs=FindOutTime(EventChannel,EventBin,PREFIX);

ExtraFPOnTimeRel=ExtraFPOnTimeAbs-SceneOnMarkerAbsAll(:,1);
ExtraFPOnTimeRelMean=nanmean(ExtraFPOnTimeRel);

ExtraFPOffTimeRel=ExtraFPOffTimeAbs-SceneOnMarkerAbsAll(:,1);
ExtraFPOffTimeRelMean=nanmean(ExtraFPOffTimeRel);

PreFPTimeRel=PreFPTimeAbs-SceneOnMarkerAbsAll(:,1);
PreFPTimeRelMean=nanmean(PreFPTimeRel);

%Then check the intertrial interval after scene offset for all scenes 

MeanPSTHSceneOff=nanmean(ResponseOff,1);
SemPSTHSceneOff=nanstd(ResponseOff,[],1)./sqrt(sum(~isnan(ResponseOff)));
%Few number correction
NumColtmp=sum(~isnan(ResponseOff),1);
RemoveCol=NumColtmp<3;
MeanPSTHSceneOff(RemoveCol)=NaN;
SemPSTHSceneOff(RemoveCol)=NaN;


MeanPSTHTimeOff=nanmean(ResponseOffTime,1);

%Loop over the block types
UniqueBlockType=unique(BlockType);
for i=1:length(UniqueBlockType)
    currBlock=BlockType==UniqueBlockType(i);
    
    ResponseOnSel=AlignedPSTH_SceneOn(currBlock,:);
    
    MeanSceneOnResponseBlock(i,:)=nanmean(ResponseOnSel,1);
    SemSceneOnResponseBlock(i,:)=nanstd(ResponseOnSel,[],1)./sqrt(sum(~isnan(ResponseOnSel)));
    
     ResponseOffSel=ResponseOff(currBlock,:);
    
    MeanSceneOffResponseBlock(i,:)=nanmean(ResponseOffSel,1);
    SemSceneOffResponseBlock(i,:)=nanstd(ResponseOffSel,[],1)./sqrt(sum(~isnan(ResponseOffSel)));
    
    
end




UniqueScene=unique(SceneName);

%
%Loop over the scenes and block types

for i=1:length(UniqueScene)
    Sel=(SceneName==UniqueScene(i));
    Sel_Fix=((BlockType==0)|(BlockType==10) )& Sel;
    Sel_Rand=((BlockType==1)|(BlockType==11)) & Sel;
    %Fix Block
    currPSTH_SceneOn_Fix=AlignedPSTH_SceneOn(Sel_Fix,:);
    currPSTHTime_SceneOn_Fix=AlignedTime_SceneOn(Sel_Fix,:);
    
    currPSTH_SceneOn_Fix_Mean(i,:)=nanmean(currPSTH_SceneOn_Fix,1);
    currPSTH_SceneOn_Fix_Sem(i,:)=nanstd(currPSTH_SceneOn_Fix,[],1)./sqrt(sum(~isnan(currPSTH_SceneOn_Fix),1));
    currPSTHTime_SceneOn_Fix_Mean(i,:)=nanmean(currPSTHTime_SceneOn_Fix,1);
    
    currPSTH_SceneOn_Fix_tmp{i}=currPSTH_SceneOn_Fix;
    %
    %Remove points of rep<3
    ValidCol=sum(~isnan(currPSTH_SceneOn_Fix),1);
    
    currPSTH_SceneOn_Fix_Mean(i,ValidCol<3)=NaN;
    currPSTH_SceneOn_Fix_Sem(i,ValidCol<3)=NaN;
    currPSTHTime_SceneOn_Fix_Mean(i,ValidCol<3)=NaN;
    
    %
    %Rand Block
     currPSTH_SceneOn_Rand=AlignedPSTH_SceneOn(Sel_Rand,:);
    currPSTHTime_SceneOn_Rand=AlignedTime_SceneOn(Sel_Rand,:);
    
    currPSTH_SceneOn_Rand_Mean(i,:)=nanmean(currPSTH_SceneOn_Rand,1);
    currPSTH_SceneOn_Rand_Sem(i,:)=nanstd(currPSTH_SceneOn_Rand,[],1)./sqrt(sum(~isnan(currPSTH_SceneOn_Rand),1));
    currPSTHTime_SceneOn_Rand_Mean(i,:)=nanmean(currPSTHTime_SceneOn_Rand,1);
   
    %Remove points of rep<3
    ValidCol=sum(~isnan(currPSTH_SceneOn_Rand),1);
    
    currPSTH_SceneOn_Rand_Mean(i,ValidCol<3)=NaN;
    currPSTH_SceneOn_Rand_Sem(i,ValidCol<3)=NaN;
    currPSTHTime_SceneOn_Rand_Mean(i,ValidCol<3)=NaN;
    
    %Next, response off
    
    
    %Fix Block
    currPSTH_SceneOff_Fix=ResponseOff(Sel_Fix,:);
    currPSTHTime_SceneOff_Fix=ResponseOffTime(Sel_Fix,:);
    
    currPSTH_SceneOff_Fix_Mean(i,:)=nanmean(currPSTH_SceneOff_Fix,1);
    currPSTH_SceneOff_Fix_Sem(i,:)=nanstd(currPSTH_SceneOff_Fix,[],1)./sqrt(sum(~isnan(currPSTH_SceneOff_Fix),1));
    currPSTHTime_SceneOff_Fix_Mean(i,:)=nanmean(currPSTHTime_SceneOff_Fix,1);
    
    %Remove points of rep<3
    ValidCol=sum(~isnan(currPSTH_SceneOff_Fix),1);
    
    currPSTH_SceneOff_Fix_Mean(i,ValidCol<3)=NaN;
    currPSTH_SceneOff_Fix_Sem(i,ValidCol<3)=NaN;
    currPSTHTime_SceneOff_Fix_Mean(i,ValidCol<3)=NaN;
    
    
    %Rand Block
     currPSTH_SceneOff_Rand=ResponseOff(Sel_Rand,:);
    currPSTHTime_SceneOff_Rand=ResponseOffTime(Sel_Rand,:);
    
    currPSTH_SceneOff_Rand_Mean(i,:)=nanmean(currPSTH_SceneOff_Rand,1);
    currPSTH_SceneOff_Rand_Sem(i,:)=nanstd(currPSTH_SceneOff_Rand,[],1)./sqrt(sum(~isnan(currPSTH_SceneOff_Rand),1));
    currPSTHTime_SceneOff_Rand_Mean(i,:)=nanmean(currPSTHTime_SceneOff_Rand,1);
   
    %Remove points of rep<3
    ValidCol=sum(~isnan(currPSTH_SceneOff_Rand),1);
    
    currPSTH_SceneOff_Rand_Mean(i,ValidCol<3)=NaN;
    currPSTH_SceneOff_Rand_Sem(i,ValidCol<3)=NaN;
    currPSTHTime_SceneOff_Rand_Mean(i,ValidCol<3)=NaN;
    
    
end
%}



%% Plots
FigureStartNum=110;
FigureIndex=1;
if ShowFigureFlag
 %  MeanPSTHSceneOn=nanmean(AlignedPSTH_SceneOn,1);
%SemPSTHSceneOn=nanstd(AlignedPSTH_SceneOn,[],1)./sqrt(sum(~isnan(AlignedPSTH_SceneOn),1));

%MeanPSTHTime=nanmean(AlignedTime_SceneOn,1);

%%%%%%%%%%%%%%%% Next figure : Average responses before scene onset
figtitlestr{FigureIndex}='BeforeSceneOnresponses';
fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[30,100, 1600,800],'Name',figtitlestr{FigureIndex});   
 subplot(2,1,1);
shadedErrorBar(MeanPSTHTime,MeanPSTHSceneOn,SemPSTHSceneOn,'lineprops',{'r'},'transparent',1,'patchSaturation',0.5);
  hold on
  yrange=ylim;
   h1=plot([SceneOnMarker(1), SceneOnMarker(1)],yrange,'--k','LineWidth',1.5);
   
   h2=plot([PreFPTimeRelMean, PreFPTimeRelMean],yrange,'--m','LineWidth',1.5);
   
    h3=plot([ExtraFPOnTimeRelMean, ExtraFPOnTimeRelMean],yrange,'--r','LineWidth',1.5);
    h4=plot([ExtraFPOffTimeRelMean, ExtraFPOffTimeRelMean],yrange,'--b','LineWidth',1.5);
 
   
   
   
   
 set( get( get( h1, 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );
 set( get( get( h2, 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );
 
  set( get( get( h3, 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );
 set( get( get( h4, 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );
 
 
 
 set(findall(gca, 'Type', 'Line'),'LineWidth',3);
  xlabel('Time from scene onset(ms)');
  ylabel('Firing Rate(Hz)');
  

  title('Before SceneOnset');
  set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
  
 subplot(2,1,2); %Scene Off
 shadedErrorBar(MeanPSTHTimeOff,MeanPSTHSceneOff,SemPSTHSceneOff,'lineprops',{'r'},'transparent',1,'patchSaturation',0.5);
  hold on
  yrange=ylim;
   h1=plot([SceneOffMarker(1), SceneOffMarker(1)],yrange,'--k','LineWidth',1.5);
   
  
   
   
 set( get( get( h1, 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );
 
 
 
 
 set(findall(gca, 'Type', 'Line'),'LineWidth',3);
  xlabel('Time from scene offset(ms)');
  ylabel('Firing Rate(Hz)');
  

  title('After SceneOffset');
  set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
  
  
 
 FigureIndex=FigureIndex+1;
FigureStartNum=FigureStartNum+1;  



%%%%%%%%%%%%%%%% Next figure : Average responses before scene onset
figtitlestr{FigureIndex}='BeforeAfterSceneOnresponsesEachBlock';
fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[30,100, 1600,800],'Name',figtitlestr{FigureIndex});   


  subplot(2,1,1); %Scene On across blocks
  
  %Loop over the block types
for i=1:length(UniqueBlockType)
    
    
    shadedErrorBar(MeanPSTHTime,MeanSceneOnResponseBlock(i,:), SemSceneOnResponseBlock(i,:),'lineprops',{ColorForEachMode(find(AssociateBlock==UniqueBlockType(i)),:)},'transparent',1,'patchSaturation',0.5);
  hold on
    
    
    
end
 
 
 yrange=ylim;
   h1=plot([SceneOnMarker(1), SceneOnMarker(1)],yrange,'--k','LineWidth',1.5);
   
   h2=plot([PreFPTimeRelMean, PreFPTimeRelMean],yrange,'--m','LineWidth',1.5);
   
    h3=plot([ExtraFPOnTimeRelMean, ExtraFPOnTimeRelMean],yrange,'--r','LineWidth',1.5);
    h4=plot([ExtraFPOffTimeRelMean, ExtraFPOffTimeRelMean],yrange,'--b','LineWidth',1.5);
 
   
   
   
   
 set( get( get( h1, 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );
 set( get( get( h2, 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );
 
  set( get( get( h3, 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );
 set( get( get( h4, 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );
 
 set(findall(gca, 'Type', 'Line'),'LineWidth',3);
 
  xlabel('Time from scene onset(ms)');
  ylabel('Firing Rate(Hz)');
  

  title('Before SceneOnset each block');
  set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
  
  
  
  subplot(2,1,2); %Scene Odd across blocks
  
  %Loop over the block types
for i=1:length(UniqueBlockType)
    
    
    shadedErrorBar(MeanPSTHTimeOff,MeanSceneOffResponseBlock(i,:), SemSceneOffResponseBlock(i,:),'lineprops',{ColorForEachMode(find(AssociateBlock==UniqueBlockType(i)),:)},'transparent',1,'patchSaturation',0.5);
  hold on
    
    
    
end
 
 
 yrange=ylim;
   h1=plot([SceneOffMarker(1), SceneOffMarker(1)],yrange,'--k','LineWidth',1.5);
   
  
   
   
   
   
 set( get( get( h1, 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );
 
 set(findall(gca, 'Type', 'Line'),'LineWidth',3);
 
  xlabel('Time from scene offset(ms)');
  ylabel('Firing Rate(Hz)');
  

  title('After SceneOffset each block');
  set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
  
 %Add data information
  %  data_msg=sprintf(' %s, SpikeChannel: %d, TrialNum: %d ',FileName,SpikeChannel,TrialNum);
    data_msg=strcat(FileName,' SpikeChannel:',num2str(SpikeChannel),';ChannelID:',num2str(ChannelID),';ClusterID:',num2str(ClusterID),';TrialNum: ',num2str(TrialNum));

    annotation('textbox', [0.17, 0.95, 0.7, 0.04], 'string', data_msg,'FontSize',20,'FontWeight','Bold','LineStyle','none');
    

 FigureIndex=FigureIndex+1;
FigureStartNum=FigureStartNum+1;  
 
 %   
%%%%%%%%%%%%%%%%%% Next figure: Scene Onset Fix Rand Comparison

 figtitlestr{FigureIndex}='SceneOnPSTH_FixRandCompare';
 fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[90,100, 1600,800],'Name',figtitlestr{FigureIndex});
 
 %Scene On Comparison
for i=1:length(UniqueScene)
   
  subplot(4,1,1);
  shadedErrorBar(currPSTHTime_SceneOn_Fix_Mean(i,:),currPSTH_SceneOn_Fix_Mean(i,:),currPSTH_SceneOn_Fix_Sem(i,:),'lineprops',{ColorForScene(i,:)},'transparent',1,'patchSaturation',0.5);
  hold on
  %
  
  %}
  if ~isnan(PreStartPlotTime)
    xlim([-PreStartPlotTime,SceneOnMarker(2)]);
  end
  %Add time marker
  if i==length(UniqueScene)
      hold on
 % plot([SceneOnMarker(1), SceneOnMarker(1)],[0,max(cellfun(@(x,y) max(x+y),SceneOn_Resp_Mean_Each,SceneOn_Resp_Sem_Each))],'-k','Color','#707276','LineWidth',1.5);
 % plot([SceneOnMarker(2), SceneOnMarker(2)],[0,max(cellfun(@(x,y) max(x+y),SceneOn_Resp_Mean_Each,SceneOn_Resp_Sem_Each))],'-k','Color','#707276','LineWidth',1.5);
   
 yrange=ylim;
 h1=plot([SceneOnMarker(1), SceneOnMarker(1)],yrange,'--k','LineWidth',1.5);
  h2=plot([SceneOnMarker(2), SceneOnMarker(2)],yrange,'--k','LineWidth',1.5);
 
 %h1=plot([SceneOnMarker(1), SceneOnMarker(1)],[0,max(currPSTH_SceneOn_Fix_Mean+currPSTH_SceneOn_Fix_Sem,[],'all')],'--k','LineWidth',1.5);
  h2=plot([SceneOnMarker(2), SceneOnMarker(2)],[0,max(currPSTH_SceneOn_Fix_Mean+currPSTH_SceneOn_Fix_Sem,[],'all')],'--k','LineWidth',1.5);
  
  set( get( get( h1, 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );
 set( get( get( h2, 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );
 
  end

  set(findall(gca, 'Type', 'Line'),'LineWidth',3);
  
  
  
  
  
  xlabel('Time from scene onset(ms)');
  ylabel('Firing Rate(Hz)');
  

  title('SceneOnset (Fix Scene Sequence)');
  set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
  
   
  subplot(4,1,3);
  shadedErrorBar(currPSTHTime_SceneOn_Rand_Mean(i,:),currPSTH_SceneOn_Rand_Mean(i,:),currPSTH_SceneOn_Rand_Sem(i,:),'lineprops',{ColorForScene(i,:)},'transparent',1,'patchSaturation',0.5);
  
   if ~isnan(PreStartPlotTime)
    xlim([-PreStartPlotTime,SceneOnMarker(2)]);
   end
  
  %Add time marker
  if i==length(UniqueScene)
      yrange=ylim;
     h3=plot([SceneOnMarker(1), SceneOnMarker(1)],yrange,'--k','LineWidth',1.5);
  h4=plot([SceneOnMarker(2), SceneOnMarker(2)],yrange,'--k','LineWidth',1.5);
 
  set( get( get( h3, 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );
 set( get( get( h4, 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );
  end
  
    
   set(findall(gca, 'Type', 'Line'),'LineWidth',3);
    xlabel('Time from scene onset(ms)');
  ylabel('Firing Rate(Hz)');
  

  title('SceneOnset (Rand Scene Sequence)');
   set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
     hold on
    
 %Response off
  subplot(4,1,2);
  shadedErrorBar(currPSTHTime_SceneOff_Fix_Mean(i,:),currPSTH_SceneOff_Fix_Mean(i,:),currPSTH_SceneOff_Fix_Sem(i,:),'lineprops',{ColorForScene(i,:)},'transparent',1,'patchSaturation',0.5);
  hold on
 
  %{
  if ~isnan(PostrStartPlotTime)
   % xlim([SceneOnMarker(1)]);
  end
%}
  xlim([-500,SceneOffMarker(2)]);
 % ylim([0,15]);
  %Add time marker
  if i==length(UniqueScene)
      hold on
 % plot([SceneOnMarker(1), SceneOnMarker(1)],[0,max(cellfun(@(x,y) max(x+y),SceneOn_Resp_Mean_Each,SceneOn_Resp_Sem_Each))],'-k','Color','#707276','LineWidth',1.5);
 % plot([SceneOnMarker(2), SceneOnMarker(2)],[0,max(cellfun(@(x,y) max(x+y),SceneOn_Resp_Mean_Each,SceneOn_Resp_Sem_Each))],'-k','Color','#707276','LineWidth',1.5);
    h1=plot([SceneOffMarker(1), SceneOffMarker(1)],yrange,'--k','LineWidth',1.5);
 % h2=plot([[SceneOffMarker(1), SceneOnMarker(2)],[0,max(currPSTH_SceneOn_Fix_Mean+currPSTH_SceneOn_Fix_Sem,[],'all')],'--k','LineWidth',1.5);
  
  set( get( get( h1, 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );
% set( get( get( h2, 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );
 
  end
  set(findall(gca, 'Type', 'Line'),'LineWidth',3);
  xlabel('Time from scene offset(ms)');
  ylabel('Firing Rate(Hz)');
  

  title('SceneOffset (Fix Scene Sequence)');
  set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
  
   
  subplot(4,1,4);
  shadedErrorBar(currPSTHTime_SceneOff_Rand_Mean(i,:),currPSTH_SceneOff_Rand_Mean(i,:),currPSTH_SceneOff_Rand_Sem(i,:),'lineprops',{ColorForScene(i,:)},'transparent',1,'patchSaturation',0.5);
  
  
  xlim([-500,SceneOffMarker(2)]);
  %Add time marker
  if i==length(UniqueScene)
      yrange=ylim;
     h3=plot([SceneOffMarker(1), SceneOffMarker(1)],yrange,'--k','LineWidth',1.5);
 % h4=plot([SceneOffMarker(2), SceneOffMarker(2)],[0,max(currPSTH_SceneOff_Rand_Mean+currPSTH_SceneOff_Rand_Sem,[],'all')],'--k','LineWidth',1.5);
 
  set( get( get( h3, 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );
 %set( get( get( h4, 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );
  end
  
    
   set(findall(gca, 'Type', 'Line'),'LineWidth',3);
    xlabel('Time from scene offset(ms)');
  ylabel('Firing Rate(Hz)');
  

  title('SceneOffset (Rand Scene Sequence)');
   set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
     hold on
    
     
     
     
     
     
end
ld=legend([string(UniqueScene)]);
%ld=legend([string(UniqueScene);'SSI']);
ld.Box='off';
ld.Position=[0.9   0.7    0.0854    0.0880];
  

 %Add data information
  %  data_msg=sprintf(' %s, SpikeChannel: %d, TrialNum: %d ',FileName,SpikeChannel,TrialNum);
    data_msg=strcat(FileName,' SpikeChannel:',num2str(SpikeChannel),';ChannelID:',num2str(ChannelID),';ClusterID:',num2str(ClusterID),';TrialNum: ',num2str(TrialNum));

    annotation('textbox', [0.17, 0.95, 0.7, 0.04], 'string', data_msg,'FontSize',20,'FontWeight','Bold','LineStyle','none')
 
FigureIndex=FigureIndex+1;
FigureStartNum=FigureStartNum+1;  

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
    if isfield(OutputFileName,'PrePostSceneSeq')
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

OutputData_New.PrePostSceneSeq.TrialType='GoodTrials';

%%%%%%%%%Save overall scene response characteristics%%%%%%%%%%%%%%%%%%%%%%%
%Organize the data
NameStr={'Empty'
    
 
};



DataLib={NaN
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

OutputData_New.PrePostSceneSeq.Task=TaskType;
OutputData_New.PrePostSceneSeq.TaskCode=TaskCode;
OutputData_New.PrePostSceneSeq.DataStamp=DataStamp;
%Store figures in common folders
ProtocolMarker='PrePostSceneSeq';
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
Task_Old=OutputData.PrePostSceneSeq.TaskCode;

     if Task_Old~=TaskCode
        %Not the same task, add the new dataset into the old one
        OutputData.PrePostSceneSeq(numel(OutputData.PrePostSceneSeq)+1)=OutputData_New.PrePostSceneSeq;
     else
         %If the same task,the same protocol,replace the old one with the new one
         
         OutputData.PrePostSceneSeq=OutputData_New.PrePostSceneSeq;
         
    
    
     end
else
    %Output the current file
   OutputData.PrePostSceneSeq=OutputData_New.PrePostSceneSeq; 
end
if OutputFlag
    cd(OutPath);
    save(OutputFileName,'OutputData');
    disp('Data Exported');
else
    disp('No data export');
end  



end %End of the loop of the spike channel


end %End of the main function

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

function [AlignedPSTH,AlignedTime]=AlignPSTH(psth,time,step,Marker);
AlignedPSTH=NaN*ones(size(psth,1),size(psth,2)*2);
AlignedTime=AlignedPSTH;
[Row,Col]=find(abs(time-Marker)<floor(step));
MaxAlignBin=max(Col);
for i=1:size(psth,1)
    fillin=psth(i,:);
    fillinTime=time(i,:);
    MarkerLoc=find(abs(fillinTime-Marker)<floor(step),1,'first');
    AlignedPSTH(i,MaxAlignBin-MarkerLoc+1:MaxAlignBin-MarkerLoc+length(fillin))=fillin;
    AlignedTime(i,MaxAlignBin-MarkerLoc+1:MaxAlignBin-MarkerLoc+length(fillinTime))=fillinTime;
 
end
%Remove unnessary cols
Col_Valid=sum(~isnan(AlignedPSTH),1);
RemoveNaNCol=Col_Valid<floor(sqrt(size(AlignedPSTH,1)));
AlignedPSTH(:,RemoveNaNCol)=[];
AlignedTime(:,RemoveNaNCol)=[];



end




