function FractalStimEffect(Data, StartTrial, EndTrial,ShowFigureFlag,OutputFlag);

ECode;




%Criterias for sexwlecting saccades
V_Threshold=40;
ContinueBin=5;
ISI_Threshold=20;

%Time bin pre saccade and time bin offset after saccade onset
PlotSacPreBuffer=200;%in ms
PlotSacPostBuffer=300;%in ms



%Color Definition for each fractals
%{
Color=[lines(6);prism(6) ];
Color(2,:)=copper(1)+0.5;
Color=[Color;Color*0.5];
%}
Color={'#648FFF','#FE6100','#009E73','#88CCEE','#CC6677','#785EF0'};
Color_RGB=[100,143,255;
           254,97,0;
           0,158,115;
          136,204,238;
          204,102,119;
          120,94,240]/255;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Unpack all the data base and task settings
ECode;
FileName=Data.FileName;
TaskType=Data.TaskType;
TaskCode=Data.TaskCode;
TrialNum=EndTrial-StartTrial+1;
%Selection=StartTrial:EndTrial;
Selection=[StartTrial:EndTrial];


%Abnormal check
CheckCodePair=[TGTCD,TGTOFF];%This code should be the same as target on for good trials
AbnormalTrials=AbnormalTrialCheck(Data,CheckCodePair);
Selection=Selection(~ismember(Selection,AbnormalTrials));
    


%%
%GetEyeData
EyeData=Data.EyeDataRaw;
EyeChannel_X_L=EyeData(Selection,1);
EyeChannel_Y_L=EyeData(Selection,2);
%EyeChannel_X_R=EyeData(Selection,3);
%EyeChannel_Y_R=EyeData(Selection,4);

Eye_Eccentricity=cellfun(@(x,y) sqrt(x.^2+y.^2),EyeChannel_X_L,EyeChannel_Y_L,'uniform',0);

EyeBinWidth=Data.EyeBinWidth;


EventChannel=Data.EventChannel(Selection,:);
EventTimeChannel=Data.EventTimeChannel(Selection,:);
EventBin=Data.EventBin(Selection,:);




%Target Information
TargetIndex=Data.TargetIndex(Selection,2);%For this task, there is only target name, but no target index

TargetName=Data.TargetIndex(Selection,3);%For this task, there is only target name, but no target index

%TargetLoc=Data.TargetLoc(Selection,:,:);


%Rearrange target into lines
TargetIndex_Line=TargetIndex;

TargetIndex_Line=reshape(TargetIndex_Line,numel(TargetIndex_Line),1);

NoNNaNSelector=~isnan(TargetIndex_Line);

TargetIndex_Line=TargetIndex_Line(~isnan(TargetIndex_Line));


TargetName_Line=TargetName;

TargetName_Line=reshape(TargetName_Line,numel(TargetName_Line),1);

TargetName_Line=TargetName_Line(~isnan(TargetName_Line));

%Rearrange target locations into lines
%% Target location

TargetLoc=Data.TargetLoc(Selection,2:end);
TargetPolar=Data.TargetPolarLoc(Selection,2:end);
TargetAngle = TargetPolar(:,1);
TargetEcc =  TargetPolar(:,2);

UniqueAngle = unique(TargetAngle);
UniqueEcc = unique(TargetEcc);
minEcc = min(TargetEcc);

TargetAngle_Line=TargetAngle;
TargetAngle_Line(TargetAngle_Line==360)=0;
TargetAngle_Line=reshape(TargetAngle_Line,numel(TargetAngle_Line),1);
TargetAngle_Line=TargetAngle_Line(NoNNaNSelector);

TargetEcc_Line=TargetEcc;
TargetEcc_Line=reshape(TargetEcc_Line,numel(TargetEcc_Line),1);
TargetEcc_Line=TargetEcc_Line(NoNNaNSelector);



UniqueAngle = unique(TargetAngle_Line);
UniqueEcc = unique(TargetEcc_Line);





%%

EventChannel=Data.EventChannel(Selection,:);
EventTimeChannel=Data.EventTimeChannel(Selection,:);
EventBin=Data.EventBin(Selection,:);


%% Information about stimulation

%Load stimulation related data

ParaLib=Data.ParaLib;
if ismember('StimType',keys(ParaLib))
    StimTime=ParaLib('StimType');
    StimTime= StimTime(Selection,:);
    StimTime= StimTime(:,2);
    UnqiueStimTime = unique(StimTime);

end

%Seperate into stim trial and sham trial
StimType=cell2mat(ReproduceFromEvent(EventChannel,[STIM1,SHAM0]))';

StimTrial = StimType == STIM1;
ControlTrial = ~StimTrial;

%%Reproduce goodbad index from event channel
%ObjectIndex=
ObjectIndex=cell2mat(ReproduceFromEvent(EventChannel,[GOODOBJ  BADOBJ]));

GoodObjectIndex=(ObjectIndex==GOODOBJ)';


%%Analysis eye trace
%Select the period between fix off and target off set
FIXOFFBin = FindOutTime(EventChannel,EventBin,FIXOFF );
SACENDBin = FIXOFFBin + 1000;
%{
FractalOffsetBin = FindOutTime(EventChannel,EventBin,TGTOFF);
FractalOnsetTime = FindOutTime(EventChannel,EventTimeChannel,TGTCD );
FractalOffsetTime = FindOutTime(EventChannel,EventTimeChannel,TGTOFF );
%}
EyeChannel_X_Selection= SelectEyeInterval(EyeChannel_X_L,FIXOFFBin,SACENDBin)';
EyeChannel_Y_Selection= SelectEyeInterval(EyeChannel_Y_L,FIXOFFBin,SACENDBin)';

EyeChannel_X_Selection=ReorganizeEye(EyeChannel_X_Selection);
EyeChannel_Y_Selection=ReorganizeEye(EyeChannel_Y_Selection);


%Analysis saccades for each trial
Saccades=SelectSaccade(EyeChannel_X_Selection,EyeChannel_Y_Selection,EyeBinWidth,V_Threshold,ContinueBin,ISI_Threshold);
SaccadeNumber=arrayfun(@(x)  x.NumOfSaccade,Saccades)';
SaccadeAngle=arrayfun(@(x)  x.SaccadeAngle,Saccades,'UniformOutput' ,0)';
SaccadeAmplitude=arrayfun(@(x)  x.SaccadeAmplitude,Saccades,'UniformOutput' ,0)';

SaccadeStartTime=arrayfun(@(x)  x.SaccadeStartTime,Saccades,'UniformOutput' ,0)';
SaccadeEndTime=arrayfun(@(x)  x.SaccadeEndTime,Saccades,'UniformOutput' ,0)';

SaccadeStartPoint=arrayfun(@(x)  x.SaccadeStartPoint,Saccades,'UniformOutput' ,0)';
SaccadeEndPoint=arrayfun(@(x)  x.SaccadeEndPoint,Saccades,'UniformOutput' ,0)';

SaccadePeakVelocity=arrayfun(@(x)  x.PeakV,Saccades,'UniformOutput' ,0)';

%Screen out the valid saccades and saccade start time above 5;

AmplitudeThreshold=minEcc/2;

TimeThreshold=50;

SaccadeAngle_Selection=cellfun(@(x,y,z) z(((x>AmplitudeThreshold)&(y>TimeThreshold))),SaccadeAmplitude,SaccadeStartTime,SaccadeAngle,'UniformOutput' ,0);
SaccadeEndPoint_Selection=cellfun(@(x,y,z) z(((x>AmplitudeThreshold)&(y>TimeThreshold)),:),SaccadeAmplitude,SaccadeStartTime,SaccadeEndPoint,'UniformOutput' ,0);
SaccadeStartPoint_Selection=cellfun(@(x,y,z) z(((x>AmplitudeThreshold)&(y>TimeThreshold)),:),SaccadeAmplitude,SaccadeStartTime,SaccadeStartPoint,'UniformOutput' ,0);

SaccadeStartTime_Selection=cellfun(@(x,y,z) z(((x>AmplitudeThreshold)&(y>TimeThreshold)),:),SaccadeAmplitude,SaccadeStartTime,SaccadeStartTime,'UniformOutput' ,0);
SaccadeEndTime_Selection=cellfun(@(x,y,z) z(((x>AmplitudeThreshold)&(y>TimeThreshold)),:),SaccadeAmplitude,SaccadeStartTime,SaccadeEndTime,'UniformOutput' ,0);


SaccadeAmplitude_Selection=cellfun(@(x,y,z) z(((x>AmplitudeThreshold)&(y>TimeThreshold)),:),SaccadeAmplitude,SaccadeStartTime,SaccadeAmplitude,'UniformOutput' ,0);

SaccadePeakVelocity_Selection = cellfun(@(x,y,z) z(((x>AmplitudeThreshold)&(y>TimeThreshold)),:),SaccadeAmplitude,SaccadeStartTime,SaccadePeakVelocity,'UniformOutput' ,0);

%Only look at the first saccade
SaccadeNumInEachTrial=cellfun(@numel,SaccadeAngle_Selection);


SaccadeAngle_First=NaN*ones(1,length(SaccadeNumInEachTrial));
SaccadeEndPoint_First=NaN*ones(length(SaccadeNumInEachTrial),2);
RT_First_Target=NaN*ones(1,length(SaccadeNumInEachTrial));
SacEnd_First_Target=NaN*ones(1,length(SaccadeNumInEachTrial));

SaccadeStartPoint_First=NaN*ones(length(SaccadeNumInEachTrial),2);

SaccadeAmplitude_First=NaN*ones(1,length(SaccadeNumInEachTrial));

SaccadeAngle_First(SaccadeNumInEachTrial>0)=cell2mat(cellfun(@(x) x(1,:),SaccadeAngle_Selection(SaccadeNumInEachTrial>0),'UniformOutput' ,0));
SaccadeEndPoint_First(SaccadeNumInEachTrial>0,:)=cell2mat(cellfun(@(x) x(1,:),SaccadeEndPoint_Selection(SaccadeNumInEachTrial>0,:),'UniformOutput' ,0));
RT_First_Target(SaccadeNumInEachTrial>0)=cell2mat(cellfun(@(x) x(1,:),SaccadeStartTime_Selection(SaccadeNumInEachTrial>0,:),'UniformOutput' ,0));
SacEnd_First_Target(SaccadeNumInEachTrial>0)=cell2mat(cellfun(@(x) x(1,:),SaccadeEndTime_Selection(SaccadeNumInEachTrial>0,:),'UniformOutput' ,0));


SaccadeStartPoint_First(SaccadeNumInEachTrial>0,:)=cell2mat(cellfun(@(x) x(1,:),SaccadeStartPoint_Selection(SaccadeNumInEachTrial>0,:),'UniformOutput' ,0));

SaccadeAmplitude_First(SaccadeNumInEachTrial>0)=cell2mat(cellfun(@(x) x(1,:),SaccadeAmplitude_Selection(SaccadeNumInEachTrial>0,:),'UniformOutput' ,0));

SaccadePeakVelocity_First(SaccadeNumInEachTrial>0)=cell2mat(cellfun(@(x) x(1,:),SaccadePeakVelocity_Selection(SaccadeNumInEachTrial>0,:),'UniformOutput' ,0));

RT_TargetToFixOff = RT_First_Target;

%StimIntervalSac = RT_TargetToFixOff -  

 for i = 1:length(UniqueAngle)
     sel_good_stim = StimTrial &  TargetAngle == UniqueAngle(i) & GoodObjectIndex;
     sel_bad_stim = StimTrial &  TargetAngle == UniqueAngle(i) & ~GoodObjectIndex;
     
     sel_good_control = ControlTrial &  TargetAngle == UniqueAngle(i) & GoodObjectIndex;
     sel_bad_control = ControlTrial &  TargetAngle == UniqueAngle(i) & ~GoodObjectIndex;
     
    RT_good_stim{i} = RT_TargetToFixOff(sel_good_stim);
    RT_bad_stim{i} = RT_TargetToFixOff(sel_bad_stim);
    
    RT_good_control{i} = RT_TargetToFixOff(sel_good_control);
    RT_bad_control{i} = RT_TargetToFixOff(sel_bad_control);
    
    
   PeakV_good_stim{i} =SaccadePeakVelocity_First(sel_good_stim);
   PeakV_good_control{i} =SaccadePeakVelocity_First(sel_good_control);
   
    PeakV_bad_stim{i} =SaccadePeakVelocity_First(sel_bad_stim);
   PeakV_bad_control{i} =SaccadePeakVelocity_First(sel_bad_control);
   
   
    
    
     RT_good_stim_mean_loc(i) = nanmean(RT_good_stim{i});
     RT_bad_stim_mean_loc(i) = nanmean(RT_bad_stim{i});
     
     RT_good_stim_sem_loc(i) = nanstd(RT_good_stim{i})/sqrt(sum(~isnan(RT_good_stim{i})));
     RT_bad_stim_sem_loc(i) = nanstd(RT_bad_stim{i})/sqrt(sum(~isnan(RT_bad_stim{i})));
     
     RT_good_control_mean_loc(i) = nanmean(RT_good_control{i});
     RT_bad_control_mean_loc(i) = nanmean(RT_bad_control{i});
     
      RT_good_control_sem_loc(i) = nanstd(RT_good_control{i})/sqrt(sum(~isnan(RT_good_control{i})));
     RT_bad_control_sem_loc(i) = nanstd(RT_bad_control{i})/sqrt(sum(~isnan(RT_bad_control{i})));
     
     
     PV_good_stim_mean_loc(i) = nanmean(PeakV_good_stim{i});
     PV_bad_stim_mean_loc(i) = nanmean(PeakV_bad_stim{i});
     
     PV_good_stim_sem_loc(i) = nanstd(PeakV_good_stim{i})/sqrt(sum(~isnan(PeakV_good_stim{i})));
     PV_bad_stim_sem_loc(i) = nanstd(PeakV_bad_stim{i})/sqrt(sum(~isnan(PeakV_bad_stim{i})));
     
     PV_good_control_mean_loc(i) = nanmean(PeakV_good_control{i});
     PV_bad_control_mean_loc(i) = nanmean(PeakV_bad_control{i});
     
      PV_good_control_sem_loc(i) = nanstd(PeakV_good_stim{i})/sqrt(sum(~isnan(PeakV_good_stim{i})));
     PV_bad_control_sem_loc(i) = nanstd(PeakV_bad_control{i})/sqrt(sum(~isnan(PeakV_bad_control{i})));
     
 end
 

 
 
%% For neural responses 
% Align on three period: STIMON, TgtOn, FixOff

StepSize=5;%in ms
BinWidth=10;%in ms

SegEvent=[TGTCD];%Event Flow of interest
%SegEvent=[ON_STIM];
%SegEventEnd=[SHOWFIXCD,FIXOFF,TGTOFF,RWDOFFCD ,NaN,EXTRAFP];%Event end
SegEventEnd=[NaN];%Event end
%SegEventEnd=[PST_STIM];%Event end
SegEventStr={'TargetOn'};
%SegEventStr={'StimOn'};
PlotStartOffset=[200];%in ms
PlotStopOffset=[200];%in mss

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

%sSpikeTime=squeeze(Data.SpikeTimeData(:,:,Selection))';
%SpikeBin=squeeze(Data.SpikeBinData(:,:,Selection))';



%%
%%%%%%%%%%%%%%%%%%%%%%%Align data on object onset
%Segregate spike data from fractal onset to fractal offset

SegmentData=EventSegTool(SegEvent,SegEventEnd,PlotStartOffset,PlotStopOffset,EventChannel,EventTimeChannel,SpikeTime,StepSize,BinWidth);

%Retrive data from the field

%Relative time bin for each event
RelativeTimeMarker=SegmentData.RelativeTimeMarker{1};%=EventTimeMarker;
%Absolute time marker for each event
AbosoluteTimeMarker=SegmentData.AbsoluteTimeMarker{1};

ObjOnTimeMarker=nanmean(RelativeTimeMarker);

%{
RwdTimeAbs=FindOutTime(EventChannel,EventBin,RWDCD);
RwdTimeRel=RwdTimeAbs-AbosoluteTimeMarker(:,1);

%Algin marker to sac-on & sac-off
RwdTimeRelSacEnd=RwdTimeRel-SacEnd_First_Target';
TgtOffTimeRelSacEnd=RelativeTimeMarker(:,2)-SacEnd_First_Target';

RwdTimeRelSacEndMean=nanmean(RwdTimeRelSacEnd);
TgtOffTimeRelSacEndMean=nanmean(TgtOffTimeRelSacEnd);
%}
%PSTH for each event
MeanPSTH=SegmentData.MeanPSTH;%=[PSTH_Event_Time',PSTH_Event'];
MeanPSTH_Time=MeanPSTH{1};
MeanPSTH_FR_FractalOn=MeanPSTH{2};

PSTH_Time_FraOn_All=nanmean(MeanPSTH_Time,1);


%% Align responses on sac-onset
%Align neural responses according to the saccade onset time
%

%Organize the spike time using the same way
%Align spike time with target onset
FractalOnsetTime = FindOutTime(EventChannel,EventTimeChannel,TGTCD );
StimOnsetTime = FindOutTime(EventChannel,EventTimeChannel,ON_STIM);
StimOffsetTime = FindOutTime(EventChannel,EventTimeChannel,OFF_STIM);
DiffStimTgt = StimOnsetTime -FractalOnsetTime;
DiffOffStimTgt = StimOffsetTime -FractalOnsetTime;

StimIntervalSac =[ nanmean(DiffStimTgt-RT_TargetToFixOff' ),nanmean(DiffOffStimTgt-RT_TargetToFixOff' )]; 

StimInterval = [nanmean(DiffStimTgt),nanmean(DiffOffStimTgt )];

AlignedSpike=SpikeAlignTool(FractalOnsetTime,SpikeTime,StepSize,BinWidth);

SpikeTimeAligned_TgtOn=AlignedSpike.SpikeTime_Aligned_Whole;


SpikeTimeAligned_TgtOn=SpikeTimeAligned_TgtOn(NoNNaNSelector,:);

%Align again on saccade onset
AlignedSpike=SpikeAlignTool(RT_First_Target',SpikeTimeAligned_TgtOn,StepSize,BinWidth);

PSTH_SacOn=AlignedSpike.PSTH_Aligned;
TimeSequence_SacOn=AlignedSpike.PSTH_Time_Aligned;
%Select interval
SelectPSTHinterval_SacOn=[-PlotSacPreBuffer,PlotSacPostBuffer];
[PSTH_SacOn,TimeSequence_SacOn]=SelectPSTHInterval(PSTH_SacOn,TimeSequence_SacOn,SelectPSTHinterval_SacOn);



TimeSequence_SacOn_mean=nanmean(TimeSequence_SacOn,1);


%Good versus Bad
PSTH_SacOn_Good_mean=nanmean(PSTH_SacOn(GoodObjectIndex,:),1);
PSTH_SacOn_Good_sem=nanstd(PSTH_SacOn(GoodObjectIndex,:),[],1)./sqrt(sum(~isnan(PSTH_SacOn(GoodObjectIndex,1)),1));

PSTH_SacOn_Bad_mean=nanmean(PSTH_SacOn(~GoodObjectIndex,:),1);
PSTH_SacOn_Bad_sem=nanstd(PSTH_SacOn(~GoodObjectIndex,:),[],1)./sqrt(sum(~isnan(PSTH_SacOn(~GoodObjectIndex,1)),1));


%% Loop over the target directions
for i = 1:length(UniqueAngle )
   
        sel = TargetAngle_Line == UniqueAngle(i);
        
        %Align on target on 
        curr_Good_Loc_stim = MeanPSTH_FR_FractalOn(sel & GoodObjectIndex & StimTrial,:);
        curr_Bad_Loc_stim = MeanPSTH_FR_FractalOn(sel & ~GoodObjectIndex & StimTrial,:);
        
        curr_Good_Loc_control = MeanPSTH_FR_FractalOn(sel & GoodObjectIndex & ControlTrial,:);
        curr_Bad_Loc_control = MeanPSTH_FR_FractalOn(sel & ~GoodObjectIndex & ControlTrial,:);
        
        
        
        Good_FracOn_Loc_Mean_stim{i} = nanmean(curr_Good_Loc_stim,1);
        Good_FracOn_Loc_Sem_stim{i} = nanstd(curr_Good_Loc_stim,[],1)./sqrt(sum(~isnan(curr_Good_Loc_stim),1));
        
        
        Bad_FracOn_Loc_Mean_stim{i} = nanmean(curr_Bad_Loc_stim,1);
        Bad_FracOn_Loc_Sem_stim{i} = nanstd(curr_Bad_Loc_stim,[],1)./sqrt(sum(~isnan(curr_Bad_Loc_stim),1));
        
        Good_FracOn_Loc_Mean_control{i} = nanmean(curr_Good_Loc_control,1);
        Good_FracOn_Loc_Sem_control{i} = nanstd(curr_Good_Loc_control,[],1)./sqrt(sum(~isnan(curr_Good_Loc_control),1));
        
        
        Bad_FracOn_Loc_Mean_control{i} = nanmean(curr_Bad_Loc_control,1);
        Bad_FracOn_Loc_Sem_control{i} = nanstd(curr_Bad_Loc_control,[],1)./sqrt(sum(~isnan(curr_Bad_Loc_control),1));
        
        %Align on saccade on
        curr_Good_sac_Loc_stim = PSTH_SacOn(sel & GoodObjectIndex & StimTrial,:);
        curr_Bad_sac_Loc_stim = PSTH_SacOn(sel & ~GoodObjectIndex & StimTrial,:);
        
        curr_Good_sac_Loc_control = PSTH_SacOn(sel & GoodObjectIndex & ControlTrial,:);
        curr_Bad_sac_Loc_control =PSTH_SacOn(sel & ~GoodObjectIndex & ControlTrial,:);
        
        Good_SacOn_Loc_Mean_stim{i} = nanmean(curr_Good_sac_Loc_stim,1);
        Good_SacOn_Loc_Sem_stim{i} = nanstd(curr_Good_sac_Loc_stim,[],1)./sqrt(sum(~isnan(curr_Good_sac_Loc_stim),1));
        
        
        Bad_SacOn_Loc_Mean_stim{i} = nanmean(curr_Bad_sac_Loc_stim,1);
        Bad_SacOn_Loc_Sem_stim{i} = nanstd(curr_Bad_sac_Loc_stim,[],1)./sqrt(sum(~isnan(curr_Bad_sac_Loc_stim),1));
        
        Good_SacOn_Loc_Mean_control{i} = nanmean(curr_Good_sac_Loc_control ,1);
        Good_SacOn_Loc_Sem_control{i} = nanstd(curr_Good_sac_Loc_control ,[],1)./sqrt(sum(~isnan(curr_Good_sac_Loc_control ),1));
        
        
        Bad_SacOn_Loc_Mean_control{i} = nanmean(curr_Bad_sac_Loc_control,1);
        Bad_SacOn_Loc_Sem_control{i} = nanstd(curr_Bad_sac_Loc_control,[],1)./sqrt(sum(~isnan(curr_Bad_sac_Loc_control),1));
        
       
        
        
        
      %  AngleEccorder{i} = [UniqueAngle(i)];
        
        
        
    
end

 if ShowFigureFlag
     
   StartNO = 60;
   FigureIndex=1;
   figtitlestr{FigureIndex}='RTDifference';
   fig{FigureIndex}=PrepareFigure(StartNO+FigureIndex-1,'w',[50,100, 1400,400],'Name',figtitlestr{FigureIndex});
   
   subplot(1,3,1)
   g=bar([RT_good_stim_mean_loc(1);RT_good_control_mean_loc(1);RT_bad_stim_mean_loc(1);RT_bad_control_mean_loc(1)],'LineWidth',3);  
   g.FaceColor = 'flat';
    %Color=rand(NumOfCondition,3)
   g.CData = [1,0,0;0.2,0,0;0,0,1;0,0,0.2];
   
   hold on
   errorbar([RT_good_stim_mean_loc(1);RT_good_control_mean_loc(1);RT_bad_stim_mean_loc(1);RT_bad_control_mean_loc(1)],[RT_good_stim_sem_loc(1);RT_good_control_sem_loc(1);RT_bad_stim_sem_loc(1);RT_bad_control_sem_loc(1)],'.k','LineWidth',3);
   
   xticks(1:4);
   xticklabels({'GoodStim','GooNonStim','BadStim','BadNonStim'});
   ylabel('ReactionTime(ms)');
   set(gca,'LineWidth',3,'FontSize',10,'FontWeight','Bold');
   
   title(sprintf('Angle%d',UniqueAngle(1)),'FontSize',15,'FontWeight','Bold');
   box off
   
    subplot(1,3,2)
   g=bar([RT_good_stim_mean_loc(2);RT_good_control_mean_loc(2);RT_bad_stim_mean_loc(2);RT_bad_control_mean_loc(2)],'LineWidth',3);  
   g.FaceColor = 'flat';
    %Color=rand(NumOfCondition,3)
   g.CData = [1,0,0;0.2,0,0;0,0,1;0,0,0.2];
   
   hold on
   errorbar([RT_good_stim_mean_loc(2);RT_good_control_mean_loc(2);RT_bad_stim_mean_loc(2);RT_bad_control_mean_loc(2)],[RT_good_stim_sem_loc(2);RT_good_control_sem_loc(2);RT_bad_stim_sem_loc(2);RT_bad_control_sem_loc(2)],'.k','LineWidth',3);
   
   xticks(1:4);
   xticklabels({'GoodStim','GoodNonStim','BadStim','BadNonStim'});
   ylabel('ReactionTime(ms)');
   set(gca,'LineWidth',3,'FontSize',10,'FontWeight','Bold');
   
   title(sprintf('Angle%d',UniqueAngle(2)),'FontSize',15,'FontWeight','Bold');
   box off
   
    subplot(1,3,3)
 
   g = bar([RT_good_stim_mean_loc(1)-RT_good_control_mean_loc(1),RT_bad_stim_mean_loc(1)-RT_bad_control_mean_loc(1),RT_good_stim_mean_loc(2)-RT_good_control_mean_loc(2),RT_bad_stim_mean_loc(2)-RT_bad_control_mean_loc(2)]);
   
   g.FaceColor = 'flat';
    %Color=rand(NumOfCondition,3)
   g.CData = [1,0,0;0,0,1;1,0,0;0,0,1];
   xticks(1:4);
   xticklabels({sprintf('Good%d',UniqueAngle(1)),sprintf('Bad%d',UniqueAngle(1)),sprintf('Good%d',UniqueAngle(2)),sprintf('Bad%d',UniqueAngle(2))});
   
   ylabel('RT,Stim-NonStim');
   set(gca,'LineWidth',3,'FontSize',10,'FontWeight','Bold');
   box off
   
   %{
    subplot(2,3,4)
   g=bar([PV_good_stim_mean_loc(1);PV_good_control_mean_loc(1);PV_bad_stim_mean_loc(1);PV_bad_control_mean_loc(1)],'LineWidth',3);  
   g.FaceColor = 'flat';
    %Color=rand(NumOfCondition,3)
   g.CData = [1,0,0;0.2,0,0;0,0,1;0,0,0.2];
   
   hold on
   errorbar([PV_good_stim_mean_loc(1);PV_good_control_mean_loc(1);PV_bad_stim_mean_loc(1);PV_bad_control_mean_loc(1)],[PV_good_stim_sem_loc(1);PV_good_control_sem_loc(1);PV_bad_stim_sem_loc(1);PV_bad_control_sem_loc(1)],'.k','LineWidth',3);
   
   xticks(1:4);
   xticklabels({'GoodStim','GooNonStim','BadStim','BadNonStim'});
   ylabel('Peak Velocity');
   set(gca,'LineWidth',3,'FontSize',10,'FontWeight','Bold');
   
   title(sprintf('Angle%d',UniqueAngle(1)),'FontSize',15,'FontWeight','Bold');
   box off
   
    subplot(2,3,5)
   g=bar([PV_good_stim_mean_loc(2);PV_good_control_mean_loc(2);PV_bad_stim_mean_loc(2);PV_bad_control_mean_loc(2)],'LineWidth',3);  
   g.FaceColor = 'flat';
    %Color=rand(NumOfCondition,3)
   g.CData = [1,0,0;0.2,0,0;0,0,1;0,0,0.2];
   
   hold on
   errorbar([PV_good_stim_mean_loc(2);PV_good_control_mean_loc(2);PV_bad_stim_mean_loc(2);PV_bad_control_mean_loc(2)],[PV_good_stim_sem_loc(2);PV_good_control_sem_loc(2);PV_bad_stim_sem_loc(2);PV_bad_control_sem_loc(2)],'.k','LineWidth',3);
   
   xticks(1:4);
   xticklabels({'GoodStim','GoodNonStim','BadStim','BadNonStim'});
   ylabel('Peak Velocity');
   set(gca,'LineWidth',3,'FontSize',10,'FontWeight','Bold');
   
   title(sprintf('Angle%d',UniqueAngle(2)),'FontSize',15,'FontWeight','Bold');
   box off
   
    subplot(2,3,6)
 
   g = bar([PV_good_stim_mean_loc(1)-PV_good_control_mean_loc(1),PV_bad_stim_mean_loc(1)-PV_bad_control_mean_loc(1),PV_good_stim_mean_loc(2)-PV_good_control_mean_loc(2),PV_bad_stim_mean_loc(2)-PV_bad_control_mean_loc(2)]);
   
   g.FaceColor = 'flat';
    %Color=rand(NumOfCondition,3)
   g.CData = [1,0,0;0,0,1;1,0,0;0,0,1];
   xticks(1:4);
   xticklabels({sprintf('Good%d',UniqueAngle(1)),sprintf('Bad%d',UniqueAngle(1)),sprintf('Good%d',UniqueAngle(2)),sprintf('Bad%d',UniqueAngle(2))});
   
   ylabel('Peak Velocity,Stim-NonStim');
   set(gca,'LineWidth',3,'FontSize',10,'FontWeight','Bold');
   box off
   
  %}
   
   StartNO = StartNO+1;
   FigureIndex = FigureIndex + 1;
   figtitlestr{FigureIndex}='PSTHbetweenStimNonStims';
   fig{FigureIndex}=PrepareFigure(StartNO+FigureIndex-1,'w',[50,100, 1400,400],'Name',figtitlestr{FigureIndex});
   subplot(2,4,1)
  % StimInterval 
  maxy = max(max([Good_FracOn_Loc_Mean_stim{1} + Good_FracOn_Loc_Sem_stim{1}]),max([Good_FracOn_Loc_Mean_control{1}+Good_FracOn_Loc_Sem_control{1}]));
   area(StimInterval,[maxy,maxy],...
    'FaceColor','k','FaceAlpha',.1,'EdgeAlpha',.1,'HandleVisibility','off');
hold on
   box off;
   
   shadedErrorBar(PSTH_Time_FraOn_All, Good_FracOn_Loc_Mean_stim{1}, Good_FracOn_Loc_Sem_stim{1},'lineprops',{[1    0    0]},'transparent',1,'patchSaturation',0.3);
   hold on
   shadedErrorBar(PSTH_Time_FraOn_All, Good_FracOn_Loc_Mean_control{1}, Good_FracOn_Loc_Sem_control{1},'lineprops',{[0    0    0]},'transparent',1,'patchSaturation',0.3);
    yrange=ylim;
    plot([0,0],yrange,'--k');
    box off
     set(findall(gca, 'Type', 'Line'),'LineWidth',3);
     title(sprintf('Angle%d',UniqueAngle(1)));
    set(gca,'FontSize',15,'FontWeight','Bold');
    set(gca,'LineWidth',3);

    xlabel('Time from Good Onset','FontSize',15,'FontWeight','Bold');
    ylabel('FR','FontSize',15,'FontWeight','Bold');
    
     subplot(2,4,2)
      % StimInterval 
  maxy = max(max([Good_SacOn_Loc_Mean_stim{1} + Good_SacOn_Loc_Sem_stim{1}]),max([Good_SacOn_Loc_Mean_control{1}+Good_SacOn_Loc_Sem_control{1}]));
   area( StimIntervalSac,[maxy,maxy],...
    'FaceColor','k','FaceAlpha',.1,'EdgeAlpha',.1,'HandleVisibility','off');
hold on
   box off;
     
    
   
   shadedErrorBar(TimeSequence_SacOn_mean, Good_SacOn_Loc_Mean_stim{1}, Good_SacOn_Loc_Sem_stim{1},'lineprops',{[1    0    0]},'transparent',1,'patchSaturation',0.3);
   hold on
   shadedErrorBar(TimeSequence_SacOn_mean, Good_SacOn_Loc_Mean_control{1}, Good_SacOn_Loc_Sem_control{1},'lineprops',{[0    0    0]},'transparent',1,'patchSaturation',0.3);
    yrange=ylim;
    plot([0,0],yrange,'--k');

    title(sprintf('Angle%d',UniqueAngle(1)));
    set(findall(gca, 'Type', 'Line'),'LineWidth',3);
    set(gca,'FontSize',15,'FontWeight','Bold');
    set(gca,'LineWidth',3);

    xlabel('Time from Sac Onset','FontSize',15,'FontWeight','Bold');
    ylabel('FR','FontSize',15,'FontWeight','Bold');
    
     subplot(2,4,3)
     
      maxy = max(max([Good_FracOn_Loc_Mean_stim{2} + Good_FracOn_Loc_Sem_stim{2}]),max([Good_FracOn_Loc_Mean_control{2}+Good_FracOn_Loc_Sem_control{2}]));
   area(StimInterval,[maxy,maxy],...
    'FaceColor','k','FaceAlpha',.1,'EdgeAlpha',.1,'HandleVisibility','off');
hold on
box off
   
   shadedErrorBar(PSTH_Time_FraOn_All, Good_FracOn_Loc_Mean_stim{2}, Good_FracOn_Loc_Sem_stim{2},'lineprops',{[1    0    0]},'transparent',1,'patchSaturation',0.3);
   hold on
   shadedErrorBar(PSTH_Time_FraOn_All, Good_FracOn_Loc_Mean_control{2}, Good_FracOn_Loc_Sem_control{2},'lineprops',{[0    0    0]},'transparent',1,'patchSaturation',0.3);
    yrange=ylim;
    plot([0,0],yrange,'--k');





title(sprintf('Angle%d',UniqueAngle(2)));
    set(findall(gca, 'Type', 'Line'),'LineWidth',3);
    set(gca,'FontSize',15,'FontWeight','Bold');
    set(gca,'LineWidth',3);

    xlabel('Time from Good Onset','FontSize',15,'FontWeight','Bold');
    ylabel('FR','FontSize',15,'FontWeight','Bold');
    
    
     subplot(2,4,4)
     
     % StimInterval 
  maxy = max(max([Good_SacOn_Loc_Mean_stim{2} + Good_SacOn_Loc_Sem_stim{2}]),max([Good_SacOn_Loc_Mean_control{2}+Good_SacOn_Loc_Sem_control{2}]));
   area( StimIntervalSac,[maxy,maxy],...
    'FaceColor','k','FaceAlpha',.1,'EdgeAlpha',.1,'HandleVisibility','off');
hold on
   box off;
     hold on
   
   shadedErrorBar(TimeSequence_SacOn_mean, Good_SacOn_Loc_Mean_stim{2}, Good_SacOn_Loc_Sem_stim{2},'lineprops',{[1    0    0]},'transparent',1,'patchSaturation',0.3);
   hold on
   shadedErrorBar(TimeSequence_SacOn_mean, Good_SacOn_Loc_Mean_control{2}, Good_SacOn_Loc_Sem_control{2},'lineprops',{[0    0    0]},'transparent',1,'patchSaturation',0.3);
    yrange=ylim;
    plot([0,0],yrange,'--k');





title(sprintf('Angle%d',UniqueAngle(2)));
    set(findall(gca, 'Type', 'Line'),'LineWidth',3);
    set(gca,'FontSize',15,'FontWeight','Bold');
    set(gca,'LineWidth',3);

    xlabel('Time from Sac Onset','FontSize',15,'FontWeight','Bold');
    ylabel('FR','FontSize',15,'FontWeight','Bold');
    
    
    
    
   
    subplot(2,4,5)
    
     % StimInterval 
  maxy = max(max([Bad_FracOn_Loc_Mean_stim{1} + Bad_FracOn_Loc_Sem_stim{1}]),max([Bad_FracOn_Loc_Mean_control{1}+Bad_FracOn_Loc_Sem_control{1}]));
   area(StimInterval,[maxy,maxy],...
    'FaceColor','k','FaceAlpha',.1,'EdgeAlpha',.1,'HandleVisibility','off');
hold on
   box off;
   
   shadedErrorBar(PSTH_Time_FraOn_All, Bad_FracOn_Loc_Mean_stim{1}, Bad_FracOn_Loc_Sem_stim{1},'lineprops',{[0    0    1]},'transparent',1,'patchSaturation',0.3);
   hold on
   shadedErrorBar(PSTH_Time_FraOn_All, Bad_FracOn_Loc_Mean_control{1}, Bad_FracOn_Loc_Sem_control{1},'lineprops',{[0    0    0]},'transparent',1,'patchSaturation',0.3);
    yrange=ylim;
    plot([0,0],yrange,'--k');






    set(findall(gca, 'Type', 'Line'),'LineWidth',3);
    set(gca,'FontSize',15,'FontWeight','Bold');
    set(gca,'LineWidth',3);

    xlabel('Time from Bad Onset','FontSize',15,'FontWeight','Bold');
    ylabel('FR','FontSize',15,'FontWeight','Bold');
    
      subplot(2,4,6)
      
        % StimInterval 
  maxy = max(max([Bad_SacOn_Loc_Mean_stim{1} + Bad_SacOn_Loc_Sem_stim{1}]),max([Bad_SacOn_Loc_Mean_control{1}+Bad_SacOn_Loc_Sem_control{1}]));
   area( StimIntervalSac,[maxy,maxy],...
    'FaceColor','k','FaceAlpha',.1,'EdgeAlpha',.1,'HandleVisibility','off');
hold on
   box off;
     hold on
     
   
   shadedErrorBar(TimeSequence_SacOn_mean, Bad_SacOn_Loc_Mean_stim{1}, Bad_SacOn_Loc_Sem_stim{1},'lineprops',{[0    0    1]},'transparent',1,'patchSaturation',0.3);
   hold on
   shadedErrorBar(TimeSequence_SacOn_mean, Bad_SacOn_Loc_Mean_control{1}, Bad_SacOn_Loc_Sem_control{1},'lineprops',{[0    0    0]},'transparent',1,'patchSaturation',0.3);
    yrange=ylim;
    plot([0,0],yrange,'--k');


    set(findall(gca, 'Type', 'Line'),'LineWidth',3);
    set(gca,'FontSize',15,'FontWeight','Bold');
    set(gca,'LineWidth',3);

    xlabel('Time from Sac Onset','FontSize',15,'FontWeight','Bold');
    ylabel('FR','FontSize',15,'FontWeight','Bold');
   
   
    
   
   
    subplot(2,4,7)
      % StimInterval 
  maxy = max(max([Bad_FracOn_Loc_Mean_stim{2} + Bad_FracOn_Loc_Sem_stim{2}]),max([Bad_FracOn_Loc_Mean_control{2}+Bad_FracOn_Loc_Sem_control{2}]));
   area(StimInterval,[maxy,maxy],...
    'FaceColor','k','FaceAlpha',.1,'EdgeAlpha',.1,'HandleVisibility','off');
hold on
   box off;
   
   shadedErrorBar(PSTH_Time_FraOn_All, Bad_FracOn_Loc_Mean_stim{2}, Bad_FracOn_Loc_Sem_stim{2},'lineprops',{[0    0    1]},'transparent',1,'patchSaturation',0.3);
   hold on
   shadedErrorBar(PSTH_Time_FraOn_All, Bad_FracOn_Loc_Mean_control{2}, Bad_FracOn_Loc_Sem_control{2},'lineprops',{[0    0    0]},'transparent',1,'patchSaturation',0.3);
    yrange=ylim;
    plot([0,0],yrange,'--k');






    set(findall(gca, 'Type', 'Line'),'LineWidth',3);
    set(gca,'FontSize',15,'FontWeight','Bold');
    set(gca,'LineWidth',3);

    xlabel('Time from Bad Onset','FontSize',15,'FontWeight','Bold');
    ylabel('FR','FontSize',15,'FontWeight','Bold');
   
    
       subplot(2,4,8)
         % StimInterval 
  maxy = max(max([Bad_SacOn_Loc_Mean_stim{2} + Bad_SacOn_Loc_Sem_stim{2}]),max([Bad_SacOn_Loc_Mean_control{2}+Bad_SacOn_Loc_Sem_control{2}]));
   area( StimIntervalSac,[maxy,maxy],...
    'FaceColor','k','FaceAlpha',.1,'EdgeAlpha',.1,'HandleVisibility','off');
hold on
   box off;
     hold on
     
   
   shadedErrorBar(TimeSequence_SacOn_mean, Bad_SacOn_Loc_Mean_stim{2}, Bad_SacOn_Loc_Sem_stim{2},'lineprops',{[0    0    1]},'transparent',1,'patchSaturation',0.3);
   hold on
   shadedErrorBar(TimeSequence_SacOn_mean, Bad_SacOn_Loc_Mean_control{2}, Bad_SacOn_Loc_Sem_control{2},'lineprops',{[0    0    0]},'transparent',1,'patchSaturation',0.3);
    yrange=ylim;
    plot([0,0],yrange,'--k');


    set(findall(gca, 'Type', 'Line'),'LineWidth',3);
    set(gca,'FontSize',15,'FontWeight','Bold');
    set(gca,'LineWidth',3);

    xlabel('Time from Sac Onset','FontSize',15,'FontWeight','Bold');
    ylabel('FR','FontSize',15,'FontWeight','Bold');
%%

    StartNO = StartNO+1;
   FigureIndex = FigureIndex + 1;
   figtitlestr{FigureIndex}='PSTHbetweenGoodBadObjects';
   fig{FigureIndex}=PrepareFigure(StartNO+FigureIndex-1,'w',[70,100, 800,300],'Name',figtitlestr{FigureIndex});
      shadedErrorBar(PSTH_Time_FraOn_All, Good_FracOn_Loc_Mean_control{1}, Good_FracOn_Loc_Sem_control{1},'lineprops',{[1    0    0]},'transparent',1,'patchSaturation',0.3);
      hold on
       shadedErrorBar(PSTH_Time_FraOn_All, Bad_FracOn_Loc_Mean_control{1}, Bad_FracOn_Loc_Sem_control{1},'lineprops',{[0    0    1]},'transparent',1,'patchSaturation',0.3);

   
  
     
 end


end