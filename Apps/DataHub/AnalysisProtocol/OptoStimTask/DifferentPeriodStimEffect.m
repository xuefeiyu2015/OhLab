function DifferentPeriodStimEffect(Data, StartTrial, EndTrial,ShowFigureFlag,OutputFlag);
ECode;

TaskType=Data.TaskType;
DataPath= Data.Path;
FileName=Data.FileName;

Selection=StartTrial:EndTrial;%Select the trials according to the pannel

%Abnormal check for equality
CheckCode=[TGTCD];
AbnormalTrials=AbnormalTrialCheck(Data,CheckCode);

Selection=Selection(~ismember(Selection,AbnormalTrials));

CheckCode=[FIXOFF];
AbnormalTrials=AbnormalTrialCheck(Data,CheckCode);

Selection=Selection(~ismember(Selection,AbnormalTrials));




%Criterias for sexwlecting saccades
V_Threshold=40;
ContinueBin=5;
ISI_Threshold=20;


%GetEyeData
EyeData=Data.EyeDataRaw;
EyeChannel_X_L=EyeData(Selection,1);
EyeChannel_Y_L=EyeData(Selection,2);
%EyeChannel_X_R=EyeData(Selection,3);
%EyeChannel_Y_R=EyeData(Selection,4);
Eye_Eccentricity=cellfun(@(x,y) sqrt(x.^2+y.^2),EyeChannel_X_L,EyeChannel_Y_L,'uniform',0);




EyeBinWidth=Data.EyeBinWidth;
EyeCalInfo=Data.EyeCalInfo;%Each row: X_Gain,X_Bias; Y_Gain; Y_bias; X=(raw+X_bias)*X_Gain
EyeCalInfo_L=EyeCalInfo(1,:);

EyeTime = 0: EyeBinWidth:max(cellfun(@numel,EyeChannel_X_L))-EyeBinWidth;

%Load event data
EventChannel=Data.EventChannel(Selection,:);
EventTimeChannel=Data.EventTimeChannel(Selection,:);
EventBin=Data.EventBin(Selection,:);


%Load stimulation related data

ParaLib=Data.ParaLib;
if ismember('StimType',keys(ParaLib))
    StimTime=ParaLib('StimType');
    StimTime= StimTime(Selection,:);
    StimTime= StimTime(:,2);
    UnqiueStimTime = unique(StimTime);

end

%Seperate into stim trial and sham trial
StimType=ReproduceFromEvent(EventChannel,[STIM1,SHAM0])';

StimTrial = StimType == STIM1;
ControlTrial = ~StimTrial;


%Load target related parameter
TargetLoc=Data.TargetLoc(Selection,2:end);
TargetPolar=Data.TargetPolarLoc(Selection,2:end);
TargetAngle = TargetPolar(:,1);
TargetEcc =  TargetPolar(:,2);

UniqueAngle = unique(TargetAngle);
UniqueEcc = unique(TargetEcc);
minEcc = min(TargetEcc);



%First analysis the behavior: compare the saccade latency or target
%acqusition time between stim and non-stim trials


%Select the period between fix off and target off set
FIXOFFBin = FindOutTime(EventChannel,EventBin,FIXOFF )-500;
SACENDBin = FIXOFFBin + 1500;


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

%Screen out the valid saccades and saccade start time above 5;

AmplitudeThreshold=minEcc/3;

TimeThreshold=0;

SaccadeAngle_Selection=cellfun(@(x,y,z) z(((x>AmplitudeThreshold)&(y>TimeThreshold))),SaccadeAmplitude,SaccadeStartTime,SaccadeAngle,'UniformOutput' ,0);
SaccadeEndPoint_Selection=cellfun(@(x,y,z) z(((x>AmplitudeThreshold)&(y>TimeThreshold)),:),SaccadeAmplitude,SaccadeStartTime,SaccadeEndPoint,'UniformOutput' ,0);
SaccadeStartPoint_Selection=cellfun(@(x,y,z) z(((x>AmplitudeThreshold)&(y>TimeThreshold)),:),SaccadeAmplitude,SaccadeStartTime,SaccadeStartPoint,'UniformOutput' ,0);

SaccadeStartTime_Selection=cellfun(@(x,y,z) z(((x>AmplitudeThreshold)&(y>TimeThreshold)),:),SaccadeAmplitude,SaccadeStartTime,SaccadeStartTime,'UniformOutput' ,0);
SaccadeEndTime_Selection=cellfun(@(x,y,z) z(((x>AmplitudeThreshold)&(y>TimeThreshold)),:),SaccadeAmplitude,SaccadeStartTime,SaccadeEndTime,'UniformOutput' ,0);


SaccadeAmplitude_Selection=cellfun(@(x,y,z) z(((x>AmplitudeThreshold)&(y>TimeThreshold)),:),SaccadeAmplitude,SaccadeStartTime,SaccadeAmplitude,'UniformOutput' ,0);
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


RT_TargetToFixOff = RT_First_Target - 500;

%% For neural responses 
% Align on three period: STIMON, TgtOn, FixOff

StepSize=1;%in ms
BinWidth=10;%in ms

SegEvent=[ON_STIM,TGTCD,FIXOFF];%Event Flow of interest
%SegEvent=[ON_STIM];
%SegEventEnd=[SHOWFIXCD,FIXOFF,TGTOFF,RWDOFFCD ,NaN,EXTRAFP];%Event end
SegEventEnd=[OFF_STIM,NaN,NaN];%Event end
%SegEventEnd=[PST_STIM];%Event end
SegEventStr={'StimOn','TargetOn','FixOff'};
%SegEventStr={'StimOn'};



PlotStartOffset=[300,300,300];%in ms
PlotStopOffset=[200,300,500];%in mss

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
    
    
    %%%%%%%%%%%%%%%%%%%%%%%Align data on each events


SegmentData=EventSegTool(SegEvent,SegEventEnd,PlotStartOffset,PlotStopOffset,EventChannel,EventTimeChannel,SpikeTime,StepSize,BinWidth);

%Retrive data from the field
%PSTH for each event
MeanPSTH=SegmentData.MeanPSTH;%=[PSTH_Event_Time',PSTH_Event'];

PSTH_Time_StimOn = MeanPSTH{1,1};
PSTH_FR_StimOn = MeanPSTH{1,2};

RelativeTimeMarkerStim = SegmentData.RelativeTimeMarker{1,1};

PSTH_Time_TgtOn = MeanPSTH{2,1};
PSTH_FR_TgtOn = MeanPSTH{2,2};


PSTH_Time_FixOff = MeanPSTH{3,1};
PSTH_FR_FixOff = MeanPSTH{3,2};
  
MaxStimOn = 0;
MaxTgtOn = 0;
MaxFixOff = 0;

%Select the raw eye trace 
%Select the eye trace during the 100pre-stim, 100 stim, 100 post stim
%period


StimStartBin=FindOutTime(EventChannel,EventBin,PRE_STIM);
StimEndBin=FindOutTime(EventChannel,EventBin,PST_STIM );

EyeChannel_X_Stim= SelectEyeInterval(EyeChannel_X_L,StimStartBin,StimEndBin)';
EyeChannel_Y_Stim= SelectEyeInterval(EyeChannel_Y_L,StimStartBin,StimEndBin)';

EyeChannel_X_Stim=ReorganizeEye(EyeChannel_X_Stim);
EyeChannel_Y_Stim=ReorganizeEye(EyeChannel_Y_Stim);

EyeTime_Stim = 0:EyeBinWidth:size(EyeChannel_X_Stim,2)-EyeBinWidth;




MarkerTime = FindOutTime(EventChannel,EventBin,ON_STIM)-StimStartBin;

Interval = [-100,nanmean(StimEndBin-FindOutTime(EventChannel,EventBin,ON_STIM))];
%Interval=[-100,200];
StimDur = [0,nanmean(FindOutTime(EventChannel,EventBin,OFF_STIM)- FindOutTime(EventChannel,EventBin,ON_STIM))];

EyeChannel_X_StimAlign = AlignEyeData(EyeChannel_X_Stim,EyeTime_Stim,MarkerTime,Interval);
EyeChannel_Y_StimAlign = AlignEyeData(EyeChannel_Y_Stim,EyeTime_Stim,MarkerTime,Interval);

%NormX = nanmean(EyeChannel_X_StimAlign(:,1:100),2);
%NormY = nanmean(EyeChannel_Y_StimAlign(:,1:100),2);

%EyeChannel_X_StimAlign_Norm = EyeChannel_X_StimAlign-NormX;
%EyeChannel_Y_StimAlign_Norm = EyeChannel_Y_StimAlign-NormY;

EyeTime_Stim_Align=Interval(1):EyeBinWidth:Interval(2)+EyeBinWidth;

%Normalized raw eye trace, separate into stim trials and control trials
EyeX_Stim = EyeChannel_X_StimAlign(StimTrial,:);
EyeY_Stim = EyeChannel_Y_StimAlign(StimTrial,:);

EyeX_Control = EyeChannel_X_StimAlign(ControlTrial,:);
EyeY_Control= EyeChannel_Y_StimAlign(ControlTrial,:);



%% Compare between different conditions
for i = 1:length(unique(StimTime))
    currStim = StimTime == UnqiueStimTime(i);
     %Eye trace 
       EyeX_Stim_loc{i} = EyeChannel_X_StimAlign(currStim & StimTrial,:);
       EyeY_Stim_loc{i} = EyeChannel_Y_StimAlign(currStim & StimTrial,:);
       
       EyeX_Control_loc{i} = EyeChannel_X_StimAlign(currStim & ControlTrial,:);
       EyeY_Control_loc{i} = EyeChannel_Y_StimAlign(currStim & ControlTrial,:);
       

    for j = 1:length(UniqueAngle)
        sel = TargetAngle == UniqueAngle(j) & currStim;
        sel_stim = sel & StimTrial;
        sel_sham = sel & ControlTrial;
        
        RT_stim_loc_period = RT_TargetToFixOff(sel_stim);
        RT_control_loc_period = RT_TargetToFixOff(sel_sham);

        RT_stim_loc_period_Mean(i,j) = nanmean(RT_stim_loc_period);
        RT_stim_loc_period_Sem(i,j) = nanstd(RT_stim_loc_period)/sqrt(sum(~isnan(RT_stim_loc_period)));

        RT_control_loc_period_Mean(i,j) = nanmean(RT_control_loc_period);
        RT_control_loc_period_Sem(i,j) = nanstd(RT_control_loc_period)/sqrt(sum(~isnan(RT_control_loc_period)));
        
      
       
        %StimOn
       PSTH_Time_stim_loc_period = PSTH_Time_StimOn(sel_stim,:);
       PSTH_FR_stim_loc_period = PSTH_FR_StimOn(sel_stim,:);
       
       PSTH_Time_control_loc_period = PSTH_Time_StimOn(sel_sham,:);
       PSTH_FR_control_loc_period = PSTH_FR_StimOn(sel_sham,:);
       
       PSTH_stim_loc_period_Mean{i,j} = nanmean(PSTH_FR_stim_loc_period);
       PSTH_stim_loc_period_Sem{i,j} = nanstd(PSTH_FR_stim_loc_period)./sqrt(sum(~isnan(PSTH_FR_stim_loc_period),1));
       
       PSTH_control_loc_period_Mean{i,j} = nanmean(PSTH_FR_control_loc_period);
       PSTH_control_loc_period_Sem{i,j} = nanstd(PSTH_FR_control_loc_period)./sqrt(sum(~isnan(PSTH_FR_control_loc_period),1));
       
       PSTH_stim_loc_period_MeanTime{i,j} = nanmean(PSTH_Time_stim_loc_period);
       PSTH_control_loc_period_MeanTime{i,j} = nanmean(PSTH_Time_control_loc_period);
       
       if MaxStimOn < max(max(PSTH_stim_loc_period_Mean{i,j} + PSTH_stim_loc_period_Sem{i,j}),max(PSTH_control_loc_period_Mean{i,j} + PSTH_control_loc_period_Sem{i,j}))
          MaxStimOn = max(max(PSTH_stim_loc_period_Mean{i,j} + PSTH_stim_loc_period_Sem{i,j}),max(PSTH_control_loc_period_Mean{i,j} + PSTH_control_loc_period_Sem{i,j}));
       end
       
       Ave_Stim_loc_period_Mean(i,j) = nanmean(PSTH_stim_loc_period_Mean{i,j});
       Ave_Control_loc_period_Mean(i,j) = nanmean(PSTH_control_loc_period_Mean{i,j});
       
       Ave_Stim_loc_period_Sem(i,j) = nanstd(PSTH_stim_loc_period_Mean{i,j})/sum(~isnan(PSTH_stim_loc_period_Mean{i,j}) );
       Ave_Control_loc_period_Sem(i,j) = nanstd(PSTH_control_loc_period_Mean{i,j})/sum(~isnan(PSTH_stim_loc_period_Mean{i,j}) );
       
       %TgtOn
       PSTH_Time_stim_loc_tgtOn = PSTH_Time_TgtOn(sel_stim,:);
       PSTH_FR_stim_loc_tgtOn = PSTH_FR_TgtOn(sel_stim,:);
       
       PSTH_Time_control_loc_tgtOn = PSTH_Time_TgtOn(sel_sham,:);
       PSTH_FR_control_loc_tgtOn =  PSTH_FR_TgtOn(sel_sham,:);
       
       PSTH_stim_loc_tgtOn_Mean{i,j} = nanmean(PSTH_FR_stim_loc_tgtOn );
       PSTH_stim_loc_tgtOn_Sem{i,j} = nanstd(PSTH_FR_stim_loc_tgtOn)./sqrt(sum(~isnan(PSTH_FR_stim_loc_tgtOn),1));
       
       PSTH_control_loc_tgtOn_Mean{i,j} = nanmean(PSTH_FR_control_loc_tgtOn);
       PSTH_control_loc_tgtOn_Sem{i,j} = nanstd(PSTH_FR_control_loc_tgtOn)./sqrt(sum(~isnan(PSTH_FR_control_loc_tgtOn),1));
       
       PSTH_stim_loc_tgtOn_MeanTime{i,j} = nanmean(PSTH_Time_control_loc_tgtOn);
       PSTH_control_loc_tgtOn_MeanTime{i,j} = nanmean(PSTH_Time_stim_loc_tgtOn);
       
       if MaxTgtOn < max(max(PSTH_stim_loc_tgtOn_Mean{i,j} + PSTH_stim_loc_tgtOn_Sem{i,j}),max(PSTH_control_loc_tgtOn_Mean{i,j} + PSTH_control_loc_tgtOn_Sem{i,j}))
          MaxTgtOn = max(max(PSTH_stim_loc_tgtOn_Mean{i,j} + PSTH_stim_loc_tgtOn_Sem{i,j}),max(PSTH_control_loc_tgtOn_Mean{i,j} + PSTH_control_loc_tgtOn_Sem{i,j}));
       end
       
       %FixOff
       PSTH_Time_stim_loc_fixOff = PSTH_Time_FixOff(sel_stim,:);
       PSTH_FR_stim_loc_fixOff = PSTH_FR_FixOff(sel_stim,:);
       
       PSTH_Time_control_loc_fixOff = PSTH_Time_FixOff(sel_sham,:);
       PSTH_FR_control_loc_fixOff =  PSTH_FR_FixOff(sel_sham,:);
       
       PSTH_stim_loc_fixOff_Mean{i,j} = nanmean(PSTH_FR_stim_loc_fixOff );
       PSTH_stim_loc_fixOff_Sem{i,j} = nanstd(PSTH_FR_stim_loc_fixOff )./sqrt(sum(~isnan(PSTH_FR_stim_loc_fixOff),1));
       
       PSTH_control_loc_fixOff_Mean{i,j} = nanmean(PSTH_FR_control_loc_fixOff);
       PSTH_control_loc_fixOff_Sem{i,j} = nanstd(PSTH_FR_control_loc_fixOff )./sqrt(sum(~isnan(PSTH_FR_control_loc_fixOff),1));
       
       PSTH_stim_loc_fixOff_MeanTime{i,j} = nanmean(PSTH_Time_control_loc_fixOff);
       PSTH_control_loc_fixOff_MeanTime{i,j} = nanmean(PSTH_Time_stim_loc_fixOff);
       
       if MaxFixOff < max(max(PSTH_stim_loc_fixOff_Mean{i,j} + PSTH_stim_loc_fixOff_Sem{i,j}),max(PSTH_control_loc_fixOff_Mean{i,j} + PSTH_control_loc_fixOff_Sem{i,j}))
          MaxFixOff = max(max(PSTH_stim_loc_fixOff_Mean{i,j} + PSTH_stim_loc_fixOff_Sem{i,j}),max(PSTH_control_loc_fixOff_Mean{i,j} + PSTH_control_loc_fixOff_Sem{i,j}));
       end
       
       

 
        
    end



end



%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% For figure plotting%%%%%%%%%%%%%%%%%
StrStimPeriod = {'DurFix','TgtOn','FixOff','SacEnd'};

if ShowFigureFlag
   StartNO = 50;
   FigureIndex=1;
   figtitlestr{FigureIndex}='RTDifference';
   fig{FigureIndex}=PrepareFigure(StartNO+FigureIndex-1,'w',[50,100, 1200,800],'Name',figtitlestr{FigureIndex});
    for i = 1:length(unique(StimTime))

      subplot(length(UnqiueStimTime),3,(i-1)*3+1);
      %Polar plot for each stim time
      %Control
        polarwitherrorbar([UniqueAngle/180*pi;UniqueAngle(1)/180*pi],[RT_control_loc_period_Mean(i,:) RT_control_loc_period_Mean(i,1)]',[RT_control_loc_period_Sem(i,:) RT_control_loc_period_Sem(i,1)]',...
            '-b',3);
        hold on
      %Stim
       polarwitherrorbar([UniqueAngle/180*pi;UniqueAngle(1)/180*pi],[RT_stim_loc_period_Mean(i,:) RT_stim_loc_period_Mean(i,1)]',[RT_stim_loc_period_Sem(i,:) RT_stim_loc_period_Sem(i,1)]',...
            '-r',3);

    set(gca,'ThetaDir' , 'counterclockwise');
    set(gca,'ThetaZeroLocation','top')
    set(gca,'FontSize',25,'FontWeight','Bold','LineWidth',3);

    thetaticks([0:45:359]);
    thetaticklabels([0:45:359]);
    title(StrStimPeriod(UnqiueStimTime(i)));
    
    subplot(length(UnqiueStimTime),3,(i-1)*3+2);
    %Raw eye trace
    plot(EyeX_Stim_loc{i}',EyeY_Stim_loc{i}','-r','LineWidth',2);
     xlim([-20,20]);
    ylim([-20,20]);
    axis equal;
    
    set(gca,'LineWidth',3,'FontSize',20,'FontWeight','Bold');
    
    subplot(length(UnqiueStimTime),3,(i-1)*3+3);
    plot(EyeX_Control_loc{i}',EyeY_Control_loc{i}','-b','LineWidth',2);
    
    xlim([-20,20]);
    ylim([-20,20]);
    axis equal;
    set(gca,'LineWidth',3,'FontSize',20,'FontWeight','Bold');


        
    end
    
    %% Figure2: PSTH for different period
    
    
    for i = 1:length(unique(StimTime))
        StartNO = StartNO+1;
        FigureIndex = FigureIndex+1;
    
        figtitlestr{FigureIndex}=sprintf('PSTHDifference_Period_%s',StrStimPeriod{UnqiueStimTime(i)});
        fig{FigureIndex}=PrepareFigure(StartNO+FigureIndex-1,'w',[70+10*i,100, 600,600],'Name',figtitlestr{FigureIndex});
    
        for j = 1:length(UniqueAngle)
            %First stim on 
            subplot(length(UniqueAngle),3,j*3-2);
            AverageRegion = nanmean(RelativeTimeMarkerStim);
            area(AverageRegion,[MaxStimOn,MaxStimOn],...
             'FaceColor','k','FaceAlpha',.1,'EdgeAlpha',.1,'HandleVisibility','off');
            hold on
            
            shadedErrorBar(PSTH_stim_loc_period_MeanTime{i,j} ,PSTH_stim_loc_period_Mean{i,j},PSTH_stim_loc_period_Sem{i,j},'lineprops',{[1,0,0]},'transparent',1,'patchSaturation',0.3)
            hold on
            shadedErrorBar(PSTH_control_loc_period_MeanTime{i,j} ,PSTH_control_loc_period_Mean{i,j},PSTH_control_loc_period_Sem{i,j},'lineprops',{[0,0,1]},'transparent',1,'patchSaturation',0.3)
             
           
            ylim([0,MaxStimOn]);
            
            set(findall(gca, 'Type', 'Line'),'LineWidth',2);
            
            if j == length(UniqueAngle)
                xlabel('Time from stim onset');
                
            end
            ylabel(sprintf('Angle:%d',UniqueAngle(j)));
            box off;
            

              %  ylabel('Firing Rate(Hz)')
            set(gca,'LineWidth',2,'FontSize',10,'FontWeight','Bold');
           % title(sprintf("Laser Power : %d mW",unique(stimPower)));
           % legend(["Sham","Stim"])
           % box off
            
            %Second target on
            subplot(length(UniqueAngle),3,j*3-1);
           
            
            shadedErrorBar(PSTH_stim_loc_tgtOn_MeanTime{i,j} ,PSTH_stim_loc_tgtOn_Mean{i,j},PSTH_stim_loc_tgtOn_Sem{i,j},'lineprops',{[1,0,0]},'transparent',1,'patchSaturation',0.3)
            hold on
            shadedErrorBar(PSTH_control_loc_tgtOn_MeanTime{i,j} ,PSTH_control_loc_tgtOn_Mean{i,j},PSTH_control_loc_tgtOn_Sem{i,j},'lineprops',{[0,0,1]},'transparent',1,'patchSaturation',0.3)
      
            yrange = ylim;
            ylim([0,MaxTgtOn]);
            set(findall(gca, 'Type', 'Line'),'LineWidth',2);
            if j == length(UniqueAngle)
                xlabel('Time from tgt onset');
                
            end
             set(gca,'LineWidth',2,'FontSize',10,'FontWeight','Bold');
            
            %Third fix off 
            subplot(length(UniqueAngle),3,j*3);
            
            shadedErrorBar(PSTH_stim_loc_fixOff_MeanTime{i,j} ,PSTH_stim_loc_fixOff_Mean{i,j},PSTH_stim_loc_fixOff_Sem{i,j},'lineprops',{[1,0,0]},'transparent',1,'patchSaturation',0.3)
            hold on
            shadedErrorBar(PSTH_control_loc_fixOff_MeanTime{i,j} ,PSTH_control_loc_fixOff_Mean{i,j},PSTH_control_loc_fixOff_Sem{i,j},'lineprops',{[0,0,1]},'transparent',1,'patchSaturation',0.3)
      
            yrange = ylim;
            ylim([0,MaxFixOff]);
            set(findall(gca, 'Type', 'Line'),'LineWidth',2);
             if j == length(UniqueAngle)
                xlabel('Time from fix offset');
                
            end
            
         set(gca,'LineWidth',2,'FontSize',10,'FontWeight','Bold');
        
        end
    end
    
    StartNO = StartNO+1;
    FigureIndex = FigureIndex+1;
        
  %% Figure3: Polar plot for stimulation under difference angle condition
   figtitlestr{FigureIndex}='ActivationDifference';
   fig{FigureIndex}=PrepareFigure(StartNO+FigureIndex-1,'w',[90,100, 800,400],'Name',figtitlestr{FigureIndex});
    for i = 1:length(unique(StimTime))

      subplot(1,length(UnqiueStimTime),i);
      %Polar plot for each stim time
      %Control
        polarwitherrorbar([UniqueAngle/180*pi;UniqueAngle(1)/180*pi],[Ave_Control_loc_period_Mean(i,:) Ave_Control_loc_period_Mean(i,1)]',[Ave_Control_loc_period_Sem(i,:) Ave_Control_loc_period_Sem(i,1)]',...
            '-b',3);
        hold on
      %Stim
       polarwitherrorbar([UniqueAngle/180*pi;UniqueAngle(1)/180*pi],[Ave_Stim_loc_period_Mean(i,:) Ave_Stim_loc_period_Mean(i,1)]',[Ave_Stim_loc_period_Sem(i,:) Ave_Stim_loc_period_Sem(i,1)]',...
            '-r',3);

    set(gca,'ThetaDir' , 'counterclockwise');
    set(gca,'ThetaZeroLocation','top')
    set(gca,'FontSize',25,'FontWeight','Bold','LineWidth',3);

    thetaticks([0:45:359]);
    thetaticklabels([0:45:359]);
    title(StrStimPeriod(UnqiueStimTime(i)));


        
    end
    


end %End of the if show figure









end %End of the spikechannel


end

function ObjectIndex=ReproduceFromEvent(event,code)
for i=1:size(event,1)
    
    ObjectIndex{i}=event(i,ismember(event(i,:),code));
    
  
    
end
numeach=cellfun(@numel,ObjectIndex);
nummax=max(numeach);
if nummax==1
    ObjectIndex=cell2mat(ObjectIndex);

    
end


end