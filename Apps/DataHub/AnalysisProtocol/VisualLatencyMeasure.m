function VisualLatencyMeasure(Data, StartTrial, EndTrial,ShowFigureFlag,OutputFlag);

ECode;
SegEvent=[TGTCD];%Event Flow of interest
%SegEventEnd=[SHOWFIXCD,FIXOFF,TGTOFF,RWDOFFCD ,NaN,EXTRAFP];%Event end
SegEventEnd=[TGTOFF];%Event end
SegEventStr={'TargetOn'};
AmpMatchFlag=0;%1:Amplitude match on ;0:Amplitude Match off


PlotStartOffset=[300];%in ms
PlotStopOffset=[200];%in ms

ExamplePlot=0;%PLots for presentations, turn off for regular use

ResponseLatencyThreshold=3;

%Parameter to define the sliding window to measure firing rates
StepSize=5;%in ms
BinWidth=15;%in ms


%Criterias for selecting saccades
V_Threshold=40;
ContinueBin=5;
ISI_Threshold=20;

%Time bin pre saccade and time bin offset after saccade onset
PlotSacPreBuffer=500;%in ms
PlotSacPostBuffer=200;%in ms



%Color Definition for each fractals
%{
Color=[lines(6);prism(6) ];
Color(2,:)=copper(1)+0.5;
Color=[Color;Color*0.5];
%}
%{
Color={'#648FFF','#FE6100','#009E73','#88CCEE','#CC6677','#785EF0'};
Color_RGB=[100,143,255;
           254,97,0;
           0,158,115;
          136,204,238;
          204,102,119;
          120,94,240]/255;

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Unpack all the data base and task settings
ECode;
FileName=Data.FileName;
TaskType=Data.TaskType;
TaskCode=Data.TaskCode;
TrialNum=EndTrial-StartTrial+1;
Selection=StartTrial:EndTrial;


%Abnormal check
CheckCodePair=[TGTCD,TGTOFF];%This code should be the same as target on for good trials
AbnormalTrials=AbnormalTrialCheck(Data,CheckCodePair);
Selection=Selection(~ismember(Selection,AbnormalTrials));
    



%Target Information


%% Target location

TargetLoc=Data.TargetLoc(Selection,2:end);
TargetPolar=Data.TargetPolarLoc(Selection,2:end);
TargetAngle = TargetPolar(:,1);
TargetEcc =  TargetPolar(:,2);

UniqueAngle = unique(TargetAngle);
UniqueEcc = unique(TargetEcc);
minEcc = min(TargetEcc);





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


%Select the period between fix off and target off set
FIXOFFBin = FindOutTime(EventChannel,EventBin,FIXOFF )-500;
SACENDBin = FIXOFFBin + 1500;
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

%Align neural responses according to the saccade onset time
%

%Organize the spike time using the same way
%Align spike time with target onset

AlignedSpike=SpikeAlignTool(FIXOFFBin,SpikeTime,StepSize,BinWidth);

SpikeTimeAligned_FixOff=AlignedSpike.SpikeTime_Aligned_Whole;




%Align again on saccade onset
AlignedSpike=SpikeAlignTool(RT_First_Target',SpikeTimeAligned_FixOff,StepSize,BinWidth);

Raster_SacOn_whole = AlignedSpike.SpikeTime_Aligned_Whole;

PSTH_SacOn=AlignedSpike.PSTH_Aligned;
TimeSequence_SacOn=AlignedSpike.PSTH_Time_Aligned;
SpikeCount_SacOn = AlignedSpike.PSTH_SpikeCount_Aligned/BinWidth*1000;
%Select interval
SelectPSTHinterval_SacOn=[-PlotSacPreBuffer,PlotSacPostBuffer];
[PSTH_SacOn,TimeSequence_SacOn]=SelectPSTHInterval(PSTH_SacOn,TimeSequence_SacOn,SelectPSTHinterval_SacOn);

[SpikeCount_SacOn,TimeSequence_SacOn]=SelectPSTHInterval(SpikeCount_SacOn,TimeSequence_SacOn,SelectPSTHinterval_SacOn);




TimeSequence_SacOn_mean=nanmean(TimeSequence_SacOn,1);
PSTH_SacOn_mean=nanmean(PSTH_SacOn,1);
PSTH_SacOn_sem=nanstd(PSTH_SacOn,[],1)./sqrt(sum(~isnan(PSTH_SacOn),1));

%Change raster
Raster_SacOn= SelectSpikeTimeInterval(Raster_SacOn_whole,SelectPSTHinterval_SacOn(1)*ones(size(Raster_SacOn_whole,1),1),SelectPSTHinterval_SacOn(end)*ones(size(Raster_SacOn_whole,1),1));


%Align Eye trace with saccade onset
%Select the eye eccentricity between target onset and target offset
Eye_Eccentricity_FixOff=SelectEyeInterval(Eye_Eccentricity,FIXOFFBin,SACENDBin)';
Eye_Eccentricity_FixOff=ReorganizeEye(Eye_Eccentricity_FixOff);



%Select the interval around saccade onset
EyeTime=0:EyeBinWidth:size(Eye_Eccentricity_FixOff,2)-1;

EyeEcc_Aligned=AlignEyeData(Eye_Eccentricity_FixOff,EyeTime,RT_First_Target,[-PlotSacPreBuffer,PlotSacPostBuffer]);



EyeEcc_Aligned_mean=nanmean(EyeEcc_Aligned,1);
EyeEcc_Aligned_sem=nanstd(EyeEcc_Aligned,[],1)./sqrt(sum(~isnan(EyeEcc_Aligned),1));



%%
%%%%%%%%%%%%%%%%%%%%%%%Align data on object onset
%Segregate spike data from fractal onset to fractal offset

SegmentData=EventSegTool(SegEvent,SegEventEnd,PlotStartOffset,PlotStopOffset,EventChannel,EventTimeChannel,SpikeTime,StepSize,BinWidth);

%Retrive data from the field
%Overall firing rate for z-score
FR_Each = SegmentData.MeanFRRaw{1};
FR_Mean = mean(FR_Each);
FR_Std = std(FR_Each);

%Relative time bin for each event
RelativeTimeMarker=SegmentData.RelativeTimeMarker{1};%=EventTimeMarker;
%Absolute time marker for each event
AbosoluteTimeMarker=SegmentData.AbsoluteTimeMarker{1};

ObjOnTimeMarker=nanmean(RelativeTimeMarker);

%% Spike Count 
SpikeCount = SegmentData.PSTH_Event_SpikeCount{1}/BinWidth*1000;

%% Raster after align
RasterTgtOn = SegmentData.RasterAfterAlign{1};

%%

RwdTimeAbs=FindOutTime(EventChannel,EventBin,RWDCD);
RwdTimeRel=RwdTimeAbs-AbosoluteTimeMarker(:,1);

%Algin marker to sac-on & sac-off
RwdTimeRelSacEnd=RwdTimeRel-SacEnd_First_Target';
TgtOffTimeRelSacEnd=RelativeTimeMarker(:,2)-SacEnd_First_Target';

RwdTimeRelSacEndMean=nanmean(RwdTimeRelSacEnd);
TgtOffTimeRelSacEndMean=nanmean(TgtOffTimeRelSacEnd);

%PSTH for each event
MeanPSTH=SegmentData.MeanPSTH;%=[PSTH_Event_Time',PSTH_Event'];
MeanPSTH_Time=MeanPSTH{1};
MeanPSTH_FR_FractalOn=MeanPSTH{2};

PSTH_Time_FraOn_All=nanmean(MeanPSTH_Time,1);

%Get mean and sem of neural responses
MeanPSTH_FR_FractalOn_mean=nanmean(MeanPSTH_FR_FractalOn,1);
MeanPSTH_FR_FractalOn_sem=nanstd(MeanPSTH_FR_FractalOn,[],1)./sqrt(sum(~isnan(MeanPSTH_FR_FractalOn),1));



%Get mean and sem of eye eccentricity
%Get eye trance distance during the interval
StartTimeBin=FindOutTime(EventChannel,EventTimeChannel,SegEvent)-PlotStartOffset;
StartTimeBin(StartTimeBin<1)=1;
EndTimeBin=FindOutTime(EventChannel,EventTimeChannel,SegEventEnd )+PlotStopOffset;
  
Eye_Eccentricity_FractalOn=SelectEyeInterval(Eye_Eccentricity,StartTimeBin,EndTimeBin)';
Eye_Eccentricity_FractalOn=ReorganizeEye(Eye_Eccentricity_FractalOn);



Eye_Eccentricity_FractalOn_mean=nanmean(Eye_Eccentricity_FractalOn,1);
Eye_Eccentricity_FractalOn_sem=nanstd(Eye_Eccentricity_FractalOn,1)./sqrt(sum(~isnan(Eye_Eccentricity_FractalOn),1));


BaseLineWindow = [-300, -100];
VisualEpoch = [90,180];
MovementEpoch = [-80,0];


FR_BaseLine = mean(MeanPSTH_FR_FractalOn(:,PSTH_Time_FraOn_All<=BaseLineWindow(2) & PSTH_Time_FraOn_All>=BaseLineWindow(1)),2);
FR_VisualEpoch = mean(MeanPSTH_FR_FractalOn(:,PSTH_Time_FraOn_All<=VisualEpoch(2) & PSTH_Time_FraOn_All>=VisualEpoch(1)),2);
FR_MovementEpoch = mean(PSTH_SacOn(:,TimeSequence_SacOn_mean<=MovementEpoch(2) & TimeSequence_SacOn_mean>=MovementEpoch(1)),2);




%% Loop over the target directions
for i = 1:length(UniqueAngle )
    for j = 1:length(UniqueEcc)
        sel = TargetAngle == UniqueAngle(i) & TargetEcc == UniqueEcc(j);
        
        
        
        %Align on target on 
        currPSTH_FR_TgtOn_Loc = MeanPSTH_FR_FractalOn(sel,:);
        TgtOn_Loc_Mean{i,j} = nanmean(currPSTH_FR_TgtOn_Loc,1);
        TgtOn_Loc_Sem{i,j} = nanstd(currPSTH_FR_TgtOn_Loc,[],1)./sqrt(sum(~isnan(currPSTH_FR_TgtOn_Loc),1));
        
        
        
       
        
        AngleEccorder{i,j} = [UniqueAngle(i),UniqueEcc(j)];

        
         %Align on sac on 
         
         curr_SacOn_Loc = PSTH_SacOn(sel,:);
         SacOn_Loc_Mean{i,j} = nanmean(curr_SacOn_Loc ,1);
         SacOn_Loc_Sem{i,j} = nanstd(curr_SacOn_Loc ,[],1)./sqrt(sum(~isnan(curr_SacOn_Loc ),1));

         %Mean
         Mean_Visual(i,j) = nanmean(FR_VisualEpoch(sel));
         Mean_Movement(i,j) = nanmean(FR_MovementEpoch(sel));

         %Raster
         Raster_TgtOn_ij{i,j} = RasterTgtOn(sel,:); 
         Raster_SacOn_ij{i,j} = Raster_SacOn(sel,:); 


         %Cell identification 
         
         [h,p] = ttest(FR_BaseLine(sel),FR_VisualEpoch(sel));
         Diff_Visual = mean(FR_VisualEpoch(sel) - FR_BaseLine(sel));
         if p<0.05
             h_Visual(i,j) = sign(Diff_Visual);
         else
             h_Visual(i,j) = 0;

         end

         [h,p] = ttest(FR_BaseLine(sel),FR_MovementEpoch(sel));
         Diff_Movement = mean(FR_MovementEpoch(sel) - FR_BaseLine(sel));

         if p<0.05
             h_Movement(i,j) = sign(Diff_Movement);
         else
             h_Movement(i,j) = 0;
         end

    end
end


if length(UniqueAngle) > 1
    RF_Angle = NaN;
    %Choose the angle(ignore the eccentricity)
    h_VisOrMov = h_Visual(:,1) ~=0|h_Movement(:,1) ~=0;
    h_VisOrMov_Sig = h_VisOrMov(h_VisOrMov==1);
    h_v_sig = h_Visual(h_VisOrMov==1);
    h_m_sig = h_Movement(h_VisOrMov==1);

    AngleSig = UniqueAngle(h_VisOrMov);
    FRVisSig = Mean_Visual(h_VisOrMov,1);
    FRMovementSig = Mean_Movement(h_VisOrMov,1);

    if length(AngleSig) >0
        [MaxVisual, index_v ]= max(abs(FRVisSig - mean(FR_BaseLine)));
        MaxVisual_Angle = AngleSig(index_v);
        h_vis_max = h_VisOrMov_Sig(index_v);

        [MaxMovement, index_m]= max(abs(FRMovementSig - mean(FR_BaseLine)));
        MaxMovement_Angle = AngleSig(index_m);
        h_move_max = h_VisOrMov_Sig(index_m);

        if h_vis_max ~=0
            RF_Angle = MaxVisual_Angle;
            RF_h_v = h_v_sig(index_v);
            RF_h_m = h_m_sig(index_v);

        else
            RF_Angle = MaxMovement_Angle;
            RF_h_v = h_v_sig(index_m);
            RF_h_m = h_m_sig(index_m);

        end



    end

    if isnan(RF_Angle)
        [MaxVisual, index_v]= max(abs(mean(Mean_Visual(:,1)) - mean(FR_BaseLine)));
         MaxVisual_Angle = UniqueAngle(index_v);
        h_vis_max = h_VisOrMov(index_v);
        RF_Angle = MaxVisual_Angle;
        RF_h_v =h_Visual(index_v);
        
        RF_h_m =h_Movement(index_v);




    end

    
else
    RF_Angle = UniqueAngle;
    RF_h_v =h_Visual;
    RF_h_m =h_Movement;
    



end


PSTH_RF_Mean = TgtOn_Loc_Mean{UniqueAngle==RF_Angle};
PSTH_RF_Sem = TgtOn_Loc_Sem{UniqueAngle==RF_Angle};

PSTH_RF_Sac_Mean = SacOn_Loc_Mean{UniqueAngle==RF_Angle};
PSTH_RF_Sac_Sem = SacOn_Loc_Sem{UniqueAngle==RF_Angle};


PSTH_RF_Mean_z = (PSTH_RF_Mean - FR_Mean)/FR_Std;
PSTH_RF_Sac_Mean_z= (PSTH_RF_Sac_Mean - FR_Mean)/FR_Std;


Raster_TgtOn_RF  = Raster_TgtOn_ij{UniqueAngle==RF_Angle,1} ;
Raster_SacOn_RF  = Raster_SacOn_ij{UniqueAngle==RF_Angle,1} ;


CellType = 'None';
CellTypeCode = 0;
if RF_h_v ~=0 & RF_h_m ==0
    CellType = 'Visual-only';
    CellTypeCode = 1;

end

if RF_h_v ==0 & RF_h_m ~=0
    CellType = 'Sac-Only';
    CellTypeCode = 2;

end

if RF_h_v ~=0 & RF_h_m ~=0
    CellType = 'Vis-Sac';
    CellTypeCode = 3;

end

%Use Leor's method to calculate the latency
EventTime_All=FindOutTime(EventChannel,EventTimeChannel,SegEvent);
EventTime = EventTime_All(TargetAngle == RF_Angle);

SpikeTime= SpikeTime(TargetAngle == RF_Angle,1:end-2);

Heading = [1:length(EventTime)]'*10000;

SpikeTime_H = reshape(Heading+SpikeTime/1000,1,numel(SpikeTime));
SpikeTime_H = SpikeTime_H(~isnan(SpikeTime_H));
EventTime_H = Heading +EventTime/1000;%in s


methodList  = {'white2017', 'yuKatz2023'};
vizBool     = false;

for iM= 1:length(methodList)
latencyTime(iM) = compute_latency(SpikeTime_H, EventTime_H, methodList{iM}, vizBool, [ ', Neuron #' num2str(ChannelID) ' - ' methodList{iM}])*1000;
end

disp(latencyTime)


%%

if ShowFigureFlag
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Figure 1: Average firing rate
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%for each location

FigIndex = 110;
for i = 1:length(UniqueEcc)
    figtitlestr{i}=sprintf('PSTH_Ecc %d,CellType:%s',UniqueEcc(i),CellType);
    fig{i}=PrepareFigure(FigIndex + i-1,'w',[50,100, 1200,800],'Name',figtitlestr{i});
    
    
    
    for j = 1:length(UniqueAngle)
        
          subplot(length(UniqueAngle),2,j*2-1);
          %Align on tgtonset
          shadedErrorBar(PSTH_Time_FraOn_All,TgtOn_Loc_Mean{j,i},TgtOn_Loc_Sem{j,i},'lineprops',{[1    0    0]},'transparent',1,'patchSaturation',0.3);
           maxyTgt = max(TgtOn_Loc_Mean{j,i}+TgtOn_Loc_Sem{j,i});
         maxySac = max(SacOn_Loc_Mean{j,i}+SacOn_Loc_Sem{j,i});
         if ~isnan(maxyTgt)
            ylim([0,max(max(maxyTgt, maxySac))]);
         
            ymax = max(max(maxyTgt, maxySac));
          hold on
          plot([0,0],[0,ymax],'--k','LineWidth',3);
          end
          
            set(findall(gca, 'Type', 'Line'),'LineWidth',3);

        set(gca,'FontSize',15,'FontWeight','Bold');
        set(gca,'LineWidth',3);
        title(sprintf("Angle: %1.1d,Ecc: %1.1d",AngleEccorder{j,i}));
       
            xlabel('Time from Tgt Onset(ms)','FontSize',15,'FontWeight','Bold');
      
        
        ylabel('Firing Rate(Hz)','FontSize',15,'FontWeight','Bold');
        
       
          
          %Align on saccade onset
          subplot(length(UniqueAngle),2,j*2);
         shadedErrorBar(TimeSequence_SacOn_mean,SacOn_Loc_Mean{j,i},SacOn_Loc_Sem{j,i},'lineprops',{[0    0    1]},'transparent',1,'patchSaturation',0.3);
         hold on
          set(findall(gca, 'Type', 'Line'),'LineWidth',3);

        set(gca,'FontSize',15,'FontWeight','Bold');
        set(gca,'LineWidth',3);
        title(sprintf("Angle: %1.1d,Ecc: %1.1d",AngleEccorder{j,i}));
       
            xlabel('Time from Sac Onset(ms)','FontSize',15,'FontWeight','Bold');
            
        
        
        ylabel('Firing Rate(Hz)','FontSize',15,'FontWeight','Bold');
        
        
         
         if ~isnan(maxyTgt)
            ylim([0,max(max(maxyTgt, maxySac))]);
         
          plot([0,0],[0,ymax],'--k','LineWidth',3);
          end
        
        
        
    end
    
    %Add title

title_msg=sprintf(' %s, SpikeChannel: %d, TrialNum: %d ',FileName,SpikeChannel,TrialNum);
a=annotation('textbox',[0.25,0.7,0.3,0.3],'String',title_msg,'FitBoxToText','on');
a.FontSize=25;
a.FontWeight='Bold';
a.LineStyle='none';
    
    
    
    
    
    
end
%% Visual-Motor responses in RF

figtitlestr{3}=sprintf('RF;Angle %d, %s',RF_Angle,CellType);
    fig{3}=PrepareFigure(FigIndex + 3-1,'w',[100,100, 600,400],'Name',figtitlestr{3});

    StartLoc = 1;
LineWidth = 2;
Barlength = 1;

 subplot(2,2,1);
 SpikeRaterPlot(Raster_TgtOn_RF,StartLoc,Barlength,'r',LineWidth);
 yrange = ylim;
 hold on
 plot([0,0],yrange,'--k','LineWidth',3);
 xlim([-PlotStartOffset,300]);

    subplot(2,2,2);
    shadedErrorBar(PSTH_Time_FraOn_All,TgtOn_Loc_Mean{UniqueAngle==RF_Angle},TgtOn_Loc_Sem{UniqueAngle==RF_Angle},'lineprops',{[1    0    0]},'transparent',1,'patchSaturation',0.3);
     set(findall(gca, 'Type', 'Line'),'LineWidth',3);

        set(gca,'FontSize',15,'FontWeight','Bold');
        set(gca,'LineWidth',3);
        xlim([-PlotStartOffset,300]);

        subplot(2,2,3);
 SpikeRaterPlot(Raster_SacOn_RF ,StartLoc,Barlength,'b',LineWidth);
 hold on
 yrange = ylim;
 plot([0,0],yrange,'--k','LineWidth',3);
 xlim([-PlotSacPreBuffer,PlotSacPostBuffer]);


    subplot(2,2,4);
     shadedErrorBar(TimeSequence_SacOn_mean,SacOn_Loc_Mean{UniqueAngle==RF_Angle},SacOn_Loc_Sem{UniqueAngle==RF_Angle},'lineprops',{[0    0    1]},'transparent',1,'patchSaturation',0.3);
      set(findall(gca, 'Type', 'Line'),'LineWidth',3);

        set(gca,'FontSize',15,'FontWeight','Bold');
        set(gca,'LineWidth',3);  
        xlim([-PlotSacPreBuffer,PlotSacPostBuffer]);
    


    


end %End of ShowFigureFlag

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

%Store figures in common folders
%{
ProtocolMarker = 'MemorySac';
FigurePath=strcat(BasicDirectory,'Figures',BasicDirectory(end),ProtocolMarker,BasicDirectory(end));%,OutputFolerName);
if ~exist(FigurePath)
    mkdir(FigurePath)
end

if ShowFigureFlag & OutputFlag
    %for i=1:length(fig)
    for i=1:1
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

%}

%


ExistFlag=0;
%Load the old file if exist
if exist(OutputFileName)

    load(OutputFileName);
    %Load OutputData into the memory
    if isfield(OutputData,'MemorySaccadeLatency')
    ExistFlag=1;
    end
end


%%%%%%%%%%%%%%%%%%%%OutputData For current analysis
OutputData_New=[];

OutputData_New.MemorySaccadeLatency.TrialType='GoodTrials';
%%%%%%%%%Save overall scene response characteristics%%%%%%%%%%%%%%%%%%%%%%%
%Organize the data
NameStr={'PSTH_Time_TgtOn','PSTH_TgtOn_RF',...
    'PSTH_Time_SacOn','PSTH_SacOn_RF',...
    'RF_h_v','RF_h_m','CellType',...
    'PSTH_TgtOn_RF_z','PSTH_SacOn_RF_z',...
    'latencyTime'
    
    
};

DataLib={PSTH_Time_FraOn_All,PSTH_RF_Mean,...
    TimeSequence_SacOn_mean,PSTH_RF_Sac_Mean,...
    RF_h_v,RF_h_m,CellTypeCode,...
    PSTH_RF_Mean_z,PSTH_RF_Sac_Mean_z,...
    latencyTime

};


DataStamp=containers.Map(NameStr,DataLib);





OutputData_New.MemorySaccadeLatency.Task=TaskType;
OutputData_New.MemorySaccadeLatency.TaskCode=TaskCode;
OutputData_New.MemorySaccadeLatency.DataStamp=DataStamp;


%Store figures
if ShowFigureFlag
   
    %OutputData_New.SceneTuning.Figures={fig99,fig100,fig101,fig102};
    %Output figures in a folder;
    FolderName=OutputFigureName;
    if ~exist(FolderName,'dir')
        mkdir(FolderName);
    end
    cd(FolderName);
    
    TotalFigure=length(fig);
    
    for j=1:TotalFigure
        OutputFigureName=strcat(OutputFigureName,figtitlestr{j},'.jpg');
        if ~isempty(fig{j})
        saveas(fig{j},OutputFigureName);
        end
        
    end

    disp('Figures have been exported to the neuron folder');

end

if ExistFlag
%Compare the old one with the new one
Task_Old=OutputData.MemorySaccadeLatency.TaskCode;

     if Task_Old~=TaskCode
        %Not the same task, add the new dataset into the old one
        OutputData.MemorySaccadeLatency(numel(OutputData.MemorySaccadeLatency)+1)=OutputData_New.MemorySaccadeLatency;
     else
         %If the same task,the same protocol,replace the old one with the new one
         
         OutputData.MemorySaccadeLatency=OutputData_New.MemorySaccadeLatency;
         
    
    
     end
else
    %Output the current file
   OutputData.MemorySaccadeLatency=OutputData_New.MemorySaccadeLatency; 
end



if OutputFlag
    cd(OutPath);
    save(OutputFileName,'OutputData');
    disp('Data Exported');
else
    disp('No data export');
end  

%}





end %End of the spike loop
end