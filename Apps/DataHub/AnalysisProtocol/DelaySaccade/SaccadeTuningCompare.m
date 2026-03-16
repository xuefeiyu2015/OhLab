function SaccadeTuningCompare(Data, StartTrial, EndTrial,ShowFigureFlag,OutputFlag);
% Comparison of the goal directed saccade during delay saccade task and
% non-goal directed saccade durign intertrial interval
% Xuefei Yu 12122022
ECode;
SegEvent=[TGTCD];%Event Flow of interest
%SegEventEnd=[SHOWFIXCD,FIXOFF,TGTOFF,RWDOFFCD ,NaN,EXTRAFP];%Event end
SegEventEnd=[TGTOFF];%Event end
SegEventStr={'TargetOn'};
AmpMatchFlag=0;%1:Amplitude match on ;0:Amplitude Match off


PlotStartOffset=[800];%in ms
PlotStopOffset=[1000];%in ms

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
PlotSacPostBuffer=400;%in ms



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

VaryEcc = 0;
MultiLoop =0;
if length(UniqueEcc) > 1
    %Use the most frequent ecc for angle  loop
    for i = 1:length(UniqueEcc)
        NumEcc(i) = sum(TargetEcc == UniqueEcc(i));
    end
    useEcc = UniqueEcc(NumEcc==max(NumEcc));
    if length(useEcc) > 1

        MultiLoop = 1;


    end
    VaryEcc=1;
    
    
else
    useEcc =UniqueEcc;
end

if FileName== 'Robin072123DS1.rst'
     useEcc =15;
end




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

%Load stimulation related data

ParaLib=Data.ParaLib;
if ismember('StimType',keys(ParaLib))
    StimTime=ParaLib('StimType');
    StimTime= StimTime(Selection,:);
    StimTime= StimTime(:,2);
    UnqiueStimTime = unique(StimTime);

    %Seperate into stim trial and sham trial
StimType=ReproduceFromEvent(EventChannel,[STIM1,SHAM0])';

StimTrial = StimType == STIM1;
ControlTrial = ~StimTrial;
else
    ControlTrial = [1:length(Selection)]';

end




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
%Select interval
SelectPSTHinterval_SacOn=[-PlotSacPreBuffer,PlotSacPostBuffer];
[PSTH_SacOn,TimeSequence_SacOn]=SelectPSTHInterval(PSTH_SacOn,TimeSequence_SacOn,SelectPSTHinterval_SacOn);



TimeSequence_SacOn_mean=nanmean(TimeSequence_SacOn,1);


%Change raster
Raster_SacOn= SelectSpikeTimeInterval(Raster_SacOn_whole,SelectPSTHinterval_SacOn(1)*ones(size(Raster_SacOn_whole,1),1),SelectPSTHinterval_SacOn(end)*ones(size(Raster_SacOn_whole,1),1));



%Align Eye trace with saccade onset
%Select the eye eccentricity between target onset and target offset
Eye_Eccentricity_FixOff=SelectEyeInterval(Eye_Eccentricity,FIXOFFBin,SACENDBin)';
Eye_Eccentricity_FixOff=ReorganizeEye(Eye_Eccentricity_FixOff);



%Select the interval around saccade onset
EyeTime=0:EyeBinWidth:size(Eye_Eccentricity_FixOff,2)-1;

EyeEcc_Aligned=AlignEyeData(Eye_Eccentricity_FixOff,EyeTime,RT_First_Target,[-PlotSacPreBuffer,PlotSacPostBuffer]);

EyeVel_Aligned= abs(diff(EyeEcc_Aligned,1,2));
EyeVel_Aligned=[EyeVel_Aligned,EyeVel_Aligned(:,end)]/EyeBinWidth*1000;


EyeEcc_Aligned_mean=nanmean(EyeEcc_Aligned,1);
EyeEcc_Aligned_sem=nanstd(EyeEcc_Aligned,[],1)./sqrt(sum(~isnan(EyeEcc_Aligned),1));

%% Align with saccade offset
AlignedSpike=SpikeAlignTool(SacEnd_First_Target',SpikeTimeAligned_FixOff,StepSize,BinWidth);

PSTH_SacOff=AlignedSpike.PSTH_Aligned;
TimeSequence_SacOff=AlignedSpike.PSTH_Time_Aligned;
%Select interval
SelectPSTHinterval_SacOff=[-PlotSacPreBuffer,PlotSacPostBuffer];
[PSTH_SacOff,TimeSequence_SacOff]=SelectPSTHInterval(PSTH_SacOff,TimeSequence_SacOff,SelectPSTHinterval_SacOff);



TimeSequence_SacOff_mean=nanmean(TimeSequence_SacOff,1);



%%
%%%%%%%%%%%%%%%%%%%%%%%Align data on object onset
%Segregate spike data from fractal onset to fractal offset

SegmentData=EventSegTool(SegEvent,SegEventEnd,PlotStartOffset,PlotStopOffset,EventChannel,EventTimeChannel,SpikeTime,StepSize,BinWidth);

FR_Each = SegmentData.MeanFRRaw{1};
FR_Mean = mean(FR_Each);
FR_Std = std(FR_Each);

%Retrive data from the field

%Relative time bin for each event
RelativeTimeMarker=SegmentData.RelativeTimeMarker{1};%=EventTimeMarker;
%Absolute time marker for each event
AbosoluteTimeMarker=SegmentData.AbsoluteTimeMarker{1};

ObjOnTimeMarker=nanmean(RelativeTimeMarker);

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


%% Raster after align
RasterTgtOn = SegmentData.RasterAfterAlign{1};



%Get mean and sem of eye eccentricity
%Get eye trance distance during the interval
StartTimeBin=FindOutTime(EventChannel,EventTimeChannel,SegEvent)-PlotStartOffset;
StartTimeBin(StartTimeBin<1)=1;
EndTimeBin=FindOutTime(EventChannel,EventTimeChannel,SegEventEnd )+PlotStopOffset;
  
Eye_Eccentricity_FractalOn=SelectEyeInterval(Eye_Eccentricity,StartTimeBin,EndTimeBin)';
Eye_Eccentricity_FractalOn=ReorganizeEye(Eye_Eccentricity_FractalOn);








BaseLineWindow = [-500, -400];
VisualEpoch = [20,150];
MovementEpoch = [-200,0];


FR_BaseLine = mean(MeanPSTH_FR_FractalOn(:,PSTH_Time_FraOn_All<=BaseLineWindow(2) & PSTH_Time_FraOn_All>=BaseLineWindow(1)),2);
FR_VisualEpoch = mean(MeanPSTH_FR_FractalOn(:,PSTH_Time_FraOn_All<=VisualEpoch(2) & PSTH_Time_FraOn_All>=VisualEpoch(1)),2);
FR_MovementEpoch = mean(PSTH_SacOn(:,TimeSequence_SacOn_mean<=MovementEpoch(2) & TimeSequence_SacOn_mean>=MovementEpoch(1)),2);



%% Loop over the target directions



for i = 1:length(UniqueAngle )
    diff_Angle = TargetAngle-UniqueAngle(i);
     AngleCount(i) = sum(diff_Angle<0.0001 & ControlTrial);

    for j = 1:length(useEcc)
   
        sel = diff_Angle<0.0001 & ControlTrial & TargetEcc == useEcc(j);

           

       
        
        
        
        %Align on target on 
        currPSTH_FR_TgtOn_Loc = MeanPSTH_FR_FractalOn(sel,:);
        TgtOn_Loc_Mean{i,j} = nanmean(currPSTH_FR_TgtOn_Loc,1);
        TgtOn_Loc_Sem{i,j} = nanstd(currPSTH_FR_TgtOn_Loc,[],1)./sqrt(sum(~isnan(currPSTH_FR_TgtOn_Loc),1));
        
        
        
       
        
        AngleEccorder{i,j} = [UniqueAngle(i),useEcc(j)];
        
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
         Diff_Visual = mean(FR_VisualEpoch(sel)) - mean(FR_BaseLine(sel));
         if p<0.05
             h_Visual(i,j) = sign(Diff_Visual);
         else
             h_Visual(i,j) = 0;

         end

         [h,p] = ttest(FR_BaseLine(sel),FR_MovementEpoch(sel));
         Diff_Movement = mean(FR_MovementEpoch(sel)) - mean(FR_BaseLine(sel));
          disp(Diff_Movement);

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
        if(sum(h_v_sig==1)>0)
            %Positive has the priority
            [MaxVisual, index_v ]= max((FRVisSig - mean(FR_BaseLine)));
        else
            [MaxVisual, index_v ]= max(abs(FRVisSig - mean(FR_BaseLine)));
        end
     % [MaxVisual, index_v ]= max((FRVisSig - mean(FR_BaseLine)));
        MaxVisual_Angle = AngleSig(index_v);
        h_vis_max = h_VisOrMov_Sig(index_v);

        if(sum(h_m_sig==1)>0)
            %Positive has the priority
            [MaxMovement, index_m]= max((FRMovementSig - mean(FR_BaseLine)));
        else
            [MaxMovement, index_m]= max(abs(FRMovementSig - mean(FR_BaseLine)));
        end

        
      %  [MaxMovement, index_m]= max((FRMovementSig - mean(FR_BaseLine)));
      
          
        MaxMovement_Angle = AngleSig(index_m);
        h_move_max = h_VisOrMov_Sig(index_m);
%
        if h_vis_max ~=0 & h_move_max ==0 %Vis-Only
            RF_Angle = [MaxVisual_Angle,MaxVisual_Angle];
            RF_h_v = h_v_sig(index_v);
            RF_h_m = h_m_sig(index_v);
        elseif h_vis_max ~=0 & h_move_max ~=0 %Vis-Move
            RF_Angle = [MaxVisual_Angle,MaxMovement_Angle];
            RF_h_v = h_v_sig(index_v);
            RF_h_m = h_m_sig(index_m);
        elseif h_vis_max ==0 & h_move_max ~=0 %Move only
            RF_Angle = [MaxMovement_Angle,MaxMovement_Angle];
            RF_h_v = h_v_sig(index_m);
            RF_h_m = h_m_sig(index_m);


        else %None
            RF_Angle = [MaxVisual_Angle,MaxVisual_Angle];
            RF_h_v = h_v_sig(index_v);
            RF_h_m = h_m_sig(index_v);

        end
%
        
    else
        [MaxVisual, index_v]= max(abs(Mean_Visual- mean(FR_BaseLine)));
        MaxVisual_Angle = UniqueAngle(index_v);
         h_vis_max = h_VisOrMov(index_v);
        
        RF_h_v =h_Visual(index_v);
        [MaxMovement, index_m]= max(abs(Mean_Movement- mean(FR_BaseLine)));
        MaxMovement_Angle = UniqueAngle(index_m);
        RF_h_m =h_Movement(index_m);

        RF_Angle = [MaxVisual_Angle,MaxMovement_Angle];


    end
%{
    if isnan(RF_Angle)
        [MaxVisual, index_v]= max(abs(mean(Mean_Visual(:,1)) - mean(FR_BaseLine)));
       % [MaxVisual, index_v]= max((mean(Mean_Visual(:,1)) - mean(FR_BaseLine)));
         MaxVisual_Angle = UniqueAngle(index_v);
        h_vis_max = h_VisOrMov(index_v);
        RF_Angle = MaxVisual_Angle;
        RF_h_v =h_Visual(index_v);
        
        RF_h_m =h_Movement(index_v);




    end
%}
    
else
    RF_Angle = [UniqueAngle,UniqueAngle];
    RF_h_v =h_Visual;
    RF_h_m =h_Movement;
    



end
%{
%Last chance try to check the average 
LeftSide = TargetAngle > 0 & TargetAngle < 180;
RightSide = TargetAngle > 180 & TargetAngle <360;

Vis_TgtOn_Left = FR_VisualEpoch(LeftSide);
Vis_TgtOn_Right = FR_VisualEpoch(RightSide);

Spon_TgtOn_Left = FR_BaseLine(LeftSide);
Spon_TgtOn_Right = FR_BaseLine(RightSide);


Diff_Visual_Left = mean(Vis_TgtOn_Left) - mean(Spon_TgtOn_Left);
Diff_Visual_Right = mean(Vis_TgtOn_Right) - mean(Spon_TgtOn_Right);

h_Visual_Left = sign(Diff_Visual_Left );   
h_Visual_Right = sign(Diff_Visual_Right );  

[h_left,p_left] = ttest(Spon_TgtOn_Left,Vis_TgtOn_Left);
[h_right,p_right] = ttest(Spon_TgtOn_Right,Vis_TgtOn_Right);


Sac_TgtOn_Left = FR_MovementEpoch(LeftSide);
Sac_TgtOn_Right = FR_MovementEpoch(RightSide);

Diff_Sac_Left = mean(Sac_TgtOn_Left) - mean(Spon_TgtOn_Left);
Diff_Sac_Right = mean(Sac_TgtOn_Right) - mean(Spon_TgtOn_Right);

h_Movement_Left = sign(Diff_Visual_Left );   
h_Movement_Right = sign(Diff_Visual_Right );   

[h_left,p_left_sac] = ttest(Spon_TgtOn_Left,Sac_TgtOn_Left);
[h_right,p_right_sac] = ttest(Spon_TgtOn_Right,Sac_TgtOn_Right);
%}







PSTH_RF_Mean = TgtOn_Loc_Mean{UniqueAngle==RF_Angle(1)};

PSTH_RF_Sem = TgtOn_Loc_Sem{UniqueAngle==RF_Angle(1)};

PSTH_RF_Sac_Mean = SacOn_Loc_Mean{UniqueAngle==RF_Angle(2)};
PSTH_RF_Sac_Sem = SacOn_Loc_Sem{UniqueAngle==RF_Angle(2)};


PSTH_RF_Mean_z = (PSTH_RF_Mean - FR_Mean)/FR_Std;
PSTH_RF_Sac_Mean_z= (PSTH_RF_Sac_Mean - FR_Mean)/FR_Std;


Raster_TgtOn_RF  = Raster_TgtOn_ij{UniqueAngle==RF_Angle(1),1} ;
Raster_SacOn_RF  = Raster_SacOn_ij{UniqueAngle==RF_Angle(2),1} ;



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

if RF_h_v ==0 & RF_h_m ==0
    CellType = 'None';
    CellTypeCode = 4;

end



useAngle = UniqueAngle(AngleCount==max(AngleCount));
if length(useAngle)>1
    useAngle = useAngle(1);
end
if VaryEcc & ~MultiLoop
    for i = 1:length(UniqueEcc )
        sel = TargetAngle == useAngle & ControlTrial & TargetEcc == UniqueEcc(i);

         %Align on target on 
        currPSTH_FR_TgtOn_Loc = MeanPSTH_FR_FractalOn(sel,:);
        TgtOn_Loc_MeanEcc{i} = nanmean(currPSTH_FR_TgtOn_Loc,1);
        TgtOn_Loc_SemEcc{i} = nanstd(currPSTH_FR_TgtOn_Loc,[],1)./sqrt(sum(~isnan(currPSTH_FR_TgtOn_Loc),1));
        
        
        
       
        
        AngleEccorder{i} = [UniqueAngle(i)];
        
         %Align on sac on 
         
         curr_SacOn_Loc = PSTH_SacOn(sel,:);
         SacOn_Loc_MeanEcc{i} = nanmean(curr_SacOn_Loc ,1);
         SacOn_Loc_SemEcc{i} = nanstd(curr_SacOn_Loc ,[],1)./sqrt(sum(~isnan(curr_SacOn_Loc ),1));
         


    end

end

%Get mean and sem of neural responses
if ~MultiLoop
    sel =  ControlTrial & TargetEcc == useEcc;
    MeanPSTH_FR_FractalOn_mean=nanmean(MeanPSTH_FR_FractalOn(sel,:),1);
    MeanPSTH_FR_FractalOn_sem=nanstd(MeanPSTH_FR_FractalOn(sel,:),[],1)./sqrt(sum(~isnan(MeanPSTH_FR_FractalOn(sel,:)),1));

    PSTH_SacOn_mean=nanmean(PSTH_SacOn(sel,:),1);
    PSTH_SacOn_sem=nanstd(PSTH_SacOn(sel,:),[],1)./sqrt(sum(~isnan(PSTH_SacOn(sel,:)),1));

    PSTH_SacOff_mean=nanmean(PSTH_SacOff,1);
    PSTH_SacOff_sem=nanstd(PSTH_SacOff,[],1)./sqrt(sum(~isnan(PSTH_SacOff),1));



    %%
%Saccade Vector Tuning
%Get average TgtOnm PreSac, DurSac,PostSac spikes
%Use spike rate

%Tgt On
PreSacDur=[0,200];%in ms
SpikeRate_TgtOn=nanmean(MeanPSTH_FR_FractalOn(:,[PSTH_Time_FraOn_All<=PreSacDur(2) & PSTH_Time_FraOn_All>=PreSacDur(1)]),2);


%PreSaccade
%PreSacDur=[-PlotSacPreBuffer,0];%in ms
PreSacDur=[-200,0];%in ms
SpikeRate_PreSac=nanmean(PSTH_SacOn(:,[TimeSequence_SacOn_mean<=PreSacDur(2) & TimeSequence_SacOn_mean>=PreSacDur(1)]),2);

%During Saccade:
SacIn=[0*ones(length(RT_First_Target),1),[SacEnd_First_Target-RT_First_Target]'];
SpikeRate_DurSac=SelectFringRateInterval(PSTH_SacOn,TimeSequence_SacOn_mean,SacIn)';


%After Saccade:
%PostSacDur=[0,PlotSacPostBuffer];%in ms
PostSacDur=[0,200];%in ms
SpikeRate_PostSac=nanmean(PSTH_SacOff(:,[TimeSequence_SacOff_mean<=PostSacDur(2) &TimeSequence_SacOff_mean>=PostSacDur(1)]),2);

 %Group responses according to saccade vectors
 %SeparationPoint=[22.5:22.5:360];
 SeparationPoint=[22.5:45:360];
 
 %For this analysis, only select those for useEcc


 
 SaccadeVectorTuning=PoloarHistGroup(SaccadeAngle_First(sel),[SpikeRate_PreSac(sel),SpikeRate_DurSac(sel),SpikeRate_PostSac(sel)],SeparationPoint);
 SaccadeVector=SaccadeVectorTuning.Center;
 FiringRate_PreDurPost_mean=SaccadeVectorTuning.DataMean;
 FiringRate_PreDurPost_sem=SaccadeVectorTuning.DataSEM;
 SacVecTuning_p=SaccadeVectorTuning.ANOVA_p;
 
 PrefVector=SaccadeVectorTuning.PrefVector;
PrefVector_Amp=SaccadeVectorTuning.PrefVector_Amp;

% TgtOn Tuning

TgtOnTuning=PoloarHistGroup(TargetAngle(sel),SpikeRate_TgtOn(sel),SeparationPoint);
 TgtVector=TgtOnTuning.Center;
 FiringRate_TgtOn_mean=TgtOnTuning.DataMean;
 FiringRate_TgtOn_sem=TgtOnTuning.DataSEM;
 TgtOnVecTuning_p=TgtOnTuning.ANOVA_p;
 
 PrefVectorTgtOn=TgtOnTuning.PrefVector;
PrefVector_TgtOnAmp=TgtOnTuning.PrefVector_Amp;


%%

else
    for i = 1:length(useEcc)
        sel =  ControlTrial & TargetEcc == useEcc(i);
        MeanPSTH_FR_FractalOn_mean{i}=nanmean(MeanPSTH_FR_FractalOn(sel,:),1);
        MeanPSTH_FR_FractalOn_sem{i}=nanstd(MeanPSTH_FR_FractalOn(sel,:),[],1)./sqrt(sum(~isnan(MeanPSTH_FR_FractalOn(sel,:)),1));

        PSTH_SacOn_mean{i} = nanmean(PSTH_SacOn(sel,:),1);
        PSTH_SacOn_sem{i} = nanstd(PSTH_SacOn(sel,:),[],1)./sqrt(sum(~isnan(PSTH_SacOn(sel,:)),1));

        PSTH_SacOff_mean{i}=nanmean(PSTH_SacOff(sel,:),1);
        PSTH_SacOff_sem{i}=nanstd(PSTH_SacOff(sel,:),[],1)./sqrt(sum(~isnan(PSTH_SacOff(sel,:)),1));


        %%
%Saccade Vector Tuning
%Get average TgtOnm PreSac, DurSac,PostSac spikes
%Use spike rate

%Tgt On
PreSacDur=[0,200];%in ms
SpikeRate_TgtOn=nanmean(MeanPSTH_FR_FractalOn(sel,[PSTH_Time_FraOn_All<=PreSacDur(2) & PSTH_Time_FraOn_All>=PreSacDur(1)]),2);


%PreSaccade
%PreSacDur=[-PlotSacPreBuffer,0];%in ms
PreSacDur=[-200,0];%in ms
SpikeRate_PreSac=nanmean(PSTH_SacOn(sel,[TimeSequence_SacOn_mean<=PreSacDur(2) & TimeSequence_SacOn_mean>=PreSacDur(1)]),2);

%During Saccade:
SacIn=[0*ones(length(RT_First_Target),1),[SacEnd_First_Target-RT_First_Target]'];
SpikeRate_DurSac_all=SelectFringRateInterval(PSTH_SacOn,TimeSequence_SacOn_mean,SacIn)';

SpikeRate_DurSac = SpikeRate_DurSac_all(sel,:);


%After Saccade:
%PostSacDur=[0,PlotSacPostBuffer];%in ms
PostSacDur=[0,200];%in ms
SpikeRate_PostSac=nanmean(PSTH_SacOff(sel,[TimeSequence_SacOff_mean<=PostSacDur(2) &TimeSequence_SacOff_mean>=PostSacDur(1)]),2);

 %Group responses according to saccade vectors
 %SeparationPoint=[22.5:22.5:360];
 SeparationPoint=[22.5:45:360];
 
 
 
 SaccadeVectorTuning=PoloarHistGroup(SaccadeAngle_First(sel),[SpikeRate_PreSac,SpikeRate_DurSac,SpikeRate_PostSac],SeparationPoint);
 SaccadeVector=SaccadeVectorTuning.Center;
 FiringRate_PreDurPost_mean{i}=SaccadeVectorTuning.DataMean;
 FiringRate_PreDurPost_sem{i}=SaccadeVectorTuning.DataSEM;
 SacVecTuning_p{i}=SaccadeVectorTuning.ANOVA_p;
 
 PrefVector{i}=SaccadeVectorTuning.PrefVector;
PrefVector_Amp{i}=SaccadeVectorTuning.PrefVector_Amp;

% TgtOn Tuning

TgtOnTuning=PoloarHistGroup(TargetAngle(sel),SpikeRate_TgtOn,SeparationPoint);
 TgtVector=TgtOnTuning.Center;
 FiringRate_TgtOn_mean{i}=TgtOnTuning.DataMean;
 FiringRate_TgtOn_sem{i}=TgtOnTuning.DataSEM;
 TgtOnVecTuning_p(i)=TgtOnTuning.ANOVA_p;
 
 PrefVectorTgtOn(i)=TgtOnTuning.PrefVector;
PrefVector_TgtOnAmp(i)=TgtOnTuning.PrefVector_Amp;




    end

end


%Plot range 
PlotTgtOn = [-800,500];
PlotSacOn = [-500,500];

sel_tgt = PSTH_Time_FraOn_All <=PlotTgtOn(2) & PSTH_Time_FraOn_All >=PlotTgtOn(1);
sel_sac = TimeSequence_SacOn_mean <=PlotSacOn(2) & TimeSequence_SacOn_mean >=PlotSacOn(1);



%%
if ShowFigureFlag

    color_seq = {'-b','-g','y'};
    color_seq2 = {'r','g','y'};

 %%%%Figure 1 Overall PSTH and 4 tunings 
 
FigIndex = 110;
FigNum = 1;
    
figtitlestr{FigNum}=sprintf('Tunings');
fig{FigNum}=PrepareFigure(FigIndex,'w',[30,100, 1200,800],'Name',figtitlestr{FigNum});

subplot(2,12,[1,2,3,4]);
if ~MultiLoop
    shadedErrorBar(PSTH_Time_FraOn_All,MeanPSTH_FR_FractalOn_mean,MeanPSTH_FR_FractalOn_sem,'lineprops',{[0.9290    0.6940    0.1250]},'transparent',1,'patchSaturation',0.3);
else
    for i = 1:length(useEcc)
        shadedErrorBar(PSTH_Time_FraOn_All,MeanPSTH_FR_FractalOn_mean{i},MeanPSTH_FR_FractalOn_sem{i},'lineprops',{[1-(i-1)/length(useEcc)    0.6940    0.1250]},'transparent',1,'patchSaturation',0.3);
    end

end

hold on

xlim([-PlotSacPreBuffer,PlotSacPostBuffer]);
yrange=ylim;
plot([0,0],yrange,'--k');
set(findall(gca, 'Type', 'Line'),'LineWidth',3);

if MultiLoop
    legend(string(useEcc));
end

set(gca,'FontSize',15,'FontWeight','Bold');
set(gca,'LineWidth',3);
xlabel('Time from TgtOn Onset','FontSize',20,'FontWeight','Bold');
ylabel('Spikes/s','FontSize',25,'FontWeight','Bold');

subplot(2,12,[5,6,7,8]);
if ~MultiLoop
    shadedErrorBar(TimeSequence_SacOn_mean,PSTH_SacOn_mean,PSTH_SacOn_sem,'lineprops',{[0.9290    0.6940    0.1250]},'transparent',1,'patchSaturation',0.3);
else
    for i = 1:length(useEcc)
        shadedErrorBar(TimeSequence_SacOn_mean,PSTH_SacOn_mean{i},PSTH_SacOn_sem{i},'lineprops',{[1-(i-1)/length(useEcc)    0.6940    0.1250 ]},'transparent',1,'patchSaturation',0.3);
        hold on

    end
end

hold on
%xlim([-PlotSacPreBuffer,PlotSacPostBuffer]);
yrange=ylim;
plot([0,0],yrange,'--k');
set(findall(gca, 'Type', 'Line'),'LineWidth',3);

set(gca,'FontSize',15,'FontWeight','Bold');
set(gca,'LineWidth',3);
xlabel('Time from SacOn Onset','FontSize',20,'FontWeight','Bold');
ylabel('Spikes/s','FontSize',25,'FontWeight','Bold');

subplot(2,12,[9,10,11,12]);
if ~MultiLoop
    shadedErrorBar(TimeSequence_SacOff_mean,PSTH_SacOff_mean,PSTH_SacOff_sem,'lineprops',{[0.9290    0.6940    0.1250]},'transparent',1,'patchSaturation',0.3);
else
    for i = 1:length(useEcc)
        shadedErrorBar(TimeSequence_SacOff_mean,PSTH_SacOff_mean{i},PSTH_SacOff_sem{i},'lineprops',{[1-(i-1)/length(useEcc)    0.6940    0.1250]},'transparent',1,'patchSaturation',0.3);
        hold on 
    end


end
hold on
xlim([-PlotSacPreBuffer,PlotSacPostBuffer]);
yrange=ylim;
plot([0,0],yrange,'--k');
set(findall(gca, 'Type', 'Line'),'LineWidth',3);



set(gca,'FontSize',15,'FontWeight','Bold');
set(gca,'LineWidth',3);
xlabel('Time from Sac Offset','FontSize',20,'FontWeight','Bold');
ylabel('Spikes/s','FontSize',25,'FontWeight','Bold');

subplot(2,12,[13,14,15]);

if ~MultiLoop
    polarwitherrorbar([TgtVector/180*pi TgtVector(1)/180*pi],[FiringRate_TgtOn_mean' FiringRate_TgtOn_mean(1)],[FiringRate_TgtOn_sem' FiringRate_TgtOn_sem(1)],'-b',3);

    hold on
    h_p = polarplot([PrefVectorTgtOn PrefVectorTgtOn]/180*pi,[0 PrefVector_TgtOnAmp],'r');
    set(h_p,'LineWidth',2.5);
    maxTgtOnFR = max(FiringRate_TgtOn_mean);
else
    

    for i = 1:length(useEcc)

        polarwitherrorbar([TgtVector/180*pi TgtVector(1)/180*pi],[FiringRate_TgtOn_mean{i}' FiringRate_TgtOn_mean{i}(1)],[FiringRate_TgtOn_sem{i}' FiringRate_TgtOn_sem{i}(1)],color_seq{i},3);
        hold on
        h_p = polarplot([PrefVectorTgtOn(i) PrefVectorTgtOn(i)]/180*pi,[0 PrefVector_TgtOnAmp(i)],color_seq2{i});
        set(h_p,'LineWidth',2.5);


        
    end
    maxTgtOnFR = max(cellfun(@max,FiringRate_TgtOn_mean));

end

set(gca,'ThetaDir' , 'counterclockwise');
set(gca,'ThetaZeroLocation','top')
set(gca,'FontSize',25,'FontWeight','Bold','LineWidth',3);

thetaticks([0:45:359]);
thetaticklabels([0:45:359]);



%thetaticks(sort(SaccadeVector));
%thetaticklabels(sort(SaccadeVector));
title('TgtOn','FontSize',25,'FontWeight','Bold');
if ~MultiLoop
pstr=sprintf('p= %2.3f',TgtOnVecTuning_p);
a=annotation('textbox',[0.1,0.015,0.1,0.1] ,'String',pstr,'FitBoxToText','on','LineStyle','none');
a.FontSize=30;
a.FontWeight='Bold';
if TgtOnVecTuning_p<0.05
    a.Color='r';
end
else
    for i = 1:length(useEcc)
        pstr=sprintf('p= %2.3f',TgtOnVecTuning_p(i));
        a=annotation('textbox',[0.1,0.05-0.05*(i-1),0.1,0.1] ,'String',pstr,'FitBoxToText','on','LineStyle','none');
        a.FontSize=30;
        a.FontWeight='Bold';
        if TgtOnVecTuning_p(i) < 0.05
            a.Color='r';
        end
    end
end

rlim([0,maxTgtOnFR*1.2]);


%% Sac vector tuning
subtitlestr={'Pre-Sac Tuning','Dur-Sac Tuning','Post-Sac Tuning'};

if ~MultiLoop

Selection_temp=~isnan(SaccadeVector);
SaccadeVector=SaccadeVector(Selection_temp);
FiringRate_PreDurPost_mean=FiringRate_PreDurPost_mean(Selection_temp,:);
FiringRate_PreDurPost_sem=FiringRate_PreDurPost_sem(Selection_temp,:);


for j=1:3
%Pre-saccade , During, Post saccade tuning
%subplot(1,3,j);
subplot(2,12,[16+(j-1)*3,17+(j-1)*3,18+(j-1)*3]);
%{
p= polarplot((SaccadeAngle_First)*pi/180',SpikeRate_PreSac,'or');
p.MarkerFaceColor='r';
p.MarkerEdgeColor='r';
%}

%hold on
if ~isempty(FiringRate_PreDurPost_mean(:,j))

    
    polarwitherrorbar([SaccadeVector/180*pi SaccadeVector(1)/180*pi],[FiringRate_PreDurPost_mean(:,j)' FiringRate_PreDurPost_mean(1,j)],[FiringRate_PreDurPost_sem(:,j)' FiringRate_PreDurPost_sem(1,j)],...
        '-b',3);
    hold on
    h_p = polarplot([PrefVector(j) PrefVector(j)]/180*pi,[0 PrefVector_Amp(j)],'r');
        set(h_p,'LineWidth',2.5);
    
set(gca,'ThetaDir' , 'counterclockwise');
set(gca,'ThetaZeroLocation','top')
set(gca,'FontSize',25,'FontWeight','Bold','LineWidth',3);

thetaticks([0:45:359]);
thetaticklabels([0:45:359]);

%thetaticks(sort(SaccadeVector));
%thetaticklabels(sort(SaccadeVector));
title(subtitlestr{j},'FontSize',25,'FontWeight','Bold');
pstr=sprintf('p= %2.3f',SacVecTuning_p(j));
a=annotation('textbox',[0.4+(j-1)*0.2,0.015,0.1,0.1] ,'String',pstr,'FitBoxToText','on','LineStyle','none');
a.FontSize=30;
a.FontWeight='Bold';
if SacVecTuning_p(j)<0.05
    a.Color='r';
end

rlim([0,max(max(FiringRate_PreDurPost_mean))*1.2]);
end
end


else

    FiringRate_PreDurPost_mean_All =FiringRate_PreDurPost_mean;
    FiringRate_PreDurPost_sem_All = FiringRate_PreDurPost_sem;

    %IF multi-loop
    for i = 1:length(useEcc)'

        Selection_temp=~isnan(SaccadeVector);
SaccadeVector=SaccadeVector(Selection_temp);
FiringRate_PreDurPost_mean=FiringRate_PreDurPost_mean_All{i}(Selection_temp,:);
FiringRate_PreDurPost_sem=FiringRate_PreDurPost_sem_All{i}(Selection_temp,:);


for j=1:3
%Pre-saccade , During, Post saccade tuning
%subplot(1,3,j);
subplot(2,12,[16+(j-1)*3,17+(j-1)*3,18+(j-1)*3]);
%{
p= polarplot((SaccadeAngle_First)*pi/180',SpikeRate_PreSac,'or');
p.MarkerFaceColor='r';
p.MarkerEdgeColor='r';
%}

%hold on
if ~isempty(FiringRate_PreDurPost_mean(:,j))

    
    polarwitherrorbar([SaccadeVector/180*pi SaccadeVector(1)/180*pi],[FiringRate_PreDurPost_mean(:,j)' FiringRate_PreDurPost_mean(1,j)],[FiringRate_PreDurPost_sem(:,j)' FiringRate_PreDurPost_sem(1,j)],...
        color_seq{i},3);
    hold on
    h_p = polarplot([PrefVector{i}(j) PrefVector{i}(j)]/180*pi,[0 PrefVector_Amp{i}(j)],color_seq2{i});
        set(h_p,'LineWidth',2.5);
    
set(gca,'ThetaDir' , 'counterclockwise');
set(gca,'ThetaZeroLocation','top')
set(gca,'FontSize',25,'FontWeight','Bold','LineWidth',3);

thetaticks([0:45:359]);
thetaticklabels([0:45:359]);

%thetaticks(sort(SaccadeVector));
%thetaticklabels(sort(SaccadeVector));
title(subtitlestr{j},'FontSize',25,'FontWeight','Bold');
pstr=sprintf('p= %2.3f',SacVecTuning_p{i}(j));
a=annotation('textbox',[0.4+(j-1)*0.2,0.05-0.05*(i-1),0.1,0.1] ,'String',pstr,'FitBoxToText','on','LineStyle','none');
a.FontSize=30;
a.FontWeight='Bold';
if SacVecTuning_p{i}(j)<0.05
    a.Color='r';
end

%rlim([0,max(max(FiringRate_PreDurPost_mean))*1.2]);
end
end







    end
end

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Figure 2: Average PSTH
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%for each location


FigIndex = FigIndex +1;
FigNum = FigNum+1;
figtitlestr{FigNum}=sprintf('PSTH_Ecc %d',useEcc);
fig{FigNum}=PrepareFigure(FigIndex,'w',[50,100, 1200,800],'Name',figtitlestr{FigNum});
    
    
   
    for j = 1:length(UniqueAngle)
        
        
          ax(j) = subplot(ceil(length(UniqueAngle)/3),3,j);

           for i = 1:length(useEcc)
            %Align on tgtonset
            shadedErrorBar(PSTH_Time_FraOn_All(sel_tgt),TgtOn_Loc_Mean{j,i}(sel_tgt),TgtOn_Loc_Sem{j,i}(sel_tgt),'lineprops',{[1    0+0.5*(i-1)    0]},'transparent',1,'patchSaturation',0.3);
            hold on
           end
            set(findall(gca, 'Type', 'Line'),'LineWidth',3);

        set(gca,'FontSize',15,'FontWeight','Bold');
        set(gca,'LineWidth',3);
       % title(sprintf("Angle: %1.1d,Ecc: %1.1d",AngleEccorder{j,i}));
        title(sprintf("Angle: %1.1d",AngleEccorder{j,1}(1)));
       
            xlabel('Time from Tgt Onset(ms)','FontSize',15,'FontWeight','Bold');
      
        
        ylabel('Firing Rate(Hz)','FontSize',15,'FontWeight','Bold');
        
       
      
        linkaxes(ax,'xy');
           
    end
    
    %Add title
    legend(string(useEcc));

title_msg=sprintf(' %s, SpikeChannel: %d, TrialNum: %d ',FileName,SpikeChannel,TrialNum);
a=annotation('textbox',[0.25,0.7,0.3,0.3],'String',title_msg,'FitBoxToText','on');
a.FontSize=25;
a.FontWeight='Bold';
a.LineStyle='none';
    
    
%% Align on saccade onset
FigNum =FigNum+1;
FigIndex =FigIndex+1;

  figtitlestr{FigNum}=sprintf('PSTH_Ecc %d',useEcc);
   fig{FigNum}=PrepareFigure(FigIndex,'w',[50,100, 1200,800],'Name',figtitlestr{FigNum});
    
    
   
    for j = 1:length(UniqueAngle)
        
        
        
          ax2(j) = subplot(ceil(length(UniqueAngle)/3),3,j);
         
        for i = 1:length(useEcc)
       
       %  
          %Align on saccade onset
          
         shadedErrorBar(TimeSequence_SacOn_mean(sel_sac ),SacOn_Loc_Mean{j,i}(sel_sac ),SacOn_Loc_Sem{j,i}(sel_sac ),'lineprops',{[0    0    1-0.3*(i-1)]},'transparent',1,'patchSaturation',0.3);
         hold on
          set(findall(gca, 'Type', 'Line'),'LineWidth',3);

        set(gca,'FontSize',15,'FontWeight','Bold');
        set(gca,'LineWidth',3);
        title(sprintf("Angle: %1.1d",AngleEccorder{j,i}(1)));
       
            xlabel('Time from Sac Onset(ms)','FontSize',15,'FontWeight','Bold');
            
        
        
        ylabel('Firing Rate(Hz)','FontSize',15,'FontWeight','Bold');
        
        
        
       %
        linkaxes(ax2,'xy');
        end
    end
    
    %Add title
legend(string(useEcc));
title_msg=sprintf(' %s, SpikeChannel: %d, TrialNum: %d ',FileName,SpikeChannel,TrialNum);
a=annotation('textbox',[0.25,0.7,0.3,0.3],'String',title_msg,'FitBoxToText','on');
a.FontSize=25;
a.FontWeight='Bold';
a.LineStyle='none';

%% If vary eccentricity
if VaryEcc & ~MultiLoop
    FigNum =FigNum+1;
    FigIndex =FigIndex+1;

    figtitlestr{FigNum}=sprintf('PSTH_VaryEcc');
    fig{FigNum}=PrepareFigure(FigIndex,'w',[70,100, 1200,800],'Name',figtitlestr{FigNum});
    
    
   
    for j = 1:length(UniqueEcc)
        
        
          ax3(j) = subplot(2,length(UniqueEcc),j);

           %Align on tgtonset
          shadedErrorBar(PSTH_Time_FraOn_All(sel_tgt),TgtOn_Loc_MeanEcc{j}(sel_tgt),TgtOn_Loc_SemEcc{j}(sel_tgt),'lineprops',{[1    0    0]},'transparent',1,'patchSaturation',0.3);
         
          
            set(findall(gca, 'Type', 'Line'),'LineWidth',3);

        set(gca,'FontSize',15,'FontWeight','Bold');
        set(gca,'LineWidth',3);
        title(sprintf("Ecc: %1.1d",UniqueEcc(j)));
       
            xlabel('Time from Tgt Onset(ms)','FontSize',15,'FontWeight','Bold');
      
        
        ylabel('Firing Rate(Hz)','FontSize',15,'FontWeight','Bold');
        
       
       
       %  
          ax3(length(UniqueEcc)+j) = subplot(2,length(UniqueEcc),length(UniqueEcc)+j);
          %Align on saccade onset
          
         shadedErrorBar(TimeSequence_SacOn_mean(sel_sac ),SacOn_Loc_MeanEcc{j}(sel_sac ),SacOn_Loc_SemEcc{j}(sel_sac ),'lineprops',{[0    0    1]},'transparent',1,'patchSaturation',0.3);
         hold on
          set(findall(gca, 'Type', 'Line'),'LineWidth',3);

        set(gca,'FontSize',15,'FontWeight','Bold');
        set(gca,'LineWidth',3);
        title(sprintf("Ecc: %1.1d",UniqueEcc(j)));
       
            xlabel('Time from Sac Onset(ms)','FontSize',15,'FontWeight','Bold');
            
        
        
        ylabel('Firing Rate(Hz)','FontSize',15,'FontWeight','Bold');
        
        
        
       %
        linkaxes(ax3,'xy');
    end
    
    %Add title

title_msg=sprintf(' %s, SpikeChannel: %d, TrialNum: %d ',FileName,SpikeChannel,TrialNum);
a=annotation('textbox',[0.25,0.7,0.3,0.3],'String',title_msg,'FitBoxToText','on');
a.FontSize=25;
a.FontWeight='Bold';
a.LineStyle='none';



end % End of if VarryEcc


 FigNum =FigNum+1;
    FigIndex =FigIndex+1;

    %% Visual-Motor responses in RF

figtitlestr{FigNum}=sprintf('RF;Angle %d,%d, %s',RF_Angle,CellType);
    fig{FigNum}=PrepareFigure(FigIndex,'w',[100,100, 600,400],'Name',figtitlestr{FigNum});

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
    shadedErrorBar(PSTH_Time_FraOn_All,TgtOn_Loc_Mean{UniqueAngle==RF_Angle(1)},TgtOn_Loc_Sem{UniqueAngle==RF_Angle(1)},'lineprops',{[1    0    0]},'transparent',1,'patchSaturation',0.3);
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
     shadedErrorBar(TimeSequence_SacOn_mean,SacOn_Loc_Mean{UniqueAngle==RF_Angle(2)},SacOn_Loc_Sem{UniqueAngle==RF_Angle(2)},'lineprops',{[0    0    1]},'transparent',1,'patchSaturation',0.3);
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



ExistFlag=0;
%Load the old file if exist
if exist(OutputFileName)

    load(OutputFileName);
    %Load OutputData into the memory
    if isfield(OutputData,'DelaySaccadeTuning')
    ExistFlag=1;
    end
end


%%%%%%%%%%%%%%%%%%%%OutputData For current analysis
OutputData_New=[];

OutputData_New.DelaySaccadeTuning.TrialType='GoodTrials';
%%%%%%%%%Save overall scene response characteristics%%%%%%%%%%%%%%%%%%%%%%%
%Organize the data



NameStr={
    'PSTH_FR_TgtOn_Time','PSTH_FR_TgtOn',...
    'PSTH_FR_SacOn_Time','PSTH_FR_SacOn',...
    'PSTH_FR_SacOff_ITI_Time','PSTH_FR_SacOff_ITI',...
    'SacVec','SacVectTuning_PreDurPost_mean','SacVectTuning_PreDurPost_sem','SacVectTuning_p','Pref_Vect','Pref_Vect_Length',...
    'TgtVec','TgtVectTuning_mean','TgtVectTuning_sem','TgtVectTuning_p','Pref_TgtVect','Pref_TgtVect_Length',...
    'SacOnLocMean','SacOnLocSem',...
    'PSTH_TgtOn_RF','PSTH_SacOn_RF',...
    'RF_h_v','RF_h_m','CellType',...
    'PSTH_TgtOn_RF_z','PSTH_SacOn_RF_z'
    
    
};

DataLib={
    PSTH_Time_FraOn_All,[MeanPSTH_FR_FractalOn_mean',MeanPSTH_FR_FractalOn_sem'],...
    TimeSequence_SacOn_mean,[PSTH_SacOn_mean',PSTH_SacOn_sem'],...
    TimeSequence_SacOff_mean,[PSTH_SacOff_mean',PSTH_SacOff_sem'],...
    SaccadeVector,FiringRate_PreDurPost_mean,FiringRate_PreDurPost_sem,SacVecTuning_p,PrefVector,PrefVector_Amp,...
    TgtVector,FiringRate_TgtOn_mean,FiringRate_TgtOn_sem,TgtOnVecTuning_p,PrefVectorTgtOn,PrefVector_TgtOnAmp,...
    SacOn_Loc_Mean,SacOn_Loc_Sem,...
    PSTH_RF_Mean,PSTH_RF_Sac_Mean,...
    RF_h_v,RF_h_m,CellTypeCode,...
    PSTH_RF_Mean_z,PSTH_RF_Sac_Mean_z

    
};



DataStamp=containers.Map(NameStr,DataLib);





OutputData_New.DelaySaccadeTuning.Task=TaskType;
OutputData_New.DelaySaccadeTuning.TaskCode=TaskCode;
OutputData_New.DelaySaccadeTuning.DataStamp=DataStamp;


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
Task_Old=OutputData.DelaySaccadeTuning.TaskCode;

     if Task_Old~=TaskCode
        %Not the same task, add the new dataset into the old one
        OutputData.DelaySaccadeTuning(numel(OutputData.DelaySaccadeTuning)+1)=OutputData_New.DelaySaccadeTuning;
     else
         %If the same task,the same protocol,replace the old one with the new one
         
         OutputData.DelaySaccadeTuning=OutputData_New.DelaySaccadeTuning;
         
    
    
     end
else
    %Output the current file
   OutputData.DelaySaccadeTuning=OutputData_New.DelaySaccadeTuning; 
end



if OutputFlag
    cd(OutPath);
    save(OutputFileName,'OutputData');
    disp('Data Exported');
else
    disp('No data export');
end  



end %End of Spike Channel loop


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
function FR=SelectFringRateInterval(FR_Raw,Time,Interval);
for i=1:size(Interval,1)
    if sum(~isnan(Interval(i,:)))==2
        FR(i)=nanmean(FR_Raw(i,Time>=Interval(i,1) & Time<=Interval(i,2)));
    else
         FR(i)=NaN;
    end
    
    
    
end




end