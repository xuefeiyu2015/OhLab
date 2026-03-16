function  SaccadeAnalysis(Data,  StartTrial, EndTrial,ShowFigureFlag,OutputFlag);

%Fuction for analysis of all kinds of saccades
%Xuefei Yu 06192020
ECode;

SegEvent=[CORRECTCD];%Event Flow of interest
%SegEventEnd=[SHOWFIXCD,FIXOFF,TGTOFF,RWDOFFCD ,NaN,EXTRAFP];%Event end
SegEventEnd=[NaN];%Event end
SegEventStr={'InterTrialInterval'};

PlotStartOffset=[-500];%in ms
PlotStopOffset=[NaN];%in ms

%Parameter to define the sliding window to measure firing rates
StepSize=5;%in ms
BinWidth=15;%in ms


%Criterias for selecting saccades
V_Threshold=40;
ContinueBin=5;
ISI_Threshold=20;

%Time bin pre saccade and time bin offset after saccade onset
PlotSacPreBuffer=300;%in ms
PlotSacPostBuffer=400;%in ms




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
Selection=StartTrial:EndTrial;


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


%% Load the spike channel
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
        %{
    if SpikeChannel<4
        SpikeTime = SpikeTime;%Correct for the delay of online channel
    end
        %}
    ChannelID=ChannelIDWhole(spk);
    ClusterID=ClusterIDWhole(spk);




%Segregate spike data from task end to trial end

SegmentData=EventSegTool(SegEvent,SegEventEnd,PlotStartOffset,PlotStopOffset,EventChannel,EventTimeChannel,SpikeTime,StepSize,BinWidth);
PSTH_ITI=SegmentData.MeanPSTH;
MeanPSTH_Time_ITI=PSTH_ITI{1};
MeanPSTH_FR_ITI=PSTH_ITI{2};

%Remove cols whose number of data is less than 10

ColDataNumThreshold=ceil(TrialNum*1/10);

Col_Num=sum(~isnan(MeanPSTH_FR_ITI),1);
ValidCol=Col_Num>ColDataNumThreshold;

MeanPSTH_Time_ITI=MeanPSTH_Time_ITI(:,ValidCol);
MeanPSTH_FR_ITI=MeanPSTH_FR_ITI(:,ValidCol);

MaxTime=max(MeanPSTH_Time_ITI,[],2);
MinTime=min(MeanPSTH_Time_ITI,[],2);


%Select the eyetrace during ITI
TaskEndBin=FindOutTime(EventChannel,EventBin,SegEvent(1));
EyeSelectionStart=ceil((MinTime-StepSize/2)/EyeBinWidth)+TaskEndBin;

EyeSelectionEnd=ceil((MaxTime-StepSize/2)/EyeBinWidth)+TaskEndBin;

EyeSelectionStart_Time=(EyeSelectionStart-1)*EyeBinWidth;



EyeChannel_X_ITI= SelectEyeInterval(EyeChannel_X_L,EyeSelectionStart,EyeSelectionEnd)';
EyeChannel_Y_ITI= SelectEyeInterval(EyeChannel_Y_L,EyeSelectionStart,EyeSelectionEnd)';

EyeChannel_X_ITI=ReorganizeEye(EyeChannel_X_ITI);
EyeChannel_Y_ITI=ReorganizeEye(EyeChannel_Y_ITI);

Eye_Eccentricity_ITI=SelectEyeInterval(Eye_Eccentricity,EyeSelectionStart,EyeSelectionEnd)';
Eye_Eccentricity_ITI=ReorganizeEye(Eye_Eccentricity_ITI);

EyeTime_ITI=0:EyeBinWidth:size(Eye_Eccentricity_ITI,2)-EyeBinWidth;

MaxTimeLength=MaxTime-MinTime;


%Scatter plot of the eye x,y and the firing rate 
%Desampling the eye trace to match firing rate;
SmoothNum=ceil(StepSize/EyeBinWidth);
EyeChannel_X_ITI_smooth=smoothdata(EyeChannel_X_ITI,SmoothNum);
EyeChannel_Y_ITI_smooth=smoothdata(EyeChannel_Y_ITI,SmoothNum);


EyeChannel_X_ITI_Desample_abs=EyeChannel_X_ITI_smooth(:,1:SmoothNum:end);
EyeChannel_Y_ITI_Desample_abs=EyeChannel_Y_ITI_smooth(:,1:SmoothNum:end);
EyeTime_ITI_Desample=EyeTime_ITI+StepSize/2;

MeanPSTH_FR_ITI_abs=MeanPSTH_FR_ITI;

ShiftTimeLib=[-200,0,100,200];

for j=1:length(ShiftTimeLib)
    ShiftTime=ShiftTimeLib(j);
    if ShiftTime>=0
        ShiftTime=ShiftTime/StepSize+1;
        %Firing Delayed
        
EyeChannel_X_ITI_Desample=EyeChannel_X_ITI_Desample_abs(:,1:end-ShiftTime+1);
EyeChannel_Y_ITI_Desample=EyeChannel_Y_ITI_Desample_abs(:,1:end-ShiftTime+1);


EyeChannel_Ecc_ITI_Desample=sqrt(EyeChannel_X_ITI_Desample.^2+EyeChannel_Y_ITI_Desample.^2);
EyeChannel_V_ITI_Desample=abs(diff(EyeChannel_Ecc_ITI_Desample,1,2)/StepSize*1000);

EyeChannel_V_ITI_Desample=[EyeChannel_V_ITI_Desample,EyeChannel_V_ITI_Desample(:,end)];



MeanPSTH_FR_ITI=MeanPSTH_FR_ITI_abs(:,ShiftTime:end);
    else
        %Eye delayed
        ShiftTime=-ShiftTime/StepSize+1;
        EyeChannel_X_ITI_Desample=EyeChannel_X_ITI_Desample_abs(:,ShiftTime:end);
EyeChannel_Y_ITI_Desample=EyeChannel_Y_ITI_Desample_abs(:,ShiftTime:end);


EyeChannel_Ecc_ITI_Desample=sqrt(EyeChannel_X_ITI_Desample.^2+EyeChannel_Y_ITI_Desample.^2);
EyeChannel_V_ITI_Desample=abs(diff(EyeChannel_Ecc_ITI_Desample,1,2)/StepSize*1000);
EyeChannel_V_ITI_Desample=[EyeChannel_V_ITI_Desample,EyeChannel_V_ITI_Desample(:,end)];


MeanPSTH_FR_ITI=MeanPSTH_FR_ITI_abs(:,1:end-ShiftTime+1);

    end


EyeChannel_X_ITI_Desample_Line=reshape(EyeChannel_X_ITI_Desample,numel(EyeChannel_X_ITI_Desample),1);
EyeChannel_Y_ITI_Desample_Line=reshape(EyeChannel_Y_ITI_Desample,numel(EyeChannel_Y_ITI_Desample),1);
EyeChannel_V_ITI_Desample_Line=reshape(EyeChannel_V_ITI_Desample,numel(EyeChannel_V_ITI_Desample),1);


MeanPSTH_FR_ITI_Line=reshape(MeanPSTH_FR_ITI,numel(MeanPSTH_FR_ITI),1);
[MeanPSTH_FR_ITI_Line_Sort,index]=sort(MeanPSTH_FR_ITI_Line,'ascend');

tmp_nan=~isnan(MeanPSTH_FR_ITI_Line_Sort);
MeanPSTH_FR_ITI_Line_Sort=MeanPSTH_FR_ITI_Line_Sort(tmp_nan);

EyeChannel_X_ITI_Desample_Line_Sort=EyeChannel_X_ITI_Desample_Line(index);
EyeChannel_Y_ITI_Desample_Line_Sort=EyeChannel_Y_ITI_Desample_Line(index);



EyeChannel_X_ITI_Desample_Line_Sort=EyeChannel_X_ITI_Desample_Line_Sort(tmp_nan);
EyeChannel_Y_ITI_Desample_Line_Sort=EyeChannel_Y_ITI_Desample_Line_Sort(tmp_nan);

EyeChannel_X_ITI_Desample_Line_Sort_Plot{j}=EyeChannel_X_ITI_Desample_Line_Sort;
EyeChannel_Y_ITI_Desample_Line_Sort_Plot{j}=EyeChannel_Y_ITI_Desample_Line_Sort;
MeanPSTH_FR_ITI_Line_Sort_Plot{j}=MeanPSTH_FR_ITI_Line_Sort;

%EyeChannel_V_ITI_Desample_Line_Sort=EyeChannel_V_ITI_Desample_Line(index);
%EyeChannel_V_ITI_Desample_Line_Sort=EyeChannel_V_ITI_Desample_Line_Sort(tmp_nan);
end


%}
%Analysis saccades during this period

Saccades=SelectSaccade(EyeChannel_X_ITI,EyeChannel_Y_ITI,EyeBinWidth,V_Threshold,ContinueBin,ISI_Threshold);

%SaccadeNumber=arrayfun(@(x)  x.NumOfSaccade,Saccades)';
SaccadeAngle=arrayfun(@(x)  x.SaccadeAngle,Saccades,'UniformOutput' ,0)';
SaccadeAmplitude=arrayfun(@(x)  x.SaccadeAmplitude,Saccades,'UniformOutput' ,0)';

SaccadeStartTime=arrayfun(@(x)  x.SaccadeStartTime,Saccades,'UniformOutput' ,0)';
SaccadeEndTime=arrayfun(@(x)  x.SaccadeEndTime,Saccades,'UniformOutput' ,0)';

SaccadeStartPoint=arrayfun(@(x)  x.SaccadeStartPoint,Saccades,'UniformOutput' ,0)';
SaccadeEndPoint=arrayfun(@(x)  x.SaccadeEndPoint,Saccades,'UniformOutput' ,0)';

%Screen out the valid saccades and saccade start time severas

AmplitudeThreshold=5;
TimeThreshold=0;
%TimeThresholdEnd=num2cell(MaxTimeLength-PlotSacPostBuffer);
%{
SaccadeAngle_ITI=fillinNaNs(cellfun(@(x,y,z,m) z(((x>AmplitudeThreshold)&(y<=m)))',SaccadeAmplitude,SaccadeEndTime,SaccadeAngle,TimeThresholdEnd,'UniformOutput' ,0));
SaccadeEndPoint_ITI=(cellfun(@(x,y,z,m) z(((x>AmplitudeThreshold)&(y<=m)),:)',SaccadeAmplitude,SaccadeEndTime,SaccadeEndPoint,TimeThresholdEnd,'UniformOutput' ,0));
SaccadeStartPoint_ITI=(cellfun(@(x,y,z,m) z(((x>AmplitudeThreshold)&(y<=m)),:)',SaccadeAmplitude,SaccadeEndTime,SaccadeStartPoint,TimeThresholdEnd,'UniformOutput' ,0));
%}

SaccadeAngle_ITI=fillinNaNs(cellfun(@(x,y,z) z(((x>AmplitudeThreshold)&(y>TimeThreshold)))',SaccadeAmplitude,SaccadeStartTime,SaccadeAngle,'UniformOutput' ,0));
SaccadeEndPoint_ITI=(cellfun(@(x,y,z) z(((x>AmplitudeThreshold)&(y>TimeThreshold)),:)',SaccadeAmplitude,SaccadeStartTime,SaccadeEndPoint,'UniformOutput' ,0));
SaccadeStartPoint_ITI=(cellfun(@(x,y,z) z(((x>AmplitudeThreshold)&(y>TimeThreshold)),:)',SaccadeAmplitude,SaccadeStartTime,SaccadeStartPoint,'UniformOutput' ,0));
%}


SaccadeNumInEachTrial=cellfun(@(x) size(x,2),SaccadeEndPoint_ITI);

SaccadeEndPoint_ITI_X=NaN*ones(length(SaccadeNumInEachTrial),max(SaccadeNumInEachTrial));
SaccadeEndPoint_ITI_Y=NaN*ones(length(SaccadeNumInEachTrial),max(SaccadeNumInEachTrial));

SaccadeStartPoint_ITI_X=NaN*ones(length(SaccadeNumInEachTrial),max(SaccadeNumInEachTrial));
SaccadeStartPoint_ITI_Y=NaN*ones(length(SaccadeNumInEachTrial),max(SaccadeNumInEachTrial));


SaccadeEndPoint_ITI_X(SaccadeNumInEachTrial>0,:)=fillinNaNs(cellfun(@(x) x(1,:),SaccadeEndPoint_ITI(SaccadeNumInEachTrial>0,:),'UniformOutput' ,0));
SaccadeEndPoint_ITI_Y(SaccadeNumInEachTrial>0,:)=fillinNaNs(cellfun(@(x) x(2,:),SaccadeEndPoint_ITI(SaccadeNumInEachTrial>0,:),'UniformOutput' ,0));

SaccadeStartPoint_ITI_X(SaccadeNumInEachTrial>0,:)=fillinNaNs(cellfun(@(x) x(1,:),SaccadeStartPoint_ITI(SaccadeNumInEachTrial>0,:),'UniformOutput' ,0));
SaccadeStartPoint_ITI_Y(SaccadeNumInEachTrial>0,:)=fillinNaNs(cellfun(@(x) x(2,:),SaccadeStartPoint_ITI(SaccadeNumInEachTrial>0,:),'UniformOutput' ,0));

%{
SaccadeStartTime_ITI=fillinNaNs(cellfun(@(x,y,z,m) z(((x>AmplitudeThreshold)&(y<=m)),:)',SaccadeAmplitude,SaccadeStartTime,SaccadeStartTime,TimeThresholdEnd,'UniformOutput' ,0));
SaccadeEndTime_ITI=fillinNaNs(cellfun(@(x,y,z,m) z(((x>AmplitudeThreshold)&(y<=m)),:)',SaccadeAmplitude,SaccadeStartTime,SaccadeEndTime,TimeThresholdEnd,'UniformOutput' ,0));


SaccadeAmplitude_ITI=fillinNaNs(cellfun(@(x,y,z,m) z(((x>AmplitudeThreshold)&(y<=m)),:)',SaccadeAmplitude,SaccadeStartTime,SaccadeAmplitude,TimeThresholdEnd,'UniformOutput' ,0));
%}
%

SaccadeStartTime_ITI=fillinNaNs(cellfun(@(x,y,z) z(((x>AmplitudeThreshold)&(y>TimeThreshold)),:)',SaccadeAmplitude,SaccadeStartTime,SaccadeStartTime,'UniformOutput' ,0));
SaccadeEndTime_ITI=fillinNaNs(cellfun(@(x,y,z) z(((x>AmplitudeThreshold)&(y>TimeThreshold)),:)',SaccadeAmplitude,SaccadeStartTime,SaccadeEndTime,'UniformOutput' ,0));


SaccadeAmplitude_ITI=fillinNaNs(cellfun(@(x,y,z) z(((x>AmplitudeThreshold)&(y>TimeThreshold)),:)',SaccadeAmplitude,SaccadeStartTime,SaccadeAmplitude,'UniformOutput' ,0));

%}
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Align spike time with saccade
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%onset%%%%%%%%%%%%%%%%%
SaccadeStartTime_Absolute_ITI=SaccadeStartTime_ITI+EyeSelectionStart_Time;
AlignedSpike=SpikeAlignTool(SaccadeStartTime_Absolute_ITI,SpikeTime,StepSize,BinWidth);

PSTH_SacOn_ITI=AlignedSpike.PSTH_Aligned;
TimeSequence_SacOn_ITI=AlignedSpike.PSTH_Time_Aligned;

%Select interval
SelectPSTHinterval_SacOnITI=[-PlotSacPreBuffer,PlotSacPostBuffer];
[PSTH_SacOn_ITI,TimeSequence_SacOn_ITI]=SelectPSTHInterval(PSTH_SacOn_ITI,TimeSequence_SacOn_ITI,SelectPSTHinterval_SacOnITI);

%Reorganize saccade amplitude
SaccadeAmplitude_ITI_Line=reshape(SaccadeAmplitude_ITI,numel(SaccadeAmplitude_ITI),1);
NaNScreener=~isnan(SaccadeAmplitude_ITI_Line);

SaccadeAmplitude_ITI_Line=SaccadeAmplitude_ITI_Line(NaNScreener);

PSTH_SacOn_ITI=PSTH_SacOn_ITI(NaNScreener,:);
TimeSequence_SacOn_ITI=TimeSequence_SacOn_ITI(NaNScreener,:);

TimeSequence_SacOn_ITI_mean=nanmean(TimeSequence_SacOn_ITI,1);
PSTH_SacOn_ITI_mean=nanmean(PSTH_SacOn_ITI,1);
PSTH_SacOn_ITI_sem=nanstd(PSTH_SacOn_ITI,[],1)./sqrt(sum(~isnan(PSTH_SacOn_ITI),1));


SaccadeStartPoint_ITI_X_Line=reshape(SaccadeStartPoint_ITI_X,numel(SaccadeStartPoint_ITI_X),1);
SaccadeStartPoint_ITI_X_Line=SaccadeStartPoint_ITI_X_Line(NaNScreener,:);

SaccadeStartPoint_ITI_Y_Line=reshape(SaccadeStartPoint_ITI_Y,numel(SaccadeStartPoint_ITI_Y),1);
SaccadeStartPoint_ITI_Y_Line=SaccadeStartPoint_ITI_Y_Line(NaNScreener,:);


%Eye trace
%Eccentricity
Eye_Eccentricity_ITI_Aligned_SacOn=AlignEyeDataToSaccades(Eye_Eccentricity_ITI,EyeTime_ITI,SaccadeStartTime_ITI,[-PlotSacPreBuffer,PlotSacPostBuffer]);
Eye_Eccentricity_ITI_Aligned_SacOn=Eye_Eccentricity_ITI_Aligned_SacOn(NaNScreener,:);

Eye_Velocity_ITI_Aligned_SacOn=abs(diff(Eye_Eccentricity_ITI_Aligned_SacOn,1,2));
Eye_Velocity_ITI_Aligned_SacOn=[Eye_Velocity_ITI_Aligned_SacOn,Eye_Velocity_ITI_Aligned_SacOn(:,end)]/EyeBinWidth*1000;

%Velocity
Eye_Velocity_ITI_Aligned_SacOn_mean=nanmean(Eye_Velocity_ITI_Aligned_SacOn,1);
Eye_Velocity_ITI_Aligned_SacOn_sem=nanstd(Eye_Velocity_ITI_Aligned_SacOn,[],1)./sqrt(sum(~isnan(Eye_Velocity_ITI_Aligned_SacOn),1));

%Relative Eccentricity from the starting point
EyeChannel_X_ITI_Aligned_SacOn=AlignEyeDataToSaccades(EyeChannel_X_ITI,EyeTime_ITI,SaccadeStartTime_ITI,[-PlotSacPreBuffer,PlotSacPostBuffer]);
EyeChannel_X_ITI_Aligned_SacOn=EyeChannel_X_ITI_Aligned_SacOn(NaNScreener,:);

EyeChannel_Y_ITI_Aligned_SacOn=AlignEyeDataToSaccades(EyeChannel_Y_ITI,EyeTime_ITI,SaccadeStartTime_ITI,[-PlotSacPreBuffer,PlotSacPostBuffer]);
EyeChannel_Y_ITI_Aligned_SacOn=EyeChannel_Y_ITI_Aligned_SacOn(NaNScreener,:);



EyeChannel_EccRel_ITI_Aligned_SacOn=sqrt((EyeChannel_X_ITI_Aligned_SacOn-SaccadeStartPoint_ITI_X_Line).^2+(EyeChannel_Y_ITI_Aligned_SacOn-SaccadeStartPoint_ITI_Y_Line).^2);
EyeChannel_EccRel_ITI_Aligned_SacOn_mean=nanmean(EyeChannel_EccRel_ITI_Aligned_SacOn,1);
EyeChannel_EccRel_ITI_Aligned_SacOn_sem=nanstd(EyeChannel_EccRel_ITI_Aligned_SacOn,[],1)./sqrt(sum(~isnan(EyeChannel_EccRel_ITI_Aligned_SacOn),1));


%I. Sort neural responses according to the saccade amplitude

GroupStr={'Highest','Medium','Lowest'};
%GroupStr={'Above30','30-20','20-10','10-5','Below 5'};

SepGroupAmp=3;
SepGroupAmpNum= SepGroupAmp;
%SepGroupAmp={[20,10,1]};
%SepGroupAmpNum=length(SepGroupAmp{1})+1;


SacOnITI_temp=SortDataAccordingToGroup(PSTH_SacOn_ITI,SaccadeAmplitude_ITI_Line,SepGroupAmp,'descend');


PSTH_SacOn_ITI_AmpGroup_mean=SacOnITI_temp.mean;
PSTH_SacOn_ITI_AmpGroup_sem=SacOnITI_temp.sem;

Amp_Group_mean=SacOnITI_temp.MarkerMean;
Amp_Group_sem=SacOnITI_temp.MarkerSem;

NumInGroupAmp=SacOnITI_temp.Num;
Amp_SortIndex=SacOnITI_temp.SortIndex;
%Eye trace
SacOn_Eyetemp=SortDataAccordingToGroup(EyeChannel_EccRel_ITI_Aligned_SacOn,SaccadeAmplitude_ITI_Line,SepGroupAmp,'descend');


EyeEccRel_SacOn_AmpGroup_mean=SacOn_Eyetemp.mean;
EyeEccRel_SacOn_AmpGroup_sem=SacOn_Eyetemp.sem;


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Align spike time with saccade offset %%%%%

SaccadeEndTime_Absolute_ITI=SaccadeEndTime_ITI+EyeSelectionStart_Time;
AlignedSpike=SpikeAlignTool(SaccadeEndTime_Absolute_ITI,SpikeTime,StepSize,BinWidth);

PSTH_SacOff_ITI=AlignedSpike.PSTH_Aligned;
TimeSequence_SacOff_ITI=AlignedSpike.PSTH_Time_Aligned;

%Select interval
SelectPSTHinterval_SacOffITI=[-PlotSacPreBuffer,PlotSacPostBuffer];
[PSTH_SacOff_ITI,TimeSequence_SacOff_ITI]=SelectPSTHInterval(PSTH_SacOff_ITI,TimeSequence_SacOff_ITI,SelectPSTHinterval_SacOffITI);


PSTH_SacOff_ITI=PSTH_SacOff_ITI(NaNScreener,:);
TimeSequence_SacOff_ITI=TimeSequence_SacOff_ITI(NaNScreener,:);

TimeSequence_SacOff_ITI_mean=nanmean(TimeSequence_SacOff_ITI,1);
PSTH_SacOff_ITI_mean=nanmean(PSTH_SacOff_ITI,1);
PSTH_SacOff_ITI_sem=nanstd(PSTH_SacOff_ITI,[],1)./sqrt(sum(~isnan(PSTH_SacOff_ITI),1));



%Eye trace
%Eccentricity
Eye_Eccentricity_ITI_Aligned_SacOff=AlignEyeDataToSaccades(Eye_Eccentricity_ITI,EyeTime_ITI,SaccadeEndTime_ITI,[-PlotSacPreBuffer,PlotSacPostBuffer]);
Eye_Eccentricity_ITI_Aligned_SacOff=Eye_Eccentricity_ITI_Aligned_SacOff(NaNScreener,:);

Eye_Velocity_ITI_Aligned_SacOff=abs(diff(Eye_Eccentricity_ITI_Aligned_SacOff,1,2));
Eye_Velocity_ITI_Aligned_SacOff=[Eye_Velocity_ITI_Aligned_SacOff,Eye_Velocity_ITI_Aligned_SacOff(:,end)]/EyeBinWidth*1000;

%Velocity
Eye_Velocity_ITI_Aligned_SacOff_mean=nanmean(Eye_Velocity_ITI_Aligned_SacOff,1);
Eye_Velocity_ITI_Aligned_SacOff_sem=nanstd(Eye_Velocity_ITI_Aligned_SacOff,[],1)./sqrt(sum(~isnan(Eye_Velocity_ITI_Aligned_SacOff),1));

%Relative Eccentricity from the starting point
EyeChannel_X_ITI_Aligned_SacOff=AlignEyeDataToSaccades(EyeChannel_X_ITI,EyeTime_ITI,SaccadeEndTime_ITI,[-PlotSacPreBuffer,PlotSacPostBuffer]);
EyeChannel_X_ITI_Aligned_SacOff=EyeChannel_X_ITI_Aligned_SacOff(NaNScreener,:);

EyeChannel_Y_ITI_Aligned_SacOff=AlignEyeDataToSaccades(EyeChannel_Y_ITI,EyeTime_ITI,SaccadeEndTime_ITI,[-PlotSacPreBuffer,PlotSacPostBuffer]);
EyeChannel_Y_ITI_Aligned_SacOff=EyeChannel_Y_ITI_Aligned_SacOff(NaNScreener,:);
%{
SaccadeEndPoint_ITI_X_Line=EyeChannel_X_ITI_Aligned_SacOff(:,ceil(PlotSacPreBuffer/EyeBinWidth)+1);
SaccadeEndPoint_ITI_Y_Line=EyeChannel_Y_ITI_Aligned_SacOff(:,ceil(PlotSacPreBuffer/EyeBinWidth)+1);
%}
SaccadeEndPoint_ITI_X_Line=reshape(SaccadeEndPoint_ITI_X,numel(SaccadeEndPoint_ITI_X),1);
SaccadeEndPoint_ITI_X_Line=SaccadeEndPoint_ITI_X_Line(NaNScreener,:);

SaccadeEndPoint_ITI_Y_Line=reshape(SaccadeEndPoint_ITI_Y,numel(SaccadeEndPoint_ITI_Y),1);
SaccadeEndPoint_ITI_Y_Line=SaccadeEndPoint_ITI_Y_Line(NaNScreener,:);
%}


EyeChannel_EccRel_ITI_Aligned_SacOff=sqrt((EyeChannel_X_ITI_Aligned_SacOff-SaccadeEndPoint_ITI_X_Line).^2+(EyeChannel_Y_ITI_Aligned_SacOff-SaccadeEndPoint_ITI_Y_Line).^2);

%EyeChannel_EccRel_ITI_Aligned_SacOff=sqrt((EyeChannel_X_ITI_Aligned_SacOff-SaccadeStartPoint_ITI_X_Line).^2+(EyeChannel_Y_ITI_Aligned_SacOff-SaccadeStartPoint_ITI_Y_Line).^2);



EyeChannel_EccRel_ITI_Aligned_SacOff_mean=nanmean(EyeChannel_EccRel_ITI_Aligned_SacOff,1);
EyeChannel_EccRel_ITI_Aligned_SacOff_sem=nanstd(EyeChannel_EccRel_ITI_Aligned_SacOff,[],1)./sqrt(sum(~isnan(EyeChannel_EccRel_ITI_Aligned_SacOff),1));

%II. Sort neural responses according to the saccade amplitude


SacOffITI_temp=SortDataAccordingToGroup(PSTH_SacOff_ITI,SaccadeAmplitude_ITI_Line,SepGroupAmp,'descend');


PSTH_SacOff_ITI_AmpGroup_mean=SacOffITI_temp.mean;
PSTH_SacOff_ITI_AmpGroup_sem=SacOffITI_temp.sem;

%Eye trace
SacOff_Eyetemp=SortDataAccordingToGroup(EyeChannel_EccRel_ITI_Aligned_SacOff,SaccadeAmplitude_ITI_Line,SepGroupAmp,'descend');


EyeEccRel_SacOff_AmpGroup_mean=SacOff_Eyetemp.mean;
EyeEccRel_SacOff_AmpGroup_sem=SacOff_Eyetemp.sem;

%%
%Saccade Vector Tuning
%Get average PreSac, DurSac,PostSac spikes
%Use spike rate

EyeSelectionStart_Time_Line=reshape(SaccadeStartTime_ITI,numel(SaccadeStartTime_ITI),1);
EyeSelectionStart_Time_Line=EyeSelectionStart_Time_Line(NaNScreener);

EyeSelectionEnd_Time_Line=reshape(SaccadeEndTime_ITI,numel(SaccadeEndTime_ITI),1);
EyeSelectionEnd_Time_Line=EyeSelectionEnd_Time_Line(NaNScreener);


%PreSaccade
%PreSacDur=[-PlotSacPreBuffer,0];%in ms
PreSacDur=[-200,0];%in ms
SpikeRate_PreSac=nanmean(PSTH_SacOn_ITI(:,[TimeSequence_SacOn_ITI_mean<=PreSacDur(2) & TimeSequence_SacOn_ITI_mean>=PreSacDur(1)]),2);

%During Saccade:
SacIn=[0*ones(length(EyeSelectionStart_Time_Line),1),EyeSelectionEnd_Time_Line-EyeSelectionStart_Time_Line];
SpikeRate_DurSac=SelectFringRateInterval(PSTH_SacOn_ITI,TimeSequence_SacOn_ITI_mean,SacIn)';


%After Saccade:
%PostSacDur=[0,PlotSacPostBuffer];%in ms
PostSacDur=[0,200];%in ms
SpikeRate_PostSac=nanmean(PSTH_SacOff_ITI(:,[TimeSequence_SacOff_ITI_mean<=PostSacDur(2) &TimeSequence_SacOff_ITI_mean>=PostSacDur(1)]),2);

 %Group responses according to saccade vectors
 %SeparationPoint=[22.5:22.5:360];
 SeparationPoint=[22.5:45:360];
 
 SaccadeAngle_ITI_Line=reshape(SaccadeAngle_ITI,numel(SaccadeAngle_ITI),1);
 SaccadeAngle_ITI_Line=SaccadeAngle_ITI_Line(NaNScreener);
 
 SaccadeVectorTuning=PoloarHistGroup(SaccadeAngle_ITI_Line,[SpikeRate_PreSac,SpikeRate_DurSac,SpikeRate_PostSac],SeparationPoint);
 SaccadeVector=SaccadeVectorTuning.Center;
 FiringRate_PreDurPost_mean=SaccadeVectorTuning.DataMean;
 FiringRate_PreDurPost_sem=SaccadeVectorTuning.DataSEM;
 SacVecTuning_p=SaccadeVectorTuning.ANOVA_p;
 
 PrefVector=SaccadeVectorTuning.PrefVector;
PrefVector_Amp=SaccadeVectorTuning.PrefVector_Amp;

%PSTH on each direction
%figure
%plot(PSTH_SacOn_ITI');
AngleIn = [0:45:359];
for i =1: length(AngleIn)
    lw = AngleIn(i) - 22.5;
    hw = AngleIn(i) + 22.5;
    if lw < 0
        lw = lw + 360;
        sel = (SaccadeAngle_ITI_Line > lw & SaccadeAngle_ITI_Line < 360)|(SaccadeAngle_ITI_Line >= 0 & SaccadeAngle_ITI_Line <= hw);
    else
        sel = SaccadeAngle_ITI_Line > lw & SaccadeAngle_ITI_Line <=hw;
    end
    %AngleCenter(i) = (lw+hw)/2;
    PSTH_SacOn_sel = PSTH_SacOn_ITI(sel,:);
    PSTH_SacOn_Mean{i} = nanmean(PSTH_SacOn_sel,1);
    PSTH_SacOn_Sem{i} = nanstd(PSTH_SacOn_sel,[],1)./sqrt(sum(~isnan(PSTH_SacOn_sel),1));

    %Eye trace
    PSTH_EyeX{i} = nanmean(EyeChannel_X_ITI_Aligned_SacOn(sel,:)-nanmean(EyeChannel_X_ITI_Aligned_SacOn(sel,1:200),2));
    PSTH_EyeY{i} = nanmean(EyeChannel_Y_ITI_Aligned_SacOn(sel,:)-nanmean(EyeChannel_Y_ITI_Aligned_SacOn(sel,1:200),2));

        
end
%{
figure

c = {'r','g','b','c','m','y','k','c'};
for i =1: length(AngleIn)
    plot(PSTH_EyeX{i}',PSTH_EyeY{i}','or','color',c{i});
    hold on
end
xlim([-30,30])
ylim([-30,30])
%}
EyeTime =[-PlotSacPreBuffer:EyeBinWidth: PlotSacPostBuffer];

%%
if ShowFigureFlag
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Figure 1: Average firing rate
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%for each fractal
figtitlestr{1}='PSTHAlignOnSacOnsetDuringITI';

fig{1}=PrepareFigure(120,'w',[50,100, 1200,800],'Name',figtitlestr{1});
title('PSTH aligned on Saccade Onset during ITI','FontSize',25,'FontWeight','Bold');
axis off;
%%%%%%%%%Align on ITI Saccade Onset%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
axes('Position',[0.1,0.75,0.25,0.15]);

%Eye trace
shadedErrorBar([-PlotSacPreBuffer:EyeBinWidth: PlotSacPostBuffer],Eye_Velocity_ITI_Aligned_SacOn_mean,Eye_Velocity_ITI_Aligned_SacOn_sem,'lineprops',{[0.9290    0.6940    0.1250]},'transparent',1,'patchSaturation',0.3);
%shadedErrorBar([-PlotSacPreBuffer:EyeBinWidth: PlotSacPostBuffer],EyeChannel_EccRel_ITI_Aligned_SacOn_mean,EyeChannel_EccRel_ITI_Aligned_SacOn_sem,'lineprops',{[0.9290    0.6940    0.1250]},'transparent',1,'patchSaturation',0.3);


hold on
%shadedErrorBar([-PlotSacPreBuffer:EyeBinWidth: PlotSacPostBuffer],Eye_Velocity_ITI_Aligned_SacOn_mean,Eye_Velocity_ITI_Aligned_SacOn_sem,'lineprops',{'#8E929A'},'transparent',1,'patchSaturation',0.3);
xlim([-PlotSacPreBuffer, PlotSacPostBuffer]);%Keep in the same scale
yrange=ylim;
xticks([]);
xticklabels([]);
plot([0,0],yrange,'--k');
set(findall(gca, 'Type', 'Line'),'LineWidth',3);

set(gca,'FontSize',15,'FontWeight','Bold');
set(gca,'LineWidth',3);
%xlabel('Time from Saccade Onset(ITI)','FontSize',25,'FontWeight','Bold');
ylabel('Eye Velocity(deg/s)','FontSize',25,'FontWeight','Bold');
%ylabel('Eye Eccentricity(deg)','FontSize',25,'FontWeight','Bold');


axes('Position',[0.1,0.4,0.25,0.3]);
shadedErrorBar(TimeSequence_SacOn_ITI_mean,PSTH_SacOn_ITI_mean,PSTH_SacOn_ITI_sem,'lineprops',{[0.9290    0.6940    0.1250]},'transparent',1,'patchSaturation',0.3);
yrange=ylim;
hold on
xlim([-PlotSacPreBuffer,PlotSacPostBuffer]);
plot([0,0],yrange,'--k');
set(findall(gca, 'Type', 'Line'),'LineWidth',3);

set(gca,'FontSize',15,'FontWeight','Bold');
set(gca,'LineWidth',3);
xlabel('Time from Saccade Onset(ITI)','FontSize',20,'FontWeight','Bold');
ylabel('Spikes/s','FontSize',25,'FontWeight','Bold');



%%%%%%%%% PSTH aligned on sac on during ITI for different . Saccade Amplitudes

%Eye Trace
axes('Position',[0.4,0.75,0.25,0.15]);
for i=1:SepGroupAmpNum
    shadedErrorBar([-PlotSacPreBuffer:EyeBinWidth: PlotSacPostBuffer],EyeEccRel_SacOn_AmpGroup_mean(i,:),EyeEccRel_SacOn_AmpGroup_sem(i,:),'lineprops',{Color{i}},'transparent',1,'patchSaturation',0.3);  
    hold on

end
%ylim([0,18]);
yrange=ylim;
h=plot([0,0],yrange,'--k');
xlim([-PlotSacPreBuffer, PlotSacPostBuffer]);%Keep in the same scale
set(findall(gca, 'Type', 'Line'),'LineWidth',3);
xticks([]);
xticklabels([]);
set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
%ylabel('Eye Ecc(degree)','FontSize',25,'FontWeight','Bold');

%Neural firing
 axes('Position',[0.4,0.4,0.25,0.3]);
 for i=1:SepGroupAmpNum
    shadedErrorBar(TimeSequence_SacOn_ITI_mean,PSTH_SacOn_ITI_AmpGroup_mean(i,:),PSTH_SacOn_ITI_AmpGroup_sem(i,:),'lineprops',{Color{i}},'transparent',1,'patchSaturation',0.3);
    hold on
 end 

 xlim([-PlotSacPreBuffer,PlotSacPostBuffer]);
 yrange=ylim;
 h=plot([0,0],yrange,'--k');
  xlim([-PlotSacPreBuffer, PlotSacPostBuffer]);%Keep in the same scale
  set(findall(gca, 'Type', 'Line'),'LineWidth',3);   
  set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
  legend(GroupStr,'Location','northwest');
% ylabel('Firing Rate (Spikes/s)','FontSize',25,'FontWeight','Bold');
xlabel('Time from Saccade Onset(ms)ITI','FontSize',20,'FontWeight','Bold');

set(get(get(h,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');

%%%%%%%%% PSTH aligned on sac off during ITI for different . Saccade Amplitudes

%Eye Trace
axes('Position',[0.7,0.75,0.25,0.15]);
for i=1:SepGroupAmpNum
    shadedErrorBar([-PlotSacPreBuffer:EyeBinWidth: PlotSacPostBuffer],EyeEccRel_SacOff_AmpGroup_mean(i,:),EyeEccRel_SacOff_AmpGroup_sem(i,:),'lineprops',{Color{i}},'transparent',1,'patchSaturation',0.3);  
    hold on

end
%ylim([0,18]);
yrange=ylim;
h=plot([0,0],yrange,'--k');
xlim([-PlotSacPreBuffer, PlotSacPostBuffer]);%Keep in the same scale
set(findall(gca, 'Type', 'Line'),'LineWidth',3);
xticks([]);
xticklabels([]);
set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
%ylabel('Eye Ecc(degree)','FontSize',25,'FontWeight','Bold');

%Neural firing
 axes('Position',[0.7,0.4,0.25,0.3]);
 for i=1:SepGroupAmpNum
    shadedErrorBar(TimeSequence_SacOff_ITI_mean,PSTH_SacOff_ITI_AmpGroup_mean(i,:),PSTH_SacOff_ITI_AmpGroup_sem(i,:),'lineprops',{Color{i}},'transparent',1,'patchSaturation',0.3);
    hold on
 end 

 xlim([-PlotSacPreBuffer,PlotSacPostBuffer]);
 yrange=ylim;
 h=plot([0,0],yrange,'--k');
  xlim([-PlotSacPreBuffer, PlotSacPostBuffer]);%Keep in the same scale
  set(findall(gca, 'Type', 'Line'),'LineWidth',3);   
  set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
  legend(GroupStr,'Location','northwest');
% ylabel('Firing Rate (Spikes/s)','FontSize',25,'FontWeight','Bold');
xlabel('Time from Saccade Offset(ms)ITI','FontSize',20,'FontWeight','Bold');

set(get(get(h,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');



%%
%Polar Plot for saccade vector tuning


subtitlestr={'Pre-Sac Tuning','Dur-Sac Tuning','Post-Sac Tuning'};

Selection_temp=~isnan(SaccadeVector);
SaccadeVector=SaccadeVector(Selection_temp);
FiringRate_PreDurPost_mean=FiringRate_PreDurPost_mean(Selection_temp,:);
FiringRate_PreDurPost_sem=FiringRate_PreDurPost_sem(Selection_temp,:);


for j=1:3
%Pre-saccade , During, Post saccade tuning
%subplot(1,3,j);
axes('Position',[0.1+(j-1)*0.3,0.05,0.2,0.2]);
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
a=annotation('textbox',[0.25+(j-1)*0.3,0.015,0.1,0.1] ,'String',pstr,'FitBoxToText','on','LineStyle','none');
a.FontSize=30;
a.FontWeight='Bold';
if SacVecTuning_p(j)<0.05
    a.Color='r';
end

rlim([0,max(max(FiringRate_PreDurPost_mean))*1.2]);
end
end
%%
%{
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Figure 2: Firing Rate relative
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%to the eye position during ITI
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figtitlestr{2}='EyePos_Firing_ITI';

fig{2}=PrepareFigure(121,'w',[70,100, 800,800],'Name',figtitlestr{2});
title('Firing Relative To Eye Position with different delays','FontSize',25,'FontWeight','Bold');

lim_c_low=min(cellfun(@min,MeanPSTH_FR_ITI_Line_Sort_Plot));
lim_c_high=min(max(cellfun(@max,MeanPSTH_FR_ITI_Line_Sort_Plot)),200);

%
for jj=1:length(ShiftTimeLib)
subplot(ceil(length(ShiftTimeLib)/3),2,jj);
%s=scatter(EyeChannel_X_ITI_Desample_Line_Sort,EyeChannel_Y_ITI_Desample_Line_Sort,[],MeanPSTH_FR_ITI_Line_Sort,'filled');
s=scatter(EyeChannel_X_ITI_Desample_Line_Sort_Plot{jj},EyeChannel_Y_ITI_Desample_Line_Sort_Plot{jj},[],MeanPSTH_FR_ITI_Line_Sort_Plot{jj},'filled');
xlim([-50,50]);
ylim([-50,50]);
axis equal;
box on;
caxis([lim_c_low,lim_c_high]);

title( sprintf('Delay %d ms',ShiftTimeLib(jj)),'FontSize',20,'FontWeight','Bold');
end
%}

%% Align on saccade onset
FigNum =2;
FigIndex =121;

  figtitlestr{FigNum}='PSTH_EachAngle';
   fig{FigNum}=PrepareFigure(FigIndex,'w',[50,100, 1200,800],'Name',figtitlestr{FigNum});
    
    %Plot range 



   
    for j = 1:length(AngleIn)
        
        
          ax2(j) = subplot(ceil(length(AngleIn)/3),3,j);
         
        
       
       %  
          %Align on saccade onset
          
         shadedErrorBar(TimeSequence_SacOn_ITI_mean,PSTH_SacOn_Mean{j},PSTH_SacOn_Sem{j},'lineprops',{[0    0    1]},'transparent',1,'patchSaturation',0.3);
         hold on
          set(findall(gca, 'Type', 'Line'),'LineWidth',3);

        set(gca,'FontSize',15,'FontWeight','Bold');
        set(gca,'LineWidth',3);
        title(sprintf("Angle: %1.1d,Ecc: %1.1d",AngleIn(j)));
       
            xlabel('Time from Sac Onset(ms)','FontSize',15,'FontWeight','Bold');
            
        
        
        ylabel('Firing Rate(Hz)','FontSize',15,'FontWeight','Bold');
        
        
        
       %
        linkaxes(ax2,'xy');
    end
    
    %Add title

title_msg=sprintf(' %s, SpikeChannel: %d, TrialNum: %d ',FileName,SpikeChannel,TrialNum);
a=annotation('textbox',[0.25,0.7,0.3,0.3],'String',title_msg,'FitBoxToText','on');
a.FontSize=25;
a.FontWeight='Bold';
a.LineStyle='none';

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
    if isfield(OutputData,'SaccadeAnalysis')
    ExistFlag=1;
    end
end


%%%%%%%%%%%%%%%%%%%%OutputData For current analysis
OutputData_New=[];

OutputData_New.SaccadeAnalysis.TrialType='GoodTrials';
%%%%%%%%%Save overall scene response characteristics%%%%%%%%%%%%%%%%%%%%%%%
%Organize the data
NameStr={
    'PSTH_FR_SacOn_ITI_Time','PSTH_FR_SacOn_ITI',...
    'PSTH_FR_SacOff_ITI_Time','PSTH_FR_SacOff_ITI',...
    'SacVec','SacVectTuning_PreDurPost_mean','SacVectTuning_PreDurPost_sem','SacVectTuning_p','Pref_Vect','Pref_Vect_Length',...
    'SacOnLocMean','SacOnLocSem'
    
    
};

DataLib={
    TimeSequence_SacOn_ITI_mean,[PSTH_SacOn_ITI_mean',PSTH_SacOn_ITI_sem'],...
    TimeSequence_SacOff_ITI_mean,[PSTH_SacOff_ITI_mean',PSTH_SacOff_ITI_sem'],...
    SaccadeVector,FiringRate_PreDurPost_mean,FiringRate_PreDurPost_sem,SacVecTuning_p,PrefVector,PrefVector_Amp,...
    PSTH_SacOn_Mean,PSTH_SacOn_Sem;

    
};

%{
NameStr={'PSTH_EyeEccRel_ITI',...
    'PSTH_FR_SacOn_ITI_Time','PSTH_FR_SacOn_ITI',...
    'PSTH_EyeEccRel_SacOn_ITI_AmpGroup_Mean','PSTH_EyeEccRel_SacOn_ITI_AmpGroup_Sem',...
    'PSTH_FR_SacOff_ITI_Time','PSTH_EyeEccRel_SacOff_ITI_AmpGroup_Mean','PSTH_EyeEccRel_SacOff_ITI_AmpGroup_Sem',...
    'Amp_Goup_Mean_SEM','SepGroupAmp',...
    'PSTH_EyeEccRel_SacOff_ITI_AmpGroup_Mean','PSTH_EyeEccRel_SacOff_ITI_AmpGroup_Sem',...
     'SacVec','SacVectTuning_PreDurPost_mean','SacVectTuning_PreDurPost_sem','SacVectTuning_p','Pref_Vect','Pref_Vect_Length'
    
    
};

DataLib={[EyeChannel_EccRel_ITI_Aligned_SacOn_mean',EyeChannel_EccRel_ITI_Aligned_SacOn_sem'],...
    TimeSequence_SacOn_ITI_mean,[PSTH_SacOn_ITI_mean',PSTH_SacOn_ITI_sem'],...
    EyeEccRel_SacOn_AmpGroup_mean,EyeEccRel_SacOn_AmpGroup_sem,...
    TimeSequence_SacOff_ITI_mean,PSTH_SacOff_ITI_AmpGroup_mean,PSTH_SacOff_ITI_AmpGroup_sem,...
    [Amp_Group_mean,Amp_Group_sem],SepGroupAmp,...
    EyeEccRel_SacOff_AmpGroup_mean,EyeEccRel_SacOff_AmpGroup_sem,...
    SaccadeVector,FiringRate_PreDurPost_mean,FiringRate_PreDurPost_sem,SacVecTuning_p,PrefVector,PrefVector_Amp

    
};
%}

DataStamp=containers.Map(NameStr,DataLib);





OutputData_New.SaccadeAnalysis.Task=TaskType;
OutputData_New.SaccadeAnalysis.TaskCode=TaskCode;
OutputData_New.SaccadeAnalysis.DataStamp=DataStamp;


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
Task_Old=OutputData.SaccadeAnalysis.TaskCode;

     if Task_Old~=TaskCode
        %Not the same task, add the new dataset into the old one
        OutputData.SaccadeAnalysis(numel(OutputData.SaccadeAnalysis)+1)=OutputData_New.SaccadeAnalysis;
     else
         %If the same task,the same protocol,replace the old one with the new one
         
         OutputData.SaccadeAnalysis=OutputData_New.SaccadeAnalysis;
         
    
    
     end
else
    %Output the current file
   OutputData.SaccadeAnalysis=OutputData_New.SaccadeAnalysis; 
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
function FR=SelectFringRateInterval(FR_Raw,Time,Interval);
for i=1:size(Interval,1)
    if sum(~isnan(Interval(i,:)))==2
        FR(i)=nanmean(FR_Raw(i,Time>=Interval(i,1) & Time<=Interval(i,2)));
    else
         FR(i)=NaN;
    end
    
    
    
end




end