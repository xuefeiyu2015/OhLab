function FreeViewEyeMovement_ExampleExport(Data, StartTrial, EndTrial,ShowFigureFlag,OutputFlag);


ECode;

TaskType=Data.TaskType;
DataPath= Data.Path;
FileName=Data.FileName;
TaskCode=Data.TaskCode;


if FileName(1)=='R'
    Monkey = 1;
else
    Monkey = 2;
end


Selection=StartTrial:EndTrial;%Select the trials according to the pannel

%Criterias for sexwlecting saccades
V_Threshold=40;
ContinueBin=5;
ISI_Threshold=20;

%User handler
PlotExample = 1;%For plotting the raw for example, not needed for every session
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


%Select the eyetrace during the "fixation"
TaskStartBin=FindOutTime(EventChannel,EventBin,SHOWFIXCD);
TaskEndBin=FindOutTime(EventChannel,EventBin,FIXOFF);

EyeChannel_X_Task= SelectEyeInterval(EyeChannel_X_L,TaskStartBin,TaskEndBin)';
EyeChannel_Y_Task= SelectEyeInterval(EyeChannel_Y_L,TaskStartBin,TaskEndBin)';

EyeChannel_X_Task=ReorganizeEye(EyeChannel_X_Task);
EyeChannel_Y_Task=ReorganizeEye(EyeChannel_Y_Task);

Eye_Eccentricity_Task=SelectEyeInterval(Eye_Eccentricity,TaskStartBin,TaskEndBin)';
Eye_Eccentricity_Task=ReorganizeEye(Eye_Eccentricity_Task);

EyeTime_Task=0:EyeBinWidth:size(Eye_Eccentricity_Task,2)-EyeBinWidth;



%Select the eye trace during the 100pre-stim, 100 stim, 100 post stim
%period


%StimStartBin=FindOutTime(EventChannel,EventBin,PRE_STIM);
StimStartBin=FindOutTime(EventChannel,EventBin,ON_STIM )-100;
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



%% SAVE THE SPIKE CHANNEL
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
    






EyeChannel_X_StimAlign = AlignEyeData(EyeChannel_X_Stim,EyeTime_Stim,MarkerTime,Interval);
EyeChannel_Y_StimAlign = AlignEyeData(EyeChannel_Y_Stim,EyeTime_Stim,MarkerTime,Interval);




NormX = nanmean(EyeChannel_X_StimAlign(:,1:100),2);
NormY = nanmean(EyeChannel_Y_StimAlign(:,1:100),2);

EyeChannel_X_StimAlign_Norm = EyeChannel_X_StimAlign-NormX;
EyeChannel_Y_StimAlign_Norm = EyeChannel_Y_StimAlign-NormY;

EyeTime_Stim_Align=Interval(1):EyeBinWidth:Interval(2)+EyeBinWidth;
%Analysis saccades during the whole "fixation" period
%{
%Save an example file for illustration of the select saccade procedural
DataPath= Data.Path;
path=DataPath;

cd(path);
EyeX = EyeChannel_X_StimAlign;
EyeY = EyeChannel_Y_StimAlign;

save('EyeX','EyeX');
save('EyeY','EyeY');
%}
Saccades=SelectSaccade(EyeChannel_X_Stim,EyeChannel_Y_Stim,EyeBinWidth,V_Threshold,ContinueBin,ISI_Threshold);

%SaccadeNumber=arrayfun(@(x)  x.NumOfSaccade,Saccades)';
SaccadeAngle=arrayfun(@(x)  x.SaccadeAngle,Saccades,'UniformOutput' ,0)';
SaccadeAmplitude=arrayfun(@(x)  x.SaccadeAmplitude,Saccades,'UniformOutput' ,0)';

SaccadeStartTime_ori=arrayfun(@(x)  x.SaccadeStartTime,Saccades,'UniformOutput' ,0)';
SaccadeEndTime_ori=arrayfun(@(x)  x.SaccadeEndTime,Saccades,'UniformOutput' ,0)';

SaccadeStartPoint=arrayfun(@(x)  x.SaccadeStartPoint,Saccades,'UniformOutput' ,0)';
SaccadeEndPoint=arrayfun(@(x)  x.SaccadeEndPoint,Saccades,'UniformOutput' ,0)';


%Screen out the valid saccades and saccade start time severas

AmplitudeThreshold=1;
TimeThreshold=0;
%TimeThresholdEnd=num2cell(MaxTimeLength-PlotSacPostBuffer);
%{
SaccadeAngle_ITI=fillinNaNs(cellfun(@(x,y,z,m) z(((x>AmplitudeThreshold)&(y<=m)))',SaccadeAmplitude,SaccadeEndTime,SaccadeAngle,TimeThresholdEnd,'UniformOutput' ,0));
SaccadeEndPoint_ITI=(cellfun(@(x,y,z,m) z(((x>AmplitudeThreshold)&(y<=m)),:)',SaccadeAmplitude,SaccadeEndTime,SaccadeEndPoint,TimeThresholdEnd,'UniformOutput' ,0));
SaccadeStartPoint_ITI=(cellfun(@(x,y,z,m) z(((x>AmplitudeThreshold)&(y<=m)),:)',SaccadeAmplitude,SaccadeEndTime,SaccadeStartPoint,TimeThresholdEnd,'UniformOutput' ,0));
%}

SaccadeAngle=fillinNaNs(cellfun(@(x,y,z) z(((x>AmplitudeThreshold)&(y>TimeThreshold)))',SaccadeAmplitude,SaccadeStartTime_ori,SaccadeAngle,'UniformOutput' ,0));
SaccadeEndPoint=(cellfun(@(x,y,z) z(((x>AmplitudeThreshold)&(y>TimeThreshold)),:)',SaccadeAmplitude,SaccadeStartTime_ori,SaccadeEndPoint,'UniformOutput' ,0));
SaccadeStartPoint=(cellfun(@(x,y,z) z(((x>AmplitudeThreshold)&(y>TimeThreshold)),:)',SaccadeAmplitude,SaccadeStartTime_ori,SaccadeStartPoint,'UniformOutput' ,0));



SaccadeNumInEachTrial=cellfun(@(x) size(x,2),SaccadeEndPoint);

SaccadeEndPoint_X=NaN*ones(length(SaccadeNumInEachTrial),max(SaccadeNumInEachTrial));
SaccadeEndPoint_Y=NaN*ones(length(SaccadeNumInEachTrial),max(SaccadeNumInEachTrial));

SaccadeStartPoint_X=NaN*ones(length(SaccadeNumInEachTrial),max(SaccadeNumInEachTrial));
SaccadeStartPoint_Y=NaN*ones(length(SaccadeNumInEachTrial),max(SaccadeNumInEachTrial));


SaccadeEndPoint_X(SaccadeNumInEachTrial>0,:)=fillinNaNs(cellfun(@(x) x(1,:),SaccadeEndPoint(SaccadeNumInEachTrial>0,:),'UniformOutput' ,0));
SaccadeEndPoint_Y(SaccadeNumInEachTrial>0,:)=fillinNaNs(cellfun(@(x) x(2,:),SaccadeEndPoint(SaccadeNumInEachTrial>0,:),'UniformOutput' ,0));

SaccadeStartPoint_X(SaccadeNumInEachTrial>0,:)=fillinNaNs(cellfun(@(x) x(1,:),SaccadeStartPoint(SaccadeNumInEachTrial>0,:),'UniformOutput' ,0));
SaccadeStartPoint_Y(SaccadeNumInEachTrial>0,:)=fillinNaNs(cellfun(@(x) x(2,:),SaccadeStartPoint(SaccadeNumInEachTrial>0,:),'UniformOutput' ,0));

%{
SaccadeStartTime_ITI=fillinNaNs(cellfun(@(x,y,z,m) z(((x>AmplitudeThreshold)&(y<=m)),:)',SaccadeAmplitude,SaccadeStartTime,SaccadeStartTime,TimeThresholdEnd,'UniformOutput' ,0));
SaccadeEndTime_ITI=fillinNaNs(cellfun(@(x,y,z,m) z(((x>AmplitudeThreshold)&(y<=m)),:)',SaccadeAmplitude,SaccadeStartTime,SaccadeEndTime,TimeThresholdEnd,'UniformOutput' ,0));


SaccadeAmplitude_ITI=fillinNaNs(cellfun(@(x,y,z,m) z(((x>AmplitudeThreshold)&(y<=m)),:)',SaccadeAmplitude,SaccadeStartTime,SaccadeAmplitude,TimeThresholdEnd,'UniformOutput' ,0));
%}
%

SaccadeStartTime=fillinNaNs(cellfun(@(x,y,z) z(((x>AmplitudeThreshold)&(y>TimeThreshold)),:)',SaccadeAmplitude,SaccadeStartTime_ori,SaccadeStartTime_ori,'UniformOutput' ,0));
SaccadeEndTime=fillinNaNs(cellfun(@(x,y,z) z(((x>AmplitudeThreshold)&(y>TimeThreshold)),:)',SaccadeAmplitude,SaccadeStartTime_ori,SaccadeEndTime_ori,'UniformOutput' ,0));


SaccadeAmplitude=fillinNaNs(cellfun(@(x,y,z) z(((x>AmplitudeThreshold)&(y>TimeThreshold)),:)',SaccadeAmplitude,SaccadeStartTime_ori,SaccadeAmplitude,'UniformOutput' ,0));



%Seperate into stim trial and sham trial
StimType=ReproduceFromEvent(EventChannel,[STIM1,SHAM0])';

StimTrial = StimType == STIM1;
ControlTrial = ~StimTrial;

NumberStimTrial = sum(StimTrial);
NumberControlTrial = sum(ControlTrial);

%Normalized raw eye trace, separate into stim trials and control trials
EyeX_Stim = EyeChannel_X_StimAlign_Norm(StimTrial,:);
EyeY_Stim = EyeChannel_Y_StimAlign_Norm(StimTrial,:);

EyeX_Control = EyeChannel_X_StimAlign_Norm(ControlTrial,:);
EyeY_Control= EyeChannel_Y_StimAlign_Norm(ControlTrial,:);




%Saccade timing

%Align eye time relative to stimulation
SaccadeStartTime = SaccadeStartTime + Interval(1);
SaccadeTime_Stim =SaccadeStartTime(StimTrial,:);
SaccadeTime_Control =SaccadeStartTime(ControlTrial,:);

SaccadeTime_Stim=[SaccadeTime_Stim,NaN*ones(size(SaccadeTime_Stim,1),2)];
SaccadeTime_Control=[SaccadeTime_Control,NaN*ones(size(SaccadeTime_Control,1),2)];


SaccadeAmplitude_Stim = SaccadeAmplitude(StimTrial,:);
SaccadeAmplitude_Control = SaccadeAmplitude(ControlTrial,:);

SaccadeAmplitude_Stim=[SaccadeAmplitude_Stim,NaN*ones(size(SaccadeAmplitude_Stim,1),2)];
SaccadeAmplitude_Control=[SaccadeAmplitude_Control,NaN*ones(size(SaccadeAmplitude_Control,1),2)];



SaccadeTime_Stim_line = reshape(SaccadeTime_Stim,1,numel(SaccadeTime_Stim));
SaccadeTime_Stim_line = SaccadeTime_Stim_line(~isnan(SaccadeTime_Stim_line));

SaccadeTime_Control_line = reshape(SaccadeTime_Control,1,numel(SaccadeTime_Control));
SaccadeTime_Control_line = SaccadeTime_Control_line(~isnan(SaccadeTime_Control_line));

SaccadeAmplitude_Stim_line = reshape(SaccadeAmplitude_Stim,1,numel(SaccadeAmplitude_Stim));
SaccadeAmplitude_Stim_line = SaccadeAmplitude_Stim_line(~isnan(SaccadeAmplitude_Stim_line));

SaccadeAmplitude_Control_line = reshape(SaccadeAmplitude_Control,1,numel(SaccadeAmplitude_Control));
SaccadeAmplitude_Control_line = SaccadeAmplitude_Control_line(~isnan(SaccadeAmplitude_Control_line));



StepSize = 5;
BinWidth = 5;
SelectIntervalTime = repmat(Interval,size(SaccadeTime_Stim,1),1);

[SaccadeRate_Stim,TimeSequence,SaccadeCount_Stim]=CalculateSpikeFiringRate(SaccadeTime_Stim,StepSize,BinWidth,SelectIntervalTime);
SaccadeCountStim_Sum = nansum(SaccadeCount_Stim);

SelectIntervalTime = repmat(Interval,size(SaccadeTime_Control,1),1);

[SaccadeRate_Control,TimeSequence,SaccadeCount_Control]=CalculateSpikeFiringRate(SaccadeTime_Control,StepSize,BinWidth,SelectIntervalTime);
SaccadeCountControl_Sum = nansum(SaccadeCount_Control);

TimeSeq = nanmean(TimeSequence);


SaccadeFreqStim_Sum = SaccadeCountStim_Sum/NumberStimTrial/5*10^3; %in saccade/s
SaccadeFreqControl_Sum = SaccadeCountControl_Sum/NumberControlTrial/5*10^3; %in saccade/s


% Select saccades during stimulation


SaccadeSelection_Stim = SaccadeTime_Stim_line >= StimDur(1) & SaccadeTime_Stim_line <= StimDur(2);
SaccadeSelection_Control = SaccadeTime_Control_line >= StimDur(1) & SaccadeTime_Control_line <= StimDur(2);

%SaccadeAngle
%Transform to relative to contralateral side 
if Monkey==2
    SaccadeAngle=360-SaccadeAngle;
    disp('Transform to relative to contralateral site for Adams');
end



SaccadeAmp_Stim = SaccadeAngle(StimTrial,:);
SaccadeAmp_Stim = reshape(SaccadeAmp_Stim,1,numel(SaccadeAmp_Stim));
SaccadeAmp_Stim = SaccadeAmp_Stim(~isnan(SaccadeAmp_Stim));




SaccadeAmp_Control = SaccadeAngle(ControlTrial,:);

SaccadeAmp_Control = reshape(SaccadeAmp_Control,1,numel(SaccadeAmp_Control));
SaccadeAmp_Control = SaccadeAmp_Control(~isnan(SaccadeAmp_Control));

%Select out the saccades during opto-stimulation
SaccadeAmp_Stim = SaccadeAmp_Stim(SaccadeSelection_Stim);%Angle
SaccadeAmp_Control = SaccadeAmp_Control(SaccadeSelection_Control);%Angle


SaccadeAmplitude_Stim_line = SaccadeAmplitude_Stim_line(SaccadeSelection_Stim);%Amplitude
SaccadeAmplitude_Control_line = SaccadeAmplitude_Control_line(SaccadeSelection_Control);%Amplitude


SeparationPoint=[22.5:22.5:360];
 
SaccadeVectorTuning=PoloarHistGroup(SaccadeAmp_Stim,ones(size(SaccadeAmp_Stim,2),1),SeparationPoint);
SaccadeVector=SaccadeVectorTuning.Center;
NumberOfAngle_Stim=SaccadeVectorTuning.DataNum/sum(SaccadeVectorTuning.DataNum)*100;%Change into proportion of saccades
NumberOfAngle_Stim_Sem=SaccadeVectorTuning.DataSEM;
 
PrefVectorNum_Stim=SaccadeVectorTuning.PrefVectorNum;
PrefVector_AmpNum_Stim=SaccadeVectorTuning.PrefVector_AmpNum;

SaccadeVectorTuning=PoloarHistGroup(SaccadeAmp_Control,ones(size(SaccadeAmp_Control,2),1),SeparationPoint);

NumberOfAngle_Control=SaccadeVectorTuning.DataNum/sum(SaccadeVectorTuning.DataNum)*100;%Change into proportion of saccades
 
PrefVectorNum_Control=SaccadeVectorTuning.PrefVectorNum;
PrefVector_AmpNum_Control=SaccadeVectorTuning.PrefVector_AmpNum;
NumberOfAngle_Control_Sem=SaccadeVectorTuning.DataSEM;



%Raw saccade vector grill plot

if Monkey==2
   SaccadeStartPoint_X=-SaccadeStartPoint_X;
SaccadeEndPoint_X=-SaccadeEndPoint_X;
    disp('reverse the sign of X for Adams');
end

SaccadeStartPoint_X_Stim=SaccadeStartPoint_X(StimTrial,:);
SaccadeStartPoint_Y_Stim=SaccadeStartPoint_Y(StimTrial,:);

SaccadeStartPoint_X_Stim_Line = reshape(SaccadeStartPoint_X_Stim,[],1);
SaccadeStartPoint_Y_Stim_Line = reshape(SaccadeStartPoint_Y_Stim,[],1);

SaccadeStartPoint_X_Control=SaccadeStartPoint_X(ControlTrial,:);
SaccadeStartPoint_Y_Control=SaccadeStartPoint_Y(ControlTrial,:);

SaccadeStartPoint_X_Control_Line = reshape(SaccadeStartPoint_X_Control,[],1);
SaccadeStartPoint_Y_Control_Line = reshape(SaccadeStartPoint_Y_Control,[],1);


SaccadeEndPoint_X_Stim=SaccadeEndPoint_X(StimTrial,:);
SaccadeEndPoint_Y_Stim=SaccadeEndPoint_Y(StimTrial,:);

SaccadeEndPoint_X_Stim_Line = reshape(SaccadeEndPoint_X_Stim,[],1);
SaccadeEndPoint_Y_Stim_Line = reshape(SaccadeEndPoint_Y_Stim,[],1);

SaccadeEndPoint_X_Control=SaccadeEndPoint_X(ControlTrial,:);
SaccadeEndPoint_Y_Control=SaccadeEndPoint_Y(ControlTrial,:);

SaccadeEndPoint_X_Control_Line = reshape(SaccadeEndPoint_X_Control,[],1);
SaccadeEndPoint_Y_Control_Line = reshape(SaccadeEndPoint_Y_Control,[],1);


%Transform and select out the raw saccade vectors during opto-stimulation


SaccadeStartPoint_X_Stim_Line = SaccadeStartPoint_X_Stim_Line(~isnan(SaccadeStartPoint_X_Stim_Line));
SaccadeStartPoint_Y_Stim_Line = SaccadeStartPoint_Y_Stim_Line(~isnan(SaccadeStartPoint_Y_Stim_Line));

SaccadeStartPoint_X_Control_Line = SaccadeStartPoint_X_Control_Line(~isnan(SaccadeStartPoint_X_Control_Line));
SaccadeStartPoint_Y_Control_Line = SaccadeStartPoint_Y_Control_Line(~isnan(SaccadeStartPoint_Y_Control_Line));

SaccadeStartPoint_X_Stim_Line = SaccadeStartPoint_X_Stim_Line(SaccadeSelection_Stim);
SaccadeStartPoint_Y_Stim_Line = SaccadeStartPoint_Y_Stim_Line(SaccadeSelection_Stim);

SaccadeStartPoint_X_Control_Line = SaccadeStartPoint_X_Control_Line(SaccadeSelection_Control);
SaccadeStartPoint_Y_Control_Line = SaccadeStartPoint_Y_Control_Line(SaccadeSelection_Control);


SaccadeEndPoint_X_Stim_Line = SaccadeEndPoint_X_Stim_Line(~isnan(SaccadeEndPoint_X_Stim_Line));
SaccadeEndPoint_Y_Stim_Line = SaccadeEndPoint_Y_Stim_Line(~isnan(SaccadeEndPoint_Y_Stim_Line));

SaccadeEndPoint_X_Control_Line = SaccadeEndPoint_X_Control_Line(~isnan(SaccadeEndPoint_X_Control_Line));
SaccadeEndPoint_Y_Control_Line = SaccadeEndPoint_Y_Control_Line(~isnan(SaccadeEndPoint_Y_Control_Line));

SaccadeEndPoint_X_Stim_Line = SaccadeEndPoint_X_Stim_Line(SaccadeSelection_Stim);
SaccadeEndPoint_Y_Stim_Line = SaccadeEndPoint_Y_Stim_Line(SaccadeSelection_Stim);

SaccadeEndPoint_X_Control_Line = SaccadeEndPoint_X_Control_Line(SaccadeSelection_Control);
SaccadeEndPoint_Y_Control_Line = SaccadeEndPoint_Y_Control_Line(SaccadeSelection_Control);





X_Comp_Control = SaccadeEndPoint_X_Control_Line - SaccadeStartPoint_X_Control_Line;
Y_Comp_Control = SaccadeEndPoint_Y_Control_Line - SaccadeStartPoint_Y_Control_Line;

X_Comp_Stim = SaccadeEndPoint_X_Stim_Line - SaccadeStartPoint_X_Stim_Line;
Y_Comp_Stim = SaccadeEndPoint_Y_Stim_Line - SaccadeStartPoint_Y_Stim_Line;









%Sorting for the arrow plot

 [val,Index_Control]= sort(Y_Comp_Control,'ascend');
 [val,Index_Stim]= sort(Y_Comp_Stim,'ascend');


X_Comp_Control = X_Comp_Control(Index_Control);
Y_Comp_Control = Y_Comp_Control(Index_Control);

X_Comp_Stim = X_Comp_Stim(Index_Stim);
Y_Comp_Stim = Y_Comp_Stim(Index_Stim);








%Compare the contrallateral saccade with ipsilateral saccades
clSAC = length(SaccadeAmp_Stim(SaccadeAmp_Stim > 0 & SaccadeAmp_Stim < 180));
ilSAC = length(SaccadeAmp_Stim(SaccadeAmp_Stim > 181 & SaccadeAmp_Stim < 359));

%
PropClSac = clSAC/(clSAC+ilSAC);
PropIlSac = 1-PropClSac;
%}
[chi_squared, df, p_prop, is_significant] = chisq_prop_test(clSAC, ilSAC, (clSAC + ilSAC)/2, (clSAC + ilSAC)/2);



%%Plotting
if ShowFigureFlag
    %{
    %Figure 1: Saccade distribution among aong the whole trial periods
    FigureStartNum=20;
    FigureIndex=1;
    figtitlestr{FigureIndex}='Distribution_Of_saccade_angle';

    fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[50,100, 600,600],'Name',figtitlestr{FigureIndex});
    %Control
    polarwitherrorbar([SaccadeVector/180*pi SaccadeVector(1)/180*pi],[NumberOfAngle_Control NumberOfAngle_Control(1)],[NumberOfAngle_Control_Sem' NumberOfAngle_Control_Sem(1)],...
    '-k',3);
hold on
h_pc = polarplot([PrefVectorNum_Control PrefVectorNum_Control]/180*pi,[0 PrefVector_AmpNum_Control],'b');

%stim

 polarwitherrorbar([SaccadeVector/180*pi SaccadeVector(1)/180*pi],[NumberOfAngle_Stim NumberOfAngle_Stim(1)],[NumberOfAngle_Stim_Sem' NumberOfAngle_Stim_Sem(1)],...
    '-r',3);
hold on
h_p = polarplot([PrefVectorNum_Stim PrefVectorNum_Stim]/180*pi,[0 PrefVector_AmpNum_Stim],'r');




    set(h_p,'LineWidth',2.5);
    set(h_pc,'LineWidth',2.5);

set(gca,'ThetaDir' , 'counterclockwise');
set(gca,'ThetaZeroLocation','top')
set(gca,'FontSize',25,'FontWeight','Bold','LineWidth',3);

thetaticks([0:45:359]);
thetaticklabels([0:45:359]);

%thetaticks(sort(SaccadeVector));
%thetaticklabels(sort(SaccadeVector));

title("Saccade Distribution",'FontSize',25,'FontWeight','Bold');


rlim([0,max(max([NumberOfAngle_Stim,NumberOfAngle_Control])*1.2)]);
    %}

%{
%Figure 2: Raw eye trace during the stimulated period
FigureIndex=FigureIndex+1;
FigureStartNum=FigureStartNum+1;

figtitlestr{FigureIndex}='RawEyetraceDuringStimulatedPeriod';
fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[70,100, 600,600],'Name',figtitlestr{FigureIndex});

EyeScalePlot = [-20,20];


subplot(2,2,1);

AverageRegion = StimDur;
area(AverageRegion ,[EyeScalePlot(2),EyeScalePlot(2)],EyeScalePlot(1),...
    'FaceColor','k','FaceAlpha',.1,'EdgeAlpha',.1,'HandleVisibility','off');
hold on


plot(repmat(EyeTime_Stim_Align(1:size(EyeX_Stim,2)),size(EyeX_Stim,1),1)',EyeX_Stim','-r','LineWidth',2);
ylim(EyeScalePlot);
xlim(Interval);
title('Stim');
ylabel('Normed X Dev(deg)');
set(gca,'LineWidth',3,'FontSize',20,'FontWeight','Bold');
box off

subplot(2,2,3);
AverageRegion = StimDur;
area(AverageRegion ,[EyeScalePlot(2),EyeScalePlot(2)],EyeScalePlot(1),...
    'FaceColor','k','FaceAlpha',.1,'EdgeAlpha',.1,'HandleVisibility','off');
hold on


plot(repmat(EyeTime_Stim_Align(1:size(EyeY_Stim,2)),size(EyeY_Stim,1),1)',EyeY_Stim','-r','LineWidth',2);
ylim(EyeScalePlot);
xlim(Interval);

ylabel('Normed Y Dev(deg)');
xlabel('Time from stimOn');
set(gca,'LineWidth',3,'FontSize',20,'FontWeight','Bold');
box off


subplot(2,2,2);
AverageRegion = StimDur;
area(AverageRegion ,[EyeScalePlot(2),EyeScalePlot(2)],EyeScalePlot(1),...
    'FaceColor','k','FaceAlpha',.1,'EdgeAlpha',.1,'HandleVisibility','off');
hold on

plot(repmat(EyeTime_Stim_Align(1:size(EyeX_Control,2)),size(EyeX_Control,1),1)',EyeX_Control','-k','LineWidth',2);
xlim(Interval);
ylim(EyeScalePlot);
title('Control');


set(gca,'LineWidth',3,'FontSize',20,'FontWeight','Bold');
box off

subplot(2,2,4);
AverageRegion = StimDur;
area(AverageRegion ,[EyeScalePlot(2),EyeScalePlot(2)],EyeScalePlot(1),...
    'FaceColor','k','FaceAlpha',.1,'EdgeAlpha',.1,'HandleVisibility','off');
hold on

plot(repmat(EyeTime_Stim_Align(1:size(EyeY_Control,2)),size(EyeY_Control,1),1)',EyeY_Control','-k','LineWidth',2);
ylim(EyeScalePlot);
xlabel('Time from shamOn');
set(gca,'LineWidth',3,'FontSize',20,'FontWeight','Bold');
box off
%}

%{]

%Figure 3: Grill plot illustration all the saccade vector
FigureIndex=1;
FigureStartNum=20;

figtitlestr{FigureIndex}='SaccadeVectorGrill';
fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[90,100, 600,600],'Name',figtitlestr{FigureIndex});
set(fig{FigureIndex}, 'PaperUnits', 'inches');
    formatFig(fig{FigureIndex}, [2 1], 'nature');
 
%subplot(1,2,1)

%quiver(SaccadeStartPoint_X_Stim_Line,SaccadeStartPoint_Y_Stim_Line,X_Comp_Stim,Y_Comp_Stim)
quiver(X_Comp_Stim,Y_Comp_Stim,'autoscale','off','color','r','LineWidth',0.5);
title('Opto-Stim');
%set(gca,'LineWidth',3,'FontSize',20,'FontWeight','Bold');
set(gca,'LineWidth',0.5,'FontSize',5);%,'FontWeight','Bold');
box off

%{
subplot(1,2,2)

%quiver(SaccadeStartPoint_X_Control_Line,SaccadeStartPoint_Y_Control_Line,X_Comp_Control,Y_Comp_Control)
quiver(X_Comp_Control,Y_Comp_Control,'autoscale','off','color','k','LineWidth',3)
title('Control')
%set(gca,'LineWidth',3,'FontSize',20,'FontWeight','Bold');
set(gca,'LineWidth',0.5,'FontSize',5);%,'FontWeight','Bold');
box off
%}
%}
%{
%Figure 4: Saccade eye trace illustrated on 2D screen
FigureIndex=1;
FigureStartNum=20;

figtitlestr{FigureIndex}='SaccadesOnScreen';
fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[100,100, 600,300],'Name',figtitlestr{FigureIndex});
subplot(1,2,1)
plot(EyeX_Stim',EyeY_Stim','-r','LineWidth',2);
title('Stim')
xlim(EyeScalePlot)
ylim(EyeScalePlot)


set(gca,'LineWidth',3,'FontSize',20,'FontWeight','Bold');

subplot(1,2,2)
plot(EyeX_Control',EyeY_Control','-k','LineWidth',2);
title('Control')
xlim(EyeScalePlot)
ylim(EyeScalePlot)

set(gca,'LineWidth',3,'FontSize',20,'FontWeight','Bold');
%}


%Figure 5: Raster for saccades
 FigureIndex = FigureIndex + 1;
FigureStartNum = FigureStartNum + 1;

%SaccadeCountStim_Sum = SaccadeFreqStim_Sum ;
%SaccadeCountControl_Sum = SaccadeFreqControl_Sum; 

figtitlestr{FigureIndex}='SaccadesRaster';
fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[100,100, 400,400],'Name',figtitlestr{FigureIndex});

set(fig{FigureIndex}, 'PaperUnits', 'inches');
    formatFig(fig{FigureIndex}, [1.7 1.7], 'nature');

StartLoc = 1;
LineWidth = 0.5;
Barlength = 1;

axes('Position',[0.15,0.83,0.8,0.16]);
AverageRegion = [0,nanmean(MarkerTime)];
area(AverageRegion,[size(MarkerTime,1)+5,size(MarkerTime,1)+5],...
    'FaceColor','k','FaceAlpha',.1,'EdgeAlpha',.1,'HandleVisibility','off');

hold on

%Raster For Stim
SpikeRaterPlot(SaccadeTime_Stim,StartLoc,Barlength,'r',LineWidth);

StartLoc = size(SaccadeTime_Stim ,1)+5;
%Raster For Control

SpikeRaterPlot(SaccadeTime_Control,StartLoc,Barlength,'k',LineWidth);
xlim(Interval);


 set(gca,'LineWidth',0.5,'FontSize',5);%,'FontWeight','Bold');
 yticks([]);
 box off
 axis off
 
%subplot(3,1,2);
axes('Position',[0.15,0.5,0.8,0.3]);

maxy = max(max(SaccadeCountStim_Sum),max(SaccadeCountControl_Sum))+2;
AverageRegion = [0,nanmean(MarkerTime)];
area(AverageRegion,[10,10],...
    'FaceColor','k','FaceAlpha',.1,'EdgeAlpha',.1,'HandleVisibility','off');

hold on

%{
plot(TimeSeq,SaccadeCountStim_Sum,'-r','LineWidth',3)
hold on
plot(TimeSeq,SaccadeCountControl_Sum,'-k','LineWidth',3)
%}
%bar(TimeSeq,SaccadeCountStim_Sum,'r')
stairs(TimeSeq,SaccadeCountStim_Sum,'r')
set(gca,'LineWidth',0.5,'FontSize',5);%,'FontWeight','Bold');
 
 box off
max_stimtime = TimeSeq(SaccadeCountStim_Sum == max(SaccadeCountStim_Sum));
maxstim = max(SaccadeCountStim_Sum)+3;
for i = 1:length(max_stimtime )
  plot([max_stimtime(i),max_stimtime(i)],[maxstim,maxstim],'v','MarkerSize',5,'Color','r','MarkerFaceColor','r');
end
ylim([0,10]);
xticks([]);

%legend({'Stim'});


%subplot(3,1,3);
axes('Position',[0.15,0.2,0.8,0.3]);
area(AverageRegion,[10,10],...
    'FaceColor','k','FaceAlpha',.1,'EdgeAlpha',.1,'HandleVisibility','off');
hold on
%bar(TimeSeq,SaccadeCountControl_Sum,'k')
stairs(TimeSeq,SaccadeCountControl_Sum,'k')



box off;
xlabel('Time from Opto-stim Onset (ms)'); 
ylabel('Number of saccades');

set(gca,'LineWidth',0.5,'FontSize',5);%,'FontWeight','Bold');

%legend({'Control'});

%% Put stim control together
FigureIndex=FigureIndex+1;
FigureStartNum=FigureStartNum+1;
figtitlestr{FigureIndex}='PSTH_Stim_Control_Together';
fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[120,100, 300,300],'Name',figtitlestr{FigureIndex});
set(fig{FigureIndex}, 'PaperUnits', 'inches');
formatFig(fig{FigureIndex}, [1.5 1.5], 'nature');
StartLoc = 1;
LineWidth = 0.5;
Barlength = 1;

%Raster
axes('Position',[0.15,0.6,0.8,0.3]);
AverageRegion = [0,nanmean(MarkerTime)];
area(AverageRegion,[size(MarkerTime,1)+5,size(MarkerTime,1)+5],...
    'FaceColor','k','FaceAlpha',.1,'EdgeAlpha',.1,'HandleVisibility','off');

hold on

%Raster For Stim
SpikeRaterPlot(SaccadeTime_Stim,StartLoc,Barlength,'r',LineWidth);

StartLoc = size(SaccadeTime_Stim ,1)+5;
%Raster For Control

SpikeRaterPlot(SaccadeTime_Control,StartLoc,Barlength,'k',LineWidth);
xlim(Interval);


 set(gca,'LineWidth',0.5,'FontSize',5);%,'FontWeight','Bold');
 yticks([]);
 box off
 axis off
 xlim([-100,200]);

%PSTH
axes('Position',[0.15,0.2,0.8,0.4]);

maxy = max(max(SaccadeCountStim_Sum),max(SaccadeCountControl_Sum))+2;
AverageRegion = [0,nanmean(MarkerTime)];
area(AverageRegion,[10,10],...
    'FaceColor','k','FaceAlpha',.1,'EdgeAlpha',.1,'HandleVisibility','off');

hold on

%{
plot(TimeSeq,SaccadeCountStim_Sum,'-r','LineWidth',3)
hold on
plot(TimeSeq,SaccadeCountControl_Sum,'-k','LineWidth',3)
%}
%bar(TimeSeq,SaccadeCountStim_Sum,'r')

stairs(TimeSeq,SaccadeCountControl_Sum,'k');
hold on
stairs(TimeSeq,SaccadeCountStim_Sum,'r');
set(gca,'LineWidth',0.5,'FontSize',5);%,'FontWeight','Bold');
 
 box off
max_stimtime = TimeSeq(SaccadeCountStim_Sum == max(SaccadeCountStim_Sum));
maxstim = max(SaccadeCountStim_Sum)+3;
for i = 1:length(max_stimtime )
  plot([max_stimtime(i),max_stimtime(i)],[maxstim,maxstim],'v','MarkerSize',5,'Color','r','MarkerFaceColor','r');
end
ylim([0,10]);
%xticks([]);
box off;
xlabel('Time from Opto-stim Onset (ms)'); 
ylabel('Number of saccades');

set(gca,'LineWidth',0.5,'FontSize',5);%,'FontWeight','Bold');
xlim([-100,200]);

%legend({'Stim'});


%% For example plot
FigureIndex=FigureIndex+1;
FigureStartNum=FigureStartNum+1;
    
    figtitlestr{FigureIndex}='SC_Eye_Vector_On_Screen_Stim';
    fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[120,100, 300,300],'Name',figtitlestr{FigureIndex});
    set(fig{FigureIndex}, 'PaperUnits', 'inches');
    formatFig(fig{FigureIndex}, [1 1], 'nature');

   % subplot(1,2,1)
    quiver(SaccadeStartPoint_X_Stim_Line,SaccadeStartPoint_Y_Stim_Line,X_Comp_Stim,Y_Comp_Stim,'AutoScale', 'off','Color','r','LineWidth',0.5);
    xlim([-40,40]);
    ylim([-40,40]);
    xticks([]);
    yticks([]);
   % axis equal;

 FigureIndex=FigureIndex+1;
FigureStartNum=FigureStartNum+1;   

figtitlestr{FigureIndex}='SC_Eye_Vector_On_Screen_Control';
    fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[120,100, 300,300],'Name',figtitlestr{FigureIndex});
    set(fig{FigureIndex}, 'PaperUnits', 'inches');
    formatFig(fig{FigureIndex}, [1 1], 'nature');


   % subplot(1,2,2)
    quiver(SaccadeStartPoint_X_Control_Line,SaccadeStartPoint_Y_Control_Line,X_Comp_Control,Y_Comp_Control,'AutoScale', 'off','Color','k','LineWidth',0.5);
    xlim([-40,40]);
    ylim([-40,40]);
    xticks([]);
    yticks([]);
   % axis equal;
    % Add scale bar
%set(gca,'Color','#C0C0C0');

%keyboard



%}
%{
%% Fig 6 Distribution of saccade amplitude 
FigureIndex=FigureIndex+1;
FigureStartNum=FigureStartNum+1;

figtitlestr{FigureIndex}='DistributionOfSaccadesAmplitude';
fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[100,100, 600,300],'Name',figtitlestr{FigureIndex});

%binEdges = min([SaccadeAmplitude_Stim_line,SaccadeAmplitude_Control_line])-2.5:2:max([SaccadeAmplitude_Stim_line,SaccadeAmplitude_Control_line])+2.5;
binEdges = -2.5:2:40+2.5;
counts_stim = histcounts(SaccadeAmplitude_Stim_line, binEdges);
counts_control = histcounts(SaccadeAmplitude_Control_line, binEdges);

binCenters = binEdges(1:end-1) + diff(binEdges) / 2; % Calculate bin centers for plotting

subplot(2,1,1);
bar(binCenters, counts_stim, 'FaceColor', 'r');
set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
%ylabel('Number of sessions');
box off
subplot(2,1,2);
bar(binCenters, counts_control, 'FaceColor', 'k');
set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
xlabel('Distribution of saccade amplitude (degree) ');
ylabel('Number of sessions');
box off
%keyboard
%}   
end %End of ifFigurePlot


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
%OutPath=strcat(BasicDirectory,'Results',BasicDirectory(end));%,OutputFolerName);
OutPath=strcat(BasicDirectory,'Results',BasicDirectory(end),'ExportFigure',BasicDirectory(end));%,OutputFolerName);
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
if TaskCode~=364
OutputFileName=sprintf('%s_%s_N%s_C%s_T%d.mat',MonkeyName,string(RecordDateOriginal),string(NeuronNum),string(ChannelNum),TaskCode);
OutputFigureName=sprintf('%s_%s_N%s_C%s_T%d',MonkeyName,string(RecordDateOriginal),string(NeuronNum),string(ChannelNum),TaskCode);
else

    OutputFileName=sprintf('%s_%s_N%s_C%s.mat',MonkeyName,string(RecordDateOriginal),string(NeuronNum),string(ChannelNum));
OutputFigureName=sprintf('%s_%s_N%s_C%s',MonkeyName,string(RecordDateOriginal),string(NeuronNum),string(ChannelNum));

end
if OutputFlag == 1
for jj=1:FigureIndex
        if ~isempty(fig{jj})
        

         OutputFigureName_Final=strcat(OutputFigureName,figtitlestr{jj},'.pdf');
         saveas(fig{jj},OutputFigureName_Final);
        end
        
            
end


    disp('Figures have been exported to the exported folder');
end
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

%{
ExistFlag=0;
%Load the old file if exist
if exist(OutputFileName)

    load(OutputFileName);
    %Load OutputData into the memory
    if isfield(OutputData,'FreeViewBehavior')
    ExistFlag=1;
    end
end


%%%%%%%%%%%%%%%%%%%%OutputData For current analysis
OutputData_New=[];

OutputData_New.FreeViewBehavior.TrialType='GoodTrials';
%%%%%%%%%Save overall scene response characteristics%%%%%%%%%%%%%%%%%%%%%%%
%Organize the data
NameStr={'PropClSac','p_prop','NumberOfAngle_Control','NumberOfAngle_Stim',...
    'TimeForSaccadeRaster','SaccadeCountStim','SaccadeCountControl'    
};

DataLib={PropClSac,p_prop,NumberOfAngle_Control,NumberOfAngle_Stim,...
    TimeSeq,SaccadeCountStim_Sum,SaccadeCountControl_Sum
};


DataStamp=containers.Map(NameStr,DataLib);





OutputData_New.FreeViewBehavior.Task=TaskType;
OutputData_New.FreeViewBehavior.TaskCode=TaskCode;
OutputData_New.FreeViewBehavior.DataStamp=DataStamp;


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
Task_Old=OutputData.FreeViewBehavior.TaskCode;

     if Task_Old~=TaskCode
        %Not the same task, add the new dataset into the old one
        OutputData.FreeViewBehavior(numel(OutputData.FreeViewBehavior)+1)=OutputData_New.FreeViewBehavior;
     else
         %If the same task,the same protocol,replace the old one with the new one
         
         OutputData.FreeViewBehavior=OutputData_New.FreeViewBehavior;
         
    
    
     end
else
    %Output the current file
   OutputData.FreeViewBehavior=OutputData_New.FreeViewBehavior; 
end



if OutputFlag
    cd(OutPath);
    save(OutputFileName,'OutputData');
    disp('Data Exported');
else
    disp('No data export');
end 
 %} 

end %End of the spike channel


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