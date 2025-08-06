function LaserMeasurement(Data, StartTrial, EndTrial,ShowFigureFlag,OutputFlag);

ECode;

TaskType=Data.TaskType;
DataPath= Data.Path;
FileName=Data.FileName;

Selection=StartTrial:EndTrial;%Select the trials according to the pannel

%Criterias for sexwlecting saccades
V_Threshold=40;
ContinueBin=5;
ISI_Threshold=20;

%Sequence information from the library
ParaLib=Data.ParaLib;

if ismember('stimPower',keys(ParaLib))
   stimPower = ParaLib('stimPower');
   stimPower=stimPower(:,2);
   stimPower=stimPower(Selection,:);
end

if ismember('preStimDur,StimDur,PostStimDur',keys(ParaLib))
   stimDurTheoretcial = ParaLib('preStimDur,StimDur,PostStimDur');
   stimDurTheoretcial=stimDurTheoretcial(:,3);
   stimDurTheoretcial=stimDurTheoretcial(Selection,:);
end

UniqueStimDur=unique(stimDurTheoretcial);

UniquePower =unique(stimPower);
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


StimStartBin=FindOutTime(EventChannel,EventBin,PRE_STIM);
StimEndBin=FindOutTime(EventChannel,EventBin,PST_STIM );

EyeChannel_X_Stim= SelectEyeInterval(EyeChannel_X_L,StimStartBin,StimEndBin)';
EyeChannel_Y_Stim= SelectEyeInterval(EyeChannel_Y_L,StimStartBin,StimEndBin)';

EyeChannel_X_Stim=ReorganizeEye(EyeChannel_X_Stim);
EyeChannel_Y_Stim=ReorganizeEye(EyeChannel_Y_Stim);

EyeTime_Stim = 0:EyeBinWidth:size(EyeChannel_X_Stim,2)-EyeBinWidth;




MarkerTime = FindOutTime(EventChannel,EventBin,ON_STIM)-StimStartBin;

Interval = [-100,nanmean(StimEndBin-FindOutTime(EventChannel,EventBin,ON_STIM))];
StimDur = [0,nanmean(FindOutTime(EventChannel,EventBin,OFF_STIM)- FindOutTime(EventChannel,EventBin,ON_STIM))];



EyeChannel_X_StimAlign = AlignEyeData(EyeChannel_X_Stim,EyeTime_Stim,MarkerTime,Interval);
EyeChannel_Y_StimAlign = AlignEyeData(EyeChannel_Y_Stim,EyeTime_Stim,MarkerTime,Interval);

NormX = nanmean(EyeChannel_X_StimAlign(:,1:100),2);
NormY = nanmean(EyeChannel_Y_StimAlign(:,1:100),2);

%EyeChannel_X_StimAlign_Norm = EyeChannel_X_StimAlign-NormX;
%EyeChannel_Y_StimAlign_Norm = EyeChannel_Y_StimAlign-NormY;

EyeChannel_X_StimAlign_Norm = EyeChannel_X_StimAlign;
EyeChannel_Y_StimAlign_Norm = EyeChannel_Y_StimAlign;

EyeTime_Stim_Align=Interval(1):EyeBinWidth:Interval(2)+EyeBinWidth;

%Seperate into stim trial and sham trial
StimType=ReproduceFromEvent(EventChannel,[STIM1,SHAM0])';

StimTrial = StimType == STIM1;
ControlTrial = ~StimTrial;


%Normalized raw eye trace, separate into stim trials and control trials
EyeX_Stim = EyeChannel_X_StimAlign_Norm(StimTrial,:);
EyeY_Stim = EyeChannel_Y_StimAlign_Norm(StimTrial,:);

EyeX_Control = EyeChannel_X_StimAlign_Norm(ControlTrial,:);
EyeY_Control= EyeChannel_Y_StimAlign_Norm(ControlTrial,:);


%%Plotting
if ShowFigureFlag
    
  FigureStartNum=20;
    FigureIndex=1;

figtitlestr{FigureIndex}='RawEyetraceDuringStimulatedPeriod';
fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[70,100, 400,300],'Name',figtitlestr{FigureIndex});

EyeScalePlot = [0,250];




AverageRegion = StimDur;
area(AverageRegion ,[EyeScalePlot(2),EyeScalePlot(2)],EyeScalePlot(1),...
    'FaceColor','k','FaceAlpha',.1,'EdgeAlpha',.1,'HandleVisibility','off');
hold on

%{
plot(repmat(EyeTime_Stim_Align(1:size(EyeX_Stim,2)),size(EyeX_Stim,1),1)',EyeX_Stim','-r','LineWidth',2);
ylim(EyeScalePlot);
xlim(Interval);
title('Stim');
ylabel('uW');
set(gca,'LineWidth',3,'FontSize',20,'FontWeight','Bold');
box off
%}



plot(repmat(EyeTime_Stim_Align(1:size(EyeY_Stim,2)),size(EyeY_Stim,1),1)',EyeY_Stim','-b','LineWidth',2);
ylim(EyeScalePlot);
xlim(Interval);
ylabel('uW');

title(sprintf("Laser Power : %d mW; Duration: %d ms",UniquePower,UniqueStimDur));



xlabel('Time from stimOn');
set(gca,'LineWidth',3,'FontSize',20,'FontWeight','Bold');
box off




end




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