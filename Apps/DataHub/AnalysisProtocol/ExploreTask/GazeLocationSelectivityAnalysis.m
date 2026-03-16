function GazeLocationSelectivityAnalysis(Data, StartTrial, EndTrial,ShowFigureFlag,OutputFlag);
%Analysis the gaze location selectivity for the location of extra fixation
%point (and others )
%by Xuefei Yu 10302020

ECode;
TaskType=Data.TaskType;
DataPath= Data.Path;
FileName=Data.FileName;

Selection=StartTrial:EndTrial;%Select the trials according to the pannel
%Abnormal check %Especially For Robin081120FORAGEXPLO4.rst
CheckCode=EXTRAFP;%This code should appear once in each trial

AbnormalTrials=AbnormalTrialCheck(Data,CheckCode);
Selection=Selection(~ismember(Selection,AbnormalTrials));

CheckCode=EXTRAFIXOFF;%This code should appear once in each trial

AbnormalTrials=AbnormalTrialCheck(Data,CheckCode);
Selection=Selection(~ismember(Selection,AbnormalTrials));

TrialNum=length(Selection);



%For segmentation of  the neural responses
%{
SegEvent=[FIXOFF,CUEON];%Event Flow of interest
SegEventEnd=[NaN,TGTCD];%Event end
SegEventStr={'BeforeExplore','CueOn'};
%}
%
%EXTRAFIXSTARTCD 
SegEvent=[EXTRAFIXSTARTCD ];%Event Flow of interest
%SegEvent=[EXTRAFP];%Event Flow of interest
SegEventEnd=[EXTRAFIXOFF];%Event end
SegEventStr={'ExtraFP'};
%}
PlotStartOffset=[200];%in ms
%PlotStartOffset=[1700];%in ms
PlotStopOffset=[200];%in ms




EventTypeFlag=[1];%1:Number same with the trial number, 2:Number same with the target number


%Spon Event
SponEvent=[SCENE_OFF];
SponPreTimeBuffer=-500;
SponEventEnd=NaN;
SponStartOffset=0;
SponStopOffset=1000;

%Parameter to define the sliding window to measure firing rates
StepSize=10;%in ms
BinWidth=20;%in ms



%GetEyeData
EyeData=Data.EyeDataRaw;
EyeChannel_X_L=EyeData(Selection,1);
EyeChannel_Y_L=EyeData(Selection,2);



EyeBinWidth=Data.EyeBinWidth;
EyeCalInfo=Data.EyeCalInfo;%Each row: X_Gain,X_Bias; Y_Gain; Y_bias; X=(raw+X_bias)*X_Gain
EyeCalInfo_L=EyeCalInfo(1,:);

%Load event data
EventChannel=Data.EventChannel(Selection,:);
EventTimeChannel=Data.EventTimeChannel(Selection,:);
EventBin=Data.EventBin(Selection,:);
%Get the information about the extra FP
%
ParaLib=Data.ParaLib;

if(ismember('ExtraFixPointID,Style,Size,Loc_X,Loc_Y',keys(ParaLib)))
ExtraFP_Info=ParaLib('ExtraFixPointID,Style,Size,Loc_X,Loc_Y');
ExtraFP_Loc=ExtraFP_Info(Selection,end-1:end);
[ExtraFP_Angle, ExtraFP_Amp]=LocToAngle(ExtraFP_Loc(:,1),ExtraFP_Loc(:,2));

else
    disp('No extra FP, change the analysis program');
    return 
    

end




%Load the spike channel
SpikeChannel=Data.SpikeChannel;
SpikeTime=squeeze(Data.SpikeTimeData(:,:,Selection))';

%{
%Temporally not used in this protocol
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
%}
%Segregate spike data from Fixation onset to target onset
SegmentData=EventSegTool(SegEvent,SegEventEnd,PlotStartOffset,PlotStopOffset,EventChannel,EventTimeChannel,SpikeTime,StepSize,BinWidth);
%Retrive data from the field
%Total average for each event
MeanFR=SegmentData.MeanFR; %[Mean,sem]%Averall mean and sem
MeanFRRaw=SegmentData.MeanFRRaw;%=AverageFiringRate_Event;

%PSTH(Spike Density Function) for each event
MeanPSTH=SegmentData.MeanPSTH;%=[PSTH_Event_Time',PSTH_Event'];
MeanPSTH_Time_ExtraFP=MeanPSTH{1,1};
MeanPSTH_FR_ExtraFP=MeanPSTH{1,2};

PSTH_Time_ExtraFP=nanmean(MeanPSTH_Time_ExtraFP,1);

%Overall mean PSTH
PSTH_FR_ExtraFP_mean=nanmean(MeanPSTH_FR_ExtraFP,1);
PSTH_FR_ExtraFP_sem=nanmean(MeanPSTH_FR_ExtraFP,1);
%Average FP
FR_Average_ExtraFP=MeanFRRaw{1,1};


UniqueAmplitude_ExtraFP=uniquetol(ExtraFP_Amp,0.001);
UniqueAngle_ExtraFP=uniquetol(ExtraFP_Angle(ExtraFP_Amp~=0),0.001);

for i=1:length(UniqueAmplitude_ExtraFP)
    currAmp=UniqueAmplitude_ExtraFP(i);
    sel_Amp=ExtraFP_Amp== currAmp;
    if sum(sel_Amp)>=3
        if currAmp==0
            %Don't need to loop over angle
            ExtraFP_FR_Amp_Sel=FR_Average_ExtraFP(sel_Amp);
             ExtraFP_FR_Amp_Mean(i,1:length(UniqueAngle_ExtraFP))=nanmean(ExtraFP_FR_Amp_Sel);
             ExtraFP_FR_Amp_Sem(i,1:length(UniqueAngle_ExtraFP))=nanstd(ExtraFP_FR_Amp_Sel)/sqrt(sum(~isnan(ExtraFP_FR_Amp_Sel)));
        else
            %Loop over the angle
            for j=1:length(UniqueAngle_ExtraFP)
                currAngle=UniqueAngle_ExtraFP(i);
                 
                    sel_Angle=sel_Amp & (ExtraFP_Angle== currAngle);
                    if sum(sel_Amp)>=3
                         ExtraFP_FR_Angle_Sel=FR_Average_ExtraFP(sel_Amp);
                
                         ExtraFP_FR_Amp_Mean(i,j)=nanmean(ExtraFP_FR_Angle_Sel);
                         ExtraFP_FR_Amp_Sem(i,j)=nanstd(ExtraFP_FR_Angle_Sel)/sqrt(sum(~isnan(ExtraFP_FR_Angle_Sel)));
                    else
                         ExtraFP_FR_Amp_Mean(i,j)=NaN;
                         ExtraFP_FR_Amp_Sem(i,j)=NaN;
                    end
            end
        end
        
        
        
        
    else
        ExtraFP_FR_Amp_Mean(i,1:length(UniqueAngle_ExtraFP))=NaN;
        ExtraFP_FR_Amp_Sem(i,1:length(UniqueAngle_ExtraFP))=NaN;
    end
 
    
end

%Center gaze responses 
ExtraFP_FR_Center_Sel=FR_Average_ExtraFP(ExtraFP_Amp== 0);
ExtraFP_FR_Center_Mean=nanmean(ExtraFP_FR_Center_Sel);
ExtraFP_FR_Center_Sem=nanstd(ExtraFP_FR_Center_Sel)/sqrt(sum(~isnan(ExtraFP_FR_Center_Sel)));

 %Group responses according to location vectors
 SeparationPoint=[22.5:45:360];
 
     
 SaccadeVectorTuning=PoloarHistGroup(ExtraFP_Angle(ExtraFP_Amp~=0),[FR_Average_ExtraFP(ExtraFP_Amp~=0)],SeparationPoint);
 SaccadeVector=SaccadeVectorTuning.Center;
 FiringRate_PreDurPost_mean=SaccadeVectorTuning.DataMean;
 FiringRate_PreDurPost_sem=SaccadeVectorTuning.DataSEM;
 SacVecTuning_p=SaccadeVectorTuning.ANOVA_p;
 
 PrefVector=SaccadeVectorTuning.PrefVector;
PrefVector_Amp=SaccadeVectorTuning.PrefVector_Amp;


if ShowFigureFlag
%Polar Plot for saccade vector tuning

figtitlestr{1}='LocationTuning_Average';

fig{1}=PrepareFigure(140,'w',[50,200, 500,500],'Name',figtitlestr{1});


subtitlestr={'ExtraFP_Location-Tuning'};

Selection_temp=~isnan(SaccadeVector);
SaccadeVector=SaccadeVector(Selection_temp);
FiringRate_PreDurPost_mean=FiringRate_PreDurPost_mean(Selection_temp,:);
FiringRate_PreDurPost_sem=FiringRate_PreDurPost_sem(Selection_temp,:);

if ~isempty(FiringRate_PreDurPost_mean)

%Pre-saccade , During, Post saccade tuning
%subplot(1,3,j);
%{
p= polarplot((SaccadeAngle_First)*pi/180',SpikeRate_PreSac,'or');
p.MarkerFaceColor='r';
p.MarkerEdgeColor='r';
%}
for j=1:1
%hold on
polarwitherrorbar([SaccadeVector/180*pi SaccadeVector(1)/180*pi],[FiringRate_PreDurPost_mean(:,j)' FiringRate_PreDurPost_mean(1,j)],[FiringRate_PreDurPost_sem(:,j)' FiringRate_PreDurPost_sem(1,j)],...
    '-b',3);
hold on
h_p = polarplot([PrefVector(j) PrefVector(j)]/180*pi,[0 PrefVector_Amp(j)],'r');
    set(h_p,'LineWidth',2.5);

set(gca,'ThetaDir' , 'counterclockwise');
set(gca,'ThetaZeroLocation','top')
set(gca,'FontSize',25,'FontWeight','Bold','LineWidth',3);

thetaticks(sort(SaccadeVector));
thetaticklabels(sort(SaccadeVector));
title(subtitlestr{j},'FontSize',25,'FontWeight','Bold');
pstr=sprintf('p= %2.3f',SacVecTuning_p(j));
a=annotation('textbox',[0.15+(j-1)/3.5,0.06,0.1,0.1] ,'String',pstr,'FitBoxToText','on','LineStyle','none');
a.FontSize=30;
a.FontWeight='Bold';
if SacVecTuning_p(j)<0.05
    a.Color='r';
end
hold on
%Add responses to the center Gaze
polarwitherrorbar([SaccadeVector/180*pi SaccadeVector(1)/180*pi],[repmat(ExtraFP_FR_Center_Mean,1,length(SaccadeVector)+1)],[repmat(ExtraFP_FR_Center_Sem,1,length(SaccadeVector)+1)],...
    '-k',3);

rlim([0,max(max(FiringRate_PreDurPost_mean))*1.2]);

%ExtraFP_FR_Center_Mean

 
end
else
    disp('No periperal extraFP');


end
end

end