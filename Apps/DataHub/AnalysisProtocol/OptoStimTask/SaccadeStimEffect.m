function SaccadeStimEffect(Data, StartTrial, EndTrial,ShowFigureFlag,OutputFlag);
ECode;

TaskType=Data.TaskType;
DataPath= Data.Path;
FileName=Data.FileName;
TaskCode=Data.TaskCode;

Selection=StartTrial:EndTrial;%Select the trials according to the pannel

if FileName(1) == 'R'
    Monkey = 1;
else
    Monkey = 2;
end

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
%Definition of RF
if Monkey==1
    %Robin
    RF = UniqueAngle<180 & UniqueAngle>=0;
    NonRF = UniqueAngle>=180 & UniqueAngle<360;
else
    RF = UniqueAngle>=180 & UniqueAngle<360;
    NonRF = UniqueAngle>=0 & UniqueAngle<180;
end

   
     %Eye trace 
      
       

    for j = 1:length(UniqueAngle)
        sel = TargetAngle == UniqueAngle(j);
        sel_stim = sel & StimTrial;
        sel_sham = sel & ControlTrial;
        
        RT_stim_loc_period = RT_TargetToFixOff(sel_stim);
        RT_control_loc_period = RT_TargetToFixOff(sel_sham);

        RT_stim_loc_period_Mean(j) = nanmean(RT_stim_loc_period);
        RT_stim_loc_period_Sem(j) = nanstd(RT_stim_loc_period)/sqrt(sum(~isnan(RT_stim_loc_period)));

        RT_control_loc_period_Mean(j) = nanmean(RT_control_loc_period);
        RT_control_loc_period_Sem(j) = nanstd(RT_control_loc_period)/sqrt(sum(~isnan(RT_control_loc_period)));

        RF_Diff_Mean(j) = RT_stim_loc_period_Mean(j)-RT_control_loc_period_Mean(j);


        EyeX_Stim_loc{j} = EyeChannel_X_StimAlign(sel_stim,:);
        EyeY_Stim_loc{j} = EyeChannel_Y_StimAlign(sel_stim,:);
       
        EyeX_Control_loc{j} = EyeChannel_X_StimAlign(sel_sham,:);
        EyeY_Control_loc{j} = EyeChannel_Y_StimAlign(sel_sham,:);

        TargetAngle_curr = UniqueAngle(j);
        RF_Angle =UniqueAngle(RF);
        NonRF_Angle = UniqueAngle(NonRF);
        

        if TargetAngle_curr == RF_Angle
            [h,p_loc(j)] = ttest2(RT_stim_loc_period,RT_control_loc_period,'tail','left'); %one-tailed t-test
            [p_loc_sr(j),h]=ranksum(RT_stim_loc_period,RT_control_loc_period,'tail','left');
        elseif TargetAngle_curr == NonRF_Angle

            [h,p_loc(j)] = ttest2(RT_stim_loc_period,RT_control_loc_period,'tail','right');
            [p_loc_sr(j),h] = ranksum(RT_stim_loc_period,RT_control_loc_period,'Tail','right');
        end

        
        
  
    end


%

   EyeX_Stim_RF = EyeX_Stim_loc{RF};
    EyeY_Stim_RF = EyeY_Stim_loc{RF};
    
   EyeEcc_Stim_RF = (EyeX_Stim_RF.^2+EyeY_Stim_RF.^2).^(1/2);
    

    EyeX_Control_RF = EyeX_Control_loc{RF};
    EyeY_Control_RF = EyeY_Control_loc{RF};

   EyeEcc_Control_RF = (EyeX_Control_RF.^2+EyeY_Control_RF.^2).^(1/2);

   
    RT_stim_RF_period_Mean = RT_stim_loc_period_Mean(RF);
    RT_control_RF_period_Mean = RT_control_loc_period_Mean(RF);

    RT_stim_RF_period_Sem = RT_stim_loc_period_Sem(RF);
    RT_control_RF_period_Sem = RT_control_loc_period_Sem(RF);


    p_RF = p_loc(RF);
    p_RF_sr = p_loc_sr(RF);

    
    EyeX_Stim_NonRF = EyeX_Stim_loc{NonRF};
    EyeY_Stim_NonRF = EyeY_Stim_loc{NonRF};

    EyeEcc_Stim_NonRF = (EyeX_Stim_NonRF.^2+EyeY_Stim_NonRF.^2).^(1/2);

    EyeX_Control_NonRF = EyeX_Control_loc{NonRF};
    EyeY_Control_NonRF = EyeY_Control_loc{NonRF};

    EyeEcc_Control_NonRF = (EyeX_Control_NonRF.^2+EyeY_Control_NonRF.^2).^(1/2);


    RT_stim_NonRF_period_Mean = RT_stim_loc_period_Mean(NonRF);
    RT_control_NonRF_period_Mean = RT_control_loc_period_Mean(NonRF);

    RT_stim_NonRF_period_Sem = RT_stim_loc_period_Sem(NonRF);
    RT_control_NonRF_period_Sem = RT_control_loc_period_Sem(NonRF);

    p_NonRF = p_loc(NonRF);
    p_NonRF_sr = p_loc_sr(NonRF);

%{
if Monkey==1
    %Robin
   % RF = UniqueAngle<180 & UniqueAngle>=0;
    EyeX_Stim_RF = EyeX_Stim_loc{RF};
    EyeY_Stim_RF = EyeY_Stim_loc{RF};

    EyeEcc_Stim_RF = (EyeX_Stim_RF.^2+EyeY_Stim_RF.^2).^(1/2);
    

    EyeX_Control_RF = EyeX_Control_loc{RF};
    EyeY_Control_RF = EyeY_Control_loc{RF};

   EyeEcc_Control_RF = (EyeX_Control_RF.^2+EyeY_Control_RF.^2).^(1/2);

   
    RT_stim_RF_period_Mean = RT_stim_loc_period_Mean(RF);
    RT_control_RF_period_Mean = RT_control_loc_period_Mean(RF);

    RT_stim_RF_period_Sem = RT_stim_loc_period_Sem(RF);
    RT_control_RF_period_Sem = RT_control_loc_period_Sem(RF);


    p_RF = p_loc(RF);

    NonRF = UniqueAngle>=180 & UniqueAngle<360;
    EyeX_Stim_NonRF = EyeX_Stim_loc{NonRF};
    EyeY_Stim_NonRF = EyeY_Stim_loc{NonRF};

    EyeEcc_Stim_NonRF = (EyeX_Stim_NonRF.^2+EyeY_Stim_NonRF.^2).^(1/2);

    EyeX_Control_NonRF = EyeX_Control_loc{NonRF};
    EyeY_Control_NonRF = EyeY_Control_loc{NonRF};

    EyeEcc_Control_NonRF = (EyeX_Control_NonRF.^2+EyeY_Control_NonRF.^2).^(1/2);


    RT_stim_NonRF_period_Mean = RT_stim_loc_period_Mean(NonRF);
    RT_control_NonRF_period_Mean = RT_control_loc_period_Mean(NonRF);

    RT_stim_NonRF_period_Sem = RT_stim_loc_period_Sem(NonRF);
    RT_control_NonRF_period_Sem = RT_control_loc_period_Sem(NonRF);

    p_NonRF = p_loc(NonRF);

else
    %Adams
   % RF = UniqueAngle>=180 & UniqueAngle<360;
    EyeX_Stim_RF = EyeX_Stim_loc{RF};
    EyeY_Stim_RF = EyeY_Stim_loc{RF};

    EyeEcc_Stim_RF = (EyeX_Stim_RF.^2+EyeY_Stim_RF.^2).^(1/2);

    EyeX_Control_RF = EyeX_Control_loc{RF};
    EyeY_Control_RF = EyeY_Control_loc{RF};

   EyeEcc_Control_RF = (EyeX_Control_RF.^2+EyeY_Control_RF.^2).^(1/2);



    RT_stim_RF_period_Mean = RT_stim_loc_period_Mean(RF);
    RT_control_RF_period_Mean = RT_control_loc_period_Mean(RF);

    RT_stim_RF_period_Sem = RT_stim_loc_period_Sem(RF);
    RT_control_RF_period_Sem = RT_control_loc_period_Sem(RF);

    p_RF = p_loc(RF);

    NonRF = UniqueAngle>=0 & UniqueAngle<180;
    EyeX_Stim_NonRF = EyeX_Stim_loc{NonRF};
    EyeY_Stim_NonRF = EyeY_Stim_loc{NonRF};

    EyeEcc_Stim_NonRF = (EyeX_Stim_NonRF.^2+EyeY_Stim_NonRF.^2).^(1/2);

     EyeX_Control_NonRF = EyeX_Control_loc{NonRF};
    EyeY_Control_NonRF = EyeY_Control_loc{NonRF};

    EyeEcc_Control_NonRF = (EyeX_Control_NonRF.^2+EyeY_Control_NonRF.^2).^(1/2);


    RT_stim_NonRF_period_Mean = RT_stim_loc_period_Mean(NonRF);
    RT_control_NonRF_period_Mean = RT_control_loc_period_Mean(NonRF);

    RT_stim_NonRF_period_Sem = RT_stim_loc_period_Sem(NonRF);
    RT_control_NonRF_period_Sem = RT_control_loc_period_Sem(NonRF);


    p_NonRF = p_loc(NonRF);


end
%}
disp('RT_stim_RF(Mean):')
disp(RT_stim_RF_period_Mean)
disp('RF_control_RF(Mean):')
disp(RT_control_RF_period_Mean)
disp('RT_stim_RF(Sem):')
disp(RT_stim_RF_period_Sem)
disp('RF_control_RF(Sem):')
disp(RT_control_RF_period_Sem)



disp('RF_stim_nonRF(Mean):')
disp(RT_stim_NonRF_period_Mean)
disp('RF_control_nonRF(Mean):')
disp(RT_control_NonRF_period_Mean)
disp('RT_stim_nonRF(Sem):')
disp(RT_stim_NonRF_period_Sem)
disp('RF_control_nonRF(Sem):')
disp(RT_control_NonRF_period_Sem)



%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% For figure plotting%%%%%%%%%%%%%%%%%

if ShowFigureFlag
   StartNO = 50;
   FigureIndex=1;
   figtitlestr{FigureIndex}='RTDifference';
   fig{FigureIndex}=PrepareFigure(StartNO+FigureIndex-1,'w',[50,100, 1200,800],'Name',figtitlestr{FigureIndex});
 
   b=bar([1,2,3,4],[RT_stim_RF_period_Mean,RT_control_RF_period_Mean,RT_stim_NonRF_period_Mean,RT_control_NonRF_period_Mean]);
   b.FaceColor = 'flat';
   b.CData=[1,0,0;0,0,0;0.3,0.1,0;0.6,0.6,0.6];
   hold on
   errorbar([1,2,3,4],[RT_stim_RF_period_Mean,RT_control_RF_period_Mean,RT_stim_NonRF_period_Mean,RT_control_NonRF_period_Mean],[RT_stim_RF_period_Sem,RT_control_RF_period_Sem,RT_stim_NonRF_period_Sem,RT_control_NonRF_period_Sem],'.k');
   xticklabels({'RF-Stim','RF-Control','NonRF-Stim','NonRF-Control'});

   legend({sprintf('pRF=%1.2f',p_RF),sprintf('pNonRF=%1.2f',p_NonRF)});
   box off


StartNO = StartNO+1;
   FigureIndex=FigureIndex+1;

   figtitlestr{FigureIndex}='RawEyeTraceDifference';
   fig{FigureIndex}=PrepareFigure(StartNO+FigureIndex-1,'w',[50,100, 1200,800],'Name',figtitlestr{FigureIndex});
   subplot(1,2,1);
   %Towards RF 
   plot(repmat(EyeTime_Stim_Align,size(EyeEcc_Control_RF,1),1)',EyeEcc_Control_RF','-k','LineWidth',1);
   hold on
   plot(repmat(EyeTime_Stim_Align,size(EyeEcc_Stim_RF,1),1)',EyeEcc_Stim_RF','-r','LineWidth',1);

   subplot(1,2,2);
   %Towards NonRF
   plot(repmat(EyeTime_Stim_Align,size(EyeEcc_Control_NonRF,1),1)',EyeEcc_Control_NonRF','-k','LineWidth',1);
   hold on
   plot(repmat(EyeTime_Stim_Align,size(EyeEcc_Stim_NonRF,1),1)',EyeEcc_Stim_NonRF','-r','LineWidth',1);

 

        
   
    
   
end %End of the if show figure


%% Output to files
%Setup output directory
Workingdirectory=pwd;
MarkerFolder='DataAnalysis';
%MarkerFolder='DataHub';
Flag=strfind(Workingdirectory,MarkerFolder);
BasicDirectory=Workingdirectory(1:Flag+length(MarkerFolder));


%%%%%%%%%%%Export data%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Set up output path
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
    if isfield(OutputData,'SaccadeStim')
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
%OutputData For current analysis
OutputData_New=[];

OutputData_New.SaccadeStim.TrialType='GoodTrials';

%%%%%%%%%Save overall scene response characteristics%%%%%%%%%%%%%%%%%%%%%%%
%Organize the data


NameStr={'RT_stim_RF_period_Mean','RT_control_RF_period_Mean','RT_stim_NonRF_period_Mean','RT_control_NonRF_period_Mean','p_RF','p_NonRF',...
    'p_RF_sr','p_NonRF_sr'
    };

DataLib={RT_stim_RF_period_Mean,RT_control_RF_period_Mean,RT_stim_NonRF_period_Mean,RT_control_NonRF_period_Mean,p_RF,p_NonRF,...
    p_RF_sr,p_NonRF_sr

};



DataStamp=containers.Map(NameStr,DataLib);





OutputData_New.SaccadeStim.Task=TaskType;
OutputData_New.SaccadeStim.TaskCode=TaskCode;
OutputData_New.SaccadeStim.DataStamp=DataStamp;

OutputData_New.SaccadeStim.StartTrial =StartTrial;
%{
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
    
    for jj=1:TotalFigure
        if ~isempty(fig{jj})
        
         OutputFigureName_Final=strcat(OutputFigureName,figtitlestr{jj},'.jpg');
         saveas(fig{jj},OutputFigureName_Final);
        end
        
            
        
    end

    disp('Figures have been exported to the neuron folder');

end

%}





if ExistFlag
%Compare the old one with the new one
Task_Old=OutputData.SaccadeStim.TaskCode;

     if Task_Old~=TaskCode
        %Not the same task, add the new dataset into the old one
        OutputData(numel(OutputData)+1)=OutputData_New;
          
     else
         %If the same task, replace the old one with the new one
         %OutputData=OutputData_New;
         OutputData.SaccadeStim=OutputData_New.SaccadeStim;
         
    
    
     end
else
    %Output the current file
   OutputData.SaccadeStim=OutputData_New.SaccadeStim; 
end



if OutputFlag
    cd(OutPath);
    save(OutputFileName,'OutputData');
    disp('Data Exported');
else
    disp('No data export');
end  





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