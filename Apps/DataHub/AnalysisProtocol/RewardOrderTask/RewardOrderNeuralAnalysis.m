function RewardOrderNeuralAnalysis(Data, StartTrial, EndTrial,ShowFigureFlag,OutputFlag);
ECode;
TaskType=Data.TaskType;
DataPath= Data.Path;
FileName=Data.FileName;

Selection=StartTrial:EndTrial;%Select the trials according to the pannel

%Criterias for sexwlecting saccades
V_Threshold=40;
ContinueBin=5;
ISI_Threshold=20;

PlotTaskInfo = 1;

SegEvent=[FIRSTTGT,SECONDTGT];%Event Flow of interest
%SegEventEnd=[SHOWFIXCD,FIXOFF,TGTOFF,RWDOFFCD ,NaN,EXTRAFP];%Event end
SegEventEnd=[NaN,NaN];%Event end
SegEventStr={'FirstObjOn','SndObjOn'};

PlotStartOffset=[500,500];%in ms
PlotStopOffset=[NaN,NaN];%in ms

%Parameter to define the sliding window to measure firing rates
StepSize=5;%in ms
BinWidth=20;%in ms

Color={'#648FFF','#FE6100','#009E73','#88CCEE','#CC6677','#785EF0'};
Color_RGB=[100,143,255;
           254,97,0;
           0,158,115;
          136,204,238;
          153,153,0;
          120,94,240]/255;



%{
%Abnormal check %Especially For Robin081120FORAGEXPLO4.rst
CheckCode=SCENE_ON;%This code should appear once in each trial
AbnormalTrials=AbnormalTrialCheck(Data,CheckCode);
Selection=Selection(~ismember(Selection,AbnormalTrials));
TrialNum=length(Selection);
%}
%Extra information : Task Type,  Reward Order
ParaLib=Data.ParaLib;
if(ismember('TaskType',keys(ParaLib)) && (ismember('CurrentRewardOrder',keys(ParaLib))))
    TaskType_Info=ParaLib('TaskType');
    TaskType=TaskType_Info(Selection,end);
    
    RewardOrderInfo = ParaLib('CurrentRewardOrder');
    RewardOrder = RewardOrderInfo(Selection,end-1:end);
    
    BlockModeInfo = ParaLib('CurrentBlockMode');
    BlockMode = BlockModeInfo(Selection, end);
    
    
else
    disp("Not enough saved paramerer!");
    return 
    
end

%Select all choice trials
ChoiceTrials = TaskType == 2;
%Selection = Selection(ChoiceTrials);


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


%Load parameters
try
    TargetName=Data.TargetIndex(Selection,2:end,2);
    TargetIndex=Data.TargetIndex(Selection,2:end,1);
catch
    TargetName=Data.TargetIndex(Selection,3);
    TargetIndex=Data.TargetIndex(Selection,2);
    
end
TargetLoc=Data.TargetLoc(Selection,2:end,:);

%TargetNum=Data.TargetNum(Selection,2:end);
TargetNum=sum(~isnan(TargetName),2);
%TargetSize=Data.TargetSize(Selection,2:end);
try
    TargetPolarAngle=Data.TargetPolarLoc(Selection,2:end,1);

    TargetPolarEccentricity=Data.TargetPolarLoc(Selection,2:end,2);
catch
    TargetPolarAngle=Data.TargetPolarLoc(Selection,2);

    TargetPolarEccentricity=Data.TargetPolarLoc(Selection,3);
    
end

RewardOrder = RewardOrder(Selection,:);
BlockMode =BlockMode(Selection,:);
TaskType =TaskType(Selection,:);

ChoiceTrial = TaskType == 2;

%% %%%Analysis the choices %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ObjectChosen_tmp=ReproduceFromEvent(EventChannel,[GGFRAC : DBFRAC])';
ObjectChosen = cellfun(@(x) x(1),ObjectChosen_tmp);

NumOfCondition = 6;
ChoicesMat = NaN*ones(NumOfCondition,NumOfCondition);
NumTrialsCon =  ChoicesMat;

RewardOrderChoice = RewardOrder(ChoiceTrials,:);
ObjectChosenChoice = ObjectChosen(ChoiceTrials,:);

%% Load the event channel
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

    %%
%%%%%%%%%%%%%%%%%%%%%%%Align data on object onset
%Segregate spike data from fractal onset to fractal offset

SegmentData=EventSegTool(SegEvent,SegEventEnd,PlotStartOffset,PlotStopOffset,EventChannel,EventTimeChannel,SpikeTime,StepSize,BinWidth);

%Retrive data from the field
%{
%Relative time bin for each event
RelativeTimeMarker=SegmentData.RelativeTimeMarker{1};%=EventTimeMarker;
%Absolute time marker for each event
AbosoluteTimeMarker=SegmentData.AbsoluteTimeMarker{1};

ObjOnTimeMarker=nanmean(RelativeTimeMarker);
%}


%PSTH for each event
MeanPSTH=SegmentData.MeanPSTH;%=[PSTH_Event_Time',PSTH_Event'];


MeanPSTH_FirstTgt_Time=MeanPSTH{1,1};
MeanPSTH_FR_FirstFractalOn=MeanPSTH{1,2};

MeanPSTH_SecondTgt_Time=MeanPSTH{2,1};
MeanPSTH_FR_SecondFractalOn=MeanPSTH{2,2};

PSTH_Time_First=nanmean(MeanPSTH_FirstTgt_Time,1);
PSTH_Time_Second=nanmean(MeanPSTH_SecondTgt_Time,1);


%Loop over the direction and condition
UniqueAngle = unique(TargetPolarAngle(:,1));

for i = 1:length(UniqueAngle)
   idx = 1;
    for j = GGFRAC:DBFRAC
        
        sel = TargetPolarAngle(:,1) == UniqueAngle(i) & RewardOrder(:,1) == j;

        PSTH_First = MeanPSTH_FR_FirstFractalOn(sel & ~ChoiceTrial,:);
        PSTH_First_Mean{i,idx} = nanmean(PSTH_First,1);
        PSTH_First_Sem{i,idx}= nanstd(PSTH_First,[],1)./sqrt(sum(~isnan(PSTH_First),1));

        PSTH_Second = MeanPSTH_FR_SecondFractalOn(sel & ~ChoiceTrial,:);
        PSTH_Second_Mean{i,idx} = nanmean(PSTH_Second,1);
        PSTH_Second_Sem{i,idx}= nanstd(PSTH_Second,[],1)./sqrt(sum(~isnan(PSTH_Second),1));

        idx=idx+1;



    end
end

%%
if ShowFigureFlag
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Figure 1:PSTH aligned for
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%first saccades


figtitlestr{1}='PSTHAlignOnFracOnset_EachLocation';

fig{1}=PrepareFigure(110,'w',[50,100, 1200,800],'Name',figtitlestr{1});

TimeInv = [-400,500];
TimeInv2 = [-400,500];
time_sel = PSTH_Time_First >= TimeInv(1) & PSTH_Time_First <= TimeInv(2) ;

time_sel2 = PSTH_Time_Second >= TimeInv2(1) & PSTH_Time_Second <= TimeInv2(2) ;
for i = 1:length(UniqueAngle)
    
    
    idx = 1;
    for j = GGFRAC:DBFRAC
        p1((i-1)*length(UniqueAngle)+1) = subplot(length(UniqueAngle),2,(i-1)*2+1);

        shadedErrorBar(PSTH_Time_First(time_sel),PSTH_First_Mean{i,idx}(time_sel),PSTH_First_Sem{i,idx}(time_sel),'lineprops',{Color_RGB(idx,:)},'transparent',1,'patchSaturation',0.3);
        hold on
   

     xlim(TimeInv );%Keep in the same scale
       set(findall(gca, 'Type', 'Line'),'LineWidth',3);

        set(gca,'FontSize',15,'FontWeight','Bold');
        set(gca,'LineWidth',3);
        title(sprintf("FirstFracOn(Angle: %1.1d)",UniqueAngle(i)));

        p2((i-1)*length(UniqueAngle)+1)=subplot(length(UniqueAngle),2,(i-1)*2+2);

        shadedErrorBar(PSTH_Time_Second(time_sel2),PSTH_Second_Mean{i,idx}(time_sel2),PSTH_Second_Sem{i,idx}(time_sel2),'lineprops',{Color_RGB(idx,:)},'transparent',1,'patchSaturation',0.3);
        hold on
    idx=idx+1;

     xlim(TimeInv );%Keep in the same scale
       set(findall(gca, 'Type', 'Line'),'LineWidth',3);

        set(gca,'FontSize',15,'FontWeight','Bold');
        set(gca,'LineWidth',3);
        title(sprintf("SecondFracOn(Angle%1.1d)",UniqueAngle(i)));


        
         
        
        

    end
    linkaxes([p1,p2],'y');
    subplot(length(UniqueAngle),2,1);
    
            
        xlabel('Time from Firstractal Onset(ms)','FontSize',15,'FontWeight','Bold');

        legend({'GG','GB','BG','BB','DG','DB'});
        
        ylabel('Firing Rate(Hz)','FontSize',15,'FontWeight','Bold');

        subplot(length(UniqueAngle),2,2);
     xlabel('Time from SecondFractal Onset(ms)','FontSize',15,'FontWeight','Bold');
            
       

        legend({'GG','GB','BG','BB','DG','DB'});
        
        ylabel('Firing Rate(Hz)','FontSize',15,'FontWeight','Bold');

end

end %End of if showfigureflag


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
ProtocolMarker='RewardOrder';

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


end %End of the spike channel loop


end % End of the function