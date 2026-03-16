function NewLearningAnalysis(Data, StartTrial, EndTrial,ShowFigureFlag,OutputFlag);
%Analysis for new learning
a=1;
ECode;
SegEvent=[SCENE_ON];%Event Flow of interest
%SegEventEnd=[SHOWFIXCD,FIXOFF,TGTOFF,RWDOFFCD ,NaN,EXTRAFP];%Event end
SegEventEnd=[SHOWFIXCD];%Event end
%SegEventEnd=[FIXOFF];%Event end
SegEventStr={'SceneOn'};


%Special request
remove_tmp=[];
SceneRemoval=[];
BlockRemoval=[];
%Remove the intermediate scene to focus on the pure scenes if needed
%SceneRemoval=[1063];


PlotStartOffset=[200];%in ms
PlotStopOffset=[700];%in ms

ResponseLatencyThreshold=3;

%Parameter to define the sliding window to measure firing rates
StepSize=5;%in ms
BinWidth=15;%in ms
%{
%Parameter seting for interval selecting plot period
StartOffset=100;%in ms
StopOffset=50;%in ms
%For scene offset
StartOffset2=500;%in ms
StopOffset2=50;%in ms
IntervalAfterMarkerStart=1000;%in ms;%This for align marker
IntervalAfterMarkerStartOff=500;%in ms;%For scene offset
%}
%Color Definition for each scenes

Color=[lines(6);prism(6);pink(6)];
%Color(2,:)=copper(1)+0.5;
Color=[Color;hsv(10)];
Color(18,:)=[0.1,0.457,0.194];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Unpack all the data base and task settings
ECode;
FileName=Data.FileName;
TaskType=Data.TaskType;
TaskCode=Data.TaskCode;
TrialNum=EndTrial-StartTrial+1;
Selection=StartTrial:EndTrial;

%Abnormal check %Especially For Robin081120FORAGEXPLO4.rst
CheckCode=SCENE_ON;%This code should appear once in each trial
AbnormalTrials=AbnormalTrialCheck(Data,CheckCode);

Selection=Selection(~ismember(Selection,AbnormalTrials));

if ~isempty(SceneRemoval)
    SceneName_tmp=Data.SceneIndex(Selection,3);
    remove_tmp=ismember(SceneName_tmp,SceneRemoval);
      
end
if ~isempty(remove_tmp)
    Selection=Selection(~remove_tmp);
end


%FixationPointIndex=Data.FixPointIndex(Selection,:);
%FixLoc=Data.FixLoc(Selection,:,:);


%TargetIndex=Data.TargetIndex(Selection,:,:);
%TargetLoc=Data.TargetLoc(Selection,:,:);

%TargetNum=Data.TargetNum(Selection,:);
%TargetSize=Data.TargetSize(Selection,:);
%TargetPolarLoc=Data.TargetPolarLoc(Selection,:,:);
%ChoiceOutcome=Data.ChoiceOutcome(Selection,:);


%GoodObjectLib=Data.GoodObjectLib;
%BadObjectLib=Data.BadObjectLib;
%SceneLib=Data.SceneLib;

%GoodObject=arrayfun(@(x)( ismember(x,GoodObjectLib(:,3))),TargetIndex(:,2:end,2));%Use real name
%GoodObjectID=TargetIndex(:,2:end,2);

%BadObject=arrayfun(@(x)( ismember(x,BadObjectLib(:,3))),TargetIndex(:,2:end,2));%Use real name
%BadObjectID=TargetIndex(:,2:end,2);
%GoodBadIndex=GoodObject*1+BadObject*2;%1:Good Object;2:Bad Object

SceneIndex=Data.SceneIndex(Selection,:);

%EyeChannel_X=Data.EyeChannel_X_L(Selection,:); %Only import the first two channels
%EyeChannel_Y=Data.EyeChannel_Y_L(Selection,:);
%EyeBinWidth=Data.EyeBinWidth;
%EyeTime=Data.EyeTime(Selection,:);

EventChannel=Data.EventChannel(Selection,:);
EventTimeChannel=Data.EventTimeChannel(Selection,:);
EventBin=Data.EventBin(Selection,:);

ChoiceOutcome=Data.ChoiceOutcome(Selection,2:end);
%GOODCHOICECD
%BADCHOICECD

%Change period automatically

StartCheck=FindOutTime(EventChannel,EventBin,SegEvent(1));
EndCheck=FindOutTime(EventChannel,EventBin,SegEventEnd(1));
EndCheck=EndCheck(:,1);
DiffCheck=EndCheck-StartCheck;
if sum(DiffCheck<100)>ceil(length(StartCheck)/2)
    SegEventEnd=[FIXOFF];%Event end;change from showfix cd to fix off for contigency
    disp('Change to FixOff for FPOn First Mode');
    
end



SpikeChannel=Data.SpikeChannel;
SpikeTime=squeeze(Data.SpikeTimeData(:,:,Selection))';
%SpikeBin=squeeze(Data.SpikeBinData(:,:,Selection))';





%Segregate spike data from scene onset to fixation onset

SegmentData=EventSegTool(SegEvent,SegEventEnd,PlotStartOffset,PlotStopOffset,EventChannel,EventTimeChannel,SpikeTime,StepSize,BinWidth);

%Retrive data from the field
%Total average for each event
MeanFR=SegmentData.MeanFR; %[Mean,sem]%Averall mean and sem
MeanFRRaw=SegmentData.MeanFRRaw{1};%=AverageFiringRate_Event;

%PSTH for each event
MeanPSTH=SegmentData.MeanPSTH;%=[PSTH_Event_Time',PSTH_Event'];
MeanPSTH_Time=MeanPSTH{1};
MeanPSTH_FR=MeanPSTH{2};


%Relative time bin for each event
RelativeTimeMarker=SegmentData.RelativeTimeMarker{1};%=EventTimeMarker;
%Absolute time marker for each event
AbosoluteTimeMarker=SegmentData.AbsoluteTimeMarker{1};

%{
%As there is only on scene onset for each trial, we don't need the
following 

%Average firing rate for each event for each incidance in each trial
MeanFREach=SegmentData.MeanFREach;%=AverageFiringRate_Event_Each;

%Average PSTH for each event for each incidance in each trial
MeanPSTHEach=SegmentData.MeanPSTHEach;%=[PSTH_Event_Each_Time',PSTH_Event_Each'];

%Event time marker for each event for each incidance in each trial
RelativeTimeMarkerEach=SegmentData.RelativeTimeMarkerEach;%=EventTimeMarkerEach;
AbosoluteTimeMarkerEach=SegmentData.AbosoluteTimeMarkerEach;%=AbosoluteTimeMarkerEach;
%}

%Get the average firing rate 200ms before scene onset as the baseline
SponEvent=[SCENE_ON];
PreTimeBuffer=200;
SponEventEnd=NaN;
StartOffset=0;
StopOffset=PreTimeBuffer;

EventTimeChannelForSpon=EventTimeChannel-PreTimeBuffer;


SegmentDataSpon=EventSegTool(SponEvent,SponEventEnd,StartOffset,StopOffset,EventChannel,EventTimeChannelForSpon,SpikeTime,StepSize,BinWidth);

Mean_Spon=SegmentDataSpon.MeanFR;
Mean_Spon_ForEachTrial=SegmentDataSpon.MeanFRRaw{1};

%Find out when the firing rate start to show scene response

[p_All,h_All]=ttest2_group(MeanPSTH_FR,Mean_Spon_ForEachTrial,'both');

p_All(1:ceil(PlotStartOffset/StepSize))=NaN;
h_All(1:ceil(PlotStartOffset/StepSize))=NaN;
h_All(h_All==0)=NaN;%For plot

%Loop over the scenes
%{
MeanPSTH_Time=MeanPSTH{1};
MeanPSTH_FR=MeanPSTH{2};
MeanFRRaw
%}
UniqueScene=unique(SceneIndex(:,3));
NumOfCondition=length(UniqueScene);
%Calculate and plot the average firing rate according to scenes

for i=1:NumOfCondition
    currentScene=UniqueScene(i);
    trialindex=SceneIndex(:,3)==currentScene;
    AverageFR_EachScene{i}=MeanFRRaw(trialindex);
    ChoiceOutcome_currScene=ChoiceOutcome(trialindex);
    
    BinaryChoice_currScene{i}=ChoiceOutcome_currScene==GOODCHOICECD;
end


%% Figures
FigureStartNum=120;
FigureIndex=1;

%If show figures 
if ShowFigureFlag
 %Figure 1: Firing rate change with trials
   figtitlestr{FigureIndex}='SceneResponseChangeWithTrial';
   fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[50,100, 2000,800],'Name',figtitlestr{FigureIndex});
   for i=1:NumOfCondition
    
    h=plot(AverageFR_EachScene{i},'-o','Color',Color(i,:),'LineWidth',3,'MarkerSize',20);
    xlabel('Trial Num');
    ylabel('FR(Hz)');
    if i<3
    set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
    set(h,'Marker','^','LineWidth',6)
    else
     set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
     set(h,'LineWidth',3)
     
    end
    
    hold on
    
    end
    
    legend(string(UniqueScene));
    
  FigureIndex=FigureIndex+1;
  FigureStartNum=FigureStartNum+1;   
  %Figure 2: Behavior change with trials
  figtitlestr{FigureIndex}='BehaviorChangeWithTrial';
   fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[50,100, 2000,800],'Name',figtitlestr{FigureIndex});
    for i=1:NumOfCondition
    
    h=plot(BinaryChoice_currScene{i}+i-1,'-','Color',Color(i,:),'LineWidth',3,'MarkerSize',20);
    xlabel('Trial Num');
    ylabel('Good Choice');
    if i<3
    set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
   % set(h,'Marker','^','LineWidth',6)
    else
     set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
   %  set(h,'LineWidth',3)
     
    end
    
    
    hold on
    
    end
    ylim([0,NumOfCondition+1]);
   FigureIndex=FigureIndex+1;
  FigureStartNum=FigureStartNum+1;   
   
   %BinaryChoice_currScene
    
    
end

end