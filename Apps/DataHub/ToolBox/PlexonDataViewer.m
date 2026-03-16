%Data Loader Temporary for laser-stim responses with V probe before virus
%injection
%09122022 by Xuefei Yu
clc;
ECode;
ParsedDataPath='/Users/yux8/WorkRelated_YXF/Program_Matlab_Local/DataAnalysis/PlexonData/ParsedData/';
Folder = '/020222/';
FileName = 'Robin020222N1T364.mat';


load(strcat(ParsedDataPath,Folder,FileName));%load su_task into workspace

SpikeTimeAll = su_task.SpikeTiming;
SpikeChannel = length(SpikeTimeAll);

Event = su_task.Event;
EventTime = su_task.EventTiming;

TrialTime=su_task.trialtimeinfo(:,2:end);
TrialEndTime=TrialTime(:,2);

TaskCode =unique(cell2mat(cellfun(@(x) x(x>=300 & x<=400), Event, 'UniformOutput', false)));




%Find out the good trials

GoodTrialsFinder=cellfun(@(x) sum(ismember(x,[CORRECTCD,ERROR])), Event)==1; %Now search for good trials among the whole dataset
GoodTrialsIndex=find(GoodTrialsFinder==1);
NumberOfGoodTrials=sum(GoodTrialsFinder);



EventChannel = cell2mat(Event(GoodTrialsIndex)')';
EventTimeChannel = cell2mat(EventTime(GoodTrialsIndex)')';


SpikeTimeGood =cellfun(@(x) x(GoodTrialsIndex),SpikeTimeAll,'uniform',0);

%Seperate into stim trial and sham trial
StimType=ReproduceFromEvent(EventChannel,[STIM1,SHAM0])';

StimTrial = StimType == STIM1;
ControlTrial = ~StimTrial;

%Laser Stim parameters--Manually constructed from Blip file
stimPower=ones(NumberOfGoodTrials,1);
stimPower(1:80) = 10;
stimPower(81:160) = 20;
stimPower(161:220) = 50;
stimPower(161:280) = 100;
stimPower(281:340) = 200;
stimPower(341:400) = 300;
stimPower(401:460) = 100;

stimDur = 200* ones(NumberOfGoodTrials,1);
stimDur(401:460) = 400;

%Chunk data according to stim-on
SegEvent = [ON_STIM];
SegEventEnd = [OFF_STIM];
SegEventStr={'StimOn'};

PlotStartOffset=[100];%in ms
PlotStopOffset=[100];%in ms

%Parameter to define the sliding window to measure firing rates
StepSize=5;%in ms
BinWidth=5;%in ms


    
%Transform spike data 
%Spike Channel
 
 SpikeTimeData=NaN*ones(max(size(SpikeTimeGood,2),1),max(1,max(max(cellfun(@numel,SpikeTimeGood))))+2,NumberOfGoodTrials);%Add one to store the spike time start and end   
  
 MaxPreLength=0;
 
 for i=1:NumberOfGoodTrials

    %Load the spiketimess
    
    for j=1:size(SpikeTimeData,1)
      
            SpikeTimeData(j,1:numel(SpikeTimeGood{j}{i}),i)=SpikeTimeGood{j}{i}';
 
   %Add spike time end 
            SpikeTimeData(j,end-1:end,i)=[0,TrialEndTime(i)];
           
  
   end
 end
 
 
for i = 1:SpikeChannel
    SpikeTime =squeeze(SpikeTimeData(i,:,:))';

SegmentData=EventSegTool(SegEvent,SegEventEnd,PlotStartOffset,PlotStopOffset,EventChannel,EventTimeChannel,SpikeTime,StepSize,BinWidth);

%Retrive data from the field
%Total average for each event
MeanFR=SegmentData.MeanFR; %[Mean,sem]%Averall mean and sem
MeanFRRaw=SegmentData.MeanFRRaw{1};%=AverageFiringRate_Event;

%Relative time bin for each event
RelativeTimeMarker=SegmentData.RelativeTimeMarker{1};%=EventTimeMarker;
%Absolute time marker for each event
AbosoluteTimeMarker=SegmentData.AbsoluteTimeMarker{1};


TimeMarker=SegmentData.RelativeTimeMarker{1}; %[Mean,sem]%Averall mean and sem

%For Raster
Raster=SegmentData.RasterAfterAlign{1};

SpikeCount = SegmentData.PSTH_Event_SpikeCount{1}/BinWidth*1000;


%PSTH for each event
MeanPSTH=SegmentData.MeanPSTH;%=[PSTH_Event_Time',PSTH_Event'];
MeanPSTH_Time=MeanPSTH{1};
MeanPSTH_FR=MeanPSTH{2};

MeanPSTH_Time_Mean = nanmean(MeanPSTH_Time);

PSTH_Stim = MeanPSTH_FR(StimTrial,:);
PSTH_Control = MeanPSTH_FR(ControlTrial,:);


PSTH_Stim_Mean = nanmean(PSTH_Stim,1);
PSTH_Stim_Sem = nanstd(PSTH_Stim,[],1)./sqrt(sum(~isnan(PSTH_Stim)));

PSTH_Control_Mean = nanmean(PSTH_Control,1);
PSTH_Control_Sem = nanstd(PSTH_Control,[],1)./sqrt(sum(~isnan(PSTH_Control)));

%Raster
Raster_Stim = Raster(StimTrial,:);
Raster_Control = Raster(ControlTrial,:);

%Spike Count

SpikeCount_Stim=SpikeCount(StimTrial,:);
SpikeCount_Control=SpikeCount(ControlTrial,:);

SpikeCountStim_Mean = nanmean(SpikeCount_Stim,1);
SpikeCountStim_Sem = nanstd(SpikeCount_Stim,[],1)./sqrt(sum(~isnan(SpikeCount_Stim)));

SpikeCountControl_Mean = nanmean(SpikeCount_Control,1);
SpikeCountControl_Sem = nanstd(SpikeCount_Control,[],1)./sqrt(sum(~isnan(SpikeCount_Control)));

%Plots 

    
FigureStartNum=90;
FigureIndex=1;
    
figtitlestr{FigureIndex}='PSTH_StimVSControl';
fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[50,100, 1000,600],'Name',figtitlestr{FigureIndex});

subplot(2,1,1)

AverageRegion = nanmean(RelativeTimeMarker);
area(AverageRegion,[max([SpikeCountControl_Mean+SpikeCountControl_Sem,SpikeCountStim_Mean+SpikeCountStim_Sem])+20,max([SpikeCountControl_Mean+SpikeCountControl_Sem,SpikeCountStim_Mean+SpikeCountStim_Sem])+20],...
    'FaceColor','k','FaceAlpha',.1,'EdgeAlpha',.1,'HandleVisibility','off');
hold on
%{
shadedErrorBar(MeanPSTH_Time_Mean,PSTH_Control_Mean,PSTH_Control_Sem,'lineprops',{[0,0,0]},'transparent',1,'patchSaturation',0.3)

shadedErrorBar(MeanPSTH_Time_Mean,PSTH_Stim_Mean,PSTH_Stim_Sem,'lineprops',{[1,0,0]},'transparent',1,'patchSaturation',0.3)
%}
shadedErrorBar(MeanPSTH_Time_Mean,SpikeCountControl_Mean,SpikeCountControl_Sem,'lineprops',{[0,0,0]},'transparent',1,'patchSaturation',0.3)

shadedErrorBar(MeanPSTH_Time_Mean,SpikeCountStim_Mean,SpikeCountStim_Sem,'lineprops',{[1,0,0]},'transparent',1,'patchSaturation',0.3)

xlim([-PlotStartOffset,PlotStopOffset+nanmean(RelativeTimeMarker(:,2))]);

%ylim([0,150]);
set(findall(gca, 'Type', 'Line'),'LineWidth',3);

ylabel('Firing Rate(Hz)')
set(gca,'LineWidth',3,'FontSize',20,'FontWeight','Bold');
title(sprintf("Laser Power : %d mW",unique(stimPower)));
legend(["Sham","Stim"])
box off


%{
%Figure2 SpikeRaterPlot
FigureStartNum = FigureStartNum+1;
FigureIndex = FigureIndex;

figtitlestr{FigureIndex}='Raster_StimVSControl';
fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[70,100, 1000,400],'Name',figtitlestr{FigureIndex});
%}
StartLoc = 1;
LineWidth = 2;
Barlength = 1;

subplot(2,1,2)
AverageRegion = nanmean(RelativeTimeMarker);
area(AverageRegion,[size(Raster,1)+5,size(Raster,1)+5],...
    'FaceColor','k','FaceAlpha',.1,'EdgeAlpha',.1,'HandleVisibility','off');

hold on

%Raster For Stim
SpikeRaterPlot(Raster_Stim,StartLoc,Barlength,'r',LineWidth);

StartLoc = size(Raster_Stim,1)+5;
%Raster For Control

SpikeRaterPlot(Raster_Control,StartLoc,Barlength,'k',LineWidth);
xlim([-PlotStartOffset,PlotStopOffset+nanmean(RelativeTimeMarker(:,2))]);
xlabel('Time from Opto-stim On');   

 set(gca,'LineWidth',3,'FontSize',20,'FontWeight','Bold');
 yticks([]);
 box off
 





end

keyboard
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

