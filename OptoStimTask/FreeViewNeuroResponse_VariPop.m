function FreeViewNeuroResponse_VariPop(Data, StartTrial, EndTrial,ShowFigureFlag,OutputFlag);
ECode;

TaskType=Data.TaskType;
DataPath= Data.Path;
FileName=Data.FileName;

Selection=StartTrial:EndTrial;%Select the trials according to the pannel


useReISI = 0;



SegEvent = [ON_STIM];
SegEventEnd = [OFF_STIM];
SegEventStr={'StimOn'};

PlotStartOffset=[100];%in ms
PlotStopOffset=[100];%in ms

%Sequence information from the library
ParaLib=Data.ParaLib;
if ismember('stimOn',keys(ParaLib))
   stimOnFlag = ParaLib('stimOn');
   TrialStimOn = unique(stimOnFlag(stimOnFlag(:,2) == 1,1));
   stimOn = [Selection',0*ones(length(Selection),1)];
   stimOn(ismember(stimOn(:,1),TrialStimOn),2) = 1;
  
   stimOn= stimOn(:,2);
%   stimOn=stimOn(Selection,:);
   
end
if ismember('stimPower',keys(ParaLib))
   stimPower_ori = ParaLib('stimPower');
   stimPower_ori=stimPower_ori(:,2);
   
  
   
end

if ismember('preStimDur,StimDur,PostStimDur',keys(ParaLib))
   stimdur_info = ParaLib('preStimDur,StimDur,PostStimDur');
   stimDur=stimdur_info(:,3);
   
   stimDur=stimDur(Selection,:);
   
end

%Manuel Modification of the power by mistake
if ismember(FileName,{'Robin061423VIDEOVIEW1.rst','Robin062223VIDEOVIEW1.rst','Robin071323VIDEOVIEW1.rst','Robin072023VIDEOVIEW1.rst','Robin072523VIDEOVIEW1.rst'})
    stimPower = RescuePower(FileName,stimPower_ori);
    
else
    stimPower = stimPower_ori;
end

 stimPower=stimPower(Selection,:);

%Parameter to define the sliding window to measure firing rates
StepSize=5;%in ms
BinWidth=5;%in ms

%{
%Event marker for segmentation
SegEvent=[PRE_STIM];%Event Flow of interest
SegEventEnd=NaN*ones(1,length(SegEvent));
%SegEventEnd(1) =SCENE_OFF;

SegEventStr={'FixOn','1TGT','2TGT','3TGT','4TGT','5TGT','LastFix'};
MarkerSizeEvent = 20*ones(1,length(SegEventStr));
MarkerSizeEvent([1,end]) = 10;

MarkerColorEvent = {'#44AA99','#FE6100','#FE6100','#FE6100','#FE6100','#FE6100','#44AA99'};
DistratorColorEvent = {'#44AA99','#648FFF','#648FFF','#648FFF','#648FFF','#648FFF','#44AA99'};

PlotStartOffset=500*ones(1,length(SegEvent));%in ms
PlotStopOffset=NaN*ones(1,length(SegEvent));%in ms
%}


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

%Load event data
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
    ChannelID=ChannelIDWhole(spk);
    ClusterID=ClusterIDWhole(spk);

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




%1/ISI
reISI = SegmentData.reISI{1};
reISIseq = SegmentData.reISIseq{1};


%PSTH for each event
MeanPSTH=SegmentData.MeanPSTH;%=[PSTH_Event_Time',PSTH_Event'];
MeanPSTH_Time=MeanPSTH{1};
MeanPSTH_FR=MeanPSTH{2};

MeanPSTH_Time_Mean = nanmean(MeanPSTH_Time);

%Seperate into stim trial and sham trial
StimType=ReproduceFromEvent(EventChannel,[STIM1,SHAM0])';

StimTrial = StimType == STIM1;
ControlTrial = ~StimTrial;

%PSTH_Stim = MeanPSTH_FR(StimTrial(stimOn==1),:);
%PSTH_Control = MeanPSTH_FR(ControlTrial(stimOn==1),:);

PSTH_Stim = MeanPSTH_FR(StimTrial,:);
PSTH_Control = MeanPSTH_FR(ControlTrial,:);

%ReISIseq
reISIseq_Stim = reISIseq(StimTrial,:);
reISIseq_Control = reISIseq(ControlTrial,:);



%Raster
Raster_Stim = Raster(StimTrial,:);
Raster_Control = Raster(ControlTrial,:);

%Spike Count

SpikeCount_Stim=SpikeCount(StimTrial,:);
SpikeCount_Control=SpikeCount(ControlTrial,:);



%Average

PSTH_Stim_Mean = nanmean(PSTH_Stim,1);
PSTH_Stim_Sem = nanstd(PSTH_Stim,[],1)./sqrt(sum(~isnan(PSTH_Stim)));

PSTH_Control_Mean = nanmean(PSTH_Control,1);
PSTH_Control_Sem = nanstd(PSTH_Control,[],1)./sqrt(sum(~isnan(PSTH_Control)));

SpikeCountStim_Mean = nanmean(SpikeCount_Stim,1);
SpikeCountStim_Sem = nanstd(SpikeCount_Stim,[],1)./sqrt(sum(~isnan(SpikeCount_Stim)));

SpikeCountControl_Mean = nanmean(SpikeCount_Control,1);
SpikeCountControl_Sem = nanstd(SpikeCount_Control,[],1)./sqrt(sum(~isnan(SpikeCount_Control)));
%}



%% Average according to the power and duration
UniquePower = uniquetol(stimPower,0.001);
UniqueDur = unique(stimDur);
FixDur = 100; %compare stimulation with different power with stimulating duration equal 100ms;
try
FixPower = uniquetol(stimPower(stimDur == 200),0.001);
catch
FixPower = uniquetol(stimPower(stimDur == 100),0.001);
end
%}
if isempty(FixPower)
    FixPower = uniquetol(stimPower(stimDur == 50),0.001);
    
end
if length(UniqueDur) ==1
    FixPower = max(UniquePower);
end
 ContinousBin = 3;
% First vary power:
maxFR_p = 0;
for i = 1:length(UniquePower)
    sel = stimPower == UniquePower(i) & stimDur == FixDur;
    if useReISI
        psth_stim_p = reISIseq(StimTrial & sel,:);
        psth_control_p = reISIseq(ControlTrial & sel,:);
    else
        psth_stim_p = SpikeCount(StimTrial & sel,:);
        psth_control_p = SpikeCount(ControlTrial & sel,:);
    end
     
    Mean_PSTH_stim_p{i,spk} = nanmean(psth_stim_p,1);
    Sem_PSTH_stim_p{i,spk} = nanstd(psth_stim_p,1)./sqrt(sum(~isnan(psth_stim_p)));
    
    Mean_PSTH_control_p{i,spk} = nanmean(psth_control_p,1);
    Sem_PSTH_control_p{i,spk} = nanstd(psth_control_p,1)./sqrt(sum(~isnan(psth_control_p)));
    
    max_p_tmp = max(max(Mean_PSTH_stim_p{i}+Sem_PSTH_stim_p{i}),max(Mean_PSTH_control_p{i}+Sem_PSTH_control_p{i}));
    if max_p_tmp > maxFR_p
        maxFR_p = max_p_tmp;
        
    end
    
   
    
    
end

%Second, vary duration
maxFR_d = 0;
for i = 1:length(UniqueDur)
    sel = stimPower == FixPower & stimDur == UniqueDur(i);
    
    
    if useReISI
        psth_stim_d = reISIseq(StimTrial & sel,:);
        psth_control_d = reISIseq(ControlTrial & sel,:);
    else
        psth_stim_d = SpikeCount(StimTrial & sel,:);
        psth_control_d = SpikeCount(ControlTrial & sel,:);
    end
    
  
     
    Mean_PSTH_stim_d{i} = nanmean(psth_stim_d,1);
    Sem_PSTH_stim_d{i} = nanstd(psth_stim_d,1)./sqrt(sum(~isnan(psth_stim_d)));
    
    Mean_PSTH_control_d{i} = nanmean(psth_control_d,1);
    Sem_PSTH_control_d{i} = nanstd(psth_control_d,1)./sqrt(sum(~isnan(psth_control_d)));
    
    max_d_tmp = max(max(Mean_PSTH_stim_d{i}+Sem_PSTH_stim_d{i}),max(Mean_PSTH_control_d{i}+Sem_PSTH_control_d{i}));
    if max_d_tmp > maxFR_d
        maxFR_d = max_d_tmp;
        
    end
    
  
    
end

end %End of the spike channel loop


NumPower = length(UniquePower);
NumDur = length(UniqueDur);


%%
%Plots 
if ShowFigureFlag
    
FigureStartNum=120;
FigureIndex=1;
    
figtitlestr{FigureIndex}='PSTH_StimVSControl_VariPower';
fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[50,100, 1200,600],'Name',figtitlestr{FigureIndex});


for i = 1: NumPower
    subplot(2,NumPower,i)
    
 %PSTH
    AverageRegion = nanmean(RelativeTimeMarker);
area(AverageRegion,[maxFR_p+20,maxFR_p+20],...
    'FaceColor','k','FaceAlpha',.1,'EdgeAlpha',.1,'HandleVisibility','off');
hold on
if useReISI == 1
    plot(MeanPSTH_Time_Mean,Mean_PSTH_control_p{i},'color',[0,0,0],'linewidth',3);
    plot(MeanPSTH_Time_Mean, Mean_PSTH_stim_p{i}, 'color',[1,0,0],'linewidth',3);
    
else
    
    shadedErrorBar(MeanPSTH_Time_Mean,Mean_PSTH_control_p{i},Sem_PSTH_control_p{i},'lineprops',{[0,0,0]},'transparent',1,'patchSaturation',0.3)
    shadedErrorBar(MeanPSTH_Time_Mean, Mean_PSTH_stim_p{i}, Sem_PSTH_stim_p{i},'lineprops',{[1,0,0]},'transparent',1,'patchSaturation',0.3)
end
xlim([-PlotStartOffset,PlotStopOffset+nanmean(RelativeTimeMarker(:,2))]);
ylim([0,maxFR_p+20]);

%ylim([0,150]);
set(findall(gca, 'Type', 'Line'),'LineWidth',3);


set(gca,'LineWidth',3,'FontSize',20,'FontWeight','Bold');

title(sprintf("Power : %1.2f mW",UniquePower(i)));
if i == 1
    ylabel('Firing Rate(Hz)');
    legend(["Sham","Stim"])
end
box off   

%Raster
subplot(2,NumPower,NumPower+i);
StartLoc = 1;
LineWidth = 2;
Barlength = 1;

AverageRegion = nanmean(RelativeTimeMarker);
area(AverageRegion,[size(Raster,1)+5,size(Raster,1)+5],...
    'FaceColor','k','FaceAlpha',.1,'EdgeAlpha',.1,'HandleVisibility','off');

hold on

%Raster For Stim
SpikeRaterPlot(Raster_Stim_p{i},StartLoc,Barlength,'r',LineWidth);

StartLoc = size(Raster_Stim_p{i},1)+5;
%Raster For Control

SpikeRaterPlot(Raster_Control_p{i},StartLoc,Barlength,'k',LineWidth);
xlim([-PlotStartOffset,PlotStopOffset+nanmean(RelativeTimeMarker(:,2))]);
if i == ceil(NumPower/2)
    xlabel('Time from Opto-stim On');   
end

 set(gca,'LineWidth',3,'FontSize',20,'FontWeight','Bold');
 yticks([]);
 box off


end

% Vary duration
FigureStartNum =FigureStartNum+1;
FigureIndex = FigureIndex+1;

figtitlestr{FigureIndex}='PSTH_StimVSControl_VariDur';
fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[70,100, 1200,600],'Name',figtitlestr{FigureIndex});


for i = 1: NumDur
    subplot(2,NumDur,i)
    
 %PSTH
    AverageRegion = [0,UniqueDur(i)];
area(AverageRegion,[maxFR_d+20,maxFR_d+20],...
    'FaceColor','k','FaceAlpha',.1,'EdgeAlpha',.1,'HandleVisibility','off');
hold on
if useReISI == 1
    plot(MeanPSTH_Time_Mean,Mean_PSTH_control_d{i},'color',[0,0,0],'linewidth',3);
    plot(MeanPSTH_Time_Mean, Mean_PSTH_stim_d{i}, 'color',[1,0,0],'linewidth',3);
    
else
    shadedErrorBar(MeanPSTH_Time_Mean,Mean_PSTH_control_d{i},Sem_PSTH_control_d{i},'lineprops',{[0,0,0]},'transparent',1,'patchSaturation',0.3)
    shadedErrorBar(MeanPSTH_Time_Mean, Mean_PSTH_stim_d{i}, Sem_PSTH_stim_d{i},'lineprops',{[1,0,0]},'transparent',1,'patchSaturation',0.3)
end
xlim([-PlotStartOffset,PlotStopOffset+UniqueDur(i)]);
ylim([0,maxFR_d+20]);

%ylim([0,150]);
set(findall(gca, 'Type', 'Line'),'LineWidth',3);


set(gca,'LineWidth',3,'FontSize',20,'FontWeight','Bold');

title(sprintf("Dur: %dms,(%1.2f mW)",UniqueDur(i),FixPower));
if i == 1
    if useReISI == 1
        ylabel('Spike Frequency Function(1/ISI)');
    else
        ylabel('Firing Rate(Hz)');
    end
    legend(["Sham","Stim"])
end
box off   

%Raster
subplot(2,NumDur,NumDur+i);
StartLoc = 1;
LineWidth = 2;
Barlength = 1;

AverageRegion = [0,UniqueDur(i)];
area(AverageRegion,[maxFR_d+20,maxFR_d+20],...
    'FaceColor','k','FaceAlpha',.1,'EdgeAlpha',.1,'HandleVisibility','off');
hold on

%Raster For Stim
SpikeRaterPlot(Raster_Stim_d{i},StartLoc,Barlength,'r',LineWidth);

StartLoc = size(Raster_Stim_d{i},1)+5;
%Raster For Control

SpikeRaterPlot(Raster_Control_d{i},StartLoc,Barlength,'k',LineWidth);
xlim([-PlotStartOffset,PlotStopOffset+UniqueDur(i)]);
if i == ceil(NumDur/2)
    xlabel('Time from Opto-stim On');   
end

 set(gca,'LineWidth',3,'FontSize',20,'FontWeight','Bold');
 yticks([]);
 box off


end





%{
%Figure2 SpikeRaterPlot
FigureStartNum = FigureStartNum+1;
FigureIndex = FigureIndex;

figtitlestr{FigureIndex}='Raster_StimVSControl';
fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[70,100, 1000,400],'Name',figtitlestr{FigureIndex});
%}

 



end %   End of if showFigureFlag

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