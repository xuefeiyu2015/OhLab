function SceneTuningFreeView(Data, StartTrial, EndTrial,ShowFigureFlag,OutputFlag)
%%%%%%%%%%%%%%%%%User controller starts here
ECode


%Event of interest
%{
SegEvent=[SCENE_ON,SHOWFIXCD,TGTCD,RWDCD ,SCENE_OFF,SCENE_OFF];%Event Flow of interest
%SegEventEnd=[SHOWFIXCD,FIXOFF,TGTOFF,RWDOFFCD ,NaN,EXTRAFP];%Event end
SegEventEnd=[SHOWFIXCD,FIXOFF,TGTOFF,RWDOFFCD ,NaN,NaN];%Event end
SegEventStr={'SceneOn','FixOn','TgtOn','RewardOn','SceneOff','Spon'};

PlotStartOffset=[100,100,100,100,100,-1000];%in ms
PlotStopOffset=[100,100,100,100,200,2000];%in ms
%}

SegEvent=[SCENE_ON];%Event Flow of interest
%SegEventEnd=[SHOWFIXCD,FIXOFF,TGTOFF,RWDOFFCD ,NaN,EXTRAFP];%Event end
%SegEventEnd=[SHOWFIXCD];%Event end
SegEventEnd=[SCENE_OFF ];%Event end
SegEventStr={'SceneOn'};



PlotStartOffset=[1000];%in ms
PlotStopOffset=[1000];%in ms

ResponseLatencyThreshold=3;

%Parameter to define the sliding window to measure firing rates
StepSize=10;%in ms
BinWidth=100;%in ms
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
Color=[Color;Color;Color];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Unpack all the data base and task settings
ECode;
FileName=Data.FileName;
TaskType=Data.TaskType;
TaskCode=Data.TaskCode;
TrialNum=EndTrial-StartTrial+1;
Selection=StartTrial:EndTrial;

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
%EventBin=Data.EventBin(Selection,:);

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

[p_All,h_All]=ttest2_group(MeanPSTH_FR,Mean_Spon_ForEachTrial,'right');

p_All(1:ceil(PlotStartOffset/StepSize))=NaN;
h_All(1:ceil(PlotStartOffset/StepSize))=NaN;
h_All(h_All==0)=NaN;%For plot

%Find out response latency

ResponseLatency=FindOutLatency(nanmean(MeanPSTH_Time),h_All,ResponseLatencyThreshold);

%Calculate response duration--the total time within particular period
%during the whole response window
ResponseDuration=nanmean(MeanPSTH_Time);
ResponseDuration=ResponseDuration(h_All==1);
ResponseDuration=ResponseDuration(ResponseDuration<=min(RelativeTimeMarker(:,2)) & ResponseDuration>=ResponseLatency);
TotalResponseDuration=length(ResponseDuration)*StepSize;
TotalEventDuration=min(RelativeTimeMarker(:,2));
PropResponseDuration=TotalResponseDuration/TotalEventDuration*100;

%PSTH mean and sem

MeanPSTH_FR_mean=nanmean(MeanPSTH_FR);
MeanPSTH_FR_sem=nanstd(MeanPSTH_FR)./sqrt(sum(~isnan(MeanPSTH_FR),1));


%First significant interval 
 %Peak Firing Rate in the first significant interval
Significant_TimePoint_temp=ResponseDuration(ResponseDuration>=ResponseLatency);
FirstIntervalEnd_Whole=find(diff(Significant_TimePoint_temp)>StepSize,1,'first');

    if isempty(FirstIntervalEnd_Whole) || FirstIntervalEnd_Whole>length(Significant_TimePoint_temp)
        FirstIntervalEnd_Whole=length(Significant_TimePoint_temp);
        
    end
    
    
FirstSigInterval_Whole=[1:FirstIntervalEnd_Whole];


if isempty(ResponseLatency) || isempty(FirstSigInterval_Whole)
    FR_Sig_Whole=NaN;
    PeakFiringTime=NaN;
    FirstIntervalDuration=NaN;
    ResponseLatency=NaN;
    PeakFiring=NaN;
else
    

FR_Sig_Whole=MeanPSTH_FR_mean(nanmean(MeanPSTH_Time)>=ResponseLatency);
    
[PeakFiring,loc]=max(FR_Sig_Whole(FirstSigInterval_Whole));
    
PeakFiringTime=Significant_TimePoint_temp(loc);
    
%First Interval Duration
FirstIntervalDuration=diff(Significant_TimePoint_temp([FirstSigInterval_Whole(1),FirstSigInterval_Whole(end)]));
end 


%Find PSTH for each scene

%#################Separate trials according to scene name
UniqueScene=unique(SceneIndex(:,3));
NumOfCondition=length(UniqueScene);
%Calculate and plot the average firing rate according to scenes

for i=1:NumOfCondition
    currentScene=UniqueScene(i);
    trialindex=SceneIndex(:,3)==currentScene;
    %For average Scene Tuning
    FR_Scene=MeanFRRaw(trialindex');
  
    FR_Scene_mean(i)=mean(FR_Scene);
    FR_Scene_sem(i)=std(FR_Scene)/sqrt(sum(trialindex));
    FR_Scene_SSE(i)=sum((FR_Scene-mean(FR_Scene)).^2); 
    
    %Measure the mean firing rate of the spon for each scene
    Spon_ForEachScene=Mean_Spon_ForEachTrial(trialindex);
    Spon_ForEachScene_Mean(i)=mean(Mean_Spon_ForEachTrial(trialindex));
    Spon_ForEachScene_SEM(i)=std(Mean_Spon_ForEachTrial(trialindex))/sqrt(sum(trialindex));
    
     FR_Spon_SSE(i)=sum((Spon_ForEachScene-mean(Spon_ForEachScene)).^2); 
    
    %p value from the spon
    [h_SceneFromSpon(i),p_SceneFromSpon(i)]=ttest2(currentScene,Spon_ForEachScene,'tail','right');
    
    
    
    
    
   %For scene-on  PSTH
    SpikeFiringRateSceneOn_mean(i,:)=mean(MeanPSTH_FR(trialindex,:),1);
    SpikeFiringRateSceneOn_sem(i,:)=std(MeanPSTH_FR(trialindex,:),[],1)/sqrt(sum(trialindex));
    try
    RelativeTimeMarkerForeachScene(i)=min(RelativeTimeMarker(trialindex,2));
    catch
        keyboard
    end
    
    SceneStr{i}=sprintf('Scene%d',UniqueScene(i));
    
    %SSE dynamics
    FR_SSE_Scene_Time(i,:)=sum((MeanPSTH_FR(trialindex,:)-mean(MeanPSTH_FR(trialindex,:),1)).^2,1); 
   
    
    
    %Find out response latency, response duration, peak firing rate, peak
    %firing time;
    
    [p_EachScene(i,:),h_EachScene(i,:)]=ttest2_group(MeanPSTH_FR(trialindex,:),Mean_Spon_ForEachTrial,'right');
    
    p_EachScene(i,1:ceil(PlotStartOffset/StepSize))=NaN;
    h_EachScene(i,1:ceil(PlotStartOffset/StepSize))=NaN;
    h_EachScene(i,h_EachScene(i,:)==0)=NaN;
    
    %Response latency
   
   
    Latency=FindOutLatency(nanmean(MeanPSTH_Time),h_EachScene(i,:),ResponseLatencyThreshold);
    if ~isempty(Latency)
        ResponseLatencyForEachScene(i)=Latency;
        
    else
       ResponseLatencyForEachScene(i)=NaN; 
    end
   
    
    %Response Duration   
    ResponseDurationScene=nanmean(MeanPSTH_Time);
    
    ResponseDurationScene=ResponseDurationScene(h_EachScene(i,:)==1);
    ResponseDurationScene=ResponseDurationScene(ResponseDurationScene<=min(RelativeTimeMarker(:,2)) & ResponseDurationScene>=ResponseLatencyForEachScene(i));
     
    TotalResponseDurationScene(i)=length(ResponseDurationScene)*StepSize;
    TotalEventDuration=min(RelativeTimeMarker(:,2));
    PropResponseDurationScene(i)=TotalResponseDurationScene(i)/TotalEventDuration*100;
    
    %Peak Firing Rate in the first significant interval
    Significant_TimePoint=ResponseDurationScene(ResponseDurationScene>=ResponseLatencyForEachScene(i));
    FirstIntervalEnd=find(diff(Significant_TimePoint)>StepSize,1,'first');
    if isempty(FirstIntervalEnd) || FirstIntervalEnd>length(Significant_TimePoint)
        FirstIntervalEnd=length(Significant_TimePoint);
        
    end
    
    
    FirstSigInterval=[1:FirstIntervalEnd];
    
    FR_Sig=SpikeFiringRateSceneOn_mean(i,nanmean(MeanPSTH_Time)>=ResponseLatencyForEachScene(i));
    
    if ~isempty(FR_Sig) && ~isempty(FirstSigInterval)
        
    [PeakFiringScene(i),loc]=max(FR_Sig(FirstSigInterval));
       
    
    
    PeakFiringTimeScene(i)=Significant_TimePoint(loc);
    
    %First Interval Duration
    FirstIntervalDurationScene(i)=diff(Significant_TimePoint([FirstSigInterval(1),FirstSigInterval(end)]));
    else
        PeakFiringScene(i)=NaN;
        PeakFiringTimeScene(i)=NaN;
        FirstIntervalDurationScene(i)=NaN;
            
    end
    
    
    %PeakFiringTimeScene(i)=ResponseDurationScene(find(ResponseDurationScene==PeakFiringScene(i),1,'first'));
    

    %SpikeFiringRateSceneOff_mean(i,:)=mean(SpikeFiringRateOff(trialindex,:),1);
    %SpikeFiringRateSceneOff_sem(i,:)=std(SpikeFiringRateOff(trialindex,:),[],1)/sqrt(sum(trialindex));
    
    

end
%Add spon
SceneStr{i+1}='SponAll';

%Statistics
%Anova for different scene responses

[pSceneDiff,tbl,stats]=anova1(MeanFRRaw,SceneIndex(:,3),'off');
%Anova between spon
[pSponDiff,tbl,stats]=anova1(Mean_Spon_ForEachTrial,SceneIndex(:,3),'off');

%Calculate scene selective index
SceneSelectiveIndex=(max(FR_Scene_mean)-min(FR_Scene_mean))/(max(FR_Scene_mean)-min(FR_Scene_mean)+2*sqrt(sum(FR_Scene_SSE)/(TrialNum-NumOfCondition)));

SponSelectiveIndex=(max(Spon_ForEachScene_Mean)-min(Spon_ForEachScene_Mean))/(max(Spon_ForEachScene_Mean)-min(Spon_ForEachScene_Mean)+2*sqrt(sum(FR_Spon_SSE)/(TrialNum-NumOfCondition)));

%SSI significance dynamics
pSceneDiffWithTime=anova1_group(MeanPSTH_FR,SceneIndex(:,3),'off');
h_sig_ssi_withtime=floor(pSceneDiffWithTime<0.05);
h_sig_ssi_withtime(h_sig_ssi_withtime==0)=NaN;


SSI_time=(max(SpikeFiringRateSceneOn_mean,[],1)-min(SpikeFiringRateSceneOn_mean,[],1))./(max(SpikeFiringRateSceneOn_mean,[],1)-min(SpikeFiringRateSceneOn_mean,[],1)+2*sqrt(sum(FR_SSE_Scene_Time,1)/(TrialNum-NumOfCondition)));

if sum(isnan(h_sig_ssi_withtime))==length(h_sig_ssi_withtime)
    %Not anyone significant
    SSI_Latency=NaN;
    PeakSSI=NaN;
    PeakSSITime=NaN;
    FirstIntervalDurationSSI=NaN;
    SSI_sigtime_Whole=NaN;
    PropSSISig=NaN;
    
else


%SSI latency
SSI_Latency=FindOutLatency(nanmean(MeanPSTH_Time),h_sig_ssi_withtime,ResponseLatencyThreshold);

%SSI peak

%SSI first interval
SSI_timepoint=nanmean(MeanPSTH_Time);



SSI_Significant_TimePoint=SSI_timepoint((SSI_timepoint>=SSI_Latency)&h_sig_ssi_withtime==1);


SSI_FirstIntervalEnd=find(diff(SSI_Significant_TimePoint)>StepSize,1,'first');

if isempty(SSI_FirstIntervalEnd) || SSI_FirstIntervalEnd>length(SSI_Significant_TimePoint)
    
   SSI_FirstIntervalEnd=length(SSI_Significant_TimePoint);
        
end
    
SSI_FirstSigInterval=[1:SSI_FirstIntervalEnd];
    
SSI_Sig=SSI_time(nanmean(MeanPSTH_Time)>=SSI_Latency);
    
[PeakSSI,loc]=max(SSI_Sig(SSI_FirstSigInterval));
    
PeakSSITime=SSI_Significant_TimePoint(loc);
    
 %SSI First Interval Duration
FirstIntervalDurationSSI=diff(SSI_Significant_TimePoint([SSI_FirstSigInterval(1),SSI_FirstSigInterval(end)]));
    
%Total SSI Interval Duration

SSI_sig_Whole=SSI_time(SSI_timepoint>=SSI_Latency &SSI_timepoint<=min(RelativeTimeMarker(:,2))&h_sig_ssi_withtime==1);

SSI_sigtime_Whole=length(SSI_sig_Whole)*StepSize;
PropSSISig=SSI_sigtime_Whole/TotalEventDuration*100;

end


%Sort the scene according to the firing rate
[FR_Scene_mean_Sorted,Index_Sort]=sort(FR_Scene_mean,'descend');
PreferScene=UniqueScene(Index_Sort(1));
SceneStrSorted=SceneStr(1:end-1);
SceneStrSorted=SceneStrSorted(Index_Sort);
SceneStrSorted=[SceneStrSorted,SceneStr(end)];

%%
%If show figures 
if ShowFigureFlag
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Bar Plot for the averagre firing rate%%%%%%%%%%%%%%%%%%%%%%%%


fig99=PrepareFigure(99,'w',[200,100, 800,600],'Name','Bar plot for response average');

%Scene response
axes('Position',[0.12,0.5,0.8,0.4])

g=bar([FR_Scene_mean_Sorted,Mean_Spon(1)],'LineWidth',3);
g.FaceColor = 'flat';
g.EdgeColor = 'flat';
ColorChosen= Color(1:NumOfCondition+1,:);

g.CData = ColorChosen([Index_Sort,NumOfCondition+1],:);

xlim([0.2,NumOfCondition+1.8]);


XData=g.XData;
yData=g.YData;

XData=XData(1:end-1);
yData=yData(1:end-1);

%Add significance level
text(XData(h_SceneFromSpon==1)-0.045,yData(h_SceneFromSpon==1)+max(FR_Scene_mean_Sorted+FR_Scene_sem(Index_Sort))*0.1,'*','FontSize',40)
text(XData(p_SceneFromSpon<0.01)-0.045,yData(p_SceneFromSpon<0.01)+max(FR_Scene_mean_Sorted+FR_Scene_sem(Index_Sort))*0.2,'*','FontSize',40)
text(XData(p_SceneFromSpon<0.001)-0.045,yData(p_SceneFromSpon<0.001)+max(FR_Scene_mean_Sorted+FR_Scene_sem(Index_Sort))*0.3,'*','FontSize',40)

ylim([0,max(FR_Scene_mean_Sorted+FR_Scene_sem(Index_Sort))*1.5]);
%}
hold on
errorbar([FR_Scene_mean_Sorted,Mean_Spon(1)],[FR_Scene_sem(Index_Sort),Mean_Spon(2)],'.k','LineWidth',3);

set(gca,'LineWidth',3)
xticks(1:NumOfCondition+1);
xticklabels(SceneStrSorted);
set(gca,'FontSize',15,'FontWeight','Bold')

ylabel(['\fontsize{20}Firing Rate to scene',newline, '\fontsize{20}(Spikes/s)']);
box off
h=gca; h.XAxis.TickLength = [0 0];

%Show scene selective index and p value
%annonation()

SSI_str=sprintf('SSI: %2.2f\np(ANOVA): %2.2e',SceneSelectiveIndex,pSceneDiff);
ssia=annotation('textbox', [0.8, 0.8, 0.7, 0.04], 'string', SSI_str,'FontSize',15,'FontWeight','Bold','LineStyle','none','FitBoxToText','on');
if pSceneDiff<0.05
    ssia.Color='r';
end

%Add data information
data_msg=sprintf(' %s, SpikeChannel: %d, TrialNum: %d ',FileName,SpikeChannel,TrialNum);
%annotation('textbox', [0.15, 0.95, 0.7, 0.04], 'string', data_msg,'FontSize',30,'FontWeight','Bold','LineStyle','none','FitBoxToText','on');
title(data_msg,'FontSize',20,'FontWeight','Bold');





%Spontaneous response
axes('Position',[0.12,0.1,0.8,0.3])

gs=bar([Spon_ForEachScene_Mean(Index_Sort),Mean_Spon(1)],'LineWidth',3);
gs.FaceColor = 'flat';
gs.EdgeColor = 'flat';

gs.CData = ColorChosen([Index_Sort,NumOfCondition+1],:);

hold on
errorbar([Spon_ForEachScene_Mean(Index_Sort),Mean_Spon(1)],[Spon_ForEachScene_SEM(Index_Sort),Mean_Spon(2)],'.k','LineWidth',3);

set(gca,'LineWidth',3)
xticks(1:NumOfCondition+1);
xticklabels(SceneStrSorted);
set(gca,'FontSize',15,'FontWeight','Bold')

ylabel(['\fontsize{20}Spontaneous Firing Rate',newline, '\fontsize{20}(Spikes/s)']);
box off
h=gca; h.XAxis.TickLength = [0 0];
xlim([0.2,NumOfCondition+1.8])

%Show spon selective index and p value
%annonation()

SponSI_str=sprintf('SponSI: %2.2f\np(ANOVA): %2.0e',SponSelectiveIndex,pSponDiff);
sponsia=annotation('textbox', [0.8, 0.4, 0.7, 0.04], 'string', SponSI_str,'FontSize',15,'FontWeight','Bold','LineStyle','none','FitBoxToText','on');
if pSponDiff<0.05
    sponsia.Color='r';
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%First plot the overall response across all scenes%%%%%%%%%%%%%%%%%%%%%%%%

fig100=PrepareFigure(100,'w',[200,100, 1000,1000],'Name','Average Scene Response');%,'Name','SceneTuning');
%axes('position',[0.15,0.2,0.65,0.7])

shadedErrorBar(nanmean(MeanPSTH_Time),nanmean(MeanPSTH_FR),nanstd(MeanPSTH_FR)./sqrt(sum(~isnan(MeanPSTH_FR),1)),'lineprops',{[1,0,0]},'transparent',1);
hold on

%Add significance level
sig=area(nanmean(MeanPSTH_Time),h_All*0.5,'FaceColor','#FFC107','EdgeColor','#FFC107','FaceAlpha',0.5,'EdgeAlpha',0.5);
xlabel('Time from scene onset(ms)');
ylabel('Firing Rate (Hz)');


%Add vertical lines to indicate event
y_max=max(nanmean(MeanPSTH_FR)+nanstd(MeanPSTH_FR)./sqrt(sum(~isnan(MeanPSTH_FR),1)));
y_min=min(nanmean(MeanPSTH_FR)-nanstd(MeanPSTH_FR)./sqrt(sum(~isnan(MeanPSTH_FR),1)));

%area(RelativeTimeMarker(:,2),20*ones(length(RelativeTimeMarker(:,1)),1),'FaceColor','k','FaceAlpha',0.1)

plot([0,0],[y_min,y_max],'--k','LineWidth',3)
plot([min(RelativeTimeMarker(:,2)),min(RelativeTimeMarker(:,2))],[y_min,y_max],'--k','LineWidth',3)

%Add vertical lines to indicate response latency
plot([ResponseLatency,ResponseLatency],[y_min,y_max],'LineStyle','--','Color','#FFC107','LineWidth',3)
%Show response latency
RLatencyStr=sprintf('Latency: %d ms',ResponseLatency);
text(ResponseLatency,y_max*1.05,RLatencyStr,'FontSize',25,'FontWeight','Bold');



set(gca,'LineWidth',3,'FontSize',30);

FRPara=sprintf('Step %dms,\nBin %dms',StepSize,BinWidth);
g=legend(FRPara,['Significantly',newline,' above spon'],'FontSize',15,'FontWeight','Bold');
g.Position=[0.87   0.7    0.08    0.0880];

%RLatencyStr=sprintf('Latency: %d ms',ResponseLatency);
%annotation('textbox',[0.6,0.2,0.5,0.5],'String',RLatencyStr,'FitBoxToText','on','LineStyle','none','FontSize',30)

set(findall(gca, 'Type', 'Line'),'LineWidth',3);

%Add data information
data_msg=sprintf(' %s, SpikeChannel: %d, TrialNum: %d ',FileName,SpikeChannel,TrialNum);
annotation('textbox', [0.15, 0.95, 0.7, 0.04], 'string', data_msg,'FontSize',30,'FontWeight','Bold','LineStyle','none','FitBoxToText','on');
%title(data_msg,'FontSize',20,'FontWeight','Bold');

%Show total response duration
DurationStr=sprintf('Response \nDuration: %dms',TotalResponseDuration);
DurationPropStr=sprintf('Duration occupy\nProportion: %2.2f%%',PropResponseDuration);

annotation('textbox', [0.81   0.58    0.08    0.0880], 'string', DurationStr,'FontSize',20,'FontWeight','Bold','LineStyle','none','FitBoxToText','on')
annotation('textbox', [0.81   0.48    0.08    0.0880], 'string', DurationPropStr,'FontSize',20,'FontWeight','Bold','LineStyle','none','FitBoxToText','on')



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Plot PSTH for each scene%%%%%%%%%%%%%%%%%%%%%%%%

%subplot(1,3,2)

fig101=PrepareFigure(101,'w',[100,200, 1000,1000],'Name','PSTH for each scene(Together)');%,'Name','SceneTuning');
fig102=PrepareFigure(102,'w',[300,200, 1000,1000],'Name','PSTH for each scene(Separate)');%,'Name','SceneTuning');


for i=1:NumOfCondition
    %figure(101)
set(groot,'CurrentFigure',fig101);
shadedErrorBar(nanmean(MeanPSTH_Time),SpikeFiringRateSceneOn_mean(i,:),SpikeFiringRateSceneOn_sem(i,:),'lineprops',{Color(i,:)},'transparent',1,'patchSaturation',0.3)
%plot(TimeSequence(i,:),SpikeFiringRateSceneOn_mean(i,:),'Color',Color(i,:))
hold on

set(groot,'CurrentFigure',fig102);
subplot(ceil(NumOfCondition/3),3,i);
shadedErrorBar(nanmean(MeanPSTH_Time),SpikeFiringRateSceneOn_mean(i,:),SpikeFiringRateSceneOn_sem(i,:),'lineprops',{Color(i,:)},'transparent',1,'patchSaturation',0.3)

hold on

ylim([-1,max(max(SpikeFiringRateSceneOn_mean+SpikeFiringRateSceneOn_sem))*1.02])
xlim([min(min(nanmean(MeanPSTH_Time)))-10,max(max(nanmean(MeanPSTH_Time)))+10])
plot([0,0],[-1,max(max(SpikeFiringRateSceneOn_mean+SpikeFiringRateSceneOn_sem))*1.02],'-.k','LineWidth',3);
plot([RelativeTimeMarkerForeachScene(i),RelativeTimeMarkerForeachScene(i)],[-1,max(max(SpikeFiringRateSceneOn_mean+SpikeFiringRateSceneOn_sem))*1.02],'--k','LineWidth',3)

%Add significance

area(nanmean(MeanPSTH_Time),-2*h_EachScene(i,:),'FaceColor',Color(i,:),'EdgeColor',Color(i,:),'FaceAlpha',0.5,'EdgeAlpha',0.5);

%Add response latency
rlg=plot([ResponseLatencyForEachScene(i),ResponseLatencyForEachScene(i)],[-1,max(max(SpikeFiringRateSceneOn_mean+SpikeFiringRateSceneOn_sem))]);
rlg.LineStyle='--';
rlg.Color=Color(i,:);

%Show response latency
RLatencyStrScene=sprintf('Latency:\n %d ms',ResponseLatencyForEachScene(i));
text(RelativeTimeMarkerForeachScene(i)/3,max(max(SpikeFiringRateSceneOn_mean+SpikeFiringRateSceneOn_sem))*0.9,RLatencyStrScene,'FontSize',15,'FontWeight','Bold');

%Show peak firing time
PTimeStrScene=sprintf('Peak Time:\n %d ms',PeakFiringTimeScene(i));
text(RelativeTimeMarkerForeachScene(i)/3,max(max(SpikeFiringRateSceneOn_mean+SpikeFiringRateSceneOn_sem))*0.7,PTimeStrScene,'FontSize',15,'FontWeight','Bold');

rlg=plot([PeakFiringTimeScene(i),PeakFiringTimeScene(i)],[-1,max(max(SpikeFiringRateSceneOn_mean+SpikeFiringRateSceneOn_sem))]);
rlg.LineStyle='--';
rlg.Color=Color(i,:);

%Show significant proportion

PropSigSceneStr=sprintf('Sig Prop:\n %2.2f%% ',PropResponseDurationScene(i));
text(RelativeTimeMarkerForeachScene(i)/3,max(max(SpikeFiringRateSceneOn_mean+SpikeFiringRateSceneOn_sem))*0.5,PropSigSceneStr,'FontSize',15,'FontWeight','Bold');





%drawaxis(gca, 'x', 0, 'movelabel', 1)
%h_EachScene(i,:)

title(SceneStr{i});

%Add xlabel and ylabel
if mod(i,3)==1
    ylabel('Spikes/s');
end
if i==NumOfCondition-1
    
    xlabel('Time from scene onset (ms)');
end


set(findall(gca, 'Type', 'Line'),'LineWidth',3);
set(gca,'LineWidth',3,'FontSize',20,'FontWeight','Bold');


end





set(groot,'CurrentFigure',fig101);
ylim([-2,max(max(SpikeFiringRateSceneOn_mean+SpikeFiringRateSceneOn_sem))*1.02])
xlim([min(min(nanmean(MeanPSTH_Time)))-10,max(max(nanmean(MeanPSTH_Time)))+10])
plot([0,0],[-2,max(max(SpikeFiringRateSceneOn_mean+SpikeFiringRateSceneOn_sem))*1.02],'-.k','LineWidth',3);
plot([min(RelativeTimeMarker(:,2)),min(RelativeTimeMarker(:,2))],[-2,max(max(SpikeFiringRateSceneOn_mean+SpikeFiringRateSceneOn_sem))*1.02],'--k','LineWidth',3)

%Add significance

area(nanmean(MeanPSTH_Time),-2*(h_sig_ssi_withtime),'FaceColor','b','EdgeColor','none','FaceAlpha',0.2,'EdgeAlpha',0.2);



%plot([mean(FixOnTime) mean(FixOnTime)],[0,max(max(SpikeFiringRateSceneOn_mean+SpikeFiringRateSceneOn_sem))*1.02],'--k');
%plot([min(TgtOnTime) min(TgtOnTime)],[0,max(max(SpikeFiringRateSceneOn_mean+SpikeFiringRateSceneOn_sem))*1.02],'--k');


%ld.String={string(UniqueScene)};


xlabel('Time from Scene Onset')
ylabel('Spikes/s')
%{
yyaxis right
plot(nanmean(MeanPSTH_Time),SSI_time,'-b');
ylabel('SSI')
ax = gca;
ax.YColor='b';
%}
%set(gca,'Color','b');
ld=legend([string(UniqueScene)]);
%ld=legend([string(UniqueScene);'SSI']);
ld.Box='off';
ld.Position=[0.16   0.7    0.0854    0.0880];

title(data_msg,'FontSize',20,'FontWeight','Bold');
set(findall(gca, 'Type', 'Line'),'LineWidth',3);

set(gca,'LineWidth',3,'FontSize',20,'FontWeight','Bold');

%{
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Plot PSTH for the SSI . %%%%%%%%%%%%%%%%%%%%%%%%
PrepareFigure(103,'w',[200,100, 1000,1000],'Name','SSI dynamics');%,'Name','SceneTuning');
axes('Position',[0.15,0.15,0.6,0.7])
%
plot(nanmean(MeanPSTH_Time),SSI_time,'-b');
hold on

%Add significance level
sig=area(nanmean(MeanPSTH_Time),h_sig_ssi_withtime*0.01,'FaceColor','b','EdgeColor','none','FaceAlpha',0.5,'EdgeAlpha',0.5);
xlabel('Time from scene onset(ms)');
ylabel('Scene selective index(SSI)');


%Add vertical lines to indicate event
y_max=max(SSI_time);
y_min=min(SSI_time);

%area(RelativeTimeMarker(:,2),20*ones(length(RelativeTimeMarker(:,1)),1),'FaceColor','k','FaceAlpha',0.1)

plot([0,0],[y_min,y_max],'--k','LineWidth',3)
plot([min(RelativeTimeMarker(:,2)),min(RelativeTimeMarker(:,2))],[y_min,y_max],'--k','LineWidth',3);



%Add vertical lines to indicate response latency
plot([SSI_Latency,SSI_Latency],[y_min,y_max],'LineStyle','--','Color','#FFC107','LineWidth',3);
%Show response latency
SSILatencyStr=sprintf('Latency: %d ms',SSI_Latency);
text(ResponseLatency,y_max*1.05,SSILatencyStr,'FontSize',25,'FontWeight','Bold');



set(gca,'LineWidth',3,'FontSize',30);

FRPara=sprintf('Step %dms,\nBin %dms',StepSize,BinWidth);
g=legend(FRPara,['Significantly',newline,' above spon'],'FontSize',15,'FontWeight','Bold');
g.Position=[0.87   0.7    0.08    0.0880];

%RLatencyStr=sprintf('Latency: %d ms',ResponseLatency);
%annotation('textbox',[0.6,0.2,0.5,0.5],'String',RLatencyStr,'FitBoxToText','on','LineStyle','none','FontSize',30)

set(findall(gca, 'Type', 'Line'),'LineWidth',3);

%Add data information
data_msg=sprintf(' %s, SpikeChannel: %d, TrialNum: %d ',FileName,SpikeChannel,TrialNum);
annotation('textbox', [0.15, 0.95, 0.7, 0.04], 'string', data_msg,'FontSize',30,'FontWeight','Bold','LineStyle','none','FitBoxToText','on');
%title(data_msg,'FontSize',20,'FontWeight','Bold');

%Show total response duration
DurationStr=sprintf('Whole Sig \nDuration: %dms',SSI_sigtime_Whole);
DurationPropStr=sprintf('Duration occupy\nProportion: %2.2f%%',PropSSISig);
PeakSSIStr=sprintf('Peak SSI:\n %2.2f%',PeakSSI);

annotation('textbox', [0.8   0.58    0.08    0.0880], 'string', DurationStr,'FontSize',20,'FontWeight','Bold','LineStyle','none','FitBoxToText','on')
annotation('textbox', [0.8   0.48    0.08    0.0880], 'string', DurationPropStr,'FontSize',20,'FontWeight','Bold','LineStyle','none','FitBoxToText','on')
%annotation('textbox', [0.8   0.38    0.08    0.0880], 'string', PeakSSIStr,'FontSize',20,'FontWeight','Bold','LineStyle','none','FitBoxToText','on')


plot([PeakSSITime,PeakSSITime],[y_min,y_max],'LineStyle','--','Color','b','LineWidth',3)
ylim([0,y_max]);
xlim([-PlotStartOffset,min(RelativeTimeMarker(:,2)+PlotStopOffset)]);

box off
%}
end %end of if show pictures


%Output to files
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
OutputFigureName=sprintf('%s_%s_N%s_C%s',MonkeyName,string(RecordDateOriginal),string(NeuronNum),string(ChannelNum));


ExistFlag=0;
%Load the old file if exist
if exist(OutputFileName)

    load(OutputFileName);
    %Load OutputData into the memory
    if isfield(OutputFileName,'SceneTuningFreeView')
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

OutputData_New.SceneTuningFreeView.TrialType='GoodTrials';

%%%%%%%%%Save overall scene response characteristics%%%%%%%%%%%%%%%%%%%%%%%
%Organize the data
NameStr={'TrialNum','AverageFR_Scene','AverageFR_Scene_EachTrial','SceneSegTimeOffset','AverageFR_Spon','AverageFR_Spon_EachTrial','SponSeg_Marker_OffSet',...
'PSTH_EachTrial','RelativeTimeMarker','AbsoluteTimeMarker' ,'p_SceneAverage_Dynamic','ResponseLatency_Average','ResponseDuration','ResponseProportion',...
'SceneIndex','UniqueScene','StepSize,BinWidth','AverageFR_EachScene','AverageFR_Spon_EachScene','p_SceneFromSpon',...
'PSTH_EachScene','RelativeTimeEnd_EachScene','p_SceneFromSpon_EachScene',...
'ResponseLatency_EachScene','ResponseDuration_EachScene','ResponseProportion_EachScene',...
'PeakFiring_EachScene','PeakFiringTimeScene','FirstIntervalDurationScene',...
'SceneSelectiveIndex','p_SSI','SponSelectiveIndex','p_SponSSI',...
'PSTH_SSI','SSI_Latency','SSI_peak','SSO_PeakTime','SSI_FirstSigInterval','SSI_SigDuration','SSI_PropSig',...
'PSTH_SceneAverage','PeakFiring_SceneAverage','PeakFiringTime_SceneAverage','FirstSigInterval','PreferScene','ResponseLatencyThreshold'

};

DataLib={TrialNum,MeanFR,MeanFRRaw,[PlotStartOffset,PlotStopOffset],Mean_Spon,Mean_Spon_ForEachTrial,[SponEvent,SponEventEnd,PreTimeBuffer,StartOffset],...
MeanPSTH,RelativeTimeMarker,AbosoluteTimeMarker,[h_All',p_All'],ResponseLatency,TotalResponseDuration,PropResponseDuration,...
SceneIndex,UniqueScene,[StepSize,BinWidth],[FR_Scene_mean',FR_Scene_sem'],[Spon_ForEachScene_Mean',Spon_ForEachScene_SEM'],[h_SceneFromSpon',p_SceneFromSpon'],...
[SpikeFiringRateSceneOn_mean',SpikeFiringRateSceneOn_sem'],RelativeTimeMarkerForeachScene,[h_EachScene',p_EachScene'],...
 ResponseLatencyForEachScene,TotalResponseDurationScene,PropResponseDurationScene,...
 PeakFiringScene,PeakFiringTimeScene,FirstIntervalDurationScene,...
 SceneSelectiveIndex,pSceneDiff,SponSelectiveIndex,pSponDiff,...
 SSI_time,SSI_Latency,PeakSSI,PeakSSITime,FirstIntervalDurationSSI,SSI_sigtime_Whole,PropSSISig,...
 [MeanPSTH_FR_mean',MeanPSTH_FR_sem'],PeakFiring,PeakFiringTime,FirstIntervalDuration,PreferScene,ResponseLatencyThreshold

};


DataStamp=containers.Map(NameStr,DataLib);





OutputData_New.SceneTuningFreeView.Task=TaskType;
OutputData_New.SceneTuningFreeView.TaskCode=TaskCode;
OutputData_New.SceneTuningFreeView.DataStamp=DataStamp;
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
    OutputFigureName1=strcat(OutputFigureName,'AverageBar.jpg');
    saveas(fig99,OutputFigureName1);
    OutputFigureName2=strcat(OutputFigureName,'AverageSceneResponse.jpg');
    saveas(fig100,OutputFigureName2);
    OutputFigureName3=strcat(OutputFigureName,'PSTH_EachScene.jpg');
    saveas(fig101,OutputFigureName3);
     OutputFigureName4=strcat(OutputFigureName,'PSTH_EachScene(Separate).jpg');
    saveas(fig102,OutputFigureName4);

end
%}






if ExistFlag
%Compare the old one with the new one
Task_Old=OutputData.SceneTuningFreeView.TaskCode;

     if Task_Old~=TaskCode
        %Not the same task, add the new dataset into the old one
        OutputData(numel(OutputData)+1)=OutputData_New;
     else
         %If the same task, replace the old one with the new one
         OutputData=OutputData_New; 
    
    
     end
else
    %Output the current file
   OutputData=OutputData_New; 
end



if OutputFlag
    cd(OutPath);
    save(OutputFileName,'OutputData');
    disp('Data Exported');
else
    disp('No data export');
end  
%{
%Variabls output to excel
                 %PreferedScene	SceneSelectiveIndex	p(SSI)	SponSelectiveIndex	p(SponSI)	ResponseLatency(Scene)	PeakFiring(Scene)	PeakTime(Scene)	DurationProportion(Scene)	FirstSigInterval(Scene)

VariablesToExcel=[PreferScene,SceneSelectiveIndex,pSceneDiff,SponSelectiveIndex,pSponDiff,ResponseLatency,PeakFiring,PeakFiringTime,PropResponseDuration,FirstIntervalDuration];
VarNum=length(VariablesToExcel);
str_header=[];
for i=1:VarNum
    str_header=[str_header,'%5.2f\t'];
end

VariablesToExcel_str=sprintf(str_header,VariablesToExcel);
disp(VariablesToExcel_str);


%Output to Excel
ExcelPath=strcat(BasicDirectory,'Results');
cd(ExcelPath);
T=readtable('DataArchieve.xlsx');

varnames = T.Properties.VariableNames;



FileNameColumn=T.File__rst_;
CurrFileLoc=find(ismember(FileNameColumn,FileName)==1);

%Find the output label location
%SceneSelectiveIndex	p(SSI)	SponSelectiveIndex	p(SponSI)	ResponseLatency(Scene)	PeakTime(Scene)	DurationProportion(Scene)	FirstSigInterval(Scene)

StartColumnLoc=find(ismember(varnames,'PreferedScene')==1);


T(CurrFileLoc,[StartColumnLoc:StartColumnLoc+VarNum-1])=num2cell(VariablesToExcel);


%writetable(T,filename,'Sheet','MyNewSheet','WriteVariableNames',false);




a=1
%}


%StartOffsetBin=StartOffset/EyeBinWidth;
%StopOffsetBin=StopOffset/EyeBinWidth;
%IntervalInBin=IntervalAfterMarkerStart/EyeBinWidth;




%{
if SaveFigFlag
disp('Saving figures... ...');
%Save Figure
%FigureName=strcat('SceneTuning_',FileName(1:end-4),'_',string(SpikeChannel));
FigureName=strcat('SceneONOffPSTH_',FileName(1:end-4),'_',string(SpikeChannel));

%saveas(103, fullfile(OutPutPath, FigureName), 'jpeg');
saveas(104, fullfile(OutPutPath, FigureName), 'jpeg');

disp('Saved Successfully!');

end

 %}      







end