function  SequenceCompareEachBlock(Data,  StartTrial, EndTrial,ShowFigureFlag,OutputFlag);
%Protocol to analysis scene tuning free view task
%Xuefei Yu 05192021
ECode;
SegEvent=[SCENE_ON];%Event Flow of interest
%SegEventEnd=[SHOWFIXCD,FIXOFF,TGTOFF,RWDOFFCD ,NaN,EXTRAFP];%Event end
%SegEventEnd=[SHOWFIXCD];%Event end
SegEventEnd=[SCENE_OFF ];%Event end
SegEventStr={'SceneOn'};



PlotStartOffset=[1500];%in ms
PlotStopOffset=[1200];%in ms

PlotRangeSceneOn=[-1500,1200];

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
Color=[Color;Color;Color];
%{
ColorForEachMode=[
             
             0.47 0.25 0.80         
             1.00 0.54 0.00 
            
             ];
%}
ColorForEachMode=[
             
             0.47 0.25 0.80         
             1.00 0.54 0.00 
             0.83 0.14 0.14  
            
             ];
AssociateBlock=[0,1,-1];
Blockstr={'FixSceneSeq','RandSceneSeq','ReverseFixSeq'};



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

BlockType=Data.BlockType(Selection,:);
BlockType=BlockType(:,2);
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
UniqueBlockType=unique(BlockType);
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
    
    RelativeTimeMarkerForeachScene(i)=min(RelativeTimeMarker(trialindex,2));
    
    
    SceneStr{i}=sprintf('Scene%d',UniqueScene(i));
    
    
    %Loop over the block type
    for j=1:length(UniqueBlockType)
        
        sel=BlockType==UniqueBlockType(j) & SceneIndex(:,3)==currentScene;
        PSTH_sceneBlock=MeanPSTH_FR(sel,:);
        PSTH_sceneBlock_Mean{i,j}=nanmean(PSTH_sceneBlock,1);
        PSTH_sceneBlock_Sem{i,j}=nanstd(PSTH_sceneBlock,[],1)./sqrt(sum(~isnan(PSTH_sceneBlock)));
        
        %Separate into first half and second half
        firsthalfindex=find(sel==1);
        firsthalf=firsthalfindex(1:floor(sum(sel)/2));
        secondhalf=firsthalfindex(floor(sum(sel)/2)+1:end);
        
        PSTH_sceneBlock_first=MeanPSTH_FR(firsthalf,:);
        PSTH_sceneBlock_second=MeanPSTH_FR(secondhalf,:);
        
        PSTH_sceneBlock_first_Mean{i,j}=nanmean(PSTH_sceneBlock_first,1);
        PSTH_sceneBlock_first_Sem{i,j}=nanstd(PSTH_sceneBlock_first,[],1)./sqrt(sum(~isnan(PSTH_sceneBlock_first)));
        
         PSTH_sceneBlock_second_Mean{i,j}=nanmean(PSTH_sceneBlock_second,1);
        PSTH_sceneBlock_second_Sem{i,j}=nanstd(PSTH_sceneBlock_second,[],1)./sqrt(sum(~isnan(PSTH_sceneBlock_second)));
        
        
    end
    
    
    
    
    
    
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


Time_PSTH_SceneOn=nanmean(MeanPSTH_Time);


%%
FigureStartNum=120;
FigureIndex=1;

%If show figures 
if ShowFigureFlag
 %%%%%%%%%%%%%%%%%%Next figure: Spike Density Function PSTH
 
 figtitlestr{FigureIndex}='SceneResponsesAcrossConditions';
 fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[50,100, 1500,800],'Name',figtitlestr{FigureIndex});   
    
 
 
    for i=1:NumOfCondition
        subplot(ceil(NumOfCondition/2),2,i);
        BlockOrderIndex=[];
        for j=1:length(UniqueBlockType)
            BlockOrderIndex(j)=find(AssociateBlock==UniqueBlockType(j));
            
            Sel=Time_PSTH_SceneOn>PlotRangeSceneOn(1) & Time_PSTH_SceneOn<PlotRangeSceneOn(2);
            
            shadedErrorBar(Time_PSTH_SceneOn(Sel),PSTH_sceneBlock_Mean{i,j}(Sel),PSTH_sceneBlock_Sem{i,j}(Sel),'lineprops',{ColorForEachMode(find(AssociateBlock==UniqueBlockType(j)),:)},'transparent',1,'patchSaturation',0.5);%PSTH for fix
 hold on
 
        end
        set(findall(gca, 'Type', 'Line'),'LineWidth',3);
 
        legend(Blockstr{BlockOrderIndex});
 title(strcat('Scene',num2str(UniqueScene(i))));       
 xlabel('Time from scene onset(ms)');

 ylabel('Firing Rates(Hz)');
            
set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
        
        
    end
    
%PSTH_sceneBlock_Mean

FigureStartNum=FigureStartNum+1;
FigureIndex=FigureIndex+1;



%%%%%%%%%%%%%%%%%%Next figure: Spike Density Function PSTH separete into
%%%%%%%%%%%%%%%%%%frist half and second half
 
 figtitlestr{FigureIndex}='SceneResponsesAcrossConditions_halfcompare';
 fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[50,100, 1500,800],'Name',figtitlestr{FigureIndex});   
    
 
 
    for i=1:NumOfCondition
        %First half
        subplot(NumOfCondition,2,1+(i-1)*2);
        BlockOrderIndex=[];
        for j=1:length(UniqueBlockType)
            BlockOrderIndex(j)=find(AssociateBlock==UniqueBlockType(j));
            
            Sel=Time_PSTH_SceneOn>PlotRangeSceneOn(1) & Time_PSTH_SceneOn<PlotRangeSceneOn(2);
            
            shadedErrorBar(Time_PSTH_SceneOn(Sel),PSTH_sceneBlock_first_Mean{i,j}(Sel),PSTH_sceneBlock_first_Sem{i,j}(Sel),'lineprops',{ColorForEachMode(find(AssociateBlock==UniqueBlockType(j)),:)},'transparent',1,'patchSaturation',0.5);%PSTH for fix
            hold on
         
        end
        set(findall(gca, 'Type', 'Line'),'LineWidth',3);
 
        legend(Blockstr{BlockOrderIndex});
 title(strcat('Scene',num2str(UniqueScene(i)),'-FirstHalf'));       
 xlabel('Time from scene onset(ms)');

 ylabel('Firing Rates(Hz)');
            
set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
        
     %Second half
      subplot(NumOfCondition,2,i*2);
        BlockOrderIndex=[];
        for j=1:length(UniqueBlockType)
            BlockOrderIndex(j)=find(AssociateBlock==UniqueBlockType(j));
            
            Sel=Time_PSTH_SceneOn>PlotRangeSceneOn(1) & Time_PSTH_SceneOn<PlotRangeSceneOn(2);
            
            
             shadedErrorBar(Time_PSTH_SceneOn(Sel),PSTH_sceneBlock_second_Mean{i,j}(Sel),PSTH_sceneBlock_second_Sem{i,j}(Sel),'lineprops',{ColorForEachMode(find(AssociateBlock==UniqueBlockType(j)),:)},'transparent',1,'patchSaturation',0.5);%PSTH for fix
 hold on
        end
        set(findall(gca, 'Type', 'Line'),'LineWidth',3);
 
        legend(Blockstr{BlockOrderIndex});
 title(strcat('Scene',num2str(UniqueScene(i)),'-SecondHalf'));       
 xlabel('Time from scene onset(ms)');

 ylabel('Firing Rates(Hz)');
            
set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');

    end
    
%PSTH_sceneBlock_Mean

FigureStartNum=FigureStartNum+1;
FigureIndex=FigureIndex+1;


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
%}




end
