function FreeViewNeuroResponse(Data, StartTrial, EndTrial,ShowFigureFlag,OutputFlag);

ECode;

TaskType=Data.TaskType;
DataPath= Data.Path;
FileName=Data.FileName;
TaskCode=Data.TaskCode;

Selection=StartTrial:EndTrial;%Select the trials according to the pannel


useReISI = 0;



SegEvent = [ON_STIM];
SegEventEnd = [OFF_STIM];
SegEventStr={'StimOn'};

PlotStartOffset=[100];%in ms
PlotStopOffset=[100];%in ms

%Sequence information from the library
ParaLib=Data.ParaLib;
%{
if ismember('stimOn',keys(ParaLib))
   stimOnFlag = ParaLib('stimOn');
   TrialStimOn = unique(stimOnFlag(stimOnFlag(:,2) == 1,1));
   stimOn = [Selection',0*ones(length(Selection),1)];
   stimOn(ismember(stimOn(:,1),TrialStimOn),2) = 1;
  
   stimOn= stimOn(:,2);
%   stimOn=stimOn(Selection,:);
   
end
%}
if ismember('stimPower',keys(ParaLib))
   stimPower_ori = ParaLib('stimPower');
   stimPower_ori=stimPower_ori(:,2);
   
  
   
end

if ismember('preStimDur,StimDur,PostStimDur',keys(ParaLib))
   stimdur_info = ParaLib('preStimDur,StimDur,PostStimDur');
   stimDur=stimdur_info(:,3);
   
  
   
end



%Manuel Modification of the power by mistake
if ismember(FileName,{'Adams041624VIDEOVIEW1.rst','Adams040324VIDEOVIEW1.rst','Robin061423VIDEOVIEW1.rst','Robin062223VIDEOVIEW1.rst','Robin071323VIDEOVIEW1.rst','Robin072023VIDEOVIEW1.rst','Robin072723VIDEOVIEW1.rst'})
    stimPower = RescuePower(FileName,stimPower_ori);
    
else
    stimPower = stimPower_ori;
end

 

%
if ismember('stimChannel',keys(ParaLib))
   stimChannel_ori = ParaLib('stimChannel');
   stimChannel_ori=stimChannel_ori(:,2);
  
   
else
    if ismember(FileName,{'Robin072523VIDEOVIEW1.rst','Robin072723VIDEOVIEW1.rst','Robin072823VIDEOVIEW1.rst'})
        if ismember(FileName,{'Robin072823VIDEOVIEW1.rst'})

            stimChannel_ori=ones(size(stimPower_ori,1),1);
            stimChannel_ori(81:160)=2;%Red light
             stimChannel_ori(561:end)=2;%Red light

        end

        if ismember(FileName,'Robin072723VIDEOVIEW1.rst')

            stimChannel_ori=ones(size(stimPower_ori,1),1);
            stimChannel_ori(1121:end)=2;%Red light
             

        end
        if ismember(FileName,{'Robin072523VIDEOVIEW1.rst'})

            stimChannel_ori=ones(size(stimPower_ori,1),1);
            stimChannel_ori(1281:end)=2;%Red light
             

        end

    else
        disp('No channel information, autofilled');
        stimChannel_ori=ones(size(stimPower_ori,1),1);

    end
   
end

%Select all stimChannel == 1
SelectChannel = find(stimChannel_ori==1);

Selection = intersect(SelectChannel,Selection);

 stimDur=stimDur(Selection,:);
 stimPower=stimPower(Selection,:);

stimChannel=stimChannel_ori(Selection);

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
        %{
    if SpikeChannel<4
        SpikeTime = SpikeTime;%Correct for the delay of online channel
    end
        %}
    ChannelID = ChannelIDWhole(spk);
    ClusterID = ClusterIDWhole(spk);

    



    
if SpikeChannel<4
    SpikeTime = SpikeTime-3;%Compensate for the online delay;
end




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


SponTimeIn = [-100,0];
%PreStim_Spon = sum(SpikeCount_Stim(:,MeanPSTH_Time_Mean > SponTimeIn(1) & MeanPSTH_Time_Mean<SponTimeIn(2)),2)/((SponTimeIn(2) - SponTimeIn(1))*1000);
%PreStim_Spon = mean(SpikeCount_Stim(:,MeanPSTH_Time_Mean > SponTimeIn(1) & MeanPSTH_Time_Mean<SponTimeIn(2)),2);
PreStim_Spon = mean(SpikeCount(:,MeanPSTH_Time_Mean > SponTimeIn(1) & MeanPSTH_Time_Mean<SponTimeIn(2)),2);

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
if isempty(FixPower)
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
if length(FixPower)>1
    FixPower = FixPower(1);
end

%Reference power for normalization
sel_ref = stimPower == FixPower & stimDur == FixDur;


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
     
    Mean_PSTH_stim_p{i} = nanmean(psth_stim_p,1);
    Sem_PSTH_stim_p{i} = nanstd(psth_stim_p,1)./sqrt(sum(~isnan(psth_stim_p)));
    
    Mean_PSTH_control_p{i} = nanmean(psth_control_p,1);
    Sem_PSTH_control_p{i} = nanstd(psth_control_p,1)./sqrt(sum(~isnan(psth_control_p)));
    
    max_p_tmp = max(max(Mean_PSTH_stim_p{i}+Sem_PSTH_stim_p{i}),max(Mean_PSTH_control_p{i}+Sem_PSTH_control_p{i}));
    if max_p_tmp > maxFR_p
        maxFR_p = max_p_tmp;
        
    end
    
    %Statistics
    
   

    ResponseLatency_p(i) = LatencyWMW(MeanPSTH_Time_Mean,psth_stim_p ,psth_control_p ,ContinousBin);

    %Determine whether there is an effect
    %Another way to define the p of the stimulation effect

    PreStim_Spon_curr =PreStim_Spon(StimTrial & sel);


  

    

    
%CompareTimeBin3 = 50;%in ms


 %p_sig(i) = CheckSignificance(CompareTimeBin3,MeanPSTH_Time_Mean,psth_stim_p,PreStim_Spon_curr);

 StimInterval = nanmean(RelativeTimeMarker(StimTrial & sel,:),1);

 Mean_Response = nanmean(psth_stim_p(:,MeanPSTH_Time_Mean > 0 & MeanPSTH_Time_Mean< StimInterval(2)),2);


 Control_curr = psth_control_p(:,MeanPSTH_Time_Mean > 0 & MeanPSTH_Time_Mean< StimInterval(2));


 Mean_Control = nanmean(Control_curr,2);

 %[h,p]= ttest(Mean_Response,PreStim_Spon_curr);
 [h,p]= ttest2(Mean_Response,Mean_Control);
 p_sig_power(i) =p;

%{
MS_Mean_power(i) = mean(Mean_Response -PreStim_Spon_curr);
MS_Sem_power(i) = std((Mean_Response-PreStim_Spon_curr))/sqrt(length(Mean_Response));
%}
 MS_Mean_power(i) = mean(Mean_Response-Mean_Control(1:length(Mean_Response)));
MS_Sem_power(i) = std((Mean_Response-Mean_Control(1:length(Mean_Response))))/sqrt(length(Mean_Response));

    %Raster;
    
    Raster_Stim_p{i} = Raster(StimTrial & sel,:);
    Raster_Control_p{i} = Raster(ControlTrial & sel,:);
    
    
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
    
    Raster_Stim_d{i} = Raster(StimTrial & sel,:);
    Raster_Control_d{i} = Raster(ControlTrial & sel,:);
    
    ResponseLatency_d(i) = LatencyWMW(MeanPSTH_Time_Mean,psth_stim_d ,psth_control_d ,ContinousBin);

     StimInterval = nanmean(RelativeTimeMarker(StimTrial & sel,:),1);

 Mean_Response = nanmean(psth_stim_d(:,MeanPSTH_Time_Mean > 0 & MeanPSTH_Time_Mean< StimInterval(2)),2);

 Control_curr = psth_control_d(:,MeanPSTH_Time_Mean > 0 & MeanPSTH_Time_Mean< StimInterval(2));


 Mean_Control = nanmean(Control_curr,2);


 %[h,p]= ttest(Mean_Response,PreStim_Spon_curr);
 [h,p]= ttest2(Mean_Response,Mean_Control);
 p_sig_dur(i) =p;

%{
MS_Mean_dur(i) = mean(Mean_Response -PreStim_Spon_curr);
MS_Sem_dur(i) = std((Mean_Response-PreStim_Spon_curr))/sqrt(length(Mean_Response));
%}
 
 MS_Mean_dur(i) = mean(Mean_Response) - mean(Mean_Control);
MS_Sem_dur(i) = std((Mean_Response- mean(Mean_Control)))/sqrt(length(Mean_Response));
 
end
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

%% Summary for varying duration and power
FigureStartNum =FigureStartNum+1;
FigureIndex = FigureIndex+1;

figtitlestr{FigureIndex}='PSTH_StimVSControl_VariPowerDur';
fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[70,100, 1200,600],'Name',figtitlestr{FigureIndex});
subplot(1,2,1)
errorbar(UniquePower,MS_Mean_power,MS_Sem_power,'-sr','LineWidth',3);
hold on
plot(linspace(min(UniquePower),max(UniquePower),100),0*ones(1,100),'--k','LineWidth',3);
box off;

xlabel('Power');
ylabel('Laser Modulation Strength');
set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');

subplot(1,2,2)
errorbar(UniqueDur,MS_Mean_dur,MS_Sem_dur,'-sr','LineWidth',3);
hold on
plot(linspace(min(UniqueDur),max(UniqueDur),100),0*ones(1,100),'--k','LineWidth',3);
box off
xlabel('Duration');
ylabel('Laser Modulation Strength');
set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');





%{
%Figure2 SpikeRaterPlot
FigureStartNum = FigureStartNum+1;
FigureIndex = FigureIndex;

figtitlestr{FigureIndex}='Raster_StimVSControl';
fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[70,100, 1000,400],'Name',figtitlestr{FigureIndex});
%}

 



end %   End of if showFigureFlag

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
    if isfield(OutputData,'FreeViewNeural')
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

OutputData_New.FreeViewNeural.TrialType='GoodTrials';

%%%%%%%%%Save overall scene response characteristics%%%%%%%%%%%%%%%%%%%%%%%
%Organize the data


NameStr={'UniquePower','UniqueDur','MS_Mean_power','MS_Mean_dur',...
    'Mean_PSTH_control_d','Mean_PSTH_control_p',...
    'SpikeChannel','ChannelID'
    };

DataLib={UniquePower,UniqueDur,MS_Mean_power,MS_Mean_dur,...
    Mean_PSTH_control_d,Mean_PSTH_control_p,...
    SpikeChannel,ChannelID



};



DataStamp=containers.Map(NameStr,DataLib);





OutputData_New.FreeViewNeural.Task=TaskType;
OutputData_New.FreeViewNeural.TaskCode=TaskCode;
OutputData_New.FreeViewNeural.DataStamp=DataStamp;
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
Task_Old=OutputData.FreeViewNeural.TaskCode;

     if Task_Old~=TaskCode
        %Not the same task, add the new dataset into the old one
        OutputData(numel(OutputData)+1)=OutputData_New;
     else
         %If the same task, replace the old one with the new one
         %OutputData=OutputData_New; 
         OutputData.FreeViewNeural=OutputData_New.FreeViewNeural; 
    
    
     end
else
    %Output the current file
   OutputData.FreeViewNeural=OutputData_New.FreeViewNeural; 
end



if OutputFlag
    cd(OutPath);
    save(OutputFileName,'OutputData');
    disp('Data Exported');
else
    disp('No data export');
end  
end %End of the if spike channel


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
function p=CheckSignificance(CompareTimeBin,time,SpikeCount_Stim,Control)
    OptoResTime = time(time > 0);

    for i = 1: length(OptoResTime)
    TimeInv_curr = [OptoResTime(i),OptoResTime(i)+CompareTimeBin];
    %SpikeCount_curr = sum(SpikeCount_Stim(:,MeanPSTH_Time_Mean>=TimeInv_curr(1) & MeanPSTH_Time_Mean < TimeInv_curr(2)),2)/CompareTimeBin;
    SpikeCount_tmp= mean(SpikeCount_Stim(:,time>=TimeInv_curr(1) & time < TimeInv_curr(2)),2);
     
    SpikeCount_Store(i) = mean(SpikeCount_tmp);

    [h,p_judge(i)] = ttest(Control,SpikeCount_tmp,'tail','both');

      
    

    end
%}

p= min(p_judge);

%Strength_sig = SpikeCount_Store(find(p_judge2<0.05,1,'first'))-mean(PreStim_Spon);




end
