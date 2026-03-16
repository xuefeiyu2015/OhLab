function  RewardResponseAnalysis(Data,  StartTrial, EndTrial,ShowFigureFlag,OutputFlag);
%Analysis reward responses and rasters; by YXF 07032020
ECode;


SegEvent=[RWDCD,CORRECTCD];%Event Flow of interest
%SegEventEnd=[SHOWFIXCD,FIXOFF,TGTOFF,RWDOFFCD ,NaN,EXTRAFP];%Event end
SegEventEnd=[RWDOFFCD,NaN ];%Event end
SegEventStr={'RewardOn','Spon'};
AmpMatchFlag=0;%1:Amplitude match on ;0:Amplitude Match off


PlotStartOffset=[200,-500];%in ms
PlotStopOffset=[100,1000];%in ms

ExamplePlot=0;%PLots for presentations, turn off for regular use

ResponseLatencyThreshold=3;

%Parameter to define the sliding window to measure firing rates
StepSize=5;%in ms
BinWidth=15;%in ms


%Criterias for selecting saccades
V_Threshold=40;
ContinueBin=5;
ISI_Threshold=20;

%Time bin pre saccade and time bin offset after saccade onset
PlotSacPreBuffer=200;%in ms
PlotSacPostBuffer=400;%in ms




%Color Definition for each fractals
%{
Color=[lines(6);prism(6) ];
Color(2,:)=copper(1)+0.5;
Color=[Color;Color*0.5];
%}
Color={'#648FFF','#FE6100','#009E73','#88CCEE','#CC6677','#785EF0'};
Color_RGB=[100,143,255;
           254,97,0;
           0,158,115;
          136,204,238;
          204,102,119;
          120,94,240]/255;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Unpack all the data base and task settings
ECode;
FileName=Data.FileName;
TaskType=Data.TaskType;
TaskCode=Data.TaskCode;

TrialNum=EndTrial-StartTrial+1;
Selection=StartTrial:EndTrial;


EventChannel=Data.EventChannel(Selection,:);
EventTimeChannel=Data.EventTimeChannel(Selection,:);
EventBin=Data.EventBin(Selection,:);

%Trials with reward 

RewardTrials=sum(EventChannel==BIG_RWD,2)>0;




SpikeChannel=Data.SpikeChannel;
SpikeTime=squeeze(Data.SpikeTimeData(:,:,Selection))';
%SpikeBin=squeeze(Data.SpikeBinData(:,:,Selection))';

PSTH_Event_Each=[];
PSTH_Event=[];
PSTH_Event_Time=[];
PSTH_Event_Each_Time=[];

EventTimeMarkerEach=[];
EventTimeMarker=[];
AverageFiringRate_Event_Mean=[];
AverageFiringRate_Event_Sem=[];
AverageFiringRate_Event_Each=[];

%Findout trials with reward


%Segregate spike data from scene onset to fixation onset

SegmentData=EventSegTool(SegEvent,SegEventEnd,PlotStartOffset,PlotStopOffset,EventChannel,EventTimeChannel,SpikeTime,StepSize,BinWidth);

%Retrive data from the field
%Total average for each event
MeanFR=SegmentData.MeanFR; %[Mean,sem]%Averall mean and sem
MeanFRRaw=SegmentData.MeanFRRaw;%=AverageFiringRate_Event;

MeanFR_Spon_Raw=MeanFRRaw{2};
%Reward response
MeanFR_Reward_Raw=MeanFRRaw{1};
%Screen out nonrewarding trials
MeanFR_Reward_Raw_Rewarding=MeanFR_Reward_Raw(RewardTrials);%Remove nonrewarding(NaN) trials
%MeanFR_Reward_Raw_NonRewarding=MeanFR_Reward_Raw(~RewardTrials);
MeanFR_Spon_Raw_Rewarding=MeanFR_Spon_Raw(RewardTrials);

[h_reward,p_reward]=ttest(MeanFR_Reward_Raw_Rewarding,MeanFR_Spon_Raw_Rewarding);

%Average 
MeanFR_Reward_Mean=nanmean(MeanFR_Reward_Raw_Rewarding);
MeanFR_Reward_Sem=nanstd(MeanFR_Reward_Raw_Rewarding)./sqrt(length(MeanFR_Reward_Raw_Rewarding));

MeanFR_Spon_Mean=nanmean(MeanFR_Spon_Raw_Rewarding);
MeanFR_Spon_Sem=nanstd(MeanFR_Spon_Raw_Rewarding)./sqrt(length(MeanFR_Spon_Raw_Rewarding));


%PSTH for each event
MeanPSTH=SegmentData.MeanPSTH;%=[PSTH_Event_Time',PSTH_Event'];
MeanPSTH_Time_reward=MeanPSTH{1,1};
MeanPSTH_FR_reward=MeanPSTH{1,2};

MeanPSTH_Time_spon=MeanPSTH{2,1};
MeanPSTH_FR_spon=MeanPSTH{2,2};
%PSTH from spike count
MeanPSTH_SpikeCount=SegmentData.PSTH_Event_SpikeCount{1};
MeanPSTH_SpikeCount_Reward=MeanPSTH_SpikeCount(RewardTrials,:);
MeanPSTH_SpikeCount_Reward_mean=sum(MeanPSTH_SpikeCount_Reward,1,'omitnan')./sum(~isnan(MeanPSTH_SpikeCount_Reward),1);

MeanPSTH_SpikeCount_NonReward=MeanPSTH_SpikeCount(~RewardTrials,:);
MeanPSTH_SpikeCount_NonReward_mean=sum(MeanPSTH_SpikeCount_NonReward,1,'omitnan')./sum(~isnan(MeanPSTH_SpikeCount_NonReward),1);


%Relative time bin for reward
RelativeTimeMarker=SegmentData.RelativeTimeMarker{1};%=EventTimeMarker;
%Absolute time marker for reward
AbosoluteTimeMarker=SegmentData.AbsoluteTimeMarker{1};


RasterAfterAlign_reward=SegmentData.RasterAfterAlign{1};

RasterAfterAlign_RewardTrials=RasterAfterAlign_reward(RewardTrials,:);
RasterAfterAlign_NonRewardTrials=RasterAfterAlign_reward(~RewardTrials,:);

MeanPSTH_FR_Reward=MeanPSTH_FR_reward(RewardTrials,:);
MeanPSTH_FR_NonReward=MeanPSTH_FR_reward(~RewardTrials,:);

MeanPSTH_FR_Reward_mean=nanmean(MeanPSTH_FR_Reward,1);
MeanPSTH_FR_Reward_sem=nanstd(MeanPSTH_FR_Reward,[],1)./sqrt(sum(~isnan(MeanPSTH_FR_Reward),1));

MeanPSTH_FR_NonReward_mean=nanmean(MeanPSTH_FR_NonReward,1);
MeanPSTH_FR_NonReward_sem=nanstd(MeanPSTH_FR_NonReward,[],1)./sqrt(sum(~isnan(MeanPSTH_FR_NonReward),1));
MeanPSTH_FR_Time=nanmean(MeanPSTH_Time_reward,1);

%PSTH for spon
MeanPSTH_FR_spon=nanmean(MeanPSTH_FR_spon,1);
SemPSTH_FR_spon=nanstd(MeanPSTH_FR_spon,1)./sqrt(sum(~isnan(MeanPSTH_FR_spon),1));
%P value for PSTH 

[p_Reward_PSTH,h_Reward_PSTH]=ttest_group(MeanPSTH_FR_Reward,MeanFR_Spon_Raw_Rewarding,'both');
[p_NonReward_PSTH,h_NonReward_PSTH]=ttest2_group(MeanPSTH_FR_NonReward,MeanFR_Spon_Raw_Rewarding,'both');
h_Reward_PSTH(h_Reward_PSTH==0)=NaN;
h_NonReward_PSTH(h_NonReward_PSTH==0)=NaN;


Barlength=2;
LineWidth=1;
if ShowFigureFlag
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Figure 1: Average firing rate
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%for each fractal
figtitlestr{1}='RasterAlignOnRewardOnset';
fig{1}=PrepareFigure(130,'w',[50,100, 1200,800],'Name',figtitlestr{1});
StartLoc=1;
%Raster
axes('Position',[0.1,0.7,0.5,0.2])
title('PSTH aligned on reward onset');
SpikeRaterPlot(RasterAfterAlign_RewardTrials,StartLoc,Barlength,'r',LineWidth);
hold on
StartLoc=sum(RewardTrials)+1;
SpikeRaterPlot(RasterAfterAlign_NonRewardTrials,StartLoc,Barlength,'b',LineWidth);
set(gca,'LineWidth',3,'FontWeight','Bold','FontSize',20);
box off
xlim([-PlotStartOffset(1),200+PlotStopOffset(1)])
hold on 
%plot([0,0],[0,length(Selection)],'--k','LineWidth',1.5);
%plot([200,200],[0,length(Selection)],'--k','LineWidth',1.5);
ylabel('Number of trials');

%PSTH
axes('Position',[0.1,0.12,0.5,0.5])
%PSTH for nonreward
shadedErrorBar(MeanPSTH_FR_Time,MeanPSTH_FR_NonReward_mean,MeanPSTH_FR_NonReward_sem,'lineprops',{[0,0,1]},'transparent',1);
hold on
%PSTH for reward
shadedErrorBar(MeanPSTH_FR_Time,MeanPSTH_FR_Reward_mean,MeanPSTH_FR_Reward_sem,'lineprops',{[1,0,0]},'transparent',1);
%Add spon
shadedErrorBar(MeanPSTH_FR_Time,MeanPSTH_FR_spon(1:length(MeanPSTH_FR_Time)),SemPSTH_FR_spon(1:length(MeanPSTH_FR_Time)),'lineprops',{[150,150,150]/255},'transparent',1);


set(findall(gca, 'Type', 'Line'),'LineWidth',3);
set(gca,'LineWidth',3,'FontWeight','Bold','FontSize',20);
box off
xlim([-PlotStartOffset(1),200+PlotStopOffset(1)])
hold on 
yrange=ylim;
plot([0,0],[0,yrange(2)],'--k','LineWidth',1.5);
plot([200,200],[0,yrange(2)],'--k','LineWidth',1.5);

%Add p value for PSTH
%Add significance level
sig_reward=area(MeanPSTH_FR_Time,h_Reward_PSTH*max(yrange(1),5)*0.5,0,'FaceColor',[1,0,0],'EdgeColor',[1,0,0],'FaceAlpha',0.5,'EdgeAlpha',0.5,'LineStyle','none', 'ShowBaseLine', 'off');
sig_nonreward=area(MeanPSTH_FR_Time,-h_NonReward_PSTH*max(yrange(1),5)*0.5,0,'FaceColor',[0,0,1],'EdgeColor',[0,0,1],'FaceAlpha',0.5,'EdgeAlpha',0.5,'LineStyle','none','ShowBaseLine', 'off');
ylim([-max(yrange(1),5)*0.8,yrange(2)]);

legend({'NoReward','Reward','Spon'})
ylabel('Firing Rate(Hz)');
xlabel('Time from reward onset');

%Add data information
data_msg=sprintf(' %s, SpikeChannel: %d, TrialNum: %d ',FileName,SpikeChannel,TrialNum);
annotation('textbox', [0.1, 0.95, 0.7, 0.04], 'string', data_msg,'FontSize',25,'FontWeight','Bold','LineStyle','none','FitBoxToText','on');

%Bar Plot
axes('Position',[0.7,0.12,0.28,0.5]);
gs=bar([MeanFR_Reward_Mean,MeanFR_Spon_Mean],'LineWidth',3);
gs.FaceColor = 'flat';
gs.EdgeColor = 'flat';
gs.CData = [255,0,0;182,185,191]/255;
hold on
errorbar([MeanFR_Reward_Mean,MeanFR_Spon_Mean],[MeanFR_Reward_Sem,MeanFR_Spon_Sem],'.k','LineWidth',3);
set(gca,'LineWidth',3)
xticks(1:2);
xticklabels(SegEventStr);
set(gca,'FontSize',15,'FontWeight','Bold')

ylabel('Average Firing Rate');
box off
h=gca; h.XAxis.TickLength = [0 0];
xlim([0.2,2+0.8]);

%Add p value
yrange=ylim;
      
ctr2 = bsxfun(@plus, gs.XData, [gs.XOffset]');
hold on
plot(ctr2(1:2), [1 1]*yrange(2)*1.1, '-k', 'LineWidth',2);
p_val_str=sprintf('p=%2.3f',p_reward);
t=text(mean(ctr2(1:2))*0.8,yrange(2)*1.2,p_val_str,'FontSize',20,'FontWeight','Bold');

if p_reward<0.05
    t.Color='r';
end

%Add spike count
axes('Position',[0.68,0.7,0.28,0.2])
bar(MeanPSTH_FR_Time,MeanPSTH_SpikeCount_Reward_mean,'FaceColor','r','EdgeColor','r');
set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold')
ylabel('Spike Count');
box off


end

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





%MarkerFolder='DataAnalysis';

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



ExistFlag=0;
%Load the old file if exist
if exist(OutputFileName)

    load(OutputFileName);
    %Load OutputData into the memory
    if isfield(OutputFileName,'RewardResponse')
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
%%%%%%%%%%%%%%%%%%%%OutputData For current analysis
OutputData_New=[];

OutputData_New.RewardResponse.TrialType='GoodTrials';

%%%%%%%%%Save overall scene response characteristics%%%%%%%%%%%%%%%%%%%%%%%
%Organize the data
NameStr={'TrialNum','AverageFR_Reward','AverageFR_Spon','p_value_RewardWithSpon',...
    'PSTH_reward','PSTH_nonreward','PSTH_spon',...
    'p_PSTH_RewardSpon','p_PSTH_NonRewardSpon',...
    'SpikeCount_PSTH_reward','SpikeCount_PSTH_nonreward'
  
};

DataLib={TrialNum,[MeanFR_Reward_Mean,MeanFR_Reward_Sem],[MeanFR_Spon_Mean,MeanFR_Spon_Sem],[h_reward,p_reward],...
    [MeanPSTH_FR_Time',MeanPSTH_FR_Reward_mean',MeanPSTH_FR_Reward_sem'],[MeanPSTH_FR_Time',MeanPSTH_FR_NonReward_mean',MeanPSTH_FR_NonReward_sem'],[MeanPSTH_FR_spon',SemPSTH_FR_spon'],...    
[h_Reward_PSTH',p_Reward_PSTH'],[h_NonReward_PSTH',p_NonReward_PSTH'],...
MeanPSTH_SpikeCount_Reward_mean',MeanPSTH_SpikeCount_NonReward_mean'

};


DataStamp=containers.Map(NameStr,DataLib);





OutputData_New.RewardResponse.Task=TaskType;
OutputData_New.RewardResponse.TaskCode=TaskCode;
OutputData_New.RewardResponse.DataStamp=DataStamp;

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
    
    for j=1:TotalFigure
        OutputFigureName=strcat(OutputFigureName,figtitlestr{j},'.jpg');
        if ~isempty(fig{j})
        saveas(fig{j},OutputFigureName);
        end
        
    end

    disp('Figures have been exported to the neuron folder');

end


if ExistFlag
%Compare the old one with the new one
Task_Old=OutputData.RewardResponse.TaskCode;

     if Task_Old~=TaskCode
        %Not the same task, add the new dataset into the old one
        OutputData.RewardResponse(numel(OutputData.RewardResponse)+1)=OutputData_New.RewardResponse;
     else
         %If the same task,the same protocol,replace the old one with the new one
         
         OutputData.RewardResponse=OutputData_New.RewardResponse;
         
    
    
     end
else
    %Output the current file
   OutputData.RewardResponse=OutputData_New.RewardResponse; 
end



if OutputFlag
    cd(OutPath);
    save(OutputFileName,'OutputData');
    disp('Data Exported');
else
    disp('No data export');
end  



      

end