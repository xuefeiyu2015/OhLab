function RandomSeqFixSeqCompare(Data, StartTrial, EndTrial,ShowFigureFlag,OutputFlag);
%This is to compare the random sequence block and the fix sequence block
%neuron response
%Xuefei Yu, 02252020
%Parameter to define the sliding window to measure firing rates
StepSize=10;%in ms
BinWidth=50;%in ms

%Parameter seting for interval selecting plot period
StartOffset=500;%in ms
StopOffset=50;%in ms
%For scene offset
StartOffset2=500;%in ms
StopOffset2=50;%in ms
IntervalAfterMarkerStart=1000;%in ms;%This for align marker
IntervalAfterMarkerStartOff=1000;%in ms;%For scene offset

%Color for fix and rand
Color_Fix=[	251,97,0]/255;
Color_Rand=[100,143,255]/255;

%Color Definition for each scenes

Color=[lines(6);prism(6);pink(6)];
%Color(2,:)=copper(1)+0.5;
Color=[Color;hsv(10)];
Color(18,:)=[0.1,0.457,0.194];


%{
%For Raster plot 
Barlength=0.5;
LineWidth=1;
%}
%Color Definition for each scenes

%ColorFix=[lines(6);prism(6);pink(6)];
%ColorFix(2,:)=copper(1)+0.5;
%ColorRand=flipud(ColorFix);
ColorFix=repmat([1,0,0],26,1);
ColorRand=repmat([0,0,1],26,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Unpack all the data base and task settings
ECode;
FileName=Data.FileName;
TaskType=Data.TaskType;
TrialNum=EndTrial-StartTrial+1;
Selection=StartTrial:EndTrial;

%FixationPointIndex=Data.FixPointIndex(Selection,:);
%FixLoc=Data.FixLoc(Selection,:,:);


TargetIndex=Data.TargetIndex(Selection,:,:);
BlockType=Data.BlockType(Selection,:,:);%0:fixed sequence;1;random sequence; 0:ScenesetA;1:SceneSetB
%BlockType(BlockType(:,3)>2,:)=999;
BlockType(BlockType(:,3)>2,3)=999;

%%
SceneIndex=Data.SceneIndex(Selection,:);

%EyeChannel_X=Data.EyeChannel_X_L(Selection,:); %Only import the first two channels
%EyeChannel_Y=Data.EyeChannel_Y_L(Selection,:);
%EyeBinWidth=Data.EyeBinWidth;
%EyeTime=Data.EyeTime(Selection,:);

EventChannel=Data.EventChannel(Selection,:);
EventTimeChannel=Data.EventTimeChannel(Selection,:);
EventBin=Data.EventBin(Selection,:);

SpikeChannel=Data.SpikeChannel;
SpikeTime=squeeze(Data.SpikeTimeData(:,:,Selection))';
SpikeBin=squeeze(Data.SpikeBinData(:,:,Selection))';




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Get the PSTH for all neurons
SegEvent=[SCENE_ON];%Event Flow of interest
%SegEventEnd=[SHOWFIXCD,FIXOFF,TGTOFF,RWDOFFCD ,NaN,EXTRAFP];%Event end
%SegEventEnd=[SHOWFIXCD];%Event end
SegEventEnd=[FIXOFF];%Event end
SegEventStr={'SceneOn'};



PlotStartOffset=[200];%in ms
PlotStopOffset=[100];%in ms

ResponseLatencyThreshold=3;

%Parameter to define the sliding window to measure firing rates
%StepSize=10;%in ms
%BinWidth=50;%in ms
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

Selection=StartTrial:EndTrial;

%Abnormal check %Especially For Robin081120FORAGEXPLO4.rst
CheckCode=SCENE_ON;%This code should appear once in each trial
AbnormalTrials=AbnormalTrialCheck(Data,CheckCode);

Selection=Selection(~ismember(Selection,AbnormalTrials));
TrialNum=length(Selection);

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




%#################Separate trials according to scene name and block type

%identify block type

UniqueBlockType=unique(BlockType(:,2:end),'row','stable');

RandFixFlag=unique(UniqueBlockType(:,1));

if length(RandFixFlag)==1
    strFile=sprintf('No random/fix comparison data in the file %s\n',FileName);
    disp(strFile);
    disp('You may need to change to Protocol Scene Tuning');
    return  
end


FixSeqTrials=(BlockType(:,2)==0)|(BlockType(:,2)==10);
RandSeqTrials=(BlockType(:,2)==1)|(BlockType(:,2)==11);


FixSeqSceneName=SceneIndex(FixSeqTrials,3);
RandSeqSceneName=SceneIndex(RandSeqTrials,3);




FixSeqCommonPart=ismember(FixSeqSceneName,unique(RandSeqSceneName));
ConmmonSceneName_FixRand=FixSeqSceneName(FixSeqCommonPart);

RandomSeqCommonPart=ismember(RandSeqSceneName,unique(FixSeqSceneName));
ConmmonSceneName_Rand=RandSeqSceneName(RandomSeqCommonPart);

if isempty(ConmmonSceneName_FixRand)
     strFile=sprintf('No overlap scenes between Fix Seq and Rand Seq %s\n',FileName);
      disp(strFile);
    disp('You may need to change to Protocol Scene Tuning');
    return
end


UniqueSceneName=unique(ConmmonSceneName_FixRand,'stable');

%%%%%%%%%%%%%%%%%%%%%%%First, overall responses%%%%%%%%%%%%%%%%%%%%%%%%%%

MeanFR_Fix_trial_temp=MeanFRRaw(FixSeqTrials);
MeanFR_Rand_trial_temp=MeanFRRaw(RandSeqTrials);


MeanFR_Fix_trial=MeanFR_Fix_trial_temp(FixSeqCommonPart);
MeanFR_Rand_trial=MeanFR_Rand_trial_temp(RandomSeqCommonPart);

AverageFR_Fix_Rand_mean=[nanmean(MeanFR_Fix_trial),nanmean(MeanFR_Rand_trial)];
AverageFR_Fix_Rand_sem=[nanstd(MeanFR_Fix_trial)/sqrt(sum(~isnan(MeanFR_Fix_trial))),nanstd(MeanFR_Rand_trial)/sqrt(sum(~isnan(MeanFR_Rand_trial)))];

[h_Rand_Fix_Whole,p_Rand_Fix_Whole]=ttest2(MeanFR_Fix_trial,MeanFR_Rand_trial);

%%%%%%%%%%%%%%%%%%%%%%%Then, PSTH difference %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

PSTH_Fix_trial_temp=MeanPSTH_FR(FixSeqTrials,:);
PSTH_Rand_trial_temp=MeanPSTH_FR(RandSeqTrials,:);


PSTH_Fix_trial=PSTH_Fix_trial_temp(FixSeqCommonPart,:);
PSTH_Rand_trial=PSTH_Rand_trial_temp(RandomSeqCommonPart,:);

PSTH_Fix_All_Mean=nanmean(PSTH_Fix_trial,1);
PSTH_Rand_All_Mean=nanmean(PSTH_Rand_trial,1);

PSTH_Fix_All_Sem=nanstd(PSTH_Fix_trial,[],1)./sqrt(sum(~isnan(PSTH_Fix_trial),1));
PSTH_Rand_All_Sem=nanstd(PSTH_Rand_trial,[],1)./sqrt(sum(~isnan(PSTH_Rand_trial),1));

[p_RandFix_All,h_RandFix_All]=ttest2_group(PSTH_Fix_trial,PSTH_Rand_trial,'both');
h_RandFix_All(h_RandFix_All==0)=NaN;


PSTH_Time=nanmean(MeanPSTH_Time,1);




%%%%%%%%%%%%%%%%%%%%%%%%% Next Selectivity
UniqueScene=unique(ConmmonSceneName_FixRand);
NumberOfUniqueScene=length(UniqueScene);

%Number of condition in each repetition
[loc,NumInRepFix,RepinSceneFix]=FindNumOfRepetition(ConmmonSceneName_FixRand);
[loc,NumInRepRand,RepinSceneRand]=FindNumOfRepetition(ConmmonSceneName_Rand);



for i=1:NumberOfUniqueScene %Loop over scene 
    
    RandTrialSelection=ConmmonSceneName_Rand==UniqueScene(i);
    FixTrialSelection=ConmmonSceneName_FixRand==UniqueScene(i);
    
    FR_Scene_Rand_Whole=MeanFR_Rand_trial(RandTrialSelection);
    FR_Scene_Fix_Whole=MeanFR_Fix_trial(FixTrialSelection);
    
    PSTH_Fix=PSTH_Fix_trial(FixTrialSelection,:);
    PSTH_Rand=PSTH_Rand_trial(RandTrialSelection,:);
    
    %Repetition
    
    FR_Scene_Rand_Whole_Rep_Cell{i}=FR_Scene_Rand_Whole';
    FR_Scene_Fix_Whole_Rep_Cell{i}=FR_Scene_Fix_Whole';
    
    PSTH_Scene_Rand_Whole_Rep_Cell{i}=PSTH_Rand;
    PSTH_Scene_Fix_Whole_Rep_Cell{i}=PSTH_Fix;
    
    
    
    %Oder in each rep
    
   % RandTrialOrder_temp=mod(temp_finder_rand,NumInRepRand');
   % RandTrialOrder_temp(RandTrialOrder_temp==0)=NumInRepRand(RandTrialOrder_temp==0);
   
   [RandTrialOrder_temp,~,~]=find(RepinSceneRand==UniqueScene(i));
   
    
    FixTrialOrder_temp=repmat(i,sum(FixTrialSelection==1),1);
    %mod(find(FixTrialSelection==1),NumberOfUniqueScene);
    %FixTrialOrder_temp(FixTrialOrder_temp==0)=6;
    RandTrialOrder{i}=RandTrialOrder_temp;
    FixTrialOrder{i}=FixTrialOrder_temp;
    
    
    %Overall SSE for rand and fix
    
    FR_Scene_mean_Fix(i)=mean(FR_Scene_Fix_Whole);
    FR_Scene_std_Fix(i)=nanstd(FR_Scene_Fix_Whole);
    FR_Scene_sem_Fix(i)=std(FR_Scene_Fix_Whole)/sqrt(sum(FixTrialSelection));
      
    FR_Scene_SSE_Fix(i)=sum((FR_Scene_Fix_Whole-mean(FR_Scene_Fix_Whole)).^2); 
    
    FR_Scene_mean_Rand(i)=mean(FR_Scene_Rand_Whole);
    FR_Scene_sem_Rand(i)=std(FR_Scene_Rand_Whole)/sqrt(sum(RandTrialSelection));
    FR_Scene_std_Rand(i)=nanstd(FR_Scene_Rand_Whole); 
    
    FR_Scene_SSE_Rand(i)=sum((FR_Scene_Rand_Whole-mean(FR_Scene_Rand_Whole)).^2); 
    
    
    %For scene-on  PSTH
    PSTH_Fix_mean(i,:)=nanmean(PSTH_Fix,1);
    PSTH_Fix_sem(i,:)=nanstd(PSTH_Fix,[],1)./sqrt(sum(~isnan(PSTH_Fix),1));
    
    PSTH_Rand_mean(i,:)=nanmean(PSTH_Rand,1);
    PSTH_Rand_sem(i,:)=nanstd(PSTH_Rand,[],1)./sqrt(sum(~isnan(PSTH_Rand),1));
    
    
    
    SceneStr{i}=sprintf('Scene%d',UniqueScene(i));
    
    %SSE dynamics
    PSTH_SSE_Fix(i,:)=sum((PSTH_Fix-mean(PSTH_Fix,1)).^2,1); 
    PSTH_SSE_Rand(i,:)=sum((PSTH_Rand-mean(PSTH_Rand,1)).^2,1); 
    
    %SSI p values
    
    
    
    
      
    RelativeTimeMarkerForeachScene(i)=min(RelativeTimeMarker(SceneIndex(:,3)==UniqueScene(i),2));
 
end


%Calculate scene selective index
SceneSelectiveIndex_Fix=(max(FR_Scene_mean_Fix)-min(FR_Scene_mean_Fix))/(max(FR_Scene_mean_Fix)-min(FR_Scene_mean_Fix)+2*sqrt(sum(FR_Scene_SSE_Fix)/(length(ConmmonSceneName_Rand)-NumberOfUniqueScene)));
SceneSelectiveIndex_Rand=(max(FR_Scene_mean_Rand)-min(FR_Scene_mean_Rand))/(max(FR_Scene_mean_Rand)-min(FR_Scene_mean_Rand)+2*sqrt(sum(FR_Scene_SSE_Rand)/(length(ConmmonSceneName_FixRand)-NumberOfUniqueScene)));

SSI_time_Fix=(max(PSTH_Fix_mean,[],1)-min(PSTH_Fix_mean,[],1))./(max(PSTH_Fix_mean,[],1)-min(PSTH_Fix_mean,[],1)+2*sqrt(sum(PSTH_SSE_Fix,1)/(length(ConmmonSceneName_Rand)-NumberOfUniqueScene)));
SSI_time_Rand=(max(PSTH_Rand_mean,[],1)-min(PSTH_Rand_mean,[],1))./(max(PSTH_Rand_mean,[],1)-min(PSTH_Rand_mean,[],1)+2*sqrt(sum(PSTH_SSE_Rand,1)/(length(ConmmonSceneName_FixRand)-NumberOfUniqueScene)));

%Statistics
%Anova for different scene responses

pSSI_fix=anova1_group(PSTH_Fix_trial,ConmmonSceneName_FixRand,'off');
pSSI_rand=anova1_group(PSTH_Rand_trial,ConmmonSceneName_Rand,'off');

h_SSI_fix=floor(pSSI_fix<0.05);
h_SSI_fix(h_SSI_fix==0)=NaN;

h_SSI_rand=floor(pSSI_rand<0.05);
h_SSI_rand(h_SSI_rand==0)=NaN;

%%
%Repetition changes

NumMaxRand=max(cell2mat(cellfun(@numel,FR_Scene_Rand_Whole_Rep_Cell,'Uniform',0)));
%Change cell into a matrix
FR_Scene_Rand_Whole_Rep=cell2mat(cellfun(@(x) [x,ones(1,NumMaxRand-numel(x))]',FR_Scene_Rand_Whole_Rep_Cell,'Uniform',0));

NumMaxFix=max(cell2mat(cellfun(@numel,FR_Scene_Fix_Whole_Rep_Cell,'Uniform',0)));
%Change cell into a matrix
FR_Scene_Fix_Whole_Rep=cell2mat(cellfun(@(x) [x,ones(1,NumMaxFix-numel(x))]', FR_Scene_Fix_Whole_Rep_Cell,'Uniform',0));

%%Calculate mean and sem over scenes


FR_Scene_Rand_Whole_Rep_Mean=nanmean(FR_Scene_Rand_Whole_Rep,2);
FR_Scene_Fix_Whole_Rep_Rep_Mean=nanmean(FR_Scene_Fix_Whole_Rep,2);

FR_Scene_Rand_Whole_Rep_Sem=nanstd(FR_Scene_Rand_Whole_Rep,[],2)/sqrt(NumberOfUniqueScene);
FR_Scene_Fix_Whole_Rep_Sem=nanstd(FR_Scene_Fix_Whole_Rep,[],2)/sqrt(NumberOfUniqueScene);


%%Loop over repetition

%Rand seq

%This time, loop over the order
for i=1:NumberOfUniqueScene
    Order_curr=i;
    FR_Order_Rand=cellfun(@(x,y) y(x==Order_curr),RandTrialOrder,FR_Scene_Rand_Whole_Rep_Cell,'uniform',0);
    FR_Order_Fix=cellfun(@(x,y) y(x==Order_curr),FixTrialOrder,FR_Scene_Fix_Whole_Rep_Cell,'uniform',0);
    
    PSTH_Order_Rand_temp=cellfun(@(x,y) y(x==Order_curr,:),RandTrialOrder,PSTH_Scene_Rand_Whole_Rep_Cell,'uniform',0);
    PSTH_Order_Fix_temp=cellfun(@(x,y) y(x==Order_curr,:),FixTrialOrder,PSTH_Scene_Fix_Whole_Rep_Cell,'uniform',0);
    
     PSTH_Order_Rand=cell2mat(cellfun(@(x) x',PSTH_Order_Rand_temp,'uniform',0))';
     PSTH_Order_Fix=cell2mat(cellfun(@(x) x',PSTH_Order_Fix_temp,'uniform',0))';
     
     
     PSTH_Order_Rand_mean(i,:)=nanmean(PSTH_Order_Rand,1);
     PSTH_Order_Rand_sem(i,:)=nanstd(PSTH_Order_Rand,[],1)./sqrt(sum(~isnan(PSTH_Order_Rand),1));
     
     PSTH_Order_Fix_mean(i,:)=nanmean(PSTH_Order_Fix,1);
     PSTH_Order_Fix_sem(i,:)=nanstd(PSTH_Order_Fix,[],1)./sqrt(sum(~isnan(PSTH_Order_Fix),1));
    
    %{
    %Normalize
    temp=num2cell([1:NumberOfUniqueScene]);
    PSTH_Order_Rand_Norm=cellfun(@(x,y) (x-FR_Scene_mean_Rand((y)))/FR_Scene_std_Rand((y)),PSTH_Order_Rand_temp,temp,'uniform',0);
    
    PSTH_Order_Rand_Norm_mat=cell2mat(cellfun(@(x) x',PSTH_Order_Rand_Norm,'uniform',0))';
    
    PSTH_Order_Rand_Norm_mean(i,:)=nanmean(PSTH_Order_Rand_Norm_mat,1);
    PSTH_Order_Rand_Norm_sem(i,:)=nanstd(PSTH_Order_Rand_Norm_mat,[],1)./sqrt(sum(~isnan(PSTH_Order_Rand_Norm_mat),1));
    
    PSTH_Order_Fix_Norm=cellfun(@(x,y) (x-FR_Scene_mean_Fix((y)))/FR_Scene_std_Fix((y)),PSTH_Order_Fix_temp,temp,'uniform',0);
    
    PSTH_Order_Fix_Norm_mat=cell2mat(cellfun(@(x) x',PSTH_Order_Fix_Norm,'uniform',0))';
    
    PSTH_Order_Fix_Norm_mean(i,:)=nanmean(PSTH_Order_Fix_Norm_mat,1);
    PSTH_Order_Fix_Norm_sem(i,:)=nanstd(PSTH_Order_Fix_Norm_mat,[],1)./sqrt(sum(~isnan(PSTH_Order_Fix_Norm_mat),1));
    %}
    
  %  FR_Scene_mean_Fix(i)=mean(FR_Scene_Fix_Whole);
  %  FR_Scene_std_Fix(i)=nanstd(FR_Scene_Fix_Whole);
  %  FR_Scene_sem_Fix(i)=std(FR_Scene_Fix_Whole)/sqrt(sum(FixTrialSelection));
    
    
    
   % PSTH_Scene_Rand_Whole_Rep_Cell{i}=PSTH_Rand;
   % PSTH_Scene_Fix_Whole_Rep_Cell{i}=PSTH_Fix;
    
    
    FR_Order_Rand_SceneLoc(i,:)=cell2mat(cellfun(@nanmean,FR_Order_Rand,'uniform',0));%Each row indicates each temporal location, each column indicate each scene
    FR_Order_Fix_SceneLoc(i,:)=cell2mat(cellfun(@nanmean,FR_Order_Fix,'uniform',0));%Each row indicates each temporal location, each column indicate each scene
    empty_cell(i)={''};
    
end

%Average over each temporal location

FR_Order_Rand_SceneLoc_mean=nanmean(FR_Order_Rand_SceneLoc,2);
FR_Order_Fix_SceneLoc_mean=nanmean(FR_Order_Fix_SceneLoc,2);

FR_Order_Rand_SceneLoc_sem=nanstd(FR_Order_Rand_SceneLoc,[],2)./sqrt(sum(~isnan(FR_Order_Rand_SceneLoc),2));
FR_Order_Fix_SceneLoc_sem=nanstd(FR_Order_Fix_SceneLoc,[],2)./sqrt(sum(~isnan(FR_Order_Fix_SceneLoc),2));





%%
%If show figures 
if ShowFigureFlag
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Bar Plot for the averagre firing rate between fixed and random sequence condition%%%%%%%%%%%%%%%%%%%%%%%%


fig104=PrepareFigure(104,'w',[200,100, 800,600],'Name','Difference between the random sequence and fixed sequence as a whole');

%Outside pannel: PSTH of fixed sequence and rand sequence

shadedErrorBar(PSTH_Time,PSTH_Fix_All_Mean,PSTH_Fix_All_Sem,'lineprops',{Color_Fix},'transparent',1);%PSTH for fix

hold on
shadedErrorBar(PSTH_Time,PSTH_Rand_All_Mean,PSTH_Rand_All_Sem,'lineprops',{Color_Rand},'transparent',1);%PSTH for rand


%Add vertical lines to indicate event
%y_max=max(PSTH_Fix_All_Mean+PSTH_Fix_All_Sem);
%y_min=min(PSTH_Fix_All_Mean+PSTH_Fix_All_Sem);
yrange=ylim;

%Add significance level
sig=area(PSTH_Time,h_RandFix_All*yrange(1),yrange(1)+2,'FaceColor','#FFC107','EdgeColor','#FFC107','FaceAlpha',0.5,'EdgeAlpha',0.5,'ShowBaseLine','off');
xlabel('Time from scene onset(ms)');
ylabel('Firing Rate (Hz)');

set(findall(gca, 'Type', 'Line'),'LineWidth',3);
set(gca,'LineWidth',3,'FontSize',20);


%area(RelativeTimeMarker(:,2),20*ones(length(RelativeTimeMarker(:,1)),1),'FaceColor','k','FaceAlpha',0.1)

plot([0,0],yrange,'--k','LineWidth',3)
plot([min(RelativeTimeMarker(:,2)),min(RelativeTimeMarker(:,2))],[yrange(1),yrange(2)/2],'--k','LineWidth',3)


set(gca,'LineWidth',3,'FontSize',30);

g=legend({'FixSeq','RandSeq'},'FontSize',15,'FontWeight','Bold','box','off');
%g.Position=[0.6   0.5    0.08    0.0880];

%RLatencyStr=sprintf('Latency: %d ms',ResponseLatency);
%annotation('textbox',[0.6,0.2,0.5,0.5],'String',RLatencyStr,'FitBoxToText','on','LineStyle','none','FontSize',30)

set(findall(gca, 'Type', 'Line'),'LineWidth',3);

%Add data information
data_msg=sprintf(' %s, SpikeChannel: %d, TrialNum: %d ',FileName,SpikeChannel,TrialNum);
annotation('textbox', [0.15, 0.95, 0.7, 0.04], 'string', data_msg,'FontSize',30,'FontWeight','Bold','LineStyle','none','FitBoxToText','on');
%title(data_msg,'FontSize',20,'FontWeight','Bold');




%Average responses difference
axes('Position',[0.6,0.6,0.1,0.2])

b1=bar([0.5,0.8],AverageFR_Fix_Rand_mean);
b1.FaceColor = 'flat';
b1.EdgeColor = 'flat';
b1.CData = [Color_Fix;Color_Rand];
b1.BarWidth=0.4;

xlim([0.25,1]);
xticks(b1.XData);

hold on
errorbar([0.5,0.8],AverageFR_Fix_Rand_mean,AverageFR_Fix_Rand_sem,'Color','k','LineStyle','none','LineWidth',3);
yrange=ylim;

%Add significance level
 

if h_Rand_Fix_Whole==1
    
    %,p_RL_sus_trans
ctr2 = bsxfun(@plus, b1.XData, [b1.XOffset]');
hold on

plot(ctr2(1:2), [1 1]*yrange(2)*1.1, '-k', 'LineWidth',2)
if p_Rand_Fix_Whole<0.05 & p_Rand_Fix_Whole>=0.01
   plot(mean(ctr2(1:2)), yrange(2)*1.2, '*k','MarkerSize',10)
elseif p_Rand_Fix_Whole<0.01 & p_Rand_Fix_Whole>=0.001
    plot(mean(ctr2(1:2))*0.9, yrange(2)*1.2, '*k','MarkerSize',10)
    plot(mean(ctr2(1:2))*1.1, yrange(2)*1.2, '*k','MarkerSize',10)
        
elseif p_Rand_Fix_Whole<0.001
    plot(mean(ctr2(1:2))*0.8, yrange(2)*1.2, '*k','MarkerSize',10)
    plot(mean(ctr2(1:2))*1, yrange(2)*1.2, '*k','MarkerSize',10)
    plot(mean(ctr2(1:2))*1.2, yrange(2)*1.2, '*k','MarkerSize',10)
    
end


end
xticklabels({'Fix','Rand'})
ylabel('Firing Rate(Hz)')

box off

set(gca,'LineWidth',3)
set(findall(gca, 'Type', 'Line'),'LineWidth',3);

set(findall(gcf,'-property','FontSize'),'FontSize',15);
set(findall(gcf,'-property','FontWeight'),'FontWeight','Bold');



end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%PSTH for fixed and random sequence condition separatly%%%%%%%%%%%%%%%%%%%%%%%%


fig105=PrepareFigure(105,'w',[200,100, 800,800],'Name','PSTH for Fix and Rand sequence condition');
fig106=PrepareFigure(106,'w',[200,100, 1000,800],'Name','PSTH for Fix and Rand sequence condition(Separate)');

%Fixed sequence
for i=1:NumberOfUniqueScene
    set(groot,'CurrentFigure',fig105);
  subplot(2,1,1);
  shadedErrorBar(PSTH_Time,PSTH_Fix_mean(i,:),PSTH_Fix_sem(i,:),'lineprops',{Color(i,:)},'transparent',1);%PSTH for fix
  set(findall(gca, 'Type', 'Line'),'LineWidth',3);
  ylim([0,max(max(max(PSTH_Fix_mean)),max(max(PSTH_Rand_mean)))]);
  hold on 
  yrange=ylim;
  
  
  set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
  title('Fix Sequence','FontSize',25);
  ylabel('Firing Rate (Hz)','FontSize',20)
  
  hold on
  subplot(2,1,2);
  shadedErrorBar(PSTH_Time,PSTH_Rand_mean(i,:),PSTH_Rand_sem(i,:),'lineprops',{Color(i,:)},'transparent',1);%PSTH for rand
  hold on
  set(findall(gca, 'Type', 'Line'),'LineWidth',3);
  ylim([0,max(max(max(PSTH_Fix_mean)),max(max(PSTH_Rand_mean)))]);
   yrange=ylim;
  
  
  set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
  title('Random Sequence','FontSize',25)
  ylabel('Firing Rate (Hz)','FontSize',20)
  
  set(groot,'CurrentFigure',fig106);
  
  subplot(ceil(NumberOfUniqueScene/2),2,i);
  shadedErrorBar(PSTH_Time,PSTH_Fix_mean(i,:),PSTH_Fix_sem(i,:),'lineprops',{Color_Fix},'transparent',1);%PSTH for fix
  hold on
  shadedErrorBar(PSTH_Time,PSTH_Rand_mean(i,:),PSTH_Rand_sem(i,:),'lineprops',{Color_Rand},'transparent',1);%PSTH for rand
  set(findall(gca, 'Type', 'Line'),'LineWidth',3);
  ylim([0,max(max(max(PSTH_Fix_mean)),max(max(PSTH_Rand_mean)))]);
   yrange=ylim;
   set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
  title(SceneStr{i},'FontSize',22);
  if mod(i,2)==1
  ylabel('Firing Rate(Hz)')
  end
  if i>NumberOfUniqueScene-3
      xlabel('Time from scene onset(ms)')

  end
  xlim([-200,RelativeTimeMarkerForeachScene(i)]);
  if i==1
      legend({'FixSeq','RandSeq'});
  end
  
  
  
  
end


set(groot,'CurrentFigure',fig105);
legend(SceneStr,'box','off');
xlabel('Time from scene onset(ms)')
%set(gca,'LineWidth',3);


%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%PSTH of SSI for fixed and random sequence condition separatly%%%%%%%%%%%%%%%%%%%%%%%%


fig107=PrepareFigure(107,'w',[200,100, 800,800],'Name','Dynamics of SSI for Fix and Rand sequence condition');

plot(PSTH_Time,SSI_time_Fix,'LineStyle','-','Color',Color_Fix);
hold on
plot(PSTH_Time,SSI_time_Rand,'LineStyle','-','Color',Color_Rand);

%Add significance level
sig_fix=area(PSTH_Time,h_SSI_fix*0.02,0,'FaceColor',Color_Fix,'EdgeColor',Color_Fix,'FaceAlpha',0.6,'EdgeAlpha',0.6,'ShowBaseLine','off');
sig_rand=area(PSTH_Time,-h_SSI_rand*0.02,0,'FaceColor',Color_Rand,'EdgeColor',Color_Rand,'FaceAlpha',0.6,'EdgeAlpha',0.6,'ShowBaseLine','off');

yrange=ylim;
ylim([-0.02,yrange(2)]);
box off

xlabel('Time from scene onset(ms)');
ylabel('Scene Selective Index');

set(findall(gca, 'Type', 'Line'),'LineWidth',3);
set(gca,'LineWidth',3,'FontSize',20);

g=legend({'Fixed Seq','Rand Seq'},'box','off','location','northwest');

%Average responses difference
axes('Position',[0.8,0.8,0.1,0.15])

b1=bar([0.5,0.8],[SceneSelectiveIndex_Fix,SceneSelectiveIndex_Rand]);
b1.FaceColor = 'flat';
b1.EdgeColor = 'flat';
b1.CData = [Color_Fix;Color_Rand];
b1.BarWidth=0.4;

xlim([0.25,1]);
xticks(b1.XData);


xticklabels({'Fix','Rand'})
ylabel('SSI')

box off

set(gca,'LineWidth',3)
set(findall(gca, 'Type', 'Line'),'LineWidth',3);

%set(findall(gcf,'-property','FontSize'),'FontSize',15);
set(findall(gcf,'-property','FontWeight'),'FontWeight','Bold');
%%
%%%%  Repetition changes for each scene

fig108=PrepareFigure(108,'w',[200,100, 800,800],'Name','Fring rate change as a function of repetition for each scene');
%{
plot(FR_Scene_Rand_Whole_Rep,'lineStyle','-','Color',Color_Rand);
hold on
plot(FR_Scene_Fix_Whole_Rep,'lineStyle','-','Color',Color_Fix);
%}
subplot(2,1,1);
shadedErrorBar(1:NumMaxRand,FR_Scene_Rand_Whole_Rep_Mean',FR_Scene_Rand_Whole_Rep_Sem','lineprops',{Color_Rand},'transparent',1);%PSTH for rand
  hold on
shadedErrorBar(1:NumMaxFix,FR_Scene_Fix_Whole_Rep_Rep_Mean,FR_Scene_Fix_Whole_Rep_Sem,'lineprops',{Color_Fix},'transparent',1);%PSTH for fix

    
set(findall(gca, 'Type', 'Line'),'LineWidth',3);
 % ylim([0,max(max(max(PSTH_Fix_mean)),max(max(PSTH_Rand_mean)))]);
 %  yrange=ylim;
  
  
  set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
 
  ylabel('Firing Rate (Hz)','FontSize',20)

xlabel('Repetition of reach scene');
legend({'Rand Seq','Fixed Seq'});

title('Firing rate as a function of repetition (avearged by scene)')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(2,1,2); %Average over the temporal order

errorbar(1:NumberOfUniqueScene,FR_Order_Rand_SceneLoc_mean,FR_Order_Rand_SceneLoc_sem,'LineStyle','-','Marker','s','MarkerSize',5,'Color',Color_Rand);
  hold on
errorbar(1:NumberOfUniqueScene,FR_Order_Fix_SceneLoc_mean,FR_Order_Fix_SceneLoc_sem,'LineStyle','-','Marker','s','MarkerSize',5,'Color',Color_Fix);%PSTH for fix

xticks([1:NumberOfUniqueScene]);
xlim([0.7,NumberOfUniqueScene+0.3]);
set(findall(gcf, '-property','LineWidth'),'LineWidth',3);    

 % ylim([0,max(max(max(PSTH_Fix_mean)),max(max(PSTH_Rand_mean)))]);
 %  yrange=ylim;
  box off
  
  set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
 
  ylabel('Firing Rate (Hz)','FontSize',20)

xlabel('Temporal order');
legend({'Rand Seq','Fixed Seq'});

title('Firing rate as a function of temporal order')

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Firing Rate as a function of both temporal order and the scene %%%%%%%%%%%%%%%%%%%%%%%%


fig109=PrepareFigure(109,'w',[200,100, 1000,500],'Name','Firing Rate as a function of both temporal order and the scene');

%subplot(1,2,1);%Fix sequence
axes('Position',[0.16,0.3,0.3,0.6]);
axis equal;
FR_Order_Rand_SceneLoc_plot=FR_Order_Rand_SceneLoc;
%FR_Order_Rand_SceneLoc_plot(isnan(FR_Order_Rand_SceneLoc_plot))=0;

m1=heatmap(FR_Order_Rand_SceneLoc_plot);
m1.MissingDataColor=[172,174,180]/255;
m1.FontSize=15;

m1.XDisplayLabels=empty_cell;
m1.YDisplayLabels=empty_cell;
title('Rand Sequence');

%contourf(FR_Order_Rand_SceneLoc);
colormap(flipud(hot));




%Add marginal1
axes('Position',[0.17,0.15,0.29,0.1]);%Average over positions for scenes
mar1=errorbar([1:NumberOfUniqueScene],nanmean(FR_Order_Rand_SceneLoc_plot,1),nanstd(FR_Order_Rand_SceneLoc_plot,[],1)./sqrt(sum(~isnan(FR_Order_Rand_SceneLoc_plot),1)));
mar1.LineWidth=3;
mar1.Marker='s';
mar1.Color=Color_Rand;
xlabel('SceneName');
ylabel('FR(Hz)');
box off
set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
xlim([0.8,NumberOfUniqueScene+0.2]);
xticks([1:NumberOfUniqueScene]);
xticklabels(UniqueScene);
%Add marginal2
axes('Position',[0.07,0.33,0.07,0.55]);%Average over scenes for positions
mar2=errorbar([1:NumberOfUniqueScene],nanmean(FR_Order_Rand_SceneLoc_plot,2),nanstd(FR_Order_Rand_SceneLoc_plot,[],2)./sqrt(sum(~isnan(FR_Order_Rand_SceneLoc_plot),2)))
mar2.LineWidth=3;
mar2.Marker='s';
mar2.Color=Color_Rand;


box off
xlim([0.8,NumberOfUniqueScene+0.2]);
xticks([1:NumberOfUniqueScene]);
view([90 90])


xlabel('Temporal Order');
ylabel('FR(Hz)');
set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');



axes('Position',[0.63,0.3,0.3,0.6]);
%subplot(1,2,2);%Fix sequence
axis equal;

FR_Order_Fix_SceneLoc_plot=FR_Order_Fix_SceneLoc;
%FR_Order_Fix_SceneLoc_plot(isnan(FR_Order_Fix_SceneLoc_plot))=0;

m2=heatmap(FR_Order_Fix_SceneLoc_plot,'title','Fixed Sequence');
m2.MissingDataColor=[172,174,180]/255;
m2.FontSize=15;


m2.XDisplayLabels=empty_cell;
m2.YDisplayLabels=empty_cell;

title('Fix Sequence');
%contourf(FR_Order_Fix_SceneLoc);
colormap(flipud(hot));




%Add marginal1
axes('Position',[0.65,0.15,0.27,0.1]);%Average over positions for scenes
mar1=errorbar([1:NumberOfUniqueScene],nanmean(FR_Order_Fix_SceneLoc_plot,1),nanstd(FR_Order_Fix_SceneLoc_plot,[],1)./sqrt(sum(~isnan(FR_Order_Fix_SceneLoc_plot),1)));
mar1.LineWidth=3;
mar1.Marker='s';
mar1.Color=Color_Fix;
xlabel('SceneName');
ylabel('FR(Hz)');
box off
set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
xlim([0.8,NumberOfUniqueScene+0.2]);
xticks([1:NumberOfUniqueScene]);
xticklabels(UniqueScene);
%Add marginal2
axes('Position',[0.55,0.33,0.07,0.55]);%Average over scenes for positions
mar2=errorbar([1:NumberOfUniqueScene],nanmean(FR_Order_Fix_SceneLoc_plot,2),nanstd(FR_Order_Fix_SceneLoc_plot,[],2)./sqrt(sum(~isnan(FR_Order_Fix_SceneLoc_plot),2)));
mar2.LineWidth=3;
mar2.Marker='s';
mar2.Color=Color_Fix;
box off
xlim([0.8,NumberOfUniqueScene+0.2]);
xticks([1:NumberOfUniqueScene]);
set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
view([90 90])


xlabel('Temporal Order');
ylabel('FR(Hz)');

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Firing Rate for the same scene in the same position %%%%%%%%%%%%%%%%%%%%%%%%

fig110=PrepareFigure(110,'w',[200,100, 1000,500],'Name','Firing Rate for the same scene in the same position');

%subplot(1,2,1);%Fix sequence
plot(diag(FR_Order_Rand_SceneLoc_plot),'Marker','s','MarkerSize',5,'LineStyle','-','Color',Color_Rand,'LineWidth',3);
hold on
plot(diag(FR_Order_Fix_SceneLoc_plot),'Marker','s','MarkerSize',5,'LineStyle','-','Color',Color_Fix,'LineWidth',3);
legend({'Rand_Seq','Fix_Seq'});
set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
xlabel('Scene in the same temporal position')
ylabel('Firing Rate(Hz)');
box off;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%PSTH separated by
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%temporal order
fig111=PrepareFigure(111,'w',[200,100, 1000,800],'Name','PSTH for different temporal order');
fig112=PrepareFigure(112,'w',[200,100, 1000,800],'Name','PSTH for Fix and Rand sequence condition(Separate)');

%Fixed sequence
for i=1:NumberOfUniqueScene
    set(groot,'CurrentFigure',fig111);
  subplot(2,1,1);
   
  shadedErrorBar(PSTH_Time,PSTH_Order_Fix_mean(i,:),PSTH_Order_Fix_sem(i,:),'lineprops',{Color(i,:)},'transparent',1,'patchSaturation',0.5);%PSTH for fix
  set(findall(gca, 'Type', 'Line'),'LineWidth',3);
  %ylim([0,max(max(max(PSTH_Order_Fix_Norm_mean)),max(max(PSTH_Order_Rand_Norm_mean)))]);
   ylim([min(min(min(PSTH_Order_Fix_mean)),min(min(PSTH_Order_Rand_mean)))*1.15,max(max(max(PSTH_Order_Fix_mean)),max(max(PSTH_Order_Rand_mean)))*1.15]);
  hold on 
  yrange=ylim;
  
  
  set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
  title('Fix Sequence','FontSize',25);
  ylabel('Normed Firing Rate','FontSize',20)
  
  hold on
  subplot(2,1,2);
  shadedErrorBar(PSTH_Time,PSTH_Order_Rand_mean(i,:),PSTH_Order_Rand_sem(i,:),'lineprops',{Color(i,:)},'transparent',1,'patchSaturation',0.5);%PSTH for rand
  hold on
  set(findall(gca, 'Type', 'Line'),'LineWidth',3);
 % ylim([0,max(max(max(PSTH_Order_Rand_Norm_mean)),max(max(PSTH_Order_Rand_Norm_mean)))]);
  ylim([min(min(min(PSTH_Order_Fix_mean)),min(min(PSTH_Order_Rand_mean)))*1.15,max(max(max(PSTH_Order_Fix_mean)),max(max(PSTH_Order_Rand_mean)))*1.15]);
   yrange=ylim;
  
  
  set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
  title('Random Sequence','FontSize',25)
  ylabel('Normed Firing Rate','FontSize',20)
  
  set(groot,'CurrentFigure',fig112);
  
  subplot(ceil(NumberOfUniqueScene/2),2,i);
  shadedErrorBar(PSTH_Time,PSTH_Order_Fix_mean(i,:),PSTH_Order_Fix_sem(i,:),'lineprops',{Color_Fix},'transparent',1,'patchSaturation',0.5);%PSTH for fix
  hold on
  shadedErrorBar(PSTH_Time,PSTH_Order_Rand_mean(i,:),PSTH_Order_Rand_sem(i,:),'lineprops',{Color_Rand},'transparent',1,'patchSaturation',0.5);%PSTH for rand
  set(findall(gca, 'Type', 'Line'),'LineWidth',3);
  ylim([min(min(min(PSTH_Order_Fix_mean)),min(min(PSTH_Order_Rand_mean)))*1.15,max(max(max(PSTH_Order_Fix_mean)),max(max(PSTH_Order_Rand_mean)))*1.15]);
   yrange=ylim;
   set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
 % title(SceneStr{i},'FontSize',22);
 titleStr{i}=sprintf('Order %d',i);
 title(titleStr{i},'FontSize',22);
  if mod(i,2)==1
  ylabel('Normed Firing Rate')
  end
  if i>NumberOfUniqueScene-3
      xlabel('Time from scene onset(ms)')

  end
  xlim([-200,RelativeTimeMarkerForeachScene(i)]);
  if i==1
      legend({'FixSeq','RandSeq'});
  end
  
  
  
  
end


set(groot,'CurrentFigure',fig111);
legend(titleStr,'box','off');
xlabel('Time from scene onset(ms)')
%set(gca,'LineWidth',3);

%}
%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Export Files
%Output to files
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
    if isfield(OutputFileName,'RandFixSeqCompare')
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

OutputData_New.RandFixSeqCompare.TrialType='GoodTrials';

%%%%%%%%%Save overall scene response characteristics%%%%%%%%%%%%%%%%%%%%%%%
%Organize the data
NameStr={'TrialNum','AverageFR_Fix_Rand','p_Rand_Fix_Whole','SceneSelectiveIndex_Fix','SceneSelectiveIndex_Rand',...
    'PSTH_SSI_Fix','PSTH_SSI_Rand','pSSI_fix','pSSI_rand','PSTH_Time','PSTH_Fix_All_Mean','PSTH_Fix_All_Sem',...
'PSTH_Rand_All_Mean','PSTH_Rand_All_Sem'

};





DataLib={TrialNum,[AverageFR_Fix_Rand_mean',AverageFR_Fix_Rand_sem'],[h_Rand_Fix_Whole',p_Rand_Fix_Whole'],SceneSelectiveIndex_Fix,SceneSelectiveIndex_Rand,...  
SSI_time_Fix,SSI_time_Rand,[h_SSI_fix',pSSI_fix'],[h_SSI_rand',pSSI_rand'],PSTH_Time,PSTH_Fix_All_Mean,PSTH_Fix_All_Sem,...
PSTH_Rand_All_Mean,PSTH_Rand_All_Sem
};


DataStamp=containers.Map(NameStr,DataLib);


OutputData_New.RandFixSeqCompare.Task=TaskType;
OutputData_New.RandFixSeqCompare.TaskCode=TaskCode;
OutputData_New.RandFixSeqCompare.DataStamp=DataStamp;
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
    OutputFigureName1=strcat(OutputFigureName,'AveragePSTHFix_Rand.jpg');
    saveas(fig104,OutputFigureName1);
    OutputFigureName2=strcat(OutputFigureName,'PSTH_EachScene_Fix_Rand.jpg');
    saveas(fig105,OutputFigureName2);
    OutputFigureName3=strcat(OutputFigureName,'PSTH_EachScene_Fix_Rand_Separate.jpg');
    saveas(fig106,OutputFigureName3);
     OutputFigureName4=strcat(OutputFigureName,'SSI_Fix_Rand.jpg');
    saveas(fig107,OutputFigureName4);
    
     OutputFigureName5=strcat(OutputFigureName,'FR_Repetition.jpg');
    saveas(fig108,OutputFigureName5);
    OutputFigureName6=strcat(OutputFigureName,'HeatMap.jpg');
    saveas(fig109,OutputFigureName6);
    OutputFigureName7=strcat(OutputFigureName,'FR_SameLoc_RandFix.jpg');
    saveas(fig110,OutputFigureName7);
     OutputFigureName8=strcat(OutputFigureName,'PSTH_Temporal_Loc.jpg');
    saveas(fig111,OutputFigureName8);
    OutputFigureName9=strcat(OutputFigureName,'PSTH_Temporal_Loc_Separate.jpg');
    saveas(fig112,OutputFigureName9);
    disp('Figures have been exported to the neuron folder');

end

%}





if ExistFlag
%Compare the old one with the new one
Task_Old=OutputData.RandFixSeqCompare.TaskCode;

     if Task_Old~=TaskCode
        %Not the same task, add the new dataset into the old one
        OutputData.RandFixSeqCompare(numel(OutputData.RandFixSeqCompare)+1)=OutputData_New.RandFixSeqCompare;
     else
         %If the same task,the same protocol,replace the old one with the new one
         
         OutputData.RandFixSeqCompare=OutputData_New.RandFixSeqCompare;
         
    
    
     end
else
    %Output the current file
   OutputData.RandFixSeqCompare=OutputData_New.RandFixSeqCompare; 
end



if OutputFlag
    cd(OutPath);
    save(OutputFileName,'OutputData');
    disp('Data Exported');
else
    disp('No data export');
end  

%Variabls output to excel
                 %PreferedScene	SceneSelectiveIndex	p(SSI)	SponSelectiveIndex	p(SponSI)	ResponseLatency(Scene)	PeakFiring(Scene)	PeakTime(Scene)	DurationProportion(Scene)	FirstSigInterval(Scene)
%{
VariablesToExcel=[PreferScene,SceneSelectiveIndex,pSceneDiff,SponSelectiveIndex,pSponDiff,ResponseLatency,PeakFiring,PeakFiringTime,PropResponseDuration,FirstIntervalDuration];
VarNum=length(VariablesToExcel);
str_header=[];
for i=1:VarNum
    str_header=[str_header,'%5.2f\t'];
end

VariablesToExcel_str=sprintf(str_header,VariablesToExcel);
disp(VariablesToExcel_str);
%}



%SceneSelectiveIndex_Fix

%{

UniqueBlockType=unique(BlockType(:,2:end),'row');
NumberOfBlock=size(UniqueBlockType,1);
UniqueSceneNameTotal=unique(SceneIndex(:,3));
AxisNum=[UniqueSceneNameTotal,[1:length(UniqueSceneNameTotal)]'];
if ishandle(130)
close 130;end
figure(130);set(130,'color','w');
set(130,'Position', [100,100, 1000,500], 'Name', 'Scene Onset Response Comparison between fixed sequence blocks and random sequence blocks');

if ishandle(131)
close 131;end
figure(131);set(131,'color','w');
set(131,'Position', [100,100, 1000,500], 'Name', 'Scene Offset Response Comparison between fixed sequence blocks and random sequence blocks');


for i=1:NumberOfBlock
    currBlockType=UniqueBlockType(i,:);
    SceneNameForCurrentBlock=SceneIndex(sum(BlockType(:,2:end)==currBlockType,2)==2,3);
 
    UniqueScene=unique(SceneNameForCurrentBlock);
    for j=1:length(UniqueScene)
        %Calculate and plot the average firing rate according to scenes

       currentScene=UniqueScene(j);
       trialindex=(sum(BlockType(:,2:end)==currBlockType,2)==2)&(SceneIndex(:,3)==currentScene);
       %{
       %For Scene Tuning
       FR_Scene=AverageFiringRateInInterval(trialindex);
  
       FR_Scene_mean(i)=mean(FR_Scene);
       FR_Scene_sem(i)=std(FR_Scene)/sqrt(sum(trialindex));
       FR_Scene_SSE(i)=sum((FR_Scene-mean(FR_Scene)).^2); 
    %}
      %For scene-on and off PSTH
       SpikeFiringRateSceneOn_mean(j,:)=mean(SpikeFiringRate(trialindex,:),1);
       SpikeFiringRateSceneOn_sem(j,:)=std(SpikeFiringRate(trialindex,:),[],1)/sqrt(sum(trialindex));

       SpikeFiringRateSceneOff_mean(j,:)=mean(SpikeFiringRateOff(trialindex,:),1);
       SpikeFiringRateSceneOff_sem(j,:)=std(SpikeFiringRateOff(trialindex,:),[],1)/sqrt(sum(trialindex));
    
    
  
    end
   
    
    %Define axis for each scene
 %  for jj=1:ceil(length(UniqueSceneNameTotal)/3)
       
  
       for  ii=1:length(UniqueScene)
           plotScene=UniqueScene(ii);
           plotAxis=AxisNum(AxisNum(:,1)==plotScene,2);
           
         %  if ii+(jj-1)*3<=length(UniqueScene)
                %Plot SceneOnset Response
           
           SceneOn_Mean=SpikeFiringRateSceneOn_mean(UniqueScene==plotScene,:);
           SceneOn_SEM=SpikeFiringRateSceneOn_sem(UniqueScene==plotScene,:);
          
           figure(130);
           subplot(ceil(length(UniqueSceneNameTotal)/3),3,plotAxis);
           
            if currBlockType(1)==0
               %Fixed sequence
            shadedErrorBar(TimeSequenceOn, SceneOn_Mean,SceneOn_SEM,'lineprops',{ColorFix(find(UniqueSceneNameTotal==plotScene),:)},'transparent',1);
           
            hold on
            plot([0,0],[min(SceneOn_Mean),max(SceneOn_Mean)],'--k')
  
           else
               %Random sequence
              shadedErrorBar(TimeSequenceOn, SceneOn_Mean,SceneOn_SEM,'lineprops',{ColorRand(find(UniqueSceneNameTotal==plotScene),:)},'transparent',1);  
              
              hold on
              plot([0,0],[min(SceneOn_Mean),max(SceneOn_Mean)],'--k')
               
            end
            if plotAxis==2 
            title_msg=sprintf(' %s, SpikeChannel: %d, TrialNum: %d ',FileName,SpikeChannel,TrialNum);
             title(title_msg);
             
            end

             set(gca,'LineWidth',3)

        set(findall(gca, 'Type', 'Line'),'LineWidth',3);

         set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
           ylabel('Spikes/s')
      xlabel('Time from Scene Onset') 
      
      
       %Plot SceneOffset Response
           
           SceneOff_Mean=SpikeFiringRateSceneOff_mean(UniqueScene==plotScene,:);
           SceneOff_SEM=SpikeFiringRateSceneOff_sem(UniqueScene==plotScene,:);
          
           figure(131);
           subplot(ceil(length(UniqueSceneNameTotal)/3),3,plotAxis);
           
            if currBlockType(1)==0
               %Fixed sequence
            shadedErrorBar(TimeSequenceOff, SceneOff_Mean,SceneOff_SEM,'lineprops',{ColorFix(find(UniqueSceneNameTotal==plotScene),:)},'transparent',1);
            
            hold on
            plot([0,0],[min(SceneOff_Mean),max(SceneOff_Mean)],'--k')
  
           else
               %Random sequence
              shadedErrorBar(TimeSequenceOff, SceneOff_Mean,SceneOff_SEM,'lineprops',{ColorRand(find(UniqueSceneNameTotal==plotScene),:)},'transparent',1);  
             
              hold on
              plot([0,0],[min(SceneOff_Mean),max(SceneOff_Mean)],'--k')
               
            end
            if  plotAxis==2 
            title_msg=sprintf(' %s, SpikeChannel: %d, TrialNum: %d ',FileName,SpikeChannel,TrialNum);
             title(title_msg);
             % legend('Fix Seq','Random Seq');
            end

             set(gca,'LineWidth',3)

        set(findall(gca, 'Type', 'Line'),'LineWidth',3);

         set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
           ylabel('Spikes/s')
      xlabel('Time from Scene Offset') 
      
      
            
      %     end
       end
    
%   end
    
    
   
   
   
    
end

   figure(130);
  l1=legend();
  legendString={'Fixed Seq','','Random Seq',''};
  uniqueSeqType=(UniqueBlockType(:,1));
  
  
  
 if length(uniqueSeqType)==2
           l1.String=legendString;
 else
    if uniqueSeqType==0
            l1.String={'Fixed Seq',''};
    else
        l1.String={'Random Seq',''};
    end
end
       
      
%}


end