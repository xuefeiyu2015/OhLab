function EventSegregationAnalysis(Data, StartTrial, EndTrial,ShowFigureFlag,OutputFlag);
%Function to segregate the firing rate aligned with different events
%Xuefei Yu 04012020
ECode;
%Parameter to define the gaussian window to measure firing rates
StepSize=5;%in ms
BinWidth=15;%in ms


%Parameter seting for interval selecting plot period


SegEvent=[SCENE_ON,SHOWFIXCD,TGTCD,RWDCD ,SCENE_OFF,SCENE_OFF];%Event Flow of interest
%SegEventEnd=[SHOWFIXCD,FIXOFF,TGTOFF,RWDOFFCD ,NaN,EXTRAFP];%Event end
SegEventEnd=[SHOWFIXCD,FIXOFF,TGTOFF,RWDOFFCD ,NaN,NaN];%Event end
SegEventStr={'SceneOn','FixOn','TgtOn','RewardOn','SceneOff','Spon'};

PlotStartOffset=[100,100,100,100,100,-1000];%in ms
PlotStopOffset=[100,100,100,100,200,2000];%in ms




%Color Definition for each events

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




EventChannel=Data.EventChannel(Selection,:);
EventTimeChannel=Data.EventTimeChannel(Selection,:);
EventBin=Data.EventBin(Selection,:);

%Automatically chanage the marker

StartCheck=FindOutTime(EventChannel,EventBin,SegEvent(SegEvent==SCENE_ON));
EndCheck=FindOutTime(EventChannel,EventBin,SegEventEnd(SegEvent==SCENE_ON));
EndCheck=EndCheck(:,1);
DiffCheck=EndCheck-StartCheck;
if sum(DiffCheck<100)>ceil(length(StartCheck)/2)
    SegEventEnd(SegEvent==SCENE_ON)=[FIXOFF];%Event end;change from showfix cd to fix off for contigency
    disp('Change to FixOff for FPOn First Mode');
    
end

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


SponFlag=find(contains(SegEventStr,'Spon')==1);

%Trials with reward 

[RewardTrials,~,~]=find(EventChannel,BIG_RWD);

%Align Data with each events
for i=1:length(SegEvent)
    
    %Align with the current event
   AlignTime=FindOutTime(EventChannel,EventTimeChannel,SegEvent(i));
   %Initialize for the mean firing rate
  AverageFiringRateInInterval_Whole=[];
  %AverageFiringRateInInterval=[];
  SpikeFiringRate_temp=[];
  TimeSequence_temp=[];
  EventTime=[];
  
  
   for j=1:size( AlignTime,2)
       AlignTime_Curr=AlignTime(:,j);
       
       if isnan(SegEventEnd(i))
          EndTime_Curr= AlignTime_Curr+PlotStopOffset(i);
          PlotStopOffset(i)=PlotStopOffset(i)/2;
           if i==SponFlag
           %SegEventEnd_Backup=FindOutTime(EventChannel,EventTimeChannel,SCENE_ON);
           %EndTime(isnan(EndTime))=SegEventEnd_Backup(isnan(EndTime));
           
           MeanFROffset=PlotStartOffset(i);
           
           end
      
           
       else
       EndTime=FindOutTime(EventChannel,EventTimeChannel,SegEventEnd(i));
      
       EndTime_Curr=EndTime(:,j);
      
       
        MeanFROffset=0; 
           
           
       end
   
     %  [mindiff,Index]=min(abs(EndTime-AlignTime_Curr),[],2,'omitnan');
     %  EndTime_Curr=EndTime(Index);
    
       
       %Remove the nans
       notanflag=~isnan(AlignTime_Curr);
       AlignTime_Curr=AlignTime_Curr(notanflag);
       
       
       
       EndTime_Curr=EndTime_Curr(notanflag);
       
       %Check whether end time is later than the start time
       
       Check=EndTime_Curr-AlignTime_Curr+MeanFROffset;
       NegativeTime=Check<0;
       
       if sum(NegativeTime)>0
           MessageStr=sprintf('Correct negative duration with the mid value for event %s, term %d',SegEventStr{i},j);
           disp(MessageStr);
           CorrectNum=sum(NegativeTime);
           MessageStr2=sprintf('Number of correction: %d',CorrectNum);
           
       
       AllNormalTime=Check(~NegativeTime);
       AllNormalTime_Sort=sort(AllNormalTime);
       Mid_AllNormalTime_Sort=AllNormalTime_Sort(floor(length(AllNormalTime_Sort)/2));
       EndTime_Curr(NegativeTime)=AlignTime_Curr(NegativeTime)+Mid_AllNormalTime_Sort;
       end
       
       
       
       
       
       
       %Select interval to compute the average firng rate 
       AverageFRTime=[AlignTime_Curr-AlignTime_Curr-MeanFROffset,EndTime_Curr-AlignTime_Curr];%in ms;%Time interval to compute average firing rate within the time bin
       %Time window to plot the PSTH      
       SelectIntervalTime=[-PlotStartOffset(i)*ones(size(AlignTime_Curr,1),1)+MeanFROffset,EndTime_Curr-AlignTime_Curr+PlotStopOffset(i)*ones(size(AlignTime_Curr,1),1)];
       
       %Align spikes according to the marker
       SpikeTimeAfterAlign=SpikeTime(notanflag,:)-AlignTime_Curr;
          
       %Using a gaussian smooth window to calculate firing rate
        [SpikeFiringRate_temp{j},TimeSequence_temp{j}]=CalculateSpikeFiringRate(SpikeTimeAfterAlign,StepSize,BinWidth,SelectIntervalTime);
        
      
        
      %Measure the average firing rate
        SpikeTimeInInterval=SelectSpikeTimeInterval(SpikeTimeAfterAlign,AverageFRTime(:,1),AverageFRTime(:,2));%Return spike count to calculate average firing rate
       SpikeCountInInterval=sum(~isnan(SpikeTimeInInterval),2);
       
       AverageFiringRateInInterval_Whole=[AverageFiringRateInInterval_Whole;SpikeCountInInterval./((AverageFRTime(:,2)-AverageFRTime(:,1))/1000)];
       
       %Store the average firing rate and PSTH
       AverageFiringRate_Event_Each{i}{j}=[SpikeCountInInterval./((AverageFRTime(:,2)-AverageFRTime(:,1))/1000)];%in Hz

       PSTH_Event_Each{i}{j}=SpikeFiringRate_temp{j};
       PSTH_Event_Each_Time{i}{j}=TimeSequence_temp{j};
       
       %Store the time stamp
        EventTime=[EventTime;AverageFRTime];
       
       EventTimeMarkerEach{i}{j}=AverageFRTime;
       
   end
   %Store the PSTH
   temp_length_max=max(cellfun(@(x) size(x,2),SpikeFiringRate_temp));
   
   PSTH_Event{i}=cell2mat(cellfun(@(x) [x,NaN*ones(size(x,1),temp_length_max-size(x,2))],SpikeFiringRate_temp,'uniform',0)');
   PSTH_Event_Time{i}=cell2mat(cellfun(@(x) [x,NaN*ones(size(x,1),temp_length_max-size(x,2))],TimeSequence_temp,'uniform',0)');
   EventTimeMarker{i}=EventTime;
   
 
    %Store the average firing rate   
    AverageFiringRate_Event{i}=AverageFiringRateInInterval_Whole;

    AverageFiringRate_Event_Mean(i)=nanmean(AverageFiringRateInInterval_Whole);
    AverageFiringRate_Event_Sem(i)=nanstd(AverageFiringRateInInterval_Whole)/sqrt(length(AverageFiringRateInInterval_Whole(~isnan(AverageFiringRateInInterval_Whole))));
     
    
    
end


%Statistics for the comparison between the average firing rate between
%different events

Spon_FR=AverageFiringRate_Event{SponFlag};%Select out spontaneous firing rate pool for comparison


for i=1:length(SegEvent)-1
    currFR=AverageFiringRate_Event{find(SegEvent==SegEvent(i))};
    %Use two sample t test
    [Sig_Event(i),p_Event(i)]=ttest2(Spon_FR, currFR);
 

end

%SceneSelectiveIndex=(max(FR_Scene_mean)-min(FR_Scene_mean))/(max(FR_Scene_mean)-min(FR_Scene_mean)+2*sqrt(sum(FR_Scene_SSE)/(TrialNum-NumOfCondition)));
%pScene_Average=anovan(AverageFiringRateInInterval',SceneIndex(:,3),'display','off');


if ShowFigureFlag

if ishandle(103)
    close 103;end
fig103=figure(103);set(103,'color','w');
set(103,'Position', [100,100, 900,600], 'Name', 'Average Firing Rate During Each event');
%axes('position',[0.15,0.7,0.65,0.23])

g=bar(AverageFiringRate_Event_Mean,'LineWidth',3);
g.FaceColor = 'flat';
g.EdgeColor = 'flat';
g.CData = Color(1:length(SegEvent),:);

XData=g.XData;
yData=g.YData;
%Add significance level
text(XData(Sig_Event==1)-0.045,yData(Sig_Event==1)+min(yData)/5,'*','FontSize',40)
text(XData(p_Event<0.01)-0.045,yData(p_Event<0.01)+min(yData)/5*1.5,'*','FontSize',40)
text(XData(p_Event<0.001)-0.045,yData(p_Event<0.001)+min(yData)/5*2,'*','FontSize',40)

ylim([0,max(yData)++min(yData)/5*3])

hold on
errorbar(AverageFiringRate_Event_Mean,AverageFiringRate_Event_Sem,'.k','LineWidth',3);


set(gca,'LineWidth',3)
xticks(1:length(SegEvent));
xticklabels(SegEventStr);
set(gca,'FontSize',20,'FontWeight','Bold')
xlabel('Event Name','FontSize',30)
ylabel(['\fontsize{30}Average Firing Rate',newline, '\fontsize{20}(Spikes/s)']);
box off
h=gca; h.XAxis.TickLength = [0 0];
%xlim([0.5,length(SegEvent)+0.5]);
%Add data information
data_msg=sprintf(' %s, SpikeChannel: %d, TrialNum: %d ',FileName,SpikeChannel,TrialNum);
annotation('textbox', [0.17, 0.9, 0.8, 0.08], 'string', data_msg,'FontSize',20,'FontWeight','Bold','LineStyle','none')

if ishandle(104)
    close 104;end
fig104=figure(104);set(104,'color','w');
set(104,'Position', [200,100, 1000,1000], 'Name', 'PSTH For Each event');
end


for i=1:length(SegEvent)
   
    PSTH_curr=PSTH_Event{i};
    PSTH=nanmean(PSTH_curr);
    PSTH_sem=nanstd(PSTH_curr)./sqrt(sum(~isnan(PSTH_curr),1));
    PSTH_Time_Curr=nanmean(PSTH_Event_Time{i});
    
    TimeMarker=EventTimeMarker{i};
    
    if ShowFigureFlag
     subplot(3,ceil(length(SegEvent)/3),i);
    shadedErrorBar(PSTH_Time_Curr,PSTH,PSTH_sem,'lineprops',{Color(i,:)},'transparent',1);
    hold on
    
    yrange=[min(PSTH-PSTH_sem)*1.1,max(PSTH+PSTH_sem)*1.1];
    if i==SponFlag
    xlim([min(TimeMarker(:,1))+PlotStartOffset(i)/2,min(TimeMarker(TimeMarker(:,2)>1,2))+PlotStopOffset(i)]);
    
    else
    xlim([min(TimeMarker(:,1))-PlotStartOffset(i),min(TimeMarker(TimeMarker(:,2)>1,2))+PlotStopOffset(i)]);
    
    
    end
    
     plot([min(TimeMarker(:,1)),min(TimeMarker(:,1))],yrange,'--k');
   
    plot([min(TimeMarker(TimeMarker(:,2)>1,2)),min(TimeMarker(TimeMarker(:,2)>1,2))],yrange,'--k');
    set(gca,'LineWidth',3)
    set(gca,'FontSize',15,'FontWeight','Bold')
    
    titleStr=sprintf('Align on %s',SegEventStr{i});
    title(titleStr,'FontSize',25);
    
    xstr=sprintf('Time from %s (ms)',SegEventStr{i});
    xlabel(xstr,'FontSize',20);
    ylabel('Firing Rate(Hz)','FontSize',20)
    end
    
    
end

if ShowFigureFlag
%Add data information
data_msg=sprintf(' %s, SpikeChannel: %d, TrialNum: %d ',FileName,SpikeChannel,TrialNum);
annotation('textbox', [0.17, 0.95, 0.7, 0.04], 'string', data_msg,'FontSize',20,'FontWeight','Bold','LineStyle','none')
end
%

%Output to command line

%Reorganize data
formatstr=[];
Output=[];
for i=1:length(SegEvent)-1
    formatstr=strcat(formatstr,'%f\t');
    Output=[ Output,Sig_Event(i),p_Event(i)];

end

formatstr=strcat(formatstr,'%f\t');
Output=[Output,nanmean(Spon_FR)];

%fprintf(formatstr,Output)

OutputToCommand=sprintf(formatstr,Output)

%Sig_Event,p_Event


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
    if isfield(OutputData,'EventSegregation')
    ExistFlag=1;
    end
end
%OutputData For current analysis
OutputData_New=[];

OutputData_New.EventSegregation.TrialType='GoodTrials';

%%%%%%%%%Save overall scene response characteristics%%%%%%%%%%%%%%%%%%%%%%%

%Organize the data
NameStr={'TrialNum','PSTH_Event','PSTH_Event_Time','EventTimeMarker','AverageFiringRate_Event','AverageFiringRate_Event_Mean','AverageFiringRate_Event_Sem',...
    'Spon_FR','p_Event','SegEvent','SegEventStr','SegEventEnd','PlotStartOffset','PlotStopOffset','StepSize','BinWidth'};

DataLib={TrialNum,PSTH_Event,PSTH_Event_Time,EventTimeMarker,AverageFiringRate_Event,AverageFiringRate_Event_Mean,AverageFiringRate_Event_Sem,...
Spon_FR,  {Sig_Event',p_Event'}  ,SegEvent,SegEventStr,SegEventEnd,PlotStartOffset,PlotStopOffset,StepSize,BinWidth

};


DataStamp=containers.Map(NameStr,DataLib);
OutputData_New.EventSegregation.Task=TaskType;
OutputData_New.EventSegregation.TaskCode=TaskCode;
OutputData_New.EventSegregation.DataStamp=DataStamp;

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
    saveas(fig103,OutputFigureName1);
    OutputFigureName2=strcat(OutputFigureName,'AverageEventResponse.jpg');
    saveas(fig104,OutputFigureName2);
   disp('Figures have been exported to the neuron folder');

end
 
if ExistFlag
%Compare the old one with the new one
Task_Old=OutputData.EventSegregation.TaskCode;

     if Task_Old~=TaskCode
        %Not the same task, add the new dataset into the old one
        OutputData.EventSegregation(numel(OutputData.EventSegregation)+1)=OutputData_New.EventSegregation;
     else
         %If the same task,the same protocol,replace the old one with the new one
         
         OutputData.EventSegregation=OutputData_New.EventSegregation;
         
    
    
     end
else
    %Output the current file
   OutputData.EventSegregation=OutputData_New.EventSegregation; 
end



if OutputFlag
    cd(OutPath);
    save(OutputFileName,'OutputData');
    disp('Data Exported');
else
    disp('No data export');
end  




end