function   Batch_FractalTuning_Foraging(Data,  StartFile, EndFile,ShowFigureFlag,OutputFlag);
%Fractal Tuning for foraging tasks
%Xuefei Yu 06072020
BatchFileName=Data.BatchFileName;
%Batch data file path
FilesName=Data.ResultFilePath(StartFile: EndFile);
RecordDate=Data.RecordDate(StartFile: EndFile);
FileNum=EndFile-StartFile+1;

RawFileName=Data.FileName;
%Color Definition for each events
Color=[lines(6);prism(6);pink(6)];
%Color(2,:)=copper(1)+0.5;
Color=[Color;hsv(10)];
Color(18,:)=[0.1,0.457,0.194];

%Load file
for i=1:length(FilesName)
  clear OutputData;
  
     load(FilesName{i});
     
     OutputDataTemp=OutputData;
     Task(i)=OutputData.FractalTuning(1).TaskCode;
     TaskName{i}=OutputData.FractalTuning(1).Task;
     
     DataLib{i}=OutputData.FractalTuning(1).DataStamp;
     EventSeg{i}=OutputData.EventSegregation(1).DataStamp;
     
end

%The overall fractal response by event segmentation
p_FractalFromSpon_Whole_temp=cellfun(@(x) x('p_Event'),EventSeg,'Uniform',0);

h_FractalFromSpon_Whole=cell2mat(cellfun(@(x) x{1}(2),p_FractalFromSpon_Whole_temp,'Uniform',0));
p_FractalFromSpon_Whole=cell2mat(cellfun(@(x) x{2}(2),p_FractalFromSpon_Whole_temp,'Uniform',0));

%Find out spon for each trial from the event segration 
AverageFR_Spon_EachTrial=cellfun(@(x) nanmean(x('Spon_FR')),EventSeg,'Uniform',0);%Abosolute spon

PSTH_Event_EachTrial=cellfun(@(x) x('PSTH_Event'),EventSeg,'Uniform',0);
PSTH_EventTime_EachTrial=cellfun(@(x) x('PSTH_Event_Time'),EventSeg,'Uniform',0);

PSTH_EventTime_Marker=cellfun(@(x) x('RelativeTimeMarker'), DataLib,'Uniform',0)';
PSTH_EventTime_Marker_mean=cell2mat(cellfun(@(x) min(x,[],1), PSTH_EventTime_Marker,'Uniform',0));

%PreEvent spon
AverageFR_Spon_PreFrac=cell2mat(cellfun(@(x) nanmean(x('AverageFR_Spon')),DataLib,'Uniform',0))';%Relative spon




PSTH_Spon_EachTrial=cell2mat(cellfun(@(x) nanmean(x{end},1)',PSTH_Event_EachTrial,'Uniform',0))';
PSTH_SponTime_EachTrial_mean=nanmean(cell2mat(cellfun(@(x) nanmean(x{end})',PSTH_EventTime_EachTrial,'Uniform',0))',1);

PSTH_Spon_EachTrial_mean=nanmean(PSTH_Spon_EachTrial,1);
PSTH_Spon_EachTrial_sem=nanstd(PSTH_Spon_EachTrial,[],1)./sqrt(sum(~isnan(PSTH_Spon_EachTrial),1));
PSTH_SponTime_EachTrial_mean=PSTH_SponTime_EachTrial_mean-1000;


%Find out fractal responsive cell
p_FractalFromSpon_PreFrac_temp=cellfun(@(x) x('p_AverageFR_Fractal'),DataLib,'Uniform',0);
h_FractalFromSpon_PreFrac=cell2mat(cellfun(@(x) x(1),p_FractalFromSpon_PreFrac_temp,'Uniform',0));
p_FractalFromSpon_PreFrac=cell2mat(cellfun(@(x) x(2),p_FractalFromSpon_PreFrac_temp,'Uniform',0));





%Use criteria 1
%FractalResponseNeuron=h_FractalFromSpon_Whole==1;

%Use criteria 2
FractalResponseNeuron=h_FractalFromSpon_PreFrac==1;
%Use criteria 3
%FractalResponseNeuron=h_FractalFromSpon_PreFrac==1 | h_FractalFromSpon_Whole==1;

%Retrieve the PSTH for the scene whole

PSTH_Fractal_Time=cellfun(@(x) x('PSTH_Time_All_mean'),DataLib,'Uniform',0);
PSTH_Fractal_FR=cellfun(@(x) x('PSTH_FractalAverage'),DataLib,'Uniform',0);
PSTH_Fractal_FR=cellfun(@(x) x(:,1)',PSTH_Fractal_FR,'Uniform',0);%Only need the average 


Max_Col_Whole=max(cell2mat(cellfun(@numel ,PSTH_Fractal_Time,'Uniform',0)));

PSTH_Fractal_Time=cell2mat(cellfun(@(x) [x,NaN*ones(1,(Max_Col_Whole-numel(x)))]',PSTH_Fractal_Time,'Uniform',0))';
PSTH_Fractal_FR=cell2mat(cellfun(@(x) [x,NaN*ones(1,(Max_Col_Whole-numel(x)))]',PSTH_Fractal_FR,'Uniform',0))';

PSTH_Fractal_Time_mean=nanmean(PSTH_Fractal_Time,1);

%Reorder PSTH according to the pvalue
[~,FractalResponseOrder]=sort(p_FractalFromSpon_PreFrac);
PSTH_Fractal_FR_Sort=PSTH_Fractal_FR(FractalResponseOrder,:);
RawFileName_Order=RawFileName(FractalResponseOrder);
%Normalize
%PSTH_Fractal_FR_Sort=PSTH_Fractal_FR_Sort-AverageFR_Spon_PreFrac;
%PSTH_Fractal_FR_Sort=(PSTH_Fractal_FR_Sort-nanmean(PSTH_Fractal_FR_Sort,2))./nanstd(PSTH_Fractal_FR_Sort,[],2);


%For fractal responsive neurons 
PSTH_Fractal_responsive_FR=PSTH_Fractal_FR(FractalResponseNeuron,:);
PSTH_Fractal_FR_mean_Fractal=nanmean(PSTH_Fractal_FR(FractalResponseNeuron,:));
PSTH_Fractal_FR_sem_Fractal=nanstd(PSTH_Fractal_FR(FractalResponseNeuron,:))./sqrt(sum(~isnan(PSTH_Fractal_FR(FractalResponseNeuron,:)),1));

PSTH_EventTime_Marker_mean_Fractal=PSTH_EventTime_Marker_mean(FractalResponseNeuron,:);

%For fractal irresponsive neurons 
PSTH_Fractal_FR_mean_NonFractal=nanmean(PSTH_Fractal_FR(~FractalResponseNeuron,:));
PSTH_Fractal_FR_sem_NonFractal=nanstd(PSTH_Fractal_FR(~FractalResponseNeuron,:))./sqrt(sum(~isnan(PSTH_Fractal_FR(~FractalResponseNeuron,:)),1));

PSTH_EventTime_Marker_mean_NonFractal=PSTH_EventTime_Marker_mean(~FractalResponseNeuron,:);

PSTH_EventTime_Marker_mean_All=mean(PSTH_EventTime_Marker_mean,1);


Num_FractalNeuron=sum(FractalResponseNeuron);
Num_NonFractalNeuron=sum(~FractalResponseNeuron);

%Overlay the eye trace
PSTH_EyeEcc=cellfun(@(x) x('Eye_Ecc_PSTH_All_mean'),DataLib,'Uniform',0);
PSTH_EyeEcc=cellfun(@(x) x(:,1)',PSTH_EyeEcc,'Uniform',0)';%Only need the mean

max_numel=max(cellfun(@numel,PSTH_EyeEcc));
PSTH_EyeEcc_Mat=cell2mat(cellfun(@(x) [x, NaN*ones(1,max_numel-length(x))],PSTH_EyeEcc,'Uniform',0));

PSTH_EyeTimeOffset=cell2mat(cellfun(@(x) x('FractalSegTimeOffset'),DataLib,'Uniform',0)');

PSTH_EyeEcc_Fractal=PSTH_EyeEcc_Mat(FractalResponseNeuron,:);
PSTH_EyeEcc_mean_Fractal=nanmean(PSTH_EyeEcc_Fractal,1);
PSTH_EyeEcc_sem_Fractal=nanstd(PSTH_EyeEcc_Fractal,[],1)./sqrt(sum(~isnan(PSTH_EyeEcc_Fractal),1));

PSTH_EyeEcc_NonFractal=PSTH_EyeEcc_Mat(~FractalResponseNeuron,:);
PSTH_EyeEcc_mean_NonFractal=nanmean(PSTH_EyeEcc_NonFractal,1);
PSTH_EyeEcc_sem_NonFractal=nanstd(PSTH_EyeEcc_NonFractal,[],1)./sqrt(sum(~isnan(PSTH_EyeEcc_NonFractal),1));



if ShowFigureFlag

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Plot the average PSTH for each type of neurons
fig1=PrepareFigure(10,'w',[100,200, 1000,800],'Name','PSTH of Fractal response neurons and NonFractal response neurons');%,'Name','SceneTuning');
subplot(2,1,1)
%Plot spon
shadedErrorBar(PSTH_SponTime_EachTrial_mean,PSTH_Spon_EachTrial_mean,PSTH_Spon_EachTrial_sem,'lineprops',{'#000000'},'transparent',1,'patchSaturation',0.2);
%Fractal Response
hold on
shadedErrorBar(PSTH_Fractal_Time_mean,PSTH_Fractal_FR_mean_Fractal,PSTH_Fractal_FR_sem_Fractal,'lineprops',{'#DC4826'},'transparent',1,'patchSaturation',0.5);
%Non Fractal response
shadedErrorBar(PSTH_Fractal_Time_mean,PSTH_Fractal_FR_mean_NonFractal,PSTH_Fractal_FR_sem_NonFractal,'lineprops',{'#648FFF'},'transparent',1,'patchSaturation',0.4);
xlim([min(PSTH_EventTime_Marker_mean_All(1))-150,min(PSTH_EventTime_Marker_mean_All(2))+200]);
%Labels
yrange=ylim;
h1=plot([PSTH_EventTime_Marker_mean_All(1),PSTH_EventTime_Marker_mean_All(1)],[yrange(1),yrange(2)],'-.k');
h2=plot([PSTH_EventTime_Marker_mean_All(2),PSTH_EventTime_Marker_mean_All(2)],[yrange(1),yrange(2)],'--k');

%Legend
Num_FractalNeuron_str=sprintf('Fractal Responsive\nneurons(N=%d)',Num_FractalNeuron);
Num_NonFractalNeuron_str=sprintf('NonFractal Responsive\nneurons(N=%d)',Num_NonFractalNeuron);

g=legend({'Spon',Num_FractalNeuron_str,Num_NonFractalNeuron_str},'color','none','Box', 'off','Position',[0.38    0.75    0.4510    0.1225]);

set(get(get(h1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(h2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');




xlabel('Time from Fractal onset (ms)');
ylabel('Fring Rate(Hz)')
title('PSTH for all neurons aligned on Fractal onset')

set(gca,'LineWidth',3)
set(findall(gca, 'Type', 'Line'),'LineWidth',3);
set(findall(gcf,'-property','FontSize'),'FontSize',20);
set(findall(gcf,'-property','FontWeight'),'FontWeight','Bold');

subplot(2,1,2) %Eye trace
%Fractal Response

shadedErrorBar([-150:length(PSTH_EyeEcc_mean_Fractal)-150-1],PSTH_EyeEcc_mean_Fractal,PSTH_EyeEcc_sem_Fractal,'lineprops',{'#DC4826'},'transparent',1,'patchSaturation',0.5);
hold on
%Non Fractal response
shadedErrorBar([-150:length(PSTH_EyeEcc_mean_NonFractal)-150-1],PSTH_EyeEcc_mean_NonFractal,PSTH_EyeEcc_sem_NonFractal,'lineprops',{'#648FFF'},'transparent',1,'patchSaturation',0.4);

%Labels
yrange=ylim;
h1=plot([PSTH_EventTime_Marker_mean_All(1),PSTH_EventTime_Marker_mean_All(1)],[yrange(1),yrange(2)],'-.k');
h2=plot([PSTH_EventTime_Marker_mean_All(2),PSTH_EventTime_Marker_mean_All(2)],[yrange(1),yrange(2)],'--k');

xlim([(PSTH_EventTime_Marker_mean_All(1))-150,(PSTH_EventTime_Marker_mean_All(2))+200]);
xlabel('Time from Fractal onset (ms)');
ylabel('Eye Eccentricity')
title('PSTH for all eyemovements aligned on Fractal onset')

set(gca,'LineWidth',3)
set(findall(gca, 'Type', 'Line'),'LineWidth',3);
set(findall(gcf,'-property','FontSize'),'FontSize',20);
set(findall(gcf,'-property','FontWeight'),'FontWeight','Bold');

%
%%Check the PSTH for each cel;

fig2=PrepareFigure(11,'w',[100,200, 1500,1000],'Name','PSTH of Each neuron');%,'Name','SceneTuning');

for i=1:size(PSTH_Fractal_FR_Sort,1)
    subplot(10,6,i)
plot(PSTH_Fractal_Time_mean,PSTH_Fractal_FR_Sort(i,:),'-r','LineWidth',3);
xlim([min(PSTH_EventTime_Marker_mean_All(1))-150,min(PSTH_EventTime_Marker_mean_All(2))+200]);
title(RawFileName_Order{i});
end
%}


end


end