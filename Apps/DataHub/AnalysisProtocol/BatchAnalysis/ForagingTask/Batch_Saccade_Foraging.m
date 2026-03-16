function  Batch_Saccade_Foraging(Data,  StartFile, EndFile,ShowFigureFlag,OutputFlag);

%Xuefei Yu 06182020
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
     Task(i)=OutputData.TargetSaccade(1).TaskCode;
     TaskName{i}=OutputData.TargetSaccade(1).Task;
     
     DataLib{i}=OutputData.TargetSaccade(1).DataStamp;
     
     
     
end

%Retrieve the PSTH for the first saccade 
PSTH_SacOn_Mean_Time=RetriveDataLib('PSTH_SacOn_Mean_Time',DataLib,0);
PSTH_SacOn_FR_Good=RetriveDataLib('PSTH_SacOn_FR_Good',DataLib,1);
PSTH_SacOn_FR_Bad=RetriveDataLib('PSTH_SacOn_FR_Bad',DataLib,1);

PSTH_SacOn_Mean_Time=nanmean(PSTH_SacOn_Mean_Time,1);
PSTH_SacOn_FR_Good_mean=nanmean(PSTH_SacOn_FR_Good,1);
PSTH_SacOn_FR_Bad_mean=nanmean(PSTH_SacOn_FR_Bad,1);

PSTH_SacOn_FR_Good_sem=nanstd(PSTH_SacOn_FR_Good,[],1)./sqrt(sum(~isnan(PSTH_SacOn_FR_Good),1));
PSTH_SacOn_FR_Bad_sem=nanstd(PSTH_SacOn_FR_Bad,[],1)./sqrt(sum(~isnan(PSTH_SacOn_FR_Bad),1));



%Retrieve the PSTH for the second saccade 

PSTH_SacOnReturn_Mean_Time=RetriveDataLib('PSTH_SacOnReturn_Mean_Time',DataLib,0);
PSTH_SacOnReturn_FR_Good=RetriveDataLib('PSTH_SacOnReturn_FR_Good',DataLib,1);
PSTH_SacOnReturn_FR_Bad=RetriveDataLib('PSTH_SacOnReturn_FR_Bad',DataLib,1);

PSTH_SacOnReturn_Mean_Time=nanmean(PSTH_SacOnReturn_Mean_Time,1);
PSTH_SacOnReturn_FR_Good_mean=nanmean(PSTH_SacOnReturn_FR_Good,1);
PSTH_SacOnReturn_FR_Bad_mean=nanmean(PSTH_SacOnReturn_FR_Bad,1);

PSTH_SacOnReturn_FR_Good_sem=nanstd(PSTH_SacOnReturn_FR_Good,[],1)./sqrt(sum(~isnan(PSTH_SacOnReturn_FR_Good),1));
PSTH_SacOnReturn_FR_Bad_sem=nanstd(PSTH_SacOnReturn_FR_Bad,[],1)./sqrt(sum(~isnan(PSTH_SacOnReturn_FR_Bad),1));



if ShowFigureFlag

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Plot the
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%average
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%PSTH for
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%all saccades
fig1=PrepareFigure(1,'w',[100,200, 1000,800],'Name','PSTH of Saccades towards and away from targets');%,'Name','SceneTuning');
subplot(1,2,1); %Saccades towards the target
%Good
shadedErrorBar(PSTH_SacOn_Mean_Time,PSTH_SacOn_FR_Good_mean,PSTH_SacOn_FR_Good_sem,'lineprops',{[1,0,0]},'transparent',1,'patchSaturation',0.5);
hold on
%Bad
shadedErrorBar(PSTH_SacOn_Mean_Time,PSTH_SacOn_FR_Bad_mean,PSTH_SacOn_FR_Bad_sem,'lineprops',{[0,0,1]},'transparent',1,'patchSaturation',0.4);

%Labels
yrange=ylim;
h1=plot([0,0],[yrange(1),yrange(2)],'-.k');


set(get(get(h1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');



xlim([min(PSTH_SacOn_Mean_Time),max(PSTH_SacOn_Mean_Time)-20]);
xlabel('Time from Sac Onset to Target (ms)');
ylabel('Fring Rate(Hz)')
title('PSTH for all neurons aligned on Sac onset')

set(gca,'LineWidth',3)
set(findall(gca, 'Type', 'Line'),'LineWidth',3);
set(findall(gcf,'-property','FontSize'),'FontSize',20);
set(findall(gcf,'-property','FontWeight'),'FontWeight','Bold');

subplot(1,2,2); %Saccades towards the target
%Good
shadedErrorBar(PSTH_SacOnReturn_Mean_Time,PSTH_SacOnReturn_FR_Good_mean,PSTH_SacOnReturn_FR_Good_sem,'lineprops',{[1,0,0]},'transparent',1,'patchSaturation',0.5);
hold on
%Bad
shadedErrorBar(PSTH_SacOnReturn_Mean_Time,PSTH_SacOnReturn_FR_Bad_mean,PSTH_SacOnReturn_FR_Bad_sem,'lineprops',{[0,0,1]},'transparent',1,'patchSaturation',0.4);

%Labels
yrange=ylim;
h1=plot([0,0],[yrange(1),yrange(2)],'-.k');

%Legend


g=legend({'Good','Bad'},'color','none','Box', 'off');%,'Position',[0.38    0.75    0.4510    0.1225]);

set(get(get(h1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');



xlim([min(PSTH_SacOn_Mean_Time),max(PSTH_SacOn_Mean_Time)-20]);
xlabel('Time from Return Sac from Target (ms)');
ylabel('Fring Rate(Hz)')
title('PSTH for all neurons aligned on Return Saccade onset')

set(gca,'LineWidth',3)
set(findall(gca, 'Type', 'Line'),'LineWidth',3);
set(findall(gcf,'-property','FontSize'),'FontSize',20);
set(findall(gcf,'-property','FontWeight'),'FontWeight','Bold');




end

end