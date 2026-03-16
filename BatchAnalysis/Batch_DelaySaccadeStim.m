function Batch_DelaySaccadeStim(Data,  StartFile, EndFile,ShowFigureFlag,OutputFlag);

BatchFileName=Data.BatchFileName;
%Batch data file path
FilesName=Data.ResultFilePath(StartFile: EndFile);
RecordDate=Data.RecordDate(StartFile: EndFile);
FileNum=EndFile-StartFile+1;
FilesName_Each=Data.FileName;


for i=1:length(FilesName)
  clear OutputData;
  
     load(FilesName{i}{1});
     
     OutputDataTemp=OutputData;
     
     
     Task(i)=OutputData.SaccadeStim(1).TaskCode;
     TaskName{i}=OutputData.SaccadeStim(1).Task;
     
     
     Data =OutputData.SaccadeStim(1).DataStamp;
    


      Monkey = FilesName_Each{i}(1);

          if Monkey == 'R'
              MonkeyIndex(i) = 1;
          else
              MonkeyIndex(i) = 2;
          end

%'RT_stim_RF_period_Mean','RT_control_RF_period_Mean','RT_stim_NonRF_period_Mean','RT_control_NonRF_period_Mean','p_RF','p_NonRF'

RT_stimRF_Mean(i) = Data('RT_stim_RF_period_Mean');
RT_controlRF_Mean(i) = Data('RT_control_RF_period_Mean');

RT_stim_NonRF_Mean(i) = Data('RT_stim_NonRF_period_Mean');
RT_control_NonRF_Mean(i) = Data('RT_control_NonRF_period_Mean');

p_RF_Each(i) = Data('p_RF')/2;
p_NonRF_Each(i) = Data('p_NonRF')/2;

p_RF_Each_sr(i) = Data('p_RF_sr')/2;
p_NonRF_Each_sr(i) = Data('p_NonRF_sr')/2;



     


     %%
 
end

[h,p_RF]=ttest(RT_stimRF_Mean,RT_controlRF_Mean,'tail','left');
[h,p_nonRF]=ttest(RT_stim_NonRF_Mean,RT_control_NonRF_Mean,'tail','right');


[p_RF_sr,h,stats]=signrank(RT_stimRF_Mean,RT_controlRF_Mean,'tail','left');
[p_nonRF_sr,h,stats]= signrank(RT_stim_NonRF_Mean,RT_control_NonRF_Mean,'tail','right');

RT_stimRF_Mean_Mean = mean(RT_stimRF_Mean);
RT_stimRF_Mean_Sem = std(RT_stimRF_Mean)/sqrt(length(RT_stimRF_Mean));

RT_controlRF_Mean_Mean = mean(RT_controlRF_Mean);
RT_controlRF_Mean_Sem = std(RT_controlRF_Mean)/sqrt(length(RT_controlRF_Mean));

RT_stimNonRF_Mean_Mean = mean(RT_stim_NonRF_Mean);
RT_stimNonRF_Mean_Sem = std(RT_stim_NonRF_Mean)/sqrt(length(RT_stim_NonRF_Mean));

RT_controlNonRF_Mean_Mean = mean(RT_control_NonRF_Mean);
RT_controlNonRF_Mean_Sem = std(RT_control_NonRF_Mean)/sqrt(length(RT_control_NonRF_Mean));

disp('RF: StimMean');
disp(RT_stimRF_Mean_Mean)

disp('RF: StimSem');
disp(RT_stimRF_Mean_Sem)

disp('RF: ControlkStimMean');
disp(RT_controlRF_Mean_Mean)

disp('RF: ControlkStimSem');
disp(RT_controlRF_Mean_Sem)

disp('RF: StimSem');
disp(RT_stimRF_Mean_Sem)

disp('NonRF: StimMean');
disp(RT_stimNonRF_Mean_Mean)

disp('NonRF: StimSem');
disp(RT_stimNonRF_Mean_Sem)

disp('NonRF: ControlMean');
disp(RT_controlNonRF_Mean_Mean)

disp('NonRF: ControlSem');
disp(RT_controlNonRF_Mean_Sem)

disp('p_RF')
disp(p_RF_sr)

disp('p_nonRF')
disp(p_nonRF_sr)

N=length(RT_stimRF_Mean);
disp('N=');
disp(N)




disp('RF: Stim-Control:');
disp(RT_stimRF_Mean_Mean-RT_controlRF_Mean_Mean);

disp('NonRF: Stim-Control:');
disp(RT_stimNonRF_Mean_Mean-RT_controlNonRF_Mean_Mean);


if ShowFigureFlag
   
%% Scatter plot
FigureStartNum=90;
FigureIndex=1;
    
figtitlestr{FigureIndex}='ScatterPlot_RF_NonRF';
fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[50,100, 1000,500],'Name',figtitlestr{FigureIndex});
set(fig{FigureIndex}, 'PaperUnits', 'inches');
formatFig(fig{FigureIndex}, [4.5 2], 'nature');


subplot(1,2,1)
scatter(RT_stimRF_Mean(MonkeyIndex==1),RT_controlRF_Mean(MonkeyIndex==1),'filled', 'Marker', '^','MarkerFaceColor',[1,0,0]);
hold on
scatter(RT_stimRF_Mean(MonkeyIndex==2),RT_controlRF_Mean(MonkeyIndex==2),'filled', 'Marker', 's','MarkerFaceColor',[1,0,0]);
%max_v = max([max(RT_stimRF_Mean),max(RT_controlRF_Mean)])+10;
%min_v = min([min(RT_stimRF_Mean),min(RT_controlRF_Mean)])-10;
plot([90:300],[90:300],'--k');
%xlim([min_v,max_v]);
%ylim([min_v,max_v]);
%xticks([floor(min_v),round((floor(min_v)+ceil(max_v))/2),ceil(max_v)]);
%yticks([floor(min_v),round((floor(min_v)+ceil(max_v))/2),ceil(max_v)]);
xlim([90,300]);
ylim([90,300]);

legend({'Monkey R','Monkey A'});

%axis equal;
xlabel('RT(Opto-Stim)');
ylabel('RT(Control)');
title('Saccades toward RF');
set(gca,'FontSize',5,'LineWidth',0.5);


subplot(1,2,2)
scatter(RT_stim_NonRF_Mean(MonkeyIndex==1),RT_control_NonRF_Mean(MonkeyIndex==1),'filled', 'Marker', '^','MarkerFaceColor',[1,0,0]);

hold on
scatter(RT_stim_NonRF_Mean(MonkeyIndex==2),RT_control_NonRF_Mean(MonkeyIndex==2),'filled', 'Marker', 's','MarkerFaceColor',[1,0,0]);
max_v = max([max(RT_stim_NonRF_Mean),max(RT_control_NonRF_Mean)])+10;
min_v = min([min(RT_stim_NonRF_Mean),min(RT_control_NonRF_Mean)])-10;
plot([90:300],[90:300],'--k');
%xlim([min_v,max_v]);
%ylim([min_v,max_v]);
%xticks([floor(min_v),round((floor(min_v)+ceil(max_v))/2),ceil(max_v)]);
%yticks([floor(min_v),round((floor(min_v)+ceil(max_v))/2),ceil(max_v)]);
xlim([90,300]);
ylim([90,300]);
%axis equal;
xlabel('RT(Opto-Stim)');
ylabel('RT(Control)');
title('Saccades away from RF');
set(gca,'FontSize',5,'LineWidth',0.5);

%% Scatter plot with significance for each case
%p_RF_Each
FigureStartNum=FigureStartNum+1;
FigureIndex=FigureIndex+1;

figtitlestr{FigureIndex}='ScatterPlot_RF_NonRF_EachSig';
fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[50,100, 1000,500],'Name',figtitlestr{FigureIndex});
set(fig{FigureIndex}, 'PaperUnits', 'inches');
formatFig(fig{FigureIndex}, [4.5 2], 'nature');


subplot(1,2,1)
Sig_RF = p_RF_Each_sr<0.05;
Sig_NonRF = p_NonRF_Each_sr<0.05;
plot([90:300],[90:300],'--k');
hold on
scatter(RT_stimRF_Mean(MonkeyIndex==1 & ~Sig_RF),RT_controlRF_Mean(MonkeyIndex==1 & ~Sig_RF), 'Marker', '^','MarkerFaceColor','w','MarkerEdgeColor','k');

scatter(RT_stimRF_Mean(MonkeyIndex==2 & ~Sig_RF),RT_controlRF_Mean(MonkeyIndex==2 & ~Sig_RF),'Marker', 's','MarkerFaceColor','w','MarkerEdgeColor','k');

scatter(RT_stimRF_Mean(MonkeyIndex==1 & Sig_RF),RT_controlRF_Mean(MonkeyIndex==1 & Sig_RF),'filled', 'Marker', '^','MarkerFaceColor',[1,0,0],'MarkerEdgeColor','k');

scatter(RT_stimRF_Mean(MonkeyIndex==2 & Sig_RF),RT_controlRF_Mean(MonkeyIndex==2 & Sig_RF),'filled', 'Marker', 's','MarkerFaceColor',[1,0,0],'MarkerEdgeColor','k');
%max_v = max([max(RT_stimRF_Mean),max(RT_controlRF_Mean)])+10;
%min_v = min([min(RT_stimRF_Mean),min(RT_controlRF_Mean)])-10;

%xlim([min_v,max_v]);
%ylim([min_v,max_v]);
%xticks([floor(min_v),round((floor(min_v)+ceil(max_v))/2),ceil(max_v)]);
%yticks([floor(min_v),round((floor(min_v)+ceil(max_v))/2),ceil(max_v)]);
xlim([90,300]);
ylim([90,300]);
yticks([100,200,300]);
xticks([100,200,300]);

legend({'Monkey R','Monkey A'});

%axis equal;
xlabel('RT(Opto-Stim)');
ylabel('RT(Control)');
title('Saccades toward RF');
set(gca,'FontSize',5,'LineWidth',0.5);
box off

subplot(1,2,2)
plot([90:300],[90:300],'--k');
hold on

scatter(RT_stim_NonRF_Mean(MonkeyIndex==1 & ~Sig_NonRF),RT_control_NonRF_Mean(MonkeyIndex==1 & ~Sig_NonRF), 'Marker', '^','MarkerFaceColor','w','MarkerEdgeColor','k');

scatter(RT_stim_NonRF_Mean(MonkeyIndex==2 & ~Sig_NonRF),RT_control_NonRF_Mean(MonkeyIndex==2 & ~Sig_NonRF),'Marker', 's','MarkerFaceColor','w','MarkerEdgeColor','k');

scatter(RT_stim_NonRF_Mean(MonkeyIndex==1 & Sig_NonRF),RT_control_NonRF_Mean(MonkeyIndex==1 & Sig_NonRF),'filled', 'Marker', '^','MarkerFaceColor',[1,0,0],'MarkerEdgeColor','k');

scatter(RT_stim_NonRF_Mean(MonkeyIndex==2 & Sig_NonRF),RT_control_NonRF_Mean(MonkeyIndex==2 & Sig_NonRF),'filled', 'Marker', 's','MarkerFaceColor',[1,0,0],'MarkerEdgeColor','k');

max_v = max([max(RT_stim_NonRF_Mean),max(RT_control_NonRF_Mean)])+10;
min_v = min([min(RT_stim_NonRF_Mean),min(RT_control_NonRF_Mean)])-10;

%xlim([min_v,max_v]);
%ylim([min_v,max_v]);
%xticks([floor(min_v),round((floor(min_v)+ceil(max_v))/2),ceil(max_v)]);
%yticks([floor(min_v),round((floor(min_v)+ceil(max_v))/2),ceil(max_v)]);
xlim([90,300]);
ylim([90,300]);
yticks([100,200,300]);
xticks([100,200,300]);
box off
%axis equal;
xlabel('RT(Opto-Stim)');
ylabel('RT(Control)');
title('Saccades away from RF');
set(gca,'FontSize',5,'LineWidth',0.5);

%% Bar plot
%% RF
FigureStartNum=FigureStartNum+1;
FigureIndex=FigureIndex+1;

figtitlestr{FigureIndex}='BarPlot_RF';
fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[70,100, 400,400],'Name',figtitlestr{FigureIndex});
set(fig{FigureIndex}, 'PaperUnits', 'inches');
formatFig(fig{FigureIndex}, [0.8 0.8], 'nature');


b = bar([1,2],[RT_stimRF_Mean_Mean,RT_controlRF_Mean_Mean]);
b.FaceColor = 'flat';
b.CData = [1,0,0;0.3,0.3,0.3];
xticklabels({'Opto-Stim','Control'});
hold on
errorbar([1,2],[RT_stimRF_Mean_Mean,RT_controlRF_Mean_Mean],[RT_stimRF_Mean_Sem,RT_controlRF_Mean_Sem],'.k');
box off
ylim([100,200]);
yticks([100,200]);
set(gca,'FontSize',5,'LineWidth',0.5);

%% NonRF
FigureStartNum=FigureStartNum+1;
FigureIndex=FigureIndex+1;

figtitlestr{FigureIndex}='BarPlot_NonRF';
fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[90,100, 400,400],'Name',figtitlestr{FigureIndex});
set(fig{FigureIndex}, 'PaperUnits', 'inches');
formatFig(fig{FigureIndex}, [0.8 0.8], 'nature');

b = bar([1,2],[RT_stimNonRF_Mean_Mean,RT_controlNonRF_Mean_Mean]);
b.FaceColor = 'flat';
b.CData = [1,0,0;0.3,0.3,0.3];
xticklabels({'Opto-Stim','Control'});
hold on
errorbar([1,2],[RT_stimNonRF_Mean_Mean,RT_controlNonRF_Mean_Mean],[RT_stimNonRF_Mean_Sem,RT_controlNonRF_Mean_Sem],'.k');
box off

ylim([100,200]);
yticks([100,200]);
set(gca,'FontSize',5,'LineWidth',0.5);




end

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
OutPath=strcat(BasicDirectory,'Results',BasicDirectory(end),'ExportFigure',BasicDirectory(end));%,OutputFolerName);
if ~exist(OutPath)
    mkdir(OutPath)
end
cd(OutPath);

%OutputFigureName='SC_VideoVeiw_Behavior';
OutputFigureName='FEF_DelayStim';

if OutputFlag == 1

for jj=1:FigureIndex
        if ~isempty(fig{jj})
        

         OutputFigureName_Final=strcat(OutputFigureName,figtitlestr{jj},'.pdf');
         saveas(fig{jj},OutputFigureName_Final);
        end
        
            
end

    disp('Figures have been exported to the exported folder');

end


end