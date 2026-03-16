function Batch_SceneFractalCorrelation_Foraging(Data,  StartFile, EndFile,ShowFigureFlag,OutputFlag);
%Function to check the response relationship between the scene and the
%fractal onset on population level
%Xuefei Yu 06162020
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
     
     DataLibFrac{i}=OutputData.FractalTuning(1).DataStamp;
     
     DataLibScene{i}=OutputData.SceneTuning(1).DataStamp;
     
     EventSeg{i}=OutputData.EventSegregation(1).DataStamp;
     
end


%Retrieve the PSTH for the scene whole

PSTH_Scene_Whole=cellfun(@(x) x('PSTH_EachTrial'),DataLibScene,'Uniform',0);
PSTH_Scene_Time=cellfun(@(x) nanmean(x{1}),PSTH_Scene_Whole,'Uniform',0);
PSTH_Scene_FR=cellfun(@(x) nanmean(x{2}),PSTH_Scene_Whole,'Uniform',0);

Max_Col_Whole=max(cell2mat(cellfun(@numel ,PSTH_Scene_Time,'Uniform',0)));

PSTH_Scene_Time=cell2mat(cellfun(@(x) [x,NaN*ones(1,(Max_Col_Whole-numel(x)))]',PSTH_Scene_Time,'Uniform',0))';
PSTH_Scene_FR=cell2mat(cellfun(@(x) [x,NaN*ones(1,(Max_Col_Whole-numel(x)))]',PSTH_Scene_FR,'Uniform',0))';

PSTH_Scene_Time_mean=nanmean(PSTH_Scene_Time);

%Retrieve the PSTH for the scene whole

PSTH_Fractal_Time=cellfun(@(x) x('PSTH_Time_All_mean'),DataLibFrac,'Uniform',0);
PSTH_Fractal_FR=cellfun(@(x) x('PSTH_FractalAverage'),DataLibFrac,'Uniform',0);
PSTH_Fractal_FR=cellfun(@(x) x(:,1)',PSTH_Fractal_FR,'Uniform',0);%Only need the average 

Max_Col_Whole=max(cell2mat(cellfun(@numel ,PSTH_Fractal_Time,'Uniform',0)));

PSTH_Fractal_Time=cell2mat(cellfun(@(x) [x,NaN*ones(1,(Max_Col_Whole-numel(x)))]',PSTH_Fractal_Time,'Uniform',0))';
PSTH_Fractal_FR=cell2mat(cellfun(@(x) [x,NaN*ones(1,(Max_Col_Whole-numel(x)))]',PSTH_Fractal_FR,'Uniform',0))';

PSTH_Fractal_Time_mean=nanmean(PSTH_Fractal_Time,1);

%Get the response proportion and examine the cell type
ResponseProportion_Scene_EachNeuron=cell2mat(cellfun(@(x) x('ResponseProportion'),DataLibScene,'Uniform',0));

Transcient_Neurons_Scene=ResponseProportion_Scene_EachNeuron<=30 &ResponseProportion_Scene_EachNeuron>0;
Sustained_Neurons_Scene=ResponseProportion_Scene_EachNeuron>30;


PSTH_Fractal_FR_Transcience=PSTH_Fractal_FR(Transcient_Neurons_Scene,:);
PSTH_Fractal_FR_Sustained=PSTH_Fractal_FR(Sustained_Neurons_Scene,:);

PSTH_Scene_FR_Transcience=PSTH_Scene_FR(Transcient_Neurons_Scene,:);
PSTH_Scene_FR_Sustained=PSTH_Scene_FR(Sustained_Neurons_Scene,:);

PSTH_Fractal_FR_Transcience_mean=nanmean(PSTH_Fractal_FR_Transcience,1);
PSTH_Fractal_FR_Transcience_sem=nanstd(PSTH_Fractal_FR_Transcience,[],1)./sqrt(sum(~isnan(PSTH_Fractal_FR_Transcience),1));

PSTH_Fractal_FR_Sustained_mean=nanmean(PSTH_Fractal_FR_Sustained,1);
PSTH_Fractal_FR_Sustained_sem=nanstd(PSTH_Fractal_FR_Sustained,[],1)./sqrt(sum(~isnan(PSTH_Fractal_FR_Sustained),1));

PSTH_Scene_FR_Transcience_mean=nanmean(PSTH_Scene_FR_Transcience,1);
PSTH_Scene_FR_Transcience_sem=nanstd(PSTH_Scene_FR_Transcience,[],1)./sqrt(sum(~isnan(PSTH_Scene_FR_Transcience),1));

PSTH_Scene_FR_Sustained_mean=nanmean(PSTH_Scene_FR_Sustained,1);
PSTH_Scene_FR_Sustained_sem=nanstd(PSTH_Scene_FR_Sustained,[],1)./sqrt(sum(~isnan(PSTH_Scene_FR_Sustained),1));


if ShowFigureFlag
    
fig1=PrepareFigure(15,'w',[100,200, 1000,800],'Name','PSTH of Scene sustain neurons and Scene trainscient neurons');%,'Name','SceneTuning');
subplot(1,2,1)
%Transcient Scene Response
hold on
shadedErrorBar(PSTH_Scene_Time_mean,PSTH_Scene_FR_Transcience_mean,PSTH_Scene_FR_Transcience_sem,'lineprops',{'#4180FB'},'transparent',1,'patchSaturation',0.5);
%Sustained Scene Response
shadedErrorBar(PSTH_Scene_Time_mean,PSTH_Scene_FR_Sustained_mean,PSTH_Scene_FR_Sustained_sem,'lineprops',{'#FB2A2A'},'transparent',1,'patchSaturation',0.4);
 xlim([min(PSTH_Scene_Time_mean),max(PSTH_Scene_Time_mean)-20]);
%Labels
yrange=ylim;
h1=plot([0,0],[yrange(1),yrange(2)],'-.k');
h2=plot([1000,1000],[yrange(1),yrange(2)],'--k');

%Legend
Num_TranscientNeuron_str=sprintf('Transcient neurons(N=%d)',sum(Transcient_Neurons_Scene));
Num_SustainedNeuron_str=sprintf('Sustained neurons(N=%d)',sum(Sustained_Neurons_Scene));

g=legend({Num_TranscientNeuron_str,Num_SustainedNeuron_str},'color','none','Box', 'off','Position',[0.2    0.8    0.4510    0.1225]);
g.Orientation='Vertical';

set(get(get(h1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(h2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');


xlabel('Time from scene onset (ms)');
ylabel('Fring Rate (Hz)')
%title('PSTH for two types of scene responsive neurons')

set(gca,'LineWidth',3)
set(findall(gca, 'Type', 'Line'),'LineWidth',3);
set(findall(gcf,'-property','FontSize'),'FontSize',20);
set(findall(gcf,'-property','FontWeight'),'FontWeight','Bold');

subplot(1,2,2)
%Transcient Scene Response
hold on
shadedErrorBar(PSTH_Fractal_Time_mean,PSTH_Fractal_FR_Transcience_mean,PSTH_Fractal_FR_Transcience_sem,'lineprops',{'#4180FB'},'transparent',1,'patchSaturation',0.5);
%Sustained Scene Response
shadedErrorBar(PSTH_Fractal_Time_mean,PSTH_Fractal_FR_Sustained_mean,PSTH_Fractal_FR_Sustained_sem,'lineprops',{'#FB2A2A'},'transparent',1,'patchSaturation',0.4);
 
%Labels
%ylim(yrange);
xlim([-150,600]);
yrange=ylim;

h1=plot([0,0],[yrange(1),yrange(2)],'-.k');
h2=plot([400,400],[yrange(1),yrange(2)],'--k');



xlabel('Time from Fractal onset (ms)');
ylabel('Fring Rate (Hz)')
%title('PSTH for two types of scene responsive neurons')

set(gca,'LineWidth',3)
set(findall(gca, 'Type', 'Line'),'LineWidth',3);
set(findall(gcf,'-property','FontSize'),'FontSize',20);
set(findall(gcf,'-property','FontWeight'),'FontWeight','Bold');
 
    
    
    
end

end
