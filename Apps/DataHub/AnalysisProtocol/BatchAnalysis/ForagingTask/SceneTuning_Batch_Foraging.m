%Batch analysis for the scene response
%Xuefei Yu
%04092020
%Change into function and merge into batch analysis
function SceneTuning_Batch_Foraging(Data, StartFile, EndFile,ShowFigureFlag,OutputFlag);

%Batch data file path
FilesName=Data.ResultFilePath(StartFile: EndFile);
RecordDate=Data.RecordDate(StartFile: EndFile);


%Color library
%Color=[lines(6);prism(6);pink(6)];
%Color(2,:)=copper(1)+0.5;
%Color=hsv(1);
%Color=[Color;lines(6)];
%Color=[Color;hsv(10)];
%Color(18,:)=[0.1,0.457,0.194];
%Color=flipud(Color);


for i=1:length(FilesName)
  clear OutputData;
  
     load(FilesName{i});
     Task(i)=OutputData.SceneTuning(1).TaskCode;
     TaskName{i}=OutputData.SceneTuning(1).Task;
     
     DataLib{i}=OutputData.SceneTuning(1).DataStamp;
     EventSeg{i}=OutputData.EventSegregation(1).DataStamp;
     
    
end

%The overall scene response by event segmentation
p_SceneFromSpon_Whole_temp=cellfun(@(x) x('p_Event'),EventSeg,'Uniform',0);

h_SceneFromSpon_Whole=cell2mat(cellfun(@(x) x{1}(1),p_SceneFromSpon_Whole_temp,'Uniform',0));
p_SceneFromSpon_Whole=cell2mat(cellfun(@(x) x{2}(1),p_SceneFromSpon_Whole_temp,'Uniform',0));

%Find out spon for each trial from the event segration 
AverageFR_Spon_EachTrial=cellfun(@(x) nanmean(x('Spon_FR')),EventSeg,'Uniform',0);


PSTH_Event_EachTrial=cellfun(@(x) x('PSTH_Event'),EventSeg,'Uniform',0);
PSTH_EventTime_EachTrial=cellfun(@(x) x('PSTH_Event_Time'),EventSeg,'Uniform',0);

PSTH_Spon_EachTrial=cell2mat(cellfun(@(x) nanmean(x{end},1)',PSTH_Event_EachTrial,'Uniform',0))';
PSTH_SponTime_EachTrial_mean=nanmean(cell2mat(cellfun(@(x) nanmean(x{end})',PSTH_EventTime_EachTrial,'Uniform',0))',1);

PSTH_Spon_EachTrial_mean=nanmean(PSTH_Spon_EachTrial,1);
PSTH_Spon_EachTrial_sem=nanstd(PSTH_Spon_EachTrial,[],1)./sqrt(sum(~isnan(PSTH_Spon_EachTrial),1));
PSTH_SponTime_EachTrial_mean=PSTH_SponTime_EachTrial_mean-1000;


%PSTH_Spon_EachTrial=cellfun(@(x) x{end},PSTH_Event_EachTrial,'Uniform',0);
%PSTH_Event','PSTH_Event_Time'



%Find out scene responsive cell
AverageFR_Scene_EachTrial=cellfun(@(x) x('AverageFR_Scene_EachTrial'),DataLib,'Uniform',0);

AverageFR_Spon_EachTrial_PreScene=cellfun(@(x) x('AverageFR_Spon_EachTrial'),DataLib,'Uniform',0);
AverageFR_Spon_Mean_PreScene=cell2mat(cellfun(@nanmean ,AverageFR_Spon_EachTrial_PreScene,'Uniform',0));
%{
[h_Scene,p_Scene]=cellfun(@(x,y) ttest2(x,y),AverageFR_Scene_EachTrial,AverageFR_Spon_EachTrial,'Uniform',0);
%First criterium: fr significantly above spon
h_Scene=cell2mat(h_Scene);
p_Scene=cell2mat(p_Scene);

%Second criterium: fr significant in any one of the scene
%p_SceneFromSpon_EachScene=cellfun(@(x) x('p_SceneFromSpon_EachScene'),DataLib,'Uniform',0);
p_SceneFromSpon=cellfun(@(x) x('p_SceneFromSpon'),DataLib,'Uniform',0);
h_Scene_c2=cell2mat(cellfun(@(x) sum(x(2,:)==1)>0,p_SceneFromSpon,'Uniform',0));
%}
%SceneResponseNeuron=(h_Scene|h_Scene_c2|h_SceneFromSpon_Whole);
SceneResponseNeuron=h_SceneFromSpon_Whole==1;

%Retrieve the PSTH for the scene whole

PSTH_Scene_Whole=cellfun(@(x) x('PSTH_EachTrial'),DataLib,'Uniform',0);
PSTH_Scene_Time=cellfun(@(x) nanmean(x{1}),PSTH_Scene_Whole,'Uniform',0);
PSTH_Scene_FR=cellfun(@(x) nanmean(x{2}),PSTH_Scene_Whole,'Uniform',0);

Max_Col_Whole=max(cell2mat(cellfun(@numel ,PSTH_Scene_Time,'Uniform',0)));

PSTH_Scene_Time=cell2mat(cellfun(@(x) [x,NaN*ones(1,(Max_Col_Whole-numel(x)))]',PSTH_Scene_Time,'Uniform',0))';
PSTH_Scene_FR=cell2mat(cellfun(@(x) [x,NaN*ones(1,(Max_Col_Whole-numel(x)))]',PSTH_Scene_FR,'Uniform',0))';


PSTH_Scene_Time_mean=nanmean(PSTH_Scene_Time);


PSTH_Scene_responsive_FR=PSTH_Scene_FR(SceneResponseNeuron,:);
PSTH_Scene_FR_mean_Scene=nanmean(PSTH_Scene_FR(SceneResponseNeuron,:));
PSTH_Scene_FR_sem_Scene=nanstd(PSTH_Scene_FR(SceneResponseNeuron,:))./sqrt(sum(~isnan(PSTH_Scene_FR(SceneResponseNeuron,:)),1));

%Spon For Scene responsive neruons
AverageFR_Spon_Mean_PreScene_Scene=AverageFR_Spon_Mean_PreScene(SceneResponseNeuron);


PSTH_Scene_FR_mean_NonScene=nanmean(PSTH_Scene_FR(~SceneResponseNeuron,:));
PSTH_Scene_FR_sem_NonScene=nanstd(PSTH_Scene_FR(~SceneResponseNeuron,:))./sqrt(sum(~isnan(PSTH_Scene_FR(~SceneResponseNeuron,:)),1));

Num_SceneNeuron=sum(SceneResponseNeuron);
Num_NonSceneNeuron=sum(~SceneResponseNeuron);

%Get the 'ResponseProportion' response latency,  'PeakFiring_SceneAverage','PeakFiringTime_SceneAverage','FirstSigInterval','PreferScene'
ResponseLatency_Scene_EachNeuron=cell2mat(cellfun(@(x) x('ResponseLatency_Average'),DataLib,'Uniform',0));
ResponseProportion_Scene_EachNeuron=cell2mat(cellfun(@(x) x('ResponseProportion'),DataLib,'Uniform',0));
AveragePeakFiring_Scene_EachNeuron=cell2mat(cellfun(@(x) x('PeakFiring_SceneAverage'),DataLib,'Uniform',0));
AveragePeakFiringTime_Scene_EachNeuron=cell2mat(cellfun(@(x) x('PeakFiringTime_SceneAverage'),DataLib,'Uniform',0));
FirstSigInterval_Scene_EachNeuron=cell2mat(cellfun(@(x) x('FirstSigInterval'),DataLib,'Uniform',0));


PreferScene_Scene_EachNeuron=cell2mat(cellfun(@(x) x('PreferScene'),DataLib,'Uniform',0));

AllUniqueScene_tested_Raw=cellfun(@(x) x('UniqueScene')',DataLib,'Uniform',0);
AllUniqueScene_tested=unique(cell2mat(AllUniqueScene_tested_Raw));
%Transform to distribution
ResponseLatency_Scene_RespNeuron=ResponseLatency_Scene_EachNeuron(SceneResponseNeuron);
ResponseProportion_Scene_RespNeuron=ResponseProportion_Scene_EachNeuron(SceneResponseNeuron);
AveragePeakFiringTime_Scene_RespNeuron_Rel=AveragePeakFiringTime_Scene_EachNeuron(SceneResponseNeuron)-ResponseLatency_Scene_EachNeuron(SceneResponseNeuron);
FirstSigInterval_Scene_RespNeuron=FirstSigInterval_Scene_EachNeuron(SceneResponseNeuron);

AveragePeakFiring_Scene_RespNeuron=AveragePeakFiring_Scene_EachNeuron(SceneResponseNeuron);

PreferScene_Scene_RespNeuron=PreferScene_Scene_EachNeuron(SceneResponseNeuron);

%%
%Separate neurons into 2 groups with different type of scene trigerred
%responses; Threshold: ResponseProportion_Scene_RespNeuron=30
Transcient_Neurons_Scene=ResponseProportion_Scene_RespNeuron<=30;
Sustained_Neurons_Scene=ResponseProportion_Scene_RespNeuron>30;

PSTH_Transcient_Scene=PSTH_Scene_responsive_FR(Transcient_Neurons_Scene,:);
PSTH_Sustained_Scene=PSTH_Scene_responsive_FR(Sustained_Neurons_Scene,:);

PSTH_Transcient_Scene_mean=nanmean(PSTH_Transcient_Scene);
PSTH_Transcient_Scene_sem=nanstd(PSTH_Transcient_Scene)./sqrt(sum(~isnan(PSTH_Transcient_Scene),1));

PSTH_Sustained_Scene_mean=nanmean(PSTH_Sustained_Scene);
PSTH_Sustained_Scene_sem=nanstd(PSTH_Sustained_Scene)./sqrt(sum(~isnan(PSTH_Sustained_Scene),1));

%p_Value for the different between the sustained and the transcient group
[p_Sus_Tr_Diff,h_Sus_Tr_Diff]=ttest2_group(PSTH_Transcient_Scene,PSTH_Sustained_Scene,'both');
%For plot
h_Sus_Tr_Diff(h_Sus_Tr_Diff==0)=NaN;

%Normalize the responses for the two type of neurons by the spon and the
%peak
%Use the prescene spon to normalize the responses
PSTH_Transcient_Scene_spon=AverageFR_Spon_Mean_PreScene_Scene(Transcient_Neurons_Scene);
PSTH_Sustained_Scene_spon=AverageFR_Spon_Mean_PreScene_Scene(Sustained_Neurons_Scene);


PSTH_Transcient_Scene_Norm=(PSTH_Transcient_Scene-PSTH_Transcient_Scene_spon')./repmat(AveragePeakFiring_Scene_RespNeuron(Transcient_Neurons_Scene)',1,size(PSTH_Transcient_Scene,2));
PSTH_Sustained_Scene_Norm=(PSTH_Sustained_Scene-PSTH_Sustained_Scene_spon')./repmat(AveragePeakFiring_Scene_RespNeuron(Sustained_Neurons_Scene)',1,size(PSTH_Sustained_Scene,2));
%PSTH_Transcient_Scene_Norm=(PSTH_Transcient_Scene-PSTH_Transcient_Scene_spon');
%PSTH_Sustained_Scene_Norm=(PSTH_Sustained_Scene-PSTH_Sustained_Scene_spon');

%p_Value for the different between the transcient and the spon
%for normalized condition
[p_Sus_Tr_Diff_Transcient_Norm,h_Sus_Tr_Diff_Transcient_Norm]=ttest2_group(PSTH_Transcient_Scene_Norm,0*ones(size(PSTH_Transcient_Scene_Norm,1),1),'both');
%For plot
h_Sus_Tr_Diff_Transcient_Norm(h_Sus_Tr_Diff_Transcient_Norm==0)=NaN;

%p_Value for the different between the sustained and the spon
%for normalized condition
[p_Sus_Tr_Diff_Sustained_Norm,h_Sus_Tr_Diff_Sustained_Norm]=ttest2_group(PSTH_Sustained_Scene_Norm,0*ones(size(PSTH_Sustained_Scene_Norm,1),1),'both');
%For plot
h_Sus_Tr_Diff_Sustained_Norm(h_Sus_Tr_Diff_Sustained_Norm==0)=NaN;




%p_Value for the different between the sustained and the transcient group
%for normalized condition
[p_Sus_Tr_Diff_Norm,h_Sus_Tr_Diff_Norm]=ttest2_group(PSTH_Transcient_Scene_Norm,PSTH_Sustained_Scene_Norm,'both');
%For plot
h_Sus_Tr_Diff_Norm(h_Sus_Tr_Diff_Norm==0)=NaN;



PSTH_Transcient_Scene_Norm_mean=nanmean(PSTH_Transcient_Scene_Norm,1);
PSTH_Transcient_Scene_Norm_sem=nanstd(PSTH_Transcient_Scene_Norm,[],1)./sqrt(sum(~isnan(PSTH_Transcient_Scene_Norm),1));

PSTH_Sustained_Scene_Norm_mean=nanmean(PSTH_Sustained_Scene_Norm,1);
PSTH_Sustained_Scene_Norm_sem=nanstd(PSTH_Sustained_Scene_Norm,[],1)./sqrt(sum(~isnan(PSTH_Sustained_Scene_Norm),1));

%%
%Properties for two different kind of scene responsive neurons
%Response latency
ResponseLatency_Transcient=ResponseLatency_Scene_RespNeuron(Transcient_Neurons_Scene);
ResponseLatency_Sustained=ResponseLatency_Scene_RespNeuron(Sustained_Neurons_Scene);

Mean_RL_Sus_Trans=[nanmean(ResponseLatency_Sustained),nanmean(ResponseLatency_Transcient)];
Sem_RL_Sus_Trans=[nanstd(ResponseLatency_Sustained)/sqrt(sum(~isnan(ResponseLatency_Sustained))),nanstd(ResponseLatency_Transcient)/sqrt(sum(~isnan(ResponseLatency_Transcient)))];

%Significance test
[h_RL_sus_trans,p_RL_sus_trans]=ttest2(ResponseLatency_Transcient,ResponseLatency_Sustained);



%First Interval Duration
FirstInterval_Transcient=FirstSigInterval_Scene_RespNeuron(Transcient_Neurons_Scene);
FirstInterval_Sustained=FirstSigInterval_Scene_RespNeuron(Sustained_Neurons_Scene);

Mean_FI_Sus_Trans=[nanmean(FirstInterval_Sustained),nanmean(FirstInterval_Transcient)];
Sem_FI_Sus_Trans=[nanstd(FirstInterval_Sustained)/sqrt(sum(~isnan(FirstInterval_Sustained))),nanstd(FirstInterval_Transcient)/sqrt(sum(~isnan(FirstInterval_Transcient)))];

%Significance test
[h_FI_sus_trans,p_FI_sus_trans]=ttest2(FirstInterval_Transcient,FirstInterval_Sustained);

%Peak Firing 

PeakFiring_Transcient=AveragePeakFiring_Scene_RespNeuron(Transcient_Neurons_Scene);
PeakFiring_Sustained=AveragePeakFiring_Scene_RespNeuron(Sustained_Neurons_Scene);

Mean_PF_Sus_Trans=[nanmean(PeakFiring_Sustained),nanmean(PeakFiring_Transcient)];
Sem_PF_Sus_Trans=[nanstd(PeakFiring_Sustained)/sqrt(sum(~isnan(PeakFiring_Sustained))),nanstd(PeakFiring_Transcient)/sqrt(sum(~isnan(PeakFiring_Transcient)))];

%Significance test
[h_PF_sus_trans,p_PF_sus_trans]=ttest2(PeakFiring_Transcient,PeakFiring_Sustained);



%Peak Firing Time Relative to response latency 

Rel_PeakFiringTime_Transcient=AveragePeakFiringTime_Scene_RespNeuron_Rel(Transcient_Neurons_Scene);
Rel_PeakFiringTime_Sustained=AveragePeakFiringTime_Scene_RespNeuron_Rel(Sustained_Neurons_Scene);

Mean_RelPFT_Sus_Trans=[nanmean(Rel_PeakFiringTime_Sustained),nanmean(Rel_PeakFiringTime_Transcient)];
Sem_RelPFT_Sus_Trans=[nanstd(Rel_PeakFiringTime_Sustained)/sqrt(sum(~isnan(Rel_PeakFiringTime_Sustained))),nanstd(Rel_PeakFiringTime_Transcient)/sqrt(sum(~isnan(Rel_PeakFiringTime_Transcient)))];

%Significance test
[h_RelPFT_sus_trans,p_RelPFT_sus_trans]=ttest2(Rel_PeakFiringTime_Transcient,Rel_PeakFiringTime_Sustained);



%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Analysis the scene selective=%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% responses .%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SSI=cell2mat(cellfun(@(x) x('SceneSelectiveIndex'),DataLib,'Uniform',0));
p_SSI=cell2mat(cellfun(@(x) x('p_SSI'),DataLib,'Uniform',0));

SponSI=cell2mat(cellfun(@(x) x('SponSelectiveIndex'),DataLib,'Uniform',0));
p_SponSI=cell2mat(cellfun(@(x) x('p_SponSSI'),DataLib,'Uniform',0));


PSTH_SSI=cellfun(@(x) x('PSTH_SSI'),DataLib,'Uniform',0);



p_SceneFromSpon_EachScene_temp=cellfun(@(x) x('p_SceneFromSpon'),DataLib,'Uniform',0);
p_SceneFromSpon_EachScene=cellfun(@(x) x(:,2),p_SceneFromSpon_EachScene_temp,'Uniform',0);
h_SceneFromSpon_EachScene=cellfun(@(x) x(:,1),p_SceneFromSpon_EachScene_temp,'Uniform',0);

NumSigScene_EachNeuron=cell2mat(cellfun(@(x) sum(x,'omitnan'),h_SceneFromSpon_EachScene,'Uniform',0));
TotalNumOfSceneTested=cell2mat(cellfun(@(x) sum(~isnan(x)),h_SceneFromSpon_EachScene,'Uniform',0));

PropSigScene_EachNeuron=NumSigScene_EachNeuron./TotalNumOfSceneTested;


%%Screen out the scene responsive neurons

SSI_Scene=SSI(SceneResponseNeuron);
p_SSI_Scene=p_SSI(SceneResponseNeuron);

PSTH_SSI_Scene=PSTH_SSI(SceneResponseNeuron);
PSTH_SSI_Scene=cell2mat(cellfun(@(x) [x,NaN*ones(1,(Max_Col_Whole-numel(x)))]',PSTH_SSI_Scene,'Uniform',0))';

%PSTH_Scene_Time=cell2mat(cellfun(@(x) [x,NaN*ones(1,(Max_Col_Whole-numel(x)))]',PSTH_Scene_Time,'Uniform',0))';


SponSI_Scene=SponSI(SceneResponseNeuron);
p_SponSI_Scene=p_SponSI(SceneResponseNeuron);

PropSigScene_Scene=PropSigScene_EachNeuron(SceneResponseNeuron);


%Define the scene selective neurons
SS_C=p_SSI_Scene<0.05;
%SS_C=PropSigScene_Scene>0&PropSigScene_Scene<1;
%SS_C=SSI_Scene==max(SSI_Scene);
%SS_C=PropSigScene_Scene<1;
%nonSS_C=PropSigScene_Scene==1;
nonSS_C=~SS_C;
%nonSS_C=SSI_Scene==min(SSI_Scene);;


PSTH_SSI_Scene_SSN=PSTH_SSI_Scene(SS_C,:);
PSTH_SSI_Scene_NonSSN=PSTH_SSI_Scene(nonSS_C,:);


Mean_PSTH_SSI_Scene_SSN=nanmean(PSTH_SSI_Scene_SSN,1);
Sem_PSTH_SSI_Scene_SSN=nanstd(PSTH_SSI_Scene_SSN,[],1)./sqrt(sum(~isnan(PSTH_SSI_Scene_SSN),1));

Mean_PSTH_SSI_Scene_NonSSN=nanmean(PSTH_SSI_Scene_NonSSN,1);
Sem_PSTH_SSI_Scene_NonSSN=nanstd(PSTH_SSI_Scene_NonSSN,[],1)./sqrt(sum(~isnan(PSTH_SSI_Scene_NonSSN),1));

[p_diff_SSN_NonSSN,h_diff_SSN_NonSSN]=ttest2_group(PSTH_SSI_Scene_SSN,PSTH_SSI_Scene_NonSSN,'right');

h_diff_SSN_NonSSN(h_diff_SSN_NonSSN==0)=NaN;

%{
figure
shadedErrorBar(PSTH_Scene_Time_mean,nanmean(PSTH_SSI_Scene_SSN,1),nanstd(PSTH_SSI_Scene_SSN,[],1)./sqrt(sum(~isnan(PSTH_SSI_Scene_SSN),1)),'lineprops',{'#DC4826'},'transparent',1,'patchSaturation',0.5);
hold on
shadedErrorBar(PSTH_Scene_Time_mean,nanmean(PSTH_SSI_Scene_NonSSN,1),nanstd(PSTH_SSI_Scene_NonSSN,[],1)./sqrt(sum(~isnan(PSTH_SSI_Scene_NonSSN),1)),'lineprops',{'#648FFF'},'transparent',1,'patchSaturation',0.4);
%}
%legend({'SSN','NonSSN'})



PropSigScene_Scene_SSN=PropSigScene_Scene(SS_C);

PropSigScene_Scene_SSN_mean=mean(PropSigScene_Scene_SSN);


%Most preferred scene for the scene responsive neurons
AllUniqueScene_tested_SceneNeuron=AllUniqueScene_tested_Raw(SceneResponseNeuron);

AllUniqueScene_tested_SceneNeuron_SSN=AllUniqueScene_tested_SceneNeuron(SS_C);


UniqueScenesTested_SceneNeuron_SSN=unique(cell2mat(AllUniqueScene_tested_SceneNeuron_SSN));%All scenes tested for those neurons


PreferScene_Scene_RespNeuron_SSN=PreferScene_Scene_RespNeuron(SS_C);


%PreferScene_Scene_RespNeuron

%PreferScene_Scene_RespNeuron
%PreferScene_Scene_RespNeuron

%Comparison between transcient neurons and sustained neurons

SSI_transcient=SSI_Scene(Transcient_Neurons_Scene);
SSI_sustained=SSI_Scene(Sustained_Neurons_Scene);

Mean_SSI_Sus_Trans=[nanmean(SSI_sustained),nanmean(SSI_transcient)];
Sem_SSI_Sus_Trans=[nanstd(SSI_sustained)/sqrt(sum(~isnan(SSI_sustained))),nanstd(SSI_transcient)/sqrt(sum(~isnan(SSI_transcient)))];

%Significance test
[h_SSI_sus_trans,p_SSI_sus_trans]=ttest2(SSI_transcient,SSI_sustained);


if ShowFigureFlag

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Plot the average PSTH for each type of neurons
fig1=PrepareFigure(1,'w',[100,200, 1000,800],'Name','PSTH of Scene response neurons and NonScene response neurons');%,'Name','SceneTuning');
%Plot spon
shadedErrorBar(PSTH_SponTime_EachTrial_mean,PSTH_Spon_EachTrial_mean,PSTH_Spon_EachTrial_sem,'lineprops',{'#000000'},'transparent',1,'patchSaturation',0.2);
%Scene Response
hold on
shadedErrorBar(PSTH_Scene_Time_mean,PSTH_Scene_FR_mean_Scene,PSTH_Scene_FR_sem_Scene,'lineprops',{'#DC4826'},'transparent',1,'patchSaturation',0.5);
%Non scene response
shadedErrorBar(PSTH_Scene_Time_mean,PSTH_Scene_FR_mean_NonScene,PSTH_Scene_FR_sem_NonScene,'lineprops',{'#648FFF'},'transparent',1,'patchSaturation',0.4);

%Labels
yrange=ylim;
h1=plot([0,0],[yrange(1),yrange(2)],'-.k');
h2=plot([1000,1000],[yrange(1),yrange(2)],'--k');

%Legend
Num_SceneNeuron_str=sprintf('Scene Responsive\nneurons(N=%d)',Num_SceneNeuron);
Num_NonSceneNeuron_str=sprintf('NonScene Responsive\nneurons(N=%d)',Num_NonSceneNeuron);

g=legend({'Spon',Num_SceneNeuron_str,Num_NonSceneNeuron_str},'color','none','Box', 'off','Position',[0.38    0.75    0.4510    0.1225]);

set(get(get(h1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(h2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');



xlim([min(PSTH_Scene_Time_mean),max(PSTH_Scene_Time_mean)-20]);
xlabel('Time from scene onset (ms)');
ylabel('Fring Rate(Hz)')
title('PSTH for all neurons aligned on scene onset')

set(gca,'LineWidth',3)
set(findall(gca, 'Type', 'Line'),'LineWidth',3);
set(findall(gcf,'-property','FontSize'),'FontSize',20);
set(findall(gcf,'-property','FontWeight'),'FontWeight','Bold');
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Histogrames for the responselatency, Proportion of
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% significant duration, Peak
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Firing Rate,
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ColorForHist={'#FFB000','#FE6100','#648FFF','#785EF0'};
fig2=PrepareFigure(2,'w',[100,200, 1000,800],'Name','Distribution of response latency for scene responsive neurons');%,'Name','SceneTuning');
subplot(2,2,1);
histogram(ResponseLatency_Scene_RespNeuron,'BinWidth',20,'FaceColor',ColorForHist{1},'EdgeColor', ColorForHist{1},'EdgeAlpha',0.5,'FaceAlpha',0.5);
set(gca,'LineWidth',3);
xlabel('Scene response latency (ms)');
ylabel('Number of neurons');
%set(gca,'FontSize',20);
%legend(sprintf('BinWidth %d',20),'FontSize',15)
box off

subplot(2,2,2);histogram(AveragePeakFiring_Scene_RespNeuron,'BinWidth',20,'FaceColor',ColorForHist{4},'EdgeColor', ColorForHist{4},'EdgeAlpha',0.5,'FaceAlpha',0.5);
set(gca,'LineWidth',3);
xlabel('Distribution of Peak firing(Hz)');
%ylabel('Number of neurons');
%legend(sprintf('BinWidth %d',20))
box off


subplot(2,2,3);histogram(ResponseProportion_Scene_RespNeuron,'BinWidth',5,'FaceColor',ColorForHist{2},'EdgeColor', ColorForHist{2},'EdgeAlpha',0.5,'FaceAlpha',0.5);
set(gca,'LineWidth',3);
xlabel(['Percentage of time with ',newline,'significant scene responses (%)']);
ylabel('Number of neurons');
%legend(sprintf('BinWidth %d',5),'FontSize',15)
box off



subplot(2,2,4);histogram(FirstSigInterval_Scene_RespNeuron,'BinWidth',50,'FaceColor',ColorForHist{3},'EdgeColor', ColorForHist{3},'EdgeAlpha',0.5,'FaceAlpha',0.5);
set(gca,'LineWidth',3);
xlabel(['First interval duration with',newline,'significant scene responses(ms)'])
%legend(sprintf('BinWidth %d',50))    
%ylabel('Number of neurons');
box off




set(findall(gca, 'Type', 'Line'),'LineWidth',3);
set(findall(gcf,'-property','FontSize'),'FontSize',20);
set(findall(gcf,'-property','FontWeight'),'FontWeight','Bold');
%fig3=PrepareFigure(3,'w',[100,200, 500,500],'Name','Distribution of prefered scene for scene responsive neurons');%,'Name','SceneTuning');
%histogram(PreferScene_Scene_RespNeuron,'BinWidth',1);


%subplot(2,2,3);histogram(AveragePeakFiringTime_Scene_RespNeuron_Rel,'BinWidth',5);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%PSTH for two type of scene responsive neurons%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%neurons%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fig3=PrepareFigure(3,'w',[100,200, 1200,1000],'Name','PSTH for two type of scene responsive neurons');%,'Name','SceneTuning');
%Raw responses
%subplot(1,2,1)
axes('Position',[0.06,0.2,0.42,0.6]);
%Transient Scene Response

shadedErrorBar(PSTH_Scene_Time_mean,PSTH_Transcient_Scene_mean,PSTH_Transcient_Scene_sem,'lineprops',{'#4180FB'},'transparent',1,'patchSaturation',0.5);
hold on
%Sustained Scene Response
shadedErrorBar(PSTH_Scene_Time_mean,PSTH_Sustained_Scene_mean,PSTH_Sustained_Scene_sem,'lineprops',{'#FB2A2A'},'transparent',1,'patchSaturation',0.4);
yrange=ylim;
%Add signifcance interval
sig=area(PSTH_Scene_Time_mean,h_Sus_Tr_Diff*max(min(PSTH_Sustained_Scene_mean),1)*0.2,'FaceColor','#FFC107','EdgeColor','#FFC107','FaceAlpha',0.5,'EdgeAlpha',0.5);

%Labels

h1=plot([0,0],[yrange(1),yrange(2)],'-.k');
h2=plot([1000,1000],[yrange(1),yrange(2)],'--k');

%Legend
Num_TranscientNeuron_str=sprintf('Transcient neurons(N=%d)',sum(Transcient_Neurons_Scene));
Num_SustainedNeuron_str=sprintf('Sustained neurons(N=%d)',sum(Sustained_Neurons_Scene));

g=legend({Num_TranscientNeuron_str,Num_SustainedNeuron_str},'color','none','Box', 'off','Position',[0.2    0.8    0.4510    0.1225]);
g.Orientation='Horizontal';

set(get(get(h1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(h2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');

xlim([min(PSTH_Scene_Time_mean),max(PSTH_Scene_Time_mean)-20]);
xlabel('Time from scene onset (ms)');
ylabel('Fring Rate (Hz)')
%title('PSTH for two types of scene responsive neurons')

set(gca,'LineWidth',3)
set(findall(gca, 'Type', 'Line'),'LineWidth',3);
set(findall(gcf,'-property','FontSize'),'FontSize',20);
set(findall(gcf,'-property','FontWeight'),'FontWeight','Bold');


%Normalized responses
%subplot(1,2,2)
axes('Position',[0.54,0.2,0.42,0.6])
%Transient Scene Response

shadedErrorBar(PSTH_Scene_Time_mean,PSTH_Transcient_Scene_Norm_mean,PSTH_Transcient_Scene_Norm_sem,'lineprops',{'#4180FB'},'transparent',1,'patchSaturation',0.5);
hold on
%Sustained Scene Response
shadedErrorBar(PSTH_Scene_Time_mean,PSTH_Sustained_Scene_Norm_mean,PSTH_Sustained_Scene_Norm_sem,'lineprops',{'#FB2A2A'},'transparent',1,'patchSaturation',0.4);
ylim([-0.2,0.7]);

yrange=ylim;

%{

%Add signifcance interval for difference between sustained and spon
sig_sus=area(PSTH_Scene_Time_mean,-h_Sus_Tr_Diff_Sustained_Norm*0.15,-0.2,'FaceColor','#FB2A2A','EdgeColor','#FB2A2A','FaceAlpha',0.5,'EdgeAlpha',0.5);

%Add signifcance interval for difference between transcient and spon

sig_tran=area(PSTH_Scene_Time_mean,-h_Sus_Tr_Diff_Transcient_Norm*0.2,-0.25,'FaceColor','#4180FB','EdgeColor','#4180FB','FaceAlpha',0.5,'EdgeAlpha',0.5);
%}

%Add signifcance interval for difference between transcient and sustained
sig_diff=area(PSTH_Scene_Time_mean,-h_Sus_Tr_Diff_Norm*0.15,-0.2,'FaceColor','#FFC107','EdgeColor','#FFC107','FaceAlpha',0.5,'EdgeAlpha',0.5,'ShowBaseLine','off');

%Labels

plot([0,0],[yrange(1),yrange(2)],'-.k');
plot([1000,1000],[yrange(1),yrange(2)],'--k');

xlim([min(PSTH_Scene_Time_mean),max(PSTH_Scene_Time_mean)-20]);
plot([min(PSTH_Scene_Time_mean),max(PSTH_Scene_Time_mean)-20],[0,0],'--k');




xlabel('Time from scene onset (ms)');
ylabel('Nomalized Response')
%title('PSTH for two types of scene responsive neurons')

set(gca,'LineWidth',3)
set(findall(gca, 'Type', 'Line'),'LineWidth',3);
set(findall(gcf,'-property','FontSize'),'FontSize',20);
set(findall(gcf,'-property','FontWeight'),'FontWeight','Bold');

a=annotation('textbox',[0.2,0.45,0.5,0.5],'String','PSTH for two types of scene responsive neurons','FitBoxToText','on','LineStyle','none');
a.FontSize=25;
a.FontWeight='Bold';


%%
%%Porperties for two types of scene responsive neurons
fig4=PrepareFigure(4,'w',[100,200, 500,500],'Name','Porperties for two types of scene responsive neurons');%,'Name','SceneTuning');

%Response Latency
subplot(2,2,1);
b1=bar(Mean_RL_Sus_Trans);
b1.FaceColor = 'flat';
b1.EdgeColor = 'flat';
b1.CData = [	251,42,42;65,128,251]/255;

hold on
errorbar(Mean_RL_Sus_Trans,Sem_RL_Sus_Trans,'Color','k','LineStyle','none','LineWidth',3);
yrange=ylim;

%Add significance level
 

if h_RL_sus_trans==1
    
    %,p_RL_sus_trans
ctr2 = bsxfun(@plus, b1.XData, [b1.XOffset]');
hold on

plot(ctr2(1:2), [1 1]*yrange(2)*1.1, '-k', 'LineWidth',2)
if p_RL_sus_trans<0.05 & p_RL_sus_trans>=0.01
   plot(mean(ctr2(1:2)), yrange(2)*1.15, '*k','MarkerSize',10)
elseif p_RL_sus_trans<0.01 & p_RL_sus_trans>=0.001
    plot(mean(ctr2(1:2))*0.9, yrange(2)*1.15, '*k','MarkerSize',10)
    plot(mean(ctr2(1:2))*1.1, yrange(2)*1.15, '*k','MarkerSize',10)
        
elseif p_RL_sus_trans<0.001
    plot(mean(ctr2(1:2))*0.8, yrange(2)*1.15, '*k','MarkerSize',10)
    plot(mean(ctr2(1:2))*1, yrange(2)*1.15, '*k','MarkerSize',10)
    plot(mean(ctr2(1:2))*1.2, yrange(2)*1.15, '*k','MarkerSize',10)
    
end


end
xticklabels({'Sus','Tran'})
ylabel('Response Latency (ms)')

box off

set(gca,'LineWidth',3)
%set(findall(gca, 'Type', 'Line'),'LineWidth',3);


%%%%%%%%%%%%%%%%%%%%%%%%%%%First interval duration
subplot(2,2,2);
b1=bar(Mean_FI_Sus_Trans);
b1.FaceColor = 'flat';
b1.EdgeColor = 'flat';
b1.CData = [	251,42,42;65,128,251]/255;

hold on
errorbar(Mean_FI_Sus_Trans,Sem_FI_Sus_Trans,'Color','k','LineStyle','none','LineWidth',3);
yrange=ylim;


%Add significance level
 

if h_FI_sus_trans==1
    
    %,p_RL_sus_trans
ctr2 = bsxfun(@plus, b1.XData, [b1.XOffset]');
hold on

plot(ctr2(1:2), [1 1]*yrange(2)*1.1, '-k', 'LineWidth',2)
if p_FI_sus_trans<0.05 & p_FI_sus_trans>=0.01
   plot(mean(ctr2(1:2)), yrange(2)*1.15, '*k','MarkerSize',10)
elseif p_FI_sus_trans<0.01 & p_FI_sus_trans>=0.001
    plot(mean(ctr2(1:2))*0.9, yrange(2)*1.15, '*k','MarkerSize',10)
    plot(mean(ctr2(1:2))*1.1, yrange(2)*1.15, '*k','MarkerSize',10)
        
elseif p_FI_sus_trans<0.001
    plot(mean(ctr2(1:2))*0.8, yrange(2)*1.15, '*k','MarkerSize',10)
    plot(mean(ctr2(1:2))*1, yrange(2)*1.15, '*k','MarkerSize',10)
    plot(mean(ctr2(1:2))*1.2, yrange(2)*1.15, '*k','MarkerSize',10)
    
end


end
xticklabels({'Sus','Tran'})
ylabel('First Interval Dur(ms)')

box off

set(gca,'LineWidth',3)





%%%%%%%%%%%%%%%%%%%%%%%%%%% Peak Firing Rate %%%%%%%%%%%%%%%%%%%%%%%%%
subplot(2,2,3);
b1=bar(Mean_PF_Sus_Trans);
b1.FaceColor = 'flat';
b1.EdgeColor = 'flat';
b1.CData = [	251,42,42;65,128,251]/255;

hold on
errorbar(Mean_PF_Sus_Trans,Sem_PF_Sus_Trans,'Color','k','LineStyle','none','LineWidth',3);
yrange=ylim;


%Add significance level
 

if h_PF_sus_trans==1
    
    %,p_RL_sus_trans
ctr2 = bsxfun(@plus, b1.XData, [b1.XOffset]');
hold on

plot(ctr2(1:2), [1 1]*yrange(2)*1.1, '-k', 'LineWidth',2)
if p_PF_sus_trans<0.05 & p_PF_sus_trans>=0.01
   plot(mean(ctr2(1:2)), yrange(2)*1.15, '*k','MarkerSize',10)
elseif p_PF_sus_trans<0.01 & p_PF_sus_trans>=0.001
    plot(mean(ctr2(1:2))*0.9, yrange(2)*1.15, '*k','MarkerSize',10)
    plot(mean(ctr2(1:2))*1.1, yrange(2)*1.15, '*k','MarkerSize',10)
        
elseif p_PF_sus_trans<0.001
    plot(mean(ctr2(1:2))*0.8, yrange(2)*1.15, '*k','MarkerSize',10)
    plot(mean(ctr2(1:2))*1, yrange(2)*1.15, '*k','MarkerSize',10)
    plot(mean(ctr2(1:2))*1.2, yrange(2)*1.15, '*k','MarkerSize',10)
    
end


end
xticklabels({'Sus','Tran'})
ylabel('Peak Firing (Hz)')

box off

set(gca,'LineWidth',3)




%%%%%%%%%%%%%%%%%%%%%%%%%%% Relative Peak Firing Time %%%%%%%%%%%%%%%%%%%%%%%%%
subplot(2,2,4);
b1=bar(Mean_RelPFT_Sus_Trans);
b1.FaceColor = 'flat';
b1.EdgeColor = 'flat';
b1.CData = [	251,42,42;65,128,251]/255;

hold on
errorbar(Mean_RelPFT_Sus_Trans,Sem_RelPFT_Sus_Trans,'Color','k','LineStyle','none','LineWidth',3);
yrange=ylim;


%Add significance level
 

if h_RelPFT_sus_trans==1
    
    %,p_RL_sus_trans
ctr2 = bsxfun(@plus, b1.XData, [b1.XOffset]');
hold on

plot(ctr2(1:2), [1 1]*yrange(2)*1.1, '-k', 'LineWidth',2)
if p_RelPFT_sus_trans<0.05 & p_RelPFT_sus_trans>=0.01
   plot(mean(ctr2(1:2)), yrange(2)*1.15, '*k','MarkerSize',10)
elseif p_RelPFT_sus_trans<0.01 & p_RelPFT_sus_trans>=0.001
    plot(mean(ctr2(1:2))*0.9, yrange(2)*1.15, '*k','MarkerSize',10)
    plot(mean(ctr2(1:2))*1.1, yrange(2)*1.15, '*k','MarkerSize',10)
        
elseif p_RelPFT_sus_trans<0.001
    plot(mean(ctr2(1:2))*0.8, yrange(2)*1.15, '*k','MarkerSize',10)
    plot(mean(ctr2(1:2))*1, yrange(2)*1.15, '*k','MarkerSize',10)
    plot(mean(ctr2(1:2))*1.2, yrange(2)*1.15, '*k','MarkerSize',10)
    
end


end
xticklabels({'Sus','Tran'})
ylabel('Rel Peak Time (ms)')

box off

set(gca,'LineWidth',3)


set(findall(gca, 'Type', 'Line'),'LineWidth',3);


set(findall(gcf,'-property','FontSize'),'FontSize',20);
set(findall(gcf,'-property','FontWeight'),'FontWeight','Bold');

%%%%%%%%%%%%%%%%%%%%%%%%%%% Criterium for selection of scene selective neurons %%%%%%%%%%%%%%%%%%%%%%%%%
%close([1,2,3,4]);
fig5=PrepareFigure(5,'w',[100,200, 800,600],'Name','Criterium for the selection of scene selective neurons');%,'Name','SceneTuning');

%Histogram for SSI
%SSI_Scene=SSI(SceneResponseNeuron);
%p_SSI_Scene=p_SSI(SceneResponseNeuron);

%For Scene Selective Index
SSI_Scene_Sig=SSI_Scene(p_SSI_Scene<0.05);
SSI_Scene_NonSig=SSI_Scene(~(p_SSI_Scene<0.05));

SSI_Range=linspace(min(SSI_Scene)*0.9,max(SSI_Scene),10);

%For Spon Selective Index
SponSI_Scene_Sig=SponSI_Scene(p_SponSI_Scene<0.05);
SponSI_Scene_NonSig=SponSI_Scene(~(p_SponSI_Scene<0.05));



%subplot(2,1,1)
axes('Position',[0.1,0.7,0.6,0.25])
h1=histogram(SSI_Scene_Sig,SSI_Range);
hold on
h2=histogram(SSI_Scene_NonSig,SSI_Range);

h3=histogram(SponSI_Scene_Sig,SSI_Range);

h4=histogram(SponSI_Scene_NonSig,SSI_Range);


Count=[h1.Values',h2.Values'];
Count_Spon=[h3.Values',h4.Values'];


hold off

%subplot(2,2,1)


b=bar([SSI_Range(1:end-1)+diff(SSI_Range);SSI_Range(1:end-1)+diff(SSI_Range)]',Count,'Stack','LineWidth',3);
b(1).FaceColor='#039E20';
b(2).FaceColor='#676968';
legend({'Sig','NonSig'},'FontSize',15,'FontWeight','Bold');

set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold')


xlabel('Scene Selective Index','FontSize',20,'FontWeight','Bold')
ylabel('Num of neurons','FontSize',20,'FontWeight','Bold')

box off

axes('Position',[0.8,0.7,0.15,0.25])
bar(PropSigScene_Scene_SSN_mean,'FaceColor','#039E20');

xticklabels('Scene Selective Neurons')

box off
%xlabel('Scene Selective Neurons','FontSize',20,'FontWeight','Bold')
set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold')
ylabel(['Average Prop',newline,' of sig Scenes'],'FontSize',20,'FontWeight','Bold')


%PropSigScene_Scene_SSN_mean



%{
axes('Position',[0.1,0.25,0.5,0.25])

b=bar([SSI_Range(1:end-1)+diff(SSI_Range);SSI_Range(1:end-1)+diff(SSI_Range)]',Count_Spon,'Stack','LineWidth',3);
b(1).FaceColor='#676968';
b(2).FaceColor='w';
legend({'Sig','NonSig'},'FontSize',15,'FontWeight','Bold');

set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold')


xlabel('Spon Selective Index','FontSize',20,'FontWeight','Bold')
ylabel('Num of neurons','FontSize',20,'FontWeight','Bold')

box off
%}


%set(findall(gca, 'Type', 'Line'),'LineWidth',3);


%set(findall(gcf,'-property','FontSize'),'FontSize',20);
%set(findall(gcf,'-property','FontWeight'),'FontWeight','Bold');



%SSI for all Scene Selective Neurons and Non Scene Selective Neurons
%subplot(2,1,2)
axes('Position',[0.1,0.1,0.6,0.5])
shadedErrorBar(PSTH_Scene_Time_mean,Mean_PSTH_SSI_Scene_NonSSN,Sem_PSTH_SSI_Scene_NonSSN,'lineprops',{'#676968'},'transparent',1,'patchSaturation',0.4);
hold on

shadedErrorBar(PSTH_Scene_Time_mean,Mean_PSTH_SSI_Scene_SSN,Sem_PSTH_SSI_Scene_SSN,'lineprops',{'#039E20'},'transparent',1,'patchSaturation',0.5);


legend({'NonSSN','SSN'})

sig_diff=area(PSTH_Scene_Time_mean,h_diff_SSN_NonSSN*0.17,0.15,'FaceColor','#039E20','EdgeColor','#039E20','FaceAlpha',0.5,'EdgeAlpha',0.5,'ShowBaseLine','off');

yrange=ylim;
h1=plot([0,0],[yrange(1),yrange(2)],'-.k');
h2=plot([1000,1000],[yrange(1),yrange(2)],'--k');



set(get(get(sig_diff,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');

set(get(get(h1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');

set(get(get(h2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');



set(findall(gca, 'Type', 'Line'),'LineWidth',3);

set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold')


xlabel('Time from Scene Onset','FontSize',20,'FontWeight','Bold')
ylabel('SSI','FontSize',20,'FontWeight','Bold')
box off


%Comparison between sustained neurons and transcient neurons
axes('Position',[0.75,0.1,0.2,0.3])


b1=bar(Mean_SSI_Sus_Trans);
b1.FaceColor = 'flat';
b1.EdgeColor = 'flat';
b1.CData = [	251,42,42;65,128,251]/255;

hold on
errorbar(Mean_SSI_Sus_Trans,Sem_SSI_Sus_Trans,'Color','k','LineStyle','none','LineWidth',3);
yrange=ylim;


%Add significance level
 

if h_SSI_sus_trans==1
    
    %,p_RL_sus_trans
ctr2 = bsxfun(@plus, b1.XData, [b1.XOffset]');
hold on

plot(ctr2(1:2), [1 1]*yrange(2)*1.1, '-k', 'LineWidth',2)
if p_SSI_sus_trans<0.05 & p_SSI_sus_trans>=0.01
   plot(mean(ctr2(1:2)), yrange(2)*1.15, '*k','MarkerSize',10)
elseif p_SSI_sus_trans<0.01 & p_SSI_sus_trans>=0.001
    plot(mean(ctr2(1:2))*0.9, yrange(2)*1.15, '*k','MarkerSize',10)
    plot(mean(ctr2(1:2))*1.1, yrange(2)*1.15, '*k','MarkerSize',10)
        
elseif p_SSI_sus_trans<0.001
    plot(mean(ctr2(1:2))*0.8, yrange(2)*1.15, '*k','MarkerSize',10)
    plot(mean(ctr2(1:2))*1, yrange(2)*1.15, '*k','MarkerSize',10)
    plot(mean(ctr2(1:2))*1.2, yrange(2)*1.15, '*k','MarkerSize',10)
    
end


end
xticklabels({'Sus','Tran'})
ylabel('SSI')

box off

%set(gca,'LineWidth',3)


set(findall(gca, 'Type', 'Line'),'LineWidth',3);

set(gca,'LineWidth',3,'FontSize',20,'FontWeight','Bold')




end

end %End of the if showfig flag











%{
%Separate p for plotting
Num_SSIN_count=[sum(PropSigScene_Scene==0),sum(PropSigScene_Scene==1)];
PropSigScene_Scene_rest=PropSigScene_Scene(PropSigScene_Scene>0 & PropSigScene_Scene<1);

edges =[0:0.1:1];
h=histogram(PropSigScene_Scene_rest,edges);
%h=histogram(PropSigScene_Scene_rest,edges,'FaceColor','g','EdgeColor', 'g','EdgeAlpha',0.5,'FaceAlpha',0.5);
counts = h.Values;

counts_all=[Num_SSIN_count(1),counts,Num_SSIN_count(2)];
Bins_all=[0,edges(2:end)-0.05,1];


b=bar(Bins_all,counts_all,'FaceColor','b');
%}



