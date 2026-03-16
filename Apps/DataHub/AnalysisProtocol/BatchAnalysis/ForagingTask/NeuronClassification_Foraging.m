
%Summary for the event segregation analysis 
%Xuefei Yu 06072020
function NeuronClassification_Foraging(Data,  StartFile, EndFile,ShowFigureFlag,OutputFlag);
BatchFileName=Data.BatchFileName;
%Batch data file path
FilesName=Data.ResultFilePath(StartFile: EndFile);
RecordDate=Data.RecordDate(StartFile: EndFile);
FileNum=EndFile-StartFile+1;

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
     Task(i)=OutputData.EventSegregation(1).TaskCode;
     TaskName{i}=OutputData.EventSegregation(1).Task;
     
     EventSeg{i}=OutputData.EventSegregation(1).DataStamp;
     
     
end
EventStr=EventSeg{1}('SegEventStr');

%The overall scene response by event segmentation
p_Marker_Whole_temp=cellfun(@(x) x('p_Event'),EventSeg,'Uniform',0);

h_EventFromSpon_Whole=cell2mat(cellfun(@(x) x{1},p_Marker_Whole_temp,'Uniform',0))';
p_EventFromSpon_Whole=cell2mat(cellfun(@(x) x{2},p_Marker_Whole_temp,'Uniform',0))';

%Find out spon for each trial from the event segration 
AverageFR_Spon_EachTrial=cell2mat(cellfun(@(x) nanmean(x('Spon_FR')),EventSeg,'Uniform',0));

NumberOfNeuronsSigEvent=sum(h_EventFromSpon_Whole);
TotalNumberSigEvent=size(h_EventFromSpon_Whole,1);

PropOfNeuronsProbEvent=NumberOfNeuronsSigEvent/TotalNumberSigEvent;

%Only take scene and object into consideration

h_EventFromSpon_SceneObj=sum(h_EventFromSpon_Whole(:,[1:2]),2);

PropSceneObj=sum(h_EventFromSpon_SceneObj==2)/size(h_EventFromSpon_SceneObj,1);

NumberOfNeuronsSceneOrObject=NumberOfNeuronsSigEvent(1:2);
NumberOfNeuronsSceneAndObject=sum(h_EventFromSpon_SceneObj==2);
%For venn plot
Venn_Scene_Object_Data=[NumberOfNeuronsSceneOrObject-NumberOfNeuronsSceneAndObject,NumberOfNeuronsSceneAndObject];
Venn_Scene_Object_Str={sprintf('Scene\nOnly'),sprintf('Object\nOnly'),'Both'};
Offset_Venn=[1.8,1,1.8];

%Average PSTH for all the events
PSTH_events_temp=cellfun(@(x) x('PSTH_Event'),EventSeg,'Uniform',0);
h_EventFromSpon_Add=[h_EventFromSpon_Whole,ones(size(h_EventFromSpon_Whole,1),1)];

PSTH_events_time_temp=cellfun(@(x) x('PSTH_Event_Time'),EventSeg,'Uniform',0);
PSTH_events_timeMarker_temp=cellfun(@(x) x('EventTimeMarker'),EventSeg,'Uniform',0);




for i=1:length(EventStr)
PSTH_events_each_neuron=cellfun(@(x) nanmean(x{i}),PSTH_events_temp,'Uniform',0);
Numel_temp=max(cellfun(@numel,PSTH_events_each_neuron));
PSTH_events_each_neuron_mat_whole=cell2mat(cellfun(@(x) [x NaN*ones(1,Numel_temp-length(x))],PSTH_events_each_neuron,'Uniform',0)');

PSTH_events_each_neuron_mat=PSTH_events_each_neuron_mat_whole(h_EventFromSpon_Add(:,i)==1,:);
PSTH_events_each_neuron_nonresp_mat=PSTH_events_each_neuron_mat_whole(h_EventFromSpon_Add(:,i)==0,:);


PSTH_events_mean{i}=nanmean(PSTH_events_each_neuron_mat,1);
PSTH_events_sem{i}=nanstd(PSTH_events_each_neuron_mat,[],1)./sqrt(sum(~isnan(PSTH_events_each_neuron_mat),1));

PSTH_events_nonresp_mean{i}=nanmean(PSTH_events_each_neuron_nonresp_mat,1);
PSTH_events_nonresp_sem{i}=nanstd(PSTH_events_each_neuron_nonresp_mat,[],1)./sqrt(sum(~isnan(PSTH_events_each_neuron_nonresp_mat),1));

%Number
NumberOfNeuronByEvent(i,:)=[sum(h_EventFromSpon_Add(:,i)==1),sum(h_EventFromSpon_Add(:,i)==0)];

%Time
PSTH_eventstime_each_neuron=cellfun(@(x) nanmean(x{i}),PSTH_events_time_temp,'Uniform',0);
PSTH_eventstime_each_neuron_mat_whole=cell2mat(cellfun(@(x) [x NaN*ones(1,Numel_temp-length(x))],PSTH_eventstime_each_neuron,'Uniform',0)');
PSTH_eventstime_mean{i}=nanmean(PSTH_eventstime_each_neuron_mat_whole,1);
%Time marker
PSTH_eventstimeMarker_each_neuron=cell2mat(cellfun(@(x) nanmean(x{i}),PSTH_events_timeMarker_temp,'Uniform',0)');

PSTH_eventstimeMarker_each_neuron_mean(i,:)=mean(PSTH_eventstimeMarker_each_neuron);

end



if ShowFigureFlag

%Figure 1 Venn diagram for scene and object responses

figtitlestr{1}='Venn for scene and object';

fig{1}=PrepareFigure(200,'w',[70,100, 400,400],'Name',figtitlestr{1});
%First Eye trace
%subplot(2,2,3)
title('Venn diagram for scene and object responses','FontSize',15);
axis off

[H,S]=venn(NumberOfNeuronsSceneOrObject,NumberOfNeuronsSceneAndObject);


%Now label each zone 
for i = 1:3
       text(S.ZoneCentroid(i,1)*Offset_Venn(i), S.ZoneCentroid(i,2), [Venn_Scene_Object_Str(i), num2str(Venn_Scene_Object_Data(i))],'FontSize',20,'FontWeight','Bold');
end

%%
%Figure 2PSTH for each events

figtitlestr{2}='PSTH_EachEvents_Neuron_Average ';

fig{2}=PrepareFigure(201,'w',[90,100, 1000,800],'Name',figtitlestr{2});

%title(BatchFileName,'FontSize',15);
axis off
for i=1:length(EventStr)
    subplot(ceil(length(EventStr)/2),2,i);
     shadedErrorBar(PSTH_eventstime_mean{i},PSTH_events_mean{i},PSTH_events_sem{i},'lineprops',{Color(i,:)},'transparent',1);
    
    hold on
    shadedErrorBar(PSTH_eventstime_mean{i},PSTH_events_nonresp_mean{i},PSTH_events_nonresp_sem{i},'lineprops',{'#9FA4B0'},'transparent',1);
     
    xlim([min(PSTH_eventstime_mean{i}),PSTH_eventstimeMarker_each_neuron_mean(i,2)]);
    xlabel(sprintf('Time from %s onset (ms)',EventStr{i}));
    ylabel('Firing Rate(Hz)');
    set(gca,'FontSize',20,'FontWeight','Bold','LineWidth',3);
    title(EventStr{i},'FontSize',25,'FontWeight','Bold');
    set(findall(gca, 'Type', 'Line'),'LineWidth',3);
    legend(sprintf('Sig: %d',NumberOfNeuronByEvent(i,1)),sprintf('NonSig: %d',NumberOfNeuronByEvent(i,2)));
end


annotation('textbox', [0.4, 0.9, 0.8, 0.08], 'string', sprintf('%s, N=%d',BatchFileName,FileNum),'FontSize',20,'FontWeight','Bold','LineStyle','none')


end


end
%{

%{
%HP
EventSepDataHP=[


1	0	0	0.075805	1	0	1	0	0	0.462636	11.595628
1	0	0	0.4804	1	0.007998	0	0.953804	0	0.051679	25.50365
0	0.481125	0	0.937296	0	0.793005	0	0.873338	0	0.682653	15.513514
0	0.744567	0	0.502797	1	0.023171	0	0.537506	0	0.334616	21.776224
0	0.357024	0	0.146532	0	0.16202	0	0.388916	0	0.53625	3.8
1	0.000025	0	0.191581	1	0.007129	1	0.000057	1	0	40.808219
1	0	1	0	1	0	1	0.002141	1	0	5.929412
1	0.000002	1	0	0	0.196908	0	0.628524	0	0.074735	40.260274
0	0.910086	0	0.138518	1	0.015448	1	0.005286	0	0.922161	8.825175
0	0.387613	1	0.004588	1	0.003816	0	0.466756	1	0.021398	17.594406
1	0	1	0	1	0	1	0	1	0	29.733945
1	0	1	0	1	0	1	0	1	0	3.736111
1	0.005429	1	0	1	0	1	0	1	0.036186	23.391608
0	0.656919	1	0	1	0	1	0.010457	0	0.642423	12.018519
0	0.620518	1	0.000758	1	0.011549	1	0.027863	0	0.058253	8.805556
1	0.007066	0	0.863632	0	0.881559	0	0.16354	1	0.008271	8
1	0.005211	1	0	1	0	1	0	1	0	15.092466
0	0.345221	0	0.579197	0	0.922332	0	0.135203	0	0.11096	16.082645
0	0.751858	0	0.253921	0	0.727284	0	0.696739	0	0.259428	14.076923
0	0.476855	0	0.118343	0	0.08949	1	0.043321	0	0.573166	34.844444
0	0.097308	0	0.387132	0	0.282213	0	0.596654	0	0.433321	11.85
1	0	1	0.000002	1	0	1	0.000022	0	0.775826	14.379167
1	0	1	0	1	0	1	0	1	0	5.811912
0	0.800716	1	0.000147	1	0.000114	1	0.00004	1	0.000655	13.227273
1	0	1	0	1	0.000871	1	0	1	0	10.186885
0	0.195924	1	0.000067	1	0.000753	0	0.571114	0	0.854772	27.121212
1	0.001046	1	0.042023	1	0.002839	0	0.362047	0	0.171416	10.352941
1	0	1	0	1	0	1	0	1	0	19.539749
1	0	1	0	1	0	0	0.278888	1	0	19.337963
1	0.000478	0	0.475791	0	0.830911	1	0.002383	1	0.000396	24.253521
1	0	1	0.001135	1	0.049023	1	0.000758	1	0.009442	5.868056
1	0.020127	1	0	1	0.046104	1	0.000639	1	0	37.652778
1	0	1	0	1	0.002099	1	0.000648	1	0.000186	29.464789
1	0	1	0	1	0	1	0	1	0.000007	14.180556
0	0.718873	0	0.259905	0	0.684782	0	0.066975	0	0.058075	0.84375
1	0	1	0	1	0	1	0.044809	1	0.049594	13.722222
0	0.134069	1	0	1	0.000263	1	0	1	0.008161	15.833333
0	0.149492	0	0.338706	0	0.971915	1	0.000066	0	0.12737	16.0625
1	0	0	0.056675	0	0.257867	1	0	0	0.060097	11.604167
0	0.842866	0	0.562104	0	0.490982	0	0.114487	0	0.804594	2.1875
1	0.01061	1	0.001353	1	0	1	0.006189	0	0.429062	42.984127
1	0	1	0	1	0	0	0.648796	1	0	5.664671
0	0.417307	0	0.185829	0	0.07613	0	0.272284	0	0.628297	6.5
1	0.012169	0	0.262141	1	0.04381	0	0.232107	0	0.552054	20.806452
1	0.037315	0	0.397168	1	0.043182	1	0.025401	0	0.302631	16.347561
1	0.000688	1	0	1	0	0	0.276402	0	0.129313	24.18239
0	0.094575	1	0	1	0	1	0.000003	1	0.001848	5.880383
0	0.203294	0	0.343236	0	0.800098	0	0.708805	0	0.25133	21.226891
0	0.665165	0	0.8221	0	0.392661	0	0.703849	0	0.846489	3.357895
0	0.062111	1	0.007205	0	0.743984	0	0.127707	0	0.520874	25.806723
1	0.00773	1	0.020555	0	0.912253	0	0.444259	1	0.045819	6.136364
1	0.016084	0	0.945221	0	0.973511	0	0.062577	0	0.213358	8.12
0	0.861987	0	0.057679	0	0.424156	0	0.383626	0	0.89127	46.861111
1	0.007111	1	0.000432	1	0.003715	1	0.046339	1	0.003271	12.327801
1	0.017906	1	0.002501	1	0.000001	1	0.001286	0	0.163817	1.551867
];

%Caudate
EventSepDataCD=[
0	0.335318	0	0.371913	0	0.414112	1	0.003651	0	0.074222	12.64
0	0.290414	1	0.000629	1	0.047144	1	0.049757	0	0.082114	17.236111
1	0.03486	0	0.195612	0	0.165193	0	0.118527	1	0.017395	2.666667
1	0	1	0	1	0	1	0.014995	1	0.000003	0.339744
1	0	1	0	1	0	0	0.082772	1	0.000012	2.116129
];
%Putamen
EventSepDataPUT=[
0	0.233067	0	0.152353	0	0.577081	0	0.638698	0	0.340881	15.458333
1	0.000445	1	0.000028	0	0.169543	0	0.054778	0	0.937622	7.398601
];

%LGN
EventSepDataLGN=[
0	0.171404	0	0.177354	0	0.076123	0	0.19672	0	0.227872	3.75641
1	0	1	0	1	0	1	0	1	0	68.352113
];


%Summary

Sum=[32	3	1	1
31	3	1	1
35	3	0	1
30	3	0	1
24	3	0	1
16.52	6.9997	11.4285	36.0543
1.52	3.3435	4.0299	32.2979
55	5	2	2];
%}

%{
EventSepData=EventSepDataHP;

summary=sum(EventSepData(:,1:2:end-1));
Spon_Ave=mean(EventSepData(:,end));
Spon_sem=std(EventSepData(:,end))/sqrt(length(EventSepData(:,end)));
output=[summary,Spon_Ave,Spon_sem,size(EventSepData,1)]'
%}

Color=[lines(6);prism(6);pink(6)];
%Color(2,:)=copper(1)+0.5;
Color=[Color;hsv(10)];
Color(18,:)=[0.1,0.457,0.194];

ColorForEachEvent=Color;

SegEventStr={'SceneOn','FixOn','TgtOn','RewardOn','SceneOff','Spon'};
AreaString={'HPC','CDt','PUT','LGN'};


if ishandle(100)
    close 100;end
figure(100);set(100,'color','w');
set(100,'Position', [100,100, 400,600], 'Name', 'Neuron classification');

%Plot the summary figure
for i=1:size(Sum,1)-3
    Sum_Prop=Sum(i,:)./Sum(end,:);
    subplot(size(Sum,1)-2,1,i);
    b=bar(Sum_Prop,'FaceColor',ColorForEachEvent(i,:),'EdgeColor',ColorForEachEvent(i,:));
    b.FaceAlpha=0.5;
    ylim([0,1]);
    box off;
    ylabel([SegEventStr{i},newline,'(Prop)']);
    set(gca,'LineWidth',3)
    set(gca,'FontSize',15,'FontWeight','Bold')
    
        set(gca,'XTick',[])
        %{
    else
        set(gca,'XTickLabels',AreaString)
        h=gca; h.XAxis.TickLength = [0 0];
       
        
    end
     %}
    
    
end
%Add spon 

subplot(size(Sum,1)-2,1,i+1);
bs=bar(Sum(end-2,:),'FaceColor',ColorForEachEvent(i+1,:),'EdgeColor',ColorForEachEvent(i+1,:));
bs.FaceAlpha=0.5;
hold on
errorbar(Sum(end-2,:),Sum(end-1,:),'Color',ColorForEachEvent(i+1,:),'Linestyle','none','LineWidth',3)
box off
ylabel(['Spon',newline,'(Hz)'],'FontSize',15,'FontWeight','Bold')
set(gca,'XTickLabels',AreaString,'LineWidth',3,'FontSize',15,'FontWeight','Bold')
h=gca; h.XAxis.TickLength = [0 0];
%}




