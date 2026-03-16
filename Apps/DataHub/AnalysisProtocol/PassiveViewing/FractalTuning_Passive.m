function FractalTuning_Passive(Data, StartTrial, EndTrial,ShowFigureFlag,OutputFlag);
%%%%%%%%%%%%%%%%%User controller starts here
%Figure output path
%OutPutPath = 'H:\ProgramesInOHLab\DataAnalysis\FiguresOutput';

ECode;


SegEvent=[TGTCD];%Event Flow of interest
%SegEventEnd=[SHOWFIXCD,FIXOFF,TGTOFF,RWDOFFCD ,NaN,EXTRAFP];%Event end
SegEventEnd=[TGTOFF];%Event end
SegEventStr={'TargetOn'};



PlotStartOffset=[200];%in ms
PlotStopOffset=[300];%in ms

ResponseLatencyThreshold=3;

%Parameter to define the sliding window to measure firing rates
StepSize=5;%in ms
BinWidth=15;%in ms
%Color Definition for each scenes

Color=[lines(6);prism(6);pink(6)];
%Color(2,:)=copper(1)+0.5;
Color=[Color;hsv(10)];
Color(18,:)=[0.1,0.457,0.194];
Color=[Color;Color*0.5];
Color=[Color;Color*0.8];
%Color Definition for each fractals
%{
Color=[lines(6);prism(6) ];
Color(2,:)=copper(1)+0.5;
Color=[Color;Color*0.5];
Color=[Color;Color*0.8];
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Unpack all the data base and task settings
ECode;
FileName=Data.FileName;
TaskType=Data.TaskType;
TaskCode=Data.TaskCode;
TrialNum=EndTrial-StartTrial+1;
Selection=StartTrial:EndTrial;

%FixationPointIndex=Data.FixPointIndex(Selection,:);
%FixLoc=Data.FixLoc(Selection,:,:);

try
TargetIndex=Data.TargetIndex(Selection,2:end,1);
TargetName=Data.TargetIndex(Selection,2:end,2);
catch
    TargetIndex=Data.TargetIndex(Selection,2:end);
    TargetName=Data.TargetIndex(Selection,2:end);
end
%TargetName=TargetIndex;%Temporally
%TargetLoc=Data.TargetLoc(Selection,:,:);

TargetNum=Data.TargetNum(Selection,2:end);

TargetLoc=Data.TargetLoc(Selection,2:end,:);
TargetPolarLoc=Data.TargetPolarLoc(Selection,2:end,:);
TargetAngle=TargetPolarLoc(:,:,1);
TargetEcc=TargetPolarLoc(:,:,2);



%GetEyeData
EyeData=Data.EyeDataRaw;
EyeChannel_X_L=EyeData(Selection,1);
EyeChannel_Y_L=EyeData(Selection,2);
%EyeChannel_X_R=EyeData(Selection,3);
%EyeChannel_Y_R=EyeData(Selection,4);

Eye_Eccentricity=cellfun(@(x,y) sqrt(x.^2+y.^2),EyeChannel_X_L,EyeChannel_Y_L,'uniform',0);

EyeBinWidth=Data.EyeBinWidth;


EventChannel=Data.EventChannel(Selection,:);
EventTimeChannel=Data.EventTimeChannel(Selection,:);
%EventBin=Data.EventBin(Selection,:);

SpikeChannel=Data.SpikeChannel;
SpikeTime=squeeze(Data.SpikeTimeData(:,:,Selection))';
%SpikeBin=squeeze(Data.SpikeBinData(:,:,Selection))';


%Segregate spike data from fractal onset to fractal offset

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

%PSTH for each event
MeanPSTH=SegmentData.MeanPSTH;%=[PSTH_Event_Time',PSTH_Event'];
MeanPSTH_Time=MeanPSTH{1};
MeanPSTH_FR=MeanPSTH{2};

%Get the average firing rate 200ms before object onset as the baseline
SponEvent=[TGTCD];
PreTimeBuffer=200;
SponEventEnd=NaN;
StartOffset=0;
StopOffset=PreTimeBuffer;

EventTimeChannelForSpon=EventTimeChannel-PreTimeBuffer;


SegmentDataSpon=EventSegTool(SponEvent,SponEventEnd,StartOffset,StopOffset,EventChannel,EventTimeChannelForSpon,SpikeTime,StepSize,BinWidth);

Mean_Spon=SegmentDataSpon.MeanFR;
Mean_Spon_ForEachTrial=SegmentDataSpon.MeanFRRaw{1};

%Get eye trance distance during the interval
StartTimeBin=FindOutTime(EventChannel,EventTimeChannel,SegEvent)-PlotStartOffset;
EndTimeBin=FindOutTime(EventChannel,EventTimeChannel,SegEventEnd )+PlotStopOffset;
  
Eye_Eccentricity_Tgt=SelectEyeInterval(Eye_Eccentricity,StartTimeBin,EndTimeBin)';
Eye_Eccentricity_Tgt=ReorganizeEye(Eye_Eccentricity_Tgt);





%Rearrange target into lines
TargetIndex_Line=TargetIndex;

TargetIndex_Line=reshape(TargetIndex_Line,numel(TargetIndex_Line),1);

NoNNaNSelector=~isnan(TargetIndex_Line);

TargetIndex_Line=TargetIndex_Line(~isnan(TargetIndex_Line));


TargetName_Line=TargetName;

TargetName_Line=reshape(TargetName_Line,numel(TargetName_Line),1);

TargetName_Line=TargetName_Line(~isnan(TargetName_Line));

%Rearrange target locations into lines
TargetAngle_Line=TargetAngle;
TargetAngle_Line(TargetAngle_Line==360)=0;
TargetAngle_Line=reshape(TargetAngle_Line,numel(TargetAngle_Line),1);
TargetAngle_Line=TargetAngle_Line(~isnan(TargetAngle_Line));

TargetEcc_Line=TargetEcc;
TargetEcc_Line=reshape(TargetEcc_Line,numel(TargetEcc_Line),1);
TargetEcc_Line=TargetEcc_Line(~isnan(TargetEcc_Line));


%Select eye trace
Eye_Eccentricity_Tgt=Eye_Eccentricity_Tgt(NoNNaNSelector,:);

%Analysis Good versus Bad versus novel
ObjType=ObjectDefine_Passive(TargetName_Line);
GoodObjectIndex=ObjType==1;
BadObjectIndex=ObjType==2;
NovelObjectIndex=ObjType==3;


%Eye eccentricity Good v.s. Bad
Eye_Eccentricity_Tgt_Good=Eye_Eccentricity_Tgt(GoodObjectIndex,:);
Eye_Eccentricity_Tgt_Bad=Eye_Eccentricity_Tgt(BadObjectIndex,:);
Eye_Eccentricity_Tgt_Novel=Eye_Eccentricity_Tgt(NovelObjectIndex,:);


Eye_Ecc_Good_mean=nanmean(Eye_Eccentricity_Tgt_Good,1);
Eye_Ecc_Bad_mean=nanmean(Eye_Eccentricity_Tgt_Bad,1);
Eye_Ecc_Novel_mean=nanmean(Eye_Eccentricity_Tgt_Novel,1);


Eye_Ecc_Good_sem=nanstd(Eye_Eccentricity_Tgt_Good)./sqrt(sum(~isnan(Eye_Eccentricity_Tgt_Good),1));
Eye_Ecc_Bad_sem=nanstd(Eye_Eccentricity_Tgt_Bad)./sqrt(sum(~isnan(Eye_Eccentricity_Tgt_Bad),1));
Eye_Ecc_Novel_sem=nanstd(Eye_Eccentricity_Tgt_Novel)./sqrt(sum(~isnan(Eye_Eccentricity_Tgt_Novel),1));



Eye_Ecc_All_mean=nanmean(Eye_Eccentricity_Tgt);
Eye_Ecc_All_sem=nanstd(Eye_Eccentricity_Tgt)./sqrt(sum(~isnan(Eye_Eccentricity_Tgt),1));



%Overall average
MeanFR_Good=MeanFRRaw(GoodObjectIndex);
MeanFR_Bad=MeanFRRaw(BadObjectIndex);
MeanFR_Novel=MeanFRRaw(NovelObjectIndex);

MeanFR_Good_mean=nanmean(MeanFR_Good);
MeanFR_Bad_mean=nanmean(MeanFR_Bad);
MeanFR_Novel_mean=nanmean(MeanFR_Novel);

MeanFR_Good_sem=nanstd(MeanFR_Good)/sqrt(length(MeanFR_Good));
MeanFR_Bad_sem=nanstd(MeanFR_Bad)/sqrt(length(MeanFR_Bad));
MeanFR_Novel_sem=nanstd(MeanFR_Novel)/sqrt(length(MeanFR_Novel));


[h_GoodvsBad,p_GoodvsBad]=ttest2(MeanFR_Good,MeanFR_Bad);
[h_NovelvsFamiliar,p_NovelvsFamiliar]=ttest2([MeanFR_Good;MeanFR_Bad],MeanFR_Novel);

%Spon for good and bad
Mean_Spon_Good=Mean_Spon_ForEachTrial(GoodObjectIndex);
Mean_Spon_Bad=Mean_Spon_ForEachTrial(BadObjectIndex);
Mean_Spon_Novel=Mean_Spon_ForEachTrial(NovelObjectIndex);

Mean_Spon_Good_mean=nanmean(Mean_Spon_ForEachTrial(GoodObjectIndex));
Mean_Spon_Bad_mean=nanmean(Mean_Spon_ForEachTrial(BadObjectIndex));
Mean_Spon_Novel_mean=nanmean(Mean_Spon_ForEachTrial(NovelObjectIndex));

Mean_Spon_Good_sem=nanstd(Mean_Spon_ForEachTrial(GoodObjectIndex))/sqrt(sum(GoodObjectIndex));
Mean_Spon_Bad_sem=nanstd(Mean_Spon_ForEachTrial(BadObjectIndex))/sqrt(sum(BadObjectIndex));
Mean_Spon_Novel_sem=nanstd(Mean_Spon_ForEachTrial(NovelObjectIndex))/sqrt(sum(NovelObjectIndex));



%Comparison between Good/Bad and Spon firing rate
[h_GoodvsSpon,p_GoodvsSpon]=ttest(MeanFR_Good,Mean_Spon_Good);
[h_BadvsSpon,p_BadvsSpon]=ttest(MeanFR_Bad,Mean_Spon_Bad);
[h_NovelvsSpon,p_NovelvsSpon]=ttest(MeanFR_Novel,Mean_Spon_Novel);

[h_AllvsSpon,p_AllvsSpon]=ttest(MeanFRRaw,Mean_Spon_ForEachTrial);

%Comparison between Good/Bad and Spon Eye movement
Spon_Eyes=Eye_Eccentricity_Tgt(:,1:PlotStartOffset);
Spon_Eyes_reshape=reshape(Spon_Eyes,1,numel(Spon_Eyes));

[p_AllvsSponEye,h_AllvsSponEye]=ttest2_group(Eye_Eccentricity_Tgt,Spon_Eyes_reshape','both');
[p_GoodvsBadEye,h_GoodvsBadEye]=ttest2_group(Eye_Eccentricity_Tgt_Good,Eye_Eccentricity_Tgt_Bad,'both');


%PSTH
PSTH_Time_Good=MeanPSTH_Time(GoodObjectIndex,:);
PSTH_FR_Good=MeanPSTH_FR(GoodObjectIndex,:);

PSTH_Time_Bad=MeanPSTH_Time(BadObjectIndex,:);
PSTH_FR_Bad=MeanPSTH_FR(BadObjectIndex,:);

PSTH_Time_Novel=MeanPSTH_Time(NovelObjectIndex,:);
PSTH_FR_Novel=MeanPSTH_FR(NovelObjectIndex,:);

PSTH_Time_All=nanmean(MeanPSTH_Time,1);
PSTH_FR_All_mean=nanmean(MeanPSTH_FR,1);
PSTH_FR_All_sem=nanstd(MeanPSTH_FR,1)./sqrt(sum(~isnan(MeanPSTH_Time),1));

PSTH_FR_Good_mean=nanmean(PSTH_FR_Good,1);
PSTH_FR_Bad_mean=nanmean(PSTH_FR_Bad,1);
PSTH_FR_Novel_mean=nanmean(PSTH_FR_Novel,1);


PSTH_FR_Good_sem=nanstd(PSTH_FR_Good,[],1)./sqrt(sum(~isnan(PSTH_FR_Good),1));
PSTH_FR_Bad_sem=nanstd(PSTH_FR_Bad,[],1)./sqrt(sum(~isnan(PSTH_FR_Bad),1));
PSTH_FR_Novel_sem=nanstd(PSTH_FR_Novel,[],1)./sqrt(sum(~isnan(PSTH_FR_Novel),1));

%PSTH p value from the spon
%[p_PSTH_Good,h_PSTH_Good]=ttest2_group(PSTH_FR_Good,Mean_Spon_ForEachTrial,'both');
[p_PSTH_All,h_PSTH_All]=ttest2_group(MeanPSTH_FR,Mean_Spon_ForEachTrial,'both');
[p_PSTH_GoodvsBad,h_PSTH_GoodvsBad]=ttest2_group(PSTH_FR_Good,PSTH_FR_Bad,'both');
[p_PSTH_NovelvsFamiliar,h_PSTH_NovelvsFamiliar]=ttest2_group(PSTH_FR_Novel,[PSTH_FR_Good;PSTH_FR_Bad],'both');

 

UniqueTargetName=unique(TargetName_Line);
NumberOfObject=length(UniqueTargetName);

%%
%%%%%%%%%%%%%Loop over the objects
for i=1:NumberOfObject
    currObject=UniqueTargetName(i);
    FR_Whole_EachObject=MeanFRRaw(TargetName_Line==currObject);
    Mean_FR_Whole_EachObject(i)=nanmean( FR_Whole_EachObject);
    SEM_FR_Whole_EachObject(i)=nanstd( FR_Whole_EachObject)/sqrt(length(FR_Whole_EachObject));
    %Spon for each object
     Spon_ForEachObject=Mean_Spon_ForEachTrial(TargetName_Line==currObject);
    %p value from the spon
    [h_ObjectFromSpon(i),p_ObjectFromSpon(i)]=ttest(FR_Whole_EachObject,Spon_ForEachObject);
    
    PSTH_FR_EachObject=MeanPSTH_FR(TargetName_Line==currObject,:);
    
    
    Mean_PSTH_FR_EachObject(i,:)=nanmean(PSTH_FR_EachObject,1);
    SEM_PSTH_FR_EachObject(i,:)=nanstd(PSTH_FR_EachObject,[],1)./sqrt(sum(~isnan(PSTH_FR_EachObject)));
    
  
    %Eye eccentricity
    Eye_Eccentricity_Tgt_ForObject=Eye_Eccentricity_Tgt(TargetName_Line==currObject,:);
    Mean_Eye_Eccentricity_ForObject(i,:)=nanmean(Eye_Eccentricity_Tgt_ForObject,1);
    SEM_Eye_Eccentricity_ForObject(i,:)=nanstd(Eye_Eccentricity_Tgt_ForObject,[],1)./sqrt(sum(~isnan(Eye_Eccentricity_Tgt_ForObject)));
    

 
    FractalStr{i}=sprintf('Fractal%d',currObject);
    
    %PSTH p value from the spon
    [p_EachFractal(i,:),h_EachScene(i,:)]=ttest2_group(PSTH_FR_EachObject,Mean_Spon_ForEachTrial,'right');
end


%For unique objects
ObjectType_Unique=ObjectDefine_Passive(UniqueTargetName);
GoodObject=ObjectType_Unique==1;
BadObject=ObjectType_Unique==2;
NovelObject=ObjectType_Unique==3;

%p Value by ANOVA
[pFractal_SSI_All,tbl,stats]=anova1(MeanFRRaw,TargetName_Line,'off');
[pFractal_SSI_Good,tbl,stats]=anova1(MeanFR_Good,TargetName_Line(GoodObjectIndex),'off');
[pFractal_SSI_Bad,tbl,stats]=anova1(MeanFR_Bad,TargetName_Line(BadObjectIndex),'off');
[pFractal_SSI_Novel,tbl,stats]=anova1(MeanFR_Novel,TargetName_Line(NovelObjectIndex),'off');

%%
%%%%%  Location Selectivity%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

UniqueTargetAngle=unique(TargetAngle_Line);
UniqueTargetEcc=unique(TargetEcc_Line);

%Location Tuning
for ii=1:length(UniqueTargetAngle)
    for jj=1:length(UniqueTargetEcc)
        Sel=(TargetAngle_Line==UniqueTargetAngle(ii)) &(TargetEcc_Line==UniqueTargetEcc(jj));
        FR_Loc=MeanFRRaw(Sel);
        FR_Loc_mean(ii,jj)=nanmean(FR_Loc);
        FR_Loc_sem(ii,jj)=nanstd(FR_Loc)/sqrt(sum(Sel));
        
        FR_Spon_Loc=Mean_Spon_ForEachTrial(Sel);
        FR_Spon_Loc_mean(ii,jj)=nanmean(FR_Spon_Loc);
        FR_Spon_Sem_mean(ii,jj)=nanstd(FR_Spon_Loc)/sqrt(sum(Sel));
        
        
    end
end


%%

if ShowFigureFlag
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Figure 1: Average firing rate
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%for each fractal
figtitlestr{1}='Fractal_Tuning';

fig{1}=PrepareFigure(104,'w',[50,100, 1200,800],'Name',figtitlestr{1});

%subplot(2,1,1)

axes('position',[0.1,0.7,0.6,0.2])


g=bar([Mean_FR_Whole_EachObject(GoodObject)';Mean_Spon_Good_mean],'LineWidth',3);  
g.FaceColor = 'flat';
%Color=rand(NumOfCondition,3)
g.CData = [Color(1:sum(GoodObject),:);0,0,0];
%ylim([0,max(Mean_FR_Whole_EachObject+SEM_FR_Whole_EachObject)]); 

hold on

errorbar([Mean_FR_Whole_EachObject(GoodObject),Mean_Spon_Good_mean],[SEM_FR_Whole_EachObject(GoodObject),Mean_Spon_Good_sem],'.k','LineWidth',3);
set(gca,'LineWidth',3)
xticks(1:sum(GoodObject)+1);
xticklabels([string(UniqueTargetName(GoodObject));'Spon_G']);
set(gca,'FontSize',15,'FontWeight','Bold')
%xlabel('Fractal Name')
ylabel('Spikes/s')

box off
text(sum(GoodObject)/10,max(Mean_FR_Whole_EachObject+SEM_FR_Whole_EachObject)*1.15,'Good Object','FontSize',15,'FontWeight','Bold')

title_msg=sprintf(' %s, SpikeChannel: %d, TrialNum: %d ',FileName,SpikeChannel,TrialNum);
text(NumberOfObject/5,max(Mean_FR_Whole_EachObject+SEM_FR_Whole_EachObject)*1.4,title_msg,'FontSize',25,'FontWeight','Bold');
%Add significance level

XData=g.XData;
yData=g.YData;

XData=XData(1:end-1);
yData=yData(1:end-1);


text(XData(h_ObjectFromSpon(GoodObject)==1)-0.045,yData(h_ObjectFromSpon(GoodObject)==1)+max(Mean_FR_Whole_EachObject+SEM_FR_Whole_EachObject)*0.1,'*','FontSize',40)
text(XData(p_ObjectFromSpon(GoodObject)<0.01)-0.045,yData(p_ObjectFromSpon(GoodObject)<0.01)+max(Mean_FR_Whole_EachObject+SEM_FR_Whole_EachObject)*0.2,'*','FontSize',40)
text(XData(p_ObjectFromSpon(GoodObject)<0.001)-0.045,yData(p_ObjectFromSpon(GoodObject)<0.001)+max(Mean_FR_Whole_EachObject+SEM_FR_Whole_EachObject)*0.3,'*','FontSize',40)


%ylim([0,max(Mean_FR_Whole_EachObject+SEM_FR_Whole_EachObject)*1.35]); 




axes('position',[0.1,0.45,0.6,0.2])
%axes('position',[0.48,0.7,0.35,0.23])
%ann=annotation('textbox',[.1 .87 .5 .1],'String',title_msg,'FitBoxToText','on');
%ann.FontSize=13;
%subplot(2,1,2)
g=bar([Mean_FR_Whole_EachObject(BadObject)';Mean_Spon_Bad_mean],'LineWidth',3);  
g.FaceColor = 'flat';
%Color=rand(NumOfCondition,3)
g.CData = [Color(1:sum(BadObject),:);0,0,0];
   

hold on

errorbar([Mean_FR_Whole_EachObject(BadObject),Mean_Spon_Bad_mean],[SEM_FR_Whole_EachObject(BadObject),Mean_Spon_Bad_sem],'.k','LineWidth',3);
set(gca,'LineWidth',3)
xticks(1:sum(BadObject)+1);
xticklabels([string(UniqueTargetName(BadObject));'Spon_B']);
set(gca,'FontSize',15,'FontWeight','Bold')
%xlabel('Fractal Name')
ylabel('Spikes/s')

box off



text(sum(BadObject)/10,max(Mean_FR_Whole_EachObject+SEM_FR_Whole_EachObject)*1.15,'Bad Object','FontSize',15,'FontWeight','Bold')

%Add significance level

XData=g.XData;
yData=g.YData;

XData=XData(1:end-1);
yData=yData(1:end-1);


text(XData(h_ObjectFromSpon(BadObject)==1)-0.045,yData(h_ObjectFromSpon(BadObject)==1)+max(Mean_FR_Whole_EachObject+SEM_FR_Whole_EachObject)*0.1,'*','FontSize',40)
text(XData(p_ObjectFromSpon(BadObject)<0.01)-0.045,yData(p_ObjectFromSpon(BadObject)<0.01)+max(Mean_FR_Whole_EachObject+SEM_FR_Whole_EachObject)*0.2,'*','FontSize',40)
text(XData(p_ObjectFromSpon(BadObject)<0.001)-0.045,yData(p_ObjectFromSpon(BadObject)<0.001)+max(Mean_FR_Whole_EachObject+SEM_FR_Whole_EachObject)*0.3,'*','FontSize',40)


%ylim([0,max(Mean_FR_Whole_EachObject+SEM_FR_Whole_EachObject)*1.35]); 
%Novel Object
axes('position',[0.1,0.2,0.6,0.2])
g=bar([Mean_FR_Whole_EachObject(NovelObject)';Mean_Spon_Novel_mean],'LineWidth',3);  
g.FaceColor = 'flat';
%Color=rand(NumOfCondition,3)
g.CData = [Color(1:sum(NovelObject),:);0,0,0];
   

hold on

errorbar([Mean_FR_Whole_EachObject(NovelObject),Mean_Spon_Novel_mean],[SEM_FR_Whole_EachObject(NovelObject),Mean_Spon_Novel_sem],'.k','LineWidth',3);
set(gca,'LineWidth',3)
xticks(1:sum(NovelObject)+1);
xticklabels([string(UniqueTargetName(NovelObject));'Spon_N']);
set(gca,'FontSize',15,'FontWeight','Bold')
xlabel('Fractal Name')
ylabel('Spikes/s')

box off



text(sum(NovelObject)/10,max(Mean_FR_Whole_EachObject+SEM_FR_Whole_EachObject)*1.15,'Novel Object','FontSize',15,'FontWeight','Bold')

%Add significance level

XData=g.XData;
yData=g.YData;

XData=XData(1:end-1);
yData=yData(1:end-1);


text(XData(h_ObjectFromSpon(NovelObject)==1)-0.045,yData(h_ObjectFromSpon(NovelObject)==1)+max(Mean_FR_Whole_EachObject+SEM_FR_Whole_EachObject)*0.1,'*','FontSize',40)
text(XData(p_ObjectFromSpon(NovelObject)<0.01)-0.045,yData(p_ObjectFromSpon(NovelObject)<0.01)+max(Mean_FR_Whole_EachObject+SEM_FR_Whole_EachObject)*0.2,'*','FontSize',40)
text(XData(p_ObjectFromSpon(NovelObject)<0.001)-0.045,yData(p_ObjectFromSpon(NovelObject)<0.001)+max(Mean_FR_Whole_EachObject+SEM_FR_Whole_EachObject)*0.3,'*','FontSize',40)



%%Good vs Bad vs Novel
axes('position',[0.75,0.2,0.25,0.3]);
g=bar([MeanFR_Good_mean,MeanFR_Bad_mean,MeanFR_Novel_mean,MeanFR(1),Mean_Spon(1)],'LineWidth',3);  
g.FaceColor = 'flat';
%Color=rand(NumOfCondition,3)
g.CData = [1,0,0;0,0,1;0.7176,0.9255,0.5529;0.7,0.4,0.15;0,0,0];
%ylim([0,max([MeanFR_Good_mean+MeanFR_Good_sem,MeanFR_Bad_mean+MeanFR_Bad_sem])]); 

hold on

errorbar([MeanFR_Good_mean,MeanFR_Bad_mean,MeanFR_Novel_mean,MeanFR(1),Mean_Spon(1)],[MeanFR_Good_sem,MeanFR_Bad_sem,MeanFR_Novel_sem,MeanFR(2),Mean_Spon(2)],'.k','LineWidth',3);
set(gca,'LineWidth',3)
xticks(1:5);
xticklabels({'Good','Bad','Novel','All','Spon'});
set(gca,'FontSize',15,'FontWeight','Bold')
xlabel('Fractal Type')
ylabel('Spikes/s')
title('Good vs Bad vs Novel');

box off

%Add significance level
XData=g.XData;
yData=g.YData;

XData=XData(1:end-1);
yData=yData(1:end-1);


h_All=[h_GoodvsSpon,h_BadvsSpon,h_NovelvsSpon,h_AllvsSpon,NaN];
p_All=[p_GoodvsSpon,p_BadvsSpon,p_NovelvsSpon,p_AllvsSpon,NaN];
yrange=ylim;
text(XData(h_All==1)-0.045,yData(h_All==1)+yrange(2)*0.1,'*','FontSize',40)
text(XData(p_All<0.01)-0.045,yData(p_All<0.01)+yrange(2)*0.2,'*','FontSize',40)
text(XData(p_All<0.001)-0.045,yData(p_All<0.001)+yrange(2)*0.3,'*','FontSize',40)


%Add significant level between good and bad


if h_GoodvsBad==1
ctr2 = bsxfun(@plus, g.XData, [g.XOffset]');
hold on
plot(ctr2(1:2), [1 1]*yrange(2)*1.35, '-k', 'LineWidth',3,'Color','r')
if p_GoodvsBad<0.05 && p_GoodvsBad>=0.01
   plot(mean(ctr2(1:2)), yrange(2)*1.4, '*k','MarkerSize',15,'Color','r','LineWidth',3);
elseif p_GoodvsBad<0.01 && p_GoodvsBad>=0.001
   plot(mean(ctr2(1:2))*0.9, yrange(2)*1.4, '*k','MarkerSize',15,'Color','r','LineWidth',3);
   plot(mean(ctr2(1:2))*1.1, yrange(2)*1.4, '*k','MarkerSize',15,'Color','r','LineWidth',3);
elseif p_GoodvsBad<0.001 
   plot(mean(ctr2(1:2))*0.9, yrange(2)*1.4, '*k','MarkerSize',15,'Color','r','LineWidth',3);
   plot(mean(ctr2(1:2))*1.1, yrange(2)*1.4, '*k','MarkerSize',15,'Color','r','LineWidth',3);
   plot(mean(ctr2(1:2))*1.3, yrange(2)*1.4, '*k','MarkerSize',15,'Color','r','LineWidth',3);
end

end

if h_NovelvsFamiliar==1
ctr2 = bsxfun(@plus, g.XData, [g.XOffset]');
hold on
%,p_NovelvsFamiliar
plot([mean(ctr2(1:2)),ctr2(3)], [1 1]*yrange(2)*1.35, '-k', 'LineWidth',3,'Color','r')
if p_NovelvsFamiliar<0.05 & p_NovelvsFamiliar>=0.01
   plot(mean(ctr2(1:3)), yrange(2)*1.4, '*k','MarkerSize',15,'Color','r','LineWidth',3);
elseif p_NovelvsFamiliar<0.01 & p_NovelvsFamiliar>=0.001
   plot(mean(ctr2(1:3))*0.9, yrange(2)*1.4, '*k','MarkerSize',15,'Color','r','LineWidth',3);
   plot(mean(ctr2(1:3))*1.1, yrange(2)*1.4, '*k','MarkerSize',15,'Color','r','LineWidth',3);
elseif p_NovelvsFamiliar<0.001 
   plot(mean(ctr2(1:3))*0.9, yrange(2)*1.4, '*k','MarkerSize',15,'Color','r','LineWidth',3);
   plot(mean(ctr2(1:3))*1.1, yrange(2)*1.4, '*k','MarkerSize',15,'Color','r','LineWidth',3);
   plot(mean(ctr2(1:3))*1.3, yrange(2)*1.4, '*k','MarkerSize',15,'Color','r','LineWidth',3);
end

end


ylim(yrange*1.45);

%%Fractal Selective Index
%axes('position',[0.75,0.65,0.2,0.2]);
%strbox_Good=sprintf('FSI-Good: %2.2f, p= %2.4f',FractalSelectiveIndex_Good,pFractal_SSI_Good);
%strbox_Bad=sprintf('FSI-Bad: %2.2f, p= %2.4f',FractalSelectiveIndex_Bad,pFractal_SSI_Bad);
%strbox_Novel=sprintf('FSI-Novel: %2.2f, p= %2.4f',FractalSelectiveIndex_Novel,pFractal_SSI_Bad);


%strbox=sprintf('FSI-All: %2.2f, p= %2.4f',FractalSelectiveIndex_All,pFractal_SSI_All);
%{
ann{1}=annotation('textbox',[0.7,0.65,0.2,0.2],'String',strbox_Good,'FitBoxToText','on');
ann{2}=annotation('textbox',[0.7,0.6,0.2,0.2],'String',strbox_Bad,'FitBoxToText','on');
ann{3}=annotation('textbox',[0.7,0.55,0.2,0.2],'String',strbox,'FitBoxToText','on');
%}
%p_temp=[pFractal_SSI_Good,pFractal_SSI_Bad,pFractal_SSI_All];
%{
for ii=1:3

ann{ii}.FontSize=25;
ann{ii}.FontWeight='Bold';
ann{ii}.LineStyle='none';

if p_temp(ii)<0.05
    ann{ii}.Color='r';
end
end

%}








%%
%%%%%%%%%%%%%%%%%%%%%%%% Figure 2 PSTH for good and
%%%%%%%%%%%%%%%%%%%%%%%% bad and novel%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figtitlestr{2}='PSTH_GoodvsBadvsNovel';

fig{2}=PrepareFigure(105,'w',[70,100, 1000,800],'Name',figtitlestr{2});
%First Eye trace
%subplot(2,2,3)
title('PSTH for all targets as a function of time','FontSize',25);
axis off

axes('Position',[0.1,0.7,0.4,0.15]);
shadedErrorBar([-PlotStartOffset:length(Eye_Ecc_Good_mean)-PlotStartOffset-1],Eye_Ecc_All_mean,Eye_Ecc_All_sem,'lineprops',{[0.9290    0.6940    0.1250]},'transparent',1,'patchSaturation',0.3);


% plot([min(TimeMarker(TimeMarker(:,2)>1,2)),min(TimeMarker(TimeMarker(:,2)>1,2))],yrange,'--k');
%xlim([-PlotStartOffset, median(TimeMarker(:,2))+PlotStopOffset]);
%xlim([-PlotStartOffset,max(min(TimeMarker(:,2))+PlotStopOffset,400)]);
xlim([-PlotStartOffset,min(TimeMarker(:,2))+PlotStopOffset]);
%xlim([-PlotStartOffset,mean(TimeMarker(:,2))+PlotStopOffset]);
%xlim([-PlotStartOffset,max(mean(TimeMarker(:,2))+PlotStopOffset,800)]);
hold on
%Add significance level
yrange=ylim;
h_AllvsSponEye(h_AllvsSponEye==0)=NaN;
sig=area(PSTH_Time_All(1):PSTH_Time_All(1)+length(h_AllvsSponEye)-1,h_AllvsSponEye*max(yrange(1),1)*0.6,0,'FaceColor','#FFC107','EdgeColor','#FFC107','FaceAlpha',0.5,'EdgeAlpha',0.5,'ShowBaseLine','off');



set(findall(gca, 'Type', 'Line'),'LineWidth',3);
ylabel(sprintf('Eye \nEccentricity'));
%ylabel('Firing Rate(Hz)')
xticks([]);
xticklabels([]);
set(gca,'LineWidth',3,'FontSize',25);






%subplot(2,2,4)
axes('Position',[0.55,0.7,0.4,0.15]);
shadedErrorBar([-PlotStartOffset:length(Eye_Ecc_Good_mean)-PlotStartOffset-1],Eye_Ecc_Good_mean,Eye_Ecc_Good_sem,'lineprops',{[1    0    0]},'transparent',1,'patchSaturation',0.3);
hold on
shadedErrorBar([-PlotStartOffset:length(Eye_Ecc_Good_mean)-PlotStartOffset-1],Eye_Ecc_Bad_mean,Eye_Ecc_Bad_sem,'lineprops',{[0    0    1]},'transparent',1,'patchSaturation',0.3);
shadedErrorBar([-PlotStartOffset:length(Eye_Ecc_Good_mean)-PlotStartOffset-1],Eye_Ecc_Novel_mean,Eye_Ecc_Novel_sem,'lineprops',{[0    1   0]},'transparent',1,'patchSaturation',0.3);

% plot([min(TimeMarker(TimeMarker(:,2)>1,2)),min(TimeMarker(TimeMarker(:,2)>1,2))],yrange,'--k');
%xlim([-PlotStartOffset,min(TimeMarker(:,2))+PlotStopOffset]);
%xlim([-PlotStartOffset, median(TimeMarker(:,2))+PlotStopOffset]);
%xlim([-PlotStartOffset,max(min(TimeMarker(:,2))+PlotStopOffset,400)]);
xlim([-PlotStartOffset,min(TimeMarker(:,2))+PlotStopOffset]);
%xlim([-PlotStartOffset,mean(TimeMarker(:,2))+PlotStopOffset]);
set(findall(gca, 'Type', 'Line'),'LineWidth',3);

%Add significance level
yrange=ylim;
h_GoodvsBadEye(h_GoodvsBadEye==0)=NaN;
sig=area(PSTH_Time_All(1):PSTH_Time_All(1)+length(h_GoodvsBadEye)-1,h_GoodvsBadEye*max(yrange(1),1)*0.8,0,'FaceColor','#FFC107','EdgeColor','#FFC107','FaceAlpha',0.5,'EdgeAlpha',0.5,'ShowBaseLine','off');







%xlabel('Time from target onset');
%ylabel('Firing Rate(Hz)')
%ylabel('Eye Eccentricity');
set(gca,'LineWidth',3,'FontSize',25);
xticks([]);
xticklabels([]);

%subplot(2,2,1)
%%%%%%%%%%%%%%%%%%%%%%PSTH for the average of all fractals
axes('Position',[0.1,0.15,0.4,0.5])
shadedErrorBar(PSTH_Time_All,PSTH_FR_All_mean,PSTH_FR_All_sem,'lineprops',{[0.9290    0.6940    0.1250]},'transparent',1,'patchSaturation',0.3);

% plot([min(TimeMarker(TimeMarker(:,2)>1,2)),min(TimeMarker(TimeMarker(:,2)>1,2))],yrange,'--k');
%xlim([-PlotStartOffset,min(TimeMarker(:,2))+PlotStopOffset]);
%xlim([-PlotStartOffset, median(TimeMarker(:,2))+PlotStopOffset]);
%xlim([-PlotStartOffset,max(min(TimeMarker(:,2))+PlotStopOffset,400)]);
xlim([-PlotStartOffset,min(TimeMarker(:,2))+PlotStopOffset]);
%xlim([-PlotStartOffset,mean(TimeMarker(:,2))+PlotStopOffset]);

set(findall(gca, 'Type', 'Line'),'LineWidth',3);
xlabel('Time from Object onset');
ylabel('Firing Rate(Hz)')

set(gca,'LineWidth',3,'FontSize',25);

yrange=ylim;
hold on
%Add significance level
h_PSTH_All(h_PSTH_All==0)=NaN;
sig=area(PSTH_Time_All,h_PSTH_All*max(yrange(1),1)*1.1,max(yrange(1),1)*1,'FaceColor','#FFC107','EdgeColor','#FFC107','FaceAlpha',0.5,'EdgeAlpha',0.5,'ShowBaseLine','off');



%subplot(2,2,2)
axes('Position',[0.55,0.15,0.4,0.5])
shadedErrorBar(PSTH_Time_All,PSTH_FR_Good_mean,PSTH_FR_Good_sem,'lineprops',{[1,0,0]},'transparent',1,'patchSaturation',0.3)
hold on
shadedErrorBar(PSTH_Time_All,PSTH_FR_Bad_mean,PSTH_FR_Bad_sem,'lineprops',{[0,0,1]},'transparent',1,'patchSaturation',0.3)
shadedErrorBar(PSTH_Time_All,PSTH_FR_Novel_mean,PSTH_FR_Novel_sem,'lineprops',{[0,1,0]},'transparent',1,'patchSaturation',0.3)


%xlim([-PlotStartOffset,min(TimeMarker(:,2))+PlotStopOffset]);
%xlim([-PlotStartOffset, median(TimeMarker(:,2))+PlotStopOffset]);
%xlim([-PlotStartOffset,max(min(TimeMarker(:,2))+PlotStopOffset,400)]);
xlim([-PlotStartOffset,min(TimeMarker(:,2))+PlotStopOffset]);
%xlim([-PlotStartOffset,mean(TimeMarker(:,2))+PlotStopOffset]);

set(findall(gca, 'Type', 'Line'),'LineWidth',3);
xlabel('Time from object onset (ms)');
%ylabel('Firing Rate(Hz)')
ldg=legend({'Good','Bad','Novel'},'Location','northwest');
set(gca,'LineWidth',3,'FontSize',25);

%Add significance level

%yrange=ylim;
%ylim(yrange*0.9);
yrange=ylim;
h_PSTH_GoodvsBad(h_PSTH_GoodvsBad==0)=NaN;
sig=area(PSTH_Time_All,h_PSTH_GoodvsBad*max(yrange(1),1)*1.2,max(yrange(1),1)*1,'FaceColor','#FFC107','EdgeColor','#FFC107','FaceAlpha',0.5,'EdgeAlpha',0.5,'ShowBaseLine','off');

set(get(get(sig,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
ylim([yrange(1),yrange(2)*1.25]);

%Novel v.s. Familiar

h_PSTH_NovelvsFamiliar(h_PSTH_NovelvsFamiliar==0)=NaN;
sig=area(PSTH_Time_All,-h_PSTH_NovelvsFamiliar*max(yrange(1),1)*1.2,min(yrange(1),1)*1,'FaceColor','#B7EC8D','EdgeColor','#B7EC8D','FaceAlpha',0.5,'EdgeAlpha',0.5,'ShowBaseLine','off');

set(get(get(sig,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
ylim([min(yrange(1))-1,yrange(2)*1.25]);

%[p_PSTH_All,h_PSTH_All]=ttest2_group(MeanPSTH_FR,Mean_Spon_ForEachTrial,'both');
%[p_PSTH_GoodvsBad,h_PSTH_GoodvsBad]=ttest2_group(PSTH_FR_Good,PSTH_FR_Bad,'both');


%%
%%%%%%%%%%%%%%%%%%%%%%%% Figure 3 & 4 PSTH for Each fractal
%%%%%%%%%%%%%%%%%%%%%%%% bad%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figtitlestr{3}='PSTH_EachFractal_Whole';
%title('PSTH for each fractals','FontSize',25,'FontWeight','Bold');

figtitlestr{4}='PSTH_EachFractal_Sep';

fig{3}=PrepareFigure(106,'w',[50,100, 1000,800],'Name',figtitlestr{3});
fig{4}=PrepareFigure(107,'w',[100,100, 1000,800],'Name',figtitlestr{4});

temp_GoodObject_mean=Mean_PSTH_FR_EachObject(GoodObject,:);
temp_BadObject_mean=Mean_PSTH_FR_EachObject(BadObject,:);
temp_NovelObject_mean=Mean_PSTH_FR_EachObject(NovelObject,:);


temp_GoodObject_sem=SEM_PSTH_FR_EachObject(GoodObject,:);
temp_BadObject_sem=SEM_PSTH_FR_EachObject(BadObject,:);
temp_NovelObject_sem=SEM_PSTH_FR_EachObject(NovelObject,:);

temp_GoodStr=cell(max([sum(GoodObject),sum(BadObject),sum(NovelObject)]),1);
temp_BadStr=temp_GoodStr;
temp_NovelStr=temp_GoodStr;

temp_GoodStr(1:sum(GoodObject))=FractalStr(GoodObject);
temp_BadStr(1:sum(BadObject))=FractalStr(BadObject);
temp_NovelStr(1:sum(NovelObject))=FractalStr(NovelObject);




for i=1:max([sum(GoodObject),sum(BadObject),sum(NovelObject)])
   
    
set(groot,'CurrentFigure',fig{3});
subplot(3,1,1);%Good Object
title('Good Object');
if i<=sum(GoodObject)
shadedErrorBar(PSTH_Time_All,temp_GoodObject_mean(i,:),temp_GoodObject_sem(i,:),'lineprops',{Color(i,:)},'transparent',1,'patchSaturation',0.3)
end
set(findall(gca, 'Type', 'Line'),'LineWidth',3);
%ylim([-0.5*min(min(Mean_PSTH_FR_EachObject)),max(max(Mean_PSTH_FR_EachObject+SEM_PSTH_FR_EachObject))*1.02]);
%ylim([-0.5*min(min(Mean_PSTH_FR_EachObject)),max(max(Mean_PSTH_FR_EachObject(ceil(PlotStartOffset):ceil(PlotStartOffset)+1000)+SEM_PSTH_FR_EachObject(ceil(PlotStartOffset):ceil(PlotStartOffset)+1000)))*1.02])

%xlim([min(min(PSTH_Time_All))-10,max(max(PSTH_Time_All))+10]);
%xlim([-PlotStartOffset,min(TimeMarker(:,2))+PlotStopOffset]);
%xlim([-PlotStartOffset, median(TimeMarker(:,2))+PlotStopOffset]);
xlim([-PlotStartOffset,max(min(TimeMarker(:,2))+PlotStopOffset,400)]);
ylabel('Spikes/s');

set(gca,'LineWidth',3,'FontSize',20,'FontWeight','Bold');
 if i==NumberOfObject/3
     legend(temp_GoodStr);
 end

hold on
subplot(3,1,2)%Bad
title('Bad Object');
if i<=sum(BadObject)
shadedErrorBar(PSTH_Time_All,temp_BadObject_mean(i,:),temp_BadObject_sem(i,:),'lineprops',{Color(i,:)},'transparent',1,'patchSaturation',0.3)
end
hold on
set(findall(gca, 'Type', 'Line'),'LineWidth',3);




%ylim([-0.5*min(min(Mean_PSTH_FR_EachObject)),max(max(Mean_PSTH_FR_EachObject(ceil(PlotStartOffset):ceil(PlotStartOffset)+1000)+SEM_PSTH_FR_EachObject(ceil(PlotStartOffset):ceil(PlotStartOffset)+1000)))*1.02])
%xlim([min(min(PSTH_Time_All))-10,max(max(PSTH_Time_All))+10]);
%xlim([-PlotStartOffset,min(TimeMarker(:,2))+PlotStopOffset]);
%xlim([-PlotStartOffset, median(TimeMarker(:,2))+PlotStopOffset]);
%xlim([-PlotStartOffset,max(min(TimeMarker(:,2))+PlotStopOffset,400)]);
xlim([-PlotStartOffset,min(TimeMarker(:,2))+PlotStopOffset]);
ylabel('Spikes/s')
xlabel('Time from fractal onset (ms)')
set(gca,'LineWidth',3,'FontSize',20,'FontWeight','Bold');
if i==NumberOfObject/3
     legend(temp_BadStr);
end

subplot(3,1,3)%Novel
title('Novel Object');
if i<=sum(NovelObject)
shadedErrorBar(PSTH_Time_All,temp_NovelObject_mean(i,:),temp_NovelObject_sem(i,:),'lineprops',{Color(i,:)},'transparent',1,'patchSaturation',0.3)
end
hold on
set(findall(gca, 'Type', 'Line'),'LineWidth',3);

xlim([-PlotStartOffset,min(TimeMarker(:,2))+PlotStopOffset]);
ylabel('Spikes/s')
xlabel('Time from fractal onset (ms)')
set(gca,'LineWidth',3,'FontSize',20,'FontWeight','Bold');
if i==NumberOfObject/3
     legend(temp_NovelStr);
end



set(groot,'CurrentFigure',fig{4});
subplot(ceil(NumberOfObject/3/4),4,i);
if i<=sum(GoodObject)
line1=shadedErrorBar(PSTH_Time_All,temp_GoodObject_mean(i,:),temp_GoodObject_sem(i,:),'lineprops',{[1,0,0]},'transparent',1,'patchSaturation',0.3);
end
hold on
if i<=sum(BadObject)
line2=shadedErrorBar(PSTH_Time_All,temp_BadObject_mean(i,:),temp_BadObject_sem(i,:),'lineprops',{[0,0,1]},'transparent',1,'patchSaturation',0.3);
end
hold on
if i<=sum(NovelObject)
line2=shadedErrorBar(PSTH_Time_All,temp_NovelObject_mean(i,:),temp_NovelObject_sem(i,:),'lineprops',{[0,1,0]},'transparent',1,'patchSaturation',0.3);
end

MaxVal=max(max(Mean_PSTH_FR_EachObject(PlotStartOffset:PlotStartOffset+100)+SEM_PSTH_FR_EachObject(PlotStartOffset:PlotStartOffset+100)))*1.02;





%ylim([-0.5*min(min(Mean_PSTH_FR_EachObject)),max(max(MaxVal+temp_GoodObject_Eye_mean+temp_GoodObject_Eye_sem))*1.02]);


%xlim([min(min(PSTH_Time_All))-10,max(max(PSTH_Time_All))+10]);
set(findall(gca, 'Type', 'Line'),'LineWidth',3);
set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');

%xlim([-PlotStartOffset,max(min(TimeMarker(:,2))+PlotStopOffset,400)]);
xlim([-PlotStartOffset,min(TimeMarker(:,2))+PlotStopOffset]);
%xlim([-PlotStartOffset, median(TimeMarker(:,2))+PlotStopOffset]);



title(sprintf('%s\n%s\n%s',temp_GoodStr{i},temp_BadStr{i},temp_NovelStr{i}),'FontSize',10,'FontWeight','Bold');




%plot([0,0],[-1,max(max(Mean_PSTH_FR_EachObject+SEM_PSTH_FR_EachObject))*1.02],'-.k','LineWidth',3);
%plot([RelativeTimeMarkerForeachScene(i),RelativeTimeMarkerForeachScene(i)],[-1,max(max(SpikeFiringRateSceneOn_mean+SpikeFiringRateSceneOn_sem))*1.02],'--k','LineWidth',3)

%Add significance

%area(nanmean(MeanPSTH_Time),-0.5*min(min(SpikeFiringRateSceneOn_mean))*h_EachScene(i,:),'FaceColor',Color(i,:),'EdgeColor',Color(i,:),'FaceAlpha',0.5,'EdgeAlpha',0.5);





%drawaxis(gca, 'x', 0, 'movelabel', 1)
%h_EachScene(i,:)

%title(SceneStr{i});

%Add xlabel and ylabel
if mod(i,4)==1
    ylabel('Spikes/s');
end
if i==NumberOfObject/3-2
    
    xlabel('Time from scene onset (ms)');
end



end



%%%%%%%%%%%%%%%%%%%%%%%% Figure 5 PSTH for each object
%%%%%%%%%%%%%%%%%%%%%%%% bad%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for j=1:ceil(NumberOfObject/16)
figtitlestr{5+j-1}=sprintf('PSTH_EachFractal_%d',j);

fig{5+j-1}=PrepareFigure(108+j,'w',[160+5*(j-1),100, 1000,800],'Name',figtitlestr{5+j-1});

NumberOfObjectInEachFigure=16;
Indicator=1;
for ii=(j-1)*NumberOfObjectInEachFigure+1:min(j*NumberOfObjectInEachFigure,NumberOfObject)
    subplot(ceil(NumberOfObjectInEachFigure/4),4,Indicator);
   
shadedErrorBar(PSTH_Time_All,Mean_PSTH_FR_EachObject(ii,:),SEM_PSTH_FR_EachObject(ii,:),'lineprops',{Color(ii,:)},'transparent',1,'patchSaturation',0.3);
set(findall(gca, 'Type', 'Line'),'LineWidth',3);
%ylim([-0.5*min(min(Mean_PSTH_FR_EachObject)),max(max(Mean_PSTH_FR_EachObject+SEM_PSTH_FR_EachObject))*1.02]);
%ylim([-0.5*min(min(Mean_PSTH_FR_EachObject)),max(max(Mean_PSTH_FR_EachObject(ceil(PlotStartOffset):ceil(PlotStartOffset)+1000)+SEM_PSTH_FR_EachObject(ceil(PlotStartOffset):ceil(PlotStartOffset)+1000)))*1.02])

%xlim([min(min(PSTH_Time_All))-10,max(max(PSTH_Time_All))+10]);
%xlim([-PlotStartOffset,min(TimeMarker(:,2))+PlotStopOffset]);
%xlim([-PlotStartOffset, median(TimeMarker(:,2))+PlotStopOffset]);
%xlim([-PlotStartOffset,max(min(TimeMarker(:,2))+PlotStopOffset,400)]);
xlim([-PlotStartOffset,min(TimeMarker(:,2))+PlotStopOffset]);
%ylabel('Spikes/s');

set(gca,'LineWidth',3,'FontSize',20,'FontWeight','Bold');
%ylim([0,60]);

hold on
%Add significance level
yrange=ylim;
h_EachScene(ii,h_EachScene(ii,:)==0)=NaN;
h_EachScene(ii,PSTH_Time_All<0)=NaN;

sig=area(PSTH_Time_All,h_EachScene(ii,:)*max(yrange(1)*0.6,3),0,'FaceColor','#FFC107','EdgeColor','#FFC107','FaceAlpha',0.5,'EdgeAlpha',0.5,'ShowBaseLine','off');



set(findall(gca, 'Type', 'Line'),'LineWidth',3);
%ylabel(sprintf('Spikes/s'));
%ylabel('Firing Rate(Hz)')
%xticks([]);
%xticklabels([]);
set(gca,'LineWidth',3,'FontSize',25);
title(FractalStr{ii},'FontSize',20,'FontWeight','Bold');
  
Indicator=Indicator+1;  
    
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Figure 6 Location
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Tuning for each object

figtitlestr{7}='PSTH_EachFractal_Sep';

fig{7}=PrepareFigure(111,'w',[50,100, 1000,800],'Name',figtitlestr{7});
for m=1:length(UniqueTargetEcc)
    subplot(length(UniqueTargetEcc),1,m)
    polarwitherrorbar([UniqueTargetAngle/180*pi; UniqueTargetAngle(1)/180*pi]',[FR_Loc_mean(:,m)' FR_Loc_mean(1,m)],[FR_Loc_sem(:,m)' FR_Loc_sem(1,m)],...
    '-b',3);
hold on
%Spon
polarwitherrorbar([UniqueTargetAngle/180*pi; UniqueTargetAngle(1)/180*pi]',[FR_Spon_Loc_mean(:,m)' FR_Spon_Loc_mean(1,m)],[FR_Spon_Sem_mean(:,m)' FR_Spon_Sem_mean(1,m)],...
    '-k',1.5);

%h_p = polarplot([PrefVector(j) PrefVector(j)]/180*pi,[0 PrefVector_Amp(j)],'r');
%    set(h_p,'LineWidth',2.5);

set(gca,'ThetaDir' , 'counterclockwise');
set(gca,'ThetaZeroLocation','top')
set(gca,'FontSize',25,'FontWeight','Bold','LineWidth',3);
    
    
end









end %End of if showfigureflag




%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Output to files%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
    if isfield(OutputFileName,'FractalTuningPassive')
    ExistFlag=1;
    end
end


%%%%%%%%%%%%%%%%%%%%OutputData For current analysis
OutputData_New=[];

OutputData_New.FractalTuningPassive.TrialType='GoodTrials';

%%%%%%%%%Save overall scene response characteristics%%%%%%%%%%%%%%%%%%%%%%%
%Organize the data
NameStr={'TrialNum','AverageFR_Fractal','p_AverageFR_Fractal','AverageFR_Fractal_Good','p_AverageFR_Fractal_Good','AverageFR_Fractal_Bad','p_AverageFR_Fractal_Bad',...
'p_GoodvsBad',...
'FractalSegTimeOffset',...
'AverageFR_Spon','AverageFR_Spon_Good','AverageFR_Spon_Bad','SponSeg_Marker_OffSet',...
'RelativeTimeMarker','AbsoluteTimeMarker' ,...
'FractalNameLine','FractalIndexLine','UniqueTarget','StepSize,BinWidth',...
'AverageFR_EachFractal','AverageFR_GoodFractal','AverageFR_BadFractal',...
'FractalSelectiveIndex_All','p_FSI_All','FractalSelectiveIndex_Good','p_FSI_Good','FractalSelectiveIndex_Bad','p_FSI_Bad',...
'PSTH_Time_All_mean','PSTH_FractalAverage','p_PSTH_All','PSTH_FractalAverage_Good','PSTH_FractalAverage_Bad','p_PSTH_GoodvsBad',...
'Eye_Ecc_PSTH_All_mean','p_Eye_Ecc_PSTH_All_mean','Eye_Ecc_PSTH_Good_mean','Eye_Ecc_PSTH_Bad_mean','p_Eye_Ecc_PSTH_Bad_mean',...

};



DataLib={TrialNum,MeanFR,[h_AllvsSpon',p_AllvsSpon'],[MeanFR_Good_mean,MeanFR_Good_sem],[h_GoodvsSpon',p_GoodvsSpon'],[MeanFR_Bad_mean,MeanFR_Bad_sem],[h_BadvsSpon',p_BadvsSpon']...
[h_GoodvsBad',p_GoodvsBad'],...
[PlotStartOffset,PlotStopOffset],...
Mean_Spon,[Mean_Spon_Good_mean,Mean_Spon_Good_sem],[Mean_Spon_Bad_mean,Mean_Spon_Bad_sem],[SponEvent,SponEventEnd,PreTimeBuffer,StartOffset],...
RelativeTimeMarker,AbosoluteTimeMarker,...
TargetName_Line,TargetIndex_Line,UniqueTargetName,[StepSize,BinWidth],...
[Mean_FR_Whole_EachObject',SEM_FR_Whole_EachObject'],[Mean_FR_Whole_EachObject(GoodObject)',SEM_FR_Whole_EachObject(GoodObject)'],[Mean_FR_Whole_EachObject(BadObject)',SEM_FR_Whole_EachObject(BadObject)'],...
FractalSelectiveIndex_All,pFractal_SSI_All,FractalSelectiveIndex_Good,pFractal_SSI_Good,FractalSelectiveIndex_Bad,pFractal_SSI_Bad,...
PSTH_Time_All,[PSTH_FR_All_mean',PSTH_FR_All_sem'],[h_PSTH_All',p_PSTH_All'],[PSTH_FR_Good_mean',PSTH_FR_Good_sem'],[PSTH_FR_Bad_mean',PSTH_FR_Bad_sem'],[h_PSTH_GoodvsBad',p_PSTH_GoodvsBad'],...
[Eye_Ecc_All_mean',Eye_Ecc_All_sem'],[h_AllvsSponEye',p_AllvsSponEye'],[Eye_Ecc_Good_mean',Eye_Ecc_Good_sem'],[Eye_Ecc_Bad_mean',Eye_Ecc_Bad_sem'],[h_GoodvsBadEye',p_GoodvsBadEye'],...

};
DataStamp=containers.Map(NameStr,DataLib);


OutputData_New.FractalTuningPassive.Task=TaskType;
OutputData_New.FractalTuningPassive.TaskCode=TaskCode;
OutputData_New.FractalTuningPassive.DataStamp=DataStamp;


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
        saveas(fig{j},OutputFigureName);
        
    end

    disp('Figures have been exported to the neuron folder');

end

if ExistFlag
%Compare the old one with the new one
Task_Old=OutputData.FractalTuningPassive.TaskCode;

     if Task_Old~=TaskCode
        %Not the same task, add the new dataset into the old one
        OutputData.FractalTuningPassive(numel(OutputData.FractalTuningPassive)+1)=OutputData_New.FractalTuningPassive;
     else
         %If the same task,the same protocol,replace the old one with the new one
         
         OutputData.FractalTuningPassive=OutputData_New.FractalTuningPassive;
         
    
    
     end
else
    %Output the current file
   OutputData.FractalTuningPassive=OutputData_New.FractalTuningPassive; 
end



if OutputFlag
    cd(OutPath);
    save(OutputFileName,'OutputData');
    disp('Data Exported');
else
    disp('No data export');
end  

%{
%Variabls output to excel
                 %PreferedScene	SceneSelectiveIndex	p(SSI)	SponSelectiveIndex	p(SponSI)	ResponseLatency(Scene)	PeakFiring(Scene)	PeakTime(Scene)	DurationProportion(Scene)	FirstSigInterval(Scene)

VariablesToExcel=[PreferScene,SceneSelectiveIndex,pSceneDiff,SponSelectiveIndex,pSponDiff,ResponseLatency,PeakFiring,PeakFiringTime,PropResponseDuration,FirstIntervalDuration];
VarNum=length(VariablesToExcel);
str_header=[];
for i=1:VarNum
    str_header=[str_header,'%5.2f\t'];
end

VariablesToExcel_str=sprintf(str_header,VariablesToExcel);
disp(VariablesToExcel_str);
%}

%{
Stats1=sprintf('SSI(Total)=%3.3f;',FractalSelectiveIndex);
txt1=text(sum(BadObject)+1,max(FR_Fractal_mean),Stats1);

set(txt1,'FontSize',15)

if pFractal_Average<0.05 && pFractal_Average>=0.01
aster='(*)';
elseif pFractal_Average<0.01 && pFractal_Average>=0.001
aster='(**)';
elseif pFractal_Average<0.001
aster='(***)';
else
aster='';
end


Stats2=sprintf('p(Anova)=%3.3f%s;',pFractal_Average,aster);
txt2=text(sum(BadObject)+1,max(FR_Fractal_mean)-max(FR_Fractal_mean)/3,Stats2);

set(txt2,'FontSize',15)
if pFractal_Average<0.05
set(txt2,'Color','r')
end
%}




%%
%{
%Separate PSTH for good and bad object
TimeSeqTgtN_meanGood=TimeSeqTgtN_mean(GoodObject,:);
FRAlignTgtN_meanGood=FRAlignTgtN_mean(GoodObject,:);
FRAlignTgtN_semGood=FRAlignTgtN_sem(GoodObject,:);
FractalGood=UniqueFractal(GoodObject);
axes('position',[0.1,0.45,0.3,0.2])
for i=1:length(FractalGood)
shadedErrorBar( TimeSeqTgtN_meanGood(i,:),FRAlignTgtN_meanGood(i,:),FRAlignTgtN_semGood(i,:),'lineprops',{Color(i,:)},'transparent',1,'patchSaturation',0.3)
%plot(TimeSequence(i,:),SpikeFiringRateSceneOn_mean(i,:),'Color',Color(i,:))
hold on
end

ylim([0,max(max(FRAlignTgtN_mean+FRAlignTgtN_sem))+10])
xlim([min(min(TimeSeqTgtN_mean))-10,max(max(TimeSeqTgtN_mean))+10])
plot([0,0],[0,max(max(FRAlignTgtN_mean+FRAlignTgtN_sem))+10],'-.k','LineWidth',3);

plot([mean(RewardTime) mean(RewardTime)],[0,max(max(FRAlignTgtN_mean+FRAlignTgtN_sem))+10],'--k');
%plot([min(TgtOnTime) min(TgtOnTime)],[0,max(max(FRAlignTgtN_mean+FRAlignTgtN_sem))+10],'--k');
text((mean(RewardTime)-50),max(max(FRAlignTgtN_mean+FRAlignTgtN_sem))+10,'Rwd On','FontSize',13,'FontWeight','Bold');
plot([mean(RewardTimeOff) mean(RewardTimeOff)],[0,max(max(FRAlignTgtN_mean+FRAlignTgtN_sem))+10],'--k');
text((mean(RewardTimeOff)-50),max(max(FRAlignTgtN_mean+FRAlignTgtN_sem))+10,'Rwd Off','FontSize',13,'FontWeight','Bold');



ld=legend(string(FractalGood));
ld.Box='off';
ld.Position=[0.385   0.3    0.0854    0.0880];

%ld.String={string(UniqueScene)};


xlabel('Time from Good Object Onset')
ylabel('Spikes/s')
set(findall(gca, 'Type', 'Line'),'LineWidth',3);

set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');


%Separate PSTH for good and bad object
TimeSeqTgtN_meanBad=TimeSeqTgtN_mean(BadObject,:);
FRAlignTgtN_meanBad=FRAlignTgtN_mean(BadObject,:);
FRAlignTgtN_semBad=FRAlignTgtN_sem(BadObject,:);
FractalBad=UniqueFractal(BadObject);
axes('position',[0.5,0.45,0.3,0.2])
for i=1:length(FractalBad)
shadedErrorBar( TimeSeqTgtN_meanBad(i,:),FRAlignTgtN_meanBad(i,:),FRAlignTgtN_semBad(i,:),'lineprops',{Color(i,:)},'transparent',1,'patchSaturation',0.3)
%plot(TimeSequence(i,:),SpikeFiringRateSceneOn_mean(i,:),'Color',Color(i,:))
hold on
end

ylim([0,max(max(FRAlignTgtN_mean+FRAlignTgtN_sem))+10])
xlim([min(min(TimeSeqTgtN_mean))-10,max(max(TimeSeqTgtN_mean))+10])
plot([0,0],[0,max(max(FRAlignTgtN_mean+FRAlignTgtN_sem))+10],'-.k','LineWidth',3);

%plot([mean(RewardTime) mean(RewardTime)],[0,max(max(FRAlignTgtN_mean+FRAlignTgtN_sem))+10],'--k');
%plot([min(TgtOnTime) min(TgtOnTime)],[0,max(max(FRAlignTgtN_mean+FRAlignTgtN_sem))+10],'--k');


ld=legend(string(FractalBad));
ld.Box='off';
ld.Position=[0.8   0.3    0.0854    0.0880];

%ld.String={string(UniqueScene)};


xlabel('Time from Bad Object Onset')
ylabel('Spikes/s')
set(findall(gca, 'Type', 'Line'),'LineWidth',3);

set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');











%FRAlignTgtNOff=

%Separate PSTH for good and bad object aligned with target offset
TimeSeqTgtN_meanGoodOff=TimeSeqTgtN_meanOff(GoodObject,:);
FRAlignTgtN_meanGoodOff=FRAlignTgtN_meanOff(GoodObject,:);
FRAlignTgtN_semGoodOff=FRAlignTgtN_semOff(GoodObject,:);


axes('position',[0.1,0.15,0.3,0.2])
for i=1:length(FractalGood)
shadedErrorBar(TimeSeqTgtN_meanGoodOff(i,:),FRAlignTgtN_meanGoodOff(i,:),FRAlignTgtN_semGoodOff(i,:),'lineprops',{Color(i,:)},'transparent',1)
end
hold on
plot([0,0],[0,max(max(FRAlignTgtN_meanGoodOff+FRAlignTgtN_semGoodOff))+10],'-.k','LineWidth',3);

plot([mean(RewardTime_RefOff) mean(RewardTime_RefOff)],[0,max(max(FRAlignTgtN_meanGoodOff+FRAlignTgtN_semGoodOff))+10],'--k');

text((mean(RewardTime_RefOff)-50),max(max(FRAlignTgtN_meanGoodOff+FRAlignTgtN_semGoodOff))+10,'Rwd On','FontSize',13,'FontWeight','Bold');
plot([mean(RewardTimeOff_RefOff) mean(RewardTimeOff_RefOff)],[0,max(max(FRAlignTgtN_meanGoodOff+FRAlignTgtN_semGoodOff))+10],'--k');
text((mean(RewardTimeOff_RefOff)-50),max(max(FRAlignTgtN_meanGoodOff+FRAlignTgtN_semGoodOff))+10,'Rwd Off','FontSize',13,'FontWeight','Bold');


ylim([0,max(max(FRAlignTgtN_meanGoodOff+FRAlignTgtN_semGoodOff))+10])
xlim([min(min(TimeSeqTgtN_meanGoodOff))-10,max(max(TimeSeqTgtN_meanGoodOff))+10])
ylabel('Spikes/s')
xlabel('Time from Good Object Offset')

set(gca,'LineWidth',3)

set(findall(gca, 'Type', 'Line'),'LineWidth',3);

set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');


%Separate PSTH for good and bad object aligned with target offset
TimeSeqTgtN_meanBadOff=TimeSeqTgtN_meanOff(BadObject,:);
FRAlignTgtN_meanBadOff=FRAlignTgtN_meanOff(BadObject,:);
FRAlignTgtN_semBadOff=FRAlignTgtN_semOff(BadObject,:);


axes('position',[0.5,0.15,0.3,0.2])
for i=1:length(FractalBad)
shadedErrorBar(TimeSeqTgtN_meanBadOff(i,:),FRAlignTgtN_meanBadOff(i,:),FRAlignTgtN_semBadOff(i,:),'lineprops',{Color(i,:)},'transparent',1)
end
hold on
plot([0,0],[0,max(max(FRAlignTgtN_meanBadOff+FRAlignTgtN_semBadOff))+10],'-.k','LineWidth',3);

%area([mean(RewardTime) mean(RewardTimeOff)],[max(max(SpikeFiringRateSceneOff_mean+SpikeFiringRateSceneOff_sem))+10,max(max(SpikeFiringRateSceneOff_mean+SpikeFiringRateSceneOff_sem))+10],'FaceColor','k','FaceAlpha',.1,'EdgeAlpha',.1,'HandleVisibility','off');

ylim([0,max(max(FRAlignTgtN_meanBadOff+FRAlignTgtN_semBadOff))+10])
xlim([min(min(TimeSeqTgtN_meanBadOff))-10,max(max(TimeSeqTgtN_meanBadOff))+10])
ylabel('Spikes/s')
xlabel('Time from Bad Object Offset')

set(gca,'LineWidth',3)

set(findall(gca, 'Type', 'Line'),'LineWidth',3);

set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');





if SaveFigFlag
disp('Saving figures... ...');
%Save Figure
%FigureName=strcat('SceneTuning_',FileName(1:end-4),'_',string(SpikeChannel));
FigureName=strcat('FractalONOffPSTH_',FileName(1:end-4),'_',string(SpikeChannel));

%saveas(103, fullfile(OutPutPath, FigureName), 'jpeg');
saveas(104, fullfile(OutPutPath, FigureName), 'jpeg');

disp('Saved Successfully!');

end

%}       







end


function organized=ReorganizeEye(data);
organized=[];
for i=1:length(data)
    data_curr=data{i};
    maxnumel=max(cellfun(@numel,data_curr));
    organized{i}=cell2mat(cellfun(@(x) [x,NaN*ones(1,maxnumel-length(x))],data_curr,'uniform',0));
    
end

 maxnumel=max(cellfun(@(x) size(x,2),organized));
  organized=cell2mat(cellfun(@(x) [x,NaN*ones(size(x,1),maxnumel-size(x,2))],organized,'uniform',0)');

end
