function ObjSaccadeAnalysis_TrialHistory(Data, StartTrial, EndTrial,ShowFigureFlag,OutputFlag);
%This needs the information from all data (including the error trials)
%by Xuefei Yu 052521

ECode;
SegEvent=[TGTCD,TGTOFF];%Event Flow of interest
%SegEventEnd=[SHOWFIXCD,FIXOFF,TGTOFF,RWDOFFCD ,NaN,EXTRAFP];%Event end
SegEventEnd=[TGTOFF,NaN];%Event end
SegEventStr={'TargetOn','TargetOff'};
AmpMatchFlag=0;%1:Amplitude match on ;0:Amplitude Match off

PlotRangeFracOn=[-500,2000];
PlotRangeFracOff=[-1000,1500];


PlotStartOffset=[500,500];%in ms
PlotStopOffset=[2000,2000];%in ms

ExamplePlot=0;%PLots for presentations, turn off for regular use

ResponseLatencyThreshold=3;

%Parameter to define the sliding window to measure firing rates
StepSize=5;%in ms
BinWidth=20;%in ms


%Criterias for sexwlecting saccades
V_Threshold=40;
ContinueBin=5;
ISI_Threshold=20;

%Time bin pre saccade and time bin offset after saccade onset
PlotSacPreBuffer=500;%in ms
PlotSacPostBuffer=2000;%in ms




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
Selection_all=StartTrial:EndTrial;




%{
if strcmp(FileName,'Robin040321FRACSAC1.rst') & SpikeChannel==22
    Selection=Selection(~ismember(Selection,[14,40]));
end
%}
%FixationPointIndex=Data.FixPointIndex(Selection,:);
%FixLoc=Data.FixLoc(Selection,:,:);
%{
try
tmp=Data.ParaLib('FlipGoodBad');
FlipGoodBad=tmp(:,2);
catch
   FlipGoodBad=0*ones(sum(Selection),1);
   disp('Auto no flipping sign');
end
%}

%This is all trials, so first sort trials into different type;

TrialType=Data.TrialType;

GoodTrialIndex=TrialType==CORRECTCD;
GoodTrial=find(GoodTrialIndex==1);

ErrorTrialIndex=TrialType==BROKEFIX | TrialType==FAILED_SAC;
%ErrorTrialIndex=TrialType==BROKEFIX;
ErrorTrial=find(ErrorTrialIndex==1);
ErrorTrial_Selection=Selection_all(ismember(Selection_all,ErrorTrial))';

Selection=Selection_all(ismember(Selection_all,GoodTrial))';


%Insert null trials for those parameters
%here is only for target index and name
Data.TargetIndex=InsertNullTrials(Data.TargetIndex);

%Target Information
TargetIndex=Data.TargetIndex(Selection,2);%For this task, there is only target name, but no target index
try
TargetName=Data.TargetIndex(Selection,3);%For this task, there is only target name, but no target index
catch
    TargetName=TargetIndex;
    disp('No targetname');
end
%TargetLoc=Data.TargetLoc(Selection,:,:);

%TargetNum=Data.TargetNum(Selection,2:end);

%TargetIndex_Error=Data.TargetIndex(ErrorTrial_Selection,2);





%{
%Analysis Good versus Bad

GoodObjectIndex_Raw=TargetIndex_Line<=4;
%%
%Flip the sign of good and bad object according to flipping index

GoodObjectIndex_Raw(FlipGoodBad==1)=~GoodObjectIndex_Raw(FlipGoodBad==1);
%}





%%
%GetEyeData
EyeData=Data.EyeDataRaw;
EyeChannel_X_L=EyeData(Selection,1);
EyeChannel_Y_L=EyeData(Selection,2);
%EyeChannel_X_R=EyeData(Selection,3);
%EyeChannel_Y_R=EyeData(Selection,4);

Eye_Eccentricity=cellfun(@(x,y) sqrt(x.^2+y.^2),EyeChannel_X_L,EyeChannel_Y_L,'uniform',0);

%EyeBinWidth=Data.EyeBinWidth;


EventChannel=Data.EventChannel(Selection,:);
EventTimeChannel=Data.EventTimeChannel(Selection,:);
EventBin=Data.EventBin(Selection,:);


EventChannelError=Data.EventChannel(ErrorTrial_Selection,:);
EventTimeChannelError=Data.EventTimeChannel(ErrorTrial_Selection,:);
EventBinError=Data.EventBin(ErrorTrial_Selection,:);

%Algin Error trials on errorbeep
%AlignMarker= ERR_BEEP  ;

SegEventError=[ERR_BEEP ];%Event Flow of interest
%SegEventEnd=[SHOWFIXCD,FIXOFF,TGTOFF,RWDOFFCD ,NaN,EXTRAFP];%Event end
SegEventEndError=[NaN];%Event end
SegEventStrError={'ErrorEnd'};



PlotStartOffsetError=[500];%in ms
PlotStopOffsetError=[2000];%in ms


%
%%Reproduce goodbad index from event channel
%ObjectIndex=
ObjectIndex=ReproduceFromEvent(EventChannel,[GOODOBJ  BADOBJ]);

GoodObjectIndex=(ObjectIndex==GOODOBJ)';
%}

%Load the spike channel
SpikeChannelWhole=Data.SpikeChannel;
SpikeTimeWhole=Data.SpikeTimeData(:,:,Selection);
ChannelIDWhole=Data.ChannelID;
ClusterIDWhole=Data.ClusterID;

SpikeTimeError_Whole=Data.SpikeTimeData(:,:,ErrorTrial_Selection);

MultiChannel=0;
if length(SpikeChannelWhole)>1
    MultiChannel=1;
    
end
for spk=1:length(SpikeChannelWhole)
    SpikeChannel=SpikeChannelWhole(spk);
    SpikeTime=squeeze(SpikeTimeWhole(spk,:,:))';
    ChannelID=ChannelIDWhole(spk);
    ClusterID=ClusterIDWhole(spk);
    
    SpikeTimeError=squeeze(SpikeTimeError_Whole(spk,:,:))';
    if size(SpikeTimeError,2)==1
        SpikeTimeError = SpikeTimeError';
    end


%For good trials
SegmentData=EventSegTool(SegEvent,SegEventEnd,PlotStartOffset,PlotStopOffset,EventChannel,EventTimeChannel,SpikeTime,StepSize,BinWidth);
%Relative time bin for each event
RelativeTimeMarker=SegmentData.RelativeTimeMarker;%=EventTimeMarker;
%Absolute time marker for each event
AbosoluteTimeMarker=SegmentData.AbsoluteTimeMarker;

ObjOnTimeMarker=nanmean(RelativeTimeMarker{1});
ObjOffTimeMarker=nanmean(RelativeTimeMarker{2});

%PSTH for each event
MeanPSTH=SegmentData.MeanPSTH;%=[PSTH_Event_Time',PSTH_Event'];
MeanPSTH_Time=MeanPSTH(:,1);
MeanPSTH_FR=MeanPSTH(:,2);

PSTH_Time_FraOn=nanmean(MeanPSTH_Time{1},1);
PSTH_FraOn_GoodTrials=MeanPSTH_FR{1};


PSTH_GoodOn_GoodTrials=PSTH_FraOn_GoodTrials(GoodObjectIndex,:);
PSTH_BadOn_GoodTrials=PSTH_FraOn_GoodTrials(~GoodObjectIndex,:);

PSTH_GoodOn_GoodTrials_Mean=nanmean(PSTH_GoodOn_GoodTrials,1);
PSTH_BadOn_GoodTrials_Mean=nanmean(PSTH_BadOn_GoodTrials,1);

PSTH_GoodOn_GoodTrials_Sem=nanstd(PSTH_GoodOn_GoodTrials,[],1)./sqrt(sum(~isnan(PSTH_GoodOn_GoodTrials),1));
PSTH_BadOn_GoodTrials_Sem=nanstd(PSTH_BadOn_GoodTrials,[],1)./sqrt(sum(~isnan(PSTH_BadOn_GoodTrials),1));


PSTH_Time_FraOff=nanmean(MeanPSTH_Time{2},1);
PSTH_FraOff_GoodTrials=MeanPSTH_FR{2};


PSTH_GoodOff_GoodTrials=PSTH_FraOff_GoodTrials(GoodObjectIndex,:);
PSTH_BadOff_GoodTrials=PSTH_FraOff_GoodTrials(~GoodObjectIndex,:);

PSTH_GoodOff_GoodTrials_Mean=nanmean(PSTH_GoodOff_GoodTrials,1);
PSTH_BadOff_GoodTrials_Mean=nanmean(PSTH_BadOff_GoodTrials,1);

PSTH_GoodOff_GoodTrials_Sem=nanstd(PSTH_GoodOff_GoodTrials,[],1)./sqrt(sum(~isnan(PSTH_GoodOff_GoodTrials),1));
PSTH_BadOff_GoodTrials_Sem=nanstd(PSTH_BadOff_GoodTrials,[],1)./sqrt(sum(~isnan(PSTH_BadOff_GoodTrials),1));

PSTH_GoodOff_GoodTrials_Mean=FewNumCorrection(PSTH_GoodOff_GoodTrials,PSTH_GoodOff_GoodTrials_Mean,3);
PSTH_GoodOff_GoodTrials_Sem=FewNumCorrection(PSTH_GoodOff_GoodTrials,PSTH_GoodOff_GoodTrials_Sem,3);

PSTH_BadOff_GoodTrials_Mean=FewNumCorrection(PSTH_BadOff_GoodTrials,PSTH_BadOff_GoodTrials_Mean,3);
PSTH_BadOff_GoodTrials_Sem=FewNumCorrection(PSTH_BadOff_GoodTrials,PSTH_BadOff_GoodTrials_Sem,3);


%For error trials

SegmentData=EventSegTool(SegEventError,SegEventEndError,PlotStartOffsetError,PlotStopOffsetError,EventChannelError,EventTimeChannelError,SpikeTimeError,StepSize,BinWidth);

%Relative time bin for each event
RelativeTimeMarker=SegmentData.RelativeTimeMarker;%=EventTimeMarker;
%Absolute time marker for each event
AbosoluteTimeMarker=SegmentData.AbsoluteTimeMarker;

ErrorTimeMarker=nanmean(RelativeTimeMarker{1});

%PSTH for each event
MeanPSTH=SegmentData.MeanPSTH;%=[PSTH_Event_Time',PSTH_Event'];
MeanPSTH_Time=MeanPSTH{1};
MeanPSTH_FR=MeanPSTH{2};


PSTH_TrialEndError_Time=nanmean(MeanPSTH_Time,1);

PSTH_TrialEndError_Mean=nanmean(MeanPSTH_FR,1);
PSTH_TrialEndError_Sem=nanstd(MeanPSTH_FR,[],1)./sqrt(sum(~isnan(MeanPSTH_FR),1));


PSTH_TrialEndError_Mean=FewNumCorrection(MeanPSTH_FR,PSTH_TrialEndError_Mean,3);
PSTH_TrialEndError_Sem=FewNumCorrection(MeanPSTH_FR,PSTH_TrialEndError_Sem,3);


%Trial History
TrialTypeLast=[NaN;TrialType(1:end-1)];
TrialTypeLast_GoodTrials=TrialTypeLast(GoodTrialIndex);

TrialTypeCurr_GoodTrials=TrialType(GoodTrialIndex);

%Separate into last trial error, last trial good, last trial bad 

LastError=~(TrialTypeLast_GoodTrials==CORRECTCD) & ~isnan(TrialTypeLast_GoodTrials);



GoodObjectIndexLast=[NaN; GoodObjectIndex(1:end-1)];

LastGood=~LastError & GoodObjectIndexLast==1;
LastBad=~LastError & GoodObjectIndexLast==0;

CurrGood=GoodObjectIndex==1;
CurrBad=GoodObjectIndex==0;

LastGoodCurrGood=LastGood & CurrGood;
LastGoodCurrBad=LastGood & CurrBad;

LastBadCurrGood=LastBad & CurrGood;
LastBadCurrBad=LastBad & CurrBad;

PSTH_FraOnHisotry_GoodTrials=PSTH_FraOn_GoodTrials(1:end,:);
PSTH_FracOn_LastError=PSTH_FraOnHisotry_GoodTrials(LastError,:);

PSTH_FracOn_LastGood=PSTH_FraOnHisotry_GoodTrials(LastGood,:);
PSTH_FracOn_LastBad=PSTH_FraOnHisotry_GoodTrials(LastBad,:);


PSTH_FracOn_LastGoodCurrGood=PSTH_FraOnHisotry_GoodTrials(LastGoodCurrGood,:);
PSTH_FracOn_LastGoodCurrBad=PSTH_FraOnHisotry_GoodTrials(LastGoodCurrBad,:);


PSTH_FracOn_LastBadCurrGood=PSTH_FraOnHisotry_GoodTrials(LastBadCurrGood,:);
PSTH_FracOn_LastBadCurrBad=PSTH_FraOnHisotry_GoodTrials(LastBadCurrBad,:);



[PSTH_FracOn_LastError_Mean,PSTH_FracOn_LastError_Sem]=StatisticsPSTH(PSTH_FracOn_LastError,3);
[PSTH_FracOn_LastGood_Mean,PSTH_FracOn_LastGood_Sem]=StatisticsPSTH(PSTH_FracOn_LastGood,3);
[PSTH_FracOn_LastBad_Mean,PSTH_FracOn_LastBad_Sem]=StatisticsPSTH(PSTH_FracOn_LastBad,3);

[PSTH_FracOn_LastGoodCurrGood_Mean,PSTH_FracOn_LastGoodCurrGood_Sem]=StatisticsPSTH(PSTH_FracOn_LastGoodCurrGood,3);
[PSTH_FracOn_LastGoodCurrBad_Mean,PSTH_FracOn_LastGoodCurrBad_Sem]=StatisticsPSTH(PSTH_FracOn_LastGoodCurrBad,3);

[PSTH_FracOn_LastBadCurrGood_Mean,PSTH_FracOn_LastBadCurrGood_Sem]=StatisticsPSTH(PSTH_FracOn_LastBadCurrGood,3);
[PSTH_FracOn_LastBadCurrBad_Mean,PSTH_FracOn_LastBadCurrBad_Sem]=StatisticsPSTH(PSTH_FracOn_LastBadCurrBad,3);



PSTH_FraOffHisotry_GoodTrials=PSTH_FraOff_GoodTrials(1:end,:);
PSTH_FracOff_LastError=PSTH_FraOffHisotry_GoodTrials(LastError,:);

PSTH_FracOff_LastGood=PSTH_FraOffHisotry_GoodTrials(LastGood,:);
PSTH_FracOff_LastBad=PSTH_FraOffHisotry_GoodTrials(LastBad,:);


PSTH_FracOff_LastGoodCurrGood=PSTH_FraOffHisotry_GoodTrials(LastGoodCurrGood,:);
PSTH_FracOff_LastGoodCurrBad=PSTH_FraOffHisotry_GoodTrials(LastGoodCurrBad,:);


PSTH_FracOff_LastBadCurrGood=PSTH_FraOffHisotry_GoodTrials(LastBadCurrGood,:);
PSTH_FracOff_LastBadCurrBad=PSTH_FraOffHisotry_GoodTrials(LastBadCurrBad,:);



[PSTH_FracOff_LastError_Mean,PSTH_FracOff_LastError_Sem]=StatisticsPSTH(PSTH_FracOff_LastError,3);
[PSTH_FracOff_LastGood_Mean,PSTH_FracOff_LastGood_Sem]=StatisticsPSTH(PSTH_FracOff_LastGood,3);
[PSTH_FracOff_LastBad_Mean,PSTH_FracOff_LastBad_Sem]=StatisticsPSTH(PSTH_FracOff_LastBad,3);

[PSTH_FracOff_LastGoodCurrGood_Mean,PSTH_FracOff_LastGoodCurrGood_Sem]=StatisticsPSTH(PSTH_FracOff_LastGoodCurrGood,3);
[PSTH_FracOff_LastGoodCurrBad_Mean,PSTH_FracOff_LastGoodCurrBad_Sem]=StatisticsPSTH(PSTH_FracOff_LastGoodCurrBad,3);

[PSTH_FracOff_LastBadCurrGood_Mean,PSTH_FracOff_LastBadCurrGood_Sem]=StatisticsPSTH(PSTH_FracOff_LastBadCurrGood,3);
[PSTH_FracOff_LastBadCurrBad_Mean,PSTH_FracOff_LastBadCurrBad_Sem]=StatisticsPSTH(PSTH_FracOff_LastBadCurrBad,3);







%%
if ShowFigureFlag
    %First figure, align on fractal on and fractal off
    figtitlestr{1}='PSTHFracOnAndFracOff';

     fig{1}=PrepareFigure(100,'w',[50,100, 1200,800],'Name',figtitlestr{1});
    
     sel=PSTH_Time_FraOn>PlotRangeFracOn(1) & PSTH_Time_FraOn<PlotRangeFracOn(2);
    subplot(2,1,1);%FracOn
    shadedErrorBar(PSTH_Time_FraOn( sel),PSTH_GoodOn_GoodTrials_Mean( sel),PSTH_GoodOn_GoodTrials_Sem( sel),'lineprops',{[1    0    0]},'transparent',1,'patchSaturation',0.3);

    hold on
     shadedErrorBar(PSTH_Time_FraOn( sel),PSTH_BadOn_GoodTrials_Mean( sel),PSTH_BadOn_GoodTrials_Sem( sel),'lineprops',{[0    0    1]},'transparent',1,'patchSaturation',0.3);
     
     set(findall(gca, 'Type', 'Line'),'LineWidth',3);
set(gca,'FontSize',15,'FontWeight','Bold');
set(gca,'LineWidth',3);

xlabel('Time from Fractal Onset(ms)','FontSize',25,'FontWeight','Bold');
ylabel('Spikes/s','FontSize',25,'FontWeight','Bold');


    subplot(2,1,2);%FracOff
    sel=PSTH_Time_FraOff>PlotRangeFracOff(1) & PSTH_Time_FraOff<PlotRangeFracOff(2);
    shadedErrorBar(PSTH_Time_FraOff(sel),PSTH_GoodOff_GoodTrials_Mean(sel),PSTH_GoodOff_GoodTrials_Sem(sel),'lineprops',{[1    0    0]},'transparent',1,'patchSaturation',0.3);

    hold on
     shadedErrorBar(PSTH_Time_FraOff(sel),PSTH_BadOff_GoodTrials_Mean(sel),PSTH_BadOff_GoodTrials_Sem(sel),'lineprops',{[0    0    1]},'transparent',1,'patchSaturation',0.3);
     
     sel=PSTH_TrialEndError_Time>PlotRangeFracOff(1) & PSTH_TrialEndError_Time<PlotRangeFracOff(2);
      shadedErrorBar(PSTH_TrialEndError_Time(sel),PSTH_TrialEndError_Mean(sel),PSTH_TrialEndError_Sem(sel),'lineprops',{[0.5    0.5    0.5]},'transparent',1,'patchSaturation',0.3);

    
    set(findall(gca, 'Type', 'Line'),'LineWidth',3);
set(gca,'FontSize',15,'FontWeight','Bold');
set(gca,'LineWidth',3);

xlabel('Time from Fractal Offset(ms)','FontSize',25,'FontWeight','Bold');
ylabel('Spikes/s','FontSize',25,'FontWeight','Bold');
    
   %Add title

title_msg=sprintf(' %s, SpikeChannel: %d, TrialNum: %d ',FileName,SpikeChannel,TrialNum);
a=annotation('textbox',[0.25,0.68,0.3,0.3],'String',title_msg,'FitBoxToText','on');
a.FontSize=25;
a.FontWeight='Bold';
a.LineStyle='none';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Figure 2, trial history

figtitlestr{2}='PSTHHistoryEffect_FracOn';

fig{2}=PrepareFigure(101,'w',[70,100, 1200,800],'Name',figtitlestr{2});
     
    subplot(3,2,1)
     sel=PSTH_Time_FraOn>PlotRangeFracOn(1) & PSTH_Time_FraOn<PlotRangeFracOn(2);
    
    shadedErrorBar(PSTH_Time_FraOn( sel),PSTH_FracOn_LastError_Mean( sel),PSTH_FracOn_LastError_Sem( sel),'lineprops',{[0.5    0.5    0.5]},'transparent',1,'patchSaturation',0.3);

    hold on
     shadedErrorBar(PSTH_Time_FraOn( sel),PSTH_FracOn_LastGood_Mean( sel),PSTH_FracOn_LastGood_Sem( sel),'lineprops',{[1    0    0]},'transparent',1,'patchSaturation',0.3);
      shadedErrorBar(PSTH_Time_FraOn( sel),PSTH_FracOn_LastBad_Mean( sel),PSTH_FracOn_LastBad_Sem( sel),'lineprops',{[0    0    1]},'transparent',1,'patchSaturation',0.3);
 
     
     set(findall(gca, 'Type', 'Line'),'LineWidth',3);
     legend({'LastError','LastGood','LastBad'});
set(gca,'FontSize',15,'FontWeight','Bold');
set(gca,'LineWidth',3);


ylabel('Spikes/s','FontSize',25,'FontWeight','Bold');


   subplot(3,2,3)
    
     shadedErrorBar(PSTH_Time_FraOn( sel),PSTH_FracOn_LastGoodCurrGood_Mean( sel),PSTH_FracOn_LastGoodCurrGood_Sem( sel),'lineprops',{[1    0    0]},'transparent',1,'patchSaturation',0.3);
      hold on
     shadedErrorBar(PSTH_Time_FraOn( sel),PSTH_FracOn_LastBadCurrGood_Mean( sel),PSTH_FracOn_LastBadCurrGood_Sem( sel),'lineprops',{[0    0    1]},'transparent',1,'patchSaturation',0.3);
 
  %
     set(findall(gca, 'Type', 'Line'),'LineWidth',3);
     legend({'LastGoodCurrGood','LastBadCurrGood'});
set(gca,'FontSize',15,'FontWeight','Bold');
set(gca,'LineWidth',3);


ylabel('Spikes/s','FontSize',25,'FontWeight','Bold');
%}
     ylabel('Spikes/s','FontSize',25,'FontWeight','Bold');
     
 subplot(3,2,5)
    
     shadedErrorBar(PSTH_Time_FraOn( sel),PSTH_FracOn_LastGoodCurrBad_Mean( sel),PSTH_FracOn_LastGoodCurrBad_Sem( sel),'lineprops',{[0.5    0    0]},'transparent',1,'patchSaturation',0.3);
      hold on
     shadedErrorBar(PSTH_Time_FraOn( sel),PSTH_FracOn_LastBadCurrBad_Mean( sel),PSTH_FracOn_LastBadCurrBad_Sem( sel),'lineprops',{[0    0    0.5]},'transparent',1,'patchSaturation',0.3);
 
    
     set(findall(gca, 'Type', 'Line'),'LineWidth',3);
     legend({'LastGoodCurrBad','LastBadCurrBad'});
set(gca,'FontSize',15,'FontWeight','Bold');
set(gca,'LineWidth',3);

xlabel('Time from Fractal Onset(ms)','FontSize',25,'FontWeight','Bold');
ylabel('Spikes/s','FontSize',25,'FontWeight','Bold');


       
    subplot(3,2,2)
     sel=PSTH_Time_FraOff>PlotRangeFracOff(1) & PSTH_Time_FraOff<PlotRangeFracOff(2);
    
    shadedErrorBar(PSTH_Time_FraOff( sel),PSTH_FracOff_LastError_Mean( sel),PSTH_FracOff_LastError_Sem( sel),'lineprops',{[0.5    0.5    0.5]},'transparent',1,'patchSaturation',0.3);

    hold on
     shadedErrorBar(PSTH_Time_FraOff( sel),PSTH_FracOff_LastGood_Mean( sel),PSTH_FracOff_LastGood_Sem( sel),'lineprops',{[1    0    0]},'transparent',1,'patchSaturation',0.3);
      shadedErrorBar(PSTH_Time_FraOff( sel),PSTH_FracOff_LastBad_Mean( sel),PSTH_FracOff_LastBad_Sem( sel),'lineprops',{[0    0    1]},'transparent',1,'patchSaturation',0.3);
 
     
     set(findall(gca, 'Type', 'Line'),'LineWidth',3);
     legend({'LastError','LastGood','LastBad'});
set(gca,'FontSize',15,'FontWeight','Bold');
set(gca,'LineWidth',3);


ylabel('Spikes/s','FontSize',25,'FontWeight','Bold');


   subplot(3,2,4)
    
     shadedErrorBar(PSTH_Time_FraOff( sel),PSTH_FracOff_LastGoodCurrGood_Mean( sel),PSTH_FracOff_LastGoodCurrGood_Sem( sel),'lineprops',{[1    0    0]},'transparent',1,'patchSaturation',0.3);
      hold on
     shadedErrorBar(PSTH_Time_FraOff( sel),PSTH_FracOff_LastBadCurrGood_Mean( sel),PSTH_FracOff_LastBadCurrGood_Sem( sel),'lineprops',{[0    0    1]},'transparent',1,'patchSaturation',0.3);
 
  %
     set(findall(gca, 'Type', 'Line'),'LineWidth',3);
     legend({'LastGoodCurrGood','LastBadCurrGood'});
set(gca,'FontSize',15,'FontWeight','Bold');
set(gca,'LineWidth',3);


ylabel('Spikes/s','FontSize',25,'FontWeight','Bold');
%}
     ylabel('Spikes/s','FontSize',25,'FontWeight','Bold');
     
 subplot(3,2,6)
    
     shadedErrorBar(PSTH_Time_FraOff( sel),PSTH_FracOff_LastGoodCurrBad_Mean( sel),PSTH_FracOff_LastGoodCurrBad_Sem( sel),'lineprops',{[0.5    0    0]},'transparent',1,'patchSaturation',0.3);
      hold on
     shadedErrorBar(PSTH_Time_FraOff( sel),PSTH_FracOff_LastBadCurrBad_Mean( sel),PSTH_FracOff_LastBadCurrBad_Sem( sel),'lineprops',{[0    0    0.5]},'transparent',1,'patchSaturation',0.3);
 
    
     set(findall(gca, 'Type', 'Line'),'LineWidth',3);
     legend({'LastGoodCurrBad','LastBadCurrBad'});
set(gca,'FontSize',15,'FontWeight','Bold');
set(gca,'LineWidth',3);

xlabel('Time from Fractal Offset(ms)','FontSize',25,'FontWeight','Bold');
ylabel('Spikes/s','FontSize',25,'FontWeight','Bold');
  
    
end


end %End of the spikechannel

end

function ObjectIndex=ReproduceFromEvent(event,code)
for i=1:size(event,1)
    
    ObjectIndex(i)=event(i,ismember(event(i,:),code));
    
  
    
end



end
function data=FewNumCorrection(mat,data_ori,the);
data=data_ori;
colnum=sum(~isnan(mat),1);
inavidcol=colnum<the;
data(inavidcol)=NaN;


end
function [mean,sem]=StatisticsPSTH(psth,the);

mean=nanmean(psth,1);
sem=nanstd(psth,[],1)./sqrt(sum(~isnan(psth),1));

mean=FewNumCorrection(psth,mean,the);
sem=FewNumCorrection(psth,sem,the);


end





