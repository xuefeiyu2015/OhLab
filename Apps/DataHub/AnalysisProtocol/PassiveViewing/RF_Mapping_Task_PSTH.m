function RF_Mapping_Task_PSTH(Data, StartTrial, EndTrial,ShowFigureFlag,OutputFlag);
ECode;


SegEvent=[TGTCD];%Event Flow of interest
%SegEventEnd=[SHOWFIXCD,FIXOFF,TGTOFF,RWDOFFCD ,NaN,EXTRAFP];%Event end
SegEventEnd=[TGTOFF];%Event end
SegEventStr={'TargetOn'};



PlotStartOffset=[200];%in ms
PlotStopOffset=[100];%in ms

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


TargetIndex=Data.TargetIndex(Selection,2:end,1);
TargetName=Data.TargetIndex(Selection,2:end,2);
%TargetLoc=Data.TargetLoc(Selection,:,:);

TargetNum=Data.TargetNum(Selection,2:end);


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

%Get the average firing rate 200ms before scene onset as the baseline
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

%Select eye trace
Eye_Eccentricity_Tgt=Eye_Eccentricity_Tgt(NoNNaNSelector,:);

%Analysis Good versus Bad
GoodObjectIndex=mod(TargetIndex_Line,2)==1;

%Eye eccentricity Good v.s. Bad
Eye_Eccentricity_Tgt_Good=Eye_Eccentricity_Tgt(GoodObjectIndex,:);
Eye_Eccentricity_Tgt_Bad=Eye_Eccentricity_Tgt(~GoodObjectIndex,:);


Eye_Ecc_Good_mean=nanmean(Eye_Eccentricity_Tgt_Good);
Eye_Ecc_Bad_mean=nanmean(Eye_Eccentricity_Tgt_Bad);

Eye_Ecc_Good_sem=nanstd(Eye_Eccentricity_Tgt_Good)./sqrt(sum(~isnan(Eye_Eccentricity_Tgt_Good),1));
Eye_Ecc_Bad_sem=nanstd(Eye_Eccentricity_Tgt_Bad)./sqrt(sum(~isnan(Eye_Eccentricity_Tgt_Bad),1));

Eye_Ecc_All_mean=nanmean(Eye_Eccentricity_Tgt);
Eye_Ecc_All_sem=nanstd(Eye_Eccentricity_Tgt)./sqrt(sum(~isnan(Eye_Eccentricity_Tgt),1));



%Overall average
MeanFR_Good=MeanFRRaw(GoodObjectIndex);
MeanFR_Bad=MeanFRRaw(~GoodObjectIndex);

MeanFR_Good_mean=nanmean(MeanFR_Good);
MeanFR_Bad_mean=nanmean(MeanFR_Bad);

MeanFR_Good_sem=nanstd(MeanFR_Good)/sqrt(length(MeanFR_Good));
MeanFR_Bad_sem=nanstd(MeanFR_Bad)/sqrt(length(MeanFR_Bad));


[h_GoodvsBad,p_GoodvsBad]=ttest2(MeanFR_Good,MeanFR_Bad);

%Spon for good and bad
Mean_Spon_Good=Mean_Spon_ForEachTrial(GoodObjectIndex);
Mean_Spon_Bad=Mean_Spon_ForEachTrial(~GoodObjectIndex);

Mean_Spon_Good_mean=nanmean(Mean_Spon_ForEachTrial(GoodObjectIndex));
Mean_Spon_Bad_mean=nanmean(Mean_Spon_ForEachTrial(~GoodObjectIndex));

Mean_Spon_Good_sem=nanstd(Mean_Spon_ForEachTrial(GoodObjectIndex))/sqrt(sum(GoodObjectIndex));
Mean_Spon_Bad_sem=nanstd(Mean_Spon_ForEachTrial(~GoodObjectIndex))/sqrt(sum(~GoodObjectIndex));



%Comparison between Good/Bad and Spon firing rate
[h_GoodvsSpon,p_GoodvsSpon]=ttest(MeanFR_Good,Mean_Spon_Good);
[h_BadvsSpon,p_BadvsSpon]=ttest(MeanFR_Bad,Mean_Spon_Bad);
[h_AllvsSpon,p_AllvsSpon]=ttest(MeanFRRaw,Mean_Spon_ForEachTrial);

%Comparison between Good/Bad and Spon Eye movement
Spon_Eyes=Eye_Eccentricity_Tgt(:,1:PlotStartOffset);
Spon_Eyes_reshape=reshape(Spon_Eyes,1,numel(Spon_Eyes));

[p_AllvsSponEye,h_AllvsSponEye]=ttest2_group(Eye_Eccentricity_Tgt,Spon_Eyes_reshape','both');
[p_GoodvsBadEye,h_GoodvsBadEye]=ttest2_group(Eye_Eccentricity_Tgt_Good,Eye_Eccentricity_Tgt_Bad,'both');


%PSTH
PSTH_Time_Good=MeanPSTH_Time(GoodObjectIndex,:);
PSTH_FR_Good=MeanPSTH_FR(GoodObjectIndex,:);

PSTH_Time_Bad=MeanPSTH_Time(~GoodObjectIndex,:);
PSTH_FR_Bad=MeanPSTH_FR(~GoodObjectIndex,:);

PSTH_Time_All=nanmean(MeanPSTH_Time,1);
PSTH_FR_All_mean=nanmean(MeanPSTH_FR,1);
PSTH_FR_All_sem=nanstd(MeanPSTH_FR,1)./sqrt(sum(~isnan(MeanPSTH_Time),1));

PSTH_FR_Good_mean=nanmean(PSTH_FR_Good,1);
PSTH_FR_Bad_mean=nanmean(PSTH_FR_Bad,1);

PSTH_FR_Good_sem=nanstd(PSTH_FR_Good,[],1)./sqrt(sum(~isnan(PSTH_FR_Good),1));
PSTH_FR_Bad_sem=nanstd(PSTH_FR_Bad,[],1)./sqrt(sum(~isnan(PSTH_FR_Bad),1));

%PSTH p value from the spon
%[p_PSTH_Good,h_PSTH_Good]=ttest2_group(PSTH_FR_Good,Mean_Spon_ForEachTrial,'both');
[p_PSTH_All,h_PSTH_All]=ttest2_group(MeanPSTH_FR,Mean_Spon_ForEachTrial,'both');
[p_PSTH_GoodvsBad,h_PSTH_GoodvsBad]=ttest2_group(PSTH_FR_Good,PSTH_FR_Bad,'both');

 

UniqueTargetName=unique(TargetName_Line);
NumberOfObject=length(UniqueTargetName);







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
