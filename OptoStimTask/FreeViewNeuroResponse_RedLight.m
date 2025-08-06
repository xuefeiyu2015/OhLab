function FreeViewNeuroResponse_RedLight(Data, StartTrial, EndTrial,ShowFigureFlag,OutputFlag);

ECode;

TaskType=Data.TaskType;
DataPath= Data.Path;
FileName=Data.FileName;
TaskCode=Data.TaskCode;

Selection=StartTrial:EndTrial;%Select the trials according to the pannel

useReISI = 0;
PlotExample = 1;
      
SegEvent = [ON_STIM];
SegEventEnd = [OFF_STIM];
SegEventStr={'StimOn'};

PlotStartOffset=[100];%in ms
PlotStopOffset=[100];%in ms

%Sequence information from the library
ParaLib=Data.ParaLib;
if ismember('stimOn',keys(ParaLib))
   stimOnFlag = ParaLib('stimOn');
   TrialStimOn = unique(stimOnFlag(stimOnFlag(:,2) == 1,1));
   stimOn = [Selection',0*ones(length(Selection),1)];
   stimOn(ismember(stimOn(:,1),TrialStimOn),2) = 1;
  
   stimOn= stimOn(:,2);
%   stimOn=stimOn(Selection,:);
   
end
if ismember('stimPower',keys(ParaLib))
   stimPower_ori = ParaLib('stimPower');
   stimPower_ori=stimPower_ori(:,2);
   
  
   
end

if ismember('preStimDur,StimDur,PostStimDur',keys(ParaLib))
   stimdur_info = ParaLib('preStimDur,StimDur,PostStimDur');
   stimDur=stimdur_info(:,3);
   
   stimDur=stimDur(Selection,:);
   
end

%Manuel Modification of the power by mistake
if ismember(FileName,{'Robin061423VIDEOVIEW1.rst','Robin062223VIDEOVIEW1.rst','Robin071323VIDEOVIEW1.rst','Robin072023VIDEOVIEW1.rst','Robin072723VIDEOVIEW1.rst'})
    stimPower = RescuePower(FileName,stimPower_ori);
    
else
    stimPower = stimPower_ori;
end

 stimPower=stimPower(Selection,:);



%Parameter to define the sliding window to measure firing rates
StepSize=5;%in ms
BinWidth=5;%in ms

%{
%Event marker for segmentation
SegEvent=[PRE_STIM];%Event Flow of interest
SegEventEnd=NaN*ones(1,length(SegEvent));
%SegEventEnd(1) =SCENE_OFF;

SegEventStr={'FixOn','1TGT','2TGT','3TGT','4TGT','5TGT','LastFix'};
MarkerSizeEvent = 20*ones(1,length(SegEventStr));
MarkerSizeEvent([1,end]) = 10;

MarkerColorEvent = {'#44AA99','#FE6100','#FE6100','#FE6100','#FE6100','#FE6100','#44AA99'};
DistratorColorEvent = {'#44AA99','#648FFF','#648FFF','#648FFF','#648FFF','#648FFF','#44AA99'};

PlotStartOffset=500*ones(1,length(SegEvent));%in ms
PlotStopOffset=NaN*ones(1,length(SegEvent));%in ms
%}

%{
%GetEyeData
EyeData=Data.EyeDataRaw;
EyeChannel_X_L=EyeData(Selection,1);
EyeChannel_Y_L=EyeData(Selection,2);
%EyeChannel_X_R=EyeData(Selection,3);
%EyeChannel_Y_R=EyeData(Selection,4);
Eye_Eccentricity=cellfun(@(x,y) sqrt(x.^2+y.^2),EyeChannel_X_L,EyeChannel_Y_L,'uniform',0);


EyeBinWidth=Data.EyeBinWidth;
EyeCalInfo=Data.EyeCalInfo;%Each row: X_Gain,X_Bias; Y_Gain; Y_bias; X=(raw+X_bias)*X_Gain
EyeCalInfo_L=EyeCalInfo(1,:);
%}
%Load event data
EventChannel=Data.EventChannel(Selection,:);
EventTimeChannel=Data.EventTimeChannel(Selection,:);
%EventBin=Data.EventBin(Selection,:);



%% Load the spike channel
SpikeChannelWhole=Data.SpikeChannel;
SpikeTimeWhole=Data.SpikeTimeData(:,:,Selection);
ChannelIDWhole=Data.ChannelID;
ClusterIDWhole=Data.ClusterID;

MultiChannel=0;
if length(SpikeChannelWhole)>1
    MultiChannel=1;
    
end

for spk=1:length(SpikeChannelWhole)
    SpikeChannel=SpikeChannelWhole(spk);
    
        SpikeTime=squeeze(SpikeTimeWhole(spk,:,:))';
        %{
    if SpikeChannel<4
        SpikeTime = SpikeTime;%Correct for the delay of online channel
    end
        %}
    ChannelID=ChannelIDWhole(spk);
    ClusterID=ClusterIDWhole(spk);

    



    
if ChannelID<4
    SpikeTime = SpikeTime-3;%Compensate for the online delay;
end

FR_All = (sum(~isnan(SpikeTime),2)-2)./SpikeTime(end)*1000;

FR_Mean = mean(FR_All);


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

Interval = nanmean(TimeMarker);

%For Raster
Raster=SegmentData.RasterAfterAlign{1};

SpikeCount = SegmentData.PSTH_Event_SpikeCount{1}/BinWidth*1000;



%PSTH for each event
MeanPSTH=SegmentData.MeanPSTH;%=[PSTH_Event_Time',PSTH_Event'];
MeanPSTH_Time=MeanPSTH{1};
MeanPSTH_FR=MeanPSTH{2};

MeanPSTH_Time_Mean = nanmean(MeanPSTH_Time);

%Seperate into stim trial and sham trial
StimType=ReproduceFromEvent(EventChannel,[STIM1,SHAM0])';

StimTrial = StimType == STIM1;
ControlTrial = ~StimTrial;

%Mean and std for z score
%Mean_FR_z = mean(MeanFRRaw(ControlTrial));
%Std_FR_z = std(MeanFRRaw(ControlTrial));

Mean_FR_z = mean(FR_All(ControlTrial));
Std_FR_z= std(FR_All(ControlTrial));






%{
if Std_FR_z <0.0000001
    Mean_FR_z = mean(MeanFRRaw);
    Std_FR_z = std(MeanFRRaw);

end
%}
%PSTH_Stim = MeanPSTH_FR(StimTrial(stimOn==1),:);
%PSTH_Control = MeanPSTH_FR(ControlTrial(stimOn==1),:);

PSTH_Stim = MeanPSTH_FR(StimTrial,:);
PSTH_Control = MeanPSTH_FR(ControlTrial,:);





%Raster
Raster_Stim = Raster(StimTrial,:);
Raster_Control = Raster(ControlTrial,:);

%1/ISI
reISI = SegmentData.reISI{1};
reISIseq = SegmentData.reISIseq{1};

reISI_Stim = reISI(StimTrial,:);
reISI_Control = reISI(ControlTrial,:);

reISI_seqStim = reISIseq(StimTrial,:);
reISI_seqControl = reISIseq(ControlTrial,:);








%Spike Count
if useReISI
        %SpikeCount_Stim = reISIseq(StimTrial ,:);
        %SpikeCount_Control = reISIseq(ControlTrial ,:);
        
        SpikeCount_Stim = reISIseq(StimTrial ,:);
        SpikeCount_Control = reISIseq(ControlTrial ,:);
else
    SpikeCount_Stim=SpikeCount(StimTrial,:);
    SpikeCount_Control=SpikeCount(ControlTrial,:);

    SpikeCount_Stim_Total = sum(SpikeCount_Stim,1);
end

%Conincide all the analysis with spike count
%Only use the spike count
PSTH_Stim = SpikeCount_Stim;
PSTH_Control = SpikeCount_Control;


%Average

PSTH_Stim_Mean = nanmean(PSTH_Stim,1);
PSTH_Stim_Sem = nanstd(PSTH_Stim,[],1)./sqrt(sum(~isnan(PSTH_Stim)));



PSTH_Control_Mean = nanmean(PSTH_Control,1);
PSTH_Control_Sem = nanstd(PSTH_Control,[],1)./sqrt(sum(~isnan(PSTH_Control)));

SpikeCountStim_Mean = nanmean(SpikeCount_Stim,1);
SpikeCountStim_Sem = nanstd(SpikeCount_Stim,[],1)./sqrt(sum(~isnan(SpikeCount_Stim)));

SpikeCountControl_Mean = nanmean(SpikeCount_Control,1);
SpikeCountControl_Sem = nanstd(SpikeCount_Control,[],1)./sqrt(sum(~isnan(SpikeCount_Control)));

%Z-scored responses

PSTH_Stim_Mean_z = nanmean((PSTH_Stim-Mean_FR_z)/Std_FR_z,1);
PSTH_Control_Mean_z = nanmean((PSTH_Control-Mean_FR_z)/Std_FR_z,1);


if ismember(FileName,{'Robin010423VIDEOVIEW2.rst'})
    closest = find(MeanPSTH_Time_Mean>=0,1,'first');

    PSTH_Stim_Mean(closest) = mean(PSTH_Stim_Mean(MeanPSTH_Time_Mean<0));
    SpikeCountStim_Mean(closest) = mean(SpikeCountStim_Mean(MeanPSTH_Time_Mean<0));

    PSTH_Stim_Mean_z(closest) =mean(PSTH_Stim_Mean_z(MeanPSTH_Time_Mean<0));

end



%{
%Smoothed psth
sigma = 5; % Standard deviation of the Gaussian kernel
kernelSize = 7; % Size of the kernel (odd number for symmetry)
x = linspace(-3*sigma, 3*sigma, kernelSize);
gaussianKernel = exp(-x.^2 / (2 * sigma^2)) / (sigma * sqrt(2 * pi));

psth_stim_pre = PSTH_Stim_Mean(MeanPSTH_Time_Mean < 0);
psth_stim_post = PSTH_Stim_Mean(MeanPSTH_Time_Mean > 0);


smoothedPSTH_post = conv(psth_stim_pre, gaussianKernel/sum(gaussianKernel), 'same');
smoothedPSTH_pre = conv(psth_stim_post, gaussianKernel/sum(gaussianKernel), 'same');


smoothedPSTH_stim = [smoothedPSTH_pre,smoothedPSTH_post];

figure
plot(smoothedPSTH_stim)
%}
%Calculate response latency
%Use stimulated trials only, compare stimulation period with pre-stimulus
%response
SponTimeIn = [-100,0];
%PreStim_Spon = sum(SpikeCount_Stim(:,MeanPSTH_Time_Mean > SponTimeIn(1) & MeanPSTH_Time_Mean<SponTimeIn(2)),2)/((SponTimeIn(2) - SponTimeIn(1))*1000);
PreStim_Spon = mean(SpikeCount_Stim(:,MeanPSTH_Time_Mean > SponTimeIn(1) & MeanPSTH_Time_Mean<SponTimeIn(2)),2);





%Search for significant time
OptoResTime = MeanPSTH_Time_Mean(MeanPSTH_Time_Mean > 0);
CompareTimeBin =10;% in ms 
ContinousBin = 2;
RespLatency = NaN;
%{
%Loop over every 1 ms 
OnsetPrecision = 1;%in 1ms

for i = 1:OnsetPrecision:OptoResTime(end)-CompareTimeBin
    TimeInv_curr = [i-1,i-1+CompareTimeBin];
    for j = 1:size(Raster_Stim,1)
        SpikeRaster_Curr(j) =  sum(Raster_Stim(j,:) >= TimeInv_curr(1) & Raster_Stim(j,:) <= TimeInv_curr(2))/CompareTimeBin*1000;

    end
    p(i) = ranksum(PreStim_Spon,SpikeRaster_Curr);
    SpikeRate_curr_Store(:,i) = SpikeRaster_Curr';

    if i >= ContinousBin
        cri = 1;
        for j = 1:ContinousBin
            if p(i-j+1) >0.05 | isnan(p(i))
                cri = 0;
            end
        end
        if cri == 1 & isnan(RespLatency)
            RespLatency = i-ContinousBin;

        end

    end




end
%}
%

%Activation significance and strength
StaInterval = 50;%in ms
%{
if ~isnan(ResponseLatency)
    TimeLine = 1:OnsetPrecision:OptoResTime(end);
    SigPeriod = (TimeLine>=RespLatency & TimeLine < RespLatency + StaInterval) & p < 0.05;
    %SpikeRate_curr_Store
    %SpikeCountIn = SpikeCount_Stim(:,MeanPSTH_Time_Mean > 0);

ActivationResponse = mean(SpikeRate_curr_Store(:,SigPeriod),2);
ActivationMean = nanmean(ActivationResponse);
ActivationStrength =  (ActivationMean - max(nanmean(PreStim_Spon),0.0001))/max(nanmean(PreStim_Spon),0.0001);

p_Modulation = ranksum(ActivationResponse,PreStim_Spon);

else
%}
    ActivationResponse = nanmean(SpikeCount_Stim(:,MeanPSTH_Time_Mean >= Interval(1) & MeanPSTH_Time_Mean<=Interval(2)),2);
    ActivationMean = nanmean(ActivationResponse);
    ActivationStrength =  (ActivationMean - max(nanmean(PreStim_Spon),0.0001))/max(nanmean(PreStim_Spon),0.0001);
  %  [h,p_Modulation] = ttest(ActivationResponse,PreStim_Spon);

%end





%Check the post stim effect
%StimMean = nanmean(SpikeCount_Stim(:,MeanPSTH_Time_Mean >= Interval(1) & MeanPSTH_Time_Mean<=Interval(2)),2);
PostStimMean = nanmean(SpikeCount_Stim(:,MeanPSTH_Time_Mean > Interval(2) & MeanPSTH_Time_Mean<=Interval(2)+100),2);
SponMean = nanmean(PreStim_Spon);

PostStimModulation = (nanmean(PostStimMean)-max(SponMean,0.001))/(max(SponMean,0.001));
[h,p_PostModulation] = ttest(PostStimMean,PreStim_Spon);

%Stim effect in the whole stimulation period

%StimModulation = (nanmean(StimMean)-max(SponMean,0.001))/(max(SponMean,0.001));

StimModulation = ActivationMean-SponMean;
StimModulation_z = StimModulation/Std_FR_z;
%p_StimModulation = ranksum(StimMean,PreStim_Spon);
[h,p_StimModulation] = ttest(ActivationResponse,PreStim_Spon);
disp(sprintf('Activation Strength is %.2f, p = %d',StimModulation ,p_StimModulation));

%Compared to control
ControlMean = nanmean(SpikeCount_Control(:,MeanPSTH_Time_Mean >= Interval(1) & MeanPSTH_Time_Mean<=Interval(2)),2);
[h,p_StimControl] = ttest2(ActivationResponse,ControlMean);

%First half,second half trials
%PSTH
FirstHalf = 1:floor(size(PSTH_Stim,1)/2);
SecondHalf = FirstHalf(2) : size(PSTH_Stim,1);
PSTH_Stim_Mean_First = nanmean(PSTH_Stim(FirstHalf,:),1);
PSTH_Stim_Sem_First = nanstd(PSTH_Stim(FirstHalf,:),[],1)./sqrt(sum(~isnan(PSTH_Stim(FirstHalf,:))));

PSTH_Stim_Mean_Second = nanmean(PSTH_Stim(SecondHalf,:),1);
PSTH_Stim_Sem_Second = nanstd(PSTH_Stim(SecondHalf,:),[],1)./sqrt(sum(~isnan(PSTH_Stim(SecondHalf,:))));


PSTH_Stim_Mean_First_z = nanmean((PSTH_Stim(FirstHalf,:)-Mean_FR_z)/Std_FR_z,1);
PSTH_Stim_Mean_Second_z = nanmean((PSTH_Stim(SecondHalf,:)-Mean_FR_z)/Std_FR_z,1);


if ismember(FileName,{'Robin010423VIDEOVIEW2.rst'})
    closest = find(MeanPSTH_Time_Mean>=0,1,'first');

    PSTH_Stim_Mean_First(closest) = mean(PSTH_Stim_Mean_First(MeanPSTH_Time_Mean<0));
    PSTH_Stim_Mean_Second(closest) = mean(PSTH_Stim_Mean_Second(MeanPSTH_Time_Mean<0));

    PSTH_Stim_Mean_First_z(closest) = mean(PSTH_Stim_Mean_First_z(MeanPSTH_Time_Mean<0));
    PSTH_Stim_Mean_Second_z(closest) =  mean(PSTH_Stim_Mean_Second_z(MeanPSTH_Time_Mean<0));

end




%Activation strength
ActivationResponse_FirstHalf = ActivationResponse(FirstHalf);
ActivationResponse_SecondHalf = ActivationResponse(SecondHalf);

PreStim_Spon_FirstHalf = PreStim_Spon(FirstHalf);
PreStim_Spon_SecondHalf = PreStim_Spon(SecondHalf);

ActivationStrength_FirstHalf =  (mean(ActivationResponse_FirstHalf) - nanmean(PreStim_Spon_FirstHalf))/Std_FR_z;%Change into z-score
ActivationStrength_SecondHalf =  (mean(ActivationResponse_SecondHalf) - nanmean(PreStim_Spon_SecondHalf))/Std_FR_z;%Change into z-score



%{
figure
plot(PSTH_Stim_Mean_norm);

%}
%Response latency
for i = 1: length(OptoResTime)
    TimeInv_curr = [OptoResTime(i),OptoResTime(i)+CompareTimeBin];
    %SpikeCount_curr = sum(SpikeCount_Stim(:,MeanPSTH_Time_Mean>=TimeInv_curr(1) & MeanPSTH_Time_Mean < TimeInv_curr(2)),2)/CompareTimeBin;
    SpikeCount_curr = mean(SpikeCount_Stim(:,MeanPSTH_Time_Mean>=TimeInv_curr(1) & MeanPSTH_Time_Mean < TimeInv_curr(2)),2);

    %
    %First judge whether it's an excitation or inhinition
    %
    if StimModulation > 0
        %Excitation 
        
        [h,p(i)] = ttest(PreStim_Spon,SpikeCount_curr,'tail','left');
    else
        [h,p(i)] = ttest(PreStim_Spon,SpikeCount_curr,'tail','right');
    end

    %}
    %}
      %  [h,p(i)] = ttest(PreStim_Spon,SpikeCount_curr,'tail','both');
    
    SpikeCount_curr_store(i) = mean(SpikeCount_curr); 
    
    if i >= ContinousBin
        cri = 1;
        for j = 1:ContinousBin
            if p(i-j+1) >0.05 | isnan(p(i))
                cri = 0;
            end
        end
        if cri == 1 & isnan(RespLatency)
            RespLatency = OptoResTime(i - ContinousBin + 1);

        end

    end

 end
%}
ResponseLatency =RespLatency;
 disp(sprintf('StimLatency is %.2f',RespLatency));

AcDur = sum(p < 0.05);
disp(sprintf('Activationn Duration is %.2f',AcDur));


 SpikeCount_curr_store_rel =  SpikeCount_curr_store - SponMean;

if (p(1)<0.05 & p(2)>0.05) ||(p(1)<0.05 & p(2)<0.05 & SpikeCount_curr_store_rel(1)>0 & SpikeCount_curr_store_rel(2)<0)
    closest = find(MeanPSTH_Time_Mean>=0,1,'first');

     PSTH_Stim_Mean(closest) = mean(PSTH_Stim_Mean(MeanPSTH_Time_Mean<0));
    SpikeCountStim_Mean(closest) = mean(SpikeCountStim_Mean(MeanPSTH_Time_Mean<0));


    PSTH_Stim_Mean_First(closest) = mean(PSTH_Stim_Mean_First(MeanPSTH_Time_Mean<0));
    PSTH_Stim_Mean_Second(closest) = mean(PSTH_Stim_Mean_Second(MeanPSTH_Time_Mean<0));

     PSTH_Stim_Mean_z(closest) =mean(PSTH_Stim_Mean_z(MeanPSTH_Time_Mean<0));

      PSTH_Stim_Mean_First_z(closest) = mean(PSTH_Stim_Mean_First_z(MeanPSTH_Time_Mean<0));
    PSTH_Stim_Mean_Second_z(closest) =  mean(PSTH_Stim_Mean_Second_z(MeanPSTH_Time_Mean<0));




    %ReCheck the post stim effect

%StimModulation = (nanmean(StimMean)-max(SponMean,0.001))/(max(SponMean,0.001));
%StimModulation = nanmean(SpikeCountStim_Mean)-SponMean;
%StimModulation_z = StimModulation/Std_FR_z;
%p_StimModulation = ranksum(StimMean,PreStim_Spon);



end

%{
%Normed responses

sel_laser = MeanPSTH_Time_Mean>=Interval(1) & MeanPSTH_Time_Mean<Interval(2);



%PeakFR = max(abs(PSTH_Stim_Mean(sel_laser)));

%Mean_Spon = mean(mean(SpikeCount_Stim(:,MeanPSTH_Time_Mean > SponTimeIn(1) & MeanPSTH_Time_Mean<SponTimeIn(2)),2));



Mean_Spon = mean(PSTH_Stim_Mean(:,MeanPSTH_Time_Mean > SponTimeIn(1) & MeanPSTH_Time_Mean<SponTimeIn(2)),2);
Mean_Spon_Control = mean(PSTH_Control_Mean(:,MeanPSTH_Time_Mean > SponTimeIn(1) & MeanPSTH_Time_Mean<SponTimeIn(2)),2);


post_sub_pstd_stim = PSTH_Stim_Mean- Mean_Spon;

post_sub_pstd_control = PSTH_Control_Mean- Mean_Spon_Control;

%PeakFR = max(abs(post_sub_pstd_stim(sel_laser)));
PeakFR = max(max(abs(post_sub_pstd_stim)),max(post_sub_pstd_control));
%
PSTH_Stim_Mean_norm = post_sub_pstd_stim/PeakFR ;
PSTH_Control_Mean_norm = post_sub_pstd_control /PeakFR ;
%}


%Norm2 divided by peak first, then subtract the spon(finally use this one)
PeakFR2 =  max(max(PSTH_Stim_Mean),max(PSTH_Control_Mean));


PSTH_Stim_Mean_div = PSTH_Stim_Mean/PeakFR2;
Spon_div = mean(PSTH_Stim_Mean_div(MeanPSTH_Time_Mean > SponTimeIn(1) & MeanPSTH_Time_Mean<SponTimeIn(2)));

PSTH_Stim_Mean_sub = PSTH_Stim_Mean_div-Spon_div;

PSTH_Stim_Mean_norm_raw = PSTH_Stim_Mean_sub;

PSTH_Control_Mean_div = PSTH_Control_Mean/PeakFR2 ;
Spon_Control_div = mean(PSTH_Control_Mean_div(MeanPSTH_Time_Mean > SponTimeIn(1) & MeanPSTH_Time_Mean<SponTimeIn(2)));
PSTH_Control_Mean_sub = PSTH_Control_Mean_div-Spon_Control_div ;

PSTH_Control_Mean_norm_raw = PSTH_Control_Mean_sub;

StimulationStrength_norm = mean(PSTH_Stim_Mean_sub(MeanPSTH_Time_Mean >= Interval(1) & MeanPSTH_Time_Mean<=Interval(2)));



PSTH_Stim_Mean_sub_pre = PSTH_Stim_Mean_sub(MeanPSTH_Time_Mean<0);
PSTH_Stim_Mean_sub_post = PSTH_Stim_Mean_sub(MeanPSTH_Time_Mean>0);

PSTH_Control_Mean_sub_pre = PSTH_Control_Mean_sub(MeanPSTH_Time_Mean<0);
PSTH_Control_Mean_sub_post = PSTH_Control_Mean_sub(MeanPSTH_Time_Mean>0);



sigma = 5; % Standard deviation of the Gaussian kernel
kernelSize = 7; % Size of the kernel (odd number for symmetry)
x = linspace(-3*sigma, 3*sigma, kernelSize);
gaussianKernel = exp(-x.^2 / (2 * sigma^2)) / (sigma * sqrt(2 * pi));

smoothedPSTH_post = conv(PSTH_Stim_Mean_sub_post, gaussianKernel/sum(gaussianKernel), 'same');
smoothedPSTH_pre = conv(PSTH_Stim_Mean_sub_pre  , gaussianKernel/sum(gaussianKernel), 'same');

psth_stim_normed2 = [smoothedPSTH_pre,smoothedPSTH_post];

smoothedPSTH_Control_post = conv(PSTH_Control_Mean_sub_post, gaussianKernel/sum(gaussianKernel), 'same');
smoothedPSTH_Control_pre = conv(PSTH_Control_Mean_sub_pre  , gaussianKernel/sum(gaussianKernel), 'same');

psth_control_normed2 = [smoothedPSTH_Control_pre,smoothedPSTH_Control_post];

%

%Another way to find latency %Not using 
CompareTimeBin2 =10;% in ms 
ContinousBin2 = 5;
CriBin = 2;

Spon_Mean=mean(PreStim_Spon);
Spon_Std = std(PreStim_Spon);
%cum = 0;

for i = 1: length(OptoResTime)
    TimeInv_curr = [OptoResTime(i),OptoResTime(i)+CompareTimeBin2];
    %SpikeCount_curr = sum(SpikeCount_Stim(:,MeanPSTH_Time_Mean>=TimeInv_curr(1) & MeanPSTH_Time_Mean < TimeInv_curr(2)),2)/CompareTimeBin;
   

    SpikeCount_curr = mean(SpikeCount_Stim(:,MeanPSTH_Time_Mean>=TimeInv_curr(1) & MeanPSTH_Time_Mean < TimeInv_curr(2)),2);

    %First judge whether it's an excitation or inhinition
    %{
    if StimModulation > 0
        %Excitation 
        
        [h,p(i)] = ttest(PreStim_Spon,SpikeCount_curr,'tail','left');
    else
        [h,p(i)] = ttest(PreStim_Spon,SpikeCount_curr,'tail','right');
    end

        %}
       % [h,p(i)] = ttest(PreStim_Spon,SpikeCount_curr,'tail','both');

    h(i) = abs(mean(SpikeCount_curr)-Spon_Mean) >3*Spon_Std;
    Time_i(i)= TimeInv_curr(1);

   

end
length_select = floor(length(h)/ContinousBin2 )*ContinousBin2 ;
h_reshape = reshape(h(1:length_select),ContinousBin2 ,length_select/ContinousBin2 );
Time_i_reshape = reshape(Time_i(1:length_select),ContinousBin2 ,length_select/ContinousBin2 );
tests_latency = sum(h_reshape,1);

if ~isempty(find(tests_latency,1,'first'))
    Latency2 = Time_i_reshape(1,find(tests_latency,1,'first'));
else
    Latency2 = NaN;
end

disp(sprintf('Latency2=%2.2f',Latency2 ));
%}

%Another way to define the p of the stimulation effect
CompareTimeBin3 = 50;%in ms
CompareTimeBin4 = 20;%in ms

for i = 1: length(OptoResTime)
    TimeInv_curr = [OptoResTime(i),OptoResTime(i)+CompareTimeBin3];
    %SpikeCount_curr = sum(SpikeCount_Stim(:,MeanPSTH_Time_Mean>=TimeInv_curr(1) & MeanPSTH_Time_Mean < TimeInv_curr(2)),2)/CompareTimeBin;
    SpikeCount_tmp= mean(SpikeCount_Stim(:,MeanPSTH_Time_Mean>=TimeInv_curr(1) & MeanPSTH_Time_Mean < TimeInv_curr(2)),2);
     
    SpikeCount_Store(i) = mean(SpikeCount_tmp);


    TimeInv_curr2 = [OptoResTime(i),OptoResTime(i)+CompareTimeBin4];
    SpikeCount_tmp2= mean(SpikeCount_Stim(:,MeanPSTH_Time_Mean>=TimeInv_curr2(1) & MeanPSTH_Time_Mean < TimeInv_curr2(2)),2);
    SpikeCount_Store2(i) = mean(SpikeCount_tmp2);



    %First judge whether it's an excitation or inhinition
    %{
    if StimModulation > 0
        %Excitation 
        
        [h,p(i)] = ttest(PreStim_Spon,SpikeCount_curr,'tail','left');
    else
        [h,p(i)] = ttest(PreStim_Spon,SpikeCount_curr,'tail','right');
    end

        %}
        [h,p_judge(i)] = ttest(PreStim_Spon,SpikeCount_tmp,'tail','both');

        [h,p_judge2(i)] = ttest(PreStim_Spon,SpikeCount_tmp2,'tail','both');




        

    
   % SpikeCount_curr_store(i) = mean(SpikeCount_curr); 
    

 end
%}

p_use_cri = min(p_judge); %Use this cri50 finally

Strength_sig = SpikeCount_Store(find(p_judge2<0.05,1,'first'))-mean(PreStim_Spon);


p_use_cri20 = min(p_judge2);




%Another criterium
for i = 1: 50: length(OptoResTime)

     TimeInv_curr = [OptoResTime(i),OptoResTime(i)+50];
    %SpikeCount_curr = sum(SpikeCount_Stim(:,MeanPSTH_Time_Mean>=TimeInv_curr(1) & MeanPSTH_Time_Mean < TimeInv_curr(2)),2)/CompareTimeBin;
    SpikeCount_tmp= mean(SpikeCount_Stim(:,MeanPSTH_Time_Mean>=TimeInv_curr(1) & MeanPSTH_Time_Mean < TimeInv_curr(2)),2);
     
    SpikeCount_Store3(i) = mean(SpikeCount_tmp);
[h,p_judge3(i)] = ttest(PreStim_Spon, SpikeCount_tmp,'tail','both');

   

end
p_whole_cri50 = min(p_judge3);
disp(sprintf('min p value for 50ms bin %2.2f',p_whole_cri50));

%keyboard

%p_whole_cri50

%Another latency
%
RespLatency = NaN;
for i = 1: length(OptoResTime)
    TimeInv_curr = [OptoResTime(i),OptoResTime(i)+CompareTimeBin];
    %SpikeCount_curr = sum(SpikeCount_Stim(:,MeanPSTH_Time_Mean>=TimeInv_curr(1) & MeanPSTH_Time_Mean < TimeInv_curr(2)),2)/CompareTimeBin;
    SpikeCount_curr = mean(SpikeCount_Stim(:,MeanPSTH_Time_Mean>=TimeInv_curr(1) & MeanPSTH_Time_Mean < TimeInv_curr(2)),2);

    %First judge whether it's an excitation or inhinition
    %
    if Strength_sig  > 0
        %Excitation 
        
        [h,p(i)] = ttest(PreStim_Spon,SpikeCount_curr,'tail','left');
    else
        [h,p(i)] = ttest(PreStim_Spon,SpikeCount_curr,'tail','right');
    end

       
       % [h,p(i)] = ttest(PreStim_Spon,SpikeCount_curr,'tail','both');
    
    SpikeCount_curr_store(i) = mean(SpikeCount_curr); 
    
    if i >= ContinousBin
        cri = 1;
        for j = 1:ContinousBin
            if p(i-j+1) >0.05 | isnan(p(i))
                cri = 0;
            end
        end
        if cri == 1 & isnan(RespLatency)
            RespLatency = OptoResTime(i - ContinousBin + 1);

        end

    end

 end

Latency3 = RespLatency;
%}
%{
ResponseLatency = LatencyWMW(MeanPSTH_Time_Mean,SpikeCount_Stim,SpikeCount_Control,ContinousBin);

%Activation significance and strength
ActivationMean = nanmean(SpikeCount_Stim(:,MeanPSTH_Time_Mean >= Interval(1) & MeanPSTH_Time_Mean<=Interval(2)),2);
ControlMean = nanmean(SpikeCount_Control(:,MeanPSTH_Time_Mean >= Interval(1) & MeanPSTH_Time_Mean<=Interval(2)),2);

p_Modulation = ranksum(ActivationMean,ControlMean);
ActivationStrength = (nanmean(ActivationMean) - nanmean(ControlMean))/nanmean(ControlMean);

disp(sprintf('Activationn Strength is %.2f, p = %d',ActivationStrength,p_Modulation));

%Activation duration
AcDur =  ActivationDuration(BinWidth,SpikeCount_Stim(:,MeanPSTH_Time_Mean>0),SpikeCount_Control(:,MeanPSTH_Time_Mean>0));

disp(sprintf('Activationn Duration is %.2f',AcDur));
%}
%%
%Plots 
if ShowFigureFlag
    
FigureStartNum=90;
FigureIndex=1;
    
figtitlestr{FigureIndex}='PSTH_StimVSControl';
fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[50,100, 1000,600],'Name',figtitlestr{FigureIndex});

subplot(2,1,1)

AverageRegion = nanmean(RelativeTimeMarker);
area(AverageRegion,[max([SpikeCountControl_Mean+SpikeCountControl_Sem,SpikeCountStim_Mean+SpikeCountStim_Sem])+20,max([SpikeCountControl_Mean+SpikeCountControl_Sem,SpikeCountStim_Mean+SpikeCountStim_Sem])+20],...
    'FaceColor','k','FaceAlpha',.1,'EdgeAlpha',.1,'HandleVisibility','off');
hold on
%{


shadedErrorBar(MeanPSTH_Time_Mean,PSTH_Control_Mean,PSTH_Control_Sem,'lineprops',{[0,0,0]},'transparent',1,'patchSaturation',0.3)

shadedErrorBar(MeanPSTH_Time_Mean,PSTH_Stim_Mean,PSTH_Stim_Sem,'lineprops',{[1,0,0]},'transparent',1,'patchSaturation',0.3)
%}
if useReISI == 1
    plot(MeanPSTH_Time_Mean,PSTH_Control_Mean,'color',[0,0,0]);

    plot(MeanPSTH_Time_Mean,PSTH_Stim_Mean,'color',[1,0,0]);
else

    shadedErrorBar(MeanPSTH_Time_Mean,SpikeCountControl_Mean,SpikeCountControl_Sem,'lineprops',{[0,0,0]},'transparent',1,'patchSaturation',0.3)
   %bar(SpikeCount_Stim_Total,'r');


   shadedErrorBar(MeanPSTH_Time_Mean,SpikeCountStim_Mean,SpikeCountStim_Sem,'lineprops',{[1,0,0]},'transparent',1,'patchSaturation',0.3)
end
xlim([-PlotStartOffset,PlotStopOffset+nanmean(RelativeTimeMarker(:,2))]);

%ylim([0,150]);
set(findall(gca, 'Type', 'Line'),'LineWidth',3);
 if useReISI == 1
    ylabel('Spike Frequency Function(1/ISI)');
 else
    ylabel('Firing Rate(Hz)')
 end
set(gca,'LineWidth',3,'FontSize',20,'FontWeight','Bold');



title(sprintf("Laser Power : %1.3f mW; Latency: %.1f ms",unique(stimPower),ResponseLatency));
legend(["Sham","Stim"])
box off


%{
%Figure2 SpikeRaterPlot
FigureStartNum = FigureStartNum+1;
FigureIndex = FigureIndex;

figtitlestr{FigureIndex}='Raster_StimVSControl';
fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[70,100, 1000,400],'Name',figtitlestr{FigureIndex});
%}
StartLoc = 1;
LineWidth = 2;
Barlength = 1;

subplot(2,1,2)
AverageRegion = nanmean(RelativeTimeMarker);
area(AverageRegion,[size(Raster,1)+5,size(Raster,1)+5],...
    'FaceColor','k','FaceAlpha',.1,'EdgeAlpha',.1,'HandleVisibility','off');

hold on

%Raster For Stim
SpikeRaterPlot(Raster_Stim,StartLoc,Barlength,'r',LineWidth);

StartLoc = size(Raster_Stim,1)+5;
%Raster For Control

SpikeRaterPlot(Raster_Control,StartLoc,Barlength,'k',LineWidth);
xlim([-PlotStartOffset,PlotStopOffset+nanmean(RelativeTimeMarker(:,2))]);
xlabel('Time from Opto-stim On');   

 set(gca,'LineWidth',3,'FontSize',20,'FontWeight','Bold');
 yticks([]);
 box off

 %% Figure 2, first half and second half
FigureStartNum = FigureStartNum + 1;
FigureIndex = FigureIndex + 1;
    
figtitlestr{FigureIndex}='PSTH_FirstHalf_SecondHalf';
fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[70,100, 1000,600],'Name',figtitlestr{FigureIndex});

shadedErrorBar(MeanPSTH_Time_Mean,PSTH_Stim_Mean_First ,PSTH_Stim_Sem_First ,'lineprops',{'#785EF0'},'transparent',1,'patchSaturation',0.3)
hold on
shadedErrorBar(MeanPSTH_Time_Mean,PSTH_Stim_Mean_Second ,PSTH_Stim_Sem_Second ,'lineprops',{'#FFB000'},'transparent',1,'patchSaturation',0.3)
set(findall(gca, 'Type', 'Line'),'LineWidth',3);
xlabel('Time from Opto-stim On');   
ylabel('Firing Rate(Hz)');
set(gca,'LineWidth',3,'FontSize',20,'FontWeight','Bold');
legend({'FirstHalf','SecondHalf'})



if PlotExample 
    %Only plot stim with color relative to excitation or inhibition
    FigureStartNum = FigureStartNum + 1;
    FigureIndex = FigureIndex + 1;
    
    figtitlestr{FigureIndex}='Optostim_Modulation';
    fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[70,100, 1000,600],'Name',figtitlestr{FigureIndex});

    if StimModulation>0 & p_StimModulation<0.05
    rcolor = [1,0,0];
    elseif StimModulation<0 & p_StimModulation<0.05
    rcolor = [0,0,1];
    else
        rcolor = [0.2,0.2,0.2];

    end


    subplot(2,1,2)

AverageRegion = nanmean(RelativeTimeMarker);
area(AverageRegion,[max([SpikeCountControl_Mean+SpikeCountControl_Sem,SpikeCountStim_Mean+SpikeCountStim_Sem])+5,max([SpikeCountControl_Mean+SpikeCountControl_Sem,SpikeCountStim_Mean+SpikeCountStim_Sem])+5],...
    'FaceColor','k','FaceAlpha',.1,'EdgeAlpha',.1,'HandleVisibility','off');
hold on
%{


shadedErrorBar(MeanPSTH_Time_Mean,PSTH_Control_Mean,PSTH_Control_Sem,'lineprops',{[0,0,0]},'transparent',1,'patchSaturation',0.3)

shadedErrorBar(MeanPSTH_Time_Mean,PSTH_Stim_Mean,PSTH_Stim_Sem,'lineprops',{[1,0,0]},'transparent',1,'patchSaturation',0.3)
%}



   %bar(SpikeCount_Stim_Total,'r');


  shadedErrorBar(MeanPSTH_Time_Mean,PSTH_Stim_Mean,PSTH_Stim_Sem,'lineprops',{rcolor},'transparent',1,'patchSaturation',0.3)
 %bar(SpikeCountStim_Mean*BinWidth/1000,'FaceColor',rcolor,'LineWidth',3);

xlim([-PlotStartOffset,PlotStopOffset+nanmean(RelativeTimeMarker(:,2))]);

%ylim([0,150]);
set(findall(gca, 'Type', 'Line'),'LineWidth',3);
 
    ylabel('Firing Rate(Hz)');

set(gca,'LineWidth',3,'FontSize',20,'FontWeight','Bold');



%title(sprintf("Laser Power : %1.3f mW; Latency: %.1f ms",unique(stimPower),ResponseLatency));
%legend(["Sham","Stim"])
box off


%{
%Figure2 SpikeRaterPlot
FigureStartNum = FigureStartNum+1;
FigureIndex = FigureIndex;

figtitlestr{FigureIndex}='Raster_StimVSControl';
fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[70,100, 1000,400],'Name',figtitlestr{FigureIndex});
%}
StartLoc = 1;
LineWidth = 2;
Barlength = 1;



subplot(2,1,1)
AverageRegion = nanmean(RelativeTimeMarker);
area(AverageRegion,[size(Raster,1)+5,size(Raster,1)+5],...
    'FaceColor','k','FaceAlpha',.1,'EdgeAlpha',.1,'HandleVisibility','off');

hold on

%Raster For Stim

SpikeRaterPlot(Raster_Stim,StartLoc,Barlength,rcolor,LineWidth);


%Raster For Control
xlim([-PlotStartOffset,PlotStopOffset+nanmean(RelativeTimeMarker(:,2))]);
xlabel('Time from Opto-stim On');   

 set(gca,'LineWidth',3,'FontSize',20,'FontWeight','Bold');
 yticks([]);
 box off
 axis off




end

end %   End of if showFigureFlag

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
OutputFigureName=sprintf('%s_%s_N%s_C%s_%CID_%NID',MonkeyName,string(RecordDateOriginal),string(NeuronNum),string(ChannelNum),string(ChannelID),string(ClusterID));





ExistFlag=0;
%Load the old file if exist
if exist(OutputFileName)

    load(OutputFileName);
    %Load OutputData into the memory
    if isfield(OutputData,'FreeViewNeuralRed')
    ExistFlag=1;
    end
end



%{
%FileName=strcat(FileName,'_RearchingTime.mat');
if SaveFlag
save(FileName,'OutputData');
disp('Seaching Time Data Exported');
else
    disp('No data export');
end
%}
%OutputData For current analysis
OutputData_New=[];

OutputData_New.FreeViewNeuralRed.TrialType='GoodTrials';

%%%%%%%%%Save overall scene response characteristics%%%%%%%%%%%%%%%%%%%%%%%
%Organize the data


NameStr={'PSTH_Time','PSTH_Stim_Mean','PSTH_Control_Mean','AcStrength','ResponseLatency',...
    'PostModulationStrength','PostModulation_p','StimModulation','StimModulation_p',...
    'PSTH_FirstHalf','PSTH_SecondHalf',...
    'AS_First','AS_Second',...
    'PSTH_Stim_Mean_z','PSTH_Control_Mean_z','PSTH_Stim_Mean_First_z','PSTH_Stim_Mean_Second_z',...
    'StimModulation_z',...
    'PSTH_Stim_Mean_norm','PSTH_Control_Mean_norm','FR_Mean','p_StimControl',...
    'AcDur','psth_stim_normed2','psth_control_normed2','p_cri50','p_cri20','Strength_sig','p_whole_cri50',...
    'Latency3','Latency2',...
    'Raster_Stim','Raster_Control','StimulationStrength_norm',...
    'ChannelID','SpikeChannel'
    };

DataLib={MeanPSTH_Time_Mean,SpikeCountStim_Mean,SpikeCountControl_Mean,ActivationStrength,ResponseLatency,...
    PostStimModulation,p_PostModulation,StimModulation,p_StimModulation,...
    PSTH_Stim_Mean_First,PSTH_Stim_Mean_Second,...
    ActivationStrength_FirstHalf,ActivationStrength_SecondHalf,...
    PSTH_Stim_Mean_z,PSTH_Control_Mean_z,PSTH_Stim_Mean_First_z,PSTH_Stim_Mean_Second_z,...
    StimModulation_z,...
    PSTH_Stim_Mean_norm_raw,PSTH_Control_Mean_norm_raw,FR_Mean,p_StimControl,...
    AcDur,psth_stim_normed2,psth_control_normed2,p_use_cri,p_use_cri20,Strength_sig,p_whole_cri50,...
    Latency3,Latency2,...
    Raster_Stim,Raster_Control,StimulationStrength_norm,...
    ChannelID,SpikeChannel

};



DataStamp=containers.Map(NameStr,DataLib);





OutputData_New.FreeViewNeuralRed.Task=TaskType;
OutputData_New.FreeViewNeuralRed.TaskCode=TaskCode;
OutputData_New.FreeViewNeuralRed.DataStamp=DataStamp;

OutputData_New.FreeViewNeuralRed.StartTrial =StartTrial;
%{
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
    
    for jj=1:TotalFigure
        if ~isempty(fig{jj})
        
         OutputFigureName_Final=strcat(OutputFigureName,figtitlestr{jj},'.jpg');
         saveas(fig{jj},OutputFigureName_Final);
        end
        
            
        
    end

    disp('Figures have been exported to the neuron folder');

end

%}





if ExistFlag
%Compare the old one with the new one
Task_Old=OutputData.FreeViewNeuralRed.TaskCode;

     if Task_Old~=TaskCode
        %Not the same task, add the new dataset into the old one
        OutputData(numel(OutputData)+1)=OutputData_New;
        %{
     elseif ismember(FileName,{'Adams062124VIDEOVIEW1'})
         %May save different task segment
         if StartTrial~=OutputData.FreeViewNeural.StartTrial

             if StartTrial > 160
                 OutputFileName=sprintf('%s_%s_N%s_C%s.mat',MonkeyName,string(RecordDateOriginal),string(2),string(ChannelNum));
                 OutputData.FreeViewNeural=OutputData_New.FreeViewNeural;

             elseif StartTrial > 360
                  OutputFileName=sprintf('%s_%s_N%s_C%s.mat',MonkeyName,string(RecordDateOriginal),string(3),string(ChannelNum));
                  OutputData.FreeViewNeural=OutputData_New.FreeViewNeural;

             end


         end
     elseif ismember(FileName,{'Robin030524VIDEOVIEW1'})
        % if StartTrial~=OutputData.FreeViewNeural.StartTrial

             if StartTrial > 80 && StartTrial <= 160
                 OutputFileName=sprintf('%s_%s_N%s_C%s.mat',MonkeyName,string(RecordDateOriginal),string(2),string(ChannelNum));
                 OutputData.FreeViewNeural=OutputData_New.FreeViewNeural;
 
             elseif StartTrial > 240 && StartTrial <= 160
                  OutputFileName=sprintf('%s_%s_N%s_C%s.mat',MonkeyName,string(RecordDateOriginal),string(3),string(ChannelNum));
                  OutputData.FreeViewNeural=OutputData_New.FreeViewNeural;

             end


       %  end
        %}     
     else
         %If the same task, replace the old one with the new one
         %OutputData=OutputData_New;
         OutputData.FreeViewNeuralRed=OutputData_New.FreeViewNeuralRed;
         
    
    
     end
else
    %Output the current file
   OutputData.FreeViewNeuralRed=OutputData_New.FreeViewNeuralRed; 
end



if OutputFlag
    cd(OutPath);
    save(OutputFileName,'OutputData');
    disp('Data Exported');
else
    disp('No data export');
end  


end %End loop of the spike channel
end
function ObjectIndex=ReproduceFromEvent(event,code)
for i=1:size(event,1)
    
    ObjectIndex{i}=event(i,ismember(event(i,:),code));
    
  
    
end
numeach=cellfun(@numel,ObjectIndex);
nummax=max(numeach);
if nummax==1
    ObjectIndex=cell2mat(ObjectIndex);

    
end


end