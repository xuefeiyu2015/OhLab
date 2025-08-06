function FreeViewNeuroResponse_Pop(Data, StartTrial, EndTrial,ShowFigureFlag,OutputFlag)
ECode;

TaskType=Data.TaskType;
DataPath= Data.Path;
FileName=Data.FileName;
TaskCode=Data.TaskCode;

Selection=StartTrial:EndTrial;%Select the trials according to the pannel
useReISI = 0;
      
SegEvent = [ON_STIM];
SegEventEnd = [OFF_STIM];
SegEventStr={'StimOn'};

PlotStartOffset=[100];%in ms
PlotStopOffset=[100];%in ms

SepChannel = 32;

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
if ismember(FileName,{'Robin061423VIDEOVIEW1.rst','Robin062223VIDEOVIEW1.rst','Robin071323VIDEOVIEW1.rst','Robin072023VIDEOVIEW1.rst','Robin072523VIDEOVIEW1.rst'})
    stimPower = RescuePower(FileName,stimPower_ori);
    
else
    stimPower = stimPower_ori;
end

 stimPower=stimPower(Selection,:);



%Parameter to define the sliding window to measure firing rates
StepSize=5;%in ms
BinWidth=5;%in ms


%Load event data
EventChannel=Data.EventChannel(Selection,:);
EventTimeChannel=Data.EventTimeChannel(Selection,:);
%EventBin=Data.EventBin(Selection,:);


%% Load the spike channel
%ExcludeChanel = [5,10,38];

%ExcludeChanel =[5,13,29];
ExcludeChanel=[];

SpikeChannelWhole=Data.SpikeChannel;

selectChannel = ~ismember(SpikeChannelWhole,ExcludeChanel);

SpikeChannelWhole = SpikeChannelWhole(selectChannel);
SpikeTimeWhole=Data.SpikeTimeData(selectChannel,:,Selection);
ChannelIDWhole=Data.ChannelID(selectChannel);
ClusterIDWhole=Data.ClusterID(selectChannel);

MultiChannel=0;
if length(SpikeChannelWhole)>1
    MultiChannel=1;
    
end

SignificanceChannel = NaN*ones(length(SpikeChannelWhole),1);
SigFastChannel = NaN*ones(length(SpikeChannelWhole),1);
SigSlowChannel = NaN*ones(length(SpikeChannelWhole),1);


ResponseLatency_store = NaN*ones(length(SpikeChannelWhole),1);
ModType =  NaN*ones(length(SpikeChannelWhole),1);

for spk=1:length(SpikeChannelWhole)
    SpikeChannel=SpikeChannelWhole(spk);
    
        SpikeTime=squeeze(SpikeTimeWhole(spk,:,:))';
        %{
    if SpikeChannel<4
        SpikeTime = SpikeTime;%Correct for the delay of online channel
    end
        %}
    SpikeChannel_curr=SpikeChannelWhole(spk);
    ClusterID=ClusterIDWhole(spk);
    ChannelID=ChannelIDWhole(spk);

    
FR_All = (sum(~isnan(SpikeTime),2)-2)./SpikeTime(end)*1000;

FR_Mean = mean(FR_All);


if SpikeChannel_curr<4
    SpikeTime = SpikeTime-3;%Compensate for the online delay;
end


SegmentData=EventSegTool(SegEvent,SegEventEnd,PlotStartOffset,PlotStopOffset,EventChannel,EventTimeChannel,SpikeTime,StepSize,BinWidth);

%Retrive data from the field
%Total average for each event
%MeanFR=SegmentData.MeanFR; %[Mean,sem]%Averall mean and sem
%MeanFRRaw=SegmentData.MeanFRRaw{1};%=AverageFiringRate_Event;

%Relative time bin for each event
RelativeTimeMarker=SegmentData.RelativeTimeMarker{1};%=EventTimeMarker;
%Absolute time marker for each event
%AbosoluteTimeMarker=SegmentData.AbsoluteTimeMarker{1};


TimeMarker=SegmentData.RelativeTimeMarker{1}; %[Mean,sem]%Averall mean and sem

Interval = nanmean(TimeMarker);

%For Raster
%Raster=SegmentData.RasterAfterAlign{1};

SpikeCount = SegmentData.PSTH_Event_SpikeCount{1}/BinWidth*1000;



%PSTH for each event
MeanPSTH=SegmentData.MeanPSTH;%=[PSTH_Event_Time',PSTH_Event'];
MeanPSTH_Time=MeanPSTH{1};
MeanPSTH_FR=MeanPSTH{2};

MeanPSTH_Time_Mean = nanmean(MeanPSTH_Time);

%Seperate into stim trial and sham trial
StimType=cell2mat(ReproduceFromEvent(EventChannel,[STIM1,SHAM0])');

StimTrial = StimType == STIM1;
ControlTrial = ~StimTrial;

%Mean and std for z score
%Mean_FR_z = mean(MeanFRRaw(ControlTrial));
%Std_FR_z = std(MeanFRRaw(ControlTrial));

%Mean_FR_z = mean(FR_All(ControlTrial));
%Std_FR_z= std(FR_All(ControlTrial));






%{
if Std_FR_z <0.0000001
    Mean_FR_z = mean(MeanFRRaw);
    Std_FR_z = std(MeanFRRaw);

end
%}
%PSTH_Stim = MeanPSTH_FR(StimTrial(stimOn==1),:);
%PSTH_Control = MeanPSTH_FR(ControlTrial(stimOn==1),:);

%PSTH_Stim = MeanPSTH_FR(StimTrial,:);
%PSTH_Control = MeanPSTH_FR(ControlTrial,:);




%Raster
%Raster_Stim = Raster(StimTrial,:);
%Raster_Control = Raster(ControlTrial,:);


%Spike Count

    SpikeCount_Stim=SpikeCount(StimTrial,:);
    SpikeCount_Control=SpikeCount(ControlTrial,:);

    %SpikeCount_Stim_Total = sum(SpikeCount_Stim,1);


%Conincide all the analysis with spike count

PSTH_Stim = SpikeCount_Stim;
PSTH_Control = SpikeCount_Control;


%Average


SpikeCountStim_Mean = nanmean(SpikeCount_Stim,1);
SpikeCountStim_Sem = nanstd(SpikeCount_Stim,[],1)./sqrt(sum(~isnan(SpikeCount_Stim)));

SpikeCountControl_Mean = nanmean(SpikeCount_Control,1);
SpikeCountControl_Sem = nanstd(SpikeCount_Control,[],1)./sqrt(sum(~isnan(SpikeCount_Control)));

%% Save for plotting
%PSTH_raw_Plot_Mean{spk} = SpikeCountStim_Mean ;
%PSTH_raw_Plot_Sem{spk} = SpikeCountStim_Sem ;

%PSTH_raw_Plot_ControlMean{spk} = SpikeCountControl_Mean ;
%PSTH_raw_Plot_ControlSem{spk} = SpikeCountControl_Sem ;
%%
PSTH_Stim_Mean = nanmean(PSTH_Stim,1);
PSTH_Stim_Sem = nanstd(PSTH_Stim,[],1)./sqrt(sum(~isnan(PSTH_Stim)));

PSTH_Control_Mean = nanmean(PSTH_Control,1);
PSTH_Control_Sem = nanstd(PSTH_Control,[],1)./sqrt(sum(~isnan(PSTH_Control)));

%Z-Score
Pre_Spon_all = SpikeCount_Stim(MeanPSTH_Time_Mean<0);

%FR_Mean_Spon = mean(Pre_Spon_all);
%FR_Std_Spon = std(Pre_Spon_all);

FR_Mean_Spon = mean(mean(SpikeCount_Stim,2));
FR_Std_Spon = std(mean(SpikeCount_Stim,2));

psth_z_scored = (SpikeCount_Stim-FR_Mean_Spon)/FR_Std_Spon;
psth_z_scored_mean = mean(psth_z_scored,1);


psth_z_pre = psth_z_scored_mean(MeanPSTH_Time_Mean<0);

psth_z_post = psth_z_scored_mean(MeanPSTH_Time_Mean>=0);

sigma = 10; % Standard deviation of the Gaussian kernel
kernelSize = 7; % Size of the kernel (odd number for symmetry)
x = linspace(-3*sigma, 3*sigma, kernelSize);
gaussianKernel = exp(-x.^2 / (2 * sigma^2)) / (sigma * sqrt(2 * pi));

smoothedPSTH_zpost = conv2(psth_z_post, gaussianKernel/sum(gaussianKernel), 'same');
smoothedPSTH_zpre = conv2(psth_z_pre, gaussianKernel/sum(gaussianKernel), 'same');

psth_z(spk,:) = [smoothedPSTH_zpre,smoothedPSTH_zpost ];


%Normed responses

sel_laser = MeanPSTH_Time_Mean>=Interval(1) & MeanPSTH_Time_Mean<Interval(2);


%Mean_Spon = mean(mean(SpikeCount_Stim(:,MeanPSTH_Time_Mean > SponTimeIn(1) & MeanPSTH_Time_Mean<SponTimeIn(2)),2));

SponTimeIn = [-100,0];

Mean_Spon = mean(PSTH_Stim_Mean(:,MeanPSTH_Time_Mean > SponTimeIn(1) & MeanPSTH_Time_Mean<SponTimeIn(2)),2);
Mean_Spon_Control = mean(PSTH_Control_Mean(:,MeanPSTH_Time_Mean > SponTimeIn(1) & MeanPSTH_Time_Mean<SponTimeIn(2)),2);


post_sub_pstd_stim = PSTH_Stim_Mean- Mean_Spon;

post_sub_pstd_control = PSTH_Control_Mean- Mean_Spon_Control;

%PeakFR = max(abs(post_sub_pstd_stim));

PeakFR = max(abs(post_sub_pstd_stim()));

%PeakFR = max(abs(post_sub_pstd_stim(sel_laser)));

PSTH_Stim_Mean_norm(spk,:) = post_sub_pstd_stim/PeakFR ;
%PSTH_Control_Mean_norm = post_sub_pstd_control /PeakFR ;

%Norm with another way
%PSTH_Stim_Mean=PSTH_Stim_Mean-PSTH_Control_Mean;

peak_W = max(PSTH_Stim_Mean);

div_psth_stim = PSTH_Stim_Mean/peak_W;
spon_Stim = mean(div_psth_stim(:,MeanPSTH_Time_Mean > SponTimeIn(1) & MeanPSTH_Time_Mean<SponTimeIn(2)),2);

PSTH_Stim_Mean_norm2_tmp = div_psth_stim -spon_Stim;

%Control



PSTH_Stim_Mean_norm2_tmp_pre = PSTH_Stim_Mean_norm2_tmp(MeanPSTH_Time_Mean < 0);
PSTH_Stim_Mean_norm2_tmp_post = PSTH_Stim_Mean_norm2_tmp(MeanPSTH_Time_Mean >0);


%Spon_after(spk) = mean(PSTH_Stim_Mean_norm2(spk,MeanPSTH_Time_Mean > SponTimeIn(1) & MeanPSTH_Time_Mean<SponTimeIn(2)),2);

%Mean_FR_z = mean(FR_All(ControlTrial));
%Std_FR_z= std(FR_All(ControlTrial));

Mean_FR_z = FR_Mean ;
Std_FR_z= std(FR_All);



%PSTH_Stim_Mean_z(spk,:) = nanmean((PSTH_Stim-Mean_FR_z)/Std_FR_z,1);
%PSTH_Control_Mean_z = nanmean((PSTH_Control-Mean_FR_z)/Std_FR_z,1);


%Smooth the PSTH for presentation 
% Define the Gaussian kernel parameters
sigma = 5; % Standard deviation of the Gaussian kernel
kernelSize = 7; % Size of the kernel (odd number for symmetry)
x = linspace(-3*sigma, 3*sigma, kernelSize);
gaussianKernel = exp(-x.^2 / (2 * sigma^2)) / (sigma * sqrt(2 * pi));

smoothedPSTH_post = conv(PSTH_Stim_Mean_norm2_tmp_post, gaussianKernel/sum(gaussianKernel), 'same');
smoothedPSTH_pre = conv(PSTH_Stim_Mean_norm2_tmp_pre , gaussianKernel/sum(gaussianKernel), 'same');



PSTH_Stim_Mean_norm2(spk,:)=[smoothedPSTH_pre,smoothedPSTH_post];
%PSTH_Stim_Mean_norm2(spk,:)=PSTH_Stim_Mean_norm2_tmp;




PSTH_Stim_Mean_norm1_tmp_pre = PSTH_Stim_Mean_norm(spk,MeanPSTH_Time_Mean < 0);
PSTH_Stim_Mean_norm1_tmp_post = PSTH_Stim_Mean_norm(spk,MeanPSTH_Time_Mean >0);

smoothedPSTH_post1 = conv(PSTH_Stim_Mean_norm1_tmp_post, gaussianKernel/sum(gaussianKernel), 'same');
smoothedPSTH_pre1 = conv(PSTH_Stim_Mean_norm1_tmp_pre , gaussianKernel/sum(gaussianKernel), 'same');

PSTH_Stim_Mean_norm1(spk,:)=[smoothedPSTH_pre1,smoothedPSTH_post1];


%Statistical test
%Search for significant time
OptoResTime = MeanPSTH_Time_Mean(MeanPSTH_Time_Mean > 0);
CompareTimeBin =10;% in ms 
ContinousBin = 3;
RespLatency = NaN;


SponTimeIn = [-100,0];
%PreStim_Spon = sum(SpikeCount_Stim(:,MeanPSTH_Time_Mean > SponTimeIn(1) & MeanPSTH_Time_Mean<SponTimeIn(2)),2)/((SponTimeIn(2) - SponTimeIn(1))*1000);
PreStim_Spon = mean(SpikeCount_Stim(:,MeanPSTH_Time_Mean > SponTimeIn(1) & MeanPSTH_Time_Mean<SponTimeIn(2)),2);

%Cri1
for i = 1: 50: length(OptoResTime)

     TimeInv_curr = [OptoResTime(i),OptoResTime(i)+50];
    %SpikeCount_curr = sum(SpikeCount_Stim(:,MeanPSTH_Time_Mean>=TimeInv_curr(1) & MeanPSTH_Time_Mean < TimeInv_curr(2)),2)/CompareTimeBin;
    SpikeCount_tmp= mean(SpikeCount_Stim(:,MeanPSTH_Time_Mean>=TimeInv_curr(1) & MeanPSTH_Time_Mean < TimeInv_curr(2)),2);

     SpikeCount_Control_tmp= mean(SpikeCount_Control(:,MeanPSTH_Time_Mean>=TimeInv_curr(1) & MeanPSTH_Time_Mean < TimeInv_curr(2)),2);

     
     
    SpikeCount_Store3(i) = mean(SpikeCount_tmp);
%[h,p_judge3(i)] = ttest(PreStim_Spon, SpikeCount_tmp,'tail','both');
[h,p_judge3(i)] = ttest2(SpikeCount_Control_tmp, SpikeCount_tmp,'tail','both');

   

end
p_whole_cri50(spk) = min(p_judge3);
%disp(sprintf('min p value for 50ms bin %2.2f',p_whole_cri50(spk)));

%Cri2
StimMean = nanmean(SpikeCount_Stim(:,MeanPSTH_Time_Mean >= Interval(1) & MeanPSTH_Time_Mean<=Interval(2)),2);
PostStimMean = nanmean(SpikeCount_Stim(:,MeanPSTH_Time_Mean > Interval(2) & MeanPSTH_Time_Mean<=Interval(2)+100),2);
SponMean = nanmean(PreStim_Spon);

ControlMean = nanmean(SpikeCount_Control(:,MeanPSTH_Time_Mean >= Interval(1) & MeanPSTH_Time_Mean<=Interval(2)),2);

[h,p_StimModulation(spk)] = ttest(StimMean,PreStim_Spon);
%[h,p_StimModulation(spk)] = ttest(StimMean,ControlMean );%Use this finally
%{
if spk == 4
    keyboard;
end
%}

StimModulation(spk) = nanmean(StimMean)-SponMean;
%StimModulation(spk) = nanmean(StimMean)-nanmean(ControlMean);




%Cri3
%PostStimModulation = (nanmean(PostStimMean)-max(SponMean,0.001))/(max(SponMean,0.001));
PostStimModulation = (nanmean(PostStimMean)-max(nanmean(ControlMean),0.001))/(max(nanmean(ControlMean),0.001));
[h,p_PostModulation(spk)] = ttest2(PostStimMean,ControlMean);

%Latency


for i = 1: length(OptoResTime)
    TimeInv_curr = [OptoResTime(i),OptoResTime(i)+CompareTimeBin];
    %SpikeCount_curr = sum(SpikeCount_Stim(:,MeanPSTH_Time_Mean>=TimeInv_curr(1) & MeanPSTH_Time_Mean < TimeInv_curr(2)),2)/CompareTimeBin;
    SpikeCount_curr = mean(SpikeCount_Stim(:,MeanPSTH_Time_Mean>=TimeInv_curr(1) & MeanPSTH_Time_Mean < TimeInv_curr(2)),2);

    SpikeCount_Control_curr = mean(SpikeCount_Control(:,MeanPSTH_Time_Mean>=TimeInv_curr(1) & MeanPSTH_Time_Mean < TimeInv_curr(2)),2);
    %First judge whether it's an excitation or inhinition
    %
    if StimModulation(spk) > 0
        %Excitation 
        
        [h,p(i)] = ttest(PreStim_Spon,SpikeCount_curr,'tail','left');
    else
        [h,p(i)] = ttest(PreStim_Spon,SpikeCount_curr,'tail','right');
    end

        %}
      % [h,p(i)] = ttest(SpikeCount_Control_curr,SpikeCount_curr,'tail','both');
    
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
  %  SpikeCount_Store2(i) = mean(SpikeCount_tmp2);



    %First judge whether it's an excitation or inhinition
    %
    if StimModulation(spk) > 0
        %Excitation 
        
        [h,p_judge(i)] = ttest(PreStim_Spon,SpikeCount_tmp,'tail','left');
        [h,p_judge2(i)] = ttest(PreStim_Spon,SpikeCount_tmp2,'tail','left');
    else
        [h,p_judge(i)] = ttest(PreStim_Spon,SpikeCount_tmp,'tail','right');
        [h,p_judge2(i)] = ttest(PreStim_Spon,SpikeCount_tmp2,'tail','right');
    end

       %{
        [h,p_judge(i)] = ttest(PreStim_Spon,SpikeCount_tmp,'tail','both');

        [h,p_judge2(i)] = ttest(PreStim_Spon,SpikeCount_tmp2,'tail','both');
       %}


p_use_cri(spk) = p_StimModulation(spk); %Needs correction; not right

%Strength_sig = SpikeCount_Store(find(p_judge2<0.05,1,'first'))-mean(PreStim_Spon);


p_use_cri20 = min(p_judge2);

%disp(sprintf('p value: %1.2f', p_StimModulation(spk)));
        

    
   % SpikeCount_curr_store(i) = mean(SpikeCount_curr); 
    

 end
%}



%Sig_m = (p_use_cri(spk)<0.05 | p_StimModulation(spk)<0.05 | p_PostModulation(spk)<0.05)&(~isnan(RespLatency));

%Sig_m = (p_StimModulation(spk)<0.05 | p_PostModulation(spk)<0.05)&(~isnan(RespLatency));

Sig_m = (p_StimModulation(spk)<0.05)&(~isnan(RespLatency));

Sig_m_fast = ResponseLatency<=20 & Sig_m;
Sig_m_slow = ResponseLatency> 20 & Sig_m;

ModType(spk) = sign(StimModulation(spk));



ResponseLatency_store(spk) = ResponseLatency;

SignificanceChannel(spk) = Sig_m==1;
SigFastExciteChannel(spk) = Sig_m_fast == 1 & ModType(spk)==1;
SigSlowExciteChannel(spk) = Sig_m_slow == 1 & ModType(spk)==1;
SigInhibitChannel(spk) = SignificanceChannel(spk) & ModType(spk)==-1;



end %End of the loop of channel number

SignificanceChannel = reshape(SignificanceChannel, length(SignificanceChannel),1);
SigFastExciteChannel = reshape(SigFastExciteChannel , length(SigFastExciteChannel),1);
SigSlowExciteChannel = reshape(SigSlowExciteChannel , length(SigSlowExciteChannel),1);
SigInhibitChannel = reshape(SigInhibitChannel, length(SigInhibitChannel),1);

%A=[ResponseLatency_store,p_whole_cri50'<0.05,p_StimModulation'<0.05,p_PostModulation'<0.05,SignificanceChannel];

SpikeChannel_Significant = SpikeChannelWhole(SignificanceChannel==1);
disp('Significant spike channel');
disp(SpikeChannel_Significant);

%SpikeChannel_Sigfast =SpikeChannelWhole(SigFastChannel == 1);

SpikeChannel_Sigfast = SpikeChannelWhole(SigFastExciteChannel==1);
disp('Significant fast excitation spike channel');
disp(SpikeChannel_Sigfast);

SpikeChannel_SigSlow = SpikeChannelWhole(SigSlowExciteChannel==1);
disp('Significant slow excitation spike channel');
disp(SpikeChannel_SigSlow);

SpikeChannel_SigInhibite = SpikeChannelWhole(SigInhibitChannel == 1);
disp('Significant slow excitation spike channel');
disp(SpikeChannel_SigInhibite);


%Plot the psth out with the heatmap
if ShowFigureFlag
    
FigureStartNum=90;
FigureIndex=1;
    
figtitlestr{FigureIndex}='PSTH_AcrossChannel';
fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[50,100, 1000,600],'Name',figtitlestr{FigureIndex});

% create blue to red colormap:

cmap = nan(128,3);

r = [0.8 0 0];

w = [0.8 0.8 0.8];

b = [0 0 0.8];

for ii = 1:3

    cmap(:,ii) = [ fliplr(linspace(w(ii), b(ii), 64)),fliplr(linspace(r(ii), w(ii), 64))];

end

cmap(65,:) = [];


%ChannelIDWhole
if max(ChannelIDWhole)>SepChannel
    %Separate into 2
    PSTH_Stim_Mean_norm_First = PSTH_Stim_Mean_norm2(ChannelIDWhole<=SepChannel,:);
    PSTH_Stim_Mean_norm_Second = PSTH_Stim_Mean_norm2(ChannelIDWhole>SepChannel,:);

    subplot(1,2,1);
    %First part of channels
    imagesc(PSTH_Stim_Mean_norm_First);

    ChannelID_first=ChannelIDWhole(ChannelIDWhole<=SepChannel);

    TimePoint = [-100,0,100,200];
for i =1:length(TimePoint)
    TimeBin_i = MeanPSTH_Time_Mean-TimePoint(i)==StepSize/2;
    try
    TimeBin(i) = find(TimeBin_i==1,1,'first');
    catch
      TimeBin(i) = NaN;
    end

end
xticks(TimeBin(~isnan(TimeBin)));
xticklabels( string(TimePoint(~isnan(TimeBin) )));
%yticks(ChannelID_first);
%yticks(ClusterIDWhole(ChannelIDWhole<=24));
set(gca, 'CLim', [-1, 1]);
%costom_color_map = [0,0,1;0,1,0;1,0,0];

%opts.cmapRwB = cmap;

colormap(cmap);
colorbar;
hold on 
%plot([0,0],[1,size(PSTH_Stim_Mean_norm,1)],'--k','LineWidth',3);
line([0, 0], [1, size(PSTH_Stim_Mean_norm_First,1)], 'Color', 'k', 'LineWidth', 2);

title(sprintf('From channel 1 to %d',SepChannel));

 subplot(1,2,2);
    %Second part of channels
    imagesc(PSTH_Stim_Mean_norm_Second);

    TimePoint = [-100,0,100,200];
for i =1:length(TimePoint)
    TimeBin_i = MeanPSTH_Time_Mean-TimePoint(i)==StepSize/2;
    try
    TimeBin(i) = find(TimeBin_i==1,1,'first');
    catch
      TimeBin(i) = NaN;
    end

end
xticks(TimeBin(~isnan(TimeBin)));
xticklabels( string(TimePoint(~isnan(TimeBin) )));
%yticks(ClusterIDWhole(ChannelIDWhole>24));
set(gca, 'CLim', [-1, 1]);
%costom_color_map = [0,0,1;0,1,0;1,0,0];

%opts.cmapRwB = cmap;

colormap(cmap);
colorbar;
hold on 
%plot([0,0],[1,size(PSTH_Stim_Mean_norm,1)],'--k','LineWidth',3);
line([0, 0], [1, size(PSTH_Stim_Mean_norm_Second,1)], 'Color', 'k', 'LineWidth', 2);
title(sprintf('From channel %d to 56',SepChannel));




else


imagesc(PSTH_Stim_Mean_norm2);
%imagesc(psth_z);


%xticks(1:50:size(PSTH_Stim_Mean_norm,2));
TimePoint = [-100,0,100,200];
for i =1:length(TimePoint)
    TimeBin_i = MeanPSTH_Time_Mean-TimePoint(i)==StepSize/2;
    try
    TimeBin(i) = find(TimeBin_i==1,1,'first');
    catch
      TimeBin(i) = NaN;
    end

end
xticks(TimeBin(~isnan(TimeBin)));
xticklabels( string(TimePoint(~isnan(TimeBin) )));
yticks(1:size(PSTH_Stim_Mean_norm,1));
%OptoChannel=3;
%yticklabels(ChannelIDWhole - OptoChannel);
set(gca, 'CLim', [-1, 1]);
%costom_color_map = [0,0,1;0,1,0;1,0,0];
% create blue to red colormap:

cmap = nan(128,3);

r = [0.8 0 0];

w = [0.8 0.8 0.8];

b = [0 0 0.8];

for ii = 1:3

    cmap(:,ii) = [ fliplr(linspace(w(ii), b(ii), 64)),fliplr(linspace(r(ii), w(ii), 64))];

end

cmap(65,:) = [];

%opts.cmapRwB = cmap;

colormap(cmap);
colorbar;
hold on 
%plot([0,0],[1,size(PSTH_Stim_Mean_norm,1)],'--k','LineWidth',3);
line([0, 0], [1, size(PSTH_Stim_Mean_norm,1)], 'Color', 'k', 'LineWidth', 2);
end 
%}
%{
for i = 1:length(SpikeChannelWhole)
    subplot(length(SpikeChannelWhole),1,i);
    plot(MeanPSTH_Time_Mean,PSTH_Stim_Mean_norm(i,:),'-r','LineWidth',3);
    box off;

end
%}

FigureStartNum=FigureStartNum+1;
FigureIndex=FigureIndex+1;

figtitlestr{FigureIndex}='PSTH_AcrossChannel_Raw';
fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[70,100, 600,1000],'Name',figtitlestr{FigureIndex});

psth = (flipud(PSTH_Stim_Mean_norm2))*10;
PannelDist = 5;
%AddInterval = 0:PannelDist:size(psth,1)*PannelDist-1;

%psth = psth+AddInterval';
%psth = flipud(psth);


if max(ChannelIDWhole)>SepChannel
psth_first = psth(flipud(ChannelIDWhole)<=SepChannel,:);
AddInterval_first = 0:PannelDist:size(psth_first,1)*PannelDist-1;
psth_first = psth_first+AddInterval_first';


psth_second = psth(flipud(ChannelIDWhole)>SepChannel,:);
AddInterval_second = 0:PannelDist:size(psth_second,1)*PannelDist-1;
psth_second = psth_second+AddInterval_second';


%SignificanceChannel_first = SignificanceChannel(ChannelIDWhole<=24);
%SignificanceChannel_second = SignificanceChannel(ChannelIDWhole>24);

SignificanceChannelfast_first = SigFastExciteChannel(ChannelIDWhole<=SepChannel);
SignificanceChannelfast_second =SigFastExciteChannel(ChannelIDWhole>SepChannel);

SignificanceChannelslow_first = SigSlowExciteChannel(ChannelIDWhole<=SepChannel);
SignificanceChannelslow_second =SigSlowExciteChannel(ChannelIDWhole>SepChannel);

SigInhibitChannel_first = SigInhibitChannel(ChannelIDWhole<=SepChannel);
SigInhibitChannel_second = SigInhibitChannel(ChannelIDWhole>SepChannel);



subplot(1,2,1);

AverageRegion = nanmean(RelativeTimeMarker);
area(AverageRegion,[max(psth_first(length(SpikeChannelWhole(ChannelIDWhole<=SepChannel)),:))+PannelDist,max(psth_first(length(SpikeChannelWhole(ChannelIDWhole<=SepChannel)),:))+PannelDist],...
    'FaceColor','k','FaceAlpha',.1,'EdgeAlpha',.1,'HandleVisibility','off');
hold on
%}
time_plot =repmat(MeanPSTH_Time_Mean,size(psth_first,1),1);
plot(time_plot',psth_first','-k','LineWidth',2);
hold on

%plot(time_plot(flipud(SigFastChannel==1),:)',psth(flipud(SigFastChannel==1),:)','-r','LineWidth',3);
plot(time_plot(flipud(SignificanceChannelfast_first==1),:)',psth_first(flipud(SignificanceChannelfast_first==1),:)','-r','LineWidth',3);
plot(time_plot(flipud(SignificanceChannelslow_first==1),:)',psth_first(flipud(SignificanceChannelslow_first==1),:)','-y','LineWidth',3);



plot(time_plot(flipud(SigInhibitChannel_first==1),:)',psth_first(flipud(SigInhibitChannel_first==1),:)','-b','LineWidth',3);


axis off
title(sprintf('Channel from 1 to %d',SepChannel));

subplot(1,2,2);

AverageRegion = nanmean(RelativeTimeMarker);
area(AverageRegion,[max(psth_second(length(SpikeChannelWhole(ChannelIDWhole>SepChannel)),:))+PannelDist,max(psth_second(length(SpikeChannelWhole(ChannelIDWhole>SepChannel)),:))+PannelDist],...
    'FaceColor','k','FaceAlpha',.1,'EdgeAlpha',.1,'HandleVisibility','off');
hold on
%}
time_plot =repmat(MeanPSTH_Time_Mean,size(psth_second,1),1);
plot(time_plot',psth_second','-k','LineWidth',2);
hold on

%plot(time_plot(flipud(SigFastChannel==1),:)',psth(flipud(SigFastChannel==1),:)','-r','LineWidth',3);
plot(time_plot(flipud(SignificanceChannelfast_second==1),:)',psth_second(flipud(SignificanceChannelfast_second==1),:)','-r','LineWidth',3);
plot(time_plot(flipud(SignificanceChannelslow_second==1),:)',psth_second(flipud(SignificanceChannelslow_second==1),:)','-y','LineWidth',3);



plot(time_plot(flipud(SigInhibitChannel_second==1),:)',psth_second(flipud(SigInhibitChannel_second==1),:)','-b','LineWidth',3);

axis off
title(sprintf('Channel from %d to 56',SepChannel+1));





else

AddInterval = 0:PannelDist:size(psth,1)*PannelDist-1;

psth = psth+AddInterval';
%psth = flipud(psth);

%
AverageRegion = nanmean(RelativeTimeMarker);
area(AverageRegion,[max(psth(length(SpikeChannelWhole),:))+PannelDist,max(psth(length(SpikeChannelWhole),:))+PannelDist],...
    'FaceColor','k','FaceAlpha',.1,'EdgeAlpha',.1,'HandleVisibility','off');
hold on
%}
time_plot =repmat(MeanPSTH_Time_Mean,size(psth,1),1);
plot(time_plot',psth','-k','LineWidth',2);
hold on

%plot(time_plot(flipud(SigFastChannel==1),:)',psth(flipud(SigFastChannel==1),:)','-r','LineWidth',3);



%plot(time_plot(flipud(SignificanceChannel==1& ModType ==1),:)',psth(flipud(SignificanceChannel==1& ModType ==1),:)','-r','LineWidth',3);
%plot(time_plot(flipud(SignificanceChannel==1& ModType ==-1),:)',psth(flipud(SignificanceChannel==1& ModType ==-1),:)','-b','LineWidth',3);



plot(time_plot(flipud(SigFastExciteChannel),:)',psth(flipud(SigFastExciteChannel),:)','-r','LineWidth',3);
plot(time_plot(flipud(SigSlowExciteChannel),:)',psth(flipud(SigSlowExciteChannel),:)','-y','LineWidth',3);
plot(time_plot(flipud(SigInhibitChannel),:)',psth(flipud(SigInhibitChannel),:)','-b','LineWidth',3);



axis off
end %End of else above 24
%{
channelorder = 1: length(SpikeChannelWhole);
yticklabels(string(channelorder));
%}
%keyboard

%{
for spk=1:length(SpikeChannelWhole)
    subplot(length(SpikeChannelWhole),1,spk);
    AverageRegion = nanmean(RelativeTimeMarker);
area(AverageRegion,[max([PSTH_raw_Plot_ControlMean{spk}+PSTH_raw_Plot_ControlSem{spk},PSTH_raw_Plot_Mean{spk}+PSTH_raw_Plot_Sem{spk}])+20,max([PSTH_raw_Plot_ControlMean{spk}+PSTH_raw_Plot_ControlSem{spk},PSTH_raw_Plot_Mean{spk}+PSTH_raw_Plot_Sem{spk}])+20],...
    'FaceColor','k','FaceAlpha',.1,'EdgeAlpha',.1,'HandleVisibility','off');
hold on
 %  shadedErrorBar(MeanPSTH_Time_Mean,PSTH_raw_Plot_ControlMean{spk},PSTH_raw_Plot_ControlSem{spk},'lineprops',{[0,0,0]},'transparent',1,'patchSaturation',0.3)
   %bar(SpikeCount_Stim_Total,'r');


   shadedErrorBar(MeanPSTH_Time_Mean,PSTH_raw_Plot_Mean{spk},PSTH_raw_Plot_Sem{spk},'lineprops',{[1,0,0]},'transparent',1,'patchSaturation',0.3)

xlim([-PlotStartOffset,PlotStopOffset+nanmean(RelativeTimeMarker(:,2))]);
box off
axis off

end
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0, 1, 1]);
set(gcf, 'Units', 'Inches');
pos = get(gcf, 'Position');
set(gcf, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Inches', 'PaperSize', [pos(3), pos(4)]);

%}
end