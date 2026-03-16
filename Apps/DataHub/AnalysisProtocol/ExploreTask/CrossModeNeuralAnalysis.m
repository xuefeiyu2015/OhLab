function CrossModeNeuralAnalysis(Data, StartTrial, EndTrial,ShowFigureFlag,OutputFlag);
%Function to Analysis the neural responses difference in different mode of
%exploration
ECode;
TaskType=Data.TaskType;
TaskCode=Data.TaskCode;
DataPath= Data.Path;
FileName=Data.FileName;

Selection=StartTrial:EndTrial;%Select the trials according to the pannel
%Abnormal check %Especially For Robin081120FORAGEXPLO4.rst
CheckCode=SCENE_ON;%This code should appear once in each trial
AbnormalTrials=AbnormalTrialCheck(Data,CheckCode);
Selection=Selection(~ismember(Selection,AbnormalTrials));
TrialNum=length(Selection);
%User defined parameter for the neural/behavioral analysis
%TargetOnOffset=100;%in ms, the time selection bin after target onset; 
%StartDelay=-50;%in ms, the time selection bin after fixation offset;

%For Plot of PSTH aligned on Each fixation on and each fixation Off
PlotFixationOn_StartOffset=100;
PlotFixationOn_StopOffset=300;%Duration of plotting

PlotFixationOff_StartOffset=300;
PlotFixationOff_StopOffset=100;%Duration of plotting

PlotSaccadeOn_StartOffset=100;
PlotSaccadeOn_StopOffset=300;

PlotSaccadeOff_StartOffset=300;
PlotSaccadeOff_StopOffset=100;


%For segmentation of  the neural responses 
SegEvent=[FIXOFF,TGTCD];%Event Flow of interest
SegEventEnd=[TGTCD,TGTOFF ];%Event end
SegEventStr={'Explore','Fractal'};

PlotStartOffset=[700,200];%in ms
%PlotStartOffset=[1700];%in ms
PlotStopOffset=[0,0];%in ms

EventTypeFlag=[2,2];%1:Number same with the trial number, 2:Number same with the target number

%Spon Event
SponEvent=[SCENE_OFF];
SponPreTimeBuffer=-500;
SponEventEnd=NaN;
SponStartOffset=0;
SponStopOffset=1000;

%Parameter to define the sliding window to measure firing rates
StepSize=5;%in ms
BinWidth=15;%in ms

%Plot Raster Flag
PlotRasterFlag=0;

%Criterias for selecting saccades
V_Threshold=40;
ContinueBin=5;
ISI_Threshold=20;

%ResponseLatencyThreshold=3;

%Some user defined parameter for visualization
%Color Definition for each saccade nums
Color={'r','g','b','m',[0.4,0.4,0],[0.28,0.78,0.89],'k',[0.6,0.1,1],[0.6,0.1,0.4],[0.3,0.5,1],[0.6,0.5,0.1],[0,0.8,0.4]};
ColorForEachMode=[0.83 0.14 0.14
             1.00 0.54 0.00
             0.47 0.25 0.80
             0.25 0.80 0.54];
AssociateBlock=[10,11,0,1];
Blockstr={'FixLoc,FixSceneSeq','FixLoc,RandSceneSeq','RandLoc,FixSceneSeq','RandLoc,RandSceneSeq'};
%Get File Information
num = regexp(Data.FileName, '\d+', 'match');
RecordDateOriginal=cell2mat(num(1));
RecordDate=datetime(RecordDateOriginal,'InputFormat','MMddyy');
%To get the monkey name
RemainFiles=erase(Data.FileName,RecordDateOriginal);
Delimiter=find(isstrprop(RemainFiles,'upper')==1);
MonkeyName=RemainFiles(Delimiter(1):Delimiter(2)-1);

%BlockSepTrial=unique(SceneLib(:,1));
BlockType=Data.BlockType(:,:);
if RecordDate<datestr('04,03,2021') && TaskCode==FORAGING_EXPLORE_TASK_ITL
    fix_tmp=BlockType(:,3)>0;
    BlockType(fix_tmp,2)=BlockType(fix_tmp,2)-10;
    
end

BlockType=BlockType(Selection,2);



%Correct early block type
%Before 03,18,2020, I only use 0 to indicate fix scene sequence and 1 to
%indicate random scene sequence, both with fixed location sequence
if RecordDate<datestr('03,17,2020')
    
   BlockType(BlockType==0)=10;
    BlockType(BlockType==1)=11;
   
    
end

%{
BlockType(BlockType==10)=11;
BlockType(BlockType==0)=1;
%}
%
BlockType(BlockType==10)=NaN;
BlockType(BlockType==0)=NaN;
%}
%BlockType(BlockType(:,3)>0,3)=999;



%GetEyeData
EyeData=Data.EyeDataRaw;
EyeChannel_X_L=EyeData(Selection,1);
EyeChannel_Y_L=EyeData(Selection,2);
%EyeChannel_X_R=EyeData(Selection,3);
%EyeChannel_Y_R=EyeData(Selection,4);


EyeBinWidth=Data.EyeBinWidth;
EyeCalInfo=Data.EyeCalInfo;%Each row: X_Gain,X_Bias; Y_Gain; Y_bias; X=(raw+X_bias)*X_Gain
EyeCalInfo_L=EyeCalInfo(1,:);

%Load event data
EventChannel=Data.EventChannel(Selection,:);
EventTimeChannel=Data.EventTimeChannel(Selection,:);
EventBin=Data.EventBin(Selection,:);
%Get the information about the extra FP
%
if RecordDate<datestr('10,28,2020') && strcmp(MonkeyName,'Robin')
%Recalibrate eyes according to the extra FP
ParaLib=Data.ParaLib;

if(ismember('ExtraFixPointID,Style,Size,Loc_X,Loc_Y',keys(ParaLib)))
ExtraFP_Info=ParaLib('ExtraFixPointID,Style,Size,Loc_X,Loc_Y');
ExtraFP_Loc=ExtraFP_Info(Selection,end-1:end);
Ref_Loc=ExtraFP_Loc;

%Get the information about the extra FP
%ExtraFPOnset=FindOutTime(EventChannel,EventTimeChannel,EXTRAFP);
ExtraEyeIn=FindOutTime(EventChannel,EventBin,EXTRAFIXSTARTCD);
ExtraEyeOff=FindOutTime(EventChannel,EventBin,EXTRAFIXOFF);

EyeInBin_X=SelectEyeInterval(EyeChannel_X_L,ExtraEyeIn,ExtraEyeOff);

EyeInBin_Y=SelectEyeInterval(EyeChannel_Y_L,ExtraEyeIn,ExtraEyeOff);


else
    
FP_Loc=Data.FixLoc(Selection,:);
Ref_Loc=FP_Loc;
%Get the information about the extra FP
%ExtraFPOnset=FindOutTime(EventChannel,EventTimeChannel,EXTRAFP);
EyeIn=FindOutTime(EventChannel,EventBin,FIXSTARTCD );
EyeOff=FindOutTime(EventChannel,EventBin,FIXOFF);

EyeInBin_X=SelectEyeInterval(EyeChannel_X_L,EyeIn,EyeOff);

EyeInBin_Y=SelectEyeInterval(EyeChannel_Y_L,EyeIn,EyeOff);
end
EyeInBin_X=EyeInBin_X{1};
EyeInBin_Y=EyeInBin_Y{1};
EyeDataForCal={EyeInBin_X,EyeInBin_Y};

EyeDataRaw={EyeChannel_X_L,EyeChannel_Y_L};
CalPlotCheck=0;%Examine the calibration results
CaledEye=EyeCalibration(EyeDataRaw,EyeDataForCal,EyeCalInfo_L,Ref_Loc,CalPlotCheck);


EyeChannel_X_L=CaledEye.EyeTrace{1};
EyeChannel_Y_L=CaledEye.EyeTrace{2};
end
%}

%Load the spike channel
SpikeChannel=Data.SpikeChannel;
SpikeTime=squeeze(Data.SpikeTimeData(:,:,Selection))';


%Get scene and fractals information
SceneName=Data.SceneIndex(Selection,3);
TargetName=Data.TargetIndex(Selection,2:end,2);
TargetIndex=Data.TargetIndex(Selection,2:end,1);
TargetLoc=Data.TargetLoc(Selection,2:end,:);

%TargetNum=Data.TargetNum(Selection,2:end);
TargetNum=sum(~isnan(TargetName),2);
TargetSize=Data.TargetSize(Selection,2:end);
TargetPolarAngle=Data.TargetPolarLoc(Selection,2:end,1);
TargetPolarEccentricity=Data.TargetPolarLoc(Selection,2:end,2);

ChoiceOutcome=Data.ChoiceOutcome(Selection,2:end);

%Select the trials at lease has one case of advance eye in

AdvanceInTargetNum=sum(EventChannel==ADVEYE_IN,2);
%AdvanceInTargetTrial=AdvanceInTargetNum>0;
AdvanceInTargetTrial=AdvanceInTargetNum==TargetNum;

strAdvanceIn=fprintf('There are %d Advance-Eye-In trials\n',sum(AdvanceInTargetTrial));


%Select Searching time
 FixOffTime=FindOutTime(EventChannel,EventTimeChannel,FIXOFF);
 %TargetOnTime=FixOffTime+300;
 AdvanceEyeInTime=FindOutTime(EventChannel,EventTimeChannel,ADVEYE_IN );
 
 TargetOnTime=FindOutTime(EventChannel,EventTimeChannel,TGTCD);
 OutScreenTime=FindOutTime(EventChannel,EventTimeChannel,EYEOUTSCREEN);
 
 FPOnTime=FindOutTime(EventChannel,EventTimeChannel,SHOWFIXCD);
 %screen out the searching trials
 FixOffTime=FixOffTime(AdvanceInTargetTrial,:);
 OutScreenTime=OutScreenTime(AdvanceInTargetTrial,:);
 FPOnTime=FPOnTime(AdvanceInTargetTrial,:);
 
 FixOffTime=ScreenOutOverScreen(FixOffTime,OutScreenTime);
 
 
 TargetOnTime=TargetOnTime(AdvanceInTargetTrial,:);
 %Remove the col
 NaNCol_Check=sum(isnan(TargetOnTime),1);

 TargetOnTime(:,NaNCol_Check==size(TargetOnTime,1))=[];
 
 SceneName=SceneName(AdvanceInTargetTrial,:);
 
 EyeChannel_X_L=EyeChannel_X_L(AdvanceInTargetTrial);
 EyeChannel_Y_L=EyeChannel_Y_L(AdvanceInTargetTrial);
 
 TargetLoc=TargetLoc(AdvanceInTargetTrial,:,:);
 TargetPolarAngle=TargetPolarAngle(AdvanceInTargetTrial,:);
 
 TargetNum=TargetNum(AdvanceInTargetTrial);
 BlockType=BlockType(AdvanceInTargetTrial);
 
 [UniqueBlockType,orderindex,reprojectorderindex]=unique(BlockType,'row');
 ValidBlocks=~isnan(UniqueBlockType);
 UniqueBlockType=UniqueBlockType(ValidBlocks);
 NumberOfBlocks=length(UniqueBlockType);
%NumberOfBlocks=max(reprojectorderindex);
 
 %If the target onset time and fixofftime dimension are not aligned 
 if size(TargetOnTime,2)~=size(FixOffTime,2)
     FixOffTime=FixOffTime(:,1:size(TargetOnTime,2));
     disp('Warning!FixOffTime Corrected!!')
 end
     
 
SearchingTime= (TargetOnTime-FixOffTime);%/1000;%Change into s for presentation


%Segregate spike data from Fixation onset to target onset
SegmentData=EventSegTool(SegEvent,SegEventEnd,PlotStartOffset,PlotStopOffset,EventChannel,EventTimeChannel,SpikeTime,StepSize,BinWidth);

%Retrive data from the field
%Total average for each event
MeanFR=SegmentData.MeanFR; %[Mean,sem]%Averall mean and sem
MeanFRRaw=SegmentData.MeanFRRaw;%=AverageFiringRate_Event;

%PSTH(Spike Density Function) for each event
MeanPSTH=SegmentData.MeanPSTH;%=[PSTH_Event_Time',PSTH_Event'];
MeanPSTH_Time=MeanPSTH{1,1};
MeanPSTH_FR=MeanPSTH{1,2};

%Raster and raw spike count PSTH
AlignedRaster=SegmentData.RasterAfterAlign{1};
PSTH_Event_SpikeCount=SegmentData.PSTH_Event_SpikeCount{1};


%Relative time bin for each event
RelativeTimeMarker=SegmentData.RelativeTimeMarker{1};%=EventTimeMarker;
%Absolute time marker for each event
AbosoluteTimeMarker=SegmentData.AbsoluteTimeMarker{1};

%Get the PSTH after object onset
MeanPSTH_Time_Fractal=MeanPSTH{2,1};
MeanPSTH_FR_Fractal=MeanPSTH{2,2};


%Reorganize parameters into line library
%Rearrange target into lines
%Target
TargetIndex_Line=TargetIndex;
TargetIndex_Line=reshape(TargetIndex_Line,numel(TargetIndex_Line),1);

NoNNaNSelector=~isnan(TargetIndex_Line);
TargetIndex_Line=TargetIndex_Line(~isnan(TargetIndex_Line));

TargetName_Line=TargetName;
TargetName_Line=reshape(TargetName_Line,numel(TargetName_Line),1);
TargetName_Line=TargetName_Line(~isnan(TargetName_Line));

%TargetNumLine
TargetNum_Line=repmat(TargetNum,1,max(TargetNum));
TargetNum_Line=reshape(TargetNum_Line,numel(TargetNum_Line),1);
TargetNum_Line=TargetNum_Line(NoNNaNSelector);


%TargetOrderLine
TargetOrder=1:max(TargetNum);
TotalNum=length(TargetNum);
TargetOrder_Mat=repmat(TargetOrder,TotalNum,1);
TargetOrder_Line=reshape(TargetOrder_Mat,numel(TargetOrder_Mat),1);
TargetOrder_Line=TargetOrder_Line(NoNNaNSelector);


%Rearrange block type into lines
BlockType_Line=repmat(BlockType,1,max(TargetNum));
BlockType_Line=reshape(BlockType_Line,numel(BlockType_Line),1);
BlockType_Line=BlockType_Line(NoNNaNSelector);

%Rearrange searching time into lines
SearchingTime_Line=reshape(SearchingTime,numel(SearchingTime),1);
SearchingTime_Line=SearchingTime_Line(NoNNaNSelector);

%Rearrange_trialIndex into lines
TrialIndex=[1:size(SearchingTime,1)]';
TrialIndex_Line=repmat(TrialIndex,1,max(TargetNum));
TrialIndex_Line=reshape(TrialIndex_Line,numel(TrialIndex_Line),1);
TrialIndex_Line=TrialIndex_Line(NoNNaNSelector);

%FR for spon
EventTimeChannelForSpon=EventTimeChannel-SponPreTimeBuffer;

SegmentDataSpon=EventSegTool(SponEvent,SponEventEnd,SponStartOffset,SponStopOffset,EventChannel,EventTimeChannelForSpon,SpikeTime,StepSize,BinWidth);

Mean_Spon=SegmentDataSpon.MeanFR;
Mean_Spon_ForEachTrial=SegmentDataSpon.MeanFRRaw{1};

Mean_Spon_ForEachTrial_Line=repmat(Mean_Spon_ForEachTrial,1,max(TargetNum));
Mean_Spon_ForEachTrial_Line=reshape(Mean_Spon_ForEachTrial_Line,numel(Mean_Spon_ForEachTrial_Line),1);
Mean_Spon_ForEachTrial_Line=Mean_Spon_ForEachTrial_Line(NoNNaNSelector);



%Select out eye trace for each target
SegEventTime=FindOutTime(EventChannel,EventTimeChannel,SegEvent(1));%For this protocol, there is only one event
SegEventEndTime=FindOutTime(EventChannel,EventTimeChannel,SegEventEnd(1));%For this protocol, there is only one event


StartTimeBin=floor((SegEventTime-PlotStartOffset(1))/EyeBinWidth);         
EndTimeBin=ceil((SegEventEndTime+PlotStopOffset(1))/EyeBinWidth);
         
Eye_Explore_X=SelectEyeInterval(EyeChannel_X_L,StartTimeBin,EndTimeBin)';


Eye_Explore_Y=SelectEyeInterval(EyeChannel_Y_L,StartTimeBin,EndTimeBin)';


%Eye_Explore_Line=ReorganizeData(Eye_Explore_X,Eye_Explore_Y);
%Eye_Explore_Line=Eye_Explore_Line(NoNNaNSelector,:);
EyeChannel_X_Explore=ReorganizeData(Eye_Explore_X);
EyeChannel_Y_Explore=ReorganizeData(Eye_Explore_Y);

EyeChannel_X_Explore=EyeChannel_X_Explore(NoNNaNSelector,:);
EyeChannel_Y_Explore=EyeChannel_Y_Explore(NoNNaNSelector,:);

%Using thresholding to select saccades

Saccades=SelectSaccade(EyeChannel_X_Explore,EyeChannel_Y_Explore,EyeBinWidth,V_Threshold,ContinueBin,ISI_Threshold);

SaccadeAngle=arrayfun(@(x)  x.SaccadeAngle,Saccades,'UniformOutput' ,0)';
SaccadeAmplitude=arrayfun(@(x)  x.SaccadeAmplitude,Saccades,'UniformOutput' ,0)';

SaccadeStartTime=arrayfun(@(x)  x.SaccadeStartTime,Saccades,'UniformOutput' ,0)';
SaccadeEndTime=arrayfun(@(x)  x.SaccadeEndTime,Saccades,'UniformOutput' ,0)';

SaccadeStartPoint=arrayfun(@(x)  x.SaccadeStartPoint,Saccades,'UniformOutput' ,0)';
SaccadeEndPoint=arrayfun(@(x)  x.SaccadeEndPoint,Saccades,'UniformOutput' ,0)';

%Screen out the valid saccades and saccade start time after fixoff

AmplitudeThreshold=1;
TimeThreshold_ori=FixOffTime-StartTimeBin; 
TimeThreshold=reshape(TimeThreshold_ori,numel(TimeThreshold_ori),1);
TimeThreshold=TimeThreshold(NoNNaNSelector);
TimeThreshold=num2cell(TimeThreshold);
%TimeThresholdEnd=num2cell(MaxTimeLength-PlotSacPostBuffer);
%{
SaccadeAngle_ITI=fillinNaNs(cellfun(@(x,y,z,m) z(((x>AmplitudeThreshold)&(y<=m)))',SaccadeAmplitude,SaccadeEndTime,SaccadeAngle,TimeThresholdEnd,'UniformOutput' ,0));
SaccadeEndPoint_ITI=(cellfun(@(x,y,z,m) z(((x>AmplitudeThreshold)&(y<=m)),:)',SaccadeAmplitude,SaccadeEndTime,SaccadeEndPoint,TimeThresholdEnd,'UniformOutput' ,0));
SaccadeStartPoint_ITI=(cellfun(@(x,y,z,m) z(((x>AmplitudeThreshold)&(y<=m)),:)',SaccadeAmplitude,SaccadeEndTime,SaccadeStartPoint,TimeThresholdEnd,'UniformOutput' ,0));
%}

SaccadeAngle_Explore=fillinNaNs(cellfun(@(x,y,z,m) z(((x>AmplitudeThreshold)&(y>m)))',SaccadeAmplitude,SaccadeStartTime,SaccadeAngle,TimeThreshold,'UniformOutput' ,0));
SaccadeEndPoint_Explore=(cellfun(@(x,y,z,m) z(((x>AmplitudeThreshold)&(y>m)),:)',SaccadeAmplitude,SaccadeStartTime,SaccadeEndPoint,TimeThreshold,'UniformOutput' ,0));
SaccadeStartPoint_Explore=(cellfun(@(x,y,z,m) z(((x>AmplitudeThreshold)&(y>m)),:)',SaccadeAmplitude,SaccadeStartTime,SaccadeStartPoint,TimeThreshold,'UniformOutput' ,0));
%}


SaccadeNumInEachTrial=cellfun(@(x) size(x,2),SaccadeEndPoint_Explore);

SaccadeEndPoint_Explore_X=NaN*ones(length(SaccadeNumInEachTrial),max(SaccadeNumInEachTrial));
SaccadeEndPoint_Explore_Y=NaN*ones(length(SaccadeNumInEachTrial),max(SaccadeNumInEachTrial));

SaccadeStartPoint_Explore_X=NaN*ones(length(SaccadeNumInEachTrial),max(SaccadeNumInEachTrial));
SaccadeStartPoint_Explore_Y=NaN*ones(length(SaccadeNumInEachTrial),max(SaccadeNumInEachTrial));


SaccadeEndPoint_Explore_X(SaccadeNumInEachTrial>0,:)=fillinNaNs(cellfun(@(x) x(1,:),SaccadeEndPoint_Explore(SaccadeNumInEachTrial>0,:),'UniformOutput' ,0));
SaccadeEndPoint_Explore_Y(SaccadeNumInEachTrial>0,:)=fillinNaNs(cellfun(@(x) x(2,:),SaccadeEndPoint_Explore(SaccadeNumInEachTrial>0,:),'UniformOutput' ,0));

SaccadeStartPoint_Explore_X(SaccadeNumInEachTrial>0,:)=fillinNaNs(cellfun(@(x) x(1,:),SaccadeStartPoint_Explore(SaccadeNumInEachTrial>0,:),'UniformOutput' ,0));
SaccadeStartPoint_Explore_Y(SaccadeNumInEachTrial>0,:)=fillinNaNs(cellfun(@(x) x(2,:),SaccadeStartPoint_Explore(SaccadeNumInEachTrial>0,:),'UniformOutput' ,0));

%{
SaccadeStartTime_Explore=fillinNaNs(cellfun(@(x,y,z,m) z(((x>AmplitudeThreshold)&(y<=m)),:)',SaccadeAmplitude,SaccadeStartTime,SaccadeStartTime,TimeThresholdEnd,'UniformOutput' ,0));
SaccadeEndTime_Explore=fillinNaNs(cellfun(@(x,y,z,m) z(((x>AmplitudeThreshold)&(y<=m)),:)',SaccadeAmplitude,SaccadeStartTime,SaccadeEndTime,TimeThresholdEnd,'UniformOutput' ,0));


SaccadeAmplitude_Explore=fillinNaNs(cellfun(@(x,y,z,m) z(((x>AmplitudeThreshold)&(y<=m)),:)',SaccadeAmplitude,SaccadeStartTime,SaccadeAmplitude,TimeThresholdEnd,'UniformOutput' ,0));
%}
%

SaccadeStartTime_Explore=fillinNaNs(cellfun(@(x,y,z,m) z(((x>AmplitudeThreshold)&(y>m)),:)',SaccadeAmplitude,SaccadeStartTime,SaccadeStartTime,TimeThreshold,'UniformOutput' ,0));
SaccadeEndTime_Explore=fillinNaNs(cellfun(@(x,y,z,m) z(((x>AmplitudeThreshold)&(y>m)),:)',SaccadeAmplitude,SaccadeStartTime,SaccadeEndTime,TimeThreshold,'UniformOutput' ,0));


SaccadeAmplitude_Explore=fillinNaNs(cellfun(@(x,y,z,m) z(((x>AmplitudeThreshold)&(y>m)),:)',SaccadeAmplitude,SaccadeStartTime,SaccadeAmplitude,TimeThreshold,'UniformOutput' ,0));

EyeChannel_Ecc_Explore=sqrt(EyeChannel_X_Explore.^2+EyeChannel_Y_Explore.^2);
EyeChannel_V_Explore=abs(diff(EyeChannel_Ecc_Explore,1,2));

SaccadeSep=[3];
ObjOrderSep=[1,2];
% Loop over the blocks
max_Val=0;
min_Val=999;
for i=1:NumberOfBlocks
    
      currBlocks=(BlockType_Line==UniqueBlockType(i));%&(SaccadeNumInEachTrial==1);%&(SearchingTime_Line<500)%& ismember(TargetNum_Line,[1])& SearchingTime_Line>1000;%& ismember(SceneName_Line,SceneName_Unique(1));%
     % RearchingTime_currBlock=SearchingTime_Line( currBlocks);
      %[num,edges]=histcounts(RearchingTime_currBlock,'BinWidth',50);
      %PSTH before searching  
      spon_currBlocks=Mean_Spon_ForEachTrial_Line(currBlocks);
     PSTH_currBlocks=MeanPSTH_FR( currBlocks,:)-spon_currBlocks;
     PSTH_time_currBlocks=MeanPSTH_Time( currBlocks,:);
     
     PSTH_currBlock_mean{i}=nanmean(PSTH_currBlocks,1);
     PSTH_currBlock_sem{i}=nanstd(PSTH_currBlocks,[],1)./sqrt(sum(~isnan(PSTH_currBlocks),1));
     PSTH_time_currBlock_mean{i}=nanmean(PSTH_time_currBlocks,1);
      
    
    %Loop over saccade conditions  
   for j=1:length(SaccadeSep)+1
          
              if j==1
                   sel=currBlocks & SaccadeNumInEachTrial<=SaccadeSep(j);
                   sac_string{j}=sprintf('Upcoming below %d sacs (N=%d)',SaccadeSep(j),sum(sel));
                   
              elseif j<=length(SaccadeSep)
                  
                   sel=currBlocks & SaccadeNumInEachTrial>SaccadeSep(j-1)& SaccadeNumInEachTrial<=SaccadeSep(j);
                   sac_string{j}=sprintf('Upcoming %d to %d sacs (N=%d)',SaccadeSep(j-1),SaccadeSep(j),sum(sel));
              else
                   sel=currBlocks & SaccadeNumInEachTrial>SaccadeSep(j-1);
                     sac_string{j}=sprintf('Upcoming over %d sacs (N=%d)',SaccadeSep(j-1),sum(sel));
              end
           
            
         
     %PSTH before searching  
     PSTH_tmp=MeanPSTH_FR(sel,:)-Mean_Spon_ForEachTrial_Line(sel);
     PSTH_time_tmp=MeanPSTH_Time(sel,:);
     
     PSTH_currBlock_SacSep_mean{i,j}=nanmean(PSTH_tmp,1);
     PSTH_currBlock_SacSep_sem{i,j}=nanstd(PSTH_tmp,[],1)./sqrt(sum(~isnan(PSTH_tmp),1));
     PSTH_time_currBlock_SacSep_mean{i,j}=nanmean(PSTH_time_tmp,1);
     
     Eye_V_Explore_currBlock_SacSep{i,j}=EyeChannel_V_Explore(sel,:);
     
     RearchingTime_SacSep{i,j}=nanmean(SearchingTime_Line( sel));
     
     max_Val_tmp=max(PSTH_currBlock_SacSep_mean{i,j}(1:(200+PlotStartOffset)/StepSize));
     min_Val_tmp=min(PSTH_currBlock_SacSep_mean{i,j}(1:(200+PlotStartOffset)/StepSize));
     if max_Val_tmp>max_Val
         max_Val=max_Val_tmp;
     end
     if min_Val_tmp<min_Val
         min_Val=min_Val_tmp;
     end
      
      
      %PSTH after researching
      
      MeanPSTH_Time_Fractal_tmp=MeanPSTH_Time_Fractal(sel,:)-Mean_Spon_ForEachTrial_Line(sel);
      MeanPSTH_PSTH_Fractal_tmp=MeanPSTH_FR_Fractal(sel,:);
      
      PSTH_currBlock_SacSep_Fractal_mean{i,j}=nanmean(MeanPSTH_PSTH_Fractal_tmp,1);
      PSTH_currBlock_SacSep_Fractal_sem{i,j}=nanstd(MeanPSTH_PSTH_Fractal_tmp,[],1)./sqrt(sum(~isnan(MeanPSTH_PSTH_Fractal_tmp),1));
     PSTH_time_currBlock_SacSep_Fractal_mean{i,j}=nanmean(MeanPSTH_Time_Fractal_tmp,1);
  end    
     %Separate into different groups
     %%
%I. Sort neural responses according to the Rearching time

%GroupStr={'Highest','Higher Medium','Mid Medium','Lowest'};



RearchingTime_currBlock=SearchingTime_Line( currBlocks);
PSTH_currBlock=MeanPSTH_FR(currBlocks,:);%-spon_currBlocks;
PSTH_Time_currBlock=MeanPSTH_Time(currBlocks,:);

%RT_First_Target=SacEnd_First_Target;
SepGroupAmp=4;
SepGroupAmpNum= SepGroupAmp;
%SepGroupAmp={[2300,1800,1100]};
%{
%SepGroupAmp={[500]};
GroupStr={'Highest','Higher Medium','Mid Medium','Lowest'};
%GroupStr={'Long','Medium','Short'};
SepGroupAmpNum=length(SepGroupAmp{1})+1;
%}
GroupStr={'Highest','Higher Medium','Mid Medium','Lowest'};
%Before Exploration
PSTH_group_temp=SortDataAccordingToGroup(PSTH_currBlock,RearchingTime_currBlock,SepGroupAmp,'descend');


PSTH_RTGroup_mean{i}=PSTH_group_temp.mean(1:SepGroupAmpNum,:);
PSTH_RTGroup_sem{i}=PSTH_group_temp.sem(1:SepGroupAmpNum,:);

RT_Group_mean{i}=PSTH_group_temp.MarkerMean;
RT_Group_sem{i}=PSTH_group_temp.MarkerSem;

NumInGroupRT=PSTH_group_temp.Num;
RT_SortIndex=PSTH_group_temp.SortIndex;

PSTH_RTGroup_Time{i}=nanmean(PSTH_Time_currBlock,1);


%After Fractal Onset   
PSTH_Fractal_currBlock=MeanPSTH_FR_Fractal(currBlocks,:);%-spon_currBlocks;
PSTH_Time__Fractal_currBlock=MeanPSTH_Time_Fractal(currBlocks,:);

PSTH_group_temp=SortDataAccordingToGroup(PSTH_Fractal_currBlock,RearchingTime_currBlock,SepGroupAmp,'descend');


PSTH_RTGroup_Fractal_mean{i}=PSTH_group_temp.mean(1:SepGroupAmpNum,:);
PSTH_RTGroup_Fractal_sem{i}=PSTH_group_temp.sem(1:SepGroupAmpNum,:);


PSTH_RTGroup_Fractal_Time{i}=nanmean(PSTH_Time__Fractal_currBlock,1);

 %Separate into different fractal orders
  for j=1:length(ObjOrderSep)+1
          if j==length(ObjOrderSep)+1
             sel=currBlocks & TargetOrder_Line>ObjOrderSep(end);
             order_string{j}=sprintf('Fractal Order> %d (N=%d)',ObjOrderSep(end),sum(sel));
          else
            sel=currBlocks & TargetOrder_Line==ObjOrderSep(j);
            order_string{j}=sprintf('Fractal Order %d (N=%d)',ObjOrderSep(j),sum(sel));
          end
     %PSTH before searching 
     PSTH_tmp=MeanPSTH_FR(sel,:);%-Mean_Spon_ForEachTrial_Line(sel);
     PSTH_time_tmp=MeanPSTH_Time(sel,:);
     
     PSTH_currBlock_OrderSep_mean{i,j}=nanmean(PSTH_tmp,1);
     PSTH_currBlock_OrderSep_sem{i,j}=nanstd(PSTH_tmp,[],1)./sqrt(sum(~isnan(PSTH_tmp),1));
     PSTH_time_currBlock_OrderSep_mean{i,j}=nanmean(PSTH_time_tmp,1);
     
     Eye_V_Explore_currBlock_OrderSep{i,j}=EyeChannel_V_Explore(sel,:);
     
     RearchingTime_OrderSep{i,j}=nanmean(SearchingTime_Line( sel));
     
     %{
     max_Val_tmp=max(PSTH_currBlock_SacSep_mean{i,j}(1:(200+PlotStartOffset)/StepSize));
     min_Val_tmp=min(PSTH_currBlock_SacSep_mean{i,j}(1:(200+PlotStartOffset)/StepSize));
     if max_Val_tmp>max_Val
         max_Val=max_Val_tmp;
     end
     if min_Val_tmp<min_Val
         min_Val=min_Val_tmp;
     end 
 %}
      %PSTH after researching
      
      MeanPSTH_Time_Fractal_tmp=MeanPSTH_Time_Fractal(sel,:);%-Mean_Spon_ForEachTrial_Line(sel);
      MeanPSTH_PSTH_Fractal_tmp=MeanPSTH_FR_Fractal(sel,:);
      
      PSTH_currBlock_OrderSep_Fractal_mean{i,j}=nanmean(MeanPSTH_PSTH_Fractal_tmp,1);
      PSTH_currBlock_OrderSep_Fractal_sem{i,j}=nanstd(MeanPSTH_PSTH_Fractal_tmp,[],1)./sqrt(sum(~isnan(MeanPSTH_PSTH_Fractal_tmp),1));
     PSTH_time_currBlock_OrderSep_Fractal_mean{i,j}=nanmean(MeanPSTH_Time_Fractal_tmp,1);
     
     
     
     
     
     
  end
      
     
end

%PlotStartOffset
if ShowFigureFlag
    FigureStartNum=140;
    FigureIndex=1;
    LineWidth=2;
    Barlength=1;
    %Figure1, Response before fixation offset for each mode
    figtitlestr{FigureIndex}='PSTH_FixOff_Each_Blocks';
    fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[30,100, 800,800],'Name',figtitlestr{FigureIndex});
 
        %subplot(ceil((length(SaccadeSep)+1)/2),2,j)
       
        for i=1:NumberOfBlocks
           % Eye_V_Explore_currBlock_SacSep{i,j}
            currBlock=UniqueBlockType(i);
            BlockOrder(i)=find(AssociateBlock==currBlock);
             shadedErrorBar(PSTH_time_currBlock_mean{i},PSTH_currBlock_mean{i},PSTH_currBlock_sem{i},'lineprops',{ColorForEachMode(AssociateBlock==currBlock,:)},'transparent',1,'patchSaturation',0.5);%PSTH for fix
             set(findall(gca, 'Type', 'Line'),'LineWidth',3);
             hold on
             % xlim([-PlotStartOffset(1),RearchingTime_SacSep{1,1}]);
               xlim([-PlotStartOffset(1),200]);
              
              
              %v_tmp=Eye_V_Explore_currBlock_SacSep{i,j};
              %x_tmp=repmat((-PlotStartOffset(1):-PlotStartOffset(1)+size(v_tmp,2)-1),size(v_tmp,1),1);
            % plot(x_tmp',yrange(2)*1.2+v_tmp','-k','color','#9D9EA1');
            %sh=plot(x_tmp',max_Val*1.2+v_tmp','-k','color',ColorForEachMode(AssociateBlock==currBlock,:));
            %Some legend off
            %for ii=1:length(h)
            %set( get( get( h(ii), 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );
            %end
 
        end
        
         legend(Blockstr{BlockOrder});
        
        
         %ylim([min_Val*0.8,max_Val*1.2+2]);   
            
       % xlim([-PlotStartOffset,RearchingTime_SacSep{i,j}]);
     %  xlim([-PlotStartOffset,RearchingTime_SacSep{1,1}]);
        xlabel('Time from FixOff(ms)');
        ylabel('Firing Rate(Hz) over spon');
       % title(sac_string{j});
         set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
               
            
    
   
     FigureStartNum=FigureStartNum+1;
     FigureIndex=FigureIndex+1;
    
    %Figure2, Response before fixation offset for different subsequent
    %number of saccades for each mode
    figtitlestr{FigureIndex}='PSTH_DifferentSacNum_Blocks';
    fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[50,100, 800,800],'Name',figtitlestr{FigureIndex});
    
     for j=1:(length(SaccadeSep)+1)
         
   
        
        subplot(ceil((length(SaccadeSep)+1)/2),2,j)
       
        for i=1:NumberOfBlocks
           % Eye_V_Explore_currBlock_SacSep{i,j}
            currBlock=UniqueBlockType(i);
            BlockOrder(i)=find(AssociateBlock==currBlock);
             shadedErrorBar(PSTH_time_currBlock_SacSep_mean{i,j},PSTH_currBlock_SacSep_mean{i,j},PSTH_currBlock_SacSep_sem{i,j},'lineprops',{ColorForEachMode(AssociateBlock==currBlock,:)},'transparent',1,'patchSaturation',0.5);%PSTH for fix
             set(findall(gca, 'Type', 'Line'),'LineWidth',3);
             hold on
            %  xlim([-PlotStartOffset(1),RearchingTime_SacSep{1,1}]);
             xlim([-PlotStartOffset(1),200]);
              
              v_tmp=Eye_V_Explore_currBlock_SacSep{i,j};
              x_tmp=repmat((-PlotStartOffset(1):-PlotStartOffset(1)+size(v_tmp,2)-1),size(v_tmp,1),1);
            % plot(x_tmp',yrange(2)*1.2+v_tmp','-k','color','#9D9EA1');
            h=plot(x_tmp',max_Val*1.2+v_tmp','-k','color',ColorForEachMode(AssociateBlock==currBlock,:));
            %Some legend off
            for ii=1:length(h)
            set( get( get( h(ii), 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );
            end
 
        end
        if j==1
         legend(Blockstr{BlockOrder});
        end
        
         ylim([min_Val*0.8,max_Val*1.2+2]);   
            
       % xlim([-PlotStartOffset,RearchingTime_SacSep{i,j}]);
     %  xlim([-PlotStartOffset,RearchingTime_SacSep{1,1}]);
        xlabel('Time from FixOff(ms)');
        ylabel('Firing Rate(Hz)');
        title(sac_string{j});
         set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
               
            
    end     
        
   
     FigureStartNum=FigureStartNum+1;
     FigureIndex=FigureIndex+1;
  %Figure2, Response before fixation offset for different subsequent
    %Searching groups for each mode
    figtitlestr{FigureIndex}='PSTH_DifferentRT_Blocks';
    fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[70,100, 800,800],'Name',figtitlestr{FigureIndex});
    
     for j=1:(SepGroupAmpNum)
         
   
        
        subplot(ceil((SepGroupAmpNum)/2),2,j)
       
        for i=1:NumberOfBlocks
           % Eye_V_Explore_currBlock_SacSep{i,j}
            currBlock=UniqueBlockType(i);
            BlockOrder(i)=find(AssociateBlock==currBlock);
             shadedErrorBar(PSTH_RTGroup_Time{i},PSTH_RTGroup_mean{i}(j,:),PSTH_RTGroup_sem{i}(j,:),'lineprops',{ColorForEachMode(AssociateBlock==currBlock,:)},'transparent',1,'patchSaturation',0.5);%PSTH for fix
             set(findall(gca, 'Type', 'Line'),'LineWidth',3);
             hold on
            %  xlim([-PlotStartOffset(1),RearchingTime_SacSep{1,1}]);
            xlim([-PlotStartOffset(1),200]);
             %{ 
              v_tmp=Eye_V_Explore_currBlock_SacSep{i,j};
              x_tmp=repmat((-PlotStartOffset:-PlotStartOffset+size(v_tmp,2)-1),size(v_tmp,1),1);
            % plot(x_tmp',yrange(2)*1.2+v_tmp','-k','color','#9D9EA1');
            h=plot(x_tmp',max_Val*1.2+v_tmp','-k','color',ColorForEachMode(AssociateBlock==currBlock,:));
            %Some legend off
            for ii=1:length(h)
            set( get( get( h(ii), 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );
            end
              %}
 
        end
        if j==1
         legend(Blockstr{BlockOrder});
        end
        
       %  ylim([min_Val*0.8,max_Val*1.2+2]);   
            
       % xlim([-PlotStartOffset,RearchingTime_SacSep{i,j}]);
     %  xlim([-PlotStartOffset,RearchingTime_SacSep{1,1}]);
        xlabel('Time from FixOff(ms)');
        ylabel('Firing Rate(Hz)');
        title(GroupStr{j});
         set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
               
            
    end     
        
   
     FigureStartNum=FigureStartNum+1;
     FigureIndex=FigureIndex+1;
     
     %Figure3, Response before fixation offset for different subsequent
    %object orders for each block
    figtitlestr{FigureIndex}='PSTH_DifferentFractalOrder_Blocks';
    fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[90,100, 800,800],'Name',figtitlestr{FigureIndex});
    
     for j=1:(length(ObjOrderSep)+1)
         
   
        
        subplot(ceil((length(ObjOrderSep)+1)/2),2,j)
       
        for i=1:NumberOfBlocks
           % Eye_V_Explore_currBlock_SacSep{i,j}
            currBlock=UniqueBlockType(i);
            BlockOrder(i)=find(AssociateBlock==currBlock);
             shadedErrorBar(PSTH_time_currBlock_OrderSep_mean{i,j},PSTH_currBlock_OrderSep_mean{i,j},PSTH_currBlock_OrderSep_sem{i,j},'lineprops',{ColorForEachMode(AssociateBlock==currBlock,:)},'transparent',1,'patchSaturation',0.5);%PSTH for fix
             set(findall(gca, 'Type', 'Line'),'LineWidth',3);
             hold on
             % xlim([-PlotStartOffset(1),RearchingTime_OrderSep{1,1}]);
              xlim([-PlotStartOffset(1),200]);
              v_tmp=Eye_V_Explore_currBlock_OrderSep{i,j};
              x_tmp=repmat((-PlotStartOffset(1):-PlotStartOffset(1)+size(v_tmp,2)-1),size(v_tmp,1),1);
            % plot(x_tmp',yrange(2)*1.2+v_tmp','-k','color','#9D9EA1');
            h=plot(x_tmp',max_Val*1.2+v_tmp','-k','color',ColorForEachMode(AssociateBlock==currBlock,:));
            %Some legend off
            for ii=1:length(h)
            set( get( get( h(ii), 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );
            end
 
        end
        if j==1
         legend(Blockstr{BlockOrder});
        end
        
       %  ylim([min_Val*0.8,max_Val*1.2+2]);   
            
       % xlim([-PlotStartOffset,RearchingTime_SacSep{i,j}]);
     %  xlim([-PlotStartOffset,RearchingTime_SacSep{1,1}]);
        xlabel('Time from FixOff(ms)');
        ylabel('Firing Rate(Hz)');
        title(order_string{j});
         set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
               
            
    end     
        
   
     FigureStartNum=FigureStartNum+1;
     FigureIndex=FigureIndex+1;
     
      %Figure4, Response after Fractal Onset for different subsequent
    %number of saccades for each mode
    figtitlestr{FigureIndex}='FractalPSTH_DifferentSacNum_Blocks';
    fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[110,100, 800,800],'Name',figtitlestr{FigureIndex});
    
     for j=1:(length(SaccadeSep)+1)
         
   
        
        subplot(ceil((length(SaccadeSep)+1)/2),2,j)
       
        for i=1:NumberOfBlocks
           % Eye_V_Explore_currBlock_SacSep{i,j}
            currBlock=UniqueBlockType(i);
            BlockOrder(i)=find(AssociateBlock==currBlock);
            
             shadedErrorBar(PSTH_time_currBlock_SacSep_Fractal_mean{i,j},PSTH_currBlock_SacSep_Fractal_mean{i,j},PSTH_currBlock_SacSep_Fractal_sem{i,j},'lineprops',{ColorForEachMode(AssociateBlock==currBlock,:)},'transparent',1,'patchSaturation',0.5);%PSTH for fix
             set(findall(gca, 'Type', 'Line'),'LineWidth',3);
           %  hold on
            % xlim([-PlotStartOffset(2),RearchingTime_SacSep{1,1}]);
              xlim([-PlotStartOffset(2),500]);
             %{ 
              v_tmp=Eye_V_Explore_currBlock_SacSep{i,j};
              x_tmp=repmat((-PlotStartOffset(2):-PlotStartOffset(2)+size(v_tmp,2)-1),size(v_tmp,1),1);
            % plot(x_tmp',yrange(2)*1.2+v_tmp','-k','color','#9D9EA1');
            h=plot(x_tmp',max_Val*1.2+v_tmp','-k','color',ColorForEachMode(AssociateBlock==currBlock,:));
            %Some legend off
            for ii=1:length(h)
            set( get( get( h(ii), 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );
            end
              %}
 
        end
        if j==1
         legend(Blockstr{BlockOrder});
        end
        
       % ylim([min_Val*0.8,max_Val*1.2+2]);   
            
       % xlim([-PlotStartOffset,RearchingTime_SacSep{i,j}]);
     %  xlim([-PlotStartOffset,RearchingTime_SacSep{1,1}]);
        xlabel('Time from ObjectOn(ms)');
        ylabel('Firing Rate(Hz)');
        title(sac_string{j});
         set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
               
            
    end     
        
   
     FigureStartNum=FigureStartNum+1;
     FigureIndex=FigureIndex+1;
     
     
   %Figure5 Response for fractals for different s
    %Searching groups for each mode
    figtitlestr{FigureIndex}='PSTH_Fractal_DifferentRT_Blocks';
    fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[130,100, 800,800],'Name',figtitlestr{FigureIndex});
    
     for j=1:(SepGroupAmpNum)
         
   
        
        subplot(ceil((SepGroupAmpNum)/2),2,j)
       
        for i=1:NumberOfBlocks
           % Eye_V_Explore_currBlock_SacSep{i,j}
            currBlock=UniqueBlockType(i);
            BlockOrder(i)=find(AssociateBlock==currBlock);
             shadedErrorBar(PSTH_RTGroup_Fractal_Time{i},PSTH_RTGroup_Fractal_mean{i}(j,:),PSTH_RTGroup_Fractal_sem{i}(j,:),'lineprops',{ColorForEachMode(AssociateBlock==currBlock,:)},'transparent',1,'patchSaturation',0.5);%PSTH for fix
             set(findall(gca, 'Type', 'Line'),'LineWidth',3);
             hold on
              xlim([-PlotStartOffset(2),500]);
             %{ 
              v_tmp=Eye_V_Explore_currBlock_SacSep{i,j};
              x_tmp=repmat((-PlotStartOffset:-PlotStartOffset+size(v_tmp,2)-1),size(v_tmp,1),1);
            % plot(x_tmp',yrange(2)*1.2+v_tmp','-k','color','#9D9EA1');
            h=plot(x_tmp',max_Val*1.2+v_tmp','-k','color',ColorForEachMode(AssociateBlock==currBlock,:));
            %Some legend off
            for ii=1:length(h)
            set( get( get( h(ii), 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );
            end
              %}
 
        end
        if j==1
         legend(Blockstr{BlockOrder});
        end
        
      %   ylim([min_Val*0.8,max_Val*1.2+2]);   
          %  ylim([0,150]); 
       % xlim([-PlotStartOffset,RearchingTime_SacSep{i,j}]);
     %  xlim([-PlotStartOffset,RearchingTime_SacSep{1,1}]);
        xlabel('Time from Fractal On(ms)');
        ylabel('Firing Rate(Hz)');
        title(GroupStr{j});
         set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
               
            
    end     
        
   
     FigureStartNum=FigureStartNum+1;
     FigureIndex=FigureIndex+1;
        
   %Figure6, Response after Fractal Onset for different subsequent
    %object orders for each block
    figtitlestr{FigureIndex}='PSTH_Fractal_DifferentFractalOrder_Blocks';
    fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[150,100, 800,800],'Name',figtitlestr{FigureIndex});
    
     for j=1:(length(ObjOrderSep)+1)
         
   
        
        subplot(ceil((length(ObjOrderSep)+1)/2),2,j)
       
        for i=1:NumberOfBlocks
           % Eye_V_Explore_currBlock_SacSep{i,j}
            currBlock=UniqueBlockType(i);
            BlockOrder(i)=find(AssociateBlock==currBlock);
             shadedErrorBar(PSTH_time_currBlock_OrderSep_Fractal_mean{i,j},PSTH_currBlock_OrderSep_Fractal_mean{i,j},PSTH_currBlock_OrderSep_Fractal_sem{i,j},'lineprops',{ColorForEachMode(AssociateBlock==currBlock,:)},'transparent',1,'patchSaturation',0.5);%PSTH for fix
             set(findall(gca, 'Type', 'Line'),'LineWidth',3);
             hold on
              xlim([-PlotStartOffset(2),500]);
              %{
              v_tmp=Eye_V_Explore_currBlock_OrderSep{i,j};
              x_tmp=repmat((-PlotStartOffset(1):-PlotStartOffset(1)+size(v_tmp,2)-1),size(v_tmp,1),1);
            % plot(x_tmp',yrange(2)*1.2+v_tmp','-k','color','#9D9EA1');
            h=plot(x_tmp',max_Val*1.2+v_tmp','-k','color',ColorForEachMode(AssociateBlock==currBlock,:));
            %Some legend off
            for ii=1:length(h)
            set( get( get( h(ii), 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );
            end
              %}
 
        end
        if j==1
         legend(Blockstr{BlockOrder});
        end
        
       %  ylim([min_Val*0.8,max_Val*1.2+2]);   
         %  ylim([0,150]); 
       % xlim([-PlotStartOffset,RearchingTime_SacSep{i,j}]);
     %  xlim([-PlotStartOffset,RearchingTime_SacSep{1,1}]);
        xlabel('Time from FractalOn(ms)');
        ylabel('Firing Rate(Hz)');
        title(order_string{j});
         set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
               
            
    end     
        
   
     FigureStartNum=FigureStartNum+1;
     FigureIndex=FigureIndex+1;
        
     
     
end   
    
%}

end
function organized=ReorganizeData(data);
organized=[];
for i=1:length(data)
    data_curr=data{i};
    maxnumel=max(cellfun(@numel,data_curr));
    organized{i}=cell2mat(cellfun(@(x) [x,NaN*ones(1,maxnumel-length(x))],data_curr,'uniform',0));
    
end

 maxnumel=max(cellfun(@(x) size(x,2),organized));
  organized=cell2mat(cellfun(@(x) [x,NaN*ones(size(x,1),maxnumel-size(x,2))],organized,'uniform',0)');


end