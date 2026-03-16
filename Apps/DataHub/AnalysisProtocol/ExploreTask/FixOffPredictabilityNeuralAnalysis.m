function FixOffPredictabilityNeuralAnalysis(Data, StartTrial, EndTrial,ShowFigureFlag,OutputFlag);
%The location predictability for responses before fix off
%Protocol to analysis the achieving fractal effect and fractal order
%Xuefei Yu 09142020
ECode;
TaskType=Data.TaskType;
DataPath= Data.Path;
FileName=Data.FileName;

Selection=StartTrial:EndTrial;%Select the trials according to the pannel
%Abnormal check %Especially For Robin081120FORAGEXPLO4.rst
CheckCode=SCENE_ON;%This code should appear once in each trial
AbnormalTrials=AbnormalTrialCheck(Data,CheckCode);
Selection=Selection(~ismember(Selection,AbnormalTrials));

if strcmp(FileName,'Robin092920FORAGEXPLO1.rst')
    %Remove bad trial for this file

CheckCodePair=[FIXOFF,SHOWFIXCD];%This code should be the same as target on for good trials
AbnormalTrials=AbnormalTrialCheck(Data,CheckCodePair);
Selection=Selection(~ismember(Selection,AbnormalTrials));
end





TrialNum=length(Selection);
%User defined parameter for the neural/behavioral analysis
%TargetOnOffset=100;%in ms, the time selection bin after target onset; 
%StartDelay=-50;%in ms, the time selection bin after fixation offset;

%For Plot of PSTH aligned on Each fixation on and each fixation Off
PlotFixationOn_StartOffset=100;
PlotFixationOn_StopOffset=300;%Duration of plotting

PlotFixationOff_StartOffset=300;
PlotFixationOff_StopOffset=500;%Duration of plotting


%For segmentation of  the neural responses
%{
SegEvent=[FIXOFF,CUEON];%Event Flow of interest
SegEventEnd=[NaN,TGTCD];%Event end
SegEventStr={'BeforeExplore','CueOn'};
%}
%
SegEvent=[FIXOFF,TGTCD];%Event Flow of interest
SegEventEnd=[NaN, TGTOFF];%Event end
SegEventStr={'BeforeExplore','CueOn'};
%}
PlotStartOffset=[1200,200];%in ms
%PlotStartOffset=[1700];%in ms
PlotStopOffset=[200,200];%in ms

%Select Average Interval
FixOffInterval=[-500,0];
FracOnInterval=[0,600];


EventTypeFlag=[2,2];%1:Number same with the trial number, 2:Number same with the target number

%Choose Num of Location response to average;
ChooseAverageLoc=3;%Average highest 3 and lowest 3;

%SceneRemove=[1054:1057,1060:1063];
SceneRemove=[];

%Spon Event
SponEvent=[SCENE_OFF];
SponPreTimeBuffer=-500;
SponEventEnd=NaN;
SponStartOffset=0;
SponStopOffset=1000;

%Parameter to define the sliding window to measure firing rates
StepSize=10;%in ms
BinWidth=20;%in ms

%Plot Raster Flag
PlotRasterFlag=0;


%Some user defined parameter for visualization
%Color Definition for each saccade nums
%Color={'#FB5441','#4173fb','#FBBB41','#0EC335',[0.28,0.78,0.89],'k',[0.6,0.1,1],[0.6,0.1,0.4],[0.3,0.5,1],[0.6,0.5,0.1],[0,0.8,0.4],'r','g','b','m'};
Color={'r',[0.3,0.8,0.0],[0.3,0.5,1],[0.7,0.3,1],[0.4,0.4,0],[0.28,0.78,0.89],'k',[0.5,0.2,1],[0.6,0.1,0.4],[0.3,0.5,1],[0.6,0.5,0.1],[0,0.8,0.4]};


ColorForEachMode=[0.83 0.14 0.14
             1.00 0.54 0.00
             0.47 0.25 0.80
             0.25 0.80 0.54];
AssociateBlock=[10,11,0,1];
Blockstr={'FixLoc,FixSceneSeq','FixLoc,RandSceneSeq','RandLoc,FixSceneSeq','RandLoc,RandSceneSeq'};


%BlockSepTrial=unique(SceneLib(:,1));
BlockType=Data.BlockType(:,:);

BlockType=BlockType(Selection,2);

%Get File Information
num = regexp(Data.FileName, '\d+', 'match');
RecordDateOriginal=cell2mat(num(1));
RecordDate=datetime(RecordDateOriginal,'InputFormat','MMddyy');

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
BlockType(BlockType==10)=NaN;
BlockType(BlockType==0)=NaN;

%BlockType(BlockType(:,3)>0,3)=999;



%GetEyeData
EyeData=Data.EyeDataRaw;
EyeChannel_X_L=EyeData(Selection,1);
EyeChannel_Y_L=EyeData(Selection,2);



EyeBinWidth=Data.EyeBinWidth;
EyeCalInfo=Data.EyeCalInfo;%Each row: X_Gain,X_Bias; Y_Gain; Y_bias; X=(raw+X_bias)*X_Gain
EyeCalInfo_L=EyeCalInfo(1,:);

%Load event data
EventChannel=Data.EventChannel(Selection,:);
EventTimeChannel=Data.EventTimeChannel(Selection,:);
EventBin=Data.EventBin(Selection,:);
%Get the information about the extra FP
%
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
 
 SceneOnTime=FindOutTime(EventChannel,EventTimeChannel,SCENE_ON);
 %screen out the searching trials
 FixOffTime=FixOffTime(AdvanceInTargetTrial,:);
 OutScreenTime=OutScreenTime(AdvanceInTargetTrial,:);
 FPOnTime=FPOnTime(AdvanceInTargetTrial,:);
 
 FixOffTime=ScreenOutOverScreen(FixOffTime,OutScreenTime);
 
 SceneOnTime=SceneOnTime(AdvanceInTargetTrial,:)-FixOffTime(:,1);%Align with the fix off
 
 
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
 ValidBlockType=~isnan(UniqueBlockType);
 UniqueBlockType=UniqueBlockType(ValidBlockType);
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
MeanPSTH_Time_FixOff=MeanPSTH{1,1};
MeanPSTH_FR_FixOff=MeanPSTH{1,2};

PSTH_Time_FixOff=nanmean(MeanPSTH_Time_FixOff,1);

MeanPSTH_Time_FracOn=MeanPSTH{2,1};
MeanPSTH_FR_FracOn=MeanPSTH{2,2};

PSTH_Time_FracOn=nanmean(MeanPSTH_Time_FracOn,1);

%Raster and raw spike count PSTH
AlignedRaster_FixOff=SegmentData.RasterAfterAlign{1};
AlignedRaster_FracOn=SegmentData.RasterAfterAlign{2};
%PSTH_Event_SpikeCount=SegmentData.PSTH_Event_SpikeCount{1};

%Relative time bin for each event
RelativeTimeMarker_FixOff=SegmentData.RelativeTimeMarker{1};%=EventTimeMarker;
RelativeTimeMarker_FracOn=SegmentData.RelativeTimeMarker{2};%=EventTimeMarker;

AbsoluteTimeMarker_FixOff=SegmentData.AbsoluteTimeMarker{1};%=EventTimeMarker;
AbsoluteTimeMarker_FixOff=AbsoluteTimeMarker_FixOff(:,1);
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

GoodObj_Line=mod(TargetIndex_Line,2)==1;
BadObj_Line=mod(TargetIndex_Line,2)==0;

TargetPolarAngle_Line=TargetPolarAngle;
TargetPolarAngle_Line=reshape(TargetPolarAngle,numel(TargetPolarAngle),1);
TargetPolarAngle_Line=TargetPolarAngle_Line(NoNNaNSelector);


%TargetOrderLine
TargetOrder=1:max(TargetNum);
TotalNum=length(TargetNum);
TargetOrder_Mat=repmat(TargetOrder,TotalNum,1);
TargetOrder_Line=reshape(TargetOrder_Mat,numel(TargetOrder_Mat),1);
TargetOrder_Line=TargetOrder_Line(NoNNaNSelector);

%SceneLine
SceneName_Mat=repmat(SceneName,1,max(TargetNum));
SceneName_Line=reshape(SceneName_Mat,numel(SceneName_Mat),1);
SceneName_Line=SceneName_Line(NoNNaNSelector);

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
%Rearrange_scene on time into lines

SceneOnTime_Line=repmat(SceneOnTime,1,max(TargetNum));
SceneOnTime_Line=reshape(SceneOnTime_Line,numel(SceneOnTime_Line),1);
SceneOnTime_Line=SceneOnTime_Line(NoNNaNSelector);
SceneOnTime_Line=SceneOnTime_Line-RelativeTimeMarker_FixOff;

%FR for spon
EventTimeChannelForSpon=EventTimeChannel-SponPreTimeBuffer;

SegmentDataSpon=EventSegTool(SponEvent,SponEventEnd,SponStartOffset,SponStopOffset,EventChannel,EventTimeChannelForSpon,SpikeTime,StepSize,BinWidth);

Mean_Spon=SegmentDataSpon.MeanFR;
Mean_Spon_ForEachTrial=SegmentDataSpon.MeanFRRaw{1};

Mean_Spon_ForEachTrial_Line=repmat(Mean_Spon_ForEachTrial,1,max(TargetNum));
Mean_Spon_ForEachTrial_Line=reshape(Mean_Spon_ForEachTrial_Line,numel(Mean_Spon_ForEachTrial_Line),1);
Mean_Spon_ForEachTrial_Line=Mean_Spon_ForEachTrial_Line(NoNNaNSelector);


%%
%ObjOrderSep=[1,2];

for i=1:NumberOfBlocks
    if ~isempty(SceneRemove)
        currBlocks=(BlockType_Line==UniqueBlockType(i))& (~ismember(SceneName_Line,SceneRemove));%&ismember(TargetOrder_Line,[1]);%;%&(SaccadeNumInEachTrial==1);%&(SearchingTime_Line<500)%& ismember(TargetNum_Line,[1])& SearchingTime_Line>1000;%& ismember(SceneName_Line,SceneName_Unique(1));%
      
    else
        currBlocks=(BlockType_Line==UniqueBlockType(i));%&(SaccadeNumInEachTrial==1);%&(SearchingTime_Line<500)%& ismember(TargetNum_Line,[1])& SearchingTime_Line>1000;%& ismember(SceneName_Line,SceneName_Unique(1));%
      
    end
    
      
    
%%Separate into different future locations
    
Spontmp=Mean_Spon_ForEachTrial_Line(currBlocks);
RearchingTime_currBlock=SearchingTime_Line( currBlocks);


MeanPSTH_FR_FixOff_currBlock=MeanPSTH_FR_FixOff(currBlocks,:);%-Spontmp;
MeanPSTH_FR_FracOn_currBlock=MeanPSTH_FR_FracOn(currBlocks,:);%-Spontmp;

%Measure the location tuning across conditions
  TargetPolarAngle_curr=TargetPolarAngle_Line(currBlocks);
  TargetPolarAngle_curr_unique{i}=unique(TargetPolarAngle_curr);
  
  TargetOrder_Line_curr=TargetOrder_Line(currBlocks);
  
  SceneOnTime_Line_curr=SceneOnTime_Line(currBlocks);
  
  for j=1:length(TargetPolarAngle_curr_unique{i})
      %Select akk
      ang=TargetPolarAngle_curr_unique{i}(j);
      sel_ang=TargetPolarAngle_curr==ang;
      
      %Only look at the First target: with both scene and before fix off responses
      %TargetPolarAngle
            sel_ang_first=TargetOrder_Line_curr==1 & sel_ang;
            SceneOnTime_Line_first(j)=nanmean(SceneOnTime_Line_curr(sel_ang_first));
      
      if sum(sel_ang)>3 %At lease three
           MeanPSTH_currBlock_FixOffLocation{i,j}=nanmean(MeanPSTH_FR_FixOff_currBlock(sel_ang,:),1);
           SemPSTH_currBlock_FixOffLocation{i,j}=nanstd(MeanPSTH_FR_FixOff_currBlock(sel_ang,:),[],1)./sqrt(sum(~isnan(MeanPSTH_FR_FixOff_currBlock(sel_ang,:)),1));
          
           MeanPSTH_currBlock_FracOnLocation{i,j}=nanmean(MeanPSTH_FR_FracOn_currBlock(sel_ang,:),1);
           SemPSTH_currBlock_FracOnLocation{i,j}=nanstd(MeanPSTH_FR_FracOn_currBlock(sel_ang,:),[],1)./sqrt(sum(~isnan(MeanPSTH_FR_FracOn_currBlock(sel_ang,:)),1));
           
           %Select Interval to average
           FixOffSel=(PSTH_Time_FixOff>=FixOffInterval(1)) &(PSTH_Time_FixOff<=FixOffInterval(2)) ;
          
           FixOff_Location_Average_Mean(i,j)=nanmean(MeanPSTH_FR_FixOff_currBlock(sel_ang,FixOffSel),'all');
           FixOff_Location_Average_Sem(i,j)=nanstd(MeanPSTH_FR_FixOff_currBlock(sel_ang,FixOffSel),0,'all')/sqrt(sum(FixOffSel));
           
           FracOnSel=(PSTH_Time_FracOn>=FracOnInterval(1)) &(PSTH_Time_FracOn<=FracOnInterval(2)) ;
          
           FracOn_Location_Average_Mean(i,j)=nanmean(MeanPSTH_FR_FracOn_currBlock(sel_ang,FracOnSel),'all');
           FracOn_Location_Average_Sem(i,j)=nanstd(MeanPSTH_FR_FracOn_currBlock(sel_ang,FracOnSel),0,'all')/sqrt(sum(FracOnSel));
         
    
      end
       if sum(sel_ang_first)>3 %At lease three
            MeanPSTH_currBlock_FixOffLocation_first{i,j}=nanmean(MeanPSTH_FR_FixOff_currBlock(sel_ang_first,:),1);
           SemPSTH_currBlock_FixOffLocation_first{i,j}=nanstd(MeanPSTH_FR_FixOff_currBlock(sel_ang_first,:),[],1)./sqrt(sum(~isnan(MeanPSTH_FR_FixOff_currBlock(sel_ang_first,:)),1));
          
           MeanPSTH_currBlock_FracOnLocation_first{i,j}=nanmean(MeanPSTH_FR_FracOn_currBlock(sel_ang_first,:),1);
           SemPSTH_currBlock_FracOnLocation_first{i,j}=nanstd(MeanPSTH_FR_FracOn_currBlock(sel_ang_first,:),[],1)./sqrt(sum(~isnan(MeanPSTH_FR_FracOn_currBlock(sel_ang_first,:)),1));
           
           %Select Interval to average
           FixOffSel=(PSTH_Time_FixOff>=FixOffInterval(1)) &(PSTH_Time_FixOff<=FixOffInterval(2)) ;
          
           FixOff_Location_Average_First_Mean(i,j)=nanmean(MeanPSTH_FR_FixOff_currBlock(sel_ang_first,FixOffSel),'all');
           FixOff_Location_Average_First_Sem(i,j)=nanstd(MeanPSTH_FR_FixOff_currBlock(sel_ang_first,FixOffSel),0,'all')/sqrt(sum(FixOffSel));
           
           FracOnSel=(PSTH_Time_FracOn>=FracOnInterval(1)) &(PSTH_Time_FracOn<=FracOnInterval(2)) ;
          
           FracOn_Location_Average_First_Mean(i,j)=nanmean(MeanPSTH_FR_FracOn_currBlock(sel_ang_first,FracOnSel),'all');
           FracOn_Location_Average_First_Sem(i,j)=nanstd(MeanPSTH_FR_FracOn_currBlock(sel_ang_first,FracOnSel),0,'all')/sqrt(sum(FracOnSel));
         
           
           
       end

      
  end
  
  %Direct correlation between response before fixoff and frac on
    FixOffSel=(PSTH_Time_FixOff>=FixOffInterval(1)) &(PSTH_Time_FixOff<=FixOffInterval(2)) ;
    FracOnSel=(PSTH_Time_FracOn>=FracOnInterval(1)) &(PSTH_Time_FracOn<=FracOnInterval(2)) ;
         
  
 FixOff_Location_Average_All{i}=nanmean(MeanPSTH_FR_FixOff_currBlock(:,FixOffSel),2);
          
 FracOn_Location_Average_All{i}=nanmean(MeanPSTH_FR_FracOn_currBlock(:,FracOnSel),2);
 
 
           
[r,p]=corrcoef(FixOff_Location_Average_All{i},FracOn_Location_Average_All{i});

r_LocCorrelation_Average{i}=r(1,2);
p_LocCorrelation_Average{i}=p(1,2);

 %Direct correlation between response before fixoff and frac on for the
 %first object
 
 FixOff_Location_Average_First{i}=nanmean(MeanPSTH_FR_FixOff_currBlock(TargetOrder_Line_curr==1,FixOffSel),2);
          
 FracOn_Location_Average_First{i}=nanmean(MeanPSTH_FR_FracOn_currBlock(TargetOrder_Line_curr==1,FracOnSel),2);
       
[r_tmp,p_tmp]=corrcoef(FixOff_Location_Average_First{i},FracOn_Location_Average_First{i});

r_LocCorrelation_Average_First{i}=r_tmp(1,2);
p_LocCorrelation_Average_First{i}=p_tmp(1,2);

%Average the PSTH according to the average responses

FracOnResp=FracOn_Location_Average_Mean(i,:);
[Val,Order]=sort(FracOnResp,'descend');

Highest=Order(1:ChooseAverageLoc);
Lowest=Order(end-ChooseAverageLoc+1:end);

HighestAngle=TargetPolarAngle_curr_unique{i}(Highest);
LowestAngle=TargetPolarAngle_curr_unique{i}(Lowest);

HighestAngle_Sel=ismember(TargetPolarAngle_curr,HighestAngle);
LowestAngle_Sel=ismember(TargetPolarAngle_curr,LowestAngle);


 MeanPSTH_currBlock_FixOff_HighestLocation{i}=nanmean(MeanPSTH_FR_FixOff_currBlock(HighestAngle_Sel,:),1);
 SemPSTH_currBlock_FixOff_HighestLocation{i}=nanstd(MeanPSTH_FR_FixOff_currBlock(HighestAngle_Sel,:),[],1)./sqrt(sum(~isnan(MeanPSTH_FR_FixOff_currBlock(HighestAngle_Sel,:)),1));
          
 MeanPSTH_currBlock_FracOn_HighestLocation{i}=nanmean(MeanPSTH_FR_FracOn_currBlock(HighestAngle_Sel,:),1);
 SemPSTH_currBlock_FracOn_HighestLocation{i}=nanstd(MeanPSTH_FR_FracOn_currBlock(HighestAngle_Sel,:),[],1)./sqrt(sum(~isnan(MeanPSTH_FR_FracOn_currBlock(HighestAngle_Sel,:)),1));
           

MeanPSTH_currBlock_FixOff_LowestLocation{i}=nanmean(MeanPSTH_FR_FixOff_currBlock(LowestAngle_Sel,:),1);
 SemPSTH_currBlock_FixOff_LowestLocation{i}=nanstd(MeanPSTH_FR_FixOff_currBlock(LowestAngle_Sel,:),[],1)./sqrt(sum(~isnan(MeanPSTH_FR_FixOff_currBlock(LowestAngle_Sel,:)),1));
          
 MeanPSTH_currBlock_FracOn_LowestLocation{i}=nanmean(MeanPSTH_FR_FracOn_currBlock(LowestAngle_Sel,:),1);
 SemPSTH_currBlock_FracOn_LowestLocation{i}=nanstd(MeanPSTH_FR_FracOn_currBlock(LowestAngle_Sel,:),[],1)./sqrt(sum(~isnan(MeanPSTH_FR_FracOn_currBlock(LowestAngle_Sel,:)),1));
           

%Average the PSTH according to the average responses for the first object

FracOnResp_First=FracOn_Location_Average_First_Mean(i,:);
[Val,Order]=sort(FracOnResp_First,'descend');

Highest_First=Order(1:ChooseAverageLoc);
Lowest_First=Order(end-ChooseAverageLoc+1:end);

HighestAngle_First=TargetPolarAngle_curr_unique{i}(Highest_First);
LowestAngle_First=TargetPolarAngle_curr_unique{i}(Lowest_First);

HighestAngle_Sel_First=ismember(TargetPolarAngle_curr,HighestAngle_First);
LowestAngle_Sel_First=ismember(TargetPolarAngle_curr,LowestAngle_First);


 MeanPSTH_currBlock_FixOff_HighestLocation_First{i}=nanmean(MeanPSTH_FR_FixOff_currBlock(HighestAngle_Sel_First,:),1);
 SemPSTH_currBlock_FixOff_HighestLocation_First{i}=nanstd(MeanPSTH_FR_FixOff_currBlock(HighestAngle_Sel_First,:),[],1)./sqrt(sum(~isnan(MeanPSTH_FR_FixOff_currBlock(HighestAngle_Sel_First,:)),1));
          
 MeanPSTH_currBlock_FracOn_HighestLocation_First{i}=nanmean(MeanPSTH_FR_FracOn_currBlock(HighestAngle_Sel_First,:),1);
 SemPSTH_currBlock_FracOn_HighestLocation_First{i}=nanstd(MeanPSTH_FR_FracOn_currBlock(HighestAngle_Sel_First,:),[],1)./sqrt(sum(~isnan(MeanPSTH_FR_FracOn_currBlock(HighestAngle_Sel_First,:)),1));
           

MeanPSTH_currBlock_FixOff_LowestLocation_First{i}=nanmean(MeanPSTH_FR_FixOff_currBlock(LowestAngle_Sel_First,:),1);
 SemPSTH_currBlock_FixOff_LowestLocation_First{i}=nanstd(MeanPSTH_FR_FixOff_currBlock(LowestAngle_Sel_First,:),[],1)./sqrt(sum(~isnan(MeanPSTH_FR_FixOff_currBlock(LowestAngle_Sel_First,:)),1));
          
 MeanPSTH_currBlock_FracOn_LowestLocation_First{i}=nanmean(MeanPSTH_FR_FracOn_currBlock(LowestAngle_Sel_First,:),1);
 SemPSTH_currBlock_FracOn_LowestLocation_First{i}=nanstd(MeanPSTH_FR_FracOn_currBlock(LowestAngle_Sel_First,:),[],1)./sqrt(sum(~isnan(MeanPSTH_FR_FracOn_currBlock(LowestAngle_Sel_First,:)),1));
           

%{

%Good Object

%After Fractal Onset   
SelGood=currBlocks&GoodObj_Line;
SelBad=currBlocks&BadObj_Line;
SponGoodtmp=Mean_Spon_ForEachTrial_Line(SelGood);
SponBadtmp=Mean_Spon_ForEachTrial_Line(SelBad);

PSTH_GoodFractal_currBlock=MeanPSTH_FR_Fractal(SelGood,:)-SponGoodtmp;
PSTH_Time_GoodFractal_currBlock=MeanPSTH_Time_Fractal(SelGood,:);

PSTH_group_temp=SortDataAccordingToGroup(PSTH_GoodFractal_currBlock,SearchingTime_Line(SelGood),SepGroupAmp,'descend');


PSTH_RTGroup_GoodFractal_mean{i}=PSTH_group_temp.mean(1:SepGroupAmpNum,:);
PSTH_RTGroup_GoodFractal_sem{i}=PSTH_group_temp.sem(1:SepGroupAmpNum,:);


PSTH_RTGroup_GoodFractal_Time{i}=nanmean(PSTH_Time_GoodFractal_currBlock,1);
%Bad Object

%After Fractal Onset   
PSTH_BadFractal_currBlock=MeanPSTH_FR_Fractal(SelBad,:)-SponBadtmp;
PSTH_Time_BadFractal_currBlock=MeanPSTH_Time_Fractal(SelBad,:);

PSTH_group_temp=SortDataAccordingToGroup(PSTH_BadFractal_currBlock,SearchingTime_Line(SelBad),SepGroupAmp,'descend');


PSTH_RTGroup_BadFractal_mean{i}=PSTH_group_temp.mean(1:SepGroupAmpNum,:);
PSTH_RTGroup_BadFractal_sem{i}=PSTH_group_temp.sem(1:SepGroupAmpNum,:);


PSTH_RTGroup_BadFractal_Time{i}=nanmean(PSTH_Time_BadFractal_currBlock,1);
%}

 
  
 
  
end

%PlotStartOffset
if ShowFigureFlag
  FigureStartNum=140;
  FigureIndex=1;
    
       %Figure1 Response for different fractals in each location  for each mode
    figtitlestr{FigureIndex}='PSTH_FractalLocation_Different_Blocks';
    fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[20,100, 800,800],'Name',figtitlestr{FigureIndex});
  %  tiledlayout(NumberOfBlocks,2);
 
     for i=1:NumberOfBlocks
    
         
  
        subplot(NumberOfBlocks,2,(i-1)*2+1)
       % ax((i-1)*2+1) = nexttile;
       
       
           % Eye_V_Explore_currBlock_SacSep{i,j}
           
            currBlock=UniqueBlockType(i);
            BlockOrder(i)=find(AssociateBlock==currBlock);
               area(FixOffInterval,[max(max(MeanPSTH_currBlock_FixOff_HighestLocation{i}+SemPSTH_currBlock_FixOff_HighestLocation{i}))*1.2,max(max(MeanPSTH_currBlock_FixOff_HighestLocation{i}+SemPSTH_currBlock_FixOff_HighestLocation{i}))*1.2],...
                    'FaceColor','k','FaceAlpha',.1,'EdgeAlpha',.1,'HandleVisibility','off');
                
            for j=1:length( TargetPolarAngle_curr_unique{i})
                if ~isempty(MeanPSTH_currBlock_FixOffLocation{i,j})
                       %Add a shaded area to indicate analysis area
               
                hold on
             shadedErrorBar(PSTH_Time_FixOff, MeanPSTH_currBlock_FixOffLocation{i,j}, SemPSTH_currBlock_FixOffLocation{i,j},'lineprops',{Color{j}},'transparent',1,'patchSaturation',0.5);%PSTH for fix
             set(findall(gca, 'Type', 'Line'),'LineWidth',3);
             hold on
              xlim([-PlotStartOffset(1),PlotStopOffset(1)*1.5]);
                end
       
            end
            if i==1
              legend( string(round(TargetPolarAngle_curr_unique{i})));
            end
       
        xlabel('Time from Fix Off(ms)');
        ylabel('Firing Rate(Hz)over spon');
       % title(GroupStr{j});
        title(Blockstr{BlockOrder(i)});
         set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
         
      %  ax(i*2) = nexttile;
        subplot(NumberOfBlocks,2,i*2)
           %Add a shaded area to indicate analysis area
       area(FracOnInterval,[max(max(MeanPSTH_currBlock_FracOn_HighestLocation{i}+SemPSTH_currBlock_FracOn_HighestLocation{i}))*1.2,max(max(MeanPSTH_currBlock_FracOn_HighestLocation{i}+SemPSTH_currBlock_FracOn_HighestLocation{i}))*1.2],...
    'FaceColor','k','FaceAlpha',.1,'EdgeAlpha',.1,'HandleVisibility','off');
hold on
          for j=1:length( TargetPolarAngle_curr_unique{i})
              if ~isempty( MeanPSTH_currBlock_FracOnLocation{i,j})
            shadedErrorBar(PSTH_Time_FracOn, MeanPSTH_currBlock_FracOnLocation{i,j}, SemPSTH_currBlock_FracOnLocation{i,j},'lineprops',{Color{j}},'transparent',1,'patchSaturation',0.5);%PSTH for fix
             set(findall(gca, 'Type', 'Line'),'LineWidth',3);
             hold on
              xlim([-PlotStartOffset(2),nanmean(RelativeTimeMarker_FracOn(:,2))]);
              end
          end
         if i==1
              legend( string(round(TargetPolarAngle_curr_unique{i})));
         end
       
        xlabel('Time from Cue on(ms)');
        ylabel('Firing Rate(Hz)over spon');
       % title(GroupStr{j});
        title(Blockstr{BlockOrder(i)});
         set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
         
         
       
     end     
      
  % linkaxes(ax,'y');
   
     FigureStartNum=FigureStartNum+1;
     FigureIndex=FigureIndex+1;  
     
 
   %Figure2 Correlation
    figtitlestr{FigureIndex}='FixOffCueOnCorrelation_Different_Blocks';
    fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[70,100, 800,600],'Name',figtitlestr{FigureIndex});
  %  tiledlayout(ceil(NumberOfBlocks/2),2);
     for i=1:NumberOfBlocks
    
        subplot(NumberOfBlocks,2,2*(i-1)+1) 
  
      %  subplot(NumberOfBlocks,2,i)
       % ax(i) = nexttile;
    
        currBlock=UniqueBlockType(i);
        BlockOrder(i)=find(AssociateBlock==currBlock);
       
       plot( FixOff_Location_Average_Mean(i,:), FracOn_Location_Average_Mean(i,:),'.r','MarkerEdgeColor',ColorForEachMode(AssociateBlock==currBlock,:),'MarkerSize',30,'LineWidth',3);
       % axis equal;
        xlim([min(min(FixOff_Location_Average_Mean(i,:)), min(FracOn_Location_Average_Mean(i,:)))*0.8,max(max(FixOff_Location_Average_Mean(i,:)), max(FracOn_Location_Average_Mean(i,:)))*1.2]);
       
       ylim([min(min(FixOff_Location_Average_Mean(i,:)), min(FracOn_Location_Average_Mean(i,:)))*0.8,max(max(FixOff_Location_Average_Mean(i,:)), max(FracOn_Location_Average_Mean(i,:)))*1.2]);
         xlabel('FixOff Response Average');
         ylabel('FracOn Response Average');
         set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
         box off;
         
         
        subplot(NumberOfBlocks,2,2*i) ;
        plot(FixOff_Location_Average_All{i}, FracOn_Location_Average_All{i},'.r','MarkerEdgeColor',ColorForEachMode(AssociateBlock==currBlock,:),'MarkerSize',30,'LineWidth',3);
       % axis equal;
       % xlim([min(min(FixOff_Location_Average_Mean(i,:)), min(FracOn_Location_Average_Mean(i,:)))*0.8,max(max(FixOff_Location_Average_Mean(i,:)), max(FracOn_Location_Average_Mean(i,:)))*1.2]);
       
      % ylim([min(min(FixOff_Location_Average_Mean(i,:)), min(FracOn_Location_Average_Mean(i,:)))*0.8,max(max(FixOff_Location_Average_Mean(i,:)), max(FracOn_Location_Average_Mean(i,:)))*1.2]);
         xlabel('FixOff Response Average');
         ylabel('FracOn Response Average');
         set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
         box off;
         
         ldgstr=sprintf('r=%2.2d\np=%2.2d',r_LocCorrelation_Average{i},p_LocCorrelation_Average{i});
         
         legend(ldgstr);

      
         
         
       
     end     
      
   %linkaxes(ax,'y');
   
     FigureStartNum=FigureStartNum+1;
     FigureIndex=FigureIndex+1;  
        
 
   %Figure3 Polar plot of the tuning curve
    figtitlestr{FigureIndex}='PolarTuningCurveCompare';
    fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[90,100, 800,800],'Name',figtitlestr{FigureIndex});
    %tiledlayout(NumberOfBlocks,2);
     for i=1:NumberOfBlocks
    
         
  
        subplot(NumberOfBlocks,2,(i-1)*2+1)
       % ax((i-1)*2+1) = nexttile;
    
        currBlock=UniqueBlockType(i);
        BlockOrder(i)=find(AssociateBlock==currBlock);
        
         angle=TargetPolarAngle_curr_unique{i}';
        polarwitherrorbar([angle/180*pi angle(1)/180*pi],[FixOff_Location_Average_Mean(i,:) FixOff_Location_Average_Mean(i,1)],[FixOff_Location_Average_Sem(i,:) FixOff_Location_Average_Sem(i,1)],...
            '-b',3);
        hold on
     %   h_p = polarplot([PrefVector(j) PrefVector(j)]/180*pi,[0 PrefVector_Amp(j)],'r');
     %   set(h_p,'LineWidth',2.5);
        
        set(gca,'ThetaDir' , 'counterclockwise');
        set(gca,'ThetaZeroLocation','top')
        set(gca,'FontSize',25,'FontWeight','Bold','LineWidth',3);
        
        thetaticks([0:45:359]);
        thetaticklabels([0:45:359]);
        title(strcat('FixOff&',Blockstr{BlockOrder(i)}),'FontSize',15,'FontWeight','Bold');
        
%        ax(i*2) = nexttile;
         subplot(NumberOfBlocks,2,i*2)
        currBlock=UniqueBlockType(i);
        BlockOrder(i)=find(AssociateBlock==currBlock);
        
         angle=TargetPolarAngle_curr_unique{i}';
        polarwitherrorbar([angle/180*pi angle(1)/180*pi],[FracOn_Location_Average_Mean(i,:) FracOn_Location_Average_Mean(i,1)],[FracOn_Location_Average_Sem(i,:) FracOn_Location_Average_Sem(i,1)],...
            '-b',3);
        hold on
     %   h_p = polarplot([PrefVector(j) PrefVector(j)]/180*pi,[0 PrefVector_Amp(j)],'r');
     %   set(h_p,'LineWidth',2.5);
        
        set(gca,'ThetaDir' , 'counterclockwise');
        set(gca,'ThetaZeroLocation','top')
        set(gca,'FontSize',25,'FontWeight','Bold','LineWidth',3);
        
        thetaticks([0:45:359]);
        thetaticklabels([0:45:359]);
       title(strcat('CueOn',Blockstr{BlockOrder(i)}),'FontSize',15,'FontWeight','Bold');
        
        
     end     
      
   
   
     FigureStartNum=FigureStartNum+1;
     FigureIndex=FigureIndex+1;     
     
     
 %Figure5 PSTH between highest and lowest response location 
    figtitlestr{FigureIndex}='PSTH_HighestLowest_LocationCompare';
    fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[110,100, 800,800],'Name',figtitlestr{FigureIndex});
    %tiledlayout(NumberOfBlocks,2);
     for i=1:NumberOfBlocks
    
         
  
        subplot(NumberOfBlocks,2,(i-1)*2+1)
       % ax((i-1)*2+1) = nexttile;
    
        currBlock=UniqueBlockType(i);
        BlockOrder(i)=find(AssociateBlock==currBlock);
        %Add a shaded area to indicate analysis area
       area(FixOffInterval,[max(max(MeanPSTH_currBlock_FixOff_HighestLocation{i}+SemPSTH_currBlock_FixOff_HighestLocation{i}))*1.2,max(max(MeanPSTH_currBlock_FixOff_HighestLocation{i}+SemPSTH_currBlock_FixOff_HighestLocation{i}))*1.2],...
    'FaceColor','k','FaceAlpha',.1,'EdgeAlpha',.1,'HandleVisibility','off');
        hold on
         shadedErrorBar(PSTH_Time_FixOff,  MeanPSTH_currBlock_FixOff_HighestLocation{i},  SemPSTH_currBlock_FixOff_HighestLocation{i},'lineprops',{Color{1}},'transparent',1,'patchSaturation',0.5);%PSTH for fix
         
         shadedErrorBar(PSTH_Time_FixOff,  MeanPSTH_currBlock_FixOff_LowestLocation{i},  SemPSTH_currBlock_FixOff_LowestLocation{i},'lineprops',{Color{2}},'transparent',1,'patchSaturation',0.5);%PSTH for fix
         box off;
         set(findall(gca, 'Type', 'Line'),'LineWidth',3);
         xlim([-PlotStartOffset(1),PlotStopOffset(1)*1.5]);
         xlabel('Time from fix off');
         ylabel('Firing Rate');
         title(Blockstr{BlockOrder(i)});
         set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
         
        
        
%        ax(i*2) = nexttile;
         subplot(NumberOfBlocks,2,i*2)
        currBlock=UniqueBlockType(i);
        BlockOrder(i)=find(AssociateBlock==currBlock);
        %Add a shaded area to indicate analysis area
        
       area(FracOnInterval,[max(max(MeanPSTH_currBlock_FracOn_HighestLocation{i}+SemPSTH_currBlock_FracOn_HighestLocation{i}))*1.2,max(max(MeanPSTH_currBlock_FracOn_HighestLocation{i}+ SemPSTH_currBlock_FracOn_HighestLocation{i}))*1.2],...
    'FaceColor','k','FaceAlpha',.1,'EdgeAlpha',.1,'HandleVisibility','off');
         hold on
        shadedErrorBar(PSTH_Time_FracOn,  MeanPSTH_currBlock_FracOn_HighestLocation{i},  SemPSTH_currBlock_FracOn_HighestLocation{i},'lineprops',{Color{1}},'transparent',1,'patchSaturation',0.5);%PSTH for fix
         shadedErrorBar(PSTH_Time_FracOn,  MeanPSTH_currBlock_FracOn_LowestLocation{i},  SemPSTH_currBlock_FracOn_LowestLocation{i},'lineprops',{Color{2}},'transparent',1,'patchSaturation',0.5);%PSTH for fix
         box off;
         set(findall(gca, 'Type', 'Line'),'LineWidth',3);
         xlim([-PlotStartOffset(2),600+PlotStopOffset(2)*1.5]);
         ylim([0,40]);
         xlabel('Time from cue on');
         ylabel('Firing Rate');
         legend({'HighestLoc','Lowest Loc'});
           title(Blockstr{BlockOrder(i)});
         set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
        
      
        
        
     end     
      
   
   
     FigureStartNum=FigureStartNum+1;
     FigureIndex=FigureIndex+1;     
         
     
     
     
       %Figure6 Response for first fractals in each location  for each
       %mode 
    figtitlestr{FigureIndex}='PSTH_FirstFractalLocation_Different_Blocks';
    fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[130,100, 800,800],'Name',figtitlestr{FigureIndex});
  %  tiledlayout(NumberOfBlocks,2);
 
     for i=1:NumberOfBlocks
    
         
  
        subplot(NumberOfBlocks,2,(i-1)*2+1)
       % ax((i-1)*2+1) = nexttile;
       
       
           % Eye_V_Explore_currBlock_SacSep{i,j}
           
            currBlock=UniqueBlockType(i);
            BlockOrder(i)=find(AssociateBlock==currBlock);
               area(FixOffInterval,[max(max(MeanPSTH_currBlock_FixOff_HighestLocation_First{i}+SemPSTH_currBlock_FixOff_HighestLocation_First{i}))*1.2,max(max(MeanPSTH_currBlock_FixOff_HighestLocation_First{i}+SemPSTH_currBlock_FixOff_HighestLocation_First{i}))*1.2],...
                    'FaceColor','k','FaceAlpha',.1,'EdgeAlpha',.1,'HandleVisibility','off');
                
            for j=1:length( TargetPolarAngle_curr_unique{i})
                if ~isempty(MeanPSTH_currBlock_FixOffLocation_first{i,j})
                       %Add a shaded area to indicate analysis area
               
                hold on
             shadedErrorBar(PSTH_Time_FixOff, MeanPSTH_currBlock_FixOffLocation_first{i,j}, SemPSTH_currBlock_FixOffLocation_first{i,j},'lineprops',{Color{j}},'transparent',1,'patchSaturation',0.5);%PSTH for fix
             set(findall(gca, 'Type', 'Line'),'LineWidth',3);
             hold on
              xlim([-PlotStartOffset(1),PlotStopOffset(1)*1.5]);
                end
                yrange=ylim;
                plot([SceneOnTime_Line_first(i),SceneOnTime_Line_first(i)],yrange,'--k');
                
                 
       
            end
            if i==1
              legend( string(round(TargetPolarAngle_curr_unique{i})));
            end
       
        xlabel('Time from First Fix Off(ms)');
        ylabel('Firing Rate(Hz)over spon');
       % title(GroupStr{j});
        title(Blockstr{BlockOrder(i)});
         set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
         
      %  ax(i*2) = nexttile;
        subplot(NumberOfBlocks,2,i*2)
           %Add a shaded area to indicate analysis area
       area(FracOnInterval,[max(max(MeanPSTH_currBlock_FracOn_HighestLocation_First{i}+SemPSTH_currBlock_FracOn_HighestLocation_First{i}))*1.2,max(max(MeanPSTH_currBlock_FracOn_HighestLocation_First{i}+SemPSTH_currBlock_FracOn_HighestLocation_First{i}))*1.2],...
    'FaceColor','k','FaceAlpha',.1,'EdgeAlpha',.1,'HandleVisibility','off');
hold on
          for j=1:length( TargetPolarAngle_curr_unique{i})
              if ~isempty( MeanPSTH_currBlock_FracOnLocation_first{i,j})
            shadedErrorBar(PSTH_Time_FracOn, MeanPSTH_currBlock_FracOnLocation_first{i,j}, SemPSTH_currBlock_FracOnLocation_first{i,j},'lineprops',{Color{j}},'transparent',1,'patchSaturation',0.5);%PSTH for fix
             set(findall(gca, 'Type', 'Line'),'LineWidth',3);
             hold on
              xlim([-PlotStartOffset(2),nanmean(RelativeTimeMarker_FracOn(:,2))]);
              end
          end
         if i==1
              legend( string(round(TargetPolarAngle_curr_unique{i})));
         end
       
        xlabel('Time from First Cue on(ms)');
        ylabel('Firing Rate(Hz)over spon');
       % title(GroupStr{j});
        title(Blockstr{BlockOrder(i)});
         set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
         
         
       
     end     
      
  % linkaxes(ax,'y');
   
     FigureStartNum=FigureStartNum+1;
     FigureIndex=FigureIndex+1;  
     
     %Figure7 Correlation for the firsty object
    figtitlestr{FigureIndex}='FirstFixOffCueOnCorrelation_Different_Blocks';
    fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[140,100, 800,600],'Name',figtitlestr{FigureIndex});
  %  tiledlayout(ceil(NumberOfBlocks/2),2);
     for i=1:NumberOfBlocks
    
        subplot(NumberOfBlocks,2,2*(i-1)+1) 
  
      %  subplot(NumberOfBlocks,2,i)
       % ax(i) = nexttile;
    
        currBlock=UniqueBlockType(i);
        BlockOrder(i)=find(AssociateBlock==currBlock);
       
       plot( FixOff_Location_Average_First_Mean(i,:), FracOn_Location_Average_First_Mean(i,:),'.r','MarkerEdgeColor',ColorForEachMode(AssociateBlock==currBlock,:),'MarkerSize',30,'LineWidth',3);
       % axis equal;
        xlim([min(min(FixOff_Location_Average_First_Mean(i,:)), min(FracOn_Location_Average_First_Mean(i,:)))*0.8,max(max(FixOff_Location_Average_First_Mean(i,:)), max(FracOn_Location_Average_First_Mean(i,:)))*1.2]);
       
       ylim([min(min(FixOff_Location_Average_First_Mean(i,:)), min(FracOn_Location_Average_First_Mean(i,:)))*0.8,max(max(FixOff_Location_Average_First_Mean(i,:)), max(FracOn_Location_Average_First_Mean(i,:)))*1.2]);
       
       
       
       xlabel('FixOff Response Average(First)');
         ylabel('FracOn Response Average(First)');
         set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
         box off;
         
         
        subplot(NumberOfBlocks,2,2*i) ;
        plot(FixOff_Location_Average_First{i}, FracOn_Location_Average_First{i},'.r','MarkerEdgeColor',ColorForEachMode(AssociateBlock==currBlock,:),'MarkerSize',30,'LineWidth',3);
       % axis equal;
       % xlim([min(min(FixOff_Location_Average_Mean(i,:)), min(FracOn_Location_Average_Mean(i,:)))*0.8,max(max(FixOff_Location_Average_Mean(i,:)), max(FracOn_Location_Average_Mean(i,:)))*1.2]);
       
      % ylim([min(min(FixOff_Location_Average_Mean(i,:)), min(FracOn_Location_Average_Mean(i,:)))*0.8,max(max(FixOff_Location_Average_Mean(i,:)), max(FracOn_Location_Average_Mean(i,:)))*1.2]);
         xlabel('FixOff Response Average(First)');
         ylabel('FracOn Response Average(First)');
         set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
         box off;
         
         ldgstr=sprintf('r=%2.2d\np=%2.2d',r_LocCorrelation_Average_First{i},p_LocCorrelation_Average_First{i});
         
         legend(ldgstr);

      
         
         
       
     end     
      
   %linkaxes(ax,'y');
   
     FigureStartNum=FigureStartNum+1;
     FigureIndex=FigureIndex+1;  
     
     %Figure8 Polar plot of the tuning curve for the first object
    figtitlestr{FigureIndex}='PolarTuningCurveCompare_First';
    fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[150,100, 800,800],'Name',figtitlestr{FigureIndex});
    %tiledlayout(NumberOfBlocks,2);
     for i=1:NumberOfBlocks
    
         
  
        subplot(NumberOfBlocks,2,(i-1)*2+1)
       % ax((i-1)*2+1) = nexttile;
    
        currBlock=UniqueBlockType(i);
        BlockOrder(i)=find(AssociateBlock==currBlock);
        
         angle=TargetPolarAngle_curr_unique{i}';
        polarwitherrorbar([angle/180*pi angle(1)/180*pi],[FixOff_Location_Average_First_Mean(i,:) FixOff_Location_Average_First_Mean(i,1)],[FixOff_Location_Average_First_Sem(i,:) FixOff_Location_Average_First_Sem(i,1)],...
            '-b',3);
        hold on
     %   h_p = polarplot([PrefVector(j) PrefVector(j)]/180*pi,[0 PrefVector_Amp(j)],'r');
     %   set(h_p,'LineWidth',2.5);
        
        set(gca,'ThetaDir' , 'counterclockwise');
        set(gca,'ThetaZeroLocation','top')
        set(gca,'FontSize',25,'FontWeight','Bold','LineWidth',3);
        
        thetaticks([0:45:359]);
        thetaticklabels([0:45:359]);
        title(strcat('First FixOff&',Blockstr{BlockOrder(i)}),'FontSize',15,'FontWeight','Bold');
        
%        ax(i*2) = nexttile;
         subplot(NumberOfBlocks,2,i*2)
        currBlock=UniqueBlockType(i);
        BlockOrder(i)=find(AssociateBlock==currBlock);
        
         angle=TargetPolarAngle_curr_unique{i}';
        polarwitherrorbar([angle/180*pi angle(1)/180*pi],[FracOn_Location_Average_First_Mean(i,:) FracOn_Location_Average_First_Mean(i,1)],[FracOn_Location_Average_First_Sem(i,:) FracOn_Location_Average_First_Sem(i,1)],...
            '-b',3);
        hold on
     %   h_p = polarplot([PrefVector(j) PrefVector(j)]/180*pi,[0 PrefVector_Amp(j)],'r');
     %   set(h_p,'LineWidth',2.5);
        
        set(gca,'ThetaDir' , 'counterclockwise');
        set(gca,'ThetaZeroLocation','top')
        set(gca,'FontSize',25,'FontWeight','Bold','LineWidth',3);
        
        thetaticks([0:45:359]);
        thetaticklabels([0:45:359]);
       title(strcat('First CueOn',Blockstr{BlockOrder(i)}),'FontSize',15,'FontWeight','Bold');
        
        
     end     
      
   
   
     FigureStartNum=FigureStartNum+1;
     FigureIndex=FigureIndex+1;     
     
     
 %Figure9 PSTH between highest and lowest response location 
    figtitlestr{FigureIndex}='PSTH_HighestLowest_LocationCompare_First';
    fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[170,100, 800,800],'Name',figtitlestr{FigureIndex});
    %tiledlayout(NumberOfBlocks,2);
     for i=1:NumberOfBlocks
    
         
  
        subplot(NumberOfBlocks,2,(i-1)*2+1)
       % ax((i-1)*2+1) = nexttile;
    
        currBlock=UniqueBlockType(i);
        BlockOrder(i)=find(AssociateBlock==currBlock);
        %Add a shaded area to indicate analysis area
       area(FixOffInterval,[max(max(MeanPSTH_currBlock_FixOff_HighestLocation_First{i}+SemPSTH_currBlock_FixOff_HighestLocation_First{i}))*1.2,max(max(MeanPSTH_currBlock_FixOff_HighestLocation_First{i}+SemPSTH_currBlock_FixOff_HighestLocation_First{i}))*1.2],...
    'FaceColor','k','FaceAlpha',.1,'EdgeAlpha',.1,'HandleVisibility','off');
        hold on
         shadedErrorBar(PSTH_Time_FixOff,  MeanPSTH_currBlock_FixOff_HighestLocation_First{i},  SemPSTH_currBlock_FixOff_HighestLocation_First{i},'lineprops',{Color{1}},'transparent',1,'patchSaturation',0.5);%PSTH for fix
         
         shadedErrorBar(PSTH_Time_FixOff,  MeanPSTH_currBlock_FixOff_LowestLocation_First{i},  SemPSTH_currBlock_FixOff_LowestLocation_First{i},'lineprops',{Color{2}},'transparent',1,'patchSaturation',0.5);%PSTH for fix
         
         box off;
         set(findall(gca, 'Type', 'Line'),'LineWidth',3);
         xlim([-PlotStartOffset(1),PlotStopOffset(1)*1.5]);
         %Add scene on line marker 
       yrange=ylim;
         plot([SceneOnTime_Line_first(i),SceneOnTime_Line_first(i)],yrange,'--k');
         
         xlabel('Time from fix off');
         ylabel('Firing Rate');
         title(Blockstr{BlockOrder(i)});
         set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
         
        
        
%        ax(i*2) = nexttile;
         subplot(NumberOfBlocks,2,i*2)
        currBlock=UniqueBlockType(i);
        BlockOrder(i)=find(AssociateBlock==currBlock);
        %Add a shaded area to indicate analysis area
        
       area(FracOnInterval,[max(max(MeanPSTH_currBlock_FracOn_HighestLocation_First{i}+SemPSTH_currBlock_FracOn_HighestLocation_First{i}))*1.2,max(max(MeanPSTH_currBlock_FracOn_HighestLocation_First{i}+ SemPSTH_currBlock_FracOn_HighestLocation_First{i}))*1.2],...
    'FaceColor','k','FaceAlpha',.1,'EdgeAlpha',.1,'HandleVisibility','off');
         hold on
        shadedErrorBar(PSTH_Time_FracOn,  MeanPSTH_currBlock_FracOn_HighestLocation_First{i},  SemPSTH_currBlock_FracOn_HighestLocation_First{i},'lineprops',{Color{1}},'transparent',1,'patchSaturation',0.5);%PSTH for fix
         shadedErrorBar(PSTH_Time_FracOn,  MeanPSTH_currBlock_FracOn_LowestLocation_First{i},  SemPSTH_currBlock_FracOn_LowestLocation_First{i},'lineprops',{Color{2}},'transparent',1,'patchSaturation',0.5);%PSTH for fix
         box off;
         set(findall(gca, 'Type', 'Line'),'LineWidth',3);
         xlim([-PlotStartOffset(2),600+PlotStopOffset(2)*1.5]);
         ylim([0,40]);
         xlabel('Time from cue on');
         ylabel('Firing Rate');
         legend({'HighestLoc','Lowest Loc'});
           title(Blockstr{BlockOrder(i)});
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