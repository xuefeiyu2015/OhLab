function  ExploreNeuralAnalysis(Data, StartTrial, EndTrial,ShowFigureFlag,OutputFlag);
%Function to analyze neural responses during exploratory behavior for both
%fix loc series mode and rand loc series mode
%by Xuefei Yu 08222020
ECode;
TaskType=Data.TaskType;
TrialNum=EndTrial-StartTrial+1;
Selection=StartTrial:EndTrial;%Select the trials according to the pannel
%Abnormal check %Especially For Robin081120FORAGEXPLO4.rst
CheckCode=SCENE_ON;%This code should appear once in each trial
AbnormalTrials=AbnormalTrialCheck(Data,CheckCode);
Selection=Selection(~ismember(Selection,AbnormalTrials));

%User defined parameter for the neural/behavioral analysis
TargetOnOffset=0;%in ms, the time selection bin after target onset; 
StartDelay=0;%in ms, the time selection bin after fixation offset;

%Determine which figure to check
PlotSpikeEventOn=1;%Former ~4 type figures 
PlotSpikeDensityOn=1;


%For segmentation of  the neural responses 
%SegEvent=[EXTRAFIXSTARTCD,SCENE_ON,SHOWFIXCD,FIXOFF,TGTCD];%Event Flow of interest
%SegEventEnd=[SHOWFIXCD,FIXOFF,TGTOFF,RWDOFFCD ,NaN,EXTRAFP];%Event end
%SegEventEnd=[EXTRAFIXOFF,SHOWFIXCD,FIXOFF,TGTCD,TGTOFF];%Event end
%SegEventStr={'ExtraFPOn','SceneOn','FixOn','Explore','ChoicePeriod'};

%PlotStartOffset=[200,200,200,200,200];%in ms
%PlotStopOffset=[200,200,200,200,300];%in ms

%EventTypeFlag=[1,1,2,2,2,2];%1:Number same with the trial number, 2:Number same with the target number

SegEvent=[FIXOFF];%Event Flow of interest

SegEventEnd=[TGTCD];%Event end; will change to cue on automatically later

SwitchToCueOn=1;

%SegEventEnd=[CUEON];%Event end
SegEventStr={'Explore'};
PlotStartOffset=[-StartDelay];%in ms
PlotStopOffset=[TargetOnOffset];%in ms

EventTypeFlag=[2];%1:Number same with the trial number, 2:Number same with the target number

%Parameter to define the sliding window to measure firing rates

BinWidth=15;%in ms

%Parameter for raster plot
Barlength=1;
LineWidth=1;
%Parameter for visualize spike event
SpikeSizePlot=10;
%Parameter for 2D histogram for spikes
%2D Density map
TwoDBinRange=[-30,30;-30,30];%Width;Height
BinDim=[2;2];%Width;Height
Param.Range=TwoDBinRange;
Param.Bin=BinDim;
%Param.FR_Threshold=100;
Param.FR_Threshold=20;
Param.GaussianSmooth=0;
Param.GaussianSigma=[1.5;1.5];


%ResponseLatencyThreshold=3;

%Some user defined parameter for visualization
%Color Definition for each scenes
Color={'r','g','b','m',[0.4,0.4,0],[0.28,0.78,0.89],'k',[0.6,0.1,1],[0.6,0.1,0.4],[0.3,0.5,1],[0.6,0.5,0.1],[0,0.8,0.4]};
ColorForEachMode=[0.83 0.14 0.14
             1.00 0.54 0.00
             0.47 0.25 0.80
             0.25 0.80 0.54];
AssociateBlock=[10,11,0,1];
Blockstr={'FixLoc,FixSceneSeq','FixLoc,RandSceneSeq','RandLoc,FixSceneSeq','RandLoc,RandSceneSeq'};

ColorFractal=[0.8500, 0.3250, 0.0980];
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

Param.OriBinWidth=EyeBinWidth;

StepSize=EyeBinWidth;%in ms %In order to organize spike count into the same interval

%Load event data
EventChannel=Data.EventChannel(Selection,:);
EventTimeChannel=Data.EventTimeChannel(Selection,:);
EventBin=Data.EventBin(Selection,:);


if SwitchToCueOn
    CueOnCheck=FindOutTime(EventChannel,EventBin,CUEON);
    if sum(isnan(CueOnCheck(:,1)))~=size(CueOnCheck,1)
        SegEventEnd=[CUEON];
    else
        disp('No Cue, TGTON still');
    end
end

%Get the information about the extra FP

%Recalibrate eyes according to the extra FP
ParaLib=Data.ParaLib;
ExtraFP_Info=ParaLib('ExtraFixPointID,Style,Size,Loc_X,Loc_Y');
ExtraFP_Loc=ExtraFP_Info(Selection,end-1:end);


%Get the information about the extra FP
%ExtraFPOnset=FindOutTime(EventChannel,EventTimeChannel,EXTRAFP);
ExtraEyeIn=FindOutTime(EventChannel,EventBin,EXTRAFIXSTARTCD);
ExtraEyeOff=FindOutTime(EventChannel,EventBin,EXTRAFIXOFF);

EyeInBin_X=SelectEyeInterval(EyeChannel_X_L,ExtraEyeIn,ExtraEyeOff);
EyeInBin_X=EyeInBin_X{1};
EyeInBin_Y=SelectEyeInterval(EyeChannel_Y_L,ExtraEyeIn,ExtraEyeOff);
EyeInBin_Y=EyeInBin_Y{1};
EyeDataForCal={EyeInBin_X,EyeInBin_Y};

EyeDataRaw={EyeChannel_X_L,EyeChannel_Y_L};
CalPlotCheck=0;%Examine the calibration results
CaledEye=EyeCalibration(EyeDataRaw,EyeDataForCal,EyeCalInfo_L,ExtraFP_Loc,CalPlotCheck);


EyeChannel_X_L=CaledEye.EyeTrace{1};
EyeChannel_Y_L=CaledEye.EyeTrace{2};

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
 %screen out the searching trials
 FixOffTime=FixOffTime(AdvanceInTargetTrial,:);
 OutScreenTime=OutScreenTime(AdvanceInTargetTrial,:);
 
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
NumberOfBlocks=max(reprojectorderindex);
 
 %If the target onset time and fixofftime dimension are not aligned 
 if size(TargetOnTime,2)~=size(FixOffTime,2)
     FixOffTime=FixOffTime(:,1:size(TargetOnTime,2));
     disp('Warning!FixOffTime Corrected!!')
 end
     
 
SearchingTime= (TargetOnTime-FixOffTime)/1000;%Change into s for presentation


%Segregate spike data from Fixation onset to target onset
SegmentData=EventSegTool(SegEvent,SegEventEnd,PlotStartOffset,PlotStopOffset,EventChannel,EventTimeChannel,SpikeTime,StepSize,BinWidth);

%Retrive data from the field
%Total average for each event
MeanFR=SegmentData.MeanFR; %[Mean,sem]%Averall mean and sem
MeanFRRaw=SegmentData.MeanFRRaw;%=AverageFiringRate_Event;

%PSTH(Spike Density Function) for each event
MeanPSTH=SegmentData.MeanPSTH;%=[PSTH_Event_Time',PSTH_Event'];
MeanPSTH_Time=MeanPSTH(:,1);
MeanPSTH_FR=MeanPSTH(:,2);

%Raster and raw spike count PSTH
AlignedRaster=SegmentData.RasterAfterAlign;
PSTH_Event_SpikeCount=SegmentData.PSTH_Event_SpikeCount;


%Relative time bin for each event
RelativeTimeMarker=SegmentData.RelativeTimeMarker;%=EventTimeMarker;
%Absolute time marker for each event
AbosoluteTimeMarker=SegmentData.AbsoluteTimeMarker;


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

%Target Location
TargetLoc_Line_X=squeeze(TargetLoc(:,:,1));
TargetLoc_Line_Y=squeeze(TargetLoc(:,:,2));

TargetLoc_Line_X=reshape(TargetLoc_Line_X,numel(TargetLoc_Line_X),1);
TargetLoc_Line_X=TargetLoc_Line_X(NoNNaNSelector);

TargetLoc_Line_Y=reshape(TargetLoc_Line_Y,numel(TargetLoc_Line_Y),1);
TargetLoc_Line_Y=TargetLoc_Line_Y(NoNNaNSelector);

TargetLoc_Line=[TargetLoc_Line_X,TargetLoc_Line_Y];

%Target Angle
TargetPolarAngle_Line=TargetPolarAngle;
TargetPolarAngle_Line=reshape(TargetPolarAngle_Line,numel(TargetPolarAngle_Line),1);
TargetPolarAngle_Line=TargetPolarAngle_Line(NoNNaNSelector);

%TargetNumLine
TargetNum_Line=repmat(TargetNum,1,max(TargetNum));
TargetNum_Line=reshape(TargetNum_Line,numel(TargetNum_Line),1);
TargetNum_Line=TargetNum_Line(NoNNaNSelector);

%Rearrange block type into lines
BlockType_Line=repmat(BlockType,1,max(TargetNum));
BlockType_Line=reshape(BlockType_Line,numel(BlockType_Line),1);
BlockType_Line=BlockType_Line(NoNNaNSelector);
%Rearrange Scene type into lines
SceneName_Line=repmat(SceneName,1,max(TargetNum));
SceneName_Line=reshape(SceneName_Line,numel(SceneName_Line),1);
SceneName_Line=SceneName_Line(NoNNaNSelector);
SceneName_Unique=unique(SceneName_Line);

%Select out eye trace for each target
StartTimeBin=floor((FixOffTime+StartDelay)/EyeBinWidth);         
EndTimeBin=ceil((TargetOnTime+TargetOnOffset)/EyeBinWidth);
         
Eye_Explore_X=SelectEyeInterval(EyeChannel_X_L,StartTimeBin,EndTimeBin)';
Eye_Explore_X_Line=ReorganizeEye(Eye_Explore_X);
Eye_Explore_X_Line=Eye_Explore_X_Line(NoNNaNSelector,:);

Eye_Explore_Y=SelectEyeInterval(EyeChannel_Y_L,StartTimeBin,EndTimeBin)';
Eye_Explore_Y_Line=ReorganizeEye(Eye_Explore_Y);
Eye_Explore_Y_Line=Eye_Explore_Y_Line(NoNNaNSelector,:);
% Loop over the blocks

for i=1:NumberOfBlocks
    %Align eye trace with spikes for each block
    
      currBlocks=(BlockType_Line==UniqueBlockType(i));%& ismember(SceneName_Line,SceneName_Unique(1));%& ismember(TargetNum_Line,[1]);
      Eye_Explore_X_Line_currBlock=Eye_Explore_X_Line(currBlocks,:);
      Eye_Explore_Y_Line_currBlock=Eye_Explore_Y_Line(currBlocks,:);
      
      TargetLoc_Line_currBlock=TargetLoc_Line(currBlocks,:);
      TargetPolarAngle_Line_currBlock=TargetPolarAngle_Line(currBlocks,:);
      
      %Rotate the eye trace relative to the target location
      [RotatedEye_X,RotatedEye_Y]=RotateEyeTrace(Eye_Explore_X_Line_currBlock,Eye_Explore_Y_Line_currBlock,TargetPolarAngle_Line_currBlock);
        %Now each rotated eye change into a row with each trial seperated by
      %a NaN
      %Select out spike information
      
      PSTH_Event_SpikeCount_currBlock=[PSTH_Event_SpikeCount{1}(currBlocks,:) NaN*ones(size(Eye_Explore_X_Line_currBlock,1),size(Eye_Explore_X_Line_currBlock,2)-size(PSTH_Event_SpikeCount{1}(currBlocks,:),2)+1)];
      
      %Reorganize spike counts
      PSTH_Event_SpikeCount_currBlock_Reorg=reshape(PSTH_Event_SpikeCount_currBlock',1,numel(PSTH_Event_SpikeCount_currBlock));
      %Thresholding
      PSTH_Event_SpikeCount_currBlock_Reorg(PSTH_Event_SpikeCount_currBlock_Reorg>0)=1;%For plotting;
      PSTH_Event_SpikeCount_currBlock_Reorg(PSTH_Event_SpikeCount_currBlock_Reorg==0)=NaN;%For plotting
      
      RotatedEye_X_Plot_Block{i}=RotatedEye_X;
      RotatedEye_Y_Plot_Block{i}=RotatedEye_Y;
      SpikeEvent_currBlock{i}=PSTH_Event_SpikeCount_currBlock_Reorg;
      
      %Re-organize into a two-dimensional map
     % Eye_X_TwoD=reshape(RotatedEye_X,sum(currBlocks),ceil(length(RotatedEye_X)/sum(currBlocks)));
     % Eye_Y_TwoD=reshape(RotatedEye_Y,sum(currBlocks),ceil(length(RotatedEye_Y)/sum(currBlocks)));
      
     %Raw eye trace and spike
      Eye_X_ori=[Eye_Explore_X_Line_currBlock,NaN*ones(size(Eye_Explore_X_Line_currBlock,1),1)];
      Eye_X_ori_Line=reshape(Eye_X_ori',1,numel(Eye_X_ori));
      Eye_Y_ori=[Eye_Explore_Y_Line_currBlock,NaN*ones(size(Eye_Explore_Y_Line_currBlock,1),1)];
      Eye_Y_ori_Line=reshape(Eye_Y_ori',1,numel(Eye_Y_ori));
      RawEye_X_Plot_Block{i}=Eye_X_ori_Line;
      RawEye_Y_Plot_Block{i}=Eye_Y_ori_Line;
      
      TargetLoc_Raw_Block{i}=unique(TargetLoc_Line_currBlock,'rows');
     
     %2D Density map rotate
      data_tmp=TwoDimHist(RotatedEye_X,RotatedEye_Y,PSTH_Event_SpikeCount_currBlock_Reorg,Param);
      
      FR_Mat_Block_Rotate{i}=data_tmp.HistData;
 
      ScreenLoc_X_Block_Rotate{i}=data_tmp.HistX;
      ScreenLoc_Y_Block_Rotate{i}=data_tmp.HistY;
       
     %2D Density map raw
     data_tmp=TwoDimHist(Eye_X_ori_Line,Eye_Y_ori_Line,PSTH_Event_SpikeCount_currBlock_Reorg,Param);
      
      FR_Mat_Block_Raw{i}=data_tmp.HistData;
 
      ScreenLoc_X_Block_Raw{i}=data_tmp.HistX;
      ScreenLoc_Y_Block_Raw{i}=data_tmp.HistY;
      
      
    %Raw eye trace and spike train for each block and each scene
    
    for j=1:length(SceneName_Unique)
        Sel=(BlockType_Line==UniqueBlockType(i)) & (SceneName_Line==SceneName_Unique(j));%& ismember(TargetNum_Line,[4]);
        Eye_Explore_X_Line_currSceneBlock=Eye_Explore_X_Line(Sel,:);
        Eye_Explore_Y_Line_currSceneBlock=Eye_Explore_Y_Line(Sel,:);
        
        TargetLoc_Line_currSceneBlock{i,j}=unique(TargetLoc_Line(Sel,:),'rows');
        TargetPolarAngle_Line_currSceneBlock{i,j}=unique(TargetPolarAngle_Line(Sel,:),'rows');
         
        
        %Select out spike information  
      PSTH_Event_SpikeCount_currSceneBlock=[PSTH_Event_SpikeCount{1}(Sel,:) NaN*ones(size(Eye_Explore_X_Line_currSceneBlock,1),size(Eye_Explore_X_Line_currSceneBlock,2)-size(PSTH_Event_SpikeCount{1}(Sel,:),2)+1)];
      %Reorganize spike counts
      PSTH_Event_SpikeCount_currSceneBlock_Reorg=reshape(PSTH_Event_SpikeCount_currSceneBlock',1,numel(PSTH_Event_SpikeCount_currSceneBlock)); 
      %Thresholding
      PSTH_Event_SpikeCount_currSceneBlock_Reorg(PSTH_Event_SpikeCount_currSceneBlock_Reorg>0)=1;%For plotting;
      PSTH_Event_SpikeCount_currSceneBlock_Reorg(PSTH_Event_SpikeCount_currSceneBlock_Reorg==0)=NaN;%For plotting
      %Reorganize eyetrace
      Eye_X_tmp=[Eye_Explore_X_Line_currSceneBlock,NaN*ones(size(Eye_Explore_X_Line_currSceneBlock,1),1)];
      Eye_Y_tmp=[Eye_Explore_Y_Line_currSceneBlock,NaN*ones(size(Eye_Explore_Y_Line_currSceneBlock,1),1)];
      
      Eye_Explore_X_Line_currSceneBlock_Reorg=reshape(Eye_X_tmp',1,numel(Eye_X_tmp));
      Eye_Explore_Y_Line_currSceneBlock_Reorg=reshape(Eye_Y_tmp',1,numel(Eye_Y_tmp));
      
      Eye_X_Plot_SceneBlock{i,j}= Eye_Explore_X_Line_currSceneBlock_Reorg;
      Eye_Y_Plot_SceneBlock{i,j}=Eye_Explore_Y_Line_currSceneBlock_Reorg;
      SpikeEvent_SceneBlock{i,j}=PSTH_Event_SpikeCount_currSceneBlock_Reorg;
      
    end
   
    %Loop over the fractal orders
    TargetNumUnique=[1,2,3,min(4,max(TargetNum))];%Only separate into Order1, Order2, Oder3,>Order3
    for j=1:length(TargetNumUnique)
        if TargetNumUnique(j)~=TargetNumUnique(end)
            TargetNumSel=TargetNumUnique(j);
        else
            TargetNumSel=[TargetNumUnique(j):max(TargetNum)];
        end
        Sel=(BlockType_Line==UniqueBlockType(i))& ismember(TargetNum_Line,TargetNumSel);
        %Temporally not rotate
        Eye_Explore_X_Line_currTgtOrdBlock=Eye_Explore_X_Line(Sel,:);
        Eye_Explore_Y_Line_currTgtOrdBlock=Eye_Explore_Y_Line(Sel,:);
        
        TargetLoc_Line_currTgtOrdBlock{i,j}=unique(TargetLoc_Line(Sel,:),'rows');
        TargetPolarAngle_Line_currTgtOrdBlock{i,j}=unique(TargetPolarAngle_Line(Sel,:),'rows');
         
        
        %Select out spike information  
      PSTH_Event_SpikeCount_currTgtOrdBlock=[PSTH_Event_SpikeCount{1}(Sel,:) NaN*ones(size(Eye_Explore_X_Line_currTgtOrdBlock,1),size(Eye_Explore_X_Line_currTgtOrdBlock,2)-size(PSTH_Event_SpikeCount{1}(Sel,:),2)+1)];
      %Reorganize spike counts
      PSTH_Event_SpikeCount_currTgtOrdBlock_Reorg=reshape(PSTH_Event_SpikeCount_currTgtOrdBlock',1,numel(PSTH_Event_SpikeCount_currTgtOrdBlock)); 
      %Thresholding
      PSTH_Event_SpikeCount_currTgtOrdBlock_Reorg(PSTH_Event_SpikeCount_currTgtOrdBlock_Reorg>0)=1;%For plotting;
      PSTH_Event_SpikeCount_currTgtOrdBlock_Reorg(PSTH_Event_SpikeCount_currTgtOrdBlock_Reorg==0)=NaN;%For plotting
      %Reorganize eyetrace
      Eye_X_tmp=[Eye_Explore_X_Line_currTgtOrdBlock,NaN*ones(size(Eye_Explore_X_Line_currTgtOrdBlock,1),1)];
      Eye_Y_tmp=[Eye_Explore_Y_Line_currTgtOrdBlock,NaN*ones(size(Eye_Explore_Y_Line_currTgtOrdBlock,1),1)];
      
      Eye_Explore_X_Line_currTgtOrdBlock_Reorg=reshape(Eye_X_tmp',1,numel(Eye_X_tmp));
      Eye_Explore_Y_Line_currTgtOrdBlock_Reorg=reshape(Eye_Y_tmp',1,numel(Eye_Y_tmp));
      
      Eye_X_Plot_TgtOrdBlock{i,j}= Eye_Explore_X_Line_currTgtOrdBlock_Reorg;
      Eye_Y_Plot_TgtOrdBlock{i,j}=Eye_Explore_Y_Line_currTgtOrdBlock_Reorg;
      SpikeEvent_TgtOrdBlock{i,j}=PSTH_Event_SpikeCount_currTgtOrdBlock_Reorg;
      
    end
   
  
end
%{
figure
for i=1:NumberOfBlocks
subplot(ceil(NumberOfBlocks/2),2,i)
       h=heatmap(FR_Mat{i});%-nanmean(FR_Mat{i},'all'));
       h.GridVisible='off';
      % xlim([0,60/BinDim(1)]);
      %ylim([0,60/BinDim(2)]);
      colormap default
%caxis([0,max(FR_Max)]);
end
%%
%}
%%
if ShowFigureFlag
    FigureStartNum=150;
    FigureIndex=1;
   if PlotSpikeEventOn
    %First figure, overlay of all the scenes for all fractals
    figtitlestr{FigureIndex}='ExploreSpikeEventPlot_Rotate';
    fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[50,100, 800,800],'Name',figtitlestr{FigureIndex});
    for i=1:NumberOfBlocks
       
        BlockOrder=find(AssociateBlock==UniqueBlockType(i));
        subplot(max(ceil(NumberOfBlocks/2),1),2,i);
        %Plot scene area
         r=25;
         ColorScene=[0.9,0.9,0.9]; 
         cplot = @(r,x0,y0) fill(r*cos(linspace(0,2*pi))+x0,r*sin(linspace(0,2*pi))+y0,[0.9,0.9,0.9]); 
        % x0 + r*cos(linspace(0,2*pi)),y0 + r*sin(linspace(0,2*pi)),'-');
         c=cplot(r,0,0);
         c.FaceColor=ColorScene;
         c.FaceAlpha=0.5;
         c.EdgeColor=ColorScene;
         c.EdgeAlpha=0.5;
        hold on
        %Overlay target location
   
          r=7.5;
           t=cplot(r,0,15);
           t.FaceColor=ColorFractal;
           t.FaceAlpha=0.1;
           t.EdgeColor=ColorFractal;
           t.EdgeAlpha=0.1;
            
        %First eye trace
         p=plot(RotatedEye_X_Plot_Block{i},RotatedEye_Y_Plot_Block{i},'-r');
  
         p.LineStyle='-';
         p.Color=ColorForEachMode(find(AssociateBlock==UniqueBlockType(i)),:);
         p.Color(4)=0.15;
         p.LineWidth=0.5;
         p.MarkerSize=0.5;
      hold on
      %Overlay spike
     % plot(RotatedEye_X_Plot_Block{i}(SpikeEvent_currBlock{i}==1),RotatedEye_Y_Plot_Block{i}(SpikeEvent_currBlock{i}==1),'or','color',[11,60,238]/255);
      s=scatter(RotatedEye_X_Plot_Block{i},RotatedEye_Y_Plot_Block{i},SpikeSizePlot*SpikeEvent_currBlock{i},[11,60,238]/255,'filled');
      xlim([-30,30]);
      ylim([-30,30]);
       %axis equal;
      axis off;
      title(Blockstr{ BlockOrder},'FontSize',15,'FontWeight','Bold');
      
  

    end
    FigureStartNum=FigureStartNum+1;
    FigureIndex=FigureIndex+1;
    
%Figure2, Raw eye trace (without roation) and spike events for each block
figtitlestr{FigureIndex}='ExploreSpikeEventPlot_Raw';
    fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[80,100, 800,800],'Name',figtitlestr{FigureIndex});
    for i=1:NumberOfBlocks
       
        BlockOrder=find(AssociateBlock==UniqueBlockType(i));
        subplot(max(ceil(NumberOfBlocks/2),1),2,i);
        %Plot scene area
         r=25;
         ColorScene=[0.9,0.9,0.9]; 
         cplot = @(r,x0,y0) fill(r*cos(linspace(0,2*pi))+x0,r*sin(linspace(0,2*pi))+y0,[0.9,0.9,0.9]); 
        % x0 + r*cos(linspace(0,2*pi)),y0 + r*sin(linspace(0,2*pi)),'-');
         c=cplot(r,0,0);
         c.FaceColor=ColorScene;
         c.FaceAlpha=0.5;
         c.EdgeColor=ColorScene;
         c.EdgeAlpha=0.5;
        hold on
        %Overlay target location
        currTgtLoc=TargetLoc_Raw_Block{i};
        for jj=1:size(currTgtLoc,1)
   
          r=7.5;
           t=cplot(r,currTgtLoc(jj,1),currTgtLoc(jj,2));
           t.FaceColor=ColorFractal;
           t.FaceAlpha=0.1;
           t.EdgeColor=ColorFractal;
           t.EdgeAlpha=0.1;
        end
            
        %First eye trace
         p=plot(RawEye_X_Plot_Block{i},RawEye_Y_Plot_Block{i},'-r');
  
         p.LineStyle='-';
         p.Color=ColorForEachMode(find(AssociateBlock==UniqueBlockType(i)),:);
         p.Color(4)=0.15;
         p.LineWidth=0.5;
         p.MarkerSize=0.5;
      hold on
      %Overlay spike
     % plot(RotatedEye_X_Plot_Block{i}(SpikeEvent_currBlock{i}==1),RotatedEye_Y_Plot_Block{i}(SpikeEvent_currBlock{i}==1),'or','color',[11,60,238]/255);
      s=scatter(RawEye_X_Plot_Block{i},RawEye_Y_Plot_Block{i},SpikeSizePlot*SpikeEvent_currBlock{i},[11,60,238]/255,'filled');
      xlim([-30,30]);
      ylim([-30,30]);
       %axis equal;
      axis off;
      title(Blockstr{ BlockOrder},'FontSize',15,'FontWeight','Bold');
      
  

    end
    FigureStartNum=FigureStartNum+1;
    FigureIndex=FigureIndex+1;
    
    
% Figure3, Separate for each scene for each mode

    for i=1:NumberOfBlocks
       BlockOrder=find(AssociateBlock==UniqueBlockType(i));
       figtitlestr{FigureIndex}=sprintf('Each Scene,%s',Blockstr{BlockOrder});
       fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[70+(i-1)*400,100, 400,700],'Name',figtitlestr{FigureIndex});
        ValidIndex=~cellfun(@isempty,Eye_X_Plot_SceneBlock(i,:));
        SceneName_Unique_currBlock=SceneName_Unique(ValidIndex');
        
        Eye_X_Plot_SceneBlock_currBlock=Eye_X_Plot_SceneBlock(i,ValidIndex);
        Eye_Y_Plot_SceneBlock_currBlock=Eye_Y_Plot_SceneBlock(i,ValidIndex);
        SpikeEvent_SceneBlock_currBlock=SpikeEvent_SceneBlock(i,ValidIndex);
        for j=1:length(SceneName_Unique_currBlock)
          %  if ~isempty(Eye_X_Plot_SceneBlock{i,j})
          subplot(max(ceil(length(SceneName_Unique_currBlock)/2),1),2,j);
        %Plot scene area
         r=25;
         ColorScene=[0.9,0.9,0.9]; 
         cplot = @(r,x0,y0) fill(r*cos(linspace(0,2*pi))+x0,r*sin(linspace(0,2*pi))+y0,[0.9,0.9,0.9]); 
        % x0 + r*cos(linspace(0,2*pi)),y0 + r*sin(linspace(0,2*pi)),'-');
         c=cplot(r,0,0);
         c.FaceColor=ColorScene;
         c.FaceAlpha=0.5;
         c.EdgeColor=ColorScene;
         c.EdgeAlpha=0.5;
        hold on
        %Overlay target location
        currTgtLoc= TargetLoc_Line_currSceneBlock{i,j};
        for jj=1:size(currTgtLoc,1)
   
          r=7.5;
           t=cplot(r,currTgtLoc(jj,1),currTgtLoc(jj,2));
           t.FaceColor=ColorFractal;
           t.FaceAlpha=0.1;
           t.EdgeColor=ColorFractal;
           t.EdgeAlpha=0.1;
        end
            
        %First eye trace
         p=plot(Eye_X_Plot_SceneBlock_currBlock{j},Eye_Y_Plot_SceneBlock_currBlock{j},'-r');
  
         p.LineStyle='-';
         p.Color=ColorForEachMode(find(AssociateBlock==UniqueBlockType(i)),:);
         p.Color(4)=0.15;
         p.LineWidth=0.5;
         p.MarkerSize=0.5;
         hold on
      %Overlay spike
     % plot(RotatedEye_X_Plot_Block{i}(SpikeEvent_currBlock{i}==1),RotatedEye_Y_Plot_Block{i}(SpikeEvent_currBlock{i}==1),'or','color',[11,60,238]/255);
      s=scatter(Eye_X_Plot_SceneBlock_currBlock{j},Eye_Y_Plot_SceneBlock_currBlock{j},SpikeSizePlot*SpikeEvent_SceneBlock_currBlock{j},[11,60,238]/255,'filled');
      xlim([-30,30]);
      ylim([-30,30]);
     % axis equal;
      axis off;
      currScene=SceneName_Unique_currBlock(j);
      title(sprintf('Scene%d',currScene),'FontSize',15,'FontWeight','Bold');   
            
          %  end
            
        end %End of loop for each scene
         FigureStartNum=FigureStartNum+1;
         FigureIndex=FigureIndex+1;
        
    end %End of loop for each block
    
% Figure3, Separate for each fractal order for each mode

    for i=1:NumberOfBlocks
       BlockOrder=find(AssociateBlock==UniqueBlockType(i));
       figtitlestr{FigureIndex}=sprintf('Fractal_Order,%s',Blockstr{BlockOrder});
       fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[90+(i-1)*400,100, 400,700],'Name',figtitlestr{FigureIndex});

        for j=1:length(TargetNumUnique)
          subplot(max(ceil(length(TargetNumUnique)/2),1),2,j);
        %Plot scene area
         r=25;
         ColorScene=[0.9,0.9,0.9]; 
         cplot = @(r,x0,y0) fill(r*cos(linspace(0,2*pi))+x0,r*sin(linspace(0,2*pi))+y0,[0.9,0.9,0.9]); 
        % x0 + r*cos(linspace(0,2*pi)),y0 + r*sin(linspace(0,2*pi)),'-');
         c=cplot(r,0,0);
         c.FaceColor=ColorScene;
         c.FaceAlpha=0.5;
         c.EdgeColor=ColorScene;
         c.EdgeAlpha=0.5;
        hold on
        %Overlay target location
        currTgtLoc= TargetLoc_Line_currTgtOrdBlock{i,j};
        for jj=1:size(currTgtLoc,1)
   
          r=7.5;
           t=cplot(r,currTgtLoc(jj,1),currTgtLoc(jj,2));
           t.FaceColor=ColorFractal;
           t.FaceAlpha=0.1;
           t.EdgeColor=ColorFractal;
           t.EdgeAlpha=0.1;
        end
            
        %First eye trace
         p=plot(Eye_X_Plot_TgtOrdBlock{i,j},Eye_Y_Plot_TgtOrdBlock{i,j},'-r');
  
         p.LineStyle='-';
         p.Color=ColorForEachMode(find(AssociateBlock==UniqueBlockType(i)),:);
         p.Color(4)=0.15;
         p.LineWidth=0.5;
         p.MarkerSize=0.5;
         hold on
      %Overlay spike
     % plot(RotatedEye_X_Plot_Block{i}(SpikeEvent_currBlock{i}==1),RotatedEye_Y_Plot_Block{i}(SpikeEvent_currBlock{i}==1),'or','color',[11,60,238]/255);
      s=scatter(Eye_X_Plot_TgtOrdBlock{i,j},Eye_Y_Plot_TgtOrdBlock{i,j},SpikeSizePlot*SpikeEvent_TgtOrdBlock{i,j},[11,60,238]/255,'filled');
      xlim([-30,30]);
      ylim([-30,30]);
     % axis equal;
      axis off;
      currTgtOrder=TargetNumUnique(j);
      if j~=length(TargetNumUnique)
         title(sprintf('Order %d',currTgtOrder),'FontSize',15,'FontWeight','Bold');  
      else
         title(sprintf('Order> %d',currTgtOrder-1),'FontSize',15,'FontWeight','Bold');  
      end
            
            
            
        end %End of loop for each scene
         FigureStartNum=FigureStartNum+1;
         FigureIndex=FigureIndex+1;
        
    end %End of loop for each block
   end %End of the if PlotSpikeEventOn
   
if PlotSpikeDensityOn   
%Figure5, 2D SpikeDenisty with Rotated eye trace for each block

figtitlestr{FigureIndex}='TwoDSpikeDensityPlot_Rotate';
fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[100,100, 800,800],'Name',figtitlestr{FigureIndex});
    for i=1:NumberOfBlocks
       
        BlockOrder=find(AssociateBlock==UniqueBlockType(i));
        subplot(max(ceil(NumberOfBlocks/2),1),2,i);
          contourf(ScreenLoc_X_Block_Rotate{i}(1:end-1)+diff(ScreenLoc_X_Block_Rotate{i})/2,ScreenLoc_Y_Block_Rotate{i}(1:end-1)+diff(ScreenLoc_Y_Block_Rotate{i})/2,FR_Mat_Block_Rotate{i});

      %   h=heatmap(FR_Mat_Block_Rotate{i},'title',Blockstr{ BlockOrder});%-nanmean(FR_Mat{i},'all'));
     %  h.GridVisible='off';
       
      % xlim([0,60/BinDim(1)]);
      %ylim([0,60/BinDim(2)]);
      colormap default;
     % h.XDisplayLabels=round(ScreenLoc_X_Block_Rotate{i}(1:end-1)+diff(ScreenLoc_X_Block_Rotate{i})/2);
     % h.YDisplayLabels=round(ScreenLoc_Y_Block_Rotate{i}(1:end-1)+diff(ScreenLoc_Y_Block_Rotate{i})/2);
 

    end
    FigureStartNum=FigureStartNum+1;
    FigureIndex=FigureIndex+1;  
    
    
    
  %  FR_Mat_Block_Raw{i}=data_tmp.HistData;
 
   %   ScreenLoc_X_Block_Raw{i}=data_tmp.HistX;
   %   ScreenLoc_Y_Block_Raw{i}=data_tmp.HistY;
   
   
  %Figure6, 2D SpikeDenisty with Raw eye trace for each block

figtitlestr{FigureIndex}='TwoDSpikeDensityPlot_Raw';
fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[120,100, 800,800],'Name',figtitlestr{FigureIndex});
    for i=1:NumberOfBlocks
       
        BlockOrder=find(AssociateBlock==UniqueBlockType(i));
        subplot(max(ceil(NumberOfBlocks/2),1),2,i);
        %{
         h=heatmap(FR_Mat_Block_Raw{i},'title',Blockstr{ BlockOrder});%-nanmean(FR_Mat{i},'all'));
       h.GridVisible='off';
       
      % xlim([0,60/BinDim(1)]);
      %ylim([0,60/BinDim(2)]);
      colormap default;
      h.XDisplayLabels=round(ScreenLoc_X_Block_Raw{i}(1:end-1)+diff(ScreenLoc_X_Block_Raw{i})/2);
      h.YDisplayLabels=round(ScreenLoc_Y_Block_Raw{i}(1:end-1)+diff(ScreenLoc_Y_Block_Raw{i})/2);
 %}
        contourf(ScreenLoc_X_Block_Raw{i}(1:end-1)+diff(ScreenLoc_X_Block_Raw{i})/2,ScreenLoc_Y_Block_Raw{i}(1:end-1)+diff(ScreenLoc_Y_Block_Raw{i})/2,FR_Mat_Block_Raw{i});

    end
    FigureStartNum=FigureStartNum+1;
    FigureIndex=FigureIndex+1;  
    
   
    
  %  FR_Mat_Block_Raw{i}=data_tmp.HistData;
 
   %   ScreenLoc_X_Block_Raw{i}=data_tmp.HistX;
   %   ScreenLoc_Y_Block_Raw{i}=data_tmp.HistY; 
   
    
end %End of the if PlotSpikeDensityOn   
    
end %End of the if ShowFigureFlag







   



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