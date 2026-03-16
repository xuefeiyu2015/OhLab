function SceneOnExploreNeuralAnalysis(Data, StartTrial, EndTrial,ShowFigureFlag,OutputFlag);
%Function to analysis the explore behavior during the scene onset period;
%Xuefei Yu 08282020
ECode;
TaskType=Data.TaskType;
TrialNum=EndTrial-StartTrial+1;
Selection=StartTrial:EndTrial;%Select the trials according to the pannel
%Abnormal check %Especially For Robin081120FORAGEXPLO4.rst
CheckCode=SCENE_ON;%This code should appear once in each trial
AbnormalTrials=AbnormalTrialCheck(Data,CheckCode);
Selection=Selection(~ismember(Selection,AbnormalTrials));

%User defined parameter for the neural/behavioral analysis
FPOnOffset=0;%in ms, the time selection bin after fixation onset; 
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

SegEvent=[SCENE_ON];%Event Flow of interest

SegEventEnd=[SHOWFIXCD];%Event end
SegEventStr={'SceneOn'};
PlotStartOffset=[StartDelay];%in ms
PlotStopOffset=[FPOnOffset];%in ms

EventTypeFlag=[1];%1:Number same with the trial number, 2:Number same with the target number

%Parameter to define the sliding window to measure firing rates

BinWidth=15;%in ms

%Parameter for raster plot
Barlength=1;
LineWidth=1;
%Parameter for visualize spike event
SpikeSizePlot=10;
%Parameter for 2D histogram for spikes
%2D Density map
TwoDBinRange=[-25,25;-25,25];%Width;Height
BinDim=[1.5;1.5];%Width;Height
Param.Range=TwoDBinRange;
Param.Bin=BinDim;
%Param.FR_Threshold=100;%ms/deg2
Param.FR_Threshold=50;%ms/degree
Param.GaussianSmooth=1;
Param.GaussianSigma=4;
Param.GaussianWindow=[13;13];
Param.OriBinWidth=EyeBinWidth;


%ResponseLatencyThreshold=3;

%Some user defined parameter for visualization
%Color Definition for each scenes
Color={'r','g','b','m',[0.4,0.4,0],[0.28,0.78,0.89],'k',[0.6,0.1,1],[0.6,0.1,0.4],[0.3,0.5,1],[0.6,0.5,0.1],[0,0.8,0.4]};
ColorForEachMode=[0.83 0.14 0.14;
             1.00 0.54 0.00;
             0.47 0.25 0.80;
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
    
FP_Loc=Data.FixLoc(Selection,2:end,:);
Ref_Loc=squeeze(FP_Loc(:,1,:));
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



 %Block type
[UniqueBlockType,orderindex,reprojectorderindex]=unique(BlockType,'row');
NumberOfBlocks=max(reprojectorderindex);
 

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

%Rearrange Scene type into lines
SceneName_Line=reshape(SceneName,numel(SceneName),1);
NoNNaNSelector=~isnan(SceneName_Line);
SceneName_Line=SceneName_Line(NoNNaNSelector);
SceneName_Unique=unique(SceneName_Line);



%Rearrange block type into lines
BlockType_Line=repmat(BlockType,1,size(SceneName,2));
BlockType_Line=reshape(BlockType_Line,numel(BlockType_Line),1);
BlockType_Line=BlockType_Line(NoNNaNSelector);




%Select out eye trace for each target
SegEventTime=FindOutTime(EventChannel,EventTimeChannel,SegEvent(1));%For this protocol, there is only one event
SegEventEndTime=FindOutTime(EventChannel,EventTimeChannel,SegEventEnd(1));%For this protocol, there is only one event


StartTimeBin=floor((SegEventTime-PlotStartOffset(1))/EyeBinWidth);         
EndTimeBin=ceil((SegEventEndTime+PlotStopOffset(1))/EyeBinWidth);
         
Eye_Explore_X=SelectEyeInterval(EyeChannel_X_L,StartTimeBin,EndTimeBin)';


Eye_Explore_Y=SelectEyeInterval(EyeChannel_Y_L,StartTimeBin,EndTimeBin)';


Eye_Explore_X_Line=ReorganizeEye(Eye_Explore_X);
Eye_Explore_X_Line=Eye_Explore_X_Line(NoNNaNSelector,:);

Eye_Explore_Y_Line=ReorganizeEye(Eye_Explore_Y);
Eye_Explore_Y_Line=Eye_Explore_Y_Line(NoNNaNSelector,:);


%Period For EyeMovemen Analysis (to remove initial visual onset responses)
% Loop over the blocks

for i=1:NumberOfBlocks
    %Align eye trace with spikes for each block
    
      currBlocks=(BlockType_Line==UniqueBlockType(i));%& ismember(SceneName_Line,SceneName_Unique(1));%& ismember(TargetNum_Line,[1]);
      Eye_Explore_X_Line_currBlock=Eye_Explore_X_Line(currBlocks,:);
      Eye_Explore_Y_Line_currBlock=Eye_Explore_Y_Line(currBlocks,:);
    
      %Select out spike information
      
      PSTH_Event_SpikeCount_currBlock=[PSTH_Event_SpikeCount{1}(currBlocks,:) NaN*ones(size(Eye_Explore_X_Line_currBlock,1),size(Eye_Explore_X_Line_currBlock,2)-size(PSTH_Event_SpikeCount{1}(currBlocks,:),2)+1)];
      
      %Reorganize spike counts
      PSTH_Event_SpikeCount_currBlock_Reorg=reshape(PSTH_Event_SpikeCount_currBlock',1,numel(PSTH_Event_SpikeCount_currBlock));
      %Thresholding
      PSTH_Event_SpikeCount_currBlock_Reorg(PSTH_Event_SpikeCount_currBlock_Reorg>0)=1;%For plotting;
      PSTH_Event_SpikeCount_currBlock_Reorg(PSTH_Event_SpikeCount_currBlock_Reorg==0)=NaN;%For plotting
      
      
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
      
     
     
    
     %2D Density map raw
     data_tmp=TwoDimHist(Eye_X_ori_Line,Eye_Y_ori_Line,PSTH_Event_SpikeCount_currBlock_Reorg,Param);
      
      FR_Mat_Block_Raw{i}=data_tmp.HistData;
 
      ScreenLoc_X_Block_Raw{i}=data_tmp.HistX;
      ScreenLoc_Y_Block_Raw{i}=data_tmp.HistY;
      
     %{ 
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
   
   
   %}
  
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
   
    
%Figure1, Raw eye trace (without roation) and spike events for each block
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
    
 %{   
% Figure3, Separate for each scene for each mode

    for i=1:NumberOfBlocks
       BlockOrder=find(AssociateBlock==UniqueBlockType(i));
       figtitlestr{FigureIndex}=sprintf('Each Scene,%s',Blockstr{BlockOrder});
       fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[70+(i-1)*400,100, 400,700],'Name',figtitlestr{FigureIndex});

        for j=1:length(SceneName_Unique)
          subplot(max(ceil(length(SceneName_Unique)/2),1),2,j);
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
         p=plot(Eye_X_Plot_SceneBlock{i,j},Eye_Y_Plot_SceneBlock{i,j},'-r');
  
         p.LineStyle='-';
         p.Color=ColorForEachMode(find(AssociateBlock==UniqueBlockType(i)),:);
         p.Color(4)=0.15;
         p.LineWidth=0.5;
         p.MarkerSize=0.5;
         hold on
      %Overlay spike
     % plot(RotatedEye_X_Plot_Block{i}(SpikeEvent_currBlock{i}==1),RotatedEye_Y_Plot_Block{i}(SpikeEvent_currBlock{i}==1),'or','color',[11,60,238]/255);
      s=scatter(Eye_X_Plot_SceneBlock{i,j},Eye_Y_Plot_SceneBlock{i,j},SpikeSizePlot*SpikeEvent_SceneBlock{i,j},[11,60,238]/255,'filled');
      xlim([-30,30]);
      ylim([-30,30]);
     % axis equal;
      axis off;
      currScene=SceneName_Unique(j);
      title(sprintf('Scene%d',currScene),'FontSize',15,'FontWeight','Bold');   
            
            
            
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
    %}
   end %End of the if PlotSpikeEventOn
   
if PlotSpikeDensityOn   

   
   
  %Figure5, 2D SpikeDenisty with Raw eye trace for each block

figtitlestr{FigureIndex}='TwoDSpikeDensityPlot_Raw';
fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[120,100, 800,800],'Name',figtitlestr{FigureIndex});
    for i=1:NumberOfBlocks
       
        BlockOrder=find(AssociateBlock==UniqueBlockType(i));
        subplot(max(ceil(NumberOfBlocks/2),1),2,i);
        
         h=heatmap(FR_Mat_Block_Raw{i},'title',Blockstr{ BlockOrder});%-nanmean(FR_Mat{i},'all'));
       h.GridVisible='off';
       
      % xlim([0,60/BinDim(1)]);
      %ylim([0,60/BinDim(2)]);
      colormap default;
      h.XDisplayLabels=round(ScreenLoc_X_Block_Raw{i}(1:end-1)+diff(ScreenLoc_X_Block_Raw{i})/2);
      h.YDisplayLabels=round(ScreenLoc_Y_Block_Raw{i}(1:end-1)+diff(ScreenLoc_Y_Block_Raw{i})/2);
 

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
