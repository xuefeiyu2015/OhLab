function  ExploreRasterCheck(Data, StartTrial, EndTrial,ShowFigureFlag,OutputFlag);
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
TargetOnOffset=50;%in ms, the time selection bin after target onset; 
StartDelay=0;%in ms, the time selection bin after fixation offset;

%For segmentation of  the neural responses 
SegEvent=[EXTRAFIXSTARTCD,SCENE_ON,SHOWFIXCD,FIXOFF,TGTCD,RWDCD  ];%Event Flow of interest
SegEventEnd=[EXTRAFIXOFF,SHOWFIXCD,FIXOFF,TGTCD,TGTOFF,];%Event end; %Automatic screening later
%SegEventEnd=[EXTRAFIXOFF,FIXOFF,FIXOFF,TGTCD,TGTOFF];%Event end
SegEventStr={'ExtraFPOn','SceneOn','FixOn','Explore','ChoicePeriod'};

PlotStartOffset=[200,200,200,100,500];%in ms
PlotStopOffset=[200,200,200,50,300];%in ms

EventTypeFlag=[1,1,2,2,2,2];%1:Number same with the trial number, 2:Number same with the target number

RasterOn=0;%Show raster, but it will be very slow to do so

%Parameter to define the sliding window to measure firing rates
StepSize=10;%in ms
BinWidth=15;%in ms

%Parameter for raster plot
Barlength=1;
LineWidth=1;


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

%Load event data
EventChannel=Data.EventChannel(Selection,:);
EventTimeChannel=Data.EventTimeChannel(Selection,:);
EventBin=Data.EventBin(Selection,:);

%Automatically chanage the marker

StartCheck=FindOutTime(EventChannel,EventBin,SegEvent(SegEvent==SCENE_ON));
EndCheck=FindOutTime(EventChannel,EventBin,SegEventEnd(SegEvent==SCENE_ON));
EndCheck=EndCheck(:,1);
DiffCheck=EndCheck-StartCheck;
if sum(DiffCheck<200)>ceil(length(StartCheck)/2)
    SegEventEnd(SegEvent==SCENE_ON)=[FIXOFF];%Event end;change from showfix cd to fix off for contigency
    disp('Change to FixOff for FPOn First Mode');
    
end

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

%Rearrange block type into lines
BlockType_Line=repmat(BlockType,1,max(TargetNum));
BlockType_Line=reshape(BlockType_Line,numel(BlockType_Line),1);
BlockType_Line=BlockType_Line(NoNNaNSelector);


% Loop over the blocks
for i=1:NumberOfBlocks
    
  
    %Organize the raster, spikecount PSTH, spike density function PSTH
    for j=1:length(AlignedRaster)
        if EventTypeFlag(j)==1 %Trial event
            currBlocks=BlockType==UniqueBlockType(i);
            
        elseif EventTypeFlag(j)==2 %Target event
            currBlocks=BlockType_Line==UniqueBlockType(i);
        end
        
        %Raster
        Raster_currEvent_Block{i,j}=AlignedRaster{j}(currBlocks,:);
        
        
        %Spike count PSTH
        PSTH_Event_SpikeCount_currEvent{i,j}=PSTH_Event_SpikeCount{j}(currBlocks,:);
        %Spike density function PSTH
        MeanPSTH_FR_currEvent{i,j}=MeanPSTH_FR{j}(currBlocks,:);
        
        %Timing 
        PSTH_Time_currEvent{i,j}=nanmean(MeanPSTH_Time{j}(currBlocks,:),1);
        RelativeTimeMarker_currEvent{i,j}=RelativeTimeMarker{j}(currBlocks,:);
        
        %{
        StartLoc=1;
        Barlength=1;
        LineWidth=1;
        SpikeRaterPlot(Raster_currEvent_Block,StartLoc,Barlength,[1,0,0],LineWidth);
        %}
    end
    
 
end

%%
FigureStartNum=110;
FigureIndex=1;
if ShowFigureFlag
   if RasterOn
    figtitlestr{FigureIndex}='Event_Raster_Check';
    fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[50,100, 2000,800],'Name',figtitlestr{FigureIndex});
   
    for i=1:length(AlignedRaster)
        
         StartLoc=1;
        for j=1:NumberOfBlocks
            %First Figure: raster
            
            subplot(1,length(AlignedRaster),i);
            NumTrials=size(Raster_currEvent_Block{j,i},1);
            currBlock=UniqueBlockType(j);
            
            %Sort Event according to the event persist time
            [Val_Sort,SortIndex]=sort(RelativeTimeMarker_currEvent{j,i}(:,2));
            
            currResponse=Raster_currEvent_Block{j,i}(SortIndex,:);
            %Event-on line
            plot([0,0],[StartLoc-1,StartLoc+NumTrials],'--k','Color','#929498','LineWidth',1);
            hold on
            %Spike raster
            SpikeRaterPlot(currResponse,StartLoc,Barlength,ColorForEachMode(AssociateBlock==currBlock,:),LineWidth);
            
            %Add event-end marker
            y_tmp=[StartLoc:StartLoc+NumTrials-1]';
            x_tmp=Val_Sort;
            plot(x_tmp,y_tmp,'>','MarkerEdgeColor','#707276','MarkerFaceColor','#707276','MarkerSize',2);
            
            xlabel('Time from event onset');
            ylabel('Trial #');
            title(SegEventStr{i});
            set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
          %  legend(Blockstr{AssociateBlock==currBlock});
        
            StartLoc=StartLoc+NumTrials;
            %{
            %Second row: Spike-Count PSTH
            subplot(3,length(AlignedRaster),i+length(AlignedRaster));
            %Event-on line
            plot([0,0],[0,max(nanmean(PSTH_Event_SpikeCount_currEvent{i,j})/StepSize*1000)],'--k','Color','#929498','LineWidth',1);
            bar(PSTH_Time_currEvent{i,j},nanmean(PSTH_Event_SpikeCount_currEvent{i,j},1)/StepSize*1000,'FaceColor',ColorForEachMode(AssociateBlock==currBlock,:));
            hold on 
            %Add event-end marker
            plot([nanmean(RelativeTimeMarker_currEvent{j,i}(:,2)), nanmean(RelativeTimeMarker_currEvent{j,i}(:,2))],[0,max(nanmean(PSTH_Event_SpikeCount_currEvent{i,j})/StepSize*1000)],'-k','Color','#707276','LineWidth',3);
             ylabel('Spike Count per s(Hz)');
            set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
            %}
        end
    end
    FigureIndex=FigureIndex+1;
    FigureStartNum=FigureStartNum+1;
  end %End of the if raster on
  
  
 %%%%%%%%%%%%%%%%%%Next figure: Spike Count PSTH 
 
 figtitlestr{FigureIndex}='SpikeCountPSTH';
 fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[80,100, 600,800],'Name',figtitlestr{FigureIndex});   
    
   for i=1:length(AlignedRaster)
        
         
        for j=1:NumberOfBlocks
            
            %Spike-Count PSTH
           
            currBlock=UniqueBlockType(j);
            BlockOrder(j)=find(AssociateBlock==currBlock);
            
            subplot(length(AlignedRaster),1,i);
            %Event-on line
            curr_Response=nanmean(PSTH_Event_SpikeCount_currEvent{j,i},1)/StepSize*1000;
            curr_Response_err=nanstd(PSTH_Event_SpikeCount_currEvent{j,i}/StepSize*1000,[],1)./sqrt(sum(~isnan(PSTH_Event_SpikeCount_currEvent{j,i}),1));
            
            
            h1=plot([0,0],[0,max(curr_Response)],'--k','Color','#929498','LineWidth',1);
            hold on
           % plot(PSTH_Time_currEvent{j,i},curr_Response,'-k','Color',ColorForEachMode(AssociateBlock==currBlock,:),'LineWidth',2);
           % bar(PSTH_Time_currEvent{j,i},nanmean(PSTH_Event_SpikeCount_currEvent{j,i},1),'FaceColor',ColorForEachMode(AssociateBlock==currBlock,:));
             shadedErrorBar(PSTH_Time_currEvent{j,i},curr_Response,curr_Response_err,'lineprops',{ColorForEachMode(AssociateBlock==currBlock,:)},'transparent',1,'patchSaturation',0.3);%PSTH
              set(findall(gca, 'Type', 'Line'),'LineWidth',3);
            hold on 
            %Add event-end marker
            h2=plot([min(RelativeTimeMarker_currEvent{j,i}(:,2)), min(RelativeTimeMarker_currEvent{j,i}(:,2))],[0,max(curr_Response)],'-k','Color','#707276','LineWidth',1.5);
            xlim([-PlotStartOffset(i),PlotStopOffset(i)+nanmean(RelativeTimeMarker_currEvent{j,i}(:,2))]);
            
            %Set ylim
           %Set ylim
            Interval=curr_Response(:,PSTH_Time_currEvent{j,i}<PlotStopOffset(i)+nanmean(RelativeTimeMarker_currEvent{j,i}(:,2)) & PlotStopOffset(i)+nanmean(RelativeTimeMarker_currEvent{j,i}(:,2))>-PlotStartOffset(i));
            Interval_max(j)=max( Interval)*1.12;
            
           %Some legend off
            set( get( get( h1, 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );
            set( get( get( h2, 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );
            if i==3
             ylabel('Spike Count per s(Hz)');
            end
            if i==length(AlignedRaster)
                xlabel('Time from event onset');
            end
             title(SegEventStr{i});
            set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
            %}
        end
        ylim([0,max( Interval_max)]);
   end
     legend(Blockstr{BlockOrder});
     
    FigureIndex=FigureIndex+1;
    FigureStartNum=FigureStartNum+1; 

%%%%%%%%%%%%%%%%%%Next figure: Spike Density Function PSTH
 
 figtitlestr{FigureIndex}='SpikeDensityFunctionPSTH';
 fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[100,100, 600,800],'Name',figtitlestr{FigureIndex});   
    
   for i=1:length(AlignedRaster)
        
         
        for j=1:NumberOfBlocks
            
            %Spike-Count PSTH
           
            currBlock=UniqueBlockType(j);
            
            subplot(length(AlignedRaster),1,i);
            %Event-on line
            curr_Response=nanmean(MeanPSTH_FR_currEvent{j,i},1);
            curr_Response_err=nanstd(MeanPSTH_FR_currEvent{j,i},[],1)./sqrt(sum(~isnan(MeanPSTH_FR_currEvent{j,i}),1));
            
            plot([0,0],[0,max(curr_Response)],'--k','Color','#929498','LineWidth',1);
            hold on
            %plot(PSTH_Time_currEvent{j,i},curr_Response,'-k','Color',ColorForEachMode(AssociateBlock==currBlock,:),'LineWidth',2);
           % bar(PSTH_Time_currEvent{j,i},nanmean(PSTH_Event_SpikeCount_currEvent{j,i},1),'FaceColor',ColorForEachMode(AssociateBlock==currBlock,:));
            shadedErrorBar(PSTH_Time_currEvent{j,i},curr_Response,curr_Response_err,'lineprops',{ColorForEachMode(AssociateBlock==currBlock,:)},'transparent',1,'patchSaturation',0.5);%PSTH for fix
             set(findall(gca, 'Type', 'Line'),'LineWidth',3);
            hold on 
            %Add event-end marker
            plot([min(RelativeTimeMarker_currEvent{j,i}(:,2)), min(RelativeTimeMarker_currEvent{j,i}(:,2))],[0,max(curr_Response)],'-k','Color','#707276','LineWidth',1.5);
            xlim([-PlotStartOffset(i),PlotStopOffset(i)+nanmean(RelativeTimeMarker_currEvent{j,i}(:,2))]);
            
            %Set ylim
            %Set ylim
            Interval=curr_Response(:,PSTH_Time_currEvent{j,i}<PlotStopOffset(i)+nanmean(RelativeTimeMarker_currEvent{j,i}(:,2)) & PlotStopOffset(i)+nanmean(RelativeTimeMarker_currEvent{j,i}(:,2))>-PlotStartOffset(i));
              Interval_max(j)=max( Interval)*1.12;
            
           
            if i==3
             ylabel('Spike Count per s(Hz)');
            end
            if i==length(AlignedRaster)
                xlabel('Time from event onset');
            end
             title(SegEventStr{i});
            set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
            %}
        end
         ylim([0,max(  Interval_max)]);
    end
    FigureIndex=FigureIndex+1;
    FigureStartNum=FigureStartNum+1; 
    
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