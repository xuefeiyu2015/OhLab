function DelaySaccadeVaryPeriodStim(Data, StartTrial, EndTrial,ShowFigureFlag,OutputFlag)
ECode;

TaskType=Data.TaskType;
TaskCode=Data.TaskCode;
DataPath= Data.Path;
FileName=Data.FileName;

Selection=StartTrial:EndTrial;%Select the trials according to the pannel

%Abnormal check for equality
%
CheckCode=[TGTCD];
AbnormalTrials=AbnormalTrialCheck(Data,CheckCode);

Selection=Selection(~ismember(Selection,AbnormalTrials));
%}
CheckCode=[FIXOFF];
AbnormalTrials=AbnormalTrialCheck(Data,CheckCode);

Selection=Selection(~ismember(Selection,AbnormalTrials));




%Criterias for sexwlecting saccades
V_Threshold=40;
ContinueBin=5;
ISI_Threshold=20;


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

EyeTime = 0: EyeBinWidth:max(cellfun(@numel,EyeChannel_X_L))-EyeBinWidth;

%Load event data
EventChannel=Data.EventChannel(Selection,:);
EventTimeChannel=Data.EventTimeChannel(Selection,:);
EventBin=Data.EventBin(Selection,:);


%Load stimulation related data

ParaLib=Data.ParaLib;
if ismember('StimType',keys(ParaLib))
    StimType=ParaLib('StimType');
    StimType= StimType(Selection,:);
    StimType= StimType(:,2);
    UnqiueStimType = unique(StimType);

end

%Sepera
% 
% te into stim trial and sham trial
%{
StimType=ReproduceFromEvent(EventChannel,[STIMT1,SHAM0,])';
StimTrial = StimType == STIM1;
ControlTrial = ~StimTrial;
%}

%{
SHAM1=196;
SHAM2=197;
SHAM3=196;
SHAM4=197;
SHAM5=196;
SHAM6=197;

STIMT1=204;
STIMT2=205;
STIMT3=206;
STIMT4=207;
STIMT5=208;
STIMT6=209;
%}



%Load target related parameter
TargetLoc=Data.TargetLoc(Selection,2:end);
TargetPolar=Data.TargetPolarLoc(Selection,2:end);
TargetAngle = TargetPolar(:,1);
TargetEcc =  TargetPolar(:,2);

UniqueAngle = unique(TargetAngle);
UniqueEcc = unique(TargetEcc);
minEcc = min(TargetEcc);



%First analysis the behavior: compare the saccade latency or target
%acqusition time between stim and non-stim trials


%Select the period between fix off and target off set
FIXOFFBin = FindOutTime(EventChannel,EventBin,FIXOFF )-500;
SACENDBin = FIXOFFBin + 1500;


EyeChannel_X_Selection= SelectEyeInterval(EyeChannel_X_L,FIXOFFBin,SACENDBin)';
EyeChannel_Y_Selection= SelectEyeInterval(EyeChannel_Y_L,FIXOFFBin,SACENDBin)';

EyeChannel_X_Selection=ReorganizeEye(EyeChannel_X_Selection);
EyeChannel_Y_Selection=ReorganizeEye(EyeChannel_Y_Selection);


%Analysis saccades for each trial
Saccades=SelectSaccade(EyeChannel_X_Selection,EyeChannel_Y_Selection,EyeBinWidth,V_Threshold,ContinueBin,ISI_Threshold);
SaccadeNumber=arrayfun(@(x)  x.NumOfSaccade,Saccades)';
SaccadeAngle=arrayfun(@(x)  x.SaccadeAngle,Saccades,'UniformOutput' ,0)';
SaccadeAmplitude=arrayfun(@(x)  x.SaccadeAmplitude,Saccades,'UniformOutput' ,0)';

SaccadeStartTime=arrayfun(@(x)  x.SaccadeStartTime,Saccades,'UniformOutput' ,0)';
SaccadeEndTime=arrayfun(@(x)  x.SaccadeEndTime,Saccades,'UniformOutput' ,0)';

SaccadeStartPoint=arrayfun(@(x)  x.SaccadeStartPoint,Saccades,'UniformOutput' ,0)';
SaccadeEndPoint=arrayfun(@(x)  x.SaccadeEndPoint,Saccades,'UniformOutput' ,0)';

%Screen out the valid saccades and saccade start time above 5;

AmplitudeThreshold=minEcc/3;

TimeThreshold=0;

SaccadeAngle_Selection=cellfun(@(x,y,z) z(((x>AmplitudeThreshold)&(y>TimeThreshold))),SaccadeAmplitude,SaccadeStartTime,SaccadeAngle,'UniformOutput' ,0);
SaccadeEndPoint_Selection=cellfun(@(x,y,z) z(((x>AmplitudeThreshold)&(y>TimeThreshold)),:),SaccadeAmplitude,SaccadeStartTime,SaccadeEndPoint,'UniformOutput' ,0);
SaccadeStartPoint_Selection=cellfun(@(x,y,z) z(((x>AmplitudeThreshold)&(y>TimeThreshold)),:),SaccadeAmplitude,SaccadeStartTime,SaccadeStartPoint,'UniformOutput' ,0);

SaccadeStartTime_Selection=cellfun(@(x,y,z) z(((x>AmplitudeThreshold)&(y>TimeThreshold)),:),SaccadeAmplitude,SaccadeStartTime,SaccadeStartTime,'UniformOutput' ,0);
SaccadeEndTime_Selection=cellfun(@(x,y,z) z(((x>AmplitudeThreshold)&(y>TimeThreshold)),:),SaccadeAmplitude,SaccadeStartTime,SaccadeEndTime,'UniformOutput' ,0);


SaccadeAmplitude_Selection=cellfun(@(x,y,z) z(((x>AmplitudeThreshold)&(y>TimeThreshold)),:),SaccadeAmplitude,SaccadeStartTime,SaccadeAmplitude,'UniformOutput' ,0);
%Only look at the first saccade
SaccadeNumInEachTrial=cellfun(@numel,SaccadeAngle_Selection);


SaccadeAngle_First=NaN*ones(1,length(SaccadeNumInEachTrial));
SaccadeEndPoint_First=NaN*ones(length(SaccadeNumInEachTrial),2);
RT_First_Target=NaN*ones(1,length(SaccadeNumInEachTrial));
SacEnd_First_Target=NaN*ones(1,length(SaccadeNumInEachTrial));

SaccadeStartPoint_First=NaN*ones(length(SaccadeNumInEachTrial),2);

SaccadeAmplitude_First=NaN*ones(1,length(SaccadeNumInEachTrial));

SaccadeAngle_First(SaccadeNumInEachTrial>0)=cell2mat(cellfun(@(x) x(1,:),SaccadeAngle_Selection(SaccadeNumInEachTrial>0),'UniformOutput' ,0));
SaccadeEndPoint_First(SaccadeNumInEachTrial>0,:)=cell2mat(cellfun(@(x) x(1,:),SaccadeEndPoint_Selection(SaccadeNumInEachTrial>0,:),'UniformOutput' ,0));
RT_First_Target(SaccadeNumInEachTrial>0)=cell2mat(cellfun(@(x) x(1,:),SaccadeStartTime_Selection(SaccadeNumInEachTrial>0,:),'UniformOutput' ,0));
SacEnd_First_Target(SaccadeNumInEachTrial>0)=cell2mat(cellfun(@(x) x(1,:),SaccadeEndTime_Selection(SaccadeNumInEachTrial>0,:),'UniformOutput' ,0));


SaccadeStartPoint_First(SaccadeNumInEachTrial>0,:)=cell2mat(cellfun(@(x) x(1,:),SaccadeStartPoint_Selection(SaccadeNumInEachTrial>0,:),'UniformOutput' ,0));

SaccadeAmplitude_First(SaccadeNumInEachTrial>0)=cell2mat(cellfun(@(x) x(1,:),SaccadeAmplitude_Selection(SaccadeNumInEachTrial>0,:),'UniformOutput' ,0));


RT_TargetToFixOff = RT_First_Target - 500;

%% For neural responses 
% Align on three period: STIMON, TgtOn, FixOff

StepSize=5;%in ms
BinWidth=10;%in ms



SegEvent=[ON_STIM,TGTCD,FIXOFF];%Event Flow of interest
%SegEvent=[ON_STIM];
%SegEventEnd=[SHOWFIXCD,FIXOFF,TGTOFF,RWDOFFCD ,NaN,EXTRAFP];%Event end
SegEventEnd=[OFF_STIM,NaN,NaN];%Event end
%SegEventEnd=[PST_STIM];%Event end
SegEventStr={'StimOn','TargetOn','FixOff'};



PlotStartOffset=[300,300,300];%in ms
PlotStopOffset=[200,300,500];%in mss

%Load the spike channel
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
    ChannelID=ChannelIDWhole(spk);
    ClusterID=ClusterIDWhole(spk);
    
    
    %%%%%%%%%%%%%%%%%%%%%%%Align data on each events


SegmentData=EventSegTool(SegEvent,SegEventEnd,PlotStartOffset,PlotStopOffset,EventChannel,EventTimeChannel,SpikeTime,StepSize,BinWidth);

%Retrive data from the field
%PSTH for each event

MeanPSTH_Each =SegmentData.MeanPSTHEach;
MeanPSTH_StimOn = MeanPSTH_Each{1,2};
TimePSTH_StimOn = MeanPSTH_Each{1,1};


MeanPSTH=SegmentData.MeanPSTH;%=[PSTH_Event_Time',PSTH_Event'];

PSTH_Time_StimOn = MeanPSTH{1,1};
PSTH_FR_StimOn = MeanPSTH{1,2};

RelativeTimeMarkerStim = SegmentData.RelativeTimeMarker{1,1};

PSTH_Time_TgtOn = MeanPSTH{2,1};
PSTH_FR_TgtOn = MeanPSTH{2,2};


PSTH_Time_FixOff = MeanPSTH{3,1};
PSTH_FR_FixOff = MeanPSTH{3,2};
  
MaxStimOn = 0;
MaxTgtOn = 0;
MaxFixOff = 0;

%Select the raw eye trace 
%Select the eye trace during the 100pre-stim, 100 stim, 100 post stim
%period


StimStartBin=FindOutTime(EventChannel,EventBin,PRE_STIM);
StimEndBin=FindOutTime(EventChannel,EventBin,PST_STIM );

EyeChannel_X_Stim= SelectEyeInterval(EyeChannel_X_L,StimStartBin,StimEndBin)';
EyeChannel_Y_Stim= SelectEyeInterval(EyeChannel_Y_L,StimStartBin,StimEndBin)';

EyeChannel_X_Stim=ReorganizeEye(EyeChannel_X_Stim);
EyeChannel_Y_Stim=ReorganizeEye(EyeChannel_Y_Stim);

EyeTime_Stim = 0:EyeBinWidth:size(EyeChannel_X_Stim,2)-EyeBinWidth;




MarkerTime = FindOutTime(EventChannel,EventBin,ON_STIM)-StimStartBin;

Interval = [-100,nanmean(StimEndBin-FindOutTime(EventChannel,EventBin,ON_STIM))];
%Interval=[-100,200];
StimDur = [0,nanmean(FindOutTime(EventChannel,EventBin,OFF_STIM)- FindOutTime(EventChannel,EventBin,ON_STIM))];

EyeChannel_X_StimAlign = AlignEyeData(EyeChannel_X_Stim,EyeTime_Stim,MarkerTime,Interval);
EyeChannel_Y_StimAlign = AlignEyeData(EyeChannel_Y_Stim,EyeTime_Stim,MarkerTime,Interval);

%NormX = nanmean(EyeChannel_X_StimAlign(:,1:100),2);
%NormY = nanmean(EyeChannel_Y_StimAlign(:,1:100),2);

%EyeChannel_X_StimAlign_Norm = EyeChannel_X_StimAlign-NormX;
%EyeChannel_Y_StimAlign_Norm = EyeChannel_Y_StimAlign-NormY;

EyeTime_Stim_Align=Interval(1):EyeBinWidth:Interval(2)+EyeBinWidth;

MaxStimOn = 0*ones(length(UnqiueStimType)-1,length(UniqueAngle));


%% Compare between different conditions
for i = 1:length(UnqiueStimType)
    currStim = StimType == UnqiueStimType(i);
     %Eye trace 
       EyeX_Stim_loc{i} = EyeChannel_X_StimAlign(currStim,:);
       EyeY_Stim_loc{i} = EyeChannel_Y_StimAlign(currStim,:);
       
     
       

    for j = 1:length(UniqueAngle)
        sel_stim = TargetAngle == UniqueAngle(j) & currStim;
       
      
        RT_stim_loc_period = RT_TargetToFixOff(sel_stim);
       

        RT_stim_loc_period_Mean(i,j) = nanmean(RT_stim_loc_period);
        RT_stim_loc_period_Sem(i,j) = nanstd(RT_stim_loc_period)/sqrt(sum(~isnan(RT_stim_loc_period)));

       
      
       
        %StimOn

        if UnqiueStimType(i) == 0
            %Control trials, separate into 6 pieces
            for ii = 1:length(UnqiueStimType)-1
                PSTH_Stim_tmp=MeanPSTH_StimOn{ii};
                PSTH_StimTime_tmp=TimePSTH_StimOn{ii};

                PSTH_Control=PSTH_Stim_tmp(sel_stim,:);
                PSTH_Control_Time = PSTH_StimTime_tmp(sel_stim,:);

                PSTH_Control_Mean{ii,j} = nanmean(PSTH_Control,1);
                PSTH_Control_Sem{ii,j} = nanstd(PSTH_Control,[],1)./sqrt(sum(~isnan(PSTH_Control),1));

                TimePSTH_Control{ii,j} = mean(PSTH_Control_Time,1);


            end
        else
                 PSTH_Stim_tmp=MeanPSTH_StimOn{UnqiueStimType(i)};
                 PSTH_StimTime_tmp=TimePSTH_StimOn{UnqiueStimType(i)};

                 PSTH_Stim=PSTH_Stim_tmp(sel_stim,:);
                 PSTH_Stim_Time = PSTH_StimTime_tmp(sel_stim,:);

                PSTH_Stim_Mean{UnqiueStimType(i),j} = nanmean(PSTH_Stim,1);
                PSTH_Stim_Sem{UnqiueStimType(i),j} = nanstd(PSTH_Stim,[],1)./sqrt(sum(~isnan(PSTH_Stim),1));

                TimePSTH_Stim{UnqiueStimType(i),j} = mean(PSTH_Stim_Time,1);

                max_tmp = max(PSTH_Stim_Mean{UnqiueStimType(i),j});
                if max_tmp>MaxStimOn(UnqiueStimType(i),j)
                    MaxStimOn(UnqiueStimType(i),j)=max_tmp;
                end
            

        end

       
       
       %{
       
       %TgtOn
       PSTH_Time_stim_loc_tgtOn = PSTH_Time_TgtOn(sel_stim,:);
       PSTH_FR_stim_loc_tgtOn = PSTH_FR_TgtOn(sel_stim,:);
       
      
       
       PSTH_stim_loc_tgtOn_Mean{i,j} = nanmean(PSTH_FR_stim_loc_tgtOn );
       PSTH_stim_loc_tgtOn_Sem{i,j} = nanstd(PSTH_FR_stim_loc_tgtOn)./sqrt(sum(~isnan(PSTH_FR_stim_loc_tgtOn),1));
       
      
       PSTH_stim_loc_tgtOn_MeanTime{i,j} = nanmean(PSTH_Time_stim_loc_tgtOn);
      
       
      
       
       %FixOff
       PSTH_Time_stim_loc_fixOff = PSTH_Time_FixOff(sel_stim,:);
       PSTH_FR_stim_loc_fixOff = PSTH_FR_FixOff(sel_stim,:);
       
       
       
       PSTH_stim_loc_fixOff_Mean{i,j} = nanmean(PSTH_FR_stim_loc_fixOff );
       PSTH_stim_loc_fixOff_Sem{i,j} = nanstd(PSTH_FR_stim_loc_fixOff )./sqrt(sum(~isnan(PSTH_FR_stim_loc_fixOff),1));
       
       
       
       PSTH_stim_loc_fixOff_MeanTime{i,j} = nanmean(PSTH_Time_stim_loc_fixOff);
      
       

       %}
        
    end



end



%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% For figure plotting%%%%%%%%%%%%%%%%%
StrStimPeriod = {'DurFix','TgtOn','FixOff','EarlyDelay','LateDelay','Excution'};
%Stim Type:       1,       2,       3,        4,        5(Before fixoff), 6(100m after fix off)          

if ShowFigureFlag
    %% Figure1: PSTH for different stimulation period
   StartNO = 50;
   FigureIndex=0;
   for j = 1:length(UniqueAngle)
       FigureIndex = FigureIndex+1;
       figtitlestr{FigureIndex}=sprintf('PSTH_Stim_Angle %dFor%sC%d',UniqueAngle(j),FileName,SpikeChannel);
    
       fig{FigureIndex}=PrepareFigure(StartNO+FigureIndex-1,'w',[50,100, 1200,800],'Name',figtitlestr{FigureIndex});

       for i =2:length(UnqiueStimType)
           subplot(ceil(length(UnqiueStimType)/2),2,i-1);

          
            AverageRegion = nanmean(RelativeTimeMarkerStim);
            area(AverageRegion,[ MaxStimOn(UnqiueStimType(i),j), MaxStimOn(UnqiueStimType(i),j)],...
             'FaceColor','k','FaceAlpha',.1,'EdgeAlpha',.3,'HandleVisibility','off');
            hold on


           shadedErrorBar(TimePSTH_Control{i-1,j} ,PSTH_Control_Mean{i-1,j},PSTH_Control_Sem{i-1,j},'lineprops',{[0,0,0]},'transparent',1,'patchSaturation',0.3);
           hold on
            shadedErrorBar( TimePSTH_Stim{i-1,j} ,PSTH_Stim_Mean{i-1,j},PSTH_Stim_Sem{i-1,j},'lineprops',{[1,0,0]},'transparent',1,'patchSaturation',0.3);
            set(findall(gca, 'Type', 'Line'),'LineWidth',2);
          %  xlim([-100,200]);
          title(StrStimPeriod(i-1));
          set(gca,'LineWidth',2,'FontSize',10,'FontWeight','Bold');
          if i-1 == length(UnqiueStimType)-1 || i-1 ==length(UnqiueStimType)
              xlabel('Time from opto-stim onset(ms)');
              ylabel('Firing Rate');

          end
          box off;


         
        end



   end
   
   
   

end %End of the if show figure


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Output to
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%files%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

%Store figures in common folders
%{
ProtocolMarker = 'MemorySac';
FigurePath=strcat(BasicDirectory,'Figures',BasicDirectory(end),ProtocolMarker,BasicDirectory(end));%,OutputFolerName);
if ~exist(FigurePath)
    mkdir(FigurePath)
end

if ShowFigureFlag & OutputFlag
    %for i=1:length(fig)
    for i=1:1
        cd(FigurePath);
        OutputFigureNameFig=[];
        titlestr=figtitlestr{i};
        FolderName=titlestr;
         if ~exist(FolderName,'dir')
          mkdir(FolderName);
         end
          cd(FolderName);
    
        OutputFigureNameFig=strcat(OutputFigureName,'.jpg');
        saveas(fig{i},OutputFigureNameFig);
        
        
    end
        
    
    disp('Figures have been exported to the common folder');
    
end

%}

%


ExistFlag=0;
%Load the old file if exist
if exist(OutputFileName)

    load(OutputFileName);
    %Load OutputData into the memory
    if isfield(OutputData,'DelaySaccadeVaryStim')
    ExistFlag=1;
    end
end


%%%%%%%%%%%%%%%%%%%%OutputData For current analysis
OutputData_New=[];

OutputData_New.DelaySaccadeVaryStim.TrialType='GoodTrials';
%%%%%%%%%Save overall scene response characteristics%%%%%%%%%%%%%%%%%%%%%%%
%Organize the data
%{
NameStr={'PSTH_Time_TgtOn','PSTH_TgtOn_RF',...
    'PSTH_Time_SacOn','PSTH_SacOn_RF',...
    'RF_h_v','RF_h_m','CellType',...
    'PSTH_TgtOn_RF_z','PSTH_SacOn_RF_z'
    
    
};

DataLib={PSTH_Time_FraOn_All,PSTH_RF_Mean,...
    TimeSequence_SacOn_mean,PSTH_RF_Sac_Mean,...
    RF_h_v,RF_h_m,CellTypeCode,...
    PSTH_RF_Mean_z,PSTH_RF_Sac_Mean_z

};
%}
NameStr={'TEST'
    
    
};
DataLib={1
};


DataStamp=containers.Map(NameStr,DataLib);





OutputData_New.DelaySaccadeVaryStim.Task=TaskType;
OutputData_New.DelaySaccadeVaryStim.TaskCode=TaskCode;
OutputData_New.DelaySaccadeVaryStim.DataStamp=DataStamp;


%Store figures
if ShowFigureFlag
   
    FolderName = 'DelaySaccadeVaryStimTask';
   
    if ~exist(FolderName,'dir')
        mkdir(FolderName);
    end
    cd(FolderName);
    
    TotalFigure=length(fig);
    %{
    for j=1:TotalFigure
        OutputFigureName=strcat(OutputFigureName,figtitlestr{j},'.jpg');
        if ~isempty(fig{j})
        saveas(fig{j},OutputFigureName);
        end
        
    end
    %}
     OutputFigureName=strcat(OutputFigureName,figtitlestr{2},'.jpg');
        if ~isempty(fig{2})
        saveas(fig{2},OutputFigureName);
        end

   % disp('Figures have been exported to the neuron folder');
   disp('Figure2 have been exported to the common folder');

end

if ExistFlag
%Compare the old one with the new one
Task_Old=OutputData.DelaySaccadeVaryStim.TaskCode;

     if Task_Old~=TaskCode
        %Not the same task, add the new dataset into the old one
        OutputData.DelaySaccadeVaryStim(numel(OutputData.DelaySaccadeVaryStim)+1)=OutputData_New.DelaySaccadeVaryStim;
     else
         %If the same task,the same protocol,replace the old one with the new one
         
         OutputData.DelaySaccadeVaryStim=OutputData_New.DelaySaccadeVaryStim;
         
    
    
     end
else
    %Output the current file
   OutputData.DelaySaccadeVaryStim=OutputData_New.DelaySaccadeVaryStim; 
end



if OutputFlag
    cd(OutPath);
    save(OutputFileName,'OutputData');
    disp('Data Exported');
else
    disp('No data export');
end  






end %End of the spikechannel


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