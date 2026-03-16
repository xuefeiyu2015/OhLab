function RTGoodBadAnalysis(Data, StartTrial, EndTrial,ShowFigureFlag,OutputFlag);
%Function to analysis neural responses and target property and saccades
%properties
%Xuefei Yu 05122020

ECode;
SegEvent=[TGTCD];%Event Flow of interest
%SegEventEnd=[SHOWFIXCD,FIXOFF,TGTOFF,RWDOFFCD ,NaN,EXTRAFP];%Event end
SegEventEnd=[TGTOFF];%Event end
SegEventStr={'TargetOn'};
AmpMatchFlag=0;%1:Amplitude match on ;0:Amplitude Match off


PlotStartOffset=[1000];%in ms
PlotStopOffset=[8000];%in ms

ExamplePlot=0;%PLots for presentations, turn off for regular use

ResponseLatencyThreshold=3;

%Parameter to define the sliding window to measure firing rates
StepSize=5;%in ms
BinWidth=15;%in ms


%Criterias for sexwlecting saccades
V_Threshold=40;
ContinueBin=5;
ISI_Threshold=20;

%Time bin pre saccade and time bin offset after saccade onset
PlotSacPreBuffer=500;%in ms
PlotSacPostBuffer=500;%in ms

PlotFracOffPreBuffer=500;%in ms
PlotFracOffPostBuffer=2000;%in ms


%Color Definition for each fractals
%{ee
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
Selection=StartTrial:EndTrial;


%Abnormal check
CheckCodePair=[TGTCD,TGTOFF];%This code should be the same as target on for good trials
AbnormalTrials=AbnormalTrialCheck(Data,CheckCodePair);
Selection=Selection(~ismember(Selection,AbnormalTrials));
    

%Remove choice trials
try
TaskType = Data.ParaLib('TaskType');
%Only maintain the sigle target trials
%SingleTarget =TaskType(TaskType(:,2) == 1,1) ;
SingleTarget =find(TaskType(:,2) == 1) ;
catch
    SingleTarget = 1:length(Selection);
end

Selection = Selection(ismember(Selection,SingleTarget));



%SpikeChannel=Data.SpikeChannel;




%FixationPointIndex=Data.FixPointIndex(Selection,:);
%FixLoc=Data.FixLoc(Selection,:,:);

try
tmp=Data.ParaLib('FlipGoodBad');
FlipGoodBad=tmp(Selection,2);
catch
   FlipGoodBad=0*ones(sum(Selection),1);
   disp('Auto no flipping sign');
end


%Target Information
TargetIndex=Data.TargetIndex(Selection,2);%For this task, there is only target name, but no target index
try
TargetName=Data.TargetIndex(Selection,3);%For this task, there is only target name, but no target index
catch
    TargetName=TargetIndex;
    disp('No targetname');
end
%TargetLoc=Data.TargetLoc(Selection,:,:);

TargetNum=Data.TargetNum(Selection,2:end);

%Rearrange target into lines
TargetIndex_Line=TargetIndex;

TargetIndex_Line=reshape(TargetIndex_Line,numel(TargetIndex_Line),1);

NoNNaNSelector=~isnan(TargetIndex_Line);

TargetIndex_Line=TargetIndex_Line(~isnan(TargetIndex_Line));


TargetName_Line=TargetName;

TargetName_Line=reshape(TargetName_Line,numel(TargetName_Line),1);

TargetName_Line=TargetName_Line(~isnan(TargetName_Line));
%Analysis Good versus Bad

GoodObjectIndex_Raw=TargetIndex_Line<=4;
%%
%Flip the sign of good and bad object according to flipping index

GoodObjectIndex_Raw(FlipGoodBad==1)=~GoodObjectIndex_Raw(FlipGoodBad==1);

%% Target location

TargetLoc=Data.TargetLoc(Selection,2:end);
%TargetPolar=Data.TargetPolarLoc(Selection,:,2:end);
if length(size(Data.TargetPolarLoc)) == 3

    TargetAngle = Data.TargetPolarLoc(Selection,2,1);

    TargetEcc =  Data.TargetPolarLoc(Selection,2,2);
else
    TargetAngle = Data.TargetPolarLoc(Selection,2);

    TargetEcc =  Data.TargetPolarLoc(Selection,3);

end

UniqueAngle = unique(TargetAngle);
UniqueEcc = unique(TargetEcc);





%%
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
EventBin=Data.EventBin(Selection,:);




%%Reproduce goodbad index from event channel
%ObjectIndex=
ObjectIndex=ReproduceFromEvent(EventChannel,[GOODOBJ  BADOBJ]);

GoodObjectIndex=(ObjectIndex==GOODOBJ)';




%Select the period between target onset and target off set
FractalOnsetBin = FindOutTime(EventChannel,EventBin,TGTCD );
FractalOffsetBin = FindOutTime(EventChannel,EventBin,TGTOFF);
FractalOnsetTime = FindOutTime(EventChannel,EventTimeChannel,TGTCD );
FractalOffsetTime = FindOutTime(EventChannel,EventTimeChannel,TGTOFF );

FixOffTime = FindOutTime(EventChannel,EventTimeChannel,FIXOFF );

FixOffTime_Rel=FixOffTime-FractalOnsetTime;

EyeChannel_X_Selection= SelectEyeInterval(EyeChannel_X_L,FractalOnsetBin ,FractalOffsetBin )';
EyeChannel_Y_Selection= SelectEyeInterval(EyeChannel_Y_L,FractalOnsetBin ,FractalOffsetBin )';

EyeChannel_X_Selection=ReorganizeEye(EyeChannel_X_Selection);
EyeChannel_Y_Selection=ReorganizeEye(EyeChannel_Y_Selection);

EyeChannel_X_Selection=EyeChannel_X_Selection(NoNNaNSelector,:);
EyeChannel_Y_Selection=EyeChannel_Y_Selection(NoNNaNSelector,:);





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

AmplitudeThreshold=3;

TimeThreshold=50;

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






%Good versus Bad raster
%Raster_SacOn_Good = Raster_SacOn(GoodObjectIndex,:);
%Raster_SacOn_Bad = Raster_SacOn(~GoodObjectIndex,:);

%Align Eye trace with saccade onset
%Select the eye eccentricity between target onset and target offset
Eye_Eccentricity_Tgt=SelectEyeInterval(Eye_Eccentricity,FractalOnsetBin,FractalOffsetBin)';
Eye_Eccentricity_Tgt=ReorganizeEye(Eye_Eccentricity_Tgt);
Eye_Eccentricity_Tgt=Eye_Eccentricity_Tgt(NoNNaNSelector,:);

%Select the interval around saccade onset
EyeTime=0:EyeBinWidth:size(Eye_Eccentricity_Tgt,2)-1;

EyeEcc_Aligned=AlignEyeData(Eye_Eccentricity_Tgt,EyeTime,RT_First_Target,[-PlotSacPreBuffer,PlotSacPostBuffer]);



EyeEcc_Aligned_mean=nanmean(EyeEcc_Aligned,1);
EyeEcc_Aligned_sem=nanstd(EyeEcc_Aligned,[],1)./sqrt(sum(~isnan(EyeEcc_Aligned),1));

%Good versus Bad

EyeEcc_Aligned_Good_mean=nanmean(EyeEcc_Aligned(GoodObjectIndex,:),1);
EyeEcc_Aligned_Good_sem=nanstd(EyeEcc_Aligned(GoodObjectIndex,:),[],1)./sqrt(sum(~isnan(EyeEcc_Aligned(GoodObjectIndex,:)),1));

EyeEcc_Aligned_Bad_mean=nanmean(EyeEcc_Aligned(~GoodObjectIndex,:),1);
EyeEcc_Aligned_Bad_sem=nanstd(EyeEcc_Aligned(~GoodObjectIndex,:),[],1)./sqrt(sum(~isnan(EyeEcc_Aligned(~GoodObjectIndex,:)),1));


%Good versus Bad raster
%Raster_TgtOn_Good = Raster_TgtOn(GoodObjectIndex,:);
%Raster_TgtOn_Bad = Raster_TgtOn(~GoodObjectIndex,:);


%Get mean and sem of eye eccentricity
%Get eye trance distance during the interval
StartTimeBin=FindOutTime(EventChannel,EventTimeChannel,SegEvent)-PlotStartOffset;
StartTimeBin(StartTimeBin<1)=1;
EndTimeBin=FindOutTime(EventChannel,EventTimeChannel,SegEventEnd )+PlotStopOffset;
  
Eye_Eccentricity_FractalOn=SelectEyeInterval(Eye_Eccentricity,StartTimeBin,EndTimeBin)';
Eye_Eccentricity_FractalOn=ReorganizeEye(Eye_Eccentricity_FractalOn);
Eye_Eccentricity_FractalOn=Eye_Eccentricity_FractalOn(NoNNaNSelector,:);





Eye_Eccentricity_FractalOn_mean=nanmean(Eye_Eccentricity_FractalOn,1);
Eye_Eccentricity_FractalOn_sem=nanstd(Eye_Eccentricity_FractalOn,1)./sqrt(sum(~isnan(Eye_Eccentricity_FractalOn),1));

%Eye eccentricity Good v.s. Bad
Eye_Eccentricity_FractalOn_Good=Eye_Eccentricity_FractalOn(GoodObjectIndex,:);
Eye_Eccentricity_FractalOn_Bad=Eye_Eccentricity_FractalOn(~GoodObjectIndex,:);


Eye_Ecc_Good_mean=nanmean(Eye_Eccentricity_FractalOn_Good,1);
Eye_Ecc_Bad_mean=nanmean(Eye_Eccentricity_FractalOn_Bad,1);

Eye_Ecc_Good_sem=nanstd(Eye_Eccentricity_FractalOn_Good,[],1)./sqrt(sum(~isnan(Eye_Eccentricity_FractalOn_Good),1));
Eye_Ecc_Bad_sem=nanstd(Eye_Eccentricity_FractalOn_Bad,[],1)./sqrt(sum(~isnan(Eye_Eccentricity_FractalOn_Bad),1));

RF_FixOff=RT_First_Target' - FixOffTime_Rel;

RT_Good=RF_FixOff(GoodObjectIndex);
RT_Bad=RF_FixOff(~GoodObjectIndex);

RT_Good_Mean = nanmean(RT_Good);
RT_Bad_Mean = nanmean(RT_Bad);

RT_Good_Sem = nanstd(RT_Good)/sqrt(length(RT_Good));
RT_Bad_Sem = nanstd(RT_Bad)/sqrt(length(RT_Good));


[h,p]=ttest2(RT_Good,RT_Bad);




%{

%% Loop over the target directions
for i = 1:length(UniqueAngle )
    for j = 1:length(UniqueEcc)
        sel = TargetAngle == UniqueAngle(i) & TargetEcc == UniqueEcc(j);
        
       

           
        



        
    end
end

%}






%%
if ShowFigureFlag
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Figure 1: Average firing rate
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%for each location
figtitlestr{1}='RT_GoodBad';

fig{1}=PrepareFigure(110,'w',[50,100, 1200,800],'Name',figtitlestr{1});

bar([1,2],[RT_Good_Mean,RT_Bad_Mean]);
hold on
errorbar([1,2],[RT_Good_Mean,RT_Bad_Mean],[RT_Good_Sem,RT_Bad_Sem],'sk');


end %End of ShowFigureFlag



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
%ChannelNum=SpikeChannel;
%ChannelNum=Data.SpikeChannel;
ChannelNum=1;


%To get the monkey name
RemainFiles=erase(FileName,RecordDateOriginal);
Delimiter=find(isstrprop(RemainFiles,'upper')==1);
MonkeyName=RemainFiles(Delimiter(1):Delimiter(2)-1);


%Output file name
OutputFileName=sprintf('%s_%s_N%s_C%s.mat',MonkeyName,string(RecordDateOriginal),string(NeuronNum),string(ChannelNum));
OutputFigureName=sprintf('%s_%s_N%s_C%s',MonkeyName,string(RecordDateOriginal),string(NeuronNum),string(ChannelNum));
%{
%Store figures in common folders
if TaskCode==OBJ_DIR_SAC_TASK  
ProtocolMarker='SatbleValue';
elseif  TaskCode==FLEXIBLE_VALUE_TASK 
    ProtocolMarker='FlexibleValue';
end
FigurePath=strcat(BasicDirectory,'Figures',BasicDirectory(end),ProtocolMarker,BasicDirectory(end));%,OutputFolerName);
if ~exist(FigurePath)
    mkdir(FigurePath)
end

if ShowFigureFlag & OutputFlag
    for i=1:length(fig)
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


ExistFlag=0;
%Load the old file if exist
if TaskCode==OBJ_DIR_SAC_TASK  
    ProtocolMarker='StableValueBehavior';
elseif  TaskCode==FLEXIBLE_VALUE_TASK 
    ProtocolMarker='FlexibleValueBehavior';
end
if exist(OutputFileName)

    load(OutputFileName);
    %Load OutputData into the memory
    if isfield(OutputData,ProtocolMarker)
    ExistFlag=1;
    end
end


%%%%%%%%%%%%%%%%%%%%OutputData For current analysis
OutputData_New=[];
if TaskCode==OBJ_DIR_SAC_TASK 
    OutputData_New.StableValueBehavior.TrialType='GoodTrials';
elseif  TaskCode==FLEXIBLE_VALUE_TASK 
    OutputData_New.FlexibleValueBehavior.TrialType='GoodTrials';
end

%%%%%%%%%Save overall scene response characteristics%%%%%%%%%%%%%%%%%%%%%%%
%Organize the data

NameStr={
'RT_Good','RT_Bad','RT_Good_Mean','RT_Bad_Mean'
  

};

DataLib={RT_Good,RT_Bad,RT_Good_Mean,RT_Bad_Mean

};


DataStamp=containers.Map(NameStr,DataLib);




if TaskCode==OBJ_DIR_SAC_TASK 
    OutputData_New.StableValueBehavior.Task=TaskType;
    OutputData_New.StableValueBehavior.TaskCode=TaskCode;
    OutputData_New.StableValueBehavior.DataStamp=DataStamp;



if ExistFlag
%Compare the old one with the new one
Task_Old=OutputData.StableValueBehavior.TaskCode;

     if Task_Old~=TaskCode
        %Not the same task, add the new dataset into the old one
        OutputData.StableValueBehavior(numel(OutputData.StableValueBehavior)+1)=OutputData_New.StableValueBehavior;
     else
         %If the same task,the same protocol,replace the old one with the new one
         
         OutputData.StableValueBehavior=OutputData_New.StableValueBehavior;
         
    
    
     end
else
    %Output the current file
   OutputData.StableValueBehavior=OutputData_New.StableValueBehavior; 
end

elseif  TaskCode==FLEXIBLE_VALUE_TASK 
        OutputData_New.FlexibleValueBehavior.Task=TaskType;
        OutputData_New.FlexibleValueBehavior.TaskCode=TaskCode;
        OutputData_New.FlexibleValueBehavior.DataStamp=DataStamp;



if ExistFlag
%Compare the old one with the new one
Task_Old=OutputData.FlexibleValueBehavior.TaskCode;

     if Task_Old~=TaskCode
        %Not the same task, add the new dataset into the old one
        OutputData.FlexibleValueBehavior(numel(OutputData.FlexibleValueBehavior)+1)=OutputData_New.FlexibleValueBehavior;
     else
         %If the same task,the same protocol,replace the old one with the new one
         
         OutputData.FlexibleValueBehavior=OutputData_New.FlexibleValueBehavior;
         
    
    
     end
else
    %Output the current file
   OutputData.FlexibleValueBehaviore=OutputData_New.FlexibleValueBehavior; 
end


end



if OutputFlag
    cd(OutPath);
    save(OutputFileName,'OutputData');
    disp('Data Exported');
else
    disp('No data export');
end  

end %End of the function

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
function FR=SelectFringRateInterval(FR_Raw,Time,Interval);
for i=1:size(Interval,1)
    if sum(~isnan(Interval(i,:)))==2
        FR(i)=nanmean(FR_Raw(i,Time>=Interval(i,1) & Time<=Interval(i,2)));
    else
         FR(i)=NaN;
    end
    
    
end




end
function ObjectIndex=ReproduceFromEvent(event,code)
for i=1:size(event,1)
    
    ObjectIndex(i)=event(i,ismember(event(i,:),code));
    
  
    
end

end