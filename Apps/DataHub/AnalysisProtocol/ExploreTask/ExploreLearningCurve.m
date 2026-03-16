function ExploreLearningCurve(Data, StartTrial, EndTrial,ShowFigureFlag,OutputFlag);
%Protocol to analysis the learning curve 
%by Xuefei Yu 10262020
%Xuefei Yu, 03,21,2020
%Unpack all the data base and task settings
ECode;
TaskType=Data.TaskType;
TaskCode=Data.TaskCode;
TrialNum=EndTrial-StartTrial+1;
Selection=StartTrial:EndTrial;%Select the trials according to the pannel
%Abnormal check %Especially For Robin081120FORAGEXPLO4.rst
CheckCode=SCENE_ON;%This code should appear once in each trial
AbnormalTrials=AbnormalTrialCheck(Data,CheckCode);

TargetOnOffset=0;%in ms, the time selection bin after target onset; 
StartDelay=0;%in ms, the time selection bin after fixation offset;

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


%BlockType(BlockType(:,3)>0,3)=999;

[UniqueBlockType,orderindex,reprojectorderindex]=unique(BlockType,'row');
NumberOfBlocks=max(reprojectorderindex);


PlotEveryNumRep=5;%Plot performance every 5 repetitions

if OutputFlag==1
h=questdlg('Output data?',...
    'Yes','No');
waitfor(h);
if strcmp(h,'No')
    OutputFlag=0;
    
end
end




%Get the searching time
%GetEyeData
EyeData=Data.EyeDataRaw;
EyeChannel_X_L=EyeData(Selection,1);
EyeChannel_Y_L=EyeData(Selection,2);
%EyeChannel_X_R=EyeData(Selection,3);
%EyeChannel_Y_R=EyeData(Selection,4);

EyeBinWidth=Data.EyeBinWidth;

%Load event data
EventChannel=Data.EventChannel(Selection,:);
EventTimeChannel=Data.EventTimeChannel(Selection,:);
EventBin=Data.EventBin(Selection,:);

try
%Recalibrate eyes according to the extra FP
EyeCalInfo=Data.EyeCalInfo;%Each row: X_Gain,X_Bias; Y_Gain; Y_bias; X=(raw+X_bias)*X_Gain
EyeCalInfo_L=EyeCalInfo(1,:);

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
catch
    disp('Skip offline eye calibration');
end




%Get scene and fractals information
SceneName=Data.SceneIndex(Selection,3);
TargetName=Data.TargetIndex(Selection,:,2);
TargetLoc=Data.TargetLoc(Selection,2:end,:);

%TargetNum=Data.TargetNum(Selection,2:end);
TargetNum=sum(~isnan(TargetName(:,2:end)),2);
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
 
 %If the target onset time and fixofftime dimension are not aligned 
 if size(TargetOnTime,2)~=size(FixOffTime,2)
     FixOffTime=FixOffTime(:,1:size(TargetOnTime,2));
     disp('Warning!FixOffTime Corrected!!')
 end
     
 
SearchingTime= (TargetOnTime-FixOffTime)/1000;%Change into s for presentation

ColorForEachMode=[0.83 0.14 0.14
             1.00 0.54 0.00
             0.47 0.25 0.80
             0.25 0.80 0.54];
AssociateBlock=[10,11,0,1];
Blockstr={'FixLoc,FixSceneSeq','FixLoc,RandSceneSeq','RandLoc,FixSceneSeq','RandLoc,RandSceneSeq'};


%Color Definition for each scenes
Color={'r','g','b','m',[0.4,0.4,0],[0.28,0.78,0.89],'k',[0.6,0.1,1],[0.6,0.1,0.4],[0.3,0.5,1],[0.6,0.5,0.1],[0,0.8,0.4],'r','g','b','m',[0.4,0.4,0],[0.28,0.78,0.89],'k',[0.6,0.1,1],[0.6,0.1,0.4],[0.3,0.5,1],[0.6,0.5,0.1],[0,0.8,0.4]};

%Setup scene library
Workingdirectory=pwd;
MarkerFolder='DataAnalysis';
Flag=strfind(Workingdirectory,MarkerFolder);
BasicDirectory=Workingdirectory(1:Flag+length(MarkerFolder));
ImageFolerName='VisualStimulus';
%ImagePath=strcat(BasicDirectory,ImageFolerName,BasicDirectory(end));

for i=1:NumberOfBlocks
    CurrentBlockType=UniqueBlockType(i,:);
    BlockLogicSelection=BlockType==CurrentBlockType;
    
    BlockOrder(i)=find(CurrentBlockType==AssociateBlock);
    
    RepStart=1;
    %{
    SceneNameCurrBlock=SceneName( BlockLogicSelection,:);
    UniqueSceneNameInCurrentBlock=unique(SceneNameCurrBlock);
  
     %Number of repetition
     RepNum=ceil(sum(BlockLogicSelection)/length(UniqueSceneNameInCurrentBlock)/PlotEveryNumRep);
     
     %Number of condition in each repetition
    [temp_loc,NumInRep,RepinScene]=FindNumOfRepetition(SceneNameCurrBlock);
    %}
    SearchingTime_currBlock=SearchingTime(BlockLogicSelection,1);
    
    RepNum=ceil(sum(BlockLogicSelection)/PlotEveryNumRep);
    
    for j=1:RepNum %Loop over the repetition
       % RepSelection=RepStart:min(RepStart+NumInRep(j)*PlotEveryNumRep-1,sum(BlockLogicSelection));  
       RepSelection=RepStart:min(RepStart+PlotEveryNumRep-1,sum(BlockLogicSelection));  
        
        SearchingTime_currBlock_rep=SearchingTime_currBlock(RepSelection);
        
        SearchingTime_currBlock_rep_mean(i,j)=nanmean(SearchingTime_currBlock_rep);
        SearchingTime_currBlock_rep_sem(i,j)=nanstd(SearchingTime_currBlock_rep)/sqrt(sum(~isnan(SearchingTime_currBlock_rep)));
        
        
       % RepStart=RepStart+NumInRep(j);
       RepStart=RepStart+PlotEveryNumRep;
   
    end
   
   
    
    
end

if ShowFigureFlag
%Bar plot for each block

if ishandle(103)
    close 103
end

titlestr{1}='SearchingLearningCurve';

fig{1}=figure(103);

set(103,'color','w','Position', [100,50 1300,800],'Name',titlestr{1})
str=strcat('Behavior analysis of ',TaskType);
title(str)
axis off


for i=1:NumberOfBlocks
   
    errorbar(SearchingTime_currBlock_rep_mean(i,:),SearchingTime_currBlock_rep_sem(i,:),'color',ColorForEachMode( BlockOrder(i),:),'LineWidth',3);
    hold on
    
end
legend(Blockstr(BlockOrder));
xlabel(sprintf('Every %d repetition',PlotEveryNumRep));
ylabel('Mean Rearching Time(s)');
set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');


end %End of the if plotfigure

%{
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Export Files
%Output to files
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
TaskShortFlag=find(TaskType=='_',1,'last');
TaskShort=TaskType(TaskShortFlag+1:end);
OutputFolerName=strcat('TrainingData_',TaskShort);
OutPath=strcat(BasicDirectory,'Results',BasicDirectory(end),OutputFolerName);
if ~exist(OutPath)
    mkdir(OutPath)
end
cd(OutPath);
FileName=Data.FileName;
FileName=FileName(1:end-4);


%Get File Information
num = regexp(FileName, '\d+', 'match');
RecordDateOriginal=cell2mat(num(1));


%To get the monkey name
RemainFiles=erase(FileName,RecordDateOriginal);
Delimiter=find(isstrprop(RemainFiles,'upper')==1);
MonkeyName=RemainFiles(Delimiter(1):Delimiter(2)-1);


%Output file name
OutputFileName=sprintf('%s_%s_%s.mat',MonkeyName,string(RecordDateOriginal),TaskShort);
%OutputFigureName=sprintf('%s_%s_N%s_C%s',MonkeyName,string(RecordDateOriginal));


%%%%%%%%%Save overall results %%%%%%%%%%%%%%%%%%%%%%%
%Organize the data
NameStr={'PlotEveryNumRep','PropGoodChoiceInNRep','PropBadChoiceInNRep','PropGoodAcceptInNRep','PropBadAcceptInNRep','PropGoodRejectInNRep','PropBadRejectInNRep'

};





DataLib={PlotEveryNumRep,PropGoodChoiceInNRep,PropBadChoiceInNRep,PropGoodAcceptInNRep,PropBadAcceptInNRep,PropGoodRejectInNRep,PropBadRejectInNRep
};


DataStamp=containers.Map(NameStr,DataLib);

OutputData.Performance.Task=TaskType;
OutputData.Performance.DataStamp=DataStamp;

if OutputFlag
   cd(OutPath);
    save(OutputFileName,'OutputData');
    disp('Data Exported');  
else
    disp('No data export');
end

%}



end