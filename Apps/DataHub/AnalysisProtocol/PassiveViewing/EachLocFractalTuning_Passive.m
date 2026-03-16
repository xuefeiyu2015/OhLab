function EachLocFractalTuning_Passive(Data, StartTrial, EndTrial,ShowFigureFlag,OutputFlag);
ECode;
SegEvent=[TGTCD];%Event Flow of interest
%SegEventEnd=[SHOWFIXCD,FIXOFF,TGTOFF,RWDOFFCD ,NaN,EXTRAFP];%Event end
SegEventEnd=[TGTOFF];%Event end
SegEventStr={'TargetOn'};
AmpMatchFlag=0;%1:Amplitude match on ;0:Amplitude Match off

FileName_raw = Data.FileName; 

if FileName_raw(1) == 'A'
    Monkey_curr = 2;%Adams
else
    Monkey_curr = 1;%Robin
end


PlotStartOffset=[400];%in ms
PlotStopOffset=[300];%in ms

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
PlotSacPreBuffer=200;%in ms
PlotSacPostBuffer=500;%in ms

PlotFracOffPreBuffer=500;%in ms
PlotFracOffPostBuffer=2000;%in ms


%Color Definition for each fractals
%{
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
    





%SpikeChannel=Data.SpikeChannel;



if strcmp(FileName,'Robin040321FRACSAC1.rst') & SpikeChannel==22
    Selection=Selection(~ismember(Selection,[14,40]));
end

%FixationPointIndex=Data.FixPointIndex(Selection,:);
%FixLoc=Data.FixLoc(Selection,:,:);

%{
try
tmp=Data.ParaLib('FlipGoodBad');
FlipGoodBad=tmp(Selection,2);
catch
   FlipGoodBad=0*ones(sum(Selection),1);
   disp('Auto no flipping sign');
end
%}

%Target Information
TargetIndex=Data.TargetIndex(Selection,2:end,1);%For this task, there is only target name, but no target index
try
TargetName=Data.TargetIndex(Selection,2:end,2);%For this task, there is only target name, but no target index
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

%Rearrange target locations into lines
%% Target location

TargetLoc=Data.TargetLoc(Selection,2:end,:);
TargetPolarLoc=Data.TargetPolarLoc(Selection,2:end,:);
TargetAngle=TargetPolarLoc(:,:,1);
TargetEcc=TargetPolarLoc(:,:,2);

TargetAngle_Line=TargetAngle;
TargetAngle_Line(TargetAngle_Line==360)=0;
TargetAngle_Line=reshape(TargetAngle_Line,numel(TargetAngle_Line),1);
TargetAngle_Line=TargetAngle_Line(NoNNaNSelector);

TargetEcc_Line=TargetEcc;
TargetEcc_Line=reshape(TargetEcc_Line,numel(TargetEcc_Line),1);
TargetEcc_Line=TargetEcc_Line(NoNNaNSelector);

%Analysis Good versus Bad

ObjType=ObjectDefine_Passive(TargetName_Line);
GoodObjectIndex=ObjType==1;
BadObjectIndex=ObjType==2;
NovelObjectIndex=ObjType==3;






UniqueAngle = unique(TargetAngle_Line);
UniqueEcc = unique(TargetEcc_Line);





%%

EventChannel=Data.EventChannel(Selection,:);
EventTimeChannel=Data.EventTimeChannel(Selection,:);
EventBin=Data.EventBin(Selection,:);




%%Reproduce goodbad index from event channel
%ObjectIndex=
%ObjectIndex=ReproduceFromEvent(EventChannel,[GOODOBJ  BADOBJ]);

%GoodObjectIndex=(ObjectIndex==GOODOBJ)';


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

%sSpikeTime=squeeze(Data.SpikeTimeData(:,:,Selection))';
%SpikeBin=squeeze(Data.SpikeBinData(:,:,Selection))';

  %find mean and std for z-score
FR_Each = (sum(~isnan(SpikeTime),2)-2)./(SpikeTime(:,end)-SpikeTime(:,end-1))*1000;
FR_Mean = mean(FR_Each);
FR_Std = std(FR_Each);


%%
%%%%%%%%%%%%%%%%%%%%%%%Align data on object onset
%Segregate spike data from fractal onset to fractal offset

SegmentData=EventSegTool(SegEvent,SegEventEnd,PlotStartOffset,PlotStopOffset,EventChannel,EventTimeChannel,SpikeTime,StepSize,BinWidth);

%Retrive data from the field

%Relative time bin for each event
RelativeTimeMarker=SegmentData.RelativeTimeMarker{1};%=EventTimeMarker;
%Absolute time marker for each event
AbosoluteTimeMarker=SegmentData.AbsoluteTimeMarker{1};

ObjOnTimeMarker=nanmean(RelativeTimeMarker);

%{
RwdTimeAbs=FindOutTime(EventChannel,EventBin,RWDCD);
RwdTimeRel=RwdTimeAbs-AbosoluteTimeMarker(:,1);

%Algin marker to sac-on & sac-off
RwdTimeRelSacEnd=RwdTimeRel-SacEnd_First_Target';
TgtOffTimeRelSacEnd=RelativeTimeMarker(:,2)-SacEnd_First_Target';

RwdTimeRelSacEndMean=nanmean(RwdTimeRelSacEnd);
TgtOffTimeRelSacEndMean=nanmean(TgtOffTimeRelSacEnd);
%}
%PSTH for each event
MeanPSTH=SegmentData.MeanPSTH;%=[PSTH_Event_Time',PSTH_Event'];
MeanPSTH_Time=MeanPSTH{1};
MeanPSTH_FR_FractalOn=MeanPSTH{2};

PSTH_Time_FraOn_All=nanmean(MeanPSTH_Time,1);


usefulfig = 0;


%% Loop over the target directions
for i = 1:length(UniqueAngle )
    for j = 1:length(UniqueEcc)
        sel = TargetAngle_Line == UniqueAngle(i) & TargetEcc_Line == UniqueEcc(j);
        
        
        ConditionCounts(i,j) = sum(sel);

        if sum(sel) > 3
            usefulfig=usefulfig+1;
        end

        

        
        %Align on target on 
        curr_GoodFracOn_Loc = MeanPSTH_FR_FractalOn(sel & GoodObjectIndex,:);
        Good_FracOn_Loc_Mean{i,j} = nanmean(curr_GoodFracOn_Loc,1);
        Good_FracOn_Loc_Sem{i,j} = nanstd(curr_GoodFracOn_Loc,[],1)./sqrt(sum(~isnan(curr_GoodFracOn_Loc),1));
        
        curr_BadFracOn_Loc = MeanPSTH_FR_FractalOn(sel & BadObjectIndex,:);
        Bad_FracOn_Loc_Mean{i,j} = nanmean(curr_BadFracOn_Loc,1);
        Bad_FracOn_Loc_Sem{i,j} = nanstd(curr_BadFracOn_Loc,[],1)./sqrt(sum(~isnan(curr_BadFracOn_Loc),1));
        
        curr_NovFracOn_Loc = MeanPSTH_FR_FractalOn(sel & NovelObjectIndex,:);
        Nov_FracOn_Loc_Mean{i,j} = nanmean(curr_NovFracOn_Loc,1);
        Nov_FracOn_Loc_Sem{i,j} = nanstd(curr_NovFracOn_Loc,[],1)./sqrt(sum(~isnan(curr_NovFracOn_Loc),1));
        
        
        
        
       
        AngleEccorder{i,j} = [UniqueAngle(i),UniqueEcc(j)];


        %Average 
        %%
       
        TimeTgt = [140,340];
        Ave_FR_TgtOn_Good = nanmean(curr_GoodFracOn_Loc(:,PSTH_Time_FraOn_All >= TimeTgt(1) & PSTH_Time_FraOn_All <= TimeTgt(2)),2);
        Ave_FR_TgtOn_Bad = nanmean(curr_BadFracOn_Loc(:,PSTH_Time_FraOn_All >= TimeTgt(1) & PSTH_Time_FraOn_All <= TimeTgt(2)),2);

        

       
    
        [h,p_Tgt(i)]= ttest2(Ave_FR_TgtOn_Good,Ave_FR_TgtOn_Bad);
       

       

        TgtGood(i) = mean(Ave_FR_TgtOn_Good);
        TgtBad(i) = mean(Ave_FR_TgtOn_Bad);

        %Time sequence comparison
        % for tgt onset
        for m = 1:length(PSTH_Time_FraOn_All)
            good_tgt=curr_GoodFracOn_Loc(:,m);
            bad_tgt =curr_BadFracOn_Loc(:,m);
            [h,p_tgttime_value_loc(i,m)]=ttest2(good_tgt,bad_tgt);

            h_tgttime_value_loc(i,m) = p_tgttime_value_loc(i,m)<0.05;

            goodbaddiff = mean(good_tgt)-mean(bad_tgt);
            goodpreferindex = goodbaddiff>0;
            
            h_tgttime_value_loc_good(i,m) = goodpreferindex & h_tgttime_value_loc(i,m);
            h_tgttime_value_loc_bad(i,m) = ~goodpreferindex & h_tgttime_value_loc(i,m);

            AUC_tgttime_loc(i,m) = rocN(good_tgt,bad_tgt);
            


        end


        

        %{

        %Raster
        Raster_TgtOn_Good{i}=Raster_TgtOn(sel & GoodObjectIndex,:);
        Raster_TgtOn_Bad{i}=Raster_TgtOn(sel & ~GoodObjectIndex,:);

     
        %}


    end
end


%Use contralateral RF for export
if Monkey_curr == 1 %Robin
    RF = UniqueAngle >= 0 & UniqueAngle <= 179.99;
else
    RF = UniqueAngle >= 180 & UniqueAngle <= 360;

end

if length(UniqueAngle)==1
    RF= 1;
end

if length(UniqueAngle)>2
    nonRF = find(abs(UniqueAngle-UniqueAngle(RF))>0,1,'first');
else
    nonRF=~RF;
end


%Time for export
ExportTime = [-400,400];
SelectTime = PSTH_Time_FraOn_All>ExportTime(1) & PSTH_Time_FraOn_All<ExportTime(2);
Time_TgtOn = PSTH_Time_FraOn_All(SelectTime);


PSTH_TgtOn_Good_RF = Good_FracOn_Loc_Mean{RF,1}(SelectTime);
PSTH_TgtOn_Bad_RF = Bad_FracOn_Loc_Mean{RF,1}(SelectTime);

PSTH_TgtOn_Good_RF_z = (PSTH_TgtOn_Good_RF - FR_Mean)/FR_Std;
PSTH_TgtOn_Bad_RF_z = (PSTH_TgtOn_Bad_RF - FR_Mean)/FR_Std;




if length(UniqueAngle) == 1
    PSTH_TgtOn_Good_nonRF = NaN*Good_FracOn_Loc_Mean{1}(SelectTime);
    PSTH_TgtOn_Bad_nonRF = NaN*Bad_FracOn_Loc_Mean{1}(SelectTime);

    

    TgtGood_nonRF = NaN;
    TgtBad_nonRF = NaN;
    

    



else
    
    PSTH_TgtOn_Good_nonRF = Good_FracOn_Loc_Mean{nonRF,1}(SelectTime);
    PSTH_TgtOn_Bad_nonRF = Bad_FracOn_Loc_Mean{nonRF,1}(SelectTime);

    

    TgtGood_nonRF = TgtGood(nonRF);
    TgtBad_nonRF = TgtBad(nonRF);
    
end



PSTH_TgtOn_Good_nonRF_z = (PSTH_TgtOn_Good_nonRF - FR_Mean)/FR_Std;
PSTH_TgtOn_Bad_nonRF_z = (PSTH_TgtOn_Bad_nonRF- FR_Mean)/FR_Std;


TgtGood_RF = TgtGood(RF);
TgtBad_RF = TgtBad(RF);







p_Tgt_RF  = p_Tgt(RF);



if p_Tgt_RF < 0.05
    if TgtGood_RF > TgtBad_RF
        PSTH_TgtOn_Prefer_RF = PSTH_TgtOn_Good_RF;
        PSTH_TgtOn_nonPrefer_RF = PSTH_TgtOn_Bad_RF;

        

        PSTH_TgtOn_Prefer_nonRF = PSTH_TgtOn_Good_nonRF;
        PSTH_TgtOn_nonPrefer_nonRF = PSTH_TgtOn_Bad_nonRF;

        

    else
        PSTH_TgtOn_Prefer_RF = PSTH_TgtOn_Bad_RF;
        PSTH_TgtOn_nonPrefer_RF = PSTH_TgtOn_Good_RF;

        

        PSTH_TgtOn_Prefer_nonRF = PSTH_TgtOn_Bad_nonRF;
        PSTH_TgtOn_nonPrefer_nonRF = PSTH_TgtOn_Good_nonRF;

        
    end

    
else
    if TgtGood_RF > TgtBad_RF
        PSTH_TgtOn_Prefer_RF = PSTH_TgtOn_Good_RF;
        PSTH_TgtOn_nonPrefer_RF = PSTH_TgtOn_Bad_RF;

        
        PSTH_TgtOn_Prefer_nonRF = PSTH_TgtOn_Good_nonRF;
        PSTH_TgtOn_nonPrefer_nonRF = PSTH_TgtOn_Bad_nonRF;

        

    else
        PSTH_TgtOn_Prefer_RF = PSTH_TgtOn_Bad_RF;
        PSTH_TgtOn_nonPrefer_RF = PSTH_TgtOn_Good_RF;

        

        PSTH_TgtOn_Prefer_nonRF = PSTH_TgtOn_Bad_nonRF;
        PSTH_TgtOn_nonPrefer_nonRF = PSTH_TgtOn_Good_nonRF;

        
    end

end




PSTH_TgtOn_Prefer_RF_z = (PSTH_TgtOn_Prefer_RF - FR_Mean)/FR_Std;

PSTH_TgtOn_nonPrefer_RF_z = (PSTH_TgtOn_nonPrefer_RF - FR_Mean)/FR_Std;


PSTH_TgtOn_Prefer_nonRF_z = (PSTH_TgtOn_Prefer_nonRF - FR_Mean)/FR_Std;
PSTH_TgtOn_nonPrefer_nonRF_z = (PSTH_TgtOn_nonPrefer_nonRF - FR_Mean)/FR_Std;



hPSTH_TgtOn_Value_RF = h_tgttime_value_loc(RF,SelectTime);

hPSTH_TgtOn_Good_RF = h_tgttime_value_loc_good(RF,SelectTime);
hPSTH_TgtOn_Bad_RF = h_tgttime_value_loc_bad(RF,SelectTime);

AUCPSTH_TgtOn_Value_RF = AUC_tgttime_loc(RF,SelectTime);

if length(UniqueAngle)==1
    hPSTH_TgtOn_Value_nonRF = NaN*h_tgttime_value_loc(RF,SelectTime);
    hPSTH_TgtOn_Good_nonRF= NaN*h_tgttime_value_loc(RF,SelectTime);
    hPSTH_TgtOn_Bad_nonRF = NaN*h_tgttime_value_loc(RF,SelectTime);

    AUCPSTH_TgtOn_Value_nonRF = NaN*h_tgttime_value_loc(RF,SelectTime);
else
    hPSTH_TgtOn_Value_nonRF = h_tgttime_value_loc(nonRF,SelectTime);
    hPSTH_TgtOn_Good_nonRF = h_tgttime_value_loc_good(nonRF,SelectTime);
    hPSTH_TgtOn_Bad_nonRF = h_tgttime_value_loc_bad(nonRF,SelectTime);

    AUCPSTH_TgtOn_Value_nonRF = AUC_tgttime_loc(nonRF,SelectTime);
end






%%
if ShowFigureFlag
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Figure 1: Average firing rate
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%for each location


figtitlestr{1}='PSTHAlignOnFracOnset_EachLocation';

fig{1}=PrepareFigure(110,'w',[50,100, 1200,800],'Name',figtitlestr{1});


%Set a layout for each angle 
if length(UniqueAngle) <= 2 | usefulfig==2
  angle = [89,270];
  subloc = [1,3];  
  layout = [1,3];
  center = 1;
    
else
angle = 0:45:359;
subloc = [2,1,4,7,8,9,6,3];
if length(UniqueEcc) == 1 || (length(UniqueEcc) == 2 && ismember(0,UniqueEcc))
    layout = [3,3];
    center = 5;
else
    totalpart = (length(UniqueEcc) - ismember(0,UniqueEcc))*length(angle)+ismember(0,UniqueEcc);
    layout = [ceil(sqrt(totalpart)),ceil(sqrt(totalpart))];
    center = (layout(1)+1)*2+1;
    
end
end

for i = 1:length(UniqueAngle )
    for j = 1:length(UniqueEcc)
        if UniqueEcc(j) == 0
            index = center;
        else
            [minValue,closestIndex] = min(abs(UniqueAngle(i)-angle));
            index = subloc(closestIndex) ;
        end
        subplot(layout(1),layout(2),index);

        if ConditionCounts(i,j)>3
        
        shadedErrorBar(PSTH_Time_FraOn_All,Good_FracOn_Loc_Mean{i,j},Good_FracOn_Loc_Sem{i,j},'lineprops',{[1    0    0]},'transparent',1,'patchSaturation',0.3);

        hold on
        shadedErrorBar(PSTH_Time_FraOn_All,Bad_FracOn_Loc_Mean{i,j},Bad_FracOn_Loc_Sem{i,j},'lineprops',{[0    0    1]},'transparent',1,'patchSaturation',0.3);
        shadedErrorBar(PSTH_Time_FraOn_All,Nov_FracOn_Loc_Mean{i,j},Nov_FracOn_Loc_Sem{i,j},'lineprops',{[0    1    0]},'transparent',1,'patchSaturation',0.3);
        end

       % yrange=ylim;
       % yrange(2)=max(MeanPSTH_FR_FractalOn_mean(PSTH_Time_FraOn_All<PlotSacPostBuffer & PSTH_Time_FraOn_All>-PlotSacPreBuffer))*1.05;
          xlim([-PlotSacPreBuffer, PlotSacPostBuffer]);%Keep in the same scale
       set(findall(gca, 'Type', 'Line'),'LineWidth',3);

        set(gca,'FontSize',15,'FontWeight','Bold');
        set(gca,'LineWidth',3);
        title(sprintf("Angle: %1.1d,Ecc: %1.1d",AngleEccorder{i,j}));
        if index == center + layout(1)
            xlabel('Time from Fractal Onset(ms)','FontSize',15,'FontWeight','Bold');
            
        end
        
        ylabel('Firing Rate(Hz)','FontSize',15,'FontWeight','Bold');
        

       
        
        
    end
end




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
ChannelNum=SpikeChannel;


%To get the monkey name
RemainFiles=erase(FileName,RecordDateOriginal);
Delimiter=find(isstrprop(RemainFiles,'upper')==1);
MonkeyName=RemainFiles(Delimiter(1):Delimiter(2)-1);


%Output file name
OutputFileName=sprintf('%s_%s_N%s_C%s.mat',MonkeyName,string(RecordDateOriginal),string(NeuronNum),string(ChannelNum));
OutputFigureName=sprintf('%s_%s_N%s_C%s',MonkeyName,string(RecordDateOriginal),string(NeuronNum),string(ChannelNum));


ExistFlag=0;
%Load the old file if exist
if exist(OutputFileName)

    load(OutputFileName);
    %Load OutputData into the memory
    if isfield(OutputData,'ObjectPassiveViewing')
    ExistFlag=1;
    end
end


%%%%%%%%%%%%%%%%%%%%OutputData For current analysis
OutputData_New=[];

OutputData_New.ObjectPassiveViewing.TrialType='GoodTrials';
%%%%%%%%%Save overall scene response characteristics%%%%%%%%%%%%%%%%%%%%%%%
%Organize the data



%Organize the data
NameStr={'PSTH_TgtOn_Good_RF','PSTH_TgtOn_Bad_RF','PSTH_TgtOn_Good_nonRF','PSTH_TgtOn_Bad_nonRF',...
    'PSTH_TgtOn_Good_RF_z','PSTH_TgtOn_Bad_RF_z',...
    'PSTH_TgtOn_Good_nonRF_z','PSTH_TgtOn_Bad_nonRF_z',...
    'TgtGood_RF','TgtBad_RF','TgtGood_nonRF','TgtBad_nonRF',...
    'PSTH_TgtOn_Prefer_RF','PSTH_TgtOn_nonPrefer_RF',...
    'PSTH_TgtOn_Prefer_nonRF','PSTH_TgtOn_nonPrefer_nonRF',...
    'Time_TgtOn','p_Tgt_RF',...
    'PSTH_TgtOn_Prefer_RF_z','PSTH_TgtOn_nonPrefer_RF_z',...
    'PSTH_TgtOn_Prefer_nonRF_z','PSTH_TgtOn_nonPrefer_nonRF_z',...
     'hPSTH_TgtOn_Value_RF','hPSTH_TgtOn_Good_RF','hPSTH_TgtOn_Bad_RF',...
    'hPSTH_TgtOn_Value_nonRF','hPSTH_TgtOn_Good_nonRF','hPSTH_TgtOn_Bad_nonRF',...
    'AUCPSTH_TgtOn_Value_RF','AUCPSTH_TgtOn_Value_nonRF',...
    'UniqueAngle'

 
};

DataLib={PSTH_TgtOn_Good_RF,PSTH_TgtOn_Bad_RF,PSTH_TgtOn_Good_nonRF,PSTH_TgtOn_Bad_nonRF,...
    PSTH_TgtOn_Good_RF_z,PSTH_TgtOn_Bad_RF_z,...
    PSTH_TgtOn_Good_nonRF_z,PSTH_TgtOn_Bad_nonRF_z,...
    TgtGood_RF,TgtBad_RF,TgtGood_nonRF,TgtBad_nonRF,...
    PSTH_TgtOn_Prefer_RF,PSTH_TgtOn_nonPrefer_RF,...
    PSTH_TgtOn_Prefer_nonRF,PSTH_TgtOn_nonPrefer_nonRF,...
    Time_TgtOn,p_Tgt_RF,...
    PSTH_TgtOn_Prefer_RF_z,PSTH_TgtOn_nonPrefer_RF_z,...
    PSTH_TgtOn_Prefer_nonRF_z,PSTH_TgtOn_nonPrefer_nonRF_z,...
    hPSTH_TgtOn_Value_RF,hPSTH_TgtOn_Good_RF,hPSTH_TgtOn_Bad_RF,...
    hPSTH_TgtOn_Value_nonRF,hPSTH_TgtOn_Good_nonRF,hPSTH_TgtOn_Bad_nonRF,...
    AUCPSTH_TgtOn_Value_RF,AUCPSTH_TgtOn_Value_nonRF,...
    UniqueAngle
};



DataStamp=containers.Map(NameStr,DataLib);





OutputData_New.ObjectPassiveViewing.Task=TaskType;
OutputData_New.ObjectPassiveViewing.TaskCode=TaskCode;
OutputData_New.ObjectPassiveViewing.DataStamp=DataStamp;

if ExistFlag
%Compare the old one with the new one
Task_Old=OutputData.ObjectPassiveViewing.TaskCode;

     if Task_Old~=TaskCode
        %Not the same task, add the new dataset into the old one
        OutputData.ObjectPassiveViewing(numel(OutputData.ObjectPassiveViewing)+1)=OutputData_New.ObjectPassiveViewingg;
     else
         %If the same task,the same protocol,replace the old one with the new one
         
         OutputData.ObjectPassiveViewing=OutputData_New.ObjectPassiveViewing;
         
    
    
     end
else
    %Output the current file
   OutputData.ObjectPassiveViewing=OutputData_New.ObjectPassiveViewing; 
end



if OutputFlag
    cd(OutPath);
    save(OutputFileName,'OutputData');
    disp('Data Exported');
else
    disp('No data export');
end  


%{
%Store figures in common folders

ProtocolMarker='PassiveView';
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
end %End of the spikechannel



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

