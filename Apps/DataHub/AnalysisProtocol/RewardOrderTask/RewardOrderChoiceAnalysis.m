function RewardOrderChoiceAnalysis(Data, StartTrial, EndTrial,ShowFigureFlag,OutputFlag);
ECode;
TaskType=Data.TaskType;
DataPath= Data.Path;
FileName=Data.FileName;

Selection=StartTrial:EndTrial;%Select the trials according to the pannel

%Criterias for sexwlecting saccades
V_Threshold=40;
ContinueBin=5;
ISI_Threshold=20;

PlotTaskInfo = 1;

%{
%Abnormal check %Especially For Robin081120FORAGEXPLO4.rst
CheckCode=SCENE_ON;%This code should appear once in each trial
AbnormalTrials=AbnormalTrialCheck(Data,CheckCode);
Selection=Selection(~ismember(Selection,AbnormalTrials));
TrialNum=length(Selection);
%}
%Extra information : Task Type,  Reward Order
ParaLib=Data.ParaLib;
if(ismember('TaskType',keys(ParaLib)) && (ismember('CurrentRewardOrder',keys(ParaLib))))
    TaskType_Info=ParaLib('TaskType');
    TaskType=TaskType_Info(Selection,end);
    
    RewardOrderInfo = ParaLib('CurrentRewardOrder');
    RewardOrder = RewardOrderInfo(Selection,end-1:end);
    
    BlockModeInfo = ParaLib('CurrentBlockMode');
    BlockMode = BlockModeInfo(Selection, end);
    
    
else
    disp("Not enough saved paramerer!");
    return 
    
end

%Select all choice trials
ChoiceTrials = TaskType == 2;
%Selection = Selection(ChoiceTrials);


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


%Load parameters
try
    TargetName=Data.TargetIndex(Selection,2:end,2);
    TargetIndex=Data.TargetIndex(Selection,2:end,1);
catch
    TargetName=Data.TargetIndex(Selection,3);
    TargetIndex=Data.TargetIndex(Selection,2);
    
end
TargetLoc=Data.TargetLoc(Selection,2:end,:);

%TargetNum=Data.TargetNum(Selection,2:end);
TargetNum=sum(~isnan(TargetName),2);
%TargetSize=Data.TargetSize(Selection,2:end);
try
    TargetPolarAngle=Data.TargetPolarLoc(Selection,2:end,1);

    TargetPolarEccentricity=Data.TargetPolarLoc(Selection,2:end,2);
catch
    TargetPolarAngle=Data.TargetPolarLoc(Selection,2);

    TargetPolarEccentricity=Data.TargetPolarLoc(Selection,3);
    
end

RewardOrder = RewardOrder(Selection,:);
BlockMode =BlockMode(Selection,:);
TaskType =TaskType(Selection,:);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Analysis the
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%choices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ObjectChosen_tmp=ReproduceFromEvent(EventChannel,[GGFRAC : DBFRAC])';
ObjectChosen = cellfun(@(x) x(1),ObjectChosen_tmp);

NumOfCondition = 6;
ChoicesMat = NaN*ones(NumOfCondition,NumOfCondition);
NumTrialsCon =  ChoicesMat;

RewardOrderChoice = RewardOrder(ChoiceTrials,:);
ObjectChosenChoice = ObjectChosen(ChoiceTrials,:);


%%Calculate the total time for each trial
StartTrialTime = FindOutTime(EventChannel,EventBin,SHOWFIXCD);
EndTrialTime = FindOutTime(EventChannel,EventBin,TGTOFF);
TotalDuration = EndTrialTime - StartTrialTime ;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Analyze the reaction time for both single objects and choice trials
%Select the period between fixation offset target reshow
offsetEye = 500;
FixOffsetBin = FindOutTime(EventChannel,EventBin,FIXOFF	)-offsetEye;%-100: offset
EndSectionBin= FindOutTime(EventChannel,EventBin,CHOICECONFIRM	);



%FractalOnsetTime = FindOutTime(EventChannel,EventTimeChannel,TGTCD );
%FractalOffsetTime = FindOutTime(EventChannel,EventTimeChannel,TGTOFF );
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


EyeChannel_X_Selection_All= SelectEyeInterval(EyeChannel_X_L,FixOffsetBin,EndSectionBin)';
EyeChannel_Y_Selection_All= SelectEyeInterval(EyeChannel_Y_L,FixOffsetBin ,EndSectionBin)';

EyeChannel_X_Selection=ReorganizeEye(EyeChannel_X_Selection_All{1});
EyeChannel_Y_Selection=ReorganizeEye(EyeChannel_Y_Selection_All{1});



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

AmplitudeThreshold=8;

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


RT_First_Target=RT_First_Target-offsetEye;

RT_First_Target=RT_First_Target';

disp("Number of trials no saccades were found: " + sum(isnan(RT_First_Target)));

%First for choice trials 
ChoiceRTMat_Mean = NaN*ones(NumOfCondition,NumOfCondition);
ChoiceRTMat_OtherMean = NaN*ones(NumOfCondition,NumOfCondition);
RTFirstChoice = RT_First_Target(ChoiceTrials);

%Then for sigle target trials
RTFirstSingleTarget = RT_First_Target(~ChoiceTrials);
RewardOrderSingle = RewardOrder(~ChoiceTrials,1);
TotalDurationSingle = TotalDuration(~ChoiceTrials);

%% Loop over all conditions for all choices
if sum(ChoiceTrials)>0
    %Looping over choice options
    for i = GGFRAC:DBFRAC
        for j = GGFRAC:DBFRAC
            if i > j
                %Choices 
                TrialsFlag = (RewardOrderChoice(:,1) == i & RewardOrderChoice(:,2) == j) | (RewardOrderChoice(:,1) == j & RewardOrderChoice(:,2) == i);
                NumTrialsCon(i-GGFRAC+1,j-GGFRAC+1) = sum(TrialsFlag);
                Chosen = ObjectChosenChoice(TrialsFlag);
              
                ChoicesMat(i-GGFRAC+1,j-GGFRAC+1) = sum(Chosen == j);
                
                
                %RT
                RT_temp = RTFirstChoice(TrialsFlag);
                RT_Choice_tmp = RT_temp(Chosen==j);
                RT_ChoiceOther_tmp = RT_temp(Chosen==i);
                
                ChoiceRTMat_Mean(i-GGFRAC+1,j-GGFRAC+1) = nanmean(RT_Choice_tmp);
                ChoiceRTMat_OtherMean(i-GGFRAC+1,j-GGFRAC+1) = nanmean(RT_ChoiceOther_tmp);
            
            end
        
        end
    end


end
ProbTopChosen = ChoicesMat./NumTrialsCon *100;

RT_SingleMean = NaN*ones(1,NumOfCondition);
RT_SingleSem = NaN*ones(1,NumOfCondition);
NumTrialSingle = NaN*ones(1,NumOfCondition);
Dur_MeanS = NaN*ones(1,NumOfCondition);
Dur_SemS = NaN*ones(1,NumOfCondition);

if(sum(~ChoiceTrials)>0)
    %For single target trials
    UniqueROSingle = unique(RewardOrderSingle);  
    
    RO = {"GG",'GB','BG','BB','2G0','2B0'};
    ROStandard = [GGFRAC:DBFRAC];
    Order = 1:NumOfCondition;
    
    for i = 1 : length(UniqueROSingle)
        select = RewardOrderSingle == UniqueROSingle(i);
        RT_tmp = RTFirstSingleTarget(select);
        
        RT_SingleMean(i) =nanmean(RT_tmp);
        RT_SingleSem(i) =nanstd(RT_tmp)/sqrt(sum(~isnan(RT_tmp)));
        
        Order(i) = find(ROStandard == UniqueROSingle(i));
        NumTrialSingle(i) = sum(select);
        
        currDur = TotalDurationSingle(select);
        Dur_MeanS(i) = nanmean(currDur);
        Dur_SemS(i) = nanstd(currDur)/sqrt(sum(select));
        
        
    end 
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Visualization%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ShowFigureFlag
    FigureStartNum=130;
    FigureIndex=1;
    
  if sum(ChoiceTrials)>0  
     
     %Choices 
    figtitlestr{FigureIndex}='ChoicesMatrix';
    fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[50,100, 1200,400],'Name',figtitlestr{FigureIndex});
    %Choices Matrix
    subplot(1,3,1)
    h=heatmap(ProbTopChosen);
    colormap default;
    h.XDisplayLabels = RO;
    h.YDisplayLabels = RO;
    h.XLabel = "Chosen Target(%Proportion)";
    h.YLabel = "Distractor Target";
    h.Title = "ChoicesMatrix";
    
    set(gca,'FontSize',15);
    
    %TrialNumber Matrix
    subplot(1,3,2)
    h=heatmap(NumTrialsCon);
    colormap default;
    h.XDisplayLabels = RO;
    h.YDisplayLabels = RO;
    h.XLabel = "Chosen Target";
    h.YLabel = "Distractor Target";
    h.Title = "TrialNumber";
    
    set(gca,'FontSize',15);
    
    %Choices Number Matrix
    subplot(1,3,3)
    h=heatmap(ChoicesMat);
    colormap default;
    h.XDisplayLabels = RO;
    h.YDisplayLabels = RO;
    h.XLabel = "Chosen Target";
    h.YLabel = "Distractor Target";
    h.Title = "RawChoicesNumber";
    
    set(gca,'FontSize',15);
    
    FigureStartNum = FigureStartNum++1;
    FigureIndex = FigureIndex+ 1;
    
    %Reaction time
    figtitlestr{FigureIndex}='RTMatrix';
    fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[70,100, 800,400],'Name',figtitlestr{FigureIndex});
    %Choose X
    subplot(1,2,1)
    h=heatmap(ChoiceRTMat_Mean);
    colormap default;
    h.XDisplayLabels = RO;
    h.YDisplayLabels = RO;
    h.XLabel = "Chosen Target";
    h.YLabel = "Distractor Target";
    h.Title = "RTOfChosen";
    
    set(gca,'FontSize',15);
    
    %TrialNumber Matrix
    subplot(1,2,2)
    h=heatmap(ChoiceRTMat_OtherMean);
    colormap default;
    h.XDisplayLabels = RO;
    h.YDisplayLabels = RO;
    h.XLabel = "Chosen Target";
    h.YLabel = "Distractor Target";
    h.Title = "RTOfOther";
    
    set(gca,'FontSize',15);
    
    
    
    FigureStartNum = FigureStartNum++1;
    FigureIndex = FigureIndex+ 1;
    
    
    
  end %End of the if choicetrials > 0   
  
  if sum(~ChoiceTrials)>0  
      %Single target
      figtitlestr{FigureIndex}='RTForSingleTarget';
      fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[100,100, 400,200],'Name',figtitlestr{FigureIndex});
      b = bar(RT_SingleMean);
      hold on
      errorbar([1:6],RT_SingleMean,RT_SingleSem,'ok');
      ylabel("RT(ms)");
      set(gca,'xTickLabel',RO(Order),'FontSize',15,'FontWeight','Bold','LineWidth',3);
      box off
      
      FigureStartNum = FigureStartNum++1;
      FigureIndex = FigureIndex+ 1;
    
      
  end
  
  if PlotTaskInfo 
      if sum(~ChoiceTrials)>0  
        figtitlestr{FigureIndex}='TaskInfoForSingleTarget';
        fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[100,120, 400,200],'Name',figtitlestr{FigureIndex});
        subplot(1,2,1)
        b = bar(NumTrialSingle);    
        ylabel("Num of trials");
        set(gca,'xTickLabel',RO(Order),'FontSize',15,'FontWeight','Bold','LineWidth',3);
        box off
        
        subplot(1,2,2)
        b = bar(Dur_MeanS); 
        hold on
        errorbar([1:6],Dur_MeanS,Dur_SemS,'ok');
        
        ylabel("Task Duration(ms)");
        set(gca,'xTickLabel',RO(Order),'FontSize',15,'FontWeight','Bold','LineWidth',3);
        box off
        
      
        FigureStartNum = FigureStartNum++1;
        FigureIndex = FigureIndex+ 1;
      end
      
      
      
  end%End of if plot task info
  
  
  
  
  
end
    
end




