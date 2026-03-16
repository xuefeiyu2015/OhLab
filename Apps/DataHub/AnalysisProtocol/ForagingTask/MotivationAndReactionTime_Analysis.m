function MotivationAndReactionTime_Analysis(Data, StartTrial, EndTrial,ShowFigureFlag,OutputFlag)
%Analysis the monkeys' motivation to start a task, and reaction time to
%select target
%Xuefei Yu 02232020

%Unpack all the data base and task settings
ECode;
FileName=Data.FileName;
TaskType=Data.TaskType;
TaskCode=Data.TaskCode;
TrialNum=EndTrial-StartTrial+1;
Selection=StartTrial:EndTrial;
TargetIndex=Data.TargetIndex(Selection,:,:);

BlockType=Data.BlockType(Selection,:,:);%0:fixed sequence;1;random sequence; 0:ScenesetA;1:SceneSetB
BlockType(BlockType(:,3)>2,3)=999;
ChoiceOutCome=Data.ChoiceOutcome(Selection,:,:);

SceneIndex=Data.SceneIndex(Selection,:);
EyeChannel_X_L= Data.EyeDataRaw(Selection,1);
EyeChannel_Y_L= Data.EyeDataRaw(Selection,2);




Eye_Eccentricity=cellfun(@(x,y) sqrt(x.^2+y.^2),EyeChannel_X_L,EyeChannel_Y_L,'uniform',0);


EyeBinWidth=Data.EyeBinWidth;
%EyeTime=Data.EyeTime(Selection,:);

EventChannel=Data.EventChannel(Selection,:);
EventTimeChannel=Data.EventTimeChannel(Selection,:);
EventBin=Data.EventBin(Selection,:);
%{
%Pre-norm the eye position
FixationOnBin=FindOutTime(EventChannel,EventBin,TGTCD);
FixationOnBin=FixationOnBin(:,1);
%}

%Color Definition for each scenes

Color=[lines(6);prism(6);pink(6)];
Color(2,:)=copper(1)+0.5;


Reward=Data.RewardAmount;
SceneReward=[SceneIndex(:,2),Reward(:,2)];
%Remove the zero reward

SceneReward(SceneReward(:,2)==0,:)=[];
SceneRewardAssociation=unique(SceneReward,'rows');

%First search for the extra fixation point
ExtraFPOnset=FindOutTime(EventChannel,EventTimeChannel,EXTRAFP);
ExtraEyeIn=FindOutTime(EventChannel,EventTimeChannel,EXTRAFIXSTARTCD);

RT_ExtraFix=ExtraEyeIn-ExtraFPOnset;

%Then search for the first fixation point onset

%Reaction time for the fixation onset
FPOnset=FindOutTime(EventChannel,EventTimeChannel,SHOWFIXCD );

RealEyeIn=FindOutTime(EventChannel,EventTimeChannel,FIXSTARTCD );
try
RT_Fix=RealEyeIn(:,2:end)-FPOnset;
catch
RT_Fix=RealEyeIn-FPOnset;
end
RT_Fix_First=RT_Fix(:,1);


%Next is the reaction time for the chosen target
%Chunk the eye data

FractalOnsetBin = FindOutTime(EventChannel,EventBin,TGTCD );
FractalOffsetBin = FindOutTime(EventChannel,EventBin,TGTOFF);

AlignTime= FindOutTime(EventChannel,EventTimeChannel,TGTCD);

%Only look at the chosen one
temp=~isnan(FractalOnsetBin);
Indices=arrayfun(@(x) find(temp( x,:), 1, 'last'), 1:size(FractalOnsetBin, 1));
ChosenFractalOnset=arrayfun(@(x,y) FractalOnsetBin(x,y),1:size(FractalOnsetBin, 1), Indices)';
ChosenFractalOffset=arrayfun(@(x,y) FractalOffsetBin(x,y),1:size(FractalOffsetBin, 1), Indices)'; 
ChosenAlignTime=arrayfun(@(x,y) AlignTime(x,y),1:size(AlignTime, 1), Indices)'; 


EyeX_in_bin=SelectEyeInterval(EyeChannel_X_L,ChosenFractalOnset ,ChosenFractalOffset);
EyeY_in_bin=SelectEyeInterval(EyeChannel_Y_L,ChosenFractalOnset ,ChosenFractalOffset);
EyeX_in_bin=EyeX_in_bin{1};
EyeY_in_bin=EyeY_in_bin{1};

EyeTime_in_Bin_Aligned=cellfun(@(x) [0:numel(x)],EyeX_in_bin,'uniform',0);




%Criterias for selecting saccades
V_Threshold=40;
ContinueBin=5;
ISI_Threshold=20;

Saccades=SelectSaccadeFromCell(EyeX_in_bin,EyeY_in_bin,EyeBinWidth,V_Threshold,ContinueBin,ISI_Threshold);

SaccadeNumber=arrayfun(@(x)  x.NumOfSaccade,Saccades)';
SaccadeAngle=arrayfun(@(x)  x.SaccadeAngle,Saccades,'UniformOutput' ,0)';
SaccadeAmplitude=arrayfun(@(x)  x.SaccadeAmplitude,Saccades,'UniformOutput' ,0)';

SaccadeStartTime=arrayfun(@(x)  x.SaccadeStartTime,Saccades,'UniformOutput' ,0)';
SaccadeEndTime=arrayfun(@(x)  x.SaccadeEndTime,Saccades,'UniformOutput' ,0)';

SaccadeStartPoint=arrayfun(@(x)  x.SaccadeStartPoint,Saccades,'UniformOutput' ,0)';
SaccadeEndPoint=arrayfun(@(x)  x.SaccadeEndPoint,Saccades,'UniformOutput' ,0)';

%Screen out the valid saccades and saccade start time above 50;

AmplitudeThreshold=5;

TimeThreshold=50;

SaccadeAngle_Selection=cellfun(@(x,y,z) z(((x>AmplitudeThreshold)&(y>TimeThreshold))),SaccadeAmplitude,SaccadeStartTime,SaccadeAngle,'UniformOutput' ,0);
SaccadeEndPoint_Selection=cellfun(@(x,y,z) z(((x>AmplitudeThreshold)&(y>TimeThreshold)),:),SaccadeAmplitude,SaccadeStartTime,SaccadeEndPoint,'UniformOutput' ,0);
SaccadeStartPoint_Selection=cellfun(@(x,y,z) z(((x>AmplitudeThreshold)&(y>TimeThreshold)),:),SaccadeAmplitude,SaccadeStartTime,SaccadeStartPoint,'UniformOutput' ,0);

SaccadeStartTime_Selection=cellfun(@(x,y,z) z(((x>AmplitudeThreshold)&(y>TimeThreshold)),:),SaccadeAmplitude,SaccadeStartTime,SaccadeStartTime,'UniformOutput' ,0);
%Only look at the first saccade
SaccadeNumInEachTrial=cellfun(@numel,SaccadeAngle_Selection);


SaccadeAngle_First=NaN*ones(1,length(SaccadeNumInEachTrial));
SaccadeEndPoint_First=NaN*ones(length(SaccadeNumInEachTrial),2);
RT_First_Target=NaN*ones(1,length(SaccadeNumInEachTrial));
SaccadeStartPoint_First=NaN*ones(length(SaccadeNumInEachTrial),2);

SaccadeAngle_First(SaccadeNumInEachTrial>0)=cell2mat(cellfun(@(x) x(1,:),SaccadeAngle_Selection(SaccadeNumInEachTrial>0),'UniformOutput' ,0));
SaccadeEndPoint_First(SaccadeNumInEachTrial>0,:)=cell2mat(cellfun(@(x) x(1,:),SaccadeEndPoint_Selection(SaccadeNumInEachTrial>0,:),'UniformOutput' ,0));
RT_First_Target(SaccadeNumInEachTrial>0)=cell2mat(cellfun(@(x) x(1,:),SaccadeStartTime_Selection(SaccadeNumInEachTrial>0,:),'UniformOutput' ,0));
SaccadeStartPoint_First(SaccadeNumInEachTrial>0,:)=cell2mat(cellfun(@(x) x(1,:),SaccadeStartPoint_Selection(SaccadeNumInEachTrial>0,:),'UniformOutput' ,0));

%}

UniqueBlockType=unique(BlockType(:,2:end),'row');
NumberOfBlock=size(UniqueBlockType,1);
%UniqueScene=unique(SceneIndex(:,2));

if ishandle(122)
close 122;end
figure(122);set(122,'color','w');
set(122,'Position', [100,100, 1000,500], 'Name', 'Before task motivation measurement for each scene');


if ishandle(123)
close 123;end
figure(123);set(123,'color','w');
set(123,'Position', [100,100, 1000,500], 'Name', 'RT to first FP for each scene');

if ishandle(124)
close 124;end
figure(124);set(124,'color','w');
set(124,'Position', [100,100, 1000,500], 'Name', 'RT to Chosen Target for each scene');



for i=1:NumberOfBlock
    currBlockType=UniqueBlockType(i,:);
    SceneIndexForCurrentBlock=SceneIndex(sum(BlockType(:,2:end)==currBlockType,2)==2,2);
    UniqueScene=unique(SceneIndexForCurrentBlock);
    clear RT_ExtraFix_CurrScene RT_Fix_CurrScene RT_First_Target_CurrScene; 
    clear RT_ExtraFix_CurrScene_Mean RT_Fix_CurrScene_Mean RT_First_Target_CurrScene_Mean;
    clear RT_ExtraFix_CurrScene_Sem RT_Fix_CurrScene_Sem RT_First_Target_CurrScene_Sem;
 
    for j=1:length(UniqueScene)
        currSceneIndex=UniqueScene(j);
        Selection=(sum(BlockType(:,2:end)==currBlockType,2)==2)&(SceneIndex(:,2)==currSceneIndex);
        
        RT_ExtraFix_CurrScene{j}=RT_ExtraFix(Selection);
        RT_Fix_CurrScene{j}=RT_Fix_First(Selection);
        RT_First_Target_CurrScene{j}=RT_First_Target(Selection);
        
        RT_ExtraFix_CurrScene_Mean(j)=nanmean(RT_ExtraFix(Selection));
        RT_Fix_CurrScene_Mean(j)=nanmean(RT_Fix_First(Selection));
        RT_First_Target_CurrScene_Mean(j)=nanmean(RT_First_Target(Selection));
        
        RT_ExtraFix_CurrScene_Sem(j)=nanstd(RT_ExtraFix(Selection))/sqrt(sum(Selection));
        RT_Fix_CurrScene_Sem(j)=nanstd(RT_Fix_First(Selection))/sqrt(sum(Selection));
        RT_First_Target_CurrScene_Sem(j)=nanstd(RT_First_Target(Selection))/sqrt(sum(Selection));
  
    end
    figure(122);
    subplot(NumberOfBlock,1,i);
    b=bar(UniqueScene,RT_ExtraFix_CurrScene_Mean);
     b.FaceColor = 'flat';
     b.CData = Color(1:length(UniqueScene),:);
  if currBlockType(1)==0
    title('Fixed Sequence Condition');
  elseif currBlockType(1)==-1
    title('Reverse Fixed Sequence Condition');
  else
      title('Random Sequence Condition');
  end
    hold on 
    errorbar(UniqueScene,RT_ExtraFix_CurrScene_Mean,RT_ExtraFix_CurrScene_Sem,'.k')
   ylabel('Reaction time to reach extra FP(ms)')
   box off
    
   figure(123);
    subplot(NumberOfBlock,1,i);
    b=bar(UniqueScene,RT_Fix_CurrScene_Mean);
     b.FaceColor = 'flat';
     b.CData = Color(1:length(UniqueScene),:);
   % b.CData = [1 0 0;0,1,0;0,0,1];
  % xticklabels({'Poor','Normal','Rich'})
  if currBlockType(1)==0
    title('Fixed Sequence Condition');
    elseif currBlockType(1)==-1
    title('Reverse Fixed Sequence Condition');
    
  else
    title('Random Sequence Condition');
  end
    hold on 
    errorbar(UniqueScene,RT_Fix_CurrScene_Mean,RT_Fix_CurrScene_Sem,'.k')
   ylabel('Reaction time to reach First FP(ms)')
   box off
   
    figure(124);
    subplot(NumberOfBlock,1,i);
    b=bar(UniqueScene,RT_First_Target_CurrScene_Mean);
     b.FaceColor = 'flat';
     b.CData = Color(1:length(UniqueScene),:);
   % b.CData = [1 0 0;0,1,0;0,0,1];
  % xticklabels({'Poor','Normal','Rich'})
  if currBlockType(1)==0
    title('Fixed Sequence Condition');
    elseif currBlockType(1)==-1
    title('Reverse Fixed Sequence Condition');
  else
    title('Random Sequence Condition');
  end
    hold on 
    errorbar(UniqueScene,RT_First_Target_CurrScene_Mean,RT_First_Target_CurrScene_Sem,'.k')
   ylabel('Reaction time to reach First Target(ms)')
   box off
   
   
    
end

if TaskCode==FORAGING_TASK_SEQ 
    
PoorSceneIndex=[THIRD_SCENE,SIXTH_SCENE];
    NormalSceneIndex=[FIRST_SCENE,FIFTH_SCENE];
    RichSceneIndex=[SECOND_SCENE,FOURTH_SCENE];
 FixedSeq=sum(BlockType(:,2:end)==[0,0],2)==2;%This analysis only for the rich-normal-poor scene trials
    RandomSeq=sum(BlockType(:,2:end)==[1,0],2)==2;

if sum(~isnan(ExtraFPOnset))==0
    disp('No extra fixation point in this block!')
else
    if sum(isnan(ExtraFPOnset))~=sum(isnan(ExtraEyeIn))
        disp('This might be not a forced fixation block')
        
        
    end
  
    
    
    
   
    
    
    ExtraFPRTForRichScene_Fixed=RT_ExtraFix((ismember(SceneIndex(:,2),RichSceneIndex))&FixedSeq);
    ExtraFPRTForPoorScene_Fixed=RT_ExtraFix((ismember(SceneIndex(:,2),PoorSceneIndex))&FixedSeq);
    ExtraFPRTForNormalScene_Fixed=RT_ExtraFix((ismember(SceneIndex(:,2),NormalSceneIndex))&FixedSeq);
    
    
    ExtraFPRTForRichScene_Random=RT_ExtraFix((ismember(SceneIndex(:,2),RichSceneIndex))&RandomSeq);
    ExtraFPRTForPoorScene_Random=RT_ExtraFix((ismember(SceneIndex(:,2),PoorSceneIndex))&RandomSeq);
    ExtraFPRTForNormalScene_Random=RT_ExtraFix((ismember(SceneIndex(:,2),NormalSceneIndex))&RandomSeq);
    
    
    
    Mean_ExtraFP_Rich_Fixed=nanmean(ExtraFPRTForRichScene_Fixed);
    Mean_ExtraFP_Poor_Fixed=nanmean(ExtraFPRTForPoorScene_Fixed);
    Mean_ExtraFP_Normal_Fixed=nanmean(ExtraFPRTForNormalScene_Fixed);
    
    SEM_ExtraFP_Rich_Fixed=nanstd(ExtraFPRTForRichScene_Fixed)/sqrt(sum(~isnan(ExtraFPRTForRichScene_Fixed)));
    SEM_ExtraFP_Poor_Fixed=nanstd(ExtraFPRTForPoorScene_Fixed)/sqrt(sum(~isnan(ExtraFPRTForPoorScene_Fixed)));
    SEM_ExtraFP_Normal_Fixed=nanstd(ExtraFPRTForNormalScene_Fixed)/sqrt(sum(~isnan(ExtraFPRTForNormalScene_Fixed)));
    
    %Statistical test
    [~,p_Rich_vs_Poor_Fixed]=ttest2(ExtraFPRTForRichScene_Fixed,ExtraFPRTForPoorScene_Fixed);
    [~,p_Rich_vs_Normal_Fixed]=ttest2(ExtraFPRTForRichScene_Fixed,ExtraFPRTForNormalScene_Fixed);
    [~,p_Normal_vs_Poor_Fixed]=ttest2(ExtraFPRTForNormalScene_Fixed,ExtraFPRTForPoorScene_Fixed);
    
    p_Fixed=[p_Rich_vs_Poor_Fixed,p_Rich_vs_Normal_Fixed,p_Normal_vs_Poor_Fixed]
    
   
    Mean_ExtraFP_Rich_Random=nanmean(ExtraFPRTForRichScene_Random);
    Mean_ExtraFP_Poor_Random=nanmean(ExtraFPRTForPoorScene_Random);
    Mean_ExtraFP_Normal_Random=nanmean(ExtraFPRTForNormalScene_Random);
    
    SEM_ExtraFP_Rich_Random=nanstd(ExtraFPRTForRichScene_Random)/sqrt(sum(~isnan(ExtraFPRTForRichScene_Random)));
    SEM_ExtraFP_Poor_Random=nanstd(ExtraFPRTForPoorScene_Random)/sqrt(sum(~isnan(ExtraFPRTForPoorScene_Random)));
    SEM_ExtraFP_Normal_Random=nanstd(ExtraFPRTForNormalScene_Random)/sqrt(sum(~isnan(ExtraFPRTForNormalScene_Random)));
    
    %Statistical test
    [~,p_Rich_vs_Poor_Random]=ttest2(ExtraFPRTForRichScene_Random,ExtraFPRTForPoorScene_Random);
    [~,p_Rich_vs_Normal_Random]=ttest2(ExtraFPRTForRichScene_Random,ExtraFPRTForNormalScene_Random);
    [~,p_Normal_vs_Poor_Random]=ttest2(ExtraFPRTForNormalScene_Random,ExtraFPRTForPoorScene_Random);
    
    p_Random=[p_Rich_vs_Poor_Random,p_Rich_vs_Normal_Random,p_Normal_vs_Poor_Random]
    
    
    
    
    if ishandle(125)
    close 125;end
    figure(125);set(125,'color','w');
    set(125,'Position', [100,100, 1000,500], 'Name', 'Before task motivation measurement');
    
   subplot(1,2,1); b=bar([Mean_ExtraFP_Poor_Fixed,Mean_ExtraFP_Normal_Fixed,Mean_ExtraFP_Rich_Fixed]);
   
   b.FaceColor = 'flat';
   b.CData = [1 0 0;0,1,0;0,0,1];
   xticklabels({'Poor','Normal','Rich'})
    title('Fixed Sequence Condition')
    hold on 
    errorbar([Mean_ExtraFP_Poor_Fixed,Mean_ExtraFP_Normal_Fixed,Mean_ExtraFP_Rich_Fixed],[SEM_ExtraFP_Poor_Fixed,SEM_ExtraFP_Normal_Fixed,SEM_ExtraFP_Rich_Fixed],'.k')
   ylabel('Reaction time to reach extra FP(ms)')
   box off
   

 subplot(1,2,2); b=bar([Mean_ExtraFP_Poor_Random,Mean_ExtraFP_Normal_Random,Mean_ExtraFP_Rich_Random]);
    b.FaceColor = 'flat';
   b.CData = [1 0 0;0,1,0;0,0,1];
   xticklabels({'Poor','Normal','Rich'})
    title('Random Sequence Condition')
    hold on 
    errorbar([Mean_ExtraFP_Poor_Random,Mean_ExtraFP_Normal_Random,Mean_ExtraFP_Rich_Random],[SEM_ExtraFP_Poor_Random,SEM_ExtraFP_Normal_Random,SEM_ExtraFP_Rich_Random],'.k')
   ylabel('Reaction time to reach extra FP(ms)')
   box off
        
   
    
end






FirstFPRTForRichScene_Fixed=RT_Fix_First((ismember(SceneIndex(:,2),RichSceneIndex))&FixedSeq);
    FirstFPRTForPoorScene_Fixed=RT_Fix_First((ismember(SceneIndex(:,2),PoorSceneIndex))&FixedSeq);
    FirstFPRTForNormalScene_Fixed=RT_Fix_First((ismember(SceneIndex(:,2),NormalSceneIndex))&FixedSeq);
    
    
    FirstFPRTForRichScene_Random=RT_Fix_First((ismember(SceneIndex(:,2),RichSceneIndex))&RandomSeq);
    FirstFPRTForPoorScene_Random=RT_Fix_First((ismember(SceneIndex(:,2),PoorSceneIndex))&RandomSeq);
    FirstFPRTForNormalScene_Random=RT_Fix_First((ismember(SceneIndex(:,2),NormalSceneIndex))&RandomSeq);
    
    
    
    Mean_FirstFP_Rich_Fixed=nanmean(FirstFPRTForRichScene_Fixed);
    Mean_FirstFP_Poor_Fixed=nanmean(FirstFPRTForPoorScene_Fixed);
    Mean_FirstFP_Normal_Fixed=nanmean(FirstFPRTForNormalScene_Fixed);
    
    SEM_FirstFP_Rich_Fixed=nanstd(FirstFPRTForRichScene_Fixed)/sqrt(sum(~isnan(FirstFPRTForRichScene_Fixed)));
    SEM_FirstFP_Poor_Fixed=nanstd(FirstFPRTForPoorScene_Fixed)/sqrt(sum(~isnan(FirstFPRTForPoorScene_Fixed)));
    SEM_FirstFP_Normal_Fixed=nanstd(FirstFPRTForNormalScene_Fixed)/sqrt(sum(~isnan(FirstFPRTForNormalScene_Fixed)));
    
    %Statistical test
    [~,p_Rich_vs_Poor_Fixed]=ttest2(FirstFPRTForRichScene_Fixed,FirstFPRTForPoorScene_Fixed);
    [~,p_Rich_vs_Normal_Fixed]=ttest2(FirstFPRTForRichScene_Fixed,FirstFPRTForNormalScene_Fixed);
    [~,p_Normal_vs_Poor_Fixed]=ttest2(FirstFPRTForNormalScene_Fixed,FirstFPRTForPoorScene_Fixed);
    
    p_Fixed_FP_First=[p_Rich_vs_Poor_Fixed,p_Rich_vs_Normal_Fixed,p_Normal_vs_Poor_Fixed]
    
   
    Mean_FirstFP_Rich_Random=nanmean(FirstFPRTForRichScene_Random);
    Mean_FirstFP_Poor_Random=nanmean(FirstFPRTForPoorScene_Random);
    Mean_FirstFP_Normal_Random=nanmean(FirstFPRTForNormalScene_Random);
    
    SEM_FirstFP_Rich_Random=nanstd(FirstFPRTForRichScene_Random)/sqrt(sum(~isnan(FirstFPRTForRichScene_Random)));
    SEM_FirstFP_Poor_Random=nanstd(FirstFPRTForPoorScene_Random)/sqrt(sum(~isnan(FirstFPRTForPoorScene_Random)));
    SEM_FirstFP_Normal_Random=nanstd(FirstFPRTForNormalScene_Random)/sqrt(sum(~isnan(FirstFPRTForNormalScene_Random)));
    
    %Statistical test
    [~,p_Rich_vs_Poor_Random]=ttest2(FirstFPRTForRichScene_Random,FirstFPRTForPoorScene_Random);
    [~,p_Rich_vs_Normal_Random]=ttest2(FirstFPRTForRichScene_Random,FirstFPRTForNormalScene_Random);
    [~,p_Normal_vs_Poor_Random]=ttest2(FirstFPRTForNormalScene_Random,FirstFPRTForPoorScene_Random);
    
    p_Random_FP_First=[p_Rich_vs_Poor_Random,p_Rich_vs_Normal_Random,p_Normal_vs_Poor_Random]


 
    if ishandle(126)
    close 126;end
    figure(126);set(126,'color','w');
    set(126,'Position', [100,100, 1000,500], 'Name', 'After Scene onset motivation measurement');
    
   subplot(1,2,1); b=bar([Mean_FirstFP_Poor_Fixed,Mean_FirstFP_Normal_Fixed,Mean_FirstFP_Rich_Fixed]);
   
   b.FaceColor = 'flat';
   b.CData = [1 0 0;0,1,0;0,0,1];
   xticklabels({'Poor','Normal','Rich'})
    title('Fixed Sequence Condition')
    hold on 
    errorbar([Mean_FirstFP_Poor_Fixed,Mean_FirstFP_Normal_Fixed,Mean_FirstFP_Rich_Fixed],[SEM_FirstFP_Poor_Fixed,SEM_FirstFP_Normal_Fixed,SEM_FirstFP_Rich_Fixed],'.k')
   ylabel('Reaction time to reach first FP(ms)')
   box off
   

 subplot(1,2,2); b=bar([Mean_FirstFP_Poor_Random,Mean_FirstFP_Normal_Random,Mean_FirstFP_Rich_Random]);
    b.FaceColor = 'flat';
   b.CData = [1 0 0;0,1,0;0,0,1];
   xticklabels({'Poor','Normal','Rich'})
    title('Random Sequence Condition')
    hold on 
    errorbar([Mean_FirstFP_Poor_Random,Mean_FirstFP_Normal_Random,Mean_FirstFP_Rich_Random],[SEM_FirstFP_Poor_Random,SEM_FirstFP_Normal_Random,SEM_FirstFP_Rich_Random],'.k')
   ylabel('Reaction time to reach first FP(ms)')
   box off



ChoiceGood=(ChoiceOutCome(:,2)==GOODCHOICECD);


    RTForRichScene_Fixed=RT_First_Target((ismember(SceneIndex(:,2),RichSceneIndex))&FixedSeq&ChoiceGood);
    RTForPoorScene_Fixed=RT_First_Target((ismember(SceneIndex(:,2),PoorSceneIndex))&FixedSeq&ChoiceGood);
    RTForNormalScene_Fixed=RT_First_Target((ismember(SceneIndex(:,2),NormalSceneIndex))&FixedSeq&ChoiceGood);
    
    
    RTForRichScene_Random=RT_First_Target((ismember(SceneIndex(:,2),RichSceneIndex))&RandomSeq&ChoiceGood);
    RTForPoorScene_Random=RT_First_Target((ismember(SceneIndex(:,2),PoorSceneIndex))&RandomSeq&ChoiceGood);
    RTForNormalScene_Random=RT_First_Target((ismember(SceneIndex(:,2),NormalSceneIndex))&RandomSeq&ChoiceGood);
    
    
    
    Mean_RT_Rich_Fixed=nanmean(RTForRichScene_Fixed);
    Mean_RT_Poor_Fixed=nanmean(RTForPoorScene_Fixed);
    Mean_RT_Normal_Fixed=nanmean(RTForNormalScene_Fixed);
    
    SEM_RT_Rich_Fixed=nanstd(RTForRichScene_Fixed)/sqrt(sum(~isnan(RTForRichScene_Fixed)));
    SEM_RT_Poor_Fixed=nanstd(RTForPoorScene_Fixed)/sqrt(sum(~isnan(RTForPoorScene_Fixed)));
    SEM_RT_Normal_Fixed=nanstd(RTForNormalScene_Fixed)/sqrt(sum(~isnan(RTForNormalScene_Fixed)));
    
    %Statistical test
    [~,p_Rich_vs_Poor_Fixed]=ttest2(RTForRichScene_Fixed,RTForPoorScene_Fixed);
    [~,p_Rich_vs_Normal_Fixed]=ttest2(RTForRichScene_Fixed,RTForNormalScene_Fixed);
    [~,p_Normal_vs_Poor_Fixed]=ttest2(RTForNormalScene_Fixed,RTForPoorScene_Fixed);
    
    p_Fixed_FP_First=[p_Rich_vs_Poor_Fixed,p_Rich_vs_Normal_Fixed,p_Normal_vs_Poor_Fixed]
    
   
    Mean_RT_Rich_Random=nanmean(RTForRichScene_Random);
    Mean_RT_Poor_Random=nanmean(RTForPoorScene_Random);
    Mean_RT_Normal_Random=nanmean(RTForNormalScene_Random);
    
    SEM_RT_Rich_Random=nanstd(RTForRichScene_Random)/sqrt(sum(~isnan(RTForRichScene_Random)));
    SEM_RT_Poor_Random=nanstd(RTForPoorScene_Random)/sqrt(sum(~isnan(RTForPoorScene_Random)));
    SEM_RT_Normal_Random=nanstd(RTForNormalScene_Random)/sqrt(sum(~isnan(RTForNormalScene_Random)));
    
    %Statistical test
    [~,p_Rich_vs_Poor_Random]=ttest2(RTForRichScene_Random,RTForPoorScene_Random);
    [~,p_Rich_vs_Normal_Random]=ttest2(RTForRichScene_Random,RTForNormalScene_Random);
    [~,p_Normal_vs_Poor_Random]=ttest2(RTForNormalScene_Random,RTForPoorScene_Random);
    
    p_Random_FP_First=[p_Rich_vs_Poor_Random,p_Rich_vs_Normal_Random,p_Normal_vs_Poor_Random]


 
    if ishandle(127)
    close 127;end
    figure(127);set(127,'color','w');
    set(127,'Position', [100,100, 1000,500], 'Name', 'RT to the chosen target');
    
   subplot(1,2,1); b=bar([Mean_RT_Poor_Fixed,Mean_RT_Normal_Fixed,Mean_RT_Rich_Fixed]);
   
   b.FaceColor = 'flat';
   b.CData = [1 0 0;0,1,0;0,0,1];
   xticklabels({'Poor','Normal','Rich'})
    title('Fixed Sequence Condition')
    hold on 
    errorbar([Mean_RT_Poor_Fixed,Mean_RT_Normal_Fixed,Mean_RT_Rich_Fixed],[SEM_RT_Poor_Fixed, SEM_RT_Normal_Fixed,SEM_RT_Rich_Fixed],'.k')
   ylabel('Reaction time to reach Chosen Target(ms)')
   box off
   

 subplot(1,2,2); b=bar([Mean_RT_Poor_Random,Mean_RT_Normal_Random,Mean_RT_Rich_Random]);
    b.FaceColor = 'flat';
   b.CData = [1 0 0;0,1,0;0,0,1];
   xticklabels({'Poor','Normal','Rich'})
    title('Random Sequence Condition')
    hold on 
    errorbar([Mean_RT_Poor_Random,Mean_RT_Normal_Random,Mean_RT_Rich_Random],[SEM_RT_Poor_Random,SEM_RT_Normal_Random,SEM_RT_Rich_Random],'.k')
   ylabel('Reaction time to reach Chosen Target(ms)')
   box off

end


end