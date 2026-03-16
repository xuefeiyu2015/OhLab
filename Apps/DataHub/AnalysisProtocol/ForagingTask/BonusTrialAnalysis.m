function BonusTrialAnalysis(Data, StartTrial, EndTrial,ShowFigureFlag,OutputFlag);
%This is the function to analysis the bonus trials
%{
//1.Trained Scene Within---Coupled Objects v.s. UnCoupled(Within)
	//2.Trained Scene Within---Coupled Objects v.s. UnCoupled(Between)
	//3.Trained Scene Within--- Coupled Objects v.s. Novel Objects
	//4.Trained Scene Within---UnCoupled Objects(Within) v.s. UnCoupled(Between)
	//5.Trained Scene Within---Novel v.s. Novel

	//6.Trained Scene Between---Trained Objects(Within) v.s. Novel Objects
	//7.Trained Scene Between---Trained Objects(Between) v.s. Novel Objects
	//8.Trained Scene Between---Trained Objects(Within) v.s. Trained(Between)
	//9.Trained Scene Between---Novel v.s. Novel

	//10.Novel Scene--Trained Objects(Within) v.s. Novel Objects
	//11.Novel Scene---Trained Objects(Between) v.s. Novel Objects
	//12.Novel Scene---Trained Objects(Between) v.s. Trained Objects(Within)
	//13.Novel Scene---Novel Objects v.s. Novel Objects

%}

%Xuefei Yu; 02182020
ECode;


BonusTrialDataRaw=Data.BonusData;
FixPointIndexBonus=BonusTrialDataRaw.FixPointIndexBonus;
FixPointLocBonus=BonusTrialDataRaw.FixPointLocBonus;
SceneIndexBonus=BonusTrialDataRaw.SceneIndexBonus;
BonusTrialCondition=BonusTrialDataRaw.BonusTrialCondition;
TargetLocBonus=BonusTrialDataRaw.TargetLocBonus;
TargetIndexBonus=BonusTrialDataRaw.TargetIndexBonus;
TargetPolarLocBonus=BonusTrialDataRaw.TargetPolarLocBonus;

EyeChannel_X_L=BonusTrialDataRaw.EyeChannel_X_L;
EyeChannel_Y_L=BonusTrialDataRaw.EyeChannel_Y_L;

EyeBinWidth=BonusTrialDataRaw.EyeBinWidth;
EyeTime=BonusTrialDataRaw.EyeTime;

EventTimeChannel=BonusTrialDataRaw.EventTimeChannel;
EventChannel=BonusTrialDataRaw.EventChannel;
EventBin=BonusTrialDataRaw.EventBin;

SpikeTimeData=BonusTrialDataRaw.SpikeTimeData;
SpikeBinData=BonusTrialDataRaw.SpikeBinData;

%{
%Select the interval between the start trial and the end trial
PreSelect=StartTrial:EndTrial;
FixPointIndexBonus=FixPointIndexBonus(PreSelect,:);
FixPointLocBonus=FixPointLocBonus(PreSelect,:);
SceneIndexBonus=SceneIndexBonus(PreSelect,:);
BonusTrialCondition=BonusTrialCondition(PreSelect,:);
TargetLocBonus=TargetLocBonus(PreSelect,:);
TargetIndexBonus=TargetIndexBonus(PreSelect,:);
TargetPolarLocBonus=TargetPolarLocBonus(PreSelect,:);

EyeChannel_X_L=EyeChannel_X_L(PreSelect,:);
EyeChannel_Y_L=EyeChannel_Y_L(PreSelect,:);

EyeTime=EyeTime(PreSelect,:);
EventTimeChannel=EventTimeChannel(PreSelect,:);
EventChannel=EventChannel(PreSelect,:);
EventBin=EventBin(PreSelect,:);
if sum(isnan(SpikeTimeData),'all')~=numel(SpikeTimeData);
SpikeTimeData=SpikeTimeData(PreSelect,:);
SpikeBinData=SpikeBinData(PreSelect,:);
end
%}
%Select out the trials which at least one pair of fractals was presented;
BonusTrialIndicator=BonusTrialDataRaw.BonusTrialNumberIndex;    

TrialIndexSelection=unique(TargetIndexBonus(:,1));
Selection=ismember(BonusTrialIndicator,TrialIndexSelection);


SceneIndexBonusUse=SceneIndexBonus(ismember(SceneIndexBonus(:,1),TrialIndexSelection),:);
BonusTrialConditionUse=BonusTrialCondition(ismember(BonusTrialCondition(:,1),TrialIndexSelection),:);
TargetPolarLocBonusUse=TargetPolarLocBonus(ismember(TargetPolarLocBonus(:,1),TrialIndexSelection),:);
TargetLocBonus=TargetLocBonus(ismember(TargetLocBonus(:,1),TrialIndexSelection),:);
TargetIndexBonusUse=TargetIndexBonus(ismember(TargetIndexBonus(:,1),TrialIndexSelection),:);




%Reconstruct the monkey's "instantaneous choice" according to eye trace
EyeChannel_X_LUse=EyeChannel_X_L(Selection,:);
EyeChannel_Y_LUse=EyeChannel_Y_L(Selection,:);

EyeTimeUse=EyeTime(Selection,:);

EventTimeChannelUse=EventTimeChannel(Selection,:);
EventChannelUse=EventChannel(Selection,:);
EventBinUse=EventBin(Selection,:);

%Chunk the eye data
FPOnsetBin = FindOutTime(EventChannelUse,EventBinUse,FIXSTARTCD);
FractalOnsetBin = FindOutTime(EventChannelUse,EventBinUse,TGTCD );
FractalOffsetBin = FindOutTime(EventChannelUse,EventBinUse,TGTOFF);

FPOnsetBin(isnan(FractalOnsetBin))=NaN;


EyeX_in_bin = ChunkData(BonusTrialIndicator(Selection),EyeChannel_X_LUse,FractalOnsetBin ,FractalOffsetBin);
EyeY_in_bin = ChunkData(BonusTrialIndicator(Selection),EyeChannel_Y_LUse,FractalOnsetBin ,FractalOffsetBin);
EyeTime_in_Bin=ChunkData(BonusTrialIndicator(Selection),EyeTimeUse,FractalOnsetBin,FractalOffsetBin);

AlignTime= FindOutTime(EventChannelUse,EventTimeChannelUse,TGTCD );
AlignTime=reshape(AlignTime',1,numel(AlignTime))';
AlignTime=AlignTime(~isnan(AlignTime));

EyeTime_in_Bin_Align=EyeTime_in_Bin-AlignTime;

EyeTime_in_Bin_Align=[EyeTime_in_Bin(:,1),EyeTime_in_Bin_Align(:,2:end)];


%Screen out the eye out of the screen
temp=EyeX_in_bin(:,2:end);
temp(temp>40)=NaN;
EyeX_in_bin(:,2:end)=temp;

EyeX_in_bin_V=diff(EyeX_in_bin(:,2:end),1,2)/EyeBinWidth*1000;
EyeX_in_bin_V=[EyeX_in_bin_V(:,1) ,EyeX_in_bin_V];



temp=EyeY_in_bin(:,2:end);
temp(temp>40)=NaN;
EyeY_in_bin(:,2:end)=temp;

EyeY_in_bin_V=diff(EyeY_in_bin(:,2:end),1,2)/EyeBinWidth*1000;
EyeY_in_bin_V=[EyeY_in_bin_V(:,1) ,EyeY_in_bin_V];

%Criterias for selecting saccades
V_Threshold=40;
ContinueBin=5;
ISI_Threshold=20;

Saccades=SelectSaccade(EyeX_in_bin(:,1),EyeX_in_bin(:,2:end),EyeY_in_bin(:,2:end),EyeTime_in_Bin_Align(:,2:end),EyeBinWidth,V_Threshold,ContinueBin,ISI_Threshold);


SaccadeNumber=arrayfun(@(x)  x.NumOfSaccade,Saccades)';
SaccadeAngle=arrayfun(@(x)  x.SaccadeAngle,Saccades,'UniformOutput' ,0)';
SaccadeAmplitude=arrayfun(@(x)  x.SaccadeAmplitude,Saccades,'UniformOutput' ,0)';

SaccadeStartTime=arrayfun(@(x)  x.SaccadeStartTime,Saccades,'UniformOutput' ,0)';
SaccadeEndTime=arrayfun(@(x)  x.SaccadeEndTime,Saccades,'UniformOutput' ,0)';

SaccadeStartPoint=arrayfun(@(x)  x.SaccadeStartPoint,Saccades,'UniformOutput' ,0)';
SaccadeEndPoint=arrayfun(@(x)  x.SaccadeEndPoint,Saccades,'UniformOutput' ,0)';

%Screen out the large saccades and saccade start time above 50;

AmplitudeThreshold=3;
%AmplitudeThreshold=0;

TimeThreshold=50;

SaccadeAngle_Selection=cellfun(@(x,y,z) z((x>AmplitudeThreshold)&(y>TimeThreshold)),SaccadeAmplitude,SaccadeStartTime,SaccadeAngle,'UniformOutput' ,0);
SaccadeEndPoint_Selection=cellfun(@(x,y,z) z(((x>AmplitudeThreshold)&(y>TimeThreshold)),:),SaccadeAmplitude,SaccadeStartTime,SaccadeEndPoint,'UniformOutput' ,0);
SaccadeStartPoint_Selection=cellfun(@(x,y,z) z(((x>AmplitudeThreshold)&(y>TimeThreshold)),:),SaccadeAmplitude,SaccadeStartTime,SaccadeStartPoint,'UniformOutput' ,0);

SaccadeStartTime_Selection=cellfun(@(x,y,z) z(((x>AmplitudeThreshold)&(y>TimeThreshold)),:),SaccadeAmplitude,SaccadeStartTime,SaccadeStartTime,'UniformOutput' ,0);
%Only look at the first saccade

SaccadeNumInEachTrial=cellfun(@numel,SaccadeAngle_Selection);

SaccadeAngle_First=cell2mat(cellfun(@(x) x(1,:),SaccadeAngle_Selection(SaccadeNumInEachTrial>0),'UniformOutput' ,0));
SaccadeEndPoint_First=cell2mat(cellfun(@(x) x(1,:),SaccadeEndPoint_Selection(SaccadeNumInEachTrial>0,:),'UniformOutput' ,0));
RT_First=cell2mat(cellfun(@(x) x(1,:),SaccadeStartTime_Selection(SaccadeNumInEachTrial>0,:),'UniformOutput' ,0));
SaccadeStartPoint_First=cell2mat(cellfun(@(x) x(1,:),SaccadeStartPoint_Selection(SaccadeNumInEachTrial>0,:),'UniformOutput' ,0));



%CheckChoice
TargetLoc1=TargetLocBonus(SaccadeNumInEachTrial>0,2:3);
TargetLoc2=TargetLocBonus(SaccadeNumInEachTrial>0,4:5);

TargetAngle1=LocToAngle(TargetLoc1(:,1),TargetLoc1(:,2));
TargetAngle2=LocToAngle(TargetLoc2(:,1),TargetLoc2(:,2));

distanceToTarget=[sqrt(sum((SaccadeEndPoint_First-TargetLoc1).^2,2)),sqrt(sum((SaccadeEndPoint_First-TargetLoc2).^2,2))];


distanceToTarget_Angle=abs([TargetAngle1,TargetAngle2]-SaccadeAngle_First);
distanceToTarget_Angle(distanceToTarget_Angle>180)=360-distanceToTarget_Angle(distanceToTarget_Angle>180);

Choice_Angle=distanceToTarget_Angle<90;

Choice_Distance=distanceToTarget==min(distanceToTarget,[],2);

Choice=Choice_Distance;


%Check the bonus trial condition


%Screen trial index


if size(BonusTrialConditionUse,1)~=size(TargetIndexBonusUse,1)
    %Restore currCon from raw scene and fractal information
    BonusTrialConditionUse=BonusTrialConditionRescue(SceneIndexBonusUse,TargetIndexBonusUse);
  
end

ContainCoupleCondition=[1,2];

ContainCoupleIndex=[BonusTrialConditionUse(:,1),ismember(BonusTrialConditionUse(:,2),ContainCoupleCondition)];
MightCoupleCondition=[7,8];
%In trial condition 7 & 8; there might be ocassionaly coupled scene
%condition


CorrectConditionSelection=ismember(BonusTrialConditionUse(:,2),MightCoupleCondition);
SceneIndexCheck=SceneIndexBonusUse(ismember(SceneIndexBonusUse(:,1),BonusTrialConditionUse(CorrectConditionSelection,1)),:);
FractalIndexCheck=TargetIndexBonusUse(ismember(TargetIndexBonusUse(:,1),BonusTrialConditionUse(CorrectConditionSelection,1)),:);

ContainCoupleIndex=BonusTrialConditionCheck(ContainCoupleIndex,SceneIndexCheck,FractalIndexCheck,BonusTrialConditionUse(CorrectConditionSelection,:));


ContainCoupleIndex_ForAnalysis=ContainCoupleIndex(SaccadeNumInEachTrial>0,:);

BonusTrialConditionSelection=BonusTrialConditionUse(SaccadeNumInEachTrial>0,:);

%Separate Good and Bad Object
TargetIndexForChoice=TargetIndexBonusUse(SaccadeNumInEachTrial>0,[2,4]);
TargetTypeForChoice=FractalTypeIdentifyBonus(TargetIndexForChoice);%1:Good;2:Bad:3:Novel


UniqueBonusCondition=unique(BonusTrialConditionSelection(:,2));

for j=1:size(UniqueBonusCondition,1)
    
     Select=BonusTrialConditionSelection(:,2)==UniqueBonusCondition(j);
    TrialNumOfCondition(j)=sum(Select);
    NumberOfChoice1(j)=sum(Choice(Select,1)==1);
    NumberOfChoice2(j)=TrialNumOfCondition(j)-NumberOfChoice1(j );
    
    
    FractalTypeCondtion=TargetTypeForChoice(Select,:);
    
 
    
    %NumberOfChoice1Good=Choice(Select&(TargetTypeForChoice==1),1);
    
    
    Choice1{j}=Choice(Select,1)==1;
    Choice2{j}=Choice(Select,2)==1;
    
    NumberOfChoice1Good(j)=sum((Choice(Select,1)==1)&sum(FractalTypeCondtion==1,2)>0);
    NumberOfTrialsGood(j)=sum(sum(FractalTypeCondtion==1,2)>0);
    
   
    
    
    NumberOfChoice1Bad(j)=sum((Choice(Select,1)==1)&sum(FractalTypeCondtion==2,2)>0);
    NumberOfTrialsBad(j)=sum(sum(FractalTypeCondtion==2,2)>0);
    
   
    RT_Selection=RT_First(Select,:);
    
    RTOfChoice1_mean(j)=mean(RT_Selection(Choice(Select,1)==1,:));
    RTOfChoice2_mean(j)=mean(RT_Selection(Choice(Select,2)==1,:));
    
    
    RTOfChoice1_Goodmean(j)=mean(RT_Selection((Choice(Select,1)==1)&( sum(FractalTypeCondtion==1,2)>0),:));
    RTOfChoice2_Goodmean(j)=mean(RT_Selection((Choice(Select,2)==1)&( sum(FractalTypeCondtion==1,2)>0),:));
    
    RTOfChoice1_Badmean(j)=mean(RT_Selection((Choice(Select,1)==1)&( sum(FractalTypeCondtion==2,2)>0),:));
    RTOfChoice2_Badmean(j)=mean(RT_Selection((Choice(Select,2)==1)&( sum(FractalTypeCondtion==2,2)>0),:));
    
    
    RTOfChoice1_Good{j}=RT_Selection((Choice(Select,1)==1)&( sum(FractalTypeCondtion==1,2)>0),:);
    RTOfChoice2_Good{j}=RT_Selection((Choice(Select,2)==1)&( sum(FractalTypeCondtion==1,2)>0),:);
    
    RTOfChoice1_Bad{j}=RT_Selection((Choice(Select,1)==1)&( sum(FractalTypeCondtion==2,2)>0),:);
    RTOfChoice2_Bad{j}=RT_Selection((Choice(Select,2)==1)&( sum(FractalTypeCondtion==2,2)>0),:);
    
   
    
    RTOfChoice1{j}=RT_Selection(Choice(Select,1)==1,:);
    RTOfChoice2{j}=RT_Selection(Choice(Select,2)==1,:);
    
    
    
    
end

ProportionOfChoices1=NumberOfChoice1./ TrialNumOfCondition;
ProportionOfChoices2=1-ProportionOfChoices1;

ProportionOfChoices1_Good=NumberOfChoice1Good./ NumberOfTrialsGood;
ProportionOfChoices2_Good=1-ProportionOfChoices1_Good;


ProportionOfChoices1_Bad=NumberOfChoice1Bad./ NumberOfTrialsBad;
ProportionOfChoices2_Bad=1-ProportionOfChoices1_Bad;


%Export file
ExportResult=1;
%OutputPath
if ExportResult
CurrPath=pwd;

OutputFileName='Results\BonusTrial';

if contains(CurrPath,OutputFileName)
    search=1;
else
   
    
  Computer= computer;
   if sum(Computer(1:2)=='PC')>0
    
   splitIndex='\';
  else
    
    splitIndex='/';
   end
Separator=find(CurrPath==splitIndex,1,'last');

CurrPath=CurrPath(1:Separator);
OutputPath=strcat(CurrPath,OutputFileName);
  if ~exist(OutputPath)
    mkdir(OutputPath);
    
  end
    
    
end




cd(OutputPath);
FileName=Data.FileName;
FileName=FileName(1:end-4);
FileName=strcat(FileName,'_Bonus.mat');


OutputData.NumberOfChoice1Good=NumberOfChoice1Good;
OutputData.NumberOfTrialsGood=NumberOfTrialsGood;

OutputData.NumberOfChoice1Bad=NumberOfChoice1Bad;
OutputData.NumberOfTrialBad=NumberOfTrialsBad;


OutputData.NumberOfChoice1=NumberOfChoice1;
OutputData.TrialNumOfCondition=TrialNumOfCondition;

OutputData.RTOfChoice1=RTOfChoice1;
OutputData.RTOfChoice2=RTOfChoice2;

OutputData.RTOfChoice1_Good=RTOfChoice1_Good;
OutputData.RTOfChoice1_Bad=RTOfChoice1_Bad;

OutputData.RTOfChoice2_Good=RTOfChoice2_Good;
OutputData.RTOfChoice2_Bad=RTOfChoice2_Bad;



save(FileName,'OutputData');
disp('Data Exported');
end


figure
set(gcf,'color','w')
subplot(3,1,1);bar([ProportionOfChoices1',ProportionOfChoices2'])
ylabel('Proportion of choices')

legend('Choice1','Choice2')
title('Choice proportion ')
subplot(3,1,2);bar([ProportionOfChoices1_Good',ProportionOfChoices2_Good'])
ylabel('Proportion of choices')

subplot(3,1,3);bar([ProportionOfChoices1_Bad',ProportionOfChoices2_Bad'])
xlabel('Trial Condition')
ylabel('Proportion of choices')


figure

set(gcf,'color','w')

subplot(3,1,1);bar([RTOfChoice1_mean',RTOfChoice2_mean'])
ylabel('RT')
title('Reaction Time')
legend('Choice1','Choice2')

subplot(3,1,2);bar([RTOfChoice1_Goodmean',RTOfChoice2_Goodmean'])
ylabel('RT')

subplot(3,1,3);bar([RTOfChoice1_Badmean',RTOfChoice2_Badmean'])
xlabel('Trial Condition')
ylabel('RT')



end