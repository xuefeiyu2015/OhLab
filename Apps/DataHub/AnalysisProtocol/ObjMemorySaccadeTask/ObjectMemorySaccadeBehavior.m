function  ObjectMemorySaccadeBehavior(Data, StartTrial, EndTrial,ShowFigureFlag,OutputFlag);
%Memory saccade behavior analysis: RT--Good Bad difference
%added by YXF 06102021
ECode;


%PlotStartOffset=[2000];%in ms
%PlotStopOffset=[5000];%in ms
%For plots
PlotRangeSample=[-500,1000];
PlotRangeTest=[-1000,3000];

%Time bin pre saccade and time bin offset after sample onset
SampleOnPreBuffer=500;%in ms
SampleOnPostBuffer=5000;%in ms



ResponseLatencyThreshold=3;


%Criterias for sexwlecting saccades
V_Threshold=40;
ContinueBin=5;
ISI_Threshold=20;

%Color Definition for each fractals
%{
Color=[lines(6);prism(6) ];
Color(2,:)=copper(1)+0.5;
Color=[Color;Color*0.5];
%}
Color={'#648FFF','#FE6100','#009E73','#88CCEE','#CC6677','#785EF0'};
Color_RGB=[
           254,97,0;
           100,143,255;
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

%SpikeChannel=Data.SpikeChannel;




%Trial outcome;
EventChannel=Data.EventChannel(Selection,:);
EventTimeChannel=Data.EventTimeChannel(Selection,:);
EventBin=Data.EventBin(Selection,:);

%GetEyeData
EyeData=Data.EyeDataRaw;
EyeChannel_X_L=EyeData(Selection,1);
EyeChannel_Y_L=EyeData(Selection,2);
%EyeChannel_X_R=EyeData(Selection,3);
%EyeChannel_Y_R=EyeData(Selection,4);

PuipleThreshold=-50;
PuipleSize=ValidPuiple(EyeData(Selection,3),PuipleThreshold);%Only valid after some date....


EyeBinWidth=Data.EyeBinWidth;



%CorrectTrials=~(sum(ismember(EventChannel,UNMATCHOBJ),2)==1);

FullySuccessfulTrials=sum(ismember(EventChannel,ADVEYE_IN),2)==1 & sum(ismember(EventChannel,ADVEYE_STAY),2)==1;

FailedGazeTrials=sum(ismember(EventChannel,ADVEYE_IN),2)==1 & sum(ismember(EventChannel,ADVEYE_STAY),2)==0;

FailedMemoryTrials=sum(ismember(EventChannel,ADVEYE_IN),2)==0 & sum(ismember(EventChannel,ADVEYE_STAY),2)==0;



disp(['Number of Fully Sucessful Memory trials=',num2str(sum(FullySuccessfulTrials))]);
disp(['Number of Failed Gaze Trials=',num2str(sum(FailedGazeTrials))]);
disp(['Number of Failed Memory Trials=',num2str(sum(FailedMemoryTrials))]);




%Target Information
TargetIndex=Data.TargetIndex(Selection,2);%For this task, there is only target name, but no target index

TargetName=Data.TargetIndex(Selection,3);%For this task, there is only target name, but no target index

%TargetLoc=Data.TargetLoc(Selection,:,:);

TargetNum=Data.TargetNum(Selection,2:end);
%{
%Separate single target trials with two target trials

TaskType_ori=Data.ParaLib('TaskType');

TaskType=TaskType_ori(Selection,2:end);

SingleTgtTrial=TaskType==1;
ChoiceTrial=TaskType>1;
%}
%SingleTgtTrial=TargetNum==1;
%ChoiceTrial=TargetNum>1;


ObjectIndex=ReproduceFromEvent(EventChannel,[GOODOBJ  BADOBJ])';
GoodObjectIndex=(ObjectIndex==GOODOBJ)';

%Reproduce the fractal presentation time

FractalOnsetTime = FindOutTime(EventChannel,EventTimeChannel,TGTCD );
FractalOffsetTime_all = FindOutTime(EventChannel,EventTimeChannel,TGTOFF );

FractalOffsetTime=FractalOffsetTime_all(:,1);

FractalPresentationDur=FractalOffsetTime-FractalOnsetTime;

%Separate into 2 interval
Short_The=100;
Long_The=300;
ShortPresent=FractalPresentationDur<=Short_The;
LongPresent=FractalPresentationDur>=Long_The;

disp(['Short trials=',num2str(sum(ShortPresent))]);
disp(['Long trials=',num2str(sum(LongPresent))]);

ShortDur=FractalPresentationDur(ShortPresent);
LongDur=FractalPresentationDur(LongPresent);


x_short=(rand(length(ShortDur),1)-0.5)*0.1+1;
x_long=(rand(length(LongDur),1)-0.5)*0.1+2;




%PuipleSize_SampleSelection_tmp= SelectEyeInterval(PuipleSize,SampleOnsetBinRange ,SampleOffsetBinRange);
%PuipleSize_SampleSelection=ReorganizeEye(PuipleSize_SampleSelection_tmp);
PuipleSize_SampleSelection=ReorganizeEye({PuipleSize});

%Align puiple data
EyeTime=0:EyeBinWidth:size(PuipleSize_SampleSelection,2)-1;

PuipleSize_Aligned=AlignEyeData(PuipleSize_SampleSelection,EyeTime,FractalOnsetTime,[-SampleOnPreBuffer,SampleOnPostBuffer]);

EyeTime_Sample=[-SampleOnPreBuffer:-SampleOnPreBuffer+size(PuipleSize_Aligned,2)-1];

SponIntervalSample=[-400,-100];
SponPreSample=EyeTime_Sample<SponIntervalSample(2) & EyeTime_Sample>SponIntervalSample(1);

PuipleSize_Aligned_SponEachTrial=nanmean(PuipleSize_Aligned(:,SponPreSample),2);
PuipleSize_AlignedOverspon=PuipleSize_Aligned-PuipleSize_Aligned_SponEachTrial;

PuipleSize_SampleAlign_Good=PuipleSize_AlignedOverspon(GoodObjectIndex,:);
PuipleSize_SampleAlign_Bad=PuipleSize_AlignedOverspon(~GoodObjectIndex,:);

[PuipleSize_SampleAlign_Good_Mean,PuipleSize_SampleAlign_Good_Sem]=StatisticsPSTH(PuipleSize_SampleAlign_Good,3);
[PuipleSize_SampleAlign_Bad_Mean,PuipleSize_SampleAlign_Bad_Sem]=StatisticsPSTH(PuipleSize_SampleAlign_Bad,3);

%Successful rate
FullySuccessfulRate_Short=sum(FullySuccessfulTrials(ShortPresent))/sum(ShortPresent);
FailedGazeRate_Short=sum(FailedGazeTrials(ShortPresent))/sum(ShortPresent);
FailedMemoryRate_Short=sum(FailedMemoryTrials(ShortPresent))/sum(ShortPresent);

FullySuccessfulRate_Long=sum(FullySuccessfulTrials(LongPresent))/sum(LongPresent);
FailedGazeRate_Long=sum(FailedGazeTrials(LongPresent))/sum(LongPresent);
FailedMemoryRate_Long=sum(FailedMemoryTrials(LongPresent))/sum(LongPresent);


SuccessfulRate_Array=[FullySuccessfulRate_Short,FailedGazeRate_Short,FailedMemoryRate_Short,FullySuccessfulRate_Long,FailedGazeRate_Long,FailedMemoryRate_Long];
SuccessfulRate_Str={'Suc_S','FailGaze_S','Failure_S','Suc_L','FailGaze_L','Failure_L'};

%Finally saccade latency

%Select the period between fixation offset target reshow
FixOffsetBin = FindOutTime(EventChannel,EventBin,FIXOFF	)-100;%-100: offset
FractalReshowBin= FindOutTime(EventChannel,EventBin,RETGT);

%FractalOnsetTime = FindOutTime(EventChannel,EventTimeChannel,TGTCD );
%FractalOffsetTime = FindOutTime(EventChannel,EventTimeChannel,TGTOFF );

EyeChannel_X_Selection_All= SelectEyeInterval(EyeChannel_X_L,FixOffsetBin,FractalReshowBin )';
EyeChannel_Y_Selection_All= SelectEyeInterval(EyeChannel_Y_L,FixOffsetBin ,FractalReshowBin)';

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

AmplitudeThreshold=3;

TimeThreshold=-100;

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


RT_First_Target=RT_First_Target-100;

RT_First_Target=RT_First_Target';

RT_First_Target_Good=RT_First_Target(GoodObjectIndex);
RT_First_Target_Bad=RT_First_Target(~GoodObjectIndex);



RT_First_Target_Good_SucShort=RT_First_Target(ShortPresent&FullySuccessfulTrials&GoodObjectIndex');
RT_First_Target_Bad_SucShort=RT_First_Target(ShortPresent&FullySuccessfulTrials&(~GoodObjectIndex'));


RT_First_Target_Good_SucLong=RT_First_Target(LongPresent&FullySuccessfulTrials&GoodObjectIndex');
RT_First_Target_Bad_SucLong=RT_First_Target(LongPresent&FullySuccessfulTrials&(~GoodObjectIndex'));


RT_First_Target_Good_SucShort_Mean=nanmean(RT_First_Target_Good_SucShort);
RT_First_Target_Good_SucShort_Sem=nanstd(RT_First_Target_Good_SucShort)/sqrt(sum(~isnan(RT_First_Target_Good_SucShort)));

RT_First_Target_Bad_SucShort_Mean=nanmean(RT_First_Target_Bad_SucShort);
RT_First_Target_Bad_SucShort_Sem=nanstd(RT_First_Target_Bad_SucShort)/sqrt(sum(~isnan(RT_First_Target_Bad_SucShort)));

RT_First_Target_Good_SucLong_Mean=nanmean(RT_First_Target_Good_SucLong);
RT_First_Target_Good_SucLong_Sem=nanstd(RT_First_Target_Good_SucLong)/sqrt(sum(~isnan(RT_First_Target_Good_SucLong)));

RT_First_Target_Bad_SucLong_Mean=nanmean(RT_First_Target_Bad_SucLong);
RT_First_Target_Bad_SucLong_Sem=nanstd(RT_First_Target_Bad_SucLong)/sqrt(sum(~isnan(RT_First_Target_Bad_SucLong)));



RT_Array_Mean=[RT_First_Target_Good_SucShort_Mean,RT_First_Target_Bad_SucShort_Mean,RT_First_Target_Good_SucLong_Mean,RT_First_Target_Bad_SucLong_Mean];
RT_Array_Sem=[RT_First_Target_Good_SucShort_Sem,RT_First_Target_Bad_SucShort_Sem,RT_First_Target_Good_SucLong_Sem,RT_First_Target_Bad_SucLong_Sem];

RT_Str={'GS','BS','GL','BL'};

%Statistics
[h,p_Short_Value]=ttest2(RT_First_Target_Good_SucShort,RT_First_Target_Bad_SucShort);
[h,p_Long_Value]=ttest2(RT_First_Target_Good_SucLong,RT_First_Target_Bad_SucLong);



if ShowFigureFlag
     figtitlestr{1}='PSTHPuipleSizeChange';

     fig{1}=PrepareFigure(1,'w',[50,100, 1200,800],'Name',figtitlestr{1});
    
     sel=EyeTime_Sample>PlotRangeSample(1) & EyeTime_Sample<PlotRangeSample(2);
    %subplot(1,2,1);%Fractal On
    axes('Position',[0.12,0.7,0.5,0.2]);
     shadedErrorBar(EyeTime_Sample( sel),PuipleSize_SampleAlign_Good_Mean( sel),PuipleSize_SampleAlign_Good_Sem( sel),'lineprops',{[1    0    0]},'transparent',1,'patchSaturation',0.3);

    hold on
     shadedErrorBar(EyeTime_Sample( sel),PuipleSize_SampleAlign_Bad_Mean( sel),PuipleSize_SampleAlign_Bad_Sem( sel),'lineprops',{[0    0    1]},'transparent',1,'patchSaturation',0.3);
     %{
     %Add marker line
     yrange=ylim;
     xlim(PlotRangeSample);
     plot([0,0],yrange,'--k');  
     plot([SampleOffMarker_Mean,SampleOffMarker_Mean],yrange,'--k');
     plot([ChoiceTgtOnMarker_Mean,  ChoiceTgtOnMarker_Mean],yrange,'--r');
     plot([ChoiceTgtOffMarker_Mean,  ChoiceTgtOffMarker_Mean],yrange,'--r');
     plot([ RwdMarker_Mean,   RwdMarker_Mean],yrange,'--y');
     %}
    
     set(findall(gca, 'Type', 'Line'),'LineWidth',3);
     set(gca,'FontSize',15,'FontWeight','Bold');
set(gca,'LineWidth',3);

xlabel('Time from Fractal Onset(ms)','FontSize',25,'FontWeight','Bold');
ylabel('PuipleSize(a.u.)','FontSize',25,'FontWeight','Bold');
title('Fractal Period');

%Two interval time distribution: 
axes('Position',[0.7,0.7,0.2,0.2]);
plot(x_short,ShortDur,'.r');
hold on
plot(x_long,LongDur,'.b');
xlim([0.5,2.5]);
ylabel('Fractal Duration(ms)');
set(gca,'FontSize',15,'FontWeight','Bold');
set(gca,'LineWidth',3);
box off;




%Correct rate
axes('Position',[0.1,0.2,0.4,0.3]);
g=bar(SuccessfulRate_Array,'LineWidth',3);

%SuccessfulRate_Str={'Suc_S','FailGaze_S','Failure_S','Suc_L','FailGaze_L','Failure_L'};

ColorChosen= Color_RGB(1:length(SuccessfulRate_Array),:);
g.FaceColor = 'flat';
g.EdgeColor = 'flat';

g.CData = ColorChosen;

xticklabels(SuccessfulRate_Str);

ylabel('SuccessfulRate');
set(gca,'FontSize',15,'FontWeight','Bold');
set(gca,'LineWidth',3);
box off;





%Saccade Latency
axes('Position',[0.55,0.2,0.4,0.3]);
g=bar(RT_Array_Mean,'LineWidth',3);

%SuccessfulRate_Str={'Suc_S','FailGaze_S','Failure_S','Suc_L','FailGaze_L','Failure_L'};

RT_Str={'GS','BS','GL','BL'};

ColorChosen=[1,0,0;0,0,1;1,0,0;0,0,1];
g.FaceColor = 'flat';
g.EdgeColor = 'flat';

g.CData = ColorChosen;

xticklabels(RT_Str);

hold on
b=errorbar(RT_Array_Mean,RT_Array_Sem,'.k','LineWidth',3);

%Add on statistics
%p_Short_Value

%if p_Short_Value<0.05
    sigline([1,2],sprintf('p=%2.3f',p_Short_Value),[],max(RT_Array_Mean+RT_Array_Sem),p_Short_Value);
%
%if p_Long_Value<0.05
    sigline([3,4],sprintf('p=%2.3f',p_Long_Value),[],max(RT_Array_Mean+RT_Array_Sem),p_Long_Value);
%end

ylabel('Reaction Time(ms)');
set(gca,'FontSize',15,'FontWeight','Bold');
set(gca,'LineWidth',3);
box off;








end



end


function organized=ReorganizeEye(data);
organized=[];

for i=1:length(data)
    data_curr=data{i};
    try
    maxnumel=max(cellfun(@numel,data_curr));
    
    organized{i}=cell2mat(cellfun(@(x) [x,NaN*ones(1,maxnumel-length(x))],data_curr,'uniform',0));
    catch
    organized{i}=data_curr;
    end
    
end


 maxnumel=max(cellfun(@(x) size(x,2),organized));
  organized=cell2mat(cellfun(@(x) [x,NaN*ones(size(x,1),maxnumel-size(x,2))],organized,'uniform',0)');

end

function output=ValidPuiple(data,the);

output=cellfun(@(x) addNaN(x,the),data,'uniform',0);


end

function output=addNaN(x,the);
output=x;
output(output<the)=NaN;

end
function data=FewNumCorrection(mat,data_ori,the);
data=data_ori;
colnum=sum(~isnan(mat),1);
inavidcol=colnum<the;
data(inavidcol)=NaN;


end

function [mean,sem]=StatisticsPSTH(psth,the);

mean=nanmean(psth,1);
sem=nanstd(psth,[],1)./sqrt(sum(~isnan(psth),1));

mean=FewNumCorrection(psth,mean,the);
sem=FewNumCorrection(psth,sem,the);


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
