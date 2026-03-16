function EyeMovementAnalysis(Data, StartTrial, EndTrial,ShowFigureFlag,OutputFlag);
close all
%Unpack all the data
ECode;
TaskType=Data.TaskType;
TrialNum=Data.TrialNum;
FixationPointIndex=Data.FixPointIndex;
FixLoc=Data.FixLoc;

TargetIndex=Data.TargetIndex;
TargetLoc=Data.TargetLoc;

GoodObject=mod(TargetIndex,floor(TargetIndex/10))<5;
BadObject=~GoodObject;

EyeChannel_X=Data.EyeChannel_X_L; %Only import the first two channels
EyeChannel_Y=Data.EyeChannel_Y_L;
EyeBinWidth=Data.EyeBinWidth;
EyeTime=Data.EyeTime;

EventChannel=Data.EventChannel;
EventBin=Data.EventBin;

StartOffset=50;%in ms
StopOffset=50;%in ms
StartOffsetBin=StartOffset/EyeBinWidth;
StopOffsetBin=StopOffset/EyeBinWidth;

StartBin=EventBin(find(EventChannel==FIXSTARTCD))-StartOffsetBin;
StopBin=EventBin(find(EventChannel==CORRECTCD  ))+StopOffsetBin;

End_Time=max(StopBin-StartBin);
TimeSeries=(0:EyeBinWidth:End_Time)-StartOffset;

EyeX_in_bin=SelectMatrixData(EyeChannel_X,StartBin,StopBin);
EyeY_in_bin=SelectMatrixData(EyeChannel_Y,StartBin,StopBin);

EyeX_Good=EyeX_in_bin(GoodObject,:);
EyeY_Good=EyeY_in_bin(GoodObject,:);

EyeX_Bad=EyeX_in_bin(BadObject,:);
EyeY_Bad=EyeY_in_bin(BadObject,:);

TimeSeriesPlotGood=repmat(TimeSeries,sum(GoodObject),1)';
TimeSeriesPlotBad=repmat(TimeSeries,sum(BadObject),1)';
if ishandle(101) 
    close 101; end
figure(101);set(101,'color','w');
set(101,'Position', [100,100 1300,500], 'Name', 'Eye Trace');
orient landscape;
axes('position',[0.07 0.57 0.5 0.35]);
plot(TimeSeriesPlotGood,EyeX_Good','-r');
 hold on;plot(TimeSeriesPlotBad,EyeX_Bad','-b');
 plot([0,0],[-20,20],'-.k','LineWidth',1.5)
    xlabel('Time relative to fix off(ms)');
    ylabel('Degree');
    title('Eye Trace Horizontal')
   
axes('position',[0.07 0.1 0.5 0.35]);    
    plot(TimeSeriesPlotGood,EyeY_Good','-r');hold on;
    plot(TimeSeriesPlotBad,EyeY_Bad','-b');
    xlabel('Time relative to fix off(ms)');
    ylabel('Degree');
    title('Eye Trace Vertical')
    
 axes('position',[0.42 0.2 0.7 0.7]);       
    plot(EyeX_Good,EyeY_Good,'.r');axis equal;
    hold on 
    plot(EyeX_Bad,EyeY_Bad,'.b');axis equal;
    xlabel('Eye X(degree)')
    ylabel('Eye Y(degree)')
    xlim([-20,20])
    ylim([-20,20])
    legend('Good Object','Bad Object')
    


keyboard





%keyboard


end