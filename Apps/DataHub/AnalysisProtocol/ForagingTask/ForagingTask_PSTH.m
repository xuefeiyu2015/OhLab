function ForagingTask_PSTH(Data, StartTrial, EndTrial,ShowFigureFlag,OutputFlag);

close all


%%%%%%%%%%%%%%%%%User controller starts here

%Parameter to define the sliding window to measure firing rates
StepSize=5;%in ms
BinWidth=15;%in ms

%Parameter seting for interval selecting plot period
StartOffset=100;%in ms
StopOffset=100;%in ms
IntervalAfterMarkerStart=1000;%in ms;%This for align marker

%For Raster plot 
Barlength=0.5;
LineWidth=1;
%Color Definition for each scenes
Color={'r','g','b','m',[0.4,0.4,0],[0.28,0.78,0.89],'k',[0.6,0.1,1],[0.6,0.1,0.4],[0.3,0.5,1],[0.6,0.5,0.1],[0,0.8,0.4]};

%Plot Scene flag
PlotScene=1;%Value equals 1 indicates figure 101 will show scene PSTH image; 
%Plot Target flag
PlotFractal=1;%Value equals 1 indicates figure 102 will show fractal PSTH image; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Unpack all the data base and task settings
ECode;
TaskType=Data.TaskType;
TrialNum=EndTrial-StartTrial+1;
Selection=StartTrial:EndTrial;

FixationPointIndex=Data.FixPointIndex(Selection,:);
FixLoc=Data.FixLoc(Selection,:,:);





TargetIndex=Data.TargetIndex(Selection,:,:);
TargetLoc=Data.TargetLoc(Selection,:,:);

TargetNum=Data.TargetNum(Selection,:);
TargetSize=Data.TargetSize(Selection,:);
TargetPolarLoc=Data.TargetPolarLoc(Selection,:,:);
ChoiceOutcome=Data.ChoiceOutcome(Selection,:);


GoodObjectLib=Data.GoodObjectLib;
BadObjectLib=Data.BadObjectLib;
SceneLib=Data.SceneLib;

GoodObject=arrayfun(@(x)( ismember(x,GoodObjectLib(:,3))),TargetIndex(:,2:end,2));%Use real name
%GoodObjectID=TargetIndex(:,2:end,2);

BadObject=arrayfun(@(x)( ismember(x,BadObjectLib(:,3))),TargetIndex(:,2:end,2));%Use real name
%BadObjectID=TargetIndex(:,2:end,2);
%GoodBadIndex=GoodObject*1+BadObject*2;%1:Good Object;2:Bad Object

SceneIndex=Data.SceneIndex(Selection,:);

%EyeChannel_X=Data.EyeChannel_X_L(Selection,:); %Only import the first two channels
%EyeChannel_Y=Data.EyeChannel_Y_L(Selection,:);
EyeBinWidth=Data.EyeBinWidth;
%EyeTime=Data.EyeTime(Selection,:);

EventChannel=Data.EventChannel(Selection,:);
EventTimeChannel=Data.EventTimeChannel(Selection,:);
EventBin=Data.EventBin(Selection,:);

SpikeTime=squeeze(Data.SpikeTimeData(:,:,Selection))';
SpikeBin=squeeze(Data.SpikeBinData(:,:,Selection))';



StartOffsetBin=StartOffset/EyeBinWidth;
StopOffsetBin=StopOffset/EyeBinWidth;
IntervalInBin=IntervalAfterMarkerStart/EyeBinWidth;


 %Align with the scene onset
AlignTime=FindOutTime(EventChannel,EventTimeChannel,SCENE_ON);
SelectIntervalTime=repmat([-StartOffset,IntervalAfterMarkerStart+StopOffset],TrialNum,1);
MarkerTime=FindOutTime(EventChannel,EventTimeChannel,TGTCD);
FirstTargetOnset=MarkerTime(:,1);
FixOnsetTime=FindOutTime(EventChannel,EventTimeChannel,SHOWFIXCD);
FixOnsetTime=FixOnsetTime(:,1);

SpikeTimeAfterAlign=SpikeTime-repmat(AlignTime,1,size(SpikeTime,2));
%EyeTimeAfterAlign=EyeTime-repmat(AlignTime,1,size(EyeTime,2));
EventTimeAfterAlign=EventTimeChannel-repmat(AlignTime,1,size(EventTimeChannel,2));
FirstTargetOnsetAfterAlign=FirstTargetOnset-AlignTime;
FixOnsetTimeAfterAlign=FixOnsetTime-AlignTime;

%Another way of align: Align with the target onset
AlignTimeWithTargetOnset=FindOutTime(EventChannel,EventTimeChannel,TGTCD);
%TargetOffRelativeToOnset=FindOutTime(EventChannel(:,5:end),EventTimeChannel(:,5:end),TGTOFF)-AlignTimeWithTargetOnset;%Start searching after visual onset
TargetOffRelativeToOnset=500+0*AlignTimeWithTargetOnset;%Start searching after visual onset
TargetOffRelativeTime=TargetOffRelativeToOnset+StopOffset;






%#####################Parameter measuring the firing rate(PSTH)
%Using a sliding window to calculate firing rate


[SpikeFiringRate,TimeSequence]=CalculateSpikeFiringRate(SpikeTimeAfterAlign,StepSize,BinWidth,SelectIntervalTime);

SpikeTimeInInterval=SelectSpikeTimeInterval(SpikeTimeAfterAlign,SelectIntervalTime(:,1),SelectIntervalTime(:,2));%For Rasterplot

%################################################
%Calculate firing rate aligned on target onset
MaxTgtTime=0;
MinTgtTime=1000;
for i=1:max(TargetNum(:,2))
    AlignTimeForTargetN=AlignTimeWithTargetOnset(:,i);
    SpikeTimeAlignTgtN=SpikeTime-repmat(AlignTimeForTargetN,1,size(SpikeTime,2));
    SelectIntervalTimeTgt=[repmat(-StartOffset,size(AlignTimeForTargetN,1),1),TargetOffRelativeTime(:,i)];
    
    
    [FRAlignTgtN{i},TimeSeqTgtN{i}]=CalculateSpikeFiringRate(SpikeTimeAlignTgtN,StepSize,BinWidth,SelectIntervalTimeTgt);
    SpikeTimeTgtNInInterval{i}=SelectSpikeTimeInterval(SpikeTimeAlignTgtN,SelectIntervalTimeTgt(:,1),SelectIntervalTimeTgt(:,2));%For Rasterplot
    
    MaxTgtTime=nanmax(MaxTgtTime, nanmax(SelectIntervalTimeTgt(:,2)));
    MinTgtTime=nanmin(MinTgtTime, nanmin(SelectIntervalTimeTgt(:,1)));
    
    
    %currentLine=GoodBadIndex(:,i);
    %GoodBadIndex;
    
end



%#################Separate trials according to scene name
UniqueScene=unique(SceneIndex(:,3));
NumOfCondition=length(UniqueScene);
%Find out the max firing rate for target onset align spikes

for i=1:NumOfCondition
    currentScene=UniqueScene(i);
    trialindex=SceneIndex(:,3)==currentScene;
    currentspikerate(i,:)=nanmean(SpikeFiringRate(trialindex,:));
    currentspikerate_sem(i,:)=nanstd(SpikeFiringRate(trialindex,:))/sqrt( sum(trialindex));
    
    currentspiketime{i}=SpikeTimeInInterval(trialindex,:);
    currentTime(i,:)=nanmean(TimeSequence(trialindex,:));
    currentFirstTargetOnsetAfterAlign(i,:)=[min(FirstTargetOnsetAfterAlign(trialindex,:)),max(FirstTargetOnsetAfterAlign(trialindex,:))];
    
    currentFixOnsetTimeAfterAlign(i,:)=[min(FixOnsetTimeAfterAlign(trialindex,:)),max(FixOnsetTimeAfterAlign(trialindex,:))];
   
    
    
    %SpikeTime for Raster
    ObjForScene_temp=cellfun(@(x)(x(trialindex,:)),SpikeTimeTgtNInInterval,'UniformOutput', false);
    ObjForScene_FiringRate_temp=cellfun(@(x)(x(trialindex,:)),FRAlignTgtN,'UniformOutput', false);
    ObjForScene_FiringRate_TimeSeq_temp=cellfun(@(x)(x(trialindex,:)),TimeSeqTgtN,'UniformOutput', false);
    
    GoodObjForScene=NaN*ones(max(TargetNum(:,2))*sum(trialindex),1000);
   
    %Select Good Object for current scene
    GoodObjCurrent=GoodObject(trialindex,:);
    
    MaxRowG=1;
    MaxColG=1;
    startrowG=1;
    %Firing rate for PSTH
    GoodObjForSceneFiringRate=GoodObjForScene;
    GoodObjForSceneFiringRateTime=GoodObjForScene;
    MaxColG_FR=1;
    
    
    
    %Select Bad Object for current scene
    
    BadObjCurrent=BadObject(trialindex,:);
    BadObjForScene=NaN*ones(max(TargetNum(:,2))*sum(trialindex),1000);
    
    MaxRowB=1;
    MaxColB=1;
    startrowB=1;
    %Firing rate for PSTH
    BadObjForSceneFiringRate=BadObjForScene;
    BadObjForSceneFiringRateTime=BadObjForScene;
    MaxColB_FR=1;
    
    
    for j=1:max(TargetNum(:,2))
        %Select spike time
        tempG=ObjForScene_temp{j}(GoodObjCurrent(:,j),:);
        GoodObjForScene(startrowG:startrowG+size(tempG,1)-1,1:size(tempG,2))=tempG;
        
        if MaxColG<size(tempG,2)
            MaxColG=size(tempG,2);
        end
        
        
        
        tempB=ObjForScene_temp{j}(BadObjCurrent(:,j),:);
        BadObjForScene(startrowB:startrowB+size(tempB,1)-1,1:size(tempB,2))=tempB;
        
        if MaxColB<size(tempB,2)
            MaxColB=size(tempB,2);
        end
       
       
        %Select spike rate
        tempG_Rate=ObjForScene_FiringRate_temp{j}(GoodObjCurrent(:,j),:);
        tempG_Rate_TimeSeq=ObjForScene_FiringRate_TimeSeq_temp{j}(GoodObjCurrent(:,j),:);
        GoodObjForSceneFiringRate(startrowG:startrowG+size(tempG_Rate,1)-1,1:size(tempG_Rate,2))=tempG_Rate;
        GoodObjForSceneFiringRateTime(startrowG:startrowG+size(tempG_Rate_TimeSeq,1)-1,1:size(tempG_Rate_TimeSeq,2))=tempG_Rate_TimeSeq;
        
        
        
        
        
        %TimeSeqTgtN{j}
         if MaxColG_FR<size(tempG_Rate,2)
            MaxColG_FR=size(tempG_Rate,2);
         end
         
        tempB_Rate=ObjForScene_FiringRate_temp{j}(BadObjCurrent(:,j),:);
        tempB_Rate_TimeSeq=ObjForScene_FiringRate_TimeSeq_temp{j}(BadObjCurrent(:,j),:);
        BadObjForSceneFiringRate(startrowB:startrowB+size(tempB_Rate,1)-1,1:size(tempB_Rate,2))=tempB_Rate;
        BadObjForSceneFiringRateTime(startrowB:startrowB+size(tempB_Rate_TimeSeq,1)-1,1:size(tempB_Rate_TimeSeq,2))=tempB_Rate_TimeSeq;
      
        
        
        if MaxColB_FR<size(tempB_Rate,2)
           MaxColB_FR=size(tempB_Rate,2);
        end
        
         startrowG=startrowG+size(tempG,1);
          MaxRowG=startrowG-1;
        
         startrowB=startrowB+size(tempB,1);
         MaxRowB=startrowB-1;
    end
    GoodObjForSceneN_SpikeTime{i}=GoodObjForScene(1:MaxRowG,1:MaxColG);
    BadObjForSceneN_SpikeTime{i}=BadObjForScene(1:MaxRowB,1:MaxColB);
    
    GoodObjForSceneN_FiringRate{i}=GoodObjForSceneFiringRate(1:MaxRowG,1:MaxColG_FR);
    BadObjForSceneN_FiringRate{i}=BadObjForSceneFiringRate(1:MaxRowB,1:MaxColB_FR);
    GoodObjForSceneN_FiringRateTime{i}= GoodObjForSceneFiringRateTime(1:MaxRowG,1:MaxColG_FR);
    BadObjForSceneN_FiringRateTime{i}=BadObjForSceneFiringRateTime(1:MaxRowB,1:MaxColB_FR);
    
   
    
    
    try
    GoodObjIndexForSceneN(i)=unique(GoodObjectLib(SceneLib(:,3)==currentScene,3));
    BadObjIndexForSceneN(i)=unique(BadObjectLib(SceneLib(:,3)==currentScene,3));
    catch
        GoodObjIndexForSceneN(i)=NaN;
        BadObjIndexForSceneN(i)=NaN;
    end
   
  
    
    
end



if(PlotScene)
%Visualize the PSTH according to the scenes
FigureColNum=4;%4Column
FigureRowNum=ceil(NumOfCondition/FigureColNum);

Color={'r','g','b','m',[0.4,0.4,0],[0.28,0.78,0.89],'k',[0.6,0.1,1],[0.6,0.1,0.4],[0.3,0.5,1],[0.6,0.5,0.1],[0,0.8,0.4],...
    [0, 0.4470, 0.7410],[0.8500, 0.3250, 0.0980],[0.9290, 0.6940, 0.1250],...
    [0.4940, 0.1840, 0.5560],[0.4660, 0.6740, 0.1880],...
    [0.3010, 0.7450, 0.9330],[0.6350, 0.0780, 0.1840],...
    [0.5678, 0.1450, 0.3330],[0.0010, 0.1450, 0.9530],...
    [0.9010, 0.5450, 0.3330],[0.5010, 0.3450, 0.5330],...
    [0.44,0.13,1],[1.0,0.2,0.77]
    };





if ishandle(101) 
    close 101; end
figure(101);set(101,'color','w');
set(101,'Position', [100,100 1300,1000], 'Name', 'PSTH Align on scene onset');
pannelcounter=1;
conditionSelector=1;
x_offset=5;
y_offset=10;
StartLoc=1;
xrange=[min(SelectIntervalTime(:,1))-x_offset,max(SelectIntervalTime(:,2))+x_offset];
yrange=[max(nanmin(nanmin(currentspikerate))-y_offset,0),nanmax(nanmax(currentspikerate))+y_offset];

while conditionSelector<=NumOfCondition

        SpikeTimePannel=currentspiketime{conditionSelector};
 
         subplot(FigureRowNum*2,FigureColNum,pannelcounter); 
        SpikeRaterPlot(SpikeTimePannel,StartLoc,Barlength,Color{conditionSelector},LineWidth);

        
        hold on
        plot([0,0],[Barlength,size(SpikeTimePannel,1)],'-.k','LineWidth',1.5);
        box off
        set(gca,'LineWidth',2)
        set(gca,'FontWeight','Bold')
        xlim(xrange);
        %Plot target on line
        plot([currentFirstTargetOnsetAfterAlign(conditionSelector,1),currentFirstTargetOnsetAfterAlign(conditionSelector,1)],yrange,'-.k','LineWidth',1)
        plot([currentFirstTargetOnsetAfterAlign(conditionSelector,2),currentFirstTargetOnsetAfterAlign(conditionSelector,2)],yrange,'-.k','LineWidth',1)
        %Plot the fixation onset line
        plot([currentFixOnsetTimeAfterAlign(conditionSelector,1), currentFixOnsetTimeAfterAlign(conditionSelector,1)],yrange,'-.k','LineWidth',1)
        plot([ currentFixOnsetTimeAfterAlign(conditionSelector,2), currentFixOnsetTimeAfterAlign(conditionSelector,2)],yrange,'-.k','LineWidth',1)
      
        
        ss=sprintf('Scene %d',UniqueScene(conditionSelector));
        title(ss);
         if (mod(pannelcounter,FigureColNum)==1)
          ylabel({'TrialNum',''},'FontSize',10,'FontWeight','Bold')
         end
 
        subplot(FigureRowNum*2,FigureColNum,pannelcounter+FigureColNum); 
        shadedErrorBar(currentTime(conditionSelector,:),currentspikerate(conditionSelector,:),currentspikerate_sem(conditionSelector,:),'lineprops',{Color{conditionSelector}},'transparent',1);
        
      
        hold on
        plot([0,0],yrange,'-.k','LineWidth',1.5);
        box off
        set(gca,'LineWidth',2)
        set(gca,'FontWeight','Bold')
        xlim(xrange);
        ylim(yrange);
        %Plot target on line
        plot([currentFirstTargetOnsetAfterAlign(conditionSelector,1),currentFirstTargetOnsetAfterAlign(conditionSelector,1)],yrange,'-.k','LineWidth',1)
        plot([currentFirstTargetOnsetAfterAlign(conditionSelector,2),currentFirstTargetOnsetAfterAlign(conditionSelector,2)],yrange,'-.k','LineWidth',1)
         %Plot the fixation onset line
        plot([currentFixOnsetTimeAfterAlign(conditionSelector,1), currentFixOnsetTimeAfterAlign(conditionSelector,1)],yrange,'-.k','LineWidth',1)
        plot([ currentFixOnsetTimeAfterAlign(conditionSelector,2), currentFixOnsetTimeAfterAlign(conditionSelector,2)],yrange,'-.k','LineWidth',1)
        if (mod((pannelcounter+FigureColNum),FigureColNum)==1)
            ylabel({'Firing Rate','(Hz)'},'FontSize',10,'FontWeight','Bold')
        end
        
        if(ceil((pannelcounter+FigureColNum)/FigureColNum)==(FigureRowNum*2))
           xlabel('Time relative to scene onset.(ms)','FontSize',10,'FontWeight','Bold')
        end
        
        pannelcounter=pannelcounter+1;
        if (mod(pannelcounter,FigureColNum)==1 && pannelcounter>1)
            
                pannelcounter=pannelcounter+FigureColNum;
                
           
        end
       
        conditionSelector=conditionSelector+1;
end
end %End of the scene plot flag

%Visualize the PSTH align on target onset

 if(PlotFractal)
FigureColNum=4;%4Column
FigureRowNum=ceil(NumOfCondition/FigureColNum);

if ishandle(102) 
    close 102; end
figure(102);set(102,'color','w');
set(102,'Position', [2,2, 1300,1000], 'Name', 'PSTH Align on fractal onset');
pannelcounter=1;
conditionSelector=1;
x_offset=5;
y_offset=10;
xrange=[MinTgtTime-x_offset,MaxTgtTime+x_offset];

MeanFiringRateG1=cellfun(@(x) nanmean(x),GoodObjForSceneN_FiringRate,'UniformOutput',false);
MaxFiringRateG2=max(cell2mat(cellfun(@(x) max(x),MeanFiringRateG1,'UniformOutput',false)));
MeanFiringRateB1=cellfun(@(x) nanmean(x),BadObjForSceneN_FiringRate,'UniformOutput',false);
MaxFiringRateB2=max(cell2mat(cellfun(@(x) max(x),MeanFiringRateB1,'UniformOutput',false)));
MinFiringRateG2=max(cell2mat(cellfun(@(x) min(x),MeanFiringRateG1,'UniformOutput',false)));
MinFiringRateB2=max(cell2mat(cellfun(@(x) min(x),MeanFiringRateB1,'UniformOutput',false)));

MaxFiringRate=max(MaxFiringRateG2,MaxFiringRateB2);
MinFiringRate=min(MinFiringRateG2,MinFiringRateB2);


yrange=[max(MinFiringRate-y_offset,0),MaxFiringRate+y_offset];

RasterWidth=1/FigureColNum*0.9;
RasterHeight=1/(FigureRowNum*6);
PSTHWidth=1/FigureColNum*0.9;
PSTHHeight=1/(FigureRowNum*5);
PlotStartLocX=0.04;
PlotStartLocY=0.98;
iter_termY=0;

while conditionSelector<=NumOfCondition
    
   

    

        GoodSpikeTimePannel=GoodObjForSceneN_SpikeTime{conditionSelector};
        BadSpikeTimePannel=BadObjForSceneN_SpikeTime{conditionSelector};
        GoodFiringRatePannel=GoodObjForSceneN_FiringRate{conditionSelector};
        BadFiringRatePannel=BadObjForSceneN_FiringRate{conditionSelector};
        GoodObj_FiringRateTimePanel= nanmean(GoodObjForSceneN_FiringRateTime{conditionSelector},1);
        BadObj_FiringRateTimePanel=nanmean(BadObjForSceneN_FiringRateTime{conditionSelector},1);
        
        currGoodObject=GoodObjIndexForSceneN(i);
        currBadObject=BadObjIndexForSceneN(i);
        %Plot good object PSTH
        if mod(pannelcounter,FigureColNum)==0
            iter_term=FigureColNum;
        else
            iter_term=mod(pannelcounter,FigureColNum);
         
        end
        
         if mod(pannelcounter,FigureColNum)==1
             %{
             if mod(pannelcounter,FigureRowNum)==0
                 iter_termY=FigureRowNum;
                 
             else
                  iter_termY=mod(pannelcounter,FigureRowNum);       
                 
             end
             %}
             if mod(pannelcounter,2*FigureColNum)==1
                 iter_termY=iter_termY+1;
                 
             end
   
         end
       
        
        axes('position',[PlotStartLocX+PSTHWidth*(iter_term-1)*1.05 PlotStartLocY-RasterHeight-(0.03+RasterHeight+PSTHHeight+RasterHeight*1.15+PSTHHeight*1.05)*(iter_termY-1)*1.15, RasterWidth RasterHeight]);

        %subplot(FigureRowNum*4,FigureColNum,pannelcounter); 
        SpikeRaterPlot(GoodSpikeTimePannel,StartLoc,Barlength,Color{conditionSelector},LineWidth);
        
        hold on
        plot([0,0],[Barlength,size(GoodSpikeTimePannel,1)],'-.k','LineWidth',1.5);
        box off
        set(gca,'LineWidth',2)
        set(gca,'FontWeight','Bold')
        xlim(xrange);
        set(gca,'xtick',[])
        %{
        %Plot target on line
        plot([currentFirstTargetOnsetAfterAlign(conditionSelector,1),currentFirstTargetOnsetAfterAlign(conditionSelector,1)],yrange,'-.k','LineWidth',1)
        plot([currentFirstTargetOnsetAfterAlign(conditionSelector,2),currentFirstTargetOnsetAfterAlign(conditionSelector,2)],yrange,'-.k','LineWidth',1)
        %Plot the fixation onset line
        plot([currentFixOnsetTimeAfterAlign(conditionSelector,1), currentFixOnsetTimeAfterAlign(conditionSelector,1)],yrange,'-.k','LineWidth',1)
        plot([ currentFixOnsetTimeAfterAlign(conditionSelector,2), currentFixOnsetTimeAfterAlign(conditionSelector,2)],yrange,'-.k','LineWidth',1)
      %}
        
        ss=sprintf('Good Fractal %d, Scene %d',currGoodObject,UniqueScene(conditionSelector));
        title(ss);
         if (mod(pannelcounter,FigureColNum)==1)
          ylabel({'TrialNum',''},'FontSize',10,'FontWeight','Bold')
         end
         axes('position',[PlotStartLocX+PSTHWidth*(iter_term-1)*1.05 PlotStartLocY-(RasterHeight+PSTHHeight)-0.01-(0.03+RasterHeight+PSTHHeight+RasterHeight*1.15+PSTHHeight*1.05)*(iter_termY-1)*1.15, PSTHWidth PSTHHeight]);
        %subplot(FigureRowNum*4,FigureColNum,pannelcounter+FigureColNum); 
       
        shadedErrorBar(GoodObj_FiringRateTimePanel,nanmean(GoodFiringRatePannel,1),nanstd(GoodFiringRatePannel,1)./sum(~isnan(GoodFiringRatePannel),1),'lineprops',{Color{conditionSelector}},'transparent',1);
        
      
        hold on
        plot([0,0],yrange,'-.k','LineWidth',1.5);
        box off
        set(gca,'LineWidth',2)
        set(gca,'FontWeight','Bold')
        xlim(xrange);
        ylim(yrange);
        %{
        %Plot target on line
        plot([currentFirstTargetOnsetAfterAlign(conditionSelector,1),currentFirstTargetOnsetAfterAlign(conditionSelector,1)],yrange,'-.k','LineWidth',1)
        plot([currentFirstTargetOnsetAfterAlign(conditionSelector,2),currentFirstTargetOnsetAfterAlign(conditionSelector,2)],yrange,'-.k','LineWidth',1)
         %Plot the fixation onset line
        plot([currentFixOnsetTimeAfterAlign(conditionSelector,1), currentFixOnsetTimeAfterAlign(conditionSelector,1)],yrange,'-.k','LineWidth',1)
        plot([ currentFixOnsetTimeAfterAlign(conditionSelector,2), currentFixOnsetTimeAfterAlign(conditionSelector,2)],yrange,'-.k','LineWidth',1)
        %}
        if (mod((pannelcounter+FigureColNum),FigureColNum)==1)
            ylabel({'Firing Rate','(Hz)'},'FontSize',10,'FontWeight','Bold')
        end
        
        
        %Plot Bad object PSTH
        
        axes('position',[PlotStartLocX+PSTHWidth*(iter_term-1)*1.05 PlotStartLocY-(0.03+RasterHeight+PSTHHeight+RasterHeight*1.15)-(0.03+RasterHeight+PSTHHeight+RasterHeight*1.15+PSTHHeight*1.05)*(iter_termY-1)*1.15, RasterWidth RasterHeight]);
        %subplot(FigureRowNum*4,FigureColNum,pannelcounter+2*FigureColNum); 
        SpikeRaterPlot(BadSpikeTimePannel,StartLoc,Barlength,Color{conditionSelector},LineWidth);
        
        hold on
        plot([0,0],[Barlength,size(BadSpikeTimePannel,1)],'-.k','LineWidth',1.5);
        box off
        set(gca,'LineWidth',2)
        set(gca,'FontWeight','Bold')
        xlim(xrange);
        set(gca,'xtick',[]);
        %{
        %Plot target on line
        plot([currentFirstTargetOnsetAfterAlign(conditionSelector,1),currentFirstTargetOnsetAfterAlign(conditionSelector,1)],yrange,'-.k','LineWidth',1)
        plot([currentFirstTargetOnsetAfterAlign(conditionSelector,2),currentFirstTargetOnsetAfterAlign(conditionSelector,2)],yrange,'-.k','LineWidth',1)
        %Plot the fixation onset line
        plot([currentFixOnsetTimeAfterAlign(conditionSelector,1), currentFixOnsetTimeAfterAlign(conditionSelector,1)],yrange,'-.k','LineWidth',1)
        plot([ currentFixOnsetTimeAfterAlign(conditionSelector,2), currentFixOnsetTimeAfterAlign(conditionSelector,2)],yrange,'-.k','LineWidth',1)
      %}
        
        ss=sprintf('Bad Fractal %d, Scene %d',currBadObject,UniqueScene(conditionSelector));
        title(ss);
         if (mod(pannelcounter,FigureColNum)==1)
          ylabel({'TrialNum',''},'FontSize',10,'FontWeight','Bold')
         end
        axes('position',[PlotStartLocX+PSTHWidth*(iter_term-1)*1.05 PlotStartLocY-(0.03+RasterHeight+PSTHHeight+RasterHeight*1.15+PSTHHeight*1.05)-(0.03+RasterHeight+PSTHHeight+RasterHeight*1.15+PSTHHeight*1.05)*(iter_termY-1)*1.15, PSTHWidth PSTHHeight]);
        %subplot(FigureRowNum*4,FigureColNum,pannelcounter+3*FigureColNum); 
      
        shadedErrorBar(BadObj_FiringRateTimePanel,nanmean(BadFiringRatePannel,1),nanstd(BadFiringRatePannel,1)./sum(~isnan(BadFiringRatePannel),1),'lineprops',{Color{conditionSelector}},'transparent',1);
       
      
        hold on
        plot([0,0],yrange,'-.k','LineWidth',1.5);
        box off
        set(gca,'LineWidth',2)
        set(gca,'FontWeight','Bold')
        xlim(xrange);
        ylim(yrange);
        
        
        if (mod((pannelcounter+FigureColNum),FigureColNum)==1)
            ylabel({'Firing Rate','(Hz)'},'FontSize',10,'FontWeight','Bold')
        end
 
        
        if(NumOfCondition-conditionSelector<=FigureColNum)
           xlabel('Time relative to fractal onset.(ms)','FontSize',10,'FontWeight','Bold')
        end
        
        pannelcounter=pannelcounter+1;
        if (mod(pannelcounter,FigureColNum)==1 && pannelcounter>1)
            
                pannelcounter=pannelcounter+FigureColNum*3;
                
           
        end
       
        conditionSelector=conditionSelector+1;
end


end %End of the target plot flag

   

    

end





