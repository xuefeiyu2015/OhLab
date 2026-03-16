function ITI_RandomSeqFixSeqCompare(Data, StartTrial, EndTrial,ShowFigureFlag,OutputFlag);

%This is to compare the random sequence block and the fix sequence block
%neuron response during ITI
%Xuefei Yu, 03032020
%Plot Interval
%Before Scene Onset
default=0;%If default ==1 , then it will plot the whole period
xlim_BeforeSceneOn=[-2000,500];
xlim_AfterSceneOn=[-500,2000];
xlim_AfterSceneOff=[-500,2000];



%Parameter to define the sliding window to measure firing rates
StepSize=10;%in ms
BinWidth=50;%in ms

%Parameter seting for interval selecting plot period
StartOffset=500;%in ms
StopOffset=50;%in ms
%For scene offset
StartOffset2=500;%in ms
StopOffset2=50;%in ms
IntervalAfterMarkerStart=1000;%in ms;%This for align marker
IntervalAfterMarkerStartOff=1000;%in ms;%For scene offset
%{
%For Raster plot 
Barlength=0.5;
LineWidth=1;
%}
%Color Definition for each scenes

%ColorFix=[lines(6);prism(6);pink(6)];
%ColorFix(2,:)=copper(1)+0.5;
%ColorRand=flipud(ColorFix);
ColorFix=repmat([1,0,0],26,1);
ColorRand=repmat([0,0,1],26,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Unpack all the data base and task settings
ECode;
FileName=Data.FileName;
TaskType=Data.TaskType;
TrialNum=EndTrial-StartTrial+1;
Selection=StartTrial:EndTrial;

%FixationPointIndex=Data.FixPointIndex(Selection,:);
%FixLoc=Data.FixLoc(Selection,:,:);


TargetIndex=Data.TargetIndex(Selection,:,:);
BlockType=Data.BlockType(Selection,:,:);%0:fixed sequence;1;random sequence; 0:ScenesetA;1:SceneSetB
BlockType(BlockType(:,2)>2,:)=999;




%TargetLoc=Data.TargetLoc(Selection,:,:);

%TargetNum=Data.TargetNum(Selection,:);
%TargetSize=Data.TargetSize(Selection,:);
%TargetPolarLoc=Data.TargetPolarLoc(Selection,:,:);
%ChoiceOutcome=Data.ChoiceOutcome(Selection,:);


GoodObjectLib=Data.GoodObjectLib;
BadObjectLib=Data.BadObjectLib;
%SceneLib=Data.SceneLib;

%GoodObject=arrayfun(@(x)( ismember(x,GoodObjectLib(:,3))),TargetIndex(:,2:end,2));%Use real name
%GoodObjectID=TargetIndex(:,2:end,2);

%BadObject=arrayfun(@(x)( ismember(x,BadObjectLib(:,3))),TargetIndex(:,2:end,2));%Use real name
%BadObjectID=TargetIndex(:,2:end,2);
%GoodBadIndex=GoodObject*1+BadObject*2;%1:Good Object;2:Bad Object

SceneIndex=Data.SceneIndex(Selection,:);

EyeChannel_X=Data.EyeChannel_X_L(Selection,:); %Only import the first two channels
EyeChannel_Y=Data.EyeChannel_Y_L(Selection,:);
EyeBinWidth=Data.EyeBinWidth;
EyeTime=Data.EyeTime(Selection,:);

EventChannel=Data.EventChannel(Selection,:);
EventTimeChannel=Data.EventTimeChannel(Selection,:);
EventBin=Data.EventBin(Selection,:);

SpikeChannel=Data.SpikeChannel;
SpikeTime=squeeze(Data.SpikeTimeData(:,:,Selection))';
SpikeBin=squeeze(Data.SpikeBinData(:,:,Selection))';

 %Align with the scene onset
AlignTime=FindOutTime(EventChannel,EventTimeChannel,SCENE_ON);

AlignTimeOff=FindOutTime(EventChannel,EventTimeChannel,SCENE_OFF);

FixOnTime=FindOutTime(EventChannel,EventTimeChannel,SHOWFIXCD )-AlignTime;
FixOnTime=FixOnTime(:,1);

TgtOnTime=FindOutTime(EventChannel,EventTimeChannel,TGTCD  )-AlignTime;
TgtOnTime=TgtOnTime(:,1);
RewardTime=FindOutTime(EventChannel,EventTimeChannel,RWDCD )-AlignTimeOff;
RewardTimeOff=FindOutTime(EventChannel,EventTimeChannel,RWDOFFCD)-AlignTimeOff;
IntervalForFR=min(AlignTimeOff-AlignTime);%in ms; Use the minimum time interval
%IntervalAfterMarkerStart=IntervalAfterMarkerStart

AverageFRTime=repmat([0,IntervalForFR],TrialNum,1);


%Findout the start and the end of the whole task period

StartTime=EyeTime(:,1);

temp=~isnan(EyeChannel_X);  
Indices=arrayfun(@(x) find(temp( x,:), 1, 'last'), 1:size(EyeChannel_X, 1));
EndTime=arrayfun(@(x,y) EyeTime(x,y),1:size(EyeTime, 1), Indices)';

%SelectIntervalTime=repmat([-StartOffset,IntervalAfterMarkerStart+StopOffset],TrialNum,1);
SelectIntervalTime=[StartTime-AlignTime,repmat([IntervalAfterMarkerStart+StopOffset],TrialNum,1)];
SelectIntervalTimeOff=[repmat([-StartOffset2],TrialNum,1),EndTime-AlignTimeOff];

SelectIntervalTimeAfter=[StartTime-AlignTime,AlignTimeOff+repmat([StopOffset],TrialNum,1)-AlignTime];
%MarkerTime=FindOutTime(EventChannel,EventTimeChannel,TGTCD);
%FirstTargetOnset=MarkerTime(:,1);
%FixOnsetTime=FindOutTime(EventChannel,EventTimeChannel,SHOWFIXCD);
%FixOnsetTime=FixOnsetTime(:,1);

SpikeTimeAfterAlign=SpikeTime-repmat(AlignTime,1,size(SpikeTime,2));
SpikeTimeAfterAlignOff=SpikeTime-repmat(AlignTimeOff,1,size(AlignTimeOff,2));

EyeTimeAfterAlign=EyeTime-repmat(AlignTime,1,size(EyeTime,2));
EventTimeAfterAlign=EventTimeChannel-repmat(AlignTime,1,size(EventTimeChannel,2));

%#####################Parameter measuring the firing rate(PSTH)
%Using a sliding window to calculate firing rate


[SpikeFiringRate,TimeSequence]=CalculateSpikeFiringRate(SpikeTimeAfterAlign,StepSize,BinWidth,SelectIntervalTime);
[SpikeFiringRateOff,TimeSequenceOff]=CalculateSpikeFiringRate(SpikeTimeAfterAlignOff,StepSize,BinWidth,SelectIntervalTimeOff);

[SpikeFiringRateAfter,TimeSequenceAfter]=CalculateSpikeFiringRate(SpikeTimeAfterAlign,StepSize,BinWidth,SelectIntervalTimeAfter);

SpikeFiringRateOn=[nan*ones(size(SpikeFiringRateOff)),SpikeFiringRate];
%SpikeFiringRateOn(2:end,:)

SpikeFiringRateOn(2:end,1:size(SpikeFiringRateOff,2))=ReShapeNaN(SpikeFiringRateOff(1:end-1,:),1);

SpikeFiringRateOff=[SpikeFiringRateOff,nan*ones(size(SpikeFiringRate))];

SpikeFiringRateOff(1:end-1,end-size(SpikeFiringRate,2)+1:end)=SpikeFiringRate(1:end-1,:);

SpikeFiringRateOff=ReShapeNaN(SpikeFiringRateOff,-1);



SpikeFiringRate=SpikeFiringRateOn;





SpikeTimeInInterval=SelectSpikeTimeInterval(SpikeTimeAfterAlign,AverageFRTime(:,1),AverageFRTime(:,2));%Return spike count to calculate average firing rate
SpikeCountInInterval=sum(~isnan(SpikeTimeInInterval),2);
AverageFiringRateInInterval=SpikeCountInInterval./(IntervalAfterMarkerStart/1000);%in Hz
TimeSequenceOn_Ori=nanmean(TimeSequence);
TimeSequenceOff_Ori=nanmean(TimeSequenceOff);

TimeSequenceAfter=nanmean(TimeSequenceAfter);

TimeSequenceOn=[TimeSequenceOff_Ori-TimeSequenceOff_Ori(end)+TimeSequenceOn_Ori(1)-StepSize TimeSequenceOn_Ori];
TimeSequenceOff=[TimeSequenceOff_Ori TimeSequenceOn_Ori-TimeSequenceOn_Ori(1)+TimeSequenceOff_Ori(end)+StepSize];

%#################Separate trials according to scene name and block type

UniqueBlockType=unique(BlockType(:,2:end),'row');
NumberOfBlock=size(UniqueBlockType,1);
UniqueSceneNameTotal=unique(SceneIndex(:,3));
AxisNum=[UniqueSceneNameTotal,[1:length(UniqueSceneNameTotal)]'];
if ishandle(130)
close 130;end
figure(130);set(130,'color','w');
set(130,'Position', [100,100, 1000,500], 'Name', 'Before Scene Onset Response Comparison between fixed sequence blocks and random sequence blocks');

if ishandle(131)
close 131;end
figure(131);set(131,'color','w');
set(131,'Position', [100,100, 1000,500], 'Name', 'After Scene Offset Response Comparison between fixed sequence blocks and random sequence blocks');

if ishandle(132)
close 132;end
figure(132);set(132,'color','w');
set(132,'Position', [100,100, 1000,500], 'Name', 'After Scene Onset Response Comparison ');


for i=1:NumberOfBlock
    currBlockType=UniqueBlockType(i,:);
    SceneNameForCurrentBlock=SceneIndex(sum(BlockType(:,2:end)==currBlockType,2)==2,3);
 
    UniqueScene=unique(SceneNameForCurrentBlock);
    for j=1:length(UniqueScene)
        %Calculate and plot the average firing rate according to scenes

       currentScene=UniqueScene(j);
       trialindex=(sum(BlockType(:,2:end)==currBlockType,2)==2)&(SceneIndex(:,3)==currentScene);
       %{
       %For Scene Tuning
       FR_Scene=AverageFiringRateInInterval(trialindex);
  
       FR_Scene_mean(i)=mean(FR_Scene);
       FR_Scene_sem(i)=std(FR_Scene)/sqrt(sum(trialindex));
       FR_Scene_SSE(i)=sum((FR_Scene-mean(FR_Scene)).^2); 
    %}
      %For scene-on and off PSTH
       SpikeFiringRateSceneOn_mean(j,:)=nanmean(SpikeFiringRate(trialindex,:),1);
       SpikeFiringRateSceneOn_sem(j,:)=nanstd(SpikeFiringRate(trialindex,:),[],1)/sqrt(sum(trialindex));

       SpikeFiringRateSceneOff_mean(j,:)=nanmean(SpikeFiringRateOff(trialindex,:),1);
       SpikeFiringRateSceneOff_sem(j,:)=nanstd(SpikeFiringRateOff(trialindex,:),[],1)/sqrt(sum(trialindex));
       
       SpikeFiringRateSceneOnAfter_mean(j,:)=nanmean(SpikeFiringRateAfter(trialindex,:),1);
       SpikeFiringRateSceneOnAfter_sem(j,:)=nanstd(SpikeFiringRateAfter(trialindex,:),[],1)/sqrt(sum(trialindex));
       
       
    
    
  
    end
   
    
    %Define axis for each scene
 %  for jj=1:ceil(length(UniqueSceneNameTotal)/3)
       
  
       for  ii=1:length(UniqueScene)
           plotScene=UniqueScene(ii);
           plotAxis=AxisNum(AxisNum(:,1)==plotScene,2);
           
         %  if ii+(jj-1)*3<=length(UniqueScene)
                %Plot SceneOnset Response
           
           SceneOn_Mean=SpikeFiringRateSceneOn_mean(UniqueScene==plotScene,:);
           SceneOn_SEM=SpikeFiringRateSceneOn_sem(UniqueScene==plotScene,:);
          
           figure(130);
           subplot(ceil(length(UniqueSceneNameTotal)/3),3,plotAxis);
           
            if currBlockType(1)==0
               %Fixed sequence
            shadedErrorBar(TimeSequenceOn, SceneOn_Mean,SceneOn_SEM,'lineprops',{ColorFix(find(UniqueSceneNameTotal==plotScene),:)},'transparent',1);
           
            hold on
            plot([0,0],[min(SceneOn_Mean),max(SceneOn_Mean)],'--k')
  
           else
               %Random sequence
              shadedErrorBar(TimeSequenceOn, SceneOn_Mean,SceneOn_SEM,'lineprops',{ColorRand(find(UniqueSceneNameTotal==plotScene),:)},'transparent',1);  
              
              hold on
              plot([0,0],[min(SceneOn_Mean),max(SceneOn_Mean)],'--k')
               
            end
            if plotAxis==2 
            title_msg=sprintf(' %s, SpikeChannel: %d, TrialNum: %d ',FileName,SpikeChannel,TrialNum);
             title(title_msg);
             
            end

             set(gca,'LineWidth',3)

        set(findall(gca, 'Type', 'Line'),'LineWidth',3);

         set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
         if plotAxis==4
           ylabel('Spikes/s')
         end
         if plotAxis==3*floor(max(AxisNum(:,2))/3)-1
      xlabel('Time from Scene Onset') 
         end
         
         if ~default
    
         xlim([xlim_BeforeSceneOn]);
         end
         
      
      
       %Plot SceneOffset Response
           
           SceneOff_Mean=SpikeFiringRateSceneOff_mean(UniqueScene==plotScene,:);
           SceneOff_SEM=SpikeFiringRateSceneOff_sem(UniqueScene==plotScene,:);
          
           figure(131);
           subplot(ceil(length(UniqueSceneNameTotal)/3),3,plotAxis);
           
            if currBlockType(1)==0
               %Fixed sequence
            shadedErrorBar(TimeSequenceOff, SceneOff_Mean,SceneOff_SEM,'lineprops',{ColorFix(find(UniqueSceneNameTotal==plotScene),:)},'transparent',1);
            
            hold on
            plot([0,0],[min(SceneOff_Mean),max(SceneOff_Mean)],'--k')
  
           else
               %Random sequence
              shadedErrorBar(TimeSequenceOff, SceneOff_Mean,SceneOff_SEM,'lineprops',{ColorRand(find(UniqueSceneNameTotal==plotScene),:)},'transparent',1);  
             
              hold on
              plot([0,0],[min(SceneOff_Mean),max(SceneOff_Mean)],'--k')
               
            end
            if  plotAxis==2 
            title_msg=sprintf(' %s, SpikeChannel: %d, TrialNum: %d ',FileName,SpikeChannel,TrialNum);
             title(title_msg);
             % legend('Fix Seq','Random Seq');
            end

             set(gca,'LineWidth',3)

        set(findall(gca, 'Type', 'Line'),'LineWidth',3);

         set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
         
           if plotAxis==4
           ylabel('Spikes/s')
           end
         if plotAxis==3*floor(max(AxisNum(:,2))/3)-1
           xlabel('Time from Scene Offset') 
         end
         if ~default
            xlim([xlim_AfterSceneOff]);     
         end
     
         
         
      
      
    %Plot After Scene onset Response
           
           SceneOnAfter_Mean=SpikeFiringRateSceneOnAfter_mean(UniqueScene==plotScene,:);
           SceneOnAfter_SEM=SpikeFiringRateSceneOnAfter_sem(UniqueScene==plotScene,:);
          
           figure(132);
           subplot(ceil(length(UniqueSceneNameTotal)/3),3,plotAxis);
           
            if currBlockType(1)==0
               %Fixed sequence
            shadedErrorBar(TimeSequenceAfter, SceneOnAfter_Mean,SceneOnAfter_SEM,'lineprops',{ColorFix(find(UniqueSceneNameTotal==plotScene),:)},'transparent',1);
            
            hold on
            plot([0,0],[min(SceneOnAfter_Mean),max(SceneOnAfter_Mean)],'--k')
  
           else
               %Random sequence
              shadedErrorBar(TimeSequenceAfter, SceneOnAfter_Mean,SceneOnAfter_SEM,'lineprops',{ColorRand(find(UniqueSceneNameTotal==plotScene),:)},'transparent',1);  
             
              hold on
              plot([0,0],[min(SceneOnAfter_Mean),max(SceneOnAfter_Mean)],'--k')
               
            end
            if  plotAxis==2 
            title_msg=sprintf(' %s, SpikeChannel: %d, TrialNum: %d ',FileName,SpikeChannel,TrialNum);
             title(title_msg);
             % legend('Fix Seq','Random Seq');
            end

             set(gca,'LineWidth',3)

        set(findall(gca, 'Type', 'Line'),'LineWidth',3);

         set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
         
          if plotAxis==4
           ylabel('Spikes/s')
           end
         if plotAxis==3*floor(max(AxisNum(:,2))/3)-1
            xlabel('Time from Scene Onset') ;  
         end
         
         if ~default
     
          xlim([xlim_AfterSceneOn]);
         end
           
     
      
      
      
      
      
      
      
            
      %     end
       end
    
%   end
    
    
   
   
   
    
end

   figure(130);
  l1=legend();
  legendString={'Fixed Seq','','Random Seq',''};
  uniqueSeqType=(UniqueBlockType(:,1));
  
  
  
 if length(uniqueSeqType)==2
           l1.String=legendString;
 else
    if uniqueSeqType==0
            l1.String={'Fixed Seq',''};
    else
        l1.String={'Random Seq',''};
    end
 end

 
 
 



end