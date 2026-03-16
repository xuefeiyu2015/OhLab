function ForagingTask_LocationSelectivity(Data, StartTrial, EndTrial,ShowFigureFlag,OutputFlag);
%This function is to analysis the location tuning measured inside the
%foraging task;
%Xuefei Yu, 11,13,2019
%close all

%%%%%%%%%%%%%%%%%User controller starts here

%Parameter to define the sliding window to measure firing rates
StepSize=20;%in ms
BinWidth=20;%in ms

%Parameter seting for interval selecting plot period
StartOffset=200;%in ms
StopOffset=50;%in ms
IntervalAfterMarkerStart=500;%in ms;%This for align marker

%Average for statistics region;
AverageRegion=[0,200];

%For Raster plot 
%Barlength=0.5;
%LineWidth=1;
%Color Definition for each location
Color={'r',[0.3,0.8,0.0],[0.3,0.5,1],[0.7,0.3,1],[0.4,0.4,0],[0.28,0.78,0.89],'k',[0.5,0.2,1],[0.6,0.1,0.4],[0.3,0.5,1],[0.6,0.5,0.1],[0,0.8,0.4]};


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Unpack all the data base and task settings
ECode;
%TaskType=Data.TaskType;
TrialNum=EndTrial-StartTrial+1;
Selection=StartTrial:EndTrial;

%FixationPointIndex=Data.FixPointIndex(Selection,:);
%FixLoc=Data.FixLoc(Selection,:,:);



TargetIndex=Data.TargetIndex(Selection,2:end,1);
TargetName=Data.TargetIndex(Selection,2:end,2);

TargetLoc=Data.TargetLoc(Selection,:,:);

TargetNum=Data.TargetNum(Selection,:);
%TargetSize=Data.TargetSize(Selection,:);
TargetPolarLoc=Data.TargetPolarLoc(Selection,:,:);
%ChoiceOutcome=Data.ChoiceOutcome(Selection,:);


%GoodObjectLib=Data.GoodObjectLib;
%BadObjectLib=Data.BadObjectLib;
%SceneLib=Data.SceneLib;

GoodObject=mod(TargetIndex,2)==1;
%GoodObjectID=TargetIndex(:,2:end,2);
BadObject=mod(TargetIndex,2)==0;
%BadObject=arrayfun(@(x)( ismember(x,BadObjectLib(:,3))),TargetIndex(:,2:end,2));%Use real name
%BadObjectID=TargetIndex(:,2:end,2);
%GoodBadIndex=GoodObject*1+BadObject*2;%1:Good Object;2:Bad Object

%SceneIndex=Data.SceneIndex(Selection,:);

%EyeChannel_X=Data.EyeChannel_X_L(Selection,:); %Only import the first two channels
%EyeChannel_Y=Data.EyeChannel_Y_L(Selection,:);
%EyeBinWidth=Data.EyeBinWidth;
%EyeTime=Data.EyeTime(Selection,:);

EventChannel=Data.EventChannel(Selection,:);
EventTimeChannel=Data.EventTimeChannel(Selection,:);
%EventBin=Data.EventBin(Selection,:);

SpikeTime=squeeze(Data.SpikeTimeData(:,:,Selection))';
%SpikeBin=squeeze(Data.SpikeBinData(:,:,Selection))';



%StartOffsetBin=StartOffset/EyeBinWidth;
%StopOffsetBin=StopOffset/EyeBinWidth;
%IntervalInBin=IntervalAfterMarkerStart/EyeBinWidth;


 
%Another way of align: Align with the target onset
AlignTimeWithTargetOnset=FindOutTime(EventChannel,EventTimeChannel,TGTCD);
%TargetOffRelativeToOnset=FindOutTime(EventChannel(:,5:end),EventTimeChannel(:,5:end),TGTOFF)-AlignTimeWithTargetOnset;%Start searching after visual onset
TargetOffRelativeToOnset=IntervalAfterMarkerStart+0*AlignTimeWithTargetOnset;%Start searching after visual onset
TargetOffRelativeTime=TargetOffRelativeToOnset+StopOffset;









%################################################
%Calculate firing rate aligned on target onset
MaxTgtTime=0;
MinTgtTime=1000;
GoodFR_AlgnTgt=[];
GoodTimeSeq_AlgnTgt=[];
TargetAngleGood=[];

BadFR_AlgnTgt=[];
BadTimeSeq_AlgnTgt=[];
TargetAngleBad=[];


for i=1:max(TargetNum(:,2))
    AlignTimeForTargetN=AlignTimeWithTargetOnset(:,i);
    SpikeTimeAlignTgtN=SpikeTime-repmat(AlignTimeForTargetN,1,size(SpikeTime,2));
    SelectIntervalTimeTgt=[repmat(-StartOffset,size(AlignTimeForTargetN,1),1),TargetOffRelativeTime(:,i)];
    
    
    [FRAlignTgtN,TimeSeqTgtN]=CalculateSpikeFiringRate(SpikeTimeAlignTgtN,StepSize,BinWidth,SelectIntervalTimeTgt);
   % SpikeTimeTgtNInInterval{i}=SelectSpikeTimeInterval(SpikeTimeAlignTgtN,SelectIntervalTimeTgt(:,1),SelectIntervalTimeTgt(:,2));%For Rasterplot
    
    MaxTgtTime=nanmax(MaxTgtTime, nanmax(SelectIntervalTimeTgt(:,2)));
    MinTgtTime=nanmin(MinTgtTime, nanmin(SelectIntervalTimeTgt(:,1)));
    
    TargetAngle=round(TargetPolarLoc( :,i+1,1));
    GoodObjectIndex=GoodObject(:,i);
    BadObjectIndex=BadObject(:,i);
    
    TargetAngleGood=[TargetAngleGood;TargetAngle(GoodObjectIndex)];
    GoodFR_AlgnTgt=[GoodFR_AlgnTgt;FRAlignTgtN(GoodObjectIndex,:)];
    GoodTimeSeq_AlgnTgt=[GoodTimeSeq_AlgnTgt;TimeSeqTgtN(GoodObjectIndex,:)];%At present, they are the same
    
    TargetAngleBad=[TargetAngleBad;TargetAngle(BadObjectIndex)];
    BadFR_AlgnTgt=[BadFR_AlgnTgt;FRAlignTgtN(BadObjectIndex,:)];
    BadTimeSeq_AlgnTgt=[BadTimeSeq_AlgnTgt;TimeSeqTgtN(BadObjectIndex,:)];%At present, they are the same
    
    %currentLine=GoodBadIndex(:,i);
    %GoodBadIndex;
    
end

%Separate trials according to the target type and location 
%First Good Object
UniqueLoc_Good=unique(TargetAngleGood);
GoodTimeForPlot=nanmean(GoodTimeSeq_AlgnTgt);
BadTimeForPlot=nanmean(BadTimeSeq_AlgnTgt);
TotalGoodNum=0;
GoodFR_AlgnTgt_InRegion_Trial=[];
GroupLabelGood=[];
for i=1:length(UniqueLoc_Good)
    Select_Trial=TargetAngleGood==UniqueLoc_Good(i);
    GoodFR_AlgnTgt_select_mean(i,:)=nanmean(GoodFR_AlgnTgt(Select_Trial,:),1);
    GoodFR_AlgnTgt_select_sem(i,:)=nanstd(GoodFR_AlgnTgt(Select_Trial,:),[],1)/sqrt(sum(Select_Trial));
    Selection_For_Trial=nanmean(GoodFR_AlgnTgt(Select_Trial,GoodTimeForPlot<AverageRegion(2)&GoodTimeForPlot>AverageRegion(1) )');
    GoodFR_AlgnTgt_InRegion_Mean(i)=nanmean(Selection_For_Trial); 
    GoodFR_AlgnTgt_InRegion_SEM(i)=nanstd(Selection_For_Trial)/sqrt(sum(Select_Trial)); 
    GoodFR_AlgnTgt_InRegion_SSE(i)=sum((Selection_For_Trial-mean(Selection_For_Trial)).^2); 
    TotalGoodNum=TotalGoodNum+sum(Select_Trial);
    GoodFR_AlgnTgt_InRegion_Trial=[GoodFR_AlgnTgt_InRegion_Trial,Selection_For_Trial];
    GroupLabelGood=[GroupLabelGood,i*ones(1,sum(Select_Trial))];
    
end
%Then Bad Object

UniqueLoc_Bad=unique(TargetAngleBad);
TotalBadNum=0;
BadFR_AlgnTgt_InRegion_Trial=[];
GroupLabelBad=[];
for i=1:length(UniqueLoc_Bad)
     Select_Trial=TargetAngleBad==UniqueLoc_Bad(i);
    
    BadFR_AlgnTgt_select_mean(i,:)=nanmean(BadFR_AlgnTgt(Select_Trial,:),1);
    BadFR_AlgnTgt_select_sem(i,:)=nanstd(BadFR_AlgnTgt(Select_Trial,:),[],1)/sqrt(sum(Select_Trial));
    Selection_For_Trial=nanmean(BadFR_AlgnTgt(Select_Trial,BadTimeForPlot<AverageRegion(2)&BadTimeForPlot>AverageRegion(1) )');
    BadFR_AlgnTgt_InRegion_Mean(i)=nanmean(Selection_For_Trial); 
    BadFR_AlgnTgt_InRegion_SEM(i)=nanstd(Selection_For_Trial)/sqrt(sum(Select_Trial)); 
    BadFR_AlgnTgt_InRegion_SSE(i)=sum((Selection_For_Trial-mean(Selection_For_Trial)).^2); 
    TotalBadNum=TotalBadNum+sum(Select_Trial);
    BadFR_AlgnTgt_InRegion_Trial=[BadFR_AlgnTgt_InRegion_Trial,Selection_For_Trial];
    GroupLabelBad=[GroupLabelBad,i*ones(1,sum(Select_Trial))];
     
end
%Comparison between Good and Bad object
GoodFR_AlgnTgt_mean=nanmean(GoodFR_AlgnTgt);
BadFR_AlgnTgt_mean=nanmean(BadFR_AlgnTgt);
GoodFR_AlgnTgt_sem=nanstd(GoodFR_AlgnTgt)/sqrt(size(GoodFR_AlgnTgt,1));
BadFR_AlgnTgt_sem=nanstd(BadFR_AlgnTgt)/sqrt(size(BadFR_AlgnTgt,1));



%Select region to compute the selectivity
GoodFiringRateInRegion=nanmean(GoodFR_AlgnTgt(:,GoodTimeForPlot<AverageRegion(2)&GoodTimeForPlot>AverageRegion(1) )');
BadFiringRateInRegion=nanmean(BadFR_AlgnTgt(:,BadTimeForPlot<AverageRegion(2)&BadTimeForPlot>AverageRegion(1) )');
ValueSelectivity=(mean(GoodFiringRateInRegion)-mean(BadFiringRateInRegion))/sqrt((std(GoodFiringRateInRegion)+std(BadFiringRateInRegion))/2);%(mu(High)-mu(Low))/sqrt((sigma(High)+sigma(Low))/2)
[~,p_ValueSelectivity]=ttest2(GoodFiringRateInRegion',BadFiringRateInRegion');
%Direction selectivity: (Rmax-Rmin)/(Rmax-Rmin+2*sqrt(SSE/(N-M))); N:Total
%condition ;M: direction condition
LocationSelectivityGood=(max(GoodFR_AlgnTgt_InRegion_Mean)-min(GoodFR_AlgnTgt_InRegion_Mean))/(max(GoodFR_AlgnTgt_InRegion_Mean)-min(GoodFR_AlgnTgt_InRegion_Mean)+2*sqrt(sum(GoodFR_AlgnTgt_InRegion_SSE)/(TotalGoodNum-length(UniqueLoc_Good))));
p_LocationSelectivityGood= anovan(GoodFR_AlgnTgt_InRegion_Trial',GroupLabelGood','display','off');
LocationSelectivityBad=(max(BadFR_AlgnTgt_InRegion_Mean)-min(BadFR_AlgnTgt_InRegion_Mean))/(max(BadFR_AlgnTgt_InRegion_Mean)-min(BadFR_AlgnTgt_InRegion_Mean)+2*sqrt(sum(BadFR_AlgnTgt_InRegion_SSE)/(TotalBadNum-length(UniqueLoc_Bad))));
p_LocationSelectivityBad= anovan(BadFR_AlgnTgt_InRegion_Trial',GroupLabelBad','display','off');
%TrialNum
%Start to plot
PlotOffset=50;
XRange=[-StartOffset-PlotOffset,IntervalAfterMarkerStart+StopOffset+PlotOffset];

%PSTH: Good Object v.s. Bad Object
if ishandle(102) 
    close 102; end
figure(102);set(102,'color','w');
set(102,'Position', [100,100 1300,500], 'Name', 'PSTH differed by object');
axes('position',[0.1,0.7,0.3,0.2])
area(AverageRegion,[max([GoodFR_AlgnTgt_mean+GoodFR_AlgnTgt_sem,BadFR_AlgnTgt_mean+BadFR_AlgnTgt_sem])+5,max([GoodFR_AlgnTgt_mean+GoodFR_AlgnTgt_sem,BadFR_AlgnTgt_mean+BadFR_AlgnTgt_sem])+5],...
    'FaceColor','k','FaceAlpha',.1,'EdgeAlpha',.1,'HandleVisibility','off');


shadedErrorBar(GoodTimeForPlot,GoodFR_AlgnTgt_mean,GoodFR_AlgnTgt_sem,'lineprops','r','transparent',1);
hold on
shadedErrorBar(BadTimeForPlot,BadFR_AlgnTgt_mean,BadFR_AlgnTgt_sem,'lineprops','b','transparent',1);
xlim(XRange);
ylim([0,max([GoodFR_AlgnTgt_mean+GoodFR_AlgnTgt_sem,BadFR_AlgnTgt_mean+BadFR_AlgnTgt_sem]+10)]);
xlabel('Time from target onset(ms)')
ylabel('Firing Rate (Hz)')

plot([0,0],[0,max([GoodFR_AlgnTgt_mean+GoodFR_AlgnTgt_sem,BadFR_AlgnTgt_mean+BadFR_AlgnTgt_sem]+5)],'-.k','LineWidth',1.5);

set(gca,'LineWidth',2)
set(gca,'FontWeight','Bold')
ldg=legend();
ldg.Position=[0.4    0.8    0.0854    0.0880];
ldg.Box='off';
ldg.String={'Good Object','Bad Object'};
str_selectivity={sprintf('Value Selectivity: %2.3f',ValueSelectivity);
                 sprintf('p_V_S: %2.3f',p_ValueSelectivity)};
anna=annotation('textbox',[.408 .5 .3 .3],'String',str_selectivity,'FitBoxToText','on');
if p_ValueSelectivity<0.05
    anna.Color='r';
end
%text(0.8,0.8,str_selectivity);



title('Comparison between Good and Bad Objects');
box off;

%PSTH: Good object in different location
axes('position',[0.1,0.4,0.3,0.2])
%Add a shaded area to indicate analysis area
area(AverageRegion,[max(max(max(GoodFR_AlgnTgt_select_mean+GoodFR_AlgnTgt_select_sem)),max(max(BadFR_AlgnTgt_select_mean+BadFR_AlgnTgt_select_sem)))+5,max(max(max(GoodFR_AlgnTgt_select_mean+GoodFR_AlgnTgt_select_sem)),max(max(BadFR_AlgnTgt_select_mean+BadFR_AlgnTgt_select_sem)))+5],...
    'FaceColor','k','FaceAlpha',.1,'EdgeAlpha',.1,'HandleVisibility','off');
for i=1:length(UniqueLoc_Good)
    shadedErrorBar(GoodTimeForPlot,GoodFR_AlgnTgt_select_mean(i,:),GoodFR_AlgnTgt_select_sem(i,:),'lineprops',{Color{i}},'transparent',1);
    hold on
end
plot([0,0],[0,max(max(max(GoodFR_AlgnTgt_select_mean+GoodFR_AlgnTgt_select_sem)),max(max(BadFR_AlgnTgt_select_mean+BadFR_AlgnTgt_select_sem)))+5],'-.k','LineWidth',1.5);
ylim([0,max(max(max(GoodFR_AlgnTgt_select_mean+GoodFR_AlgnTgt_select_sem)),max(max(BadFR_AlgnTgt_select_mean+BadFR_AlgnTgt_select_sem)))+10]);
xlim(XRange);
box off;
ldg=legend();
xlabel('Time from target onset(ms)')
ylabel('Firing Rate (Hz)')  
set(gca,'LineWidth',2)
set(gca,'FontWeight','Bold')
ldg.Position=[0.4     0.57    0.0854    0.0880];
ldg.Box='off';

ldg.String=string(UniqueLoc_Good)';
ldg.LineWidth=15;
tit=title('Good Object in Different Location');
tit.Color='r';
str_selectivity={sprintf('Location Selectivity(Good): %2.3f',LocationSelectivityGood);
                 sprintf('p_L_S_G_O_O_D: %2.5f',p_LocationSelectivityGood)};
anna=annotation('textbox',[.408 .32 .2 .2],'String',str_selectivity,'FitBoxToText','on');
if p_LocationSelectivityGood<0.05
    anna.Color='r';
end


%PSTH: Bad object in different location
axes('position',[0.1,0.1,0.3,0.2])
%Add a shaded area to indicate analysis area
area(AverageRegion,[max(max(max(GoodFR_AlgnTgt_select_mean+GoodFR_AlgnTgt_select_sem)),max(max(BadFR_AlgnTgt_select_mean+BadFR_AlgnTgt_select_sem)))+5,max(max(max(GoodFR_AlgnTgt_select_mean+GoodFR_AlgnTgt_select_sem)),max(max(BadFR_AlgnTgt_select_mean+BadFR_AlgnTgt_select_sem)))+5],...
    'FaceColor','k','FaceAlpha',.1,'EdgeAlpha',.1,'HandleVisibility','off');
for i=1:length(UniqueLoc_Bad)
    shadedErrorBar(BadTimeForPlot,BadFR_AlgnTgt_select_mean(i,:),BadFR_AlgnTgt_select_sem(i,:),'lineprops',{Color{i}},'transparent',1);
    hold on
end
plot([0,0],[0,max(max(max(GoodFR_AlgnTgt_select_mean+GoodFR_AlgnTgt_select_sem)),max(max(BadFR_AlgnTgt_select_mean+BadFR_AlgnTgt_select_sem)))+5],'-.k','LineWidth',1.5);

ldg=legend();
xlim(XRange);
ylim([0,max(max(max(GoodFR_AlgnTgt_select_mean+GoodFR_AlgnTgt_select_sem)),max(max(BadFR_AlgnTgt_select_mean+BadFR_AlgnTgt_select_sem)))+10]);



xlabel('Time from target onset(ms)')
ylabel('Firing Rate (Hz)')   
set(gca,'LineWidth',2)
set(gca,'FontWeight','Bold')
ldg.Position=[0.4     0.25    0.0854    0.0880];
ldg.Box='off';
ldg.String=string(UniqueLoc_Good);
tit=title('Bad Object in Different Location');
tit.Color='b';
str_selectivity={sprintf('Location Selectivity(Bad): %2.3f',LocationSelectivityBad);
                 sprintf('p_L_S_B_A_D: %2.5f',p_LocationSelectivityBad)};
anna=annotation('textbox',[.408 0 0.2 0.2],'String',str_selectivity,'FitBoxToText','on');
if p_LocationSelectivityBad<0.05
    anna.Color='r';
end
box off


%PSTH: Location Tuning Polar Plot
axes('position',[0.45,0.2,0.6,0.6])
LineWidth=3;
UniqueLoc_GoodPolar=UniqueLoc_Good+90;
if sum(UniqueLoc_GoodPolar>=360)
    UniqueLoc_GoodPolar(UniqueLoc_GoodPolar>=360)=UniqueLoc_GoodPolar(UniqueLoc_GoodPolar>=360)-360;%Make zero pointing towards upward
    
end

polarwitherrorbar([UniqueLoc_GoodPolar'/180*pi, UniqueLoc_GoodPolar(1)/180*pi], [GoodFR_AlgnTgt_InRegion_Mean, GoodFR_AlgnTgt_InRegion_Mean(1)],[GoodFR_AlgnTgt_InRegion_SEM,GoodFR_AlgnTgt_InRegion_SEM(1)],...
    '-r',LineWidth);
hold on
UniqueLoc_BadPolar=UniqueLoc_Bad+90;
if sum(UniqueLoc_BadPolar>=360)
    UniqueLoc_BadPolar(UniqueLoc_BadPolar>=360)=UniqueLoc_BadPolar(UniqueLoc_BadPolar>=360)-360;%Make zero pointing towards upward
    
end

polarwitherrorbar([UniqueLoc_BadPolar'/180*pi, UniqueLoc_BadPolar(1)/180*pi], [BadFR_AlgnTgt_InRegion_Mean, BadFR_AlgnTgt_InRegion_Mean(1)],[BadFR_AlgnTgt_InRegion_SEM,BadFR_AlgnTgt_InRegion_SEM(1)],...
    '-b',LineWidth);
title('Polar Plot For Location Selectivity')
ShowLabel=[0:30:360];
thetaticklabels([ShowLabel(end-3:end-1),ShowLabel(1:end-4)]);
set(gca,'fontsize',15);

%set(gca,'XTickLabel',string([UniqueLoc_Bad'/180*pi, UniqueLoc_Bad(1)/180*pi]));
%set(findall(gcf,'fontsize',10),'fontsize',15);




















