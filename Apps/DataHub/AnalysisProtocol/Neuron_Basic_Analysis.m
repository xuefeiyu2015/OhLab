function Neuron_Basic_Analysis(Data, StartTrial, EndTrial,ShowFigureFlag,OutputFlag);
%Analysis Neuronal Reponse Stability;
%YXF 11110219
close all
SpikeTime_ori=squeeze(Data.SpikeTimeData(1,:,StartTrial:EndTrial))';
%EyeTime=Data.EyeTime;





ISI=[];
for i=1:size(SpikeTime_ori,1)
    Duration=(SpikeTime_ori(i,end)-SpikeTime_ori(i,end-1));
    SpikeTime=SpikeTime_ori(i,1:end-2);
    SpikeCount_Trial=length(SpikeTime(~isnan(SpikeTime)));
    FiringRate_Trial(i)=SpikeCount_Trial/(Duration/1000);
    ISI=[ISI diff(SpikeTime(~isnan(SpikeTime)))];
    
end

ISID=diff(SpikeTime_ori(:,1:end-2),1,2);
FiringRate_Trial=1000./nanmean(ISID,2);
ISI=reshape(ISID,numel(ISID),1);
ISI=ISI(~isnan(ISI));




figure
set(gcf,'color','w');
plot([1:size(SpikeTime_ori,1)],FiringRate_Trial,'-or');
xlabel('Trial Number')
ylabel('Average Firing Rate(Hz)')
title('Firing Rate Changes With Trials')
figure
set(gcf,'color','w')
%center=[1:2:max(ISI)];
BinWidth=1;
SepPoint=100;
Edge=[0:BinWidth:SepPoint,max(ISI)+1];
%[N,BIN] =histc(ISI,center);
[Num,EdgeO]=histcounts(ISI,Edge);
x=[BinWidth/2:BinWidth:SepPoint,SepPoint];

bar(x,Num,'LineWidth',3,'FaceColor','b','EdgeColor','b');

%bar(center,N,'histc','facecolor','o') ;
xlabel('Inter Spike Interval(ms)');
ylabel('Number of spikes');
title('Distribution of InterSpikeInterval');
%Number of spikes with ISI<1ms
SpikeWithinRefactoryPeriod=sum(ISI<1);
Proprotion=SpikeWithinRefactoryPeriod/length(ISI)*100;%In %
Shown=sprintf('P(ISI<1ms) : %2.2f  %%',Proprotion);
xlim([-BinWidth/2,x(end)+BinWidth]);
xloc=xlim;
xloc=xloc(2)/3;
yloc=ylim;
yloc=yloc(2)*3/4;
text(xloc,yloc,Shown,'FontWeight','Bold','FontSize',25);
set(gca,'LineWidth',3,'FontWeight','Bold','FontSize',20);
%set(gca,'xscale','log')


end