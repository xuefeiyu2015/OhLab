%Batch analysis for the learning of foraging explore task
%Xuefei Yu; 03242020
Path='/Users/yux8/WorkRelated_YXF/Program_Matlab_Local/DataAnalysis/Results/ForagingExploreTask';
cd(Path);
clear all

FilesName =cellstr( uigetfile('*.mat', 'Multiselect', 'on'));
for i=1:length(FilesName)
  
     load(FilesName{i});
     FileName=FilesName{i};
     ProportionOfAdvanceTrials(i)=OutputData.ProportionOfAdvanceTrials;
     FirstSearchingTime_Mean=OutputData.FirstSearchingTime_Mean/1000;%Changed from ms to s from 05172020
     
     FirstSearchingTime_SEM=OutputData.FirstSearchingTime_SEM/1000;
     BlocType=OutputData.BlocType;
     TrialNumInEachBlock=OutputData.TrialNumInEachBlock;
     
     
     
    FirstSearchingTime_Mean_All(i)=sum(FirstSearchingTime_Mean(BlocType==10|BlocType==11).*TrialNumInEachBlock(BlocType==10|BlocType==11))/sum(TrialNumInEachBlock(BlocType==10|BlocType==11));
    %Rough
    FirstSearchingTime_SEM_All(i)=sum(FirstSearchingTime_SEM(BlocType==10|BlocType==11).*TrialNumInEachBlock(BlocType==10|BlocType==11))/sum(TrialNumInEachBlock(BlocType==10|BlocType==11));
    
     ExploreSearchingTimeMean=sum(FirstSearchingTime_Mean(BlocType==0|BlocType==1).*TrialNumInEachBlock(BlocType==0|BlocType==1))/sum(TrialNumInEachBlock(BlocType==0|BlocType==1));
    ExploreSearchingTime_SEM=sum(FirstSearchingTime_SEM(BlocType==0|BlocType==1).*TrialNumInEachBlock(BlocType==0|BlocType==1))/sum(TrialNumInEachBlock(BlocType==0|BlocType==1));
    
end

if ishandle(100) 
    close(100); end
figure(100);set(100,'color','w');
set(100,'Position', [100,100 800,500], 'Name', 'Mean Exploration Time Training');
Plot_Y=[FirstSearchingTime_Mean_All];
Plor_Err=[FirstSearchingTime_SEM_All];
Plot_X=1:length(Plot_Y);
errorbar(Plot_Y,Plor_Err,'-s','Color',[0.83 0.14 0.14],'LineWidth',5);
hold on 
errorbar(length(Plot_Y)+1,ExploreSearchingTimeMean,ExploreSearchingTime_SEM,'-s','Color','b','LineWidth',5);
xlim([0.5,length(Plot_Y)+1.5]);
box off
for i=1:length(Plot_Y)
    strDay{i}=sprintf('Day %d',i);
  
end

strDay{i+1}='Rand Loc';
xticklabels(strDay);
ylabel('Average Exploration Time (s)');
title('Average Exploration Time during Training');
set(gca,'FontWeight','Bold');
set(gca,'FontSize',25);
set(gca,'LineWidth',3)


if ishandle(101) 
    close(101); end
figure(101);set(101,'color','w');
set(101,'Position', [100,100 800,500], 'Name', 'Porportion of Trigger TargetShowup Rate');

plot(ProportionOfAdvanceTrials,'-s','Color',[0.83 0.14 0.14],'LineWidth',5,'MarkerSize',20,'MarkerFaceColor',[0.83 0.14 0.14]);
xticklabels({strDay{1:end-1}});
xticks(Plot_X);
xlim([0.5,length(ProportionOfAdvanceTrials)+0.5]);
box off
ylabel('Object Locating Success Rate.');
title('Suceess rate during training');
set(gca,'FontWeight','Bold');
set(gca,'FontSize',25);
set(gca,'LineWidth',3)
ylim([0.4,1]);


keyboard