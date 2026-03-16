function Batch_FreeViewBehavior(Data,  StartFile, EndFile,ShowFigureFlag,OutputFlag);
BatchFileName=Data.BatchFileName;
%Batch data file path
FilesName=Data.ResultFilePath(StartFile: EndFile);
RecordDate=Data.RecordDate(StartFile: EndFile);
FileNum=EndFile-StartFile+1;
FilesName_Each=Data.FileName;

UseAll = 1;%1: Use all the recorded sites; 0: only use significant sites 

%{
%Color Definition for each events
Color=[lines(6);prism(6);pink(6)];
%Color(2,:)=copper(1)+0.5;
Color=[Color;hsv(10)];
Color(18,:)=[0.1,0.457,0.194];
%}
%Load file
for i=1:length(FilesName)
  clear OutputData;
  
     load(FilesName{i}{1});
     
     OutputDataTemp=OutputData;
     
     
     Task(i)=OutputData.FreeViewBehavior(1).TaskCode;
     TaskName{i}=OutputData.FreeViewBehavior(1).Task;
     
     
     Data =OutputData.FreeViewBehavior.DataStamp;
    


      Monkey = FilesName_Each{i}(1);

          if Monkey == 'R'
              MonkeyIndex(i) = 1;
          else
              MonkeyIndex(i) = 2;
          end


     PropClSac(i) = Data('PropClSac');
     p_prop(i) = Data('p_prop');

     NumberOfAngle_Control(i,:) = Data('NumberOfAngle_Control');
     NumberOfAngle_Stim(i,:) = Data('NumberOfAngle_Stim');

     SaccadeCountStim_Sum(i,:) = Data('SaccadeCountStim');
     SaccadeCountControl_Sum(i,:) = Data('SaccadeCountControl');

  %   theta_stim{i}=Data('theta_Stim');
   %  amp_stim{i}=Data('Amp_Stim');
   %  theta_control{i}=Data('theta_Control');
    % amp_control{i}=Data('Amp_Control');


      Theta_Stim{i} = Data('theta_Stim');
      Theta_Control{i} = Data('theta_Control');



      NSac_Total_Stim(i) = Data('NSac_Total_Stim');
      NSac_Total_Control(i) = Data('NSac_Total_Control');
      SacFreq_Total_Stim(i) = Data('SacFreq_Total_Stim');
      SacFreq_Total_Control(i) = Data('SacFreq_Total_Control');


      SaccadeFreqStim_Sum(i,:) = Data('SaccadeFreqStim_Sum');
      SaccadeFreqControl_Sum(i,:)=Data('SaccadeFreqControl_Sum');




     
     
     TimeSeq = Data('TimeForSaccadeRaster');

     Data =OutputData.FreeViewNeural.DataStamp;


     SM_p(i) = Data('StimModulation_p');

     %%
 
end

    

%Significant case 
Total = SM_p < 0.05 |p_prop < 0.05;
Sig = p_prop < 0.05;
NumSig = sum(Sig);
NumTotal = sum(Total);
disp(sprintf('Total: %d',NumTotal));
disp(sprintf('SigNum: %d',NumSig));





if UseAll
    PropCS_total = PropClSac;
    PropCS_Sig = PropClSac;
    Sig = p_prop<1;
    Angle_Stim = Theta_Stim;
    Angle_Control = Theta_Control;
    SaccadeFreqStim_Sum=SaccadeFreqStim_Sum;
    SaccadeFreqControl_Sum=SaccadeFreqControl_Sum;





    
else
    PropCS_total = PropClSac(Total);
    PropCS_Sig = PropClSac(Sig);

    Angle_Stim = Theta_Stim(Sig);
    Angle_Control = Theta_Control(Sig);


     NSac_Total_Stim = NSac_Total_Stim(Sig);
      NSac_Total_Control = NSac_Total_Control(Sig);
      SacFreq_Total_Stim = SacFreq_Total_Stim(Sig);
      SacFreq_Total_Control = SacFreq_Total_Control(Sig);

      SaccadeFreqStim_Sum = SaccadeFreqStim_Sum(Sig,:);
      SaccadeFreqControl_Sum = SaccadeFreqControl_Sum(Sig,:);
      MonkeyIndex=MonkeyIndex(Sig);

      NumberOfAngle_Control=NumberOfAngle_Control(Sig,:);
      NumberOfAngle_Stim=NumberOfAngle_Stim(Sig,:);

      SaccadeCountStim_Sum=SaccadeCountStim_Sum(Sig,:);
      SaccadeCountControl_Sum=SaccadeCountControl_Sum(Sig,:);


  %  theta_stim=cell2mat(theta_stim(Sig))+90;
  %  amp_stim=cell2mat(amp_stim(Sig));

  %  MonkeyIndex_Sig = MonkeyIndex(Sig);

end

TotalNumSaccade_Stim = sum(NSac_Total_Stim);
TotalNumSaccade_Control = sum(NSac_Total_Control);

SacFreq_Total_Stim_Mean = mean(SacFreq_Total_Stim);
SacFreq_Total_Control_Mean = mean(SacFreq_Total_Control);

[p,h] = ttest(SacFreq_Total_Stim,SacFreq_Total_Control);

SaccadeFreqStim_Sum_R=SaccadeFreqStim_Sum(MonkeyIndex == 1,:);
SaccadeFreqStim_Sum_A=SaccadeFreqStim_Sum(MonkeyIndex == 2,:);

SaccadeFreqControl_Sum_R=SaccadeFreqControl_Sum(MonkeyIndex == 1,:);
SaccadeFreqControl_Sum_A=SaccadeFreqControl_Sum(MonkeyIndex == 2,:);

SaccadeFreqStim_Sum_R_Mean = mean(SaccadeFreqStim_Sum_R,1);
SaccadeFreqStim_Sum_A_Mean = mean(SaccadeFreqStim_Sum_A,1);

SaccadeFreqControl_Sum_R_Mean = mean(SaccadeFreqControl_Sum_R,1);
SaccadeFreqControl_Sum_A_Mean = mean(SaccadeFreqControl_Sum_A,1);

%keyboard
%{
 DataPath = '/Users/yux8/WorkRelated_YXF/Presentation_Summary/Optogenetics/Optogenetics_paper_preparation/NatureNeuroscienceSubmission/DataShare/ExportData';
    
    fileName1 = fullfile(DataPath, 'FEF_SaccadeFreq_OpticalStim_MonkeyR.mat');
    tmp1 =SaccadeFreqStim_Sum_R_Mean;
    save(fileName1, 'tmp1');
    fileName2 = fullfile(DataPath, 'FEF_SaccadeFreq_Control_MonkeyR.mat');
    tmp2 =SaccadeFreqControl_Sum_R_Mean;
    save(fileName2, 'tmp2');

    fileName3 = fullfile(DataPath, 'FEF_SaccadeFreq_OpticalStim_MonkeyA.mat');
    tmp3 =SaccadeFreqStim_Sum_A_Mean;
    save(fileName3, 'tmp3');
    fileName4 = fullfile(DataPath, 'FEF_SaccadeFreq_Control_MonkeyA.mat');
    tmp4 =SaccadeFreqControl_Sum_A_Mean;
    save(fileName4, 'tmp4');
%}

%keyboard
%{
figure
plot(SaccadeFreqStim_Sum_R_Mean)
hold on
plot(SaccadeFreqStim_Sum_A_Mean)
keyboard
%}
%keyboard

%{
[x, y] = pol2cart(theta_stim, amp_stim/max(amp_stim));
%{
figure
plot(x,y,'or');

figure
polarplot(theta_stim,amp_stim,'.');
%}

% Define number of bins
numBins = 100;

% Create bin edges based on the data ranges
xEdges = linspace(min(x), max(x), numBins + 1);  % +1 because edges are boundaries
yEdges = linspace(min(y), max(y), numBins + 1);

% Create 2D histogram
counts = histcounts2(x, y, xEdges, yEdges);
% Smooth the histogram counts using Gaussian filtering
countsSmoothed = imgaussfilt(counts, 2);  % Adjust the 2 for different levels of smoothing
% Get the bin centers
binCentersX = (xEdges(1:end-1) + xEdges(2:end)) / 2;
binCentersY = (yEdges(1:end-1) + yEdges(2:end)) / 2;

% Create a meshgrid of bin centers for plotting
[XGrid, YGrid] = meshgrid(binCentersX, binCentersY);


contourf(XGrid, YGrid,countsSmoothed)
%polarhistogram(theta_stim,50);
%}
Angle_Stim = cell2mat(Angle_Stim);
Angle_Control = cell2mat(Angle_Control);

[h,p_distributionDifference]=kstest2(Angle_Stim,Angle_Control);


BinNum = 5;
edges = 0:0.1:1;
[num_total,edge] = histcounts(PropCS_total,edges);
[num_sig,edge] = histcounts(PropCS_Sig,edges);

%{
NumAngleControl = sum(NumberOfAngle_Control(Sig,:),1);

NumAngleStim = sum(NumberOfAngle_Stim(Sig,:),1);
%}
NumAngleControl = sum(NumberOfAngle_Control,1);

NumAngleStim = sum(NumberOfAngle_Stim,1);

NumAngleControl_R = sum(NumberOfAngle_Control(MonkeyIndex==1,:),1);

NumAngleStim_R = sum(NumberOfAngle_Stim(MonkeyIndex==1,:),1);

NumAngleControl_A = sum(NumberOfAngle_Control(MonkeyIndex==2,:),1);

NumAngleStim_A = sum(NumberOfAngle_Stim(MonkeyIndex==2,:),1);

%{
NumAngleControl_R = sum(NumberOfAngle_Control(Sig&MonkeyIndex==1,:),1);

NumAngleStim_R = sum(NumberOfAngle_Stim(Sig&MonkeyIndex==1,:),1);

NumAngleControl_A = sum(NumberOfAngle_Control(Sig&MonkeyIndex==2,:),1);

NumAngleStim_A = sum(NumberOfAngle_Stim(Sig&MonkeyIndex==2,:),1);
%}
NumberOfAngle_Control_Sem_R = NaN * NumAngleControl_R';
NumberOfAngle_Stim_Sem_R = NaN * NumAngleStim_R';

NumberOfAngle_Control_Sem_A = NaN * NumAngleControl_A';
NumberOfAngle_Stim_Sem_A = NaN * NumAngleStim_A';

%Change to proportion to put two monkeys data together
PropAngleControl_R = (NumAngleControl_R)/sum(NumAngleControl_R);

PropAngleStim_R =(NumAngleStim_R)/sum(NumAngleStim_R);

PropAngleControl_A = (NumAngleControl_A)/sum(NumAngleControl_A);

PropAngleStim_A =(NumAngleStim_A)/sum(NumAngleStim_A);

PropAngleStim_All=NumAngleStim/sum(NumAngleStim);
PropAngleStim_Control=NumAngleControl/sum(NumAngleControl);
%{
%Export Data
 DataPath = '/Users/yux8/WorkRelated_YXF/Presentation_Summary/Optogenetics/Optogenetics_paper_preparation/NatureNeuroscienceSubmission/DataShare/ExportData';
    
 fileName5 = fullfile(DataPath, 'FEF_SaccadeAngle_OpticalStim_MonkeyR.mat');
 tmp5 =PropAngleStim_R;
 save(fileName5, 'tmp5');

 fileName6 = fullfile(DataPath, 'FEF_SaccadeAngle_OpticalStim_MonkeyA.mat');
 tmp6 =PropAngleStim_A;
 save(fileName6, 'tmp6');

 fileName7 = fullfile(DataPath, 'FEF_SaccadeAngle_Control_MonkeyR.mat');
 tmp7 =PropAngleControl_R;
 save(fileName7, 'tmp7');

 fileName8 = fullfile(DataPath, 'FEF_SaccadeAngle_Control_MonkeyA.mat');
 tmp8 =PropAngleControl_A;
 save(fileName8, 'tmp8');

 %keyboard
%}




SaccadeVector = [33.7500000000000	56.2500000000000	78.7500000000000	101.250000000000	123.750000000000	146.250000000000	168.750000000000	191.250000000000	213.750000000000	236.250000000000	258.750000000000	281.250000000000	303.750000000000	326.250000000000	348.750000000000	11.2500000000000];


NumberOfAngle_Control_Sem = NaN*NumAngleControl';
NumberOfAngle_Stim_Sem = NaN*NumAngleControl';
%{
SaccadeCountStim_Sum_all =sum(SaccadeCountStim_Sum(Sig,:),1);
SaccadeCountControl_Sum_all =sum(SaccadeCountControl_Sum(Sig,:),1);
%}

SaccadeCountStim_Sum_all =sum(SaccadeCountStim_Sum,1);
SaccadeCountControl_Sum_all =sum(SaccadeCountControl_Sum,1);

SaccadeCountStim_Sum_R =sum(SaccadeCountStim_Sum( MonkeyIndex==1,:),1);
SaccadeCountControl_Sum_R =sum(SaccadeCountControl_Sum( MonkeyIndex==1,:),1);

SaccadeCountStim_Sum_A =sum(SaccadeCountStim_Sum(MonkeyIndex==2,:),1);
SaccadeCountControl_Sum_A =sum(SaccadeCountControl_Sum( MonkeyIndex==2,:),1);

SaccadeCountStim_Prop_R =sum(SaccadeCountStim_Sum( MonkeyIndex==1,:),1)/sum(SaccadeCountStim_Sum( MonkeyIndex==1,:),'all');
SaccadeCountControl_Prop_R =sum(SaccadeCountControl_Sum( MonkeyIndex==1,:),1)/sum(SaccadeCountStim_Sum( MonkeyIndex==1,:),'all');

SaccadeCountStim_Prop_A =sum(SaccadeCountStim_Sum( MonkeyIndex==2,:),1)/sum(SaccadeCountStim_Sum(MonkeyIndex==2,:),'all');
SaccadeCountControl_Prop_A =sum(SaccadeCountControl_Sum( MonkeyIndex==2,:),1)/sum(SaccadeCountControl_Sum( MonkeyIndex==2,:),'all');


%{

SaccadeCountStim_Sum_R =sum(SaccadeCountStim_Sum(Sig & MonkeyIndex==1,:),1);
SaccadeCountControl_Sum_R =sum(SaccadeCountControl_Sum(Sig & MonkeyIndex==1,:),1);

SaccadeCountStim_Sum_A =sum(SaccadeCountStim_Sum(Sig & MonkeyIndex==2,:),1);
SaccadeCountControl_Sum_A =sum(SaccadeCountControl_Sum(Sig & MonkeyIndex==2,:),1);

SaccadeCountStim_Prop_R =sum(SaccadeCountStim_Sum(Sig & MonkeyIndex==1,:),1)/sum(SaccadeCountStim_Sum(Sig & MonkeyIndex==1,:),'all');
SaccadeCountControl_Prop_R =sum(SaccadeCountControl_Sum(Sig & MonkeyIndex==1,:),1)/sum(SaccadeCountStim_Sum(Sig & MonkeyIndex==1,:),'all');

SaccadeCountStim_Prop_A =sum(SaccadeCountStim_Sum(Sig & MonkeyIndex==2,:),1)/sum(SaccadeCountStim_Sum(Sig & MonkeyIndex==2,:),'all');
SaccadeCountControl_Prop_A =sum(SaccadeCountControl_Sum(Sig & MonkeyIndex==2,:),1)/sum(SaccadeCountControl_Sum(Sig & MonkeyIndex==2,:),'all');
%}

TotalSaccades_A_PSTH = sum(SaccadeCountStim_Sum( MonkeyIndex==1,:),'all') + sum(SaccadeCountControl_Sum(MonkeyIndex==1,:),'all');
TotalSaccades_R_PSTH = sum(SaccadeCountStim_Sum( MonkeyIndex==2,:),'all') + sum(SaccadeCountControl_Sum(MonkeyIndex==2,:),'all');

%keyboard
%%Plotting
if ShowFigureFlag
    %Figure 1: Saccade distribution among aong the whole trial periods
    FigureStartNum=20;
    FigureIndex=1;
    figtitlestr{FigureIndex}='Distribution_Of_saccade_angle';

    fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[50,100, 600,600],'Name',figtitlestr{FigureIndex});
%
   % set(fig{FigureIndex}, 'PaperUnits', 'inches');
   % formatFig(fig{FigureIndex}, [1.5 1.5], 'nature');
%}
    %Control
   % pax = polaraxes;
   % polarplot(pax, [SaccadeVector/180*pi SaccadeVector(1)/180*pi],[NumAngleControl NumAngleControl(1)], '-r', 'LineWidth', 2);
  %  hold on
  % h = polarwitherrorbar([SaccadeVector/180*pi SaccadeVector(1)/180*pi],[PropAngleStim_Control PropAngleStim_Control(1)],[NumberOfAngle_Control_Sem' NumberOfAngle_Control_Sem(1)],...
  %  '-k',3);
%
     h = polarwitherrorbar([SaccadeVector/180*pi SaccadeVector(1)/180*pi],[PropAngleControl_R PropAngleControl_R(1)],[NumberOfAngle_Control_Sem_R' NumberOfAngle_Control_Sem_R(1)],...
    '-k',3);
     hold on
     h = polarwitherrorbar([SaccadeVector/180*pi SaccadeVector(1)/180*pi],[PropAngleControl_A PropAngleControl_A(1)],[NumberOfAngle_Control_Sem_A' NumberOfAngle_Control_Sem_A(1)],...
    '-b',3);
%}
   
%h.MarkerSize = 0.5;

%stim
hold on
%s = polarwitherrorbar([SaccadeVector/180*pi SaccadeVector(1)/180*pi],[PropAngleStim_All PropAngleStim_All(1)],[NumberOfAngle_Stim_Sem' NumberOfAngle_Stim_Sem(1)],...
%    '-r',3);
%
s = polarwitherrorbar([SaccadeVector/180*pi SaccadeVector(1)/180*pi],[PropAngleStim_R PropAngleStim_R(1)],[NumberOfAngle_Stim_Sem_R' NumberOfAngle_Stim_Sem_R(1)],...
    '-r',3);
%s.MarkerSize = 0.5;
hold on
s = polarwitherrorbar([SaccadeVector/180*pi SaccadeVector(1)/180*pi],[PropAngleStim_A PropAngleStim_A(1)],[NumberOfAngle_Stim_Sem_A' NumberOfAngle_Stim_Sem_A(1)],...
    '-m',3);

%}



set(gca,'ThetaDir' , 'counterclockwise');
set(gca,'ThetaZeroLocation','top')
set(gca,'FontSize',30,'FontWeight','Bold','LineWidth',3);
%set(gca,'FontSize',30,'LineWidth',0.5);


thetaticks([0:45:359]);
thetaticklabels([0:45:359]);
rlim([0,0.8]);

%thetaticks(sort(SaccadeVector));
%thetaticklabels(sort(SaccadeVector));

%title("Saccade Distribution",'FontSize',25,'FontWeight','Bold');


%rlim([0,max(max([max(PropAngleStim_R,PropAngleStim_A),max(PropAngleControl_R,PropAngleControl_A)])*1.2)]);

%% Figure 5: Raster for saccades
SaccadeCountStim_Prop_R = SaccadeFreqStim_Sum_R_Mean;
SaccadeCountStim_Prop_A = SaccadeFreqStim_Sum_A_Mean;

SaccadeCountControl_Prop_R = SaccadeFreqControl_Sum_R_Mean;
SaccadeCountControl_Prop_A = SaccadeFreqControl_Sum_A_Mean;

FigureIndex=FigureIndex+1;
FigureStartNum=FigureStartNum+1;

figtitlestr{FigureIndex}='SaccadesRaster';
fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[100,100, 600,300],'Name',figtitlestr{FigureIndex});

set(fig{FigureIndex}, 'PaperUnits', 'inches');
    formatFig(fig{FigureIndex}, [1.5 1.5], 'nature');



%subplot(2,1,1);
axes('Position',[0.15,0.57,0.8,0.35]);
%maxy = max(max(SaccadeCountStim_Sum_all),max(SaccadeCountControl_Sum_all))+2;
maxy_R = max(max(SaccadeCountStim_Prop_R),max(SaccadeCountControl_Prop_R))+0.01;
maxy_A = max(max(SaccadeCountStim_Prop_A),max(SaccadeCountControl_Prop_A))+0.01;
maxy = max(maxy_R ,maxy_A);

AverageRegion = [0,100];
area(AverageRegion,[maxy,maxy],...
    'FaceColor','k','FaceAlpha',.1,'EdgeAlpha',.1,'HandleVisibility','off');

hold on

%{
plot(TimeSeq,SaccadeCountStim_Sum,'-r','LineWidth',3)
hold on
plot(TimeSeq,SaccadeCountControl_Sum,'-k','LineWidth',3)
%}
bar(TimeSeq,SaccadeCountStim_Sum_all,'r')
%stairs(TimeSeq,SaccadeCountStim_Prop_R,'r');
hold on
%stairs(TimeSeq,SaccadeCountStim_Prop_A,'m');
%set(gca,'LineWidth',3,'FontSize',20,'FontWeight','Bold');
set(gca,'LineWidth',0.5,'FontSize',5);%,'FontWeight','Bold');
 
 box off
max_stimtime = TimeSeq(SaccadeCountStim_Sum_all == max(SaccadeCountStim_Sum_all));
max_stimtime_R = TimeSeq(SaccadeCountStim_Prop_R == max(SaccadeCountStim_Prop_R));
maxstim = max(SaccadeCountStim_Sum_all)+20;
maxstim_R = max(SaccadeCountStim_Prop_R)+0.01;

max_stimtime_A = TimeSeq(SaccadeCountStim_Prop_A == max(SaccadeCountStim_Prop_A));

maxstim_A = max(SaccadeCountStim_Prop_A)+0.01;
%
for i = 1:length(max_stimtime )
  plot([max_stimtime(i),max_stimtime(i)],[maxstim,maxstim],'v','MarkerSize',5,'Color','r','MarkerFaceColor','r');
end
ylim([0,maxstim]);
%}
%{
for i = 1:length(max_stimtime_R )
  plot([max_stimtime_R(i),max_stimtime_R(i)],[max(maxstim_R,maxstim_A),max(maxstim_R,maxstim_A)],'v','MarkerSize',5,'Color','r','MarkerFaceColor','r');
end


for i = 1:length(max_stimtime_A )
  plot([max_stimtime_A(i),max_stimtime_A(i)],[max(maxstim_R,maxstim_A),max(maxstim_R,maxstim_A)],'v','MarkerSize',5,'Color','m','MarkerFaceColor','m');
end

ylim([0,max(maxstim_R,maxstim_A)]);
%}
%legend({'Stim'});
xticks([]);
xlim([-100,200]);

%subplot(2,1,2);
axes('Position',[0.15,0.2,0.8,0.35])
maxy = max(max(SaccadeCountStim_Sum_all),max(SaccadeCountControl_Sum_all))+2;
%maxy_R = max(max(SaccadeCountStim_Prop_R),max(SaccadeCountControl_Prop_R))+0.01;
%maxy_A = max(max(SaccadeCountStim_Prop_A),max(SaccadeCountControl_Prop_A))+0.01;
%maxy = max(maxy_R,maxy_A);
AverageRegion = [0,100];
area(AverageRegion,[maxy,maxy],...
    'FaceColor','k','FaceAlpha',.1,'EdgeAlpha',.1,'HandleVisibility','off');

hold on

bar(TimeSeq,SaccadeCountControl_Sum_all,'k')
%stairs(TimeSeq,SaccadeCountControl_Prop_R,'k');
%stairs(TimeSeq,SaccadeCountControl_Prop_A,'b');

%ylim([0,max(maxstim_R,maxstim_A)]);
xlim([-100,200]);
box off;
xlabel('Time from Opto-stim Onset (ms)'); 
ylabel('Number of saccades');

%set(gca,'LineWidth',3,'FontSize',20,'FontWeight','Bold');
set(gca,'LineWidth',0.5,'FontSize',5);

%legend({'Control'});

%% Figure 3 plot stim control together
%Use saccade frequency
SaccadeCountStim_Prop_R = SaccadeFreqStim_Sum_R_Mean;
SaccadeCountStim_Prop_A = SaccadeFreqStim_Sum_A_Mean;

SaccadeCountControl_Prop_R = SaccadeFreqControl_Sum_R_Mean;
SaccadeCountControl_Prop_A = SaccadeFreqControl_Sum_A_Mean;


FigureIndex=FigureIndex+1;
FigureStartNum=FigureStartNum+1;

figtitlestr{FigureIndex}='SaccadesRasterTogether';
fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[100,100, 600,300],'Name',figtitlestr{FigureIndex});

set(fig{FigureIndex}, 'PaperUnits', 'inches');
    formatFig(fig{FigureIndex}, [1.5 1.5], 'nature');

maxy_R = max(max(SaccadeCountStim_Prop_R),max(SaccadeCountControl_Prop_R))+0.01;
maxy_A = max(max(SaccadeCountStim_Prop_A),max(SaccadeCountControl_Prop_A))+0.01;
maxy = max(maxy_R ,maxy_A);

AverageRegion = [0,100];
area(AverageRegion,[maxy,maxy],...
    'FaceColor','k','FaceAlpha',.1,'EdgeAlpha',.1,'HandleVisibility','off');

hold on

%Control

AverageRegion = [0,100];
area(AverageRegion,[maxy,maxy],...
    'FaceColor','k','FaceAlpha',.1,'EdgeAlpha',.1,'HandleVisibility','off');

hold on

%bar(TimeSeq,SaccadeCountControl_Sum_all,'k')
stairs(TimeSeq,SaccadeCountControl_Prop_R,'k');
stairs(TimeSeq,SaccadeCountControl_Prop_A,'b');




%Stim
stairs(TimeSeq,SaccadeCountStim_Prop_R,'r');
hold on
stairs(TimeSeq,SaccadeCountStim_Prop_A,'m');
%set(gca,'LineWidth',3,'FontSize',20,'FontWeight','Bold');
set(gca,'LineWidth',0.5,'FontSize',5);%,'FontWeight','Bold');
 
 box off
%max_stimtime = TimeSeq(SaccadeCountStim_Sum_all == max(SaccadeCountStim_Sum_all));
max_stimtime_R = TimeSeq(SaccadeCountStim_Prop_R == max(SaccadeCountStim_Prop_R));
%maxstim = max(SaccadeCountStim_Sum_all)+20;
maxstim_R = max(SaccadeCountStim_Prop_R)+0.01;

max_stimtime_A = TimeSeq(SaccadeCountStim_Prop_A == max(SaccadeCountStim_Prop_A));

maxstim_A = max(SaccadeCountStim_Prop_A)+0.01;
%{
for i = 1:length(max_stimtime )
  plot([max_stimtime(i),max_stimtime(i)],[maxstim,maxstim],'v','MarkerSize',5,'Color','r','MarkerFaceColor','r');
end
ylim([0,maxstim]);
%}
for i = 1:length(max_stimtime_R )
  plot([max_stimtime_R(i),max_stimtime_R(i)],[max(maxstim_R,maxstim_A),max(maxstim_R,maxstim_A)],'v','MarkerSize',5,'Color','r','MarkerFaceColor','r');
end


for i = 1:length(max_stimtime_A )
  plot([max_stimtime_A(i),max_stimtime_A(i)],[max(maxstim_R,maxstim_A),max(maxstim_R,maxstim_A)],'v','MarkerSize',5,'Color','m','MarkerFaceColor','m');
end

ylim([0,max(maxstim_R,maxstim_A)]);

box off;
xlabel('Time from Opto-stim Onset (ms)'); 
ylabel('Saccade Frequency (Saccades/s)');

%set(gca,'LineWidth',3,'FontSize',20,'FontWeight','Bold');
set(gca,'LineWidth',0.5,'FontSize',5);
xlim([-100,200]);
%keyboard
%legend({'Stim'});
%xticks([]);

%keyboard
%% Figure 3 saccade bias distribution
%{
FigureStartNum=FigureStartNum+1;
    FigureIndex=FigureIndex+1;
    figtitlestr{FigureIndex}='SaccadeBiasPopulation';

    fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[50,100, 600,600],'Name',figtitlestr{FigureIndex});
bar(edge(2:end),num_total,'FaceColor','w','LineWidth',3);
hold on
bar(edge(2:end),num_sig,'FaceColor','r');
box off

xlabel('Proportion of contralateral saccades')
ylabel('Number of cases');
set(gca,'LineWidth',3,'FontWeight','Bold','FontSize',15);
end
%}
%% Output to files
%Setup output directory
Workingdirectory=pwd;
MarkerFolder='DataAnalysis';
%MarkerFolder='DataHub';
Flag=strfind(Workingdirectory,MarkerFolder);
BasicDirectory=Workingdirectory(1:Flag+length(MarkerFolder));


%%%%%%%%%%%Export data%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Set up output path
%OutputFolerName='ForagingTask';
OutPath=strcat(BasicDirectory,'Results',BasicDirectory(end),'ExportFigure',BasicDirectory(end));%,OutputFolerName);
if ~exist(OutPath)
    mkdir(OutPath)
end
cd(OutPath);
%{
FileName=Data.FileName;
FileName=FileName(1:end-4);


%Get File Information
num = regexp(FileName, '\d+', 'match');
RecordDateOriginal=cell2mat(num(1));

NeuronNum=cell2mat(num(2));
ChannelNum=SpikeChannel;


%To get the monkey name
RemainFiles=erase(FileName,RecordDateOriginal);
Delimiter=find(isstrprop(RemainFiles,'upper')==1);
MonkeyName=RemainFiles(Delimiter(1):Delimiter(2)-1);
%}

%Output file name
%OutputFileName=sprintf('%s_%s_N%s_C%s.mat',MonkeyName,string(RecordDateOriginal),string(NeuronNum),string(ChannelNum));
%OutputFigureName=sprintf('%s_%s_N%s_C%s_%CID_%NID',MonkeyName,string(RecordDateOriginal),string(NeuronNum),string(ChannelNum),string(ChannelID),string(ClusterID));

%OutputFigureName='SC_VideoVeiw_Behavior';
%OutputFigureName='SC_VideoVeiw_Behavior_Old';
OutputFigureName='FEF_VideoView_Behavior';


for jj=1:FigureIndex
        if ~isempty(fig{jj})
        

         OutputFigureName_Final=strcat(OutputFigureName,figtitlestr{jj},'.pdf');
         saveas(fig{jj},OutputFigureName_Final);
        end
        
            
end

    disp('Figures have been exported to the exported folder');
%keyboard
%{
ExistFlag=0;
%Load the old file if exist
if exist(OutputFileName)

    load(OutputFileName);
    %Load OutputData into the memory
    if isfield(OutputData,'FreeViewNeural')
    ExistFlag=1;
    end
end



%{
%FileName=strcat(FileName,'_RearchingTime.mat');
if SaveFlag
save(FileName,'OutputData');
disp('Seaching Time Data Exported');
else
    disp('No data export');
end
%}
%OutputData For current analysis
OutputData_New=[];

OutputData_New.FreeViewNeural.TrialType='GoodTrials';

%%%%%%%%%Save overall scene response characteristics%%%%%%%%%%%%%%%%%%%%%%%
%Organize the data


NameStr={'PSTH_Time','PSTH_Stim_Mean','PSTH_Control_Mean','AcStrength','ResponseLatency',...
    'PostModulationStrength','PostModulation_p','StimModulation','StimModulation_p',...
    'PSTH_FirstHalf','PSTH_SecondHalf',...
    'AS_First','AS_Second',...
    'PSTH_Stim_Mean_z','PSTH_Control_Mean_z','PSTH_Stim_Mean_First_z','PSTH_Stim_Mean_Second_z',...
    'StimModulation_z',...
    'PSTH_Stim_Mean_norm','PSTH_Control_Mean_norm','FR_Mean','p_StimControl',...
    'AcDur','psth_stim_normed2','psth_control_normed2','p_cri50','p_cri20','Strength_sig','p_whole_cri50',...
    'Latency3','Latency2'...
    'Raster_Stim','Raster_Control'
    };

DataLib={MeanPSTH_Time_Mean,SpikeCountStim_Mean,SpikeCountControl_Mean,ActivationStrength,ResponseLatency,...
    PostStimModulation,p_PostModulation,StimModulation,p_StimModulation,...
    PSTH_Stim_Mean_First,PSTH_Stim_Mean_Second,...
    ActivationStrength_FirstHalf,ActivationStrength_SecondHalf,...
    PSTH_Stim_Mean_z,PSTH_Control_Mean_z,PSTH_Stim_Mean_First_z,PSTH_Stim_Mean_Second_z,...
    StimModulation_z,...
    PSTH_Stim_Mean_norm_raw,PSTH_Control_Mean_norm_raw,FR_Mean,p_StimControl,...
    AcDur,psth_stim_normed2,psth_control_normed2,p_use_cri,p_use_cri20,Strength_sig,p_whole_cri50,...
    Latency3,Latency2,...
    Raster_Stim,Raster_Control



};



DataStamp=containers.Map(NameStr,DataLib);





OutputData_New.FreeViewNeural.Task=TaskType;
OutputData_New.FreeViewNeural.TaskCode=TaskCode;
OutputData_New.FreeViewNeural.DataStamp=DataStamp;

OutputData_New.FreeViewNeural.StartTrial =StartTrial;
%{
%Store figures
if ShowFigureFlag
   
     %OutputData_New.SceneTuning.Figures={fig99,fig100,fig101,fig102};
    %Output figures in a folder;
    FolderName=OutputFigureName;
    if ~exist(FolderName,'dir')
        mkdir(FolderName);
    end
    cd(FolderName);
    
    TotalFigure=length(fig);
    
    for jj=1:TotalFigure
        if ~isempty(fig{jj})
        
         OutputFigureName_Final=strcat(OutputFigureName,figtitlestr{jj},'.jpg');
         saveas(fig{jj},OutputFigureName_Final);
        end
        
            
        
    end

    disp('Figures have been exported to the neuron folder');

end

%}





if ExistFlag
%Compare the old one with the new one
Task_Old=OutputData.FreeViewNeural.TaskCode;

     if Task_Old~=TaskCode
        %Not the same task, add the new dataset into the old one
        OutputData(numel(OutputData)+1)=OutputData_New;
     elseif ismember(FileName,{'Adams062124VIDEOVIEW1'})
         %May save different task segment
         if StartTrial~=OutputData.FreeViewNeural.StartTrial

             if StartTrial > 160
                 OutputFileName=sprintf('%s_%s_N%s_C%s.mat',MonkeyName,string(RecordDateOriginal),string(NeuronNum+1),string(ChannelNum));
                 OutputData.FreeViewNeural=OutputData_New.FreeViewNeural;

             elseif StartTrial > 360
                  OutputFileName=sprintf('%s_%s_N%s_C%s.mat',MonkeyName,string(RecordDateOriginal),string(NeuronNum+2),string(ChannelNum));
                  OutputData.FreeViewNeural=OutputData_New.FreeViewNeural;

             end


         end

          
     else
         %If the same task, replace the old one with the new one
         %OutputData=OutputData_New;
         OutputData.FreeViewNeural=OutputData_New.FreeViewNeural;
         
    
    
     end
else
    %Output the current file
   OutputData.FreeViewNeural=OutputData_New.FreeViewNeural; 
end



if OutputFlag
    cd(OutPath);
    save(OutputFileName,'OutputData');
    disp('Data Exported');
else
    disp('No data export');
end  
%}

end

