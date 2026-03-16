function Batch_VideoFreeView_Neural(Data,  StartFile, EndFile,ShowFigureFlag,OutputFlag);
BatchFileName=Data.BatchFileName;
%Batch data file path
FilesName=Data.ResultFilePath(StartFile: EndFile);
RecordDate=Data.RecordDate(StartFile: EndFile);
FileNum=EndFile-StartFile+1;

FilesName_Each=Data.FileName;
ChannelID=Data.ChannelNumber;
%{
%Color Definition for each events
Color=[lines(6);prism(6);pink(6)];
%Color(2,:)=copper(1)+0.5;
Color=[Color;hsv(10)];
Color(18,:)=[0.1,0.457,0.194];
%}
%Load file
index = 1;

ChannelFull=[];
for i=1:length(FilesName)
    ChannelFull=[ChannelFull,ChannelID(i,1):ChannelID(i,2)];
  clear OutputData;
  for j = 1:length(FilesName{i})
      
     load(FilesName{i}{j});

     FilesNameAll{index}=FilesName_Each{i};
     
     OutputDataTemp=OutputData;
     
     Task(index)=OutputData.FreeViewNeural(1).TaskCode;
     TaskName{index}=OutputData.FreeViewNeural(1).Task;
     
     Data =OutputData.FreeViewNeural.DataStamp;

     %PSTH_Stim_Mean{index} = Data('PSTH_Stim_Mean');
     %PSTH_Stim_Mean{index} = Data('PSTH_Stim_Mean_norm');%Use the z scored ones
     %PSTH_Control_Mean{index} = Data('PSTH_Control_Mean_norm');%Use the z scored ones

     PSTH_Stim_Mean{index} = Data('PSTH_Stim_Mean');%Use the z scored ones
     PSTH_Control_Mean{index} = Data('PSTH_Control_Mean');%Use the z scored ones
     
     AS(index) = Data('AcStrength');
     Latency(index) = Data('ResponseLatency');
     p_AS(index) = Data('Ac_p');
     PMS(index) = Data('PostModulationStrength');
     PMS_p(index) = Data('PostModulation_p');
     SM(index) = Data('StimModulation');
     SM_p(index) = Data('StimModulation_p');
     %PSTH_First{index} = Data('PSTH_FirstHalf');
     %PSTH_Second{index} = Data('PSTH_SecondHalf');
     PSTH_First{index} = Data('PSTH_Stim_Mean_First_z');
     PSTH_Second{index} = Data('PSTH_Stim_Mean_Second_z');


     AS_First(index) = Data('AS_First');
     AS_Second(index) = Data('AS_Second');
     
     index = index +1;
      
      

  end
     
     
     
end

PSTH_Time = Data('PSTH_Time');
%{
Ex = AS > 0 & p_AS < 0.05;
In = AS < 0 & p_AS < 0.05;
%}
%
Ex = SM > 0 & SM_p < 0.05;
In = SM < 0 & SM_p < 0.05;
%}
PSTH = cell2mat(PSTH_Stim_Mean'); 
PSTH_control = cell2mat(PSTH_Control_Mean');





%{

figure

subplot(1,5,1);
screened_psth = PSTH_Diff(Ex_Fast_Sig ,:);

peak=max(screened_psth,[],2);
plot((PSTH_Diff(Ex_Fast_Sig ,:)./repmat(peak,1,size(screened_psth,2)))');

subplot(1,5,2);
plot(PSTH_Diff(Ex_Slow_Sig ,:)');

subplot(1,5,3);
plot(PSTH_Diff(In_Fast_Sig ,:)');

subplot(1,5,4);
plot(PSTH_Diff(In_Slow_Sig ,:)');

subplot(1,5,5);
plot(PSTH_Diff(Non_Sig_m ,:)');






keyboard
%}

%Ex_PSTH_Mean = median(PSTH(Ex,:),1) ;
%In_PSTH_Mean = median(PSTH(In,:),1) ;

PSTH()

Ex_PSTH_Mean = mean(PSTH(Ex,:),1) ;
In_PSTH_Mean = mean(PSTH(In,:),1) ;

Ex_PSTH_Sem = nanstd(PSTH(Ex,:),[],1)/sqrt(sum(Ex)) ;
In_PSTH_Sem = nanstd(PSTH(In,:),[],1)/sqrt(sum(Ex)) ;

Ex_PSTH_Mean_Control = mean(PSTH_control(Ex,:),1) ;
In_PSTH_Mean_Control = mean(PSTH_control(In,:),1) ;

Ex_PSTH_Sem_Control = nanstd(PSTH_control(Ex,:),[],1)/sqrt(sum(Ex)) ;
In_PSTH_Sem_Control = nanstd(PSTH_control(In,:),[],1)/sqrt(sum(Ex)) ;



N_Ex = sum(Ex);
N_In = sum(In);
disp(sprintf('Ex n = %d',N_Ex));
disp(sprintf('N_In = %d',N_In));

Latency_Ex = Latency(Ex);
Latency_In = Latency(In);
%{
figure
subplot(1,2,1);
plot(PSTH(Ex,:)');
subplot(1,2,2);
plot(PSTH(In,:)');
%}
First_PSTH = cell2mat(PSTH_First');
Second_PSTH = cell2mat(PSTH_Second');

First_Ex = First_PSTH(Ex,:);
Second_Ex = Second_PSTH(Ex,:);

First_In = First_PSTH(In,:);
Second_In = Second_PSTH(In,:);

First_Ex_Mean = mean(First_Ex,1);
Second_Ex_Mean = mean(Second_Ex,1);

First_In_Mean = mean(First_In,1);
Second_In_Mean = mean(Second_In,1);

First_Ex_Sem = nanstd(First_Ex,[],1)/sqrt(sum(Ex)) ;
Second_Ex_Sem = nanstd(Second_Ex,[],1)/sqrt(sum(Ex)) ;

First_In_Sem = nanstd(First_In,[],1)/sqrt(sum(In)) ;
Second_In_Sem = nanstd(Second_In,[],1)/sqrt(sum(In)) ;

AS_First_Ex_Mean = mean(AS_First(Ex));
AS_First_In_Mean = mean(AS_First(In));

AS_Second_Ex_Mean = mean(AS_Second(Ex));
AS_Second_In_Mean = mean(AS_Second(In));

AS_First_Ex_Sem = std(AS_First(Ex))/sqrt(sum(Ex));
AS_First_In_Sem = std(AS_First(In))/sqrt(sum(In));


AS_Second_Ex_Mean = mean(AS_Second(Ex));
AS_Second_In_Mean = mean(AS_Second(In));

AS_Second_Ex_Sem = std(AS_Second(Ex))/sqrt(sum(Ex));
AS_Second_In_Sem = std(AS_Second(In))/sqrt(sum(In));

[h,p_ex] = ttest(AS_First(Ex),AS_Second(Ex))
[h,p_in] = ttest(AS_First(In),AS_Second(In))




if ShowFigureFlag
    
FigureStartNum=100;
FigureIndex=1;
    
figtitlestr{FigureIndex}='PSTH_Stim_Population';
fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[50,100, 500,600],'Name',figtitlestr{FigureIndex});
subplot(2,1,1)
AverageRegion = [0,100];
area(AverageRegion,[1.5,1.5],...
    'FaceColor','k','FaceAlpha',.1,'EdgeAlpha',.1,'HandleVisibility','off');
hold on
shadedErrorBar(PSTH_Time,Ex_PSTH_Mean_Control,Ex_PSTH_Sem_Control,'lineprops',{[0,0,0]},'transparent',1,'patchSaturation',0.3)


 shadedErrorBar(PSTH_Time,Ex_PSTH_Mean,Ex_PSTH_Sem,'lineprops',{[1,0,0]},'transparent',1,'patchSaturation',0.3)

 set(findall(gca, 'Type', 'Line'),'LineWidth',3);
 set(gca,'LineWidth',3,'FontSize',20,'FontWeight','Bold');
%ylim([0,160]);
%ylim([-1,1.5]);
box off
 subplot(2,1,2)
AverageRegion = [0,100];
area(AverageRegion,[1.5,1.5],...
    'FaceColor','k','FaceAlpha',.1,'EdgeAlpha',.1,'HandleVisibility','off');
hold on
    shadedErrorBar(PSTH_Time,In_PSTH_Mean_Control,In_PSTH_Sem_Control,'lineprops',{[0,0,0]},'transparent',1,'patchSaturation',0.3)
 shadedErrorBar(PSTH_Time,In_PSTH_Mean,In_PSTH_Sem,'lineprops',{[0,0,1]},'transparent',1,'patchSaturation',0.3)

 set(findall(gca, 'Type', 'Line'),'LineWidth',3);
%ylim([-1,1.5])
set(gca,'LineWidth',3,'FontSize',20,'FontWeight','Bold');
box off


%% 
FigureStartNum=FigureStartNum+1;
FigureIndex=FigureIndex+1;
    
figtitlestr{FigureIndex}='ModulationStrength';
fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[50,100, 500,600],'Name',figtitlestr{FigureIndex});
histogram(SM(SM_p<0.05),'BinWidth',10,'FaceColor','r','LineWidth',3);
set(gca,'LineWidth',3,'FontSize',20,'FontWeight','Bold');
box off;
xlim([-50,250]);
%% 

FigureStartNum=FigureStartNum+1;
FigureIndex=FigureIndex+1;
    
figtitlestr{FigureIndex}='ResponseLatency';
fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[50,100, 500,600],'Name',figtitlestr{FigureIndex});
subplot(2,1,1)
histogram(Latency_Ex,'BinWidth',5,'FaceColor','r','LineWidth',3);
set(gca,'LineWidth',3,'FontSize',20,'FontWeight','Bold');

box off;
xlim([-5,100]);


subplot(2,1,2)
histogram(Latency_In,'BinWidth',5,'FaceColor','b','LineWidth',3);
set(gca,'LineWidth',3,'FontSize',20,'FontWeight','Bold');

box off;
xlim([-5,100]);

%% 
FigureStartNum=FigureStartNum+1;
FigureIndex=FigureIndex+1;
    
figtitlestr{FigureIndex}='FirstHalfSecondHalf';
fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[50,100, 500,600],'Name',figtitlestr{FigureIndex});

subplot(2,1,1)
AverageRegion = [0,100];
area(AverageRegion,[1.5,1.5],...
    'FaceColor','k','FaceAlpha',.1,'EdgeAlpha',.1,'HandleVisibility','off');
hold on

 shadedErrorBar(PSTH_Time,First_Ex_Mean,First_Ex_Sem,'lineprops',{[1,0,0]},'transparent',1,'patchSaturation',0.3)
 hold on
 shadedErrorBar(PSTH_Time,Second_Ex_Mean,Second_Ex_Sem,'lineprops',{'#FFB000'},'transparent',1,'patchSaturation',0.3)

 set(findall(gca, 'Type', 'Line'),'LineWidth',3);
 set(gca,'LineWidth',3,'FontSize',20,'FontWeight','Bold');
ylim([-1,1.5]);
legend({'FirstHalf','LaterHalf'});
box off
 subplot(2,1,2)
AverageRegion = [0,100];
area(AverageRegion,[1.5,1.5],...
    'FaceColor','k','FaceAlpha',.1,'EdgeAlpha',.1,'HandleVisibility','off');
hold on

 shadedErrorBar(PSTH_Time,First_In_Mean,First_In_Sem,'lineprops',{'#6543FF'},'transparent',1,'patchSaturation',0.3)
 hold on
 shadedErrorBar(PSTH_Time,Second_In_Mean,Second_In_Sem,'lineprops',{'#648FFF'},'transparent',1,'patchSaturation',0.3)

 set(findall(gca, 'Type', 'Line'),'LineWidth',3);
ylim([-1,1.5]);
legend({'FirstHalf','LaterHalf'});
set(gca,'LineWidth',3,'FontSize',20,'FontWeight','Bold');
box off

%%
FigureStartNum=FigureStartNum+1;
FigureIndex=FigureIndex+1;
    
figtitlestr{FigureIndex}='FirstHalfSecondHalf';
fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[50,100, 100,200],'Name',figtitlestr{FigureIndex});

subplot(2,1,1)
bar([1,2],[AS_First_Ex_Mean,AS_Second_Ex_Mean],'red','LineWidth',3);
hold on 
errorbar([1,2],[AS_First_Ex_Mean,AS_Second_Ex_Mean],[AS_First_Ex_Sem,AS_Second_Ex_Sem],'.r','LineWidth',3);
set(gca,'LineWidth',3);
xlabel('Excitation');
box off

subplot(2,1,2)
bar([1,2],[AS_First_In_Mean,AS_Second_In_Mean],'b');
hold on 
errorbar([1,2],[AS_First_In_Mean,AS_Second_In_Mean],[AS_First_In_Sem,AS_Second_In_Sem],'.b','LineWidth',3);
set(gca,'LineWidth',3);
xlabel('Inhibition');
box off

end
end