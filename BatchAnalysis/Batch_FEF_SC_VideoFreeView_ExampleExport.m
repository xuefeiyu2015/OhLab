function Batch_FEF_SC_VideoFreeView_ExampleExport(Data,  StartFile, EndFile,ShowFigureFlag,OutputFlag);

BatchFileName=Data.BatchFileName;
%Batch data file path
FilesName=Data.ResultFilePath(StartFile: EndFile);
RecordDate=Data.RecordDate(StartFile: EndFile);
FileNum=EndFile-StartFile+1;
%{
%Color Definition for each events
Color=[lines(6);prism(6);pink(6)];
%Color(2,:)=copper(1)+0.5;
Color=[Color;hsv(10)];
Color(18,:)=[0.1,0.457,0.194];
%}
%Load file


FilesName_Each=Data.FileName(StartFile: EndFile);
ChannelID=Data.ChannelNumber(StartFile: EndFile,:);

index =1;

ChannelFull=[];

FileIndex=[];
AbsoluteIndex = [];


for i=1:length(FilesName)
  clear OutputData;
   ChannelFull=[ChannelFull,ChannelID(i,1):ChannelID(i,2)];

    for j = 1:length(FilesName{i})

        %For Video Free View
 
     load(FilesName{i}{j});
  
     
     OutputDataTemp=OutputData;

     FilesNameAll{index}=FilesName_Each{i};

     AbsoluteIndex(index)=index;

     if index == 1
         FileIndex(index) = 1;

     else
         if ~strcmp(FilesNameAll{index},FilesNameAll{index-1})
             FileIndex(index) = FileIndex(index-1)+1;
         else
             FileIndex(index) = FileIndex(index-1);

         end
     end





      Data =OutputData.FreeViewNeural.DataStamp;
%
    PSTH_Stim_Norm{index} = Data('PSTH_Stim_Mean_norm');
     PSTH_Control_Norm{index} = Data('PSTH_Control_Mean_norm');
     PSTH_Stim_Mean_z{index} = Data('PSTH_Stim_Mean_z');

   %   PSTH_Stim_Mean_z{index} = Data('PSTH_Stim_Mean_z');%Use the z scored ones
   %  PSTH_Control_Mean_z{index} = Data('PSTH_Control_Mean_z');%Use the z scored ones
%}
  %    PSTH_Stim_Norm{index} = Data('PSTH_Stim_Mean_z');
  %   PSTH_Control_Norm{index} = Data('PSTH_Control_Mean_z');

      PSTH_Stim_Mean{index} = Data('PSTH_Stim_Mean');
     PSTH_Control_Mean{index} = Data('PSTH_Control_Mean');

     OpticalStim_Mean(index) = mean(PSTH_Stim_Mean{index},'all');
     ControlStim_Mean(index) = mean(PSTH_Control_Mean{index},'all');

     FR_Mean_Video(index) = Data('FR_Mean');

     p_StimControl(index) = Data('p_StimControl');

     p_Cri20(index) = Data('p_cri20');
     p_Cri50(index) = Data('p_cri50');
     Latency2(index) = Data('Latency2');

     PSTH_Stim_Norm2{index} = Data('psth_stim_normed2');


   

     
     AS(index) = Data('AcDur');
     Latency(index) = Data('ResponseLatency');
    % p_AS(index) = Data('Ac_p');
     PMS(index) = Data('PostModulationStrength');
     PMS_p(index) = Data('PostModulation_p');
     SM(index) = Data('StimModulation');
     SM_p(index) = Data('StimModulation_p');
     %PSTH_First{index} = Data('PSTH_FirstHalf');
     %PSTH_Second{index} = Data('PSTH_SecondHalf');
     PSTH_First{index} = Data('PSTH_Stim_Mean_First_z');
     PSTH_Second{index} = Data('PSTH_Stim_Mean_Second_z');

     p_whole_cri50{index} = Data('p_whole_cri50');



%      SM_norm(index) = Data('StimulationStrength_norm');


     AS_First(index) = Data('AS_First');
     AS_Second(index) = Data('AS_Second');

     PSTH_Opto_Time = Data('PSTH_Time');

 

     %Memory saccade
     if isfield(OutputData,'MemorySaccade')
     
     MemorySaccadeIndex(index) = 1;
     
     Data =OutputData.MemorySaccade(1).DataStamp;

     %PSTH_Tgt_RF{index} = Data('PSTH_TgtOn_RF');
     %PSTH_Sac_RF{index} = Data('PSTH_SacOn_RF');
     
    % maxFR =  max(max(Data('PSTH_TgtOn_RF_z')),max(Data('PSTH_SacOn_RF_z')));
     PSTH_Tgt_RF{index} = Data('PSTH_TgtOn_RF_z');%/maxFR;
     PSTH_Sac_RF{index} = Data('PSTH_SacOn_RF_z');%/maxFR;

     Mean_FR_MS(index) = mean([nanmean( Data('PSTH_TgtOn_RF'),'all'),nanmean(Data('PSTH_SacOn_RF'),'all')]);
      Mean_FR_MS_z(index) = mean([nanmean( Data('PSTH_TgtOn_RF_z'),'all'),nanmean(Data('PSTH_SacOn_RF_z'),'all')]);


       PSTH_Tgt_RF_norm{index} = Data('PSTH_RF_Mean_Norm');%/maxFR;
     PSTH_Sac_RF_norm{index} = Data('PSTH_RF_Sac_Mean_Norm');%/maxFR;

      




     
try
    RF_h_v(index) = Data('RF_h_v');
    RF_h_m(index) = Data('RF_h_m');
catch
    tmp = Data('RF_h_v');
    RF_h_v(index) =tmp(1);
    tmp = Data('RF_h_m');
    RF_h_m(index)=tmp(1);
end

    
     %
     if RF_h_v(index)==1 & (RF_h_m(index)==0|RF_h_m(index)==-1)

         CellType(index) = 1;%Visual-Only

     elseif RF_h_v(index)==1 & RF_h_m(index)==1
         CellType(index) = 2;% Visual-Motor

      elseif (RF_h_v(index)==0|RF_h_v(index)==-1) & RF_h_m(index)==1
         CellType(index) = 3;% Motor-Only

      elseif RF_h_v(index)==0 & RF_h_m(index)==0

         CellType(index) = 5;%Null Type 
     else
         CellType(index) = 4; %Other type
     end
    
     

     PSTH_Time_TgtOn = Data('PSTH_Time_TgtOn');
     PSTH_Time_SacOn = Data('PSTH_Time_SacOn');
     else

         MemorySaccadeIndex(index) = 0;

     PSTH_Tgt_RF{index} = NaN*ones(1,140);
     PSTH_Sac_RF{index} = NaN*ones(1,140);


     

     RF_h_v(index) = NaN;
     RF_h_m(index) = NaN;

     CellType(index) =NaN;%1:visual;2:motor;3:vis-motor

   




     end %End of memory saccade

     %% Delay saccade
%{   
      if isfield(OutputData,'DelaySaccadeTuning')

          DelaySaccadeIndex(index) = 1;

          Data =OutputData.DelaySaccadeTuning(1).DataStamp;

          PSTH_TgtTime_DS = Data('PSTH_FR_TgtOn_Time');
          PSTH_SacTime_DS = Data('PSTH_FR_SacOn_Time');


         
          Prefer_SaccadeVector(index,:) = Data('Pref_Vect');%PreDurPost
         
          Prefer_TgtVector(index) = Data('Pref_TgtVect');

          p_SacV_prefer(index,:) = Data('SacVectTuning_p');

          p_TgtV_prefer(index) = Data('TgtVectTuning_p');

          TuningStrength(index) = Data('Pref_TgtVect_Length');

          SacTuningStrength(index,:) = Data('Pref_Vect_Length');

          VisTuning_Mean(index,:)=Data('TgtVectTuning_mean')/max(Data('TgtVectTuning_mean'));
          VisTuning_Sem(index,:)=Data('TgtVectTuning_sem')/max(Data('TgtVectTuning_mean'));


          tmp = Data('SacVectTuning_PreDurPost_mean');

          SacTuning_Mean(index,:) = tmp(:,1)'/max(tmp(:,1)');

          tmp = Data('SacVectTuning_PreDurPost_sem');

          SacTuning_Sem(index,:) = tmp(:,1)'/max(tmp(:,1)');
          



psth_tgt_ds(index,:) = Data('PSTH_TgtOn_RF_z');
psth_sac_ds(index,:) = Data('PSTH_SacOn_RF_z');


          






      else
          DelaySaccadeIndex(index) = 0;

           Prefer_SaccadeVector(index,:) =NaN*ones(1,3);%PreDurPost
          Prefer_TgtVector(index) =NaN;

          p_SacV_prefer(index,:) = NaN*ones(1,3);

          p_TgtV_prefer(index) = NaN;

          TuningStrength(index)= NaN;

          SacTuningStrength(index,:)=NaN*ones(1,3);

          psth_tgt_ds(index,:) = NaN*ones(1,817);
          psth_sac_ds(index,:) = NaN*ones(1,180);

      end %End of delay saccade 
%}

     index = index+1;   
    end %End of for file name
     
end

%Plot controller
PlotOptoEffect = 1;
PlotMemorySaccade = 1;
PlotDelaySaccade = 0;


%Analysis controller 
AnalysisOS = 1;% For Video Free View
AnalysisMS = 1;% For memory saccade
AnalysisDS = 0;%For Delay saccade

%Selection 
%ManuelOut
  A=[OpticalStim_Mean',ControlStim_Mean',FR_Mean_Video'];

CellIndex = 1:length(FR_Mean_Video);
%ManuelCheck=[550 469	478	477	473	479	425	408	353	351	320	239	240	217	213	201	178	163	160	139	141	68	43	48 115 132 137 157 181 200 221 496 516 529 538 581 623 97 131 159 315 325 329 343 350 400 427 435];
%ManuelCheck=[6,18,20:37];
ManuelCheck=[];
ManuelOut = ~ismember(CellIndex,ManuelCheck);


Cri1 = OpticalStim_Mean>1 & abs(Mean_FR_MS_z)>0& Mean_FR_MS>1;



ScreenedIndex = CellIndex(Cri1 & ManuelOut);

ScreenCell = ismember(CellIndex,ScreenedIndex);

ScreenCellNo = CellIndex(ScreenCell);


%% For VideoFreeView

PSTH = cell2mat(PSTH_Stim_Mean(ScreenCell)'); 
PSTH_control = cell2mat(PSTH_Control_Mean(ScreenCell)');

PSTH_norm = cell2mat(PSTH_Stim_Norm(ScreenCell)'); 
PSTH_control_norm = cell2mat(PSTH_Control_Norm(ScreenCell)');


PSTH_norm2 =cell2mat(PSTH_Stim_Norm2(ScreenCell)');



PSTH_z = cell2mat(PSTH_Stim_Mean_z(ScreenCell)');


p_whole_cri50 = cell2mat(p_whole_cri50);

%maxnum_ms = max(cellfun(@numel,PSTH_Tgt_RF));

%% For Memory saccade

PSTH_Tgt_RF = cell2mat(cellfun(@(x) x(40:140), PSTH_Tgt_RF(ScreenCell),'UniformOutput',false)');
PSTH_Sac_RF = cell2mat(PSTH_Sac_RF(ScreenCell)');

PSTH_Time_TgtOn = PSTH_Time_TgtOn(40:140);
%PSTH_Time_SacOn = PSTH_Time_SacOn(1:140);

 PSTH_Tgt_RF_norm=cell2mat(cellfun(@(x) x(40:140), PSTH_Tgt_RF_norm(ScreenCell),'UniformOutput',false)');
 PSTH_Sac_RF_norm = cell2mat(cellfun(@(x) x(1:140), PSTH_Sac_RF_norm(ScreenCell),'UniformOutput',false)');
%{
p_Cri20(index) = Data('p_cri50');
     p_Cri50(index) = Data('p_Cri50');
     Latency2(index) = Data('Latency2');

     PSTH_Stim_Norm2(index) = Data('psth_stim_normed2');

%}

%Other information

%SM,SM_p,Latency,ChannelNO, CellType,FR_Mean_Video

%All_Info = [SM',SM_p',Latency',ChannelFull', CellType',FR_Mean_Video',MemorySaccadeIndex',DelaySaccadeIndex',Prefer_TgtVector',p_TgtV_prefer', TuningStrength',AS',p_Cri20',p_Cri50',Latency2',AbsoluteIndex',p_whole_cri50',PMS_p'];
All_Info = [SM',SM_p',Latency',ChannelFull', CellType',FR_Mean_Video',MemorySaccadeIndex',AS',p_Cri20',p_Cri50',Latency2',AbsoluteIndex',p_whole_cri50',PMS_p'];

All_Info_Screened = All_Info(ScreenCell,:);

FilesNameScreened =FilesNameAll(ScreenCell);

SM = All_Info_Screened(:,1);
SM_p= All_Info_Screened(:,2);

Latency = All_Info_Screened(:,3);
ChannelScreened = All_Info_Screened(:,4);
CellType = All_Info_Screened(:,5);

% CellType(find(ChannelScreened==21,1,'first'))=2;

FR_Mean_Video = All_Info_Screened(:,6);

MemorySaccadeIndex = All_Info_Screened(:,7);

%DelaySaccadeIndex = All_Info_Screened(:,8);

%Prefer_TgtVector = All_Info_Screened(:,9);

%p_TgtV_prefer = All_Info_Screened(:,10);

%TuningStrength =  All_Info_Screened(:,11);

%{
AS =  All_Info_Screened(:,12);

p_Cri20 = All_Info_Screened(:,13);

p_Cri50 = All_Info_Screened(:,14);

Latency2 =  All_Info_Screened(:,15);

AbsoluteIndex = All_Info_Screened(:,16);

p_whole_cri50= All_Info_Screened(:,17);


PMS_p = All_Info_Screened(:,18);
%}

AS =  All_Info_Screened(:,8);

p_Cri20 = All_Info_Screened(:,9);

p_Cri50 = All_Info_Screened(:,10);

Latency2 =  All_Info_Screened(:,11);

AbsoluteIndex = All_Info_Screened(:,12);

p_whole_cri50= All_Info_Screened(:,13);


PMS_p = All_Info_Screened(:,14);

%SM_norm=SM_norm(ScreenCell)';




%Prefer_SaccadeVector =  Prefer_SaccadeVector(ScreenCell,:);

%p_SacV_prefer =  p_SacV_prefer(ScreenCell,:);
%SacTuningStrength_All =  SacTuningStrength(ScreenCell,:);

%% Delay saccade
%{
psth_tgt_ds=psth_tgt_ds(ScreenCell,:);
psth_sac_ds=psth_sac_ds(ScreenCell,:);

%}
   

%%
%Separate cells into different modulation groups according to laser
%response

%Proportion of SC neurons which is responsive to FEF optical stimulation
%Sig_m = (p_whole_cri50<0.05 | SM_p<0.05 | PMS_p<0.05)&(~isnan(Latency));
Sig_m = SM_p<0.05 &(~isnan(Latency));

FastLat = 20;

FastResponse = Latency < FastLat;
SlowResponse = ~FastResponse;
%FastResponse(4)=1;
%Sig_m(4)=1;
%SlowResponse(4)=1;
SlowResponse(10)=0;
FastResponse(10)=1;





Fast_Ex_Sel = Sig_m & FastResponse & SM>0;
Slow_Ex_Sel = Sig_m & SlowResponse & SM>0;
Fast_In_Sel = Sig_m & FastResponse & SM<0;
Slow_In_Sel = Sig_m & SlowResponse & SM<0;



psth_use = PSTH_norm2;


%Memory saccade
PSTH_Time_TgtOn=PSTH_Time_TgtOn-50;
PSTH_TgtOn_use = PSTH_Tgt_RF_norm*10;
PSTH_TgtOn_use(5,:)=PSTH_TgtOn_use(5,:)*8;
PSTH_SacOn_use = PSTH_Sac_RF_norm*10;




psth_v = (flipud(PSTH_TgtOn_use));
psth_s = (flipud(PSTH_SacOn_use));
PannelDist = 5;


AddInterval_v = 0:PannelDist:size(psth_v,1)*PannelDist-1;


psth_v = psth_v+AddInterval_v';




AddInterval_s = 0:PannelDist:size(psth_s,1)*PannelDist-1;
psth_s = psth_s+AddInterval_s';










%Select neurons 
 %PSTH_Tgt_RF
 %PSTH_Sac_RF
 %CellType

StepSize = 5;

 %For the plots
 if ShowFigureFlag == 1

     FigureStartNum=90;
FigureIndex=1;
    
figtitlestr{FigureIndex}='PSTH_AcrossChannel';
fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[50,100, 1200,600],'Name',figtitlestr{FigureIndex});

%% PSTH heatmap
subplot(1,4,1);
imagesc(psth_use);
%imagesc(psth_z);


%xticks(1:50:size(PSTH_Stim_Mean_norm,2));
TimePoint = [-100,0,100,200];
for i =1:length(TimePoint)
    TimeBin_i = PSTH_Opto_Time-TimePoint(i)==StepSize/2;
    try
    TimeBin(i) = find(TimeBin_i==1,1,'first');
    catch
      TimeBin(i) = NaN;
    end

end
xticks(TimeBin(~isnan(TimeBin)));
xticklabels( string(TimePoint(~isnan(TimeBin) )));
yticks(1:size(psth_use,1));
set(gca, 'CLim', [-1, 1],'FontSize',3,'LineWidth',0.3);
% create blue to red colormap:

cmap = nan(128,3);

r = [0.8 0 0];

w = [0.8 0.8 0.8];

b = [0 0 0.8];

for ii = 1:3

    cmap(:,ii) = [ fliplr(linspace(w(ii), b(ii), 64)),fliplr(linspace(r(ii), w(ii), 64))];

end

cmap(65,:) = [];
colormap(cmap);
%colorbar;

%% PSTH raw FR
subplot(1,4,2);
psth = (flipud(psth_use))*10;
PannelDist = 5;

AddInterval = 0:PannelDist:size(psth,1)*PannelDist-1;

psth = psth+AddInterval';
%psth = flipud(psth);

SpikeChannelWhole = 1:size(psth,1);
%
AverageRegion = [0,100];
area(AverageRegion,[max(psth(length(SpikeChannelWhole),:))+PannelDist,max(psth(length(SpikeChannelWhole),:))+PannelDist],...
    'FaceColor','k','FaceAlpha',.1,'EdgeAlpha',.1,'HandleVisibility','off');
hold on
%}
time_plot =repmat(PSTH_Opto_Time,size(psth,1),1);
plot(time_plot',psth','-k','LineWidth',2);
hold on

%plot(time_plot(flipud(SigFastChannel==1),:)',psth(flipud(SigFastChannel==1),:)','-r','LineWidth',3);



%plot(time_plot(flipud(SignificanceChannel==1& ModType ==1),:)',psth(flipud(SignificanceChannel==1& ModType ==1),:)','-r','LineWidth',3);
%plot(time_plot(flipud(SignificanceChannel==1& ModType ==-1),:)',psth(flipud(SignificanceChannel==1& ModType ==-1),:)','-b','LineWidth',3);



plot(time_plot(flipud(Fast_Ex_Sel),:)',psth(flipud(Fast_Ex_Sel),:)','-r','LineWidth',3);
plot(time_plot(flipud(Slow_Ex_Sel),:)',psth(flipud(Slow_Ex_Sel),:)','-y','LineWidth',3);
plot(time_plot(flipud(Fast_In_Sel),:)',psth(flipud(Fast_In_Sel),:)','-b','LineWidth',3);
plot(time_plot(flipud(Slow_In_Sel),:)',psth(flipud(Slow_In_Sel),:)','-g','LineWidth',3);

ylim([-2,max(psth,[],'all')+1]);
xlim([-100,200]);
yallrange = ylim;
ylim([-5,100]);
%axis off

%% Memory saccade 






%TimePlot = [-200,600];
%TimePlotSac = [-600,400];
%{
PSTH_TgtOn_use = psth_tgt_ds(:,PSTH_TgtTime_DS>TimePlot(1)&PSTH_TgtTime_DS<TimePlot(2));
PSTH_SacOn_use = psth_sac_ds(:,PSTH_Time_SacOn>TimePlotSac(1)&PSTH_Time_SacOn<TimePlotSac(2));

PSTH_Time_TgtOn=PSTH_TgtTime_DS(PSTH_TgtTime_DS>TimePlot(1)&PSTH_TgtTime_DS<TimePlot(2));
PSTH_Time_SacOn=PSTH_SacTime_DS(PSTH_Time_SacOn>TimePlotSac(1)&PSTH_Time_SacOn<TimePlotSac(2));
%}
%
%Use normalized responses

%MaxFR = max(max(PSTH_Tgt_RF(:,PSTH_Time_TgtOn>0 & PSTH_Time_TgtOn<100),[],2),max(PSTH_Sac_RF(:,PSTH_Time_SacOn>-100 & PSTH_Time_SacOn<0),[],2));






%psth_v=[psth_v,NaN*ones(size(psth_v,1),BlankInterval),psth_s];

% Tgt onset
subplot(1,4,3);
AverageRegion =[0,100];
area(AverageRegion,[max(psth_v(length(SpikeChannelWhole),:))+PannelDist,max(psth_v(length(SpikeChannelWhole),:))+PannelDist],...
    'FaceColor','k','FaceAlpha',.1,'EdgeAlpha',.1,'HandleVisibility','off');
hold on
%{
AverageRegion_Sac = PSTH_Time_TgtOn(end)+10+1-PSTH_Time_SacOn(1)+[0,10];
area(AverageRegion_Sac,[max(psth_v(length(SpikeChannelWhole),:))+PannelDist,max(psth_v(length(SpikeChannelWhole),:))+PannelDist],...
    'FaceColor','k','FaceAlpha',.1,'EdgeAlpha',.1,'HandleVisibility','off');
%}





time_plot =repmat(PSTH_Time_TgtOn,size(psth_v,1),1);
plot(time_plot',psth_v','-k','LineWidth',2);
plot(time_plot(flipud(Fast_Ex_Sel),:)',psth_v(flipud(Fast_Ex_Sel),:)','-r','LineWidth',3);
plot(time_plot(flipud(Slow_Ex_Sel),:)',psth_v(flipud(Slow_Ex_Sel),:)','-y','LineWidth',3);
plot(time_plot(flipud(Fast_In_Sel),:)',psth_v(flipud(Fast_In_Sel),:)','-b','LineWidth',3);
plot(time_plot(flipud(Slow_In_Sel),:)',psth_v(flipud(Slow_In_Sel),:)','-g','LineWidth',3);

%hold on
ylim([-2,max(max(psth_s,[],'all'),max(psth_v,[],'all'))+1]);

%ylim(yallrange);
ylim([-5,100]);
%axis off



% saccade onset
subplot(1,4,4);
AverageRegion =[0,10];
area(AverageRegion,[max(psth_s(length(SpikeChannelWhole),:))+PannelDist,max(psth_s(length(SpikeChannelWhole),:))+PannelDist],...
    'FaceColor','k','FaceAlpha',.1,'EdgeAlpha',.1,'HandleVisibility','off');
hold on
%{
AverageRegion_Sac = PSTH_Time_TgtOn(end)+10+1-PSTH_Time_SacOn(1)+[0,10];
area(AverageRegion_Sac,[max(psth_v(length(SpikeChannelWhole),:))+PannelDist,max(psth_v(length(SpikeChannelWhole),:))+PannelDist],...
    'FaceColor','k','FaceAlpha',.1,'EdgeAlpha',.1,'HandleVisibility','off');


%}



time_plot =repmat(PSTH_Time_SacOn,size(psth_s,1),1);
plot(time_plot',psth_s','-k','LineWidth',2);
plot(time_plot(flipud(Fast_Ex_Sel),:)',psth_s(flipud(Fast_Ex_Sel),:)','-r','LineWidth',3);
plot(time_plot(flipud(Slow_Ex_Sel),:)',psth_s(flipud(Slow_Ex_Sel),:)','-y','LineWidth',3);
plot(time_plot(flipud(Fast_In_Sel),:)',psth_s(flipud(Fast_In_Sel),:)','-b','LineWidth',3);
plot(time_plot(flipud(Slow_In_Sel),:)',psth_s(flipud(Slow_In_Sel),:)','-g','LineWidth',3);

%hold on
ylim([-2,max(max(psth_s,[],'all'),max(psth_v,[],'all'))+1]);

%ylim(yallrange);
ylim([-5,100]);

 end
  


end
