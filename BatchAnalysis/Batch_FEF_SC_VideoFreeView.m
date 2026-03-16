function Batch_FEF_SC_VideoFreeView(Data,  StartFile, EndFile,ShowFigureFlag,OutputFlag);
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
     %Latency(index) = Data('ResponseLatency');
    Latency(index) = Data('LatencyC');
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


try
      SM_norm(index) = Data('StimulationStrength_norm');
catch
    keyboard
end


     AS_First(index) = Data('AS_First');
     AS_Second(index) = Data('AS_Second');

     PSTH_Opto_Time = Data('PSTH_Time');

 

     %Memory saccade
     if isfield(OutputData,'MemorySaccade')
     
     MemorySaccadeIndex(index) = 1;
     
     Data =OutputData.MemorySaccade(1).DataStamp;

     %PSTH_Tgt_RF{index} = Data('PSTH_TgtOn_RF');
     %PSTH_Sac_RF{index} = Data('PSTH_SacOn_RF');
     

     PSTH_Tgt_RF{index} = Data('PSTH_TgtOn_RF_z');
     PSTH_Sac_RF{index} = Data('PSTH_SacOn_RF_z');

     Mean_FR_MS(index) = mean([nanmean(Data('PSTH_TgtOn_RF'),'all'),nanmean(Data('PSTH_SacOn_RF'),'all')]);

     
try
    RF_h_v(index) = Data('RF_h_v');
    RF_h_m(index) = Data('RF_h_m');
catch
    tmp = Data('RF_h_v');
    RF_h_v(index) =tmp(1);
    tmp = Data('RF_h_m');
    RF_h_m(index)=tmp(1);
end

     %CellType(index) = Data('CellType');%1:visual;2:motor;3:vis-motor
     %Reasign cell_type
     %{
     if RF_h_v(index)==1 & RF_h_m(index)==0

         CellType(index) = 1;%Visual-Only

     elseif RF_h_v(index)==1 & RF_h_m(index)==1
         CellType(index) = 2;% Visual-Motor

      elseif RF_h_v(index)==0 & RF_h_m(index)==1
         CellType(index) = 3;% Motor-Only

      elseif RF_h_v(index)==0 & RF_h_m(index)==0

         CellType(index) = 5;%Null Type 
     else
         CellType(index) = 4; %Other type
     end
     %}
       %CellType(index) = Data('CellType');%1:visual;2:motor;3:vis-motor
     %Reasign cell_type
     %
     if RF_h_v(index)==1 & (RF_h_m(index)~=1)

         CellType(index) = 1;%Visual-Only

     elseif RF_h_v(index)==1 & RF_h_m(index)==1
         CellType(index) = 2;% Visual-Motor

      elseif (RF_h_v(index)~=1) & RF_h_m(index)==1
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

         % PSTH_TgtTime_DS = Data('PSTH_FR_TgtOn_Time');
         % PSTH_SacTime_DS = Data('PSTH_FR_SacOn_Time');

         
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
          






          






      else
          DelaySaccadeIndex(index) = 0;

           Prefer_SaccadeVector(index,:) =NaN*ones(1,3);%PreDurPost
          Prefer_TgtVector(index) =NaN;

          p_SacV_prefer(index,:) = NaN*ones(1,3);

          p_TgtV_prefer(index) = NaN;

          TuningStrength(index)= NaN;

          SacTuningStrength(index,:)=NaN*ones(1,3);

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
ManuelCheck=[];
ManuelOut = ~ismember(CellIndex,ManuelCheck);


Cri1 = OpticalStim_Mean>1 & CellType<5;



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

 CellType(find(ChannelScreened==21,1,'first'))=2;

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

SM_norm=SM_norm(ScreenCell)';




%Prefer_SaccadeVector =  Prefer_SaccadeVector(ScreenCell,:);

%p_SacV_prefer =  p_SacV_prefer(ScreenCell,:);
%SacTuningStrength_All =  SacTuningStrength(ScreenCell,:);





   

%%
%Separate cells into different modulation groups according to laser
%response

%Proportion of SC neurons which is responsive to FEF optical stimulation
%Sig_m = (p_whole_cri50<0.05 | SM_p<0.05 | PMS_p<0.05)&(~isnan(Latency));
%Sig_m = SM_p<0.05 &(~isnan(Latency));
Sig_m = SM_p<0.05 
Prop_sig = sum(Sig_m)/length(Sig_m);
disp(sprintf('Proportion of SC neurons which can be modulated: %1.2f(%d)',Prop_sig,length(Sig_m)));

FastLat = 10;
Latency_Sig = Latency(Sig_m);

SM_norm_sig = SM_norm(Sig_m);

SM_sig =SM(Sig_m);

Fast_Activate = Latency_Sig < FastLat;
Prop_FastA = sum(Fast_Activate)/sum(Sig_m);

%SM_Fast = SM(Sig_m & Latency<FastLat);
Prop_Exc = sum(SM_norm_sig>0)/length(SM_norm_sig);
Prop_Inhi = 1-Prop_Exc;

Fast_Act_Sel = Sig_m & Latency<FastLat;
Fast_Ex_Sel = Fast_Act_Sel & SM>0;
In_Sel = Sig_m & SM<=0;

Slow_Exc_Sel = Sig_m & Latency>FastLat & SM>0;


disp(sprintf('Proportion of fast activation: %1.2f(%d/%d)',Prop_FastA,sum(Fast_Activate),sum(Sig_m)));
disp(sprintf('Proportion of excitation: %1.2f(%d/%d)',Prop_Exc,(sum(SM_norm_sig>0)),length(SM_norm_sig)));
disp(sprintf('Proportion of inhibition: %1.2f(%d/%d)',Prop_Inhi,(sum(SM_norm_sig<=0)),length(SM_norm_sig)));

PSTH_Use = PSTH_norm;

PSTH_FastEx = PSTH_Use(Fast_Ex_Sel,:);
PSTH_In = PSTH_Use(In_Sel,:);

PSTH_FastEx_Mean = nanmean(PSTH_FastEx,1);
PSTH_In_Mean = nanmean(PSTH_In,1);

PSTH_FastEx_Sem = nanstd(PSTH_FastEx,[],1)/sqrt(size(PSTH_FastEx,1));
PSTH_In_Sem = nanstd(PSTH_In,[],1)/sqrt(size(PSTH_In,1));



%figure
%plot(PSTH_Opto_Time,PSTH_FastEx_Mean,'-r');


%% For analysis of memory saccade

%Select neurons 


if AnalysisMS

    MemorySaccade_Sel = MemorySaccadeIndex ==1;

    CellTypeUnique = unique(CellType);

   % CellTypeUnique = CellTypeUnique(~isnan(CellTypeUnique) & CellTypeUnique<=4);
    CellTypeUnique = CellTypeUnique(~isnan(CellTypeUnique));% & CellTypeUnique<=4);

  
    %For Fast Modulated cell

    CellType_Op = CellType(Fast_Act_Sel & ~isnan(CellType));
    %CellTypeUnique = unique(CellType_Op);

    Latency_Op = Latency(Fast_Act_Sel  & ~isnan(CellType));

    ChannelScreened_curr = ChannelScreened(Fast_Act_Sel  & ~isnan(CellType),:);

for i = 1:length(CellTypeUnique)
    sel = CellType_Op == CellTypeUnique(i);
    NumOfTypes(i) = sum(CellType_Op == CellTypeUnique(i));
     %NumOfTypes_Total(i) = sum(CellType== CellTypeUnique(i));
    
   
    

    PSTH_tgt_c_tmp = PSTH_Tgt_RF(CellType_Op==i,:);
    PSTH_sac_c_tmp = PSTH_Sac_RF(CellType_Op==i,:);

    PSTH_Tgt_Mean{i} = nanmean(PSTH_tgt_c_tmp,1);
    PSTH_Sac_Mean{i} = nanmean(PSTH_sac_c_tmp,1);

    PSTH_Tgt_Sem{i,j} = nanstd( PSTH_tgt_c_tmp ,[],1)/sqrt(sum(sel));
    PSTH_Sac_Sem{i,j} = nanstd( PSTH_tgt_c_tmp ,[],1)/sqrt(sum(sel));
%{
    figure
    subplot(1,2,1)
    plot(PSTH_tgt_c_tmp')
    
    subplot(1,2,2)
    plot(PSTH_sac_c_tmp')


    disp(Latency_Op(sel));
    disp(ChannelScreened_curr(sel));
%}
    %Total number of neurons recorded for each type
    sel_total =  CellType(~isnan(CellType))==CellTypeUnique(i);
    Total_Cell_EachType(i) = round(sum(sel_total));
   % disp(sprintf('Total cells for current type: %d',sum(sel_total )));
disp(sprintf('Cell Type %d N= %d(%d)',CellTypeUnique(i),NumOfTypes(i),Total_Cell_EachType(i)));

end


  TypeFreq =  NumOfTypes./Total_Cell_EachType;
disp(sprintf('N=%d in total',sum(Total_Cell_EachType)));


%Statistics:
 Class1_total = Total_Cell_EachType(1,1:3);
 Class1_group =NumOfTypes(1,1:3);
 Class1_others =Class1_total -Class1_group;
 observed = [Class1_group; Class1_others];

 

% Use built-in function if available
[p_ms_class1,Q] = chi2test(observed);



%keyboard
 


  %For slow modulation cells
  

  %  CellType_Op_slow = CellType(Slow_Act_Sel   & ~isnan(CellType));
    %CellTypeUnique = unique(CellType_Op);

 %   Latency_Op_slow = Latency(Slow_Act_Sel   & ~isnan(CellType));

 %   ChannelScreened_curr_slow = ChannelScreened(Slow_Act_Sel   & ~isnan(CellType),:);
%disp('For slow modulation cells:')
%for i = 1:length(CellTypeUnique)
 %   sel = CellType_Op_slow == CellTypeUnique(i);
 %   NumOfTypes_slow(i) = sum(CellType_Op_slow == CellTypeUnique(i));
 %   disp(sprintf('Cell Type %d N = %d',CellTypeUnique(i),NumOfTypes_slow(i)));
%{
    PSTH_tgt_c_tmp = PSTH_Tgt_RF(CellType_Op==i,:);
    PSTH_sac_c_tmp = PSTH_Sac_RF(CellType_Op==i,:);

    PSTH_Tgt_Mean{i} = nanmean(PSTH_tgt_c_tmp,1);
    PSTH_Sac_Mean{i} = nanmean(PSTH_sac_c_tmp,1);

    PSTH_Tgt_Sem{i,j} = nanstd( PSTH_tgt_c_tmp ,[],1)/sqrt(sum(sel));
    PSTH_Sac_Sem{i,j} = nanstd( PSTH_tgt_c_tmp ,[],1)/sqrt(sum(sel));
%}
%{
    figure
    subplot(1,2,1)
    plot(PSTH_tgt_c_tmp')
    
    subplot(1,2,2)
    plot(PSTH_sac_c_tmp')


    disp(Latency_Op(sel));
    disp(ChannelScreened_curr(sel));
%}
    %Total number of neurons recorded for each type
  %  sel_total =  CellType(~isnan(CellType))==CellTypeUnique(i);
  %  Total_Cell_EachType_slow = sum(sel_total);
   % disp(sprintf('Total cells for current type: %d',sum(sel_total )));


%end


%  TypeFreq_slow =  NumOfTypes_slow./Total_Cell_EachType;

 
end 

 
%end 


%keyboard %Temporary stop here

%% Delay saccacde preference analysis

%figure
if AnalysisDS

    DelaySaccade_Sel = DelaySaccadeIndex ==1;

   % ScreenCellNo_DS=ScreenCellNo(DelaySaccade_Sel );

     for i = 1:length(UniqueMol)

         sel = ModuType == UniqueMol(i) & DelaySaccade_Sel;

         %Significant visual preference 

         vis_sig = p_TgtV_prefer < 0.05;

         Prefer_Vis{i} = Prefer_TgtVector(sel & vis_sig);

         Prefer_VisStrength{i} = TuningStrength(sel & vis_sig);

         ContraPreference(i) =  sum(Prefer_Vis{i} < 180 & Prefer_Vis{i} >= 0);
         IpsiPreference(i) =  sum(~(Prefer_Vis{i} < 180 & Prefer_Vis{i} >= 0));



         VisTuning_Mean_Each{i} = VisTuning_Mean(sel,:);
         VisTuning_Sem_Each{i} = VisTuning_Sem(sel,:);
         VisPre_All{i} = Prefer_TgtVector(sel);
         VisPre_p_All{i} = p_TgtV_prefer (sel);

          VisTuningStrength{i} = TuningStrength(sel);


          ScreenCellNo_DS_Group{i}=ScreenCellNo(sel);

         



       



         presac_sig = p_SacV_prefer(:,1) < 0.05;

         Prefer_preSac{i} = Prefer_SaccadeVector(sel & presac_sig,1);

         Prefer_preSacStrength{i} =SacTuningStrength(sel & vis_sig,1);


         ContraPreference_preSac(i) =  sum(Prefer_preSac{i} < 180 & Prefer_preSac{i} >= 0);
         IpsiPreference_preSac(i) =  sum(~(Prefer_preSac{i}  < 180 & Prefer_preSac{i}  >= 0));


         SacTuning_Mean_Each{i} = SacTuning_Mean(sel,:);
         SacTuning_Sem_Each{i} = SacTuning_Sem(sel,:);
         SacPre_All{i} = Prefer_SaccadeVector(sel,1);
         SacPre_p_All{i} = p_SacV_prefer(sel,1);


         SacTuningStrength_Sac{i} = SacTuningStrength_All(sel,1);






         dursac_sig = p_SacV_prefer(:,2) < 0.05;

         Prefer_durSac{i} = Prefer_SaccadeVector(sel & dursac_sig,2);

         Prefer_durSacStrength{i} =SacTuningStrength(sel & vis_sig,2);


         ContraPreference_durSac(i) =  sum(Prefer_durSac{i} < 180 & Prefer_durSac{i} >= 0);
         IpsiPreference_durSac(i) =  sum(~(Prefer_durSac{i}  < 180 & Prefer_durSac{i}  >= 0));


         postsac_sig = p_SacV_prefer(:,3) < 0.05;

         Prefer_postSac{i} = Prefer_SaccadeVector(sel & postsac_sig,3);

         Prefer_postSacStrength{i} =SacTuningStrength(sel & vis_sig,3);


         ContraPreference_postSac(i) =  sum(Prefer_postSac{i} < 180 & Prefer_postSac{i} >= 0);
         IpsiPreference_postSac(i) =  sum(~(Prefer_postSac{i}  < 180 & Prefer_postSac{i}  >= 0));


         

        NumberOfNeuronsInGroup(i) =sum(sel & (vis_sig|presac_sig));





         %{

         subplot(1,4,i);
         edges = 0:45:360;

         histogram(Prefer_postSac{i},edges);

         hold on
         plot([180,180],[0,10],'--k');
          plot([90,90],[0,10],'--k');
           plot([270,270],[0,10],'--k');
         %}




     end

     disp('DelaySaccade:')

NumberOfNeuronsInGroup




end %End of Analysis DS






%% Plot groups of opto-genetic modulation

if ShowFigureFlag

    %Color order for different groups of neurons

    ColorOrder = [1,0,0;
                  0.9098    0.7020    0.1373;
                  0,0,1;
                  0.3569    0.3451    0.3451];



    FigureStartNum=100;
    FigureIndex=1;

    if PlotOptoEffect

%% Figure 1: Distribution of Modultion Type
    FigureStartNum=FigureStartNum+1;
    FigureIndex=FigureIndex+1;

     figtitlestr{FigureIndex}='PSTH_OptoModulation_Type';
    fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[50,100, 800,800],'Name',figtitlestr{FigureIndex});
    set(fig{FigureIndex}, 'PaperUnits', 'inches');
%formatFig(fig{FigureIndex}, [1.5 1.5], 'nature');
    formatFig(fig{FigureIndex}, [2 2], 'nature');


    %Edge= -1:0.1:1;
    Edge= min(SM)-5:5:max(SM)+5;
    %SM_norm
    %SM
    %{
    histogram(SM_norm,Edge,'FaceColor',[0.5,0.5,0.5]);
    hold on
    histogram(SM_norm_sig(SM_norm_sig>0),Edge,'FaceColor','r');
    histogram(SM_norm_sig(SM_norm_sig<0),Edge,'FaceColor','b');
    xlabel('Distribution of Opto-modulation Index');
    ylabel('Number of neurons');
    legend({'NoEffect','Exc','Inh'});
    box off;

    set(gca,'LineWidth',0.5,'FontSize',5);
    %}
    histogram(SM,Edge,'FaceColor',[0.5,0.5,0.5]);
    hold on
    histogram(SM_sig(SM_sig>0),Edge,'FaceColor','r');
    histogram(SM_sig(SM_sig<0),Edge,'FaceColor','b');
    xlabel('Distribution of Opto-modulation Index');
    ylabel('Number of neurons');
    legend({'NoEffect','Exc','Inh'});
    box off;

    set(gca,'LineWidth',0.5,'FontSize',5);
    
%% %% Figure 2: Distribution of Latency
 FigureStartNum=FigureStartNum+1;
    FigureIndex=FigureIndex+1;

     figtitlestr{FigureIndex}='PSTH_OptoModulation_Latency';
    fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[50,100, 800,800],'Name',figtitlestr{FigureIndex});
    set(fig{FigureIndex}, 'PaperUnits', 'inches');
%formatFig(fig{FigureIndex}, [1.5 1.5], 'nature');
    formatFig(fig{FigureIndex}, [2 2], 'nature');

subplot(2,1,1)
    Edge= [0:5:100];
    %SM_norm
    histogram(Latency_Sig(SM_norm_sig>0& Latency_Sig>10),Edge,'FaceColor','y');
    hold on
    histogram(Latency_Sig(SM_norm_sig>0 & Latency_Sig<10),Edge,'FaceColor','r');
   % histogram(Latency_Sig(SM_norm_sig<0),Edge,'FaceColor','b');
  %  xlabel('Distribution of Opto-modulation Latency (ms)');
    ylabel('Number of neurons');
    legend({'SlowExc','FastExc'});
    box off;
ylim([0,40]);
    set(gca,'LineWidth',0.5,'FontSize',5);

    subplot(2,1,2)
       Edge= [0:5:100];
    %SM_nor
    histogram(Latency_Sig(SM_norm_sig<0),Edge,'FaceColor','b');
   ylim([0,40]);
    xlabel('Distribution of Opto-modulation Latency (ms)');
    ylabel('Number of neurons');
    legend({'Inh'});
    box off;

    set(gca,'LineWidth',0.5,'FontSize',5);

%% %% Figure 3: Porportion of Exc & Inhinition
 FigureStartNum=FigureStartNum+1;
    FigureIndex=FigureIndex+1;



    ModuTotal = length(SM_norm_sig);
    TotalNum = length(SM_norm);
    ExcTotal = sum(SM_norm_sig>0);
    InhTotal = sum(SM_norm_sig<0);
    Prop_Modu = ModuTotal/TotalNum;
    Prop_Exc = ExcTotal/ModuTotal;
    Prop_Inh = InhTotal/ModuTotal;

     figtitlestr{FigureIndex}='PSTH_OptoModulation_Proportion';
    fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[50,100, 300,400],'Name',figtitlestr{FigureIndex});
    set(fig{FigureIndex}, 'PaperUnits', 'inches');
%formatFig(fig{FigureIndex}, [1.5 1.5], 'nature');
    formatFig(fig{FigureIndex}, [0.5 0.6], 'nature');


    h = bar([1,2,3],[Prop_Modu,Prop_Exc,Prop_Inh]);
    h.FaceColor = 'flat';
    h.CData = [0.5,0.5,0.5;1,0,0;0,0,1];
    xticks([1,2,3]);
    xticklabels({'Mod','Exc','Inh'});
    ylabel('Proportion');
    set(gca,'LineWidth',0.5,'FontSize',5);
    box off




%Latency_SigLatency_Sig

%{
%% Figure 2

    figtitlestr{FigureIndex}='PSTH_OptoModulation_Type';
    fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[50,100, 800,800],'Name',figtitlestr{FigureIndex});

    

    for i = 1:length(UniqueMol)
        %subplot(1,length(UniqueMol),i);
        subplot(2,2,i);
        AverageRegion = [0,100];
        area(AverageRegion,[1,1],...
        'FaceColor','k','FaceAlpha',.1,'EdgeAlpha',.1,'HandleVisibility','off');
        hold on

        AverageRegion = [0,100];
        area(AverageRegion,[-1,-1],...
        'FaceColor','k','FaceAlpha',.1,'EdgeAlpha',.1,'HandleVisibility','off');

        shadedErrorBar(PSTH_Opto_Time,PSTH_Stim_Mean_Group(i,:),PSTH_Stim_Sem_Group(i,:),'lineprops',{ColorOrder(i,:)},'transparent',1,'patchSaturation',0.3)
        set(findall(gca, 'Type', 'Line'),'LineWidth',3);
        
        legend(sprintf('N=%d', ModuTypeCount(i)));
        
        title(ModulString(i));
        

        set(gca,'LineWidth',3,'FontSize',20,'FontWeight','Bold');
        xlabel('Time from Opto onset(ms)');
        ylabel('Normalized Responses');

        box off;
        

    end


    FigureStartNum = FigureStartNum+1;
    FigureIndex = FigureIndex+1;

    figtitlestr{FigureIndex}='Latency_Dist';
    fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[50,100, 600,400],'Name',figtitlestr{FigureIndex});

    %Latency
    subplot(1,2,1);

    b = bar(Latency_Mean_Group(1:3),'FaceColor','flat');

    b.CData = ColorOrder(1:3,:);

    xticks(UniqueMol(1:3));
    xticklabels({'Fast Ex','Slow Ex','In'});
    
    ylabel('Latency(ms)');


    hold on
    errorbar(UniqueMol(1:3),Latency_Mean_Group(1:3),Latency_Sem_Group(1:3),'ok','LineWidth',3);
    
    box off
    set(gca,'LineWidth',3,'FontSize',20,'FontWeight','Bold');

    %Strength
    subplot(1,2,2);

    b = bar(MS_Mean_Group(1:3),'FaceColor','flat');

    b.CData = ColorOrder(1:3,:);

    xticks(UniqueMol(1:3));
    xticklabels({'Fast Ex','Slow Ex','In'});
    
    ylabel('Strength');


    hold on
    errorbar(UniqueMol(1:3),MS_Mean_Group(1:3),MS_Sem_Group(1:3),'ok','LineWidth',3);
    
    box off
    set(gca,'LineWidth',3,'FontSize',20,'FontWeight','Bold');
%}
%{

%% NonSigGroup
 FigureStartNum = FigureStartNum+1;
    FigureIndex = FigureIndex+1;

  figtitlestr{FigureIndex}='PSTH_NonSigSep';
   fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[50,100, 1200,400],'Name',figtitlestr{FigureIndex});

   ColorOrder_NonSig=[1,0,0;0,0,1;0.5,0,0;0,0,0.5];

    %PSTH_NonSigStim_Mean_Group

    subplot(1,2,1);

     AverageRegion = [0,100];
        area(AverageRegion,[1,1],...
        'FaceColor','k','FaceAlpha',.1,'EdgeAlpha',.1,'HandleVisibility','off');
        hold on

        AverageRegion = [0,100];
        area(AverageRegion,[-1,-1],...
        'FaceColor','k','FaceAlpha',.1,'EdgeAlpha',.1,'HandleVisibility','off');

        shadedErrorBar(PSTH_Opto_Time,PSTH_NonSigStim_Mean_Group(1,:),PSTH_NonSigStim_Sem_Group(1,:),'lineprops',{ColorOrder_NonSig(1,:)},'transparent',1,'patchSaturation',0.3);
        hold on
        shadedErrorBar(PSTH_Opto_Time,PSTH_NonSigStim_Mean_Group(2,:),PSTH_NonSigStim_Sem_Group(2,:),'lineprops',{ColorOrder_NonSig(2,:)},'transparent',1,'patchSaturation',0.3);

        set(findall(gca, 'Type', 'Line'),'LineWidth',3);

        
        %legend(sprintf('N=%d', NonSigCount(i)));
        legend({sprintf('Ex;N=%d',NonSigCount(1)),sprintf('In;N=%d',NonSigCount(2))});
        
       % title(ModulString(i));
        

        set(gca,'LineWidth',3,'FontSize',20,'FontWeight','Bold');
        xlabel('Time from Opto onset(ms)');
        ylabel('Normalized Responses');

        box off;


    subplot(1,2,2);

       AverageRegion = [0,100];
        area(AverageRegion,[1,1],...
        'FaceColor','k','FaceAlpha',.1,'EdgeAlpha',.1,'HandleVisibility','off');
        hold on

        AverageRegion = [0,100];
        area(AverageRegion,[-1,-1],...
        'FaceColor','k','FaceAlpha',.1,'EdgeAlpha',.1,'HandleVisibility','off');

        shadedErrorBar(PSTH_Opto_Time,PSTH_NonSigStim_Mean_Group(3,:),PSTH_NonSigStim_Sem_Group(3,:),'lineprops',{ColorOrder_NonSig(3,:)},'transparent',1,'patchSaturation',0.3);
        hold on
        shadedErrorBar(PSTH_Opto_Time,PSTH_NonSigStim_Mean_Group(4,:),PSTH_NonSigStim_Sem_Group(4,:),'lineprops',{ColorOrder_NonSig(4,:)},'transparent',1,'patchSaturation',0.3);

        set(findall(gca, 'Type', 'Line'),'LineWidth',3);

        
        %legend(sprintf('N=%d', NonSigCount(i)));
        legend({sprintf('Ex;N=%d',NonSigCount(3)),sprintf('In;N=%d',NonSigCount(4))});
        
       % title(ModulString(i));
        

        set(gca,'LineWidth',3,'FontSize',20,'FontWeight','Bold');
        xlabel('Time from Opto onset(ms)');
        ylabel('Normalized Responses');

        box off;
%}
%{
        %% Latency Distribution
         FigureStartNum = FigureStartNum+1;
    FigureIndex = FigureIndex+1;

        figtitlestr{FigureIndex}='Latency_Distribution';
   fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[50,100, 600,400],'Name',figtitlestr{FigureIndex});
   Edges =[0:5:200];

   subplot(2,1,1);
   histogram(Latency_Ex,Edges,'Facecolor','r','LineWidth',3);
   set(gca, 'YScale', 'log');
   xlabel('Latency');
   ylabel('Number of neurons');
   set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
   box off



   subplot(2,1,2);
   histogram(Latency_In,Edges,'Facecolor','b','LineWidth',3);
   set(gca, 'YScale', 'log');
   xlabel('Latency');
   ylabel('Number of neurons');
   set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
    box off
%}
%{
%% Bar plot 
    FigureStartNum = FigureStartNum+1;
    FigureIndex = FigureIndex+1;

        figtitlestr{FigureIndex}='Modulation';
   fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[50,100, 600,400],'Name',figtitlestr{FigureIndex});

   piechart([ExcitationNum,InhibitionNum,NonSigNum],["Excitation","Inhibition","No-Effect"]);
%{

    ExcitationNum = sum(Excitation & Sig_m);
InhibitionNum = sum(Inhibition & Sig_m);
NonSigNum = sum(~Sig_m);

%}
%}
%{
%% Modulation Strength
%{
MS_Mean_Group(i) = mean(abs(SM(sel)));
    MS_Sem_Group(i) = std(abs(SM(sel)))/sqrt(sum(sel));
%}

FigureStartNum = FigureStartNum+1;
FigureIndex = FigureIndex+1;

 figtitlestr{FigureIndex}='ModulationStrength';
   fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[50,100, 600,400],'Name',figtitlestr{FigureIndex});
        
   b = bar([1,2,3],MS_Mean_Group(1:3));



%}  
   
    end %End of if plot opto effect

%% Memory saccade
  if  PlotMemorySaccade 

%%For memory saccade plot

%Bar plot
FigureStartNum = FigureStartNum+1;
    FigureIndex = FigureIndex+1;

     MS_TypeOrder = [0.8824    0.7451    0.4157;
                        0.2510    0.6902    0.6510;
                        0.6471    0.5020    0.9020];

    figtitlestr{FigureIndex}='PSTH_MemorySac_Type';
    fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[50,100, 500,600],'Name',figtitlestr{FigureIndex});
 set(fig{FigureIndex}, 'PaperUnits', 'inches');
%formatFig(fig{FigureIndex}, [1.5 1.5], 'nature');
    formatFig(fig{FigureIndex}, [2 2], 'nature');


    % b = bar(NumOfTypes(1:3),'FaceColor','r');
    b = bar(TypeFreq(1:3),'FaceColor','r');
    

   %  b.FaceColor = 'Flat';
   %  b.CData = MS_TypeOrder;

    

     box off
     xticks([1,2,3]);
%     xticks(UniqueMol(1:4));
    xticklabels({'Vis','Vis-Mov','Mov'});
    ylabel('Proportion of neurons');
   % legend({'Vis','Vis-Mov','Mov'});

     set(gca,'LineWidth',0.5,'FontSize',5);%,'FontWeight','Bold');
%keyboard

 %Population psth
%{
 for i = 1: length(UniqueMol)
    FigureStartNum = FigureStartNum+1;
    FigureIndex = FigureIndex+1;

    figtitlestr{FigureIndex}= sprintf('PSTH_MemorySac_Type:%d',UniqueMol(i));
    fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[50,100, 1200,600],'Name',figtitlestr{FigureIndex});

    for j = 1:length(CellTypeUnique)
        subplot(length(CellTypeUnique),2,2*(j-1)+1);

        %Target on 


        shadedErrorBar(PSTH_Time_TgtOn ,PSTH_Tgt_Mean{i,j},PSTH_Tgt_Sem{i,j},'lineprops',{MS_TypeOrder(j,:)},'transparent',1,'patchSaturation',0.3)
        set(findall(gca, 'Type', 'Line'),'LineWidth',3);
   %      ylim([-0.5,3]);

        subplot(length(CellTypeUnique),2,2*(j-1)+2);

        %Saccade on 

        shadedErrorBar(PSTH_Time_SacOn ,PSTH_Sac_Mean{i,j},PSTH_Sac_Sem{i,j},'lineprops',{MS_TypeOrder(j,:)},'transparent',1,'patchSaturation',0.3)
        set(findall(gca, 'Type', 'Line'),'LineWidth',3);


      %  ylim([-0.5,3]);


%PSTH_Tgt_Mean{i,j}


    end

 end
%}

     %{
%Bar plot for fast activated neuron only
 FigureStartNum = FigureStartNum+1;
    FigureIndex = FigureIndex+1;

     figtitlestr{FigureIndex}='PSTH_MemorySac_ProjectionNeuron';
    fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[50,100, 500,600],'Name',figtitlestr{FigureIndex});
     b = bar(TypeFreq(1,:));

     
%         b.FaceColor = 'Flat';
        % b.CData = MS_TypeOrder(1,:);
     

     box off
     xticks(UniqueMol(1:3));
    xticklabels({'Vis-Only','Vis-Motor','Motor'});
    ylabel('Proportion');
  %  legend({'Vis-Only','Vis-Motor','Motor'});

     set(gca,'LineWidth',3,'FontSize',20,'FontWeight','Bold');


 


     %}


 end %End of the if plot memory saccade

 %% 
 if PlotDelaySaccade

     FigureStartNum = FigureStartNum+1;
    FigureIndex = FigureIndex+1;


    figtitlestr{FigureIndex}='DelaySaccadeVectorDist';
    fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[50,100, 500,600],'Name',figtitlestr{FigureIndex});
   
    TitleStr = {'ExFast','ExSlow','In','NoEffect'};

    BinSize = 30;

     for i = 1: length(UniqueMol)
         subplot(2,length(UniqueMol),i);

      %  h =  polarplot([zeros(size(Prefer_Vis{i})),Prefer_Vis{i}]',[zeros(size(Prefer_Vis{i})),Prefer_VisStrength{i}]');
      edges=deg2rad([0:BinSize:360]);

     %  polarhistogram(deg2rad(Prefer_Vis{i}),edges);

     h=polarhistogram(deg2rad(Prefer_Vis{i}),edges,'Normalization','probability');
     h.FaceColor = ColorOrder(i,:);

     rlim([0,0.37]);


  %   hold on
%     yrange = ylim;
 %    plot([180,180],[yrange(1) yrange(2)],'--k','LineWidth',3);
     

        title(sprintf('Vis,%s',TitleStr{i}));
        set(gca,'ThetaDir' , 'counterclockwise');
set(gca,'ThetaZeroLocation','top')
set(gca,'FontSize',1,'FontWeight','Bold','LineWidth',3);




         subplot(2,length(UniqueMol),length(UniqueMol)+i);

      %  h =  polarplot([zeros(size(Prefer_Vis{i})),Prefer_Vis{i}]',[zeros(size(Prefer_Vis{i})),Prefer_VisStrength{i}]');
      edges=deg2rad([0:BinSize:360]);

       h=polarhistogram(deg2rad(Prefer_preSac{i}),edges,'Normalization','probability');
       h.FaceColor = ColorOrder(i,:);
    %   hold on
   %    yrange = ylim;
    %    plot([180,180],[yrange(1) yrange(2)],'--k','LineWidth',3);
 rlim([0,0.25]);
    
 
        title(sprintf('PreSac,%s',TitleStr{i}));
        set(gca,'ThetaDir' , 'counterclockwise');
set(gca,'ThetaZeroLocation','top')
set(gca,'FontSize',1,'FontWeight','Bold','LineWidth',3);

     

       
     end

%%

        FigureStartNum = FigureStartNum+1;
    FigureIndex = FigureIndex+1;


    figtitlestr{FigureIndex}='MotorTuningAll';
    fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[50,100, 600,400],'Name',figtitlestr{FigureIndex});

    TitleStr = {'ExFast','ExSlow','In','NoEffect'};
    SaccadeVector= [45,    90,   135,   180,   225,   270,   315, 0];

     for i = 1: length(UniqueMol)
      % VisTuning_Mean_Each{i}
         

         Group= i;
         TotalNum = size( VisTuning_Mean_Each{Group},1);
         

        subplot(2,4,i);
         %For visual motor cell

    for ii = 0:TotalNum-1

             p = VisPre_p_All{Group}(ii+1);
        psac = SacPre_p_All{Group}(ii+1);

%{
          % polarwitherrorbar([SaccadeVector/180*pi SaccadeVector(1)/180*pi],[VisTuning_Mean_Each{Group}(ii+1,:) VisTuning_Mean_Each{Group}(ii+1,1)],[VisTuning_Sem_Each{Group}(ii+1,:) VisTuning_Sem_Each{Group}(ii+1,1)],...
       % '-b',3);
       polarplot([SaccadeVector/180*pi SaccadeVector(1)/180*pi],[VisTuning_Mean_Each{Group}(ii+1,:) VisTuning_Mean_Each{Group}(ii+1,1)]);
    hold on
    h_p = polarplot([VisPre_All{Group}(ii+1) VisPre_All{Group}(ii+1)]/180*pi,[0 1],'r');

        set(h_p,'LineWidth',1);

%}   


      %  polarwitherrorbar([SaccadeVector/180*pi SaccadeVector(1)/180*pi],[SacTuning_Mean_Each{Group}(ii+1,:) SacTuning_Mean_Each{Group}(ii+1,1)],[SacTuning_Sem_Each{Group}(ii+1,:) SacTuning_Sem_Each{Group}(ii+1,1)],...
      %  '-g',3);
      if p<0.05 & psac<0.05
         h_p = polarplot([SaccadeVector/180*pi SaccadeVector(1)/180*pi],[SacTuning_Mean_Each{Group}(ii+1,:) SacTuning_Mean_Each{Group}(ii+1,1)]);
        hold on
        h_p.Color=[0.8,0.8,0.8];
   
        set(h_p,'LineWidth',1);
        
      end
        

    end
     set(gca,'ThetaDir' , 'counterclockwise');
set(gca,'ThetaZeroLocation','top')

    for ii = 0:TotalNum-1
            p = VisPre_p_All{Group}(ii+1);
        psac = SacPre_p_All{Group}(ii+1);
          if p<0.05 & psac<0.05
           h_p = polarplot([SacPre_All{Group}(ii+1) SacPre_All{Group}(ii+1)]/180*pi,[0 1],'b');
           set(h_p,'LineWidth',2);
           hold on
          end

    end
     set(gca,'ThetaDir' , 'counterclockwise');
set(gca,'ThetaZeroLocation','top')
set(gca,'FontSize',1,'FontWeight','Bold','LineWidth',3);

      %For motor cell
      subplot(2,4,4+i);

    for ii = 0:TotalNum-1

             p = VisPre_p_All{Group}(ii+1);
        psac = SacPre_p_All{Group}(ii+1);


      if p>0.05 & psac<0.05
         h_p = polarplot([SaccadeVector/180*pi SaccadeVector(1)/180*pi],[SacTuning_Mean_Each{Group}(ii+1,:) SacTuning_Mean_Each{Group}(ii+1,1)]);
        hold on
         h_p.Color=[0.8,0.8,0.8];
   
        set(h_p,'LineWidth',1);
      end
        

    end
    set(gca,'ThetaDir' , 'counterclockwise');
set(gca,'ThetaZeroLocation','top')

    for ii = 0:TotalNum-1
            p = VisPre_p_All{Group}(ii+1);
        psac = SacPre_p_All{Group}(ii+1);
          if p>0.05 & psac<0.05
           h_p = polarplot([SacPre_All{Group}(ii+1) SacPre_All{Group}(ii+1)]/180*pi,[0 1],'b');
           set(h_p,'LineWidth',2);
           hold on
          end

    end
    set(gca,'ThetaDir' , 'counterclockwise');
set(gca,'ThetaZeroLocation','top')
set(gca,'FontSize',1,'FontWeight','Bold','LineWidth',3);


     end %End of if UniqueModu


  %%
   FigureStartNum = FigureStartNum+1;
    FigureIndex = FigureIndex+1;


    figtitlestr{FigureIndex}='DelaySaccadeVectorDist_Two groups';
    fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[50,100, 500,600],'Name',figtitlestr{FigureIndex});
   
    TitleStr = {'ExFast','Others'};

    BinSize = 30;

    ExtraFast=Prefer_Vis{1};
    Others = [Prefer_Vis{2}',Prefer_Vis{3}',Prefer_Vis{4}'];

    subplot(2,2,1);
     edges=deg2rad([0:BinSize:360]);
      h=polarhistogram(deg2rad(ExtraFast),edges,'Normalization','probability');
     h.FaceColor = [1,0,0];
      title(sprintf('Vis,%s',TitleStr{1}));
        set(gca,'ThetaDir' , 'counterclockwise');
set(gca,'ThetaZeroLocation','top')
set(gca,'FontSize',1,'FontWeight','Bold','LineWidth',3);
rlim([0,0.3]);


     subplot(2,2,2);
     edges=deg2rad([0:BinSize:360]);
      h=polarhistogram(deg2rad(Others),edges,'Normalization','probability');
     h.FaceColor = [1,0,0];
      title(sprintf('PreSac,%s','Others'));
        set(gca,'ThetaDir' , 'counterclockwise');
set(gca,'ThetaZeroLocation','top')
set(gca,'FontSize',1,'FontWeight','Bold','LineWidth',3);
rlim([0,0.3]);

ExtraFastSac=Prefer_preSac{1};
OthersSac = [Prefer_preSac{2}',Prefer_preSac{3}',Prefer_preSac{4}'];

 subplot(2,2,3);
     edges=deg2rad([0:BinSize:360]);
      h=polarhistogram(deg2rad(ExtraFastSac),edges,'Normalization','probability');
     h.FaceColor = [0,0,1];
      title(sprintf('PreSac,%s',TitleStr{1}));
        set(gca,'ThetaDir' , 'counterclockwise');
set(gca,'ThetaZeroLocation','top')
set(gca,'FontSize',1,'FontWeight','Bold','LineWidth',3);
rlim([0,0.3]);


     subplot(2,2,4);
     edges=deg2rad([0:BinSize:360]);
      h=polarhistogram(deg2rad(OthersSac),edges,'Normalization','probability');
     h.FaceColor = [0,0,1];
      title(sprintf('PreSac,%s','Others'));
        set(gca,'ThetaDir' , 'counterclockwise');
set(gca,'ThetaZeroLocation','top')
set(gca,'FontSize',1,'FontWeight','Bold','LineWidth',3);

rlim([0,0.3]);

     %rlim([0,0.37]);

end %End of if plot delay saccade


end %End of if figureplot

%% Output to files
%Setup output directory
Workingdirectory=pwd;
MarkerFolder='DataAnalysis';
%MarkerFolder='DataHub';
Flag=strfind(Workingdirectory,MarkerFolder);
BasicDirectory=Workingdirectory(1:Flag+length(MarkerFolder));

if OutputFlag == 1
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
OutputFigureName='FEF_SC_Opto_';


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

end %End if the OutputFlag == 1
end %End of the function

function checkGroup(name,data)
        for ii = 0:length(name )-1

        fig_index = floor(ii/16)+1;
        subplot_index = mod(ii,16)+1;
        figure(fig_index);
        subplot(4,4,subplot_index);
    

       plot(data(ii+1,:));
       % index=index+1;
       title(name(ii+1))

        end


 end




