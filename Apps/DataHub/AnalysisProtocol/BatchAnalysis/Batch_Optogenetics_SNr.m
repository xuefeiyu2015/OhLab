function Batch_Optogenetics_SNr(Data,  StartFile, EndFile,ShowFigureFlag,OutputFlag);
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


FilesName_Each=Data.FileName;
ChannelID=Data.ChannelNumber;

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
%     p_AS(index) = Data('Ac_p');
     PMS(index) = Data('PostModulationStrength');
     PMS_p(index) = Data('PostModulation_p');
     SM(index) = Data('StimModulation');
     SM_p(index) = Data('StimModulation_p');
     %PSTH_First{index} = Data('PSTH_FirstHalf');
     %PSTH_Second{index} = Data('PSTH_SecondHalf');
     PSTH_First{index} = Data('PSTH_Stim_Mean_First_z');
     PSTH_Second{index} = Data('PSTH_Stim_Mean_Second_z');

     p_whole_cri50{index} = Data('p_whole_cri50');



     


     AS_First(index) = Data('AS_First');
     AS_Second(index) = Data('AS_Second');

     PSTH_Opto_Time = Data('PSTH_Time');

 

     %Memory saccade
     if isfield(OutputData,'MemorySaccade')
     
     MemorySaccadeIndex(index) = 1;
     
     Data =OutputData.MemorySaccade(1).DataStamp;

     %PSTH_Tgt_RF{index} = Data('PSTH_TgtOn_RF');
     %PSTH_Sac_RF{index} = Data('PSTH_SacOn_RF');

     Mean_FR_MS(index) = mean([nanmean(Data('PSTH_TgtOn_RF'),'all'),nanmean(Data('PSTH_SacOn_RF'),'all')]);
     

     PSTH_Tgt_RF{index} = Data('PSTH_TgtOn_RF_z');
     PSTH_Sac_RF{index} = Data('PSTH_SacOn_RF_z');


     
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
     %}
     %{
     if abs(RF_h_v(index))==1 & (RF_h_m(index)==0)

         CellType(index) = 1;%Visual-Only

     elseif abs(RF_h_v(index))==1 & abs(RF_h_m(index))==1
         CellType(index) = 2;% Visual-Motor

      elseif RF_h_v(index)==0 & abs(RF_h_m(index))==1
         CellType(index) = 3;% Motor-Only

      elseif RF_h_v(index)==0 & RF_h_m(index)==0

         CellType(index) = 5;%Null Type 
     else
         CellType(index) = 4; %Other type
     end
     %}
     PSTH_Time_TgtOn = Data('PSTH_Time_TgtOn');
     PSTH_Time_SacOn = Data('PSTH_Time_SacOn');
     else

         MemorySaccadeIndex(index) = 0;

     PSTH_Tgt_RF{index} = NaN*ones(1,140);
     PSTH_Sac_RF{index} = NaN*ones(1,140);


     

     RF_h_v(index) = NaN;
     RF_h_m(index) = NaN;

     CellType(index) =NaN;%1:visual;2:motor;3:vis-motor

   PSTH_Time_TgtOn = NaN*ones(1,140);
    PSTH_Time_SacOn = NaN*ones(1,140);

     
         Mean_FR_MS(index)=NaN;


     end %End of memory saccade

     %% Delay saccade
      if isfield(OutputData,'DelaySaccadeTuning')

          DelaySaccadeIndex(index) = 1;

          Data =OutputData.DelaySaccadeTuning(1).DataStamp;

         % PSTH_TgtTime_DS = Data('PSTH_FR_TgtOn_Time');
         % PSTH_SacTime_DS = Data('PSTH_FR_SacOn_Time');

         try
          Prefer_SaccadeVector(index,:) = Data('Pref_Vect');%PreDurPost
         catch
             keyboard
         end
         
          Prefer_TgtVector(index) = Data('Pref_TgtVect');

          p_SacV_prefer(index,:) = Data('SacVectTuning_p');

          p_TgtV_prefer(index) = Data('TgtVectTuning_p');

          TuningStrength(index) = Data('Pref_TgtVect_Length');

          SacTuningStrength(index,:) = Data('Pref_Vect_Length');

          VisTuning_Mean(index,:)=Data('TgtVectTuning_mean');%/max(Data('TgtVectTuning_mean'));
          VisTuning_Sem(index,:)=Data('TgtVectTuning_sem');%/max(Data('TgtVectTuning_mean'));


          tmp = Data('SacVectTuning_PreDurPost_mean');

          SacTuning_Mean(index,:) = tmp(:,1)';%/max(tmp(:,1)');

          tmp = Data('SacVectTuning_PreDurPost_sem');

          SacTuning_Sem(index,:) = tmp(:,1)';%/max(tmp(:,1)');
          






          %{'TgtVectTuning_mean','TgtVectTuning_sem'
         % 'SacVectTuning_PreDurPost_mean','SacVectTuning_PreDurPost_sem'
          %}







      else
          DelaySaccadeIndex(index) = 0;

           Prefer_SaccadeVector(index,:) =NaN*ones(1,3);%PreDurPost
          Prefer_TgtVector(index) =NaN;

          p_SacV_prefer(index,:) = NaN*ones(1,3);

          p_TgtV_prefer(index) = NaN;

          TuningStrength(index)= NaN;

          SacTuningStrength(index,:)=NaN*ones(1,3);

          VisTuning_Mean(index,:)=NaN*ones(1,8);
           VisTuning_Sem(index,:)=NaN*ones(1,8);
           SacTuning_Mean(index,:)=NaN*ones(1,8);
           SacTuning_Sem(index,:)=NaN*ones(1,8);

      end %End of delay saccade 


      %% %Stable Value
     if isfield(OutputData,'StableValue')

          StableValueIndex(index) = 1;

          Data =OutputData.StableValue(1).DataStamp;


         p_Tgt_RF_stable(index) = Data('p_Tgt_RF');
         if isempty(Data('p_Tgt_nonRF'))
              p_Tgt_nonRF_stable(index) = NaN;
         else
             
             p_Tgt_nonRF_stable(index) = Data('p_Tgt_nonRF');
         end


      %  PreTgtGood_RF(index) = Data('PreTgtGood_RF');
      %  PreTgtBad_RF(index) = Data('PreTgtBad_RF');
        p_Sac_RF_stable(index) = Data('p_Sac_RF');
        if isempty(Data('p_Sac_nonRF'))
            p_Sac_nonRF_stable(index)=NaN;

        else
            p_Sac_nonRF_stable(index) = Data('p_Sac_nonRF');
        end

        

        PSTH_TgtOn_Prefer_RF_stable{index} = Data('PSTH_TgtOn_Prefer_RF_z');
        PSTH_TgtOn_nonPrefer_RF_stable{index} = Data('PSTH_TgtOn_nonPrefer_RF_z');
        PSTH_TgtOn_Prefer_nonRF_stable{index} = Data('PSTH_TgtOn_Prefer_nonRF_z');
        PSTH_TgtOn_nonPrefer_nonRF_stable{index} = Data('PSTH_TgtOn_nonPrefer_nonRF_z');

        PSTH_SacOn_Prefer_RF_stable{index} = Data('PSTH_SacOn_Prefer_RF_z');
        PSTH_SacOn_nonPrefer_RF_stable{index} = Data('PSTH_SacOn_nonPrefer_RF_z');
        PSTH_SacOn_Prefer_nonRF_stable{index} = Data('PSTH_SacOn_Prefer_nonRF_z');
        PSTH_SacOn_nonPrefer_nonRF_stable{index} = Data('PSTH_SacOn_nonPrefer_nonRF_z');



        PSTH_TgtOn_Good_RF_stable{index} = Data('PSTH_TgtOn_Good_RF_z');
        PSTH_TgtOn_Bad_RF_stable{index} = Data('PSTH_TgtOn_Bad_RF_z');
        PSTH_TgtOn_Good_nonRF_stable{index} = Data('PSTH_TgtOn_Good_nonRF_z');
        PSTH_TgtOn_Bad_nonRF_stable{index} = Data('PSTH_TgtOn_Bad_nonRF_z');

        PSTH_SacOn_Good_RF_stable{index} = Data('PSTH_SacOn_Good_RF_z');
        PSTH_SacOn_Bad_RF_stable{index} = Data('PSTH_SacOn_Bad_RF_z');
        PSTH_SacOn_Good_nonRF_stable{index} = Data('PSTH_SacOn_Good_nonRF_z');
        PSTH_SacOn_Bad_nonRF_stable{index} = Data('PSTH_SacOn_Bad_nonRF_z');



        PSTH_Time_TgtOn_stable = Data('Time_TgtOn');
        PSTH_Time_SacOn_stable = Data('Time_SacOn');

        TgtGood_RF_Stable(index) = Data('TgtGood_RF');
        TgtBad_RF_Stable(index) = Data('TgtBad_RF');

         TgtGood_nonRF_Stable(index) = Data('TgtGood_nonRF');
        TgtBad_nonRF_Stable(index) = Data('TgtBad_nonRF');

        SacGood_RF_Stable(index) = Data('SacGood_RF');
        SacBad_RF_Stable(index) = Data('SacBad_RF');

         SacGood_nonRF_Stable(index) = Data('SacGood_nonRF');
        SacBad_nonRF_Stable(index) = Data('SacBad_nonRF');

        
       AUC_Stable_RF(index) = Data('TgtAUC_RF');
       pAUC_Stable_RF(index) = Data('p_AUC_Tgt_RF');

       if isempty(Data('TgtAUC_nonRF'))
           AUC_Stable_nonRF(index) =NaN;
           pAUC_Stable_nonRF(index)=NaN;

       else
            AUC_Stable_nonRF(index) = Data('TgtAUC_nonRF');
             pAUC_Stable_nonRF(index) = Data('p_AUC_Tgt_nonRF ');
       end

        
       
  

p_Response_interval_RF(index,:)=Data('p_response_interval_RF');
TgtAUC_interval_RF(index,:) = Data('TgtAUC_interval_RF');

if isempty(Data('p_response_interval_nonRF'))
    p_Response_interval_nonRF(index,:) = NaN*ones(1,4);
    TgtAUC_interval_nonRF(index,:) = NaN*ones(1,4);

else
    p_Response_interval_nonRF(index,:) = Data('p_response_interval_nonRF');
    TgtAUC_interval_nonRF(index,:)=Data('TgtAUC_interval_nonRF');

end




   




        Value_Vis_RF(index) = sign(TgtGood_RF_Stable(index)-TgtBad_RF_Stable(index));

        Value_Vis_nonRF(index) = sign(TgtGood_nonRF_Stable(index)-TgtBad_nonRF_Stable(index));

        Value_Sac_RF(index) = sign(SacGood_RF_Stable(index)-SacBad_RF_Stable(index));

        Value_Sac_nonRF(index) = sign(SacGood_nonRF_Stable(index)-SacBad_nonRF_Stable(index));




     else
          StableValueIndex(index) = 0;

 p_Sac_nonRF_stable(index) = NaN;

           p_Tgt_RF_stable(index) = NaN;


      %  PreTgtGood_RF(index) = Data('PreTgtGood_RF');
      %  PreTgtBad_RF(index) = Data('PreTgtBad_RF');
        p_Sac_RF_stable(index) = NaN;

        

        PSTH_TgtOn_Prefer_RF_stable{index} = NaN;
        PSTH_TgtOn_nonPrefer_RF_stable{index} = NaN;
        PSTH_TgtOn_Prefer_nonRF_stable{index} = NaN;
        PSTH_TgtOn_nonPrefer_nonRF_stable{index} = NaN;

        PSTH_SacOn_Prefer_RF_stable{index} = NaN;
        PSTH_SacOn_nonPrefer_RF_stable{index} = NaN;
        PSTH_SacOn_Prefer_nonRF_stable{index} = NaN;
        PSTH_SacOn_nonPrefer_nonRF_stable{index} = NaN;



        PSTH_TgtOn_Good_RF_stable{index} =NaN;
        PSTH_TgtOn_Bad_RF_stable{index} = NaN;
        PSTH_TgtOn_Good_nonRF_stable{index} = NaN;
        PSTH_TgtOn_Bad_nonRF_stable{index} = NaN;

        PSTH_SacOn_Good_RF_stable{index} = NaN;
        PSTH_SacOn_Bad_RF_stable{index} = NaN;
        PSTH_SacOn_Good_nonRF_stable{index} = NaN;
        PSTH_SacOn_Bad_nonRF_stable{index} = NaN;


       

        TgtGood_RF_Stable(index) = NaN;
        TgtBad_RF_Stable(index) = NaN;

         TgtGood_nonRF_Stable(index) = NaN;
        TgtBad_nonRF_Stable(index) = NaN;

        SacGood_RF_Stable(index) = NaN;
        SacBad_RF_Stable(index) = NaN;

         SacGood_nonRF_Stable(index) = NaN;
        SacBad_nonRF_Stable(index) = NaN;

        Value_Vis_RF(index) = NaN;

        Value_Vis_nonRF(index) = NaN;

        Value_Sac_RF(index) = NaN;

        Value_Sac_nonRF(index) = NaN;

         AUC_Stable_RF(index) = NaN;
        AUC_Stable_nonRF(index) =NaN;

        pAUC_Stable_RF(index) = NaN;
        pAUC_Stable_nonRF(index) = NaN;

         p_Response_interval_nonRF(index,:) = NaN*ones(1,4);
         TgtAUC_interval_nonRF(index,:) = NaN*ones(1,4);

p_Tgt_nonRF_stable(index) = NaN;


     end


      %% %Flexible Value
     if isfield(OutputData,'FlexibleValue')

          FlexibleValueIndex(index) = 1;

          Data =OutputData.FlexibleValue(1).DataStamp;


         p_Tgt_RF_flexible(index) = Data('p_Tgt_RF');


      %  PreTgtGood_RF(index) = Data('PreTgtGood_RF');
      %  PreTgtBad_RF(index) = Data('PreTgtBad_RF');
        p_Sac_RF_flexible(index) = Data('p_Sac_RF');

        

        PSTH_TgtOn_Prefer_RF_flexible{index} = Data('PSTH_TgtOn_Prefer_RF_z');
        PSTH_TgtOn_nonPrefer_RF_flexible{index} = Data('PSTH_TgtOn_nonPrefer_RF_z');
        PSTH_TgtOn_Prefer_nonRF_flexible{index} = Data('PSTH_TgtOn_Prefer_nonRF_z');
        PSTH_TgtOn_nonPrefer_nonRF_flexible{index} = Data('PSTH_TgtOn_nonPrefer_nonRF_z');

        PSTH_SacOn_Prefer_RF_flexible{index} = Data('PSTH_SacOn_Prefer_RF_z');
        PSTH_SacOn_nonPrefer_RF_flexible{index} = Data('PSTH_SacOn_nonPrefer_RF_z');
        PSTH_SacOn_Prefer_nonRF_flexible{index} = Data('PSTH_SacOn_Prefer_nonRF_z');
        PSTH_SacOn_nonPrefer_nonRF_flexible{index} = Data('PSTH_SacOn_nonPrefer_nonRF_z');



        PSTH_TgtOn_Good_RF_flexible{index} = Data('PSTH_TgtOn_Good_RF_z');
        PSTH_TgtOn_Bad_RF_flexible{index} = Data('PSTH_TgtOn_Bad_RF_z');
        PSTH_TgtOn_Good_nonRF_flexible{index} = Data('PSTH_TgtOn_Good_nonRF_z');
        PSTH_TgtOn_Bad_nonRF_flexible{index} = Data('PSTH_TgtOn_Bad_nonRF_z');

        PSTH_SacOn_Good_RF_flexible{index} = Data('PSTH_SacOn_Good_RF_z');
        PSTH_SacOn_Bad_RF_flexible{index} = Data('PSTH_SacOn_Bad_RF_z');
        PSTH_SacOn_Good_nonRF_flexible{index} = Data('PSTH_SacOn_Good_nonRF_z');
        PSTH_SacOn_Bad_nonRF_flexible{index} = Data('PSTH_SacOn_Bad_nonRF_z');



        PSTH_Time_TgtOn_flexible = Data('Time_TgtOn');
        PSTH_Time_SacOn_flexible = Data('Time_SacOn');

        TgtGood_RF_flexible(index) = Data('TgtGood_RF');
        TgtBad_RF_flexible(index) = Data('TgtBad_RF');

         TgtGood_nonRF_flexible(index) = Data('TgtGood_nonRF');
        TgtBad_nonRF_flexible(index) = Data('TgtBad_nonRF');

        SacGood_RF_flexible(index) = Data('SacGood_RF');
        SacBad_RF_flexible(index) = Data('SacBad_RF');

         SacGood_nonRF_flexible(index) = Data('SacGood_nonRF');
        SacBad_nonRF_flexible(index) = Data('SacBad_nonRF');

        Value_Vis_RF_flexible(index) = sign(TgtGood_RF_flexible(index)-TgtBad_RF_flexible(index));

        Value_Vis_nonRF_flexible(index) = sign(TgtGood_nonRF_flexible(index)-TgtBad_nonRF_flexible(index));

        Value_Sac_RF_flexible(index) = sign(SacGood_RF_flexible(index)-SacBad_RF_flexible(index));

        Value_Sac_nonRF_flexible(index) = sign(SacGood_nonRF_flexible(index)-SacBad_nonRF_flexible(index));

        





      






     else
          FlexibleValueIndex(index) = 0;



           p_Tgt_RF_flexible(index) = NaN;


      %  PreTgtGood_RF(index) = Data('PreTgtGood_RF');
      %  PreTgtBad_RF(index) = Data('PreTgtBad_RF');
        p_Sac_RF_flexible(index) = NaN;

        

        PSTH_TgtOn_Prefer_RF_flexible{index} = NaN;
        PSTH_TgtOn_nonPrefer_RF_flexible{index} = NaN;
        PSTH_TgtOn_Prefer_nonRF_flexible{index} = NaN;
        PSTH_TgtOn_nonPrefer_nonRF_flexible{index} = NaN;

        PSTH_SacOn_Prefer_RF_flexible{index} = NaN;
        PSTH_SacOn_nonPrefer_RF_flexible{index} = NaN;
        PSTH_SacOn_Prefer_nonRF_flexible{index} = NaN;
        PSTH_SacOn_nonPrefer_nonRF_flexible{index} = NaN;



        PSTH_TgtOn_Good_RF_flexible{index} =NaN;
        PSTH_TgtOn_Bad_RF_flexible{index} = NaN;
        PSTH_TgtOn_Good_nonRF_flexible{index} = NaN;
        PSTH_TgtOn_Bad_nonRF_flexible{index} = NaN;

        PSTH_SacOn_Good_RF_flexible{index} = NaN;
        PSTH_SacOn_Bad_RF_flexible{index} = NaN;
        PSTH_SacOn_Good_nonRF_flexible{index} = NaN;
        PSTH_SacOn_Bad_nonRF_flexible{index} = NaN;


       

        TgtGood_RF_flexible(index) = NaN;
        TgtBad_RF_flexible(index) = NaN;

         TgtGood_nonRF_flexible(index) = NaN;
        TgtBad_nonRF_flexible(index) = NaN;

        SacGood_RF_flexible(index) = NaN;
        SacBad_RF_flexible(index) = NaN;

         SacGood_nonRF_flexible(index) = NaN;
        SacBad_nonRF_flexible(index) = NaN;

        Value_Vis_RF_flexible(index) = NaN;

        Value_Vis_nonRF_flexible(index) = NaN;

        Value_Sac_RF_flexible(index) = NaN;

        Value_Sac_nonRF_flexible(index) = NaN;



     end


%% For one DR task

if isfield(OutputData,'OneDR')

    OneDRIndex(index) = 1;

     Data =OutputData.OneDR.DataStamp;

     Tgt_Good_RF(index) = Data('TgtGoodRF');
     Tgt_Bad_RF(index) = Data('TgtBad_RF');   
     p_Tgt_RF(index) = Data('p_Tgt_RF');






     PreTgtGood_RF(index) = Data('PreTgtGood_RF');
     PreTgtBad_RF(index) = Data('PreTgtBad_RF');
     p_preTgt_RF(index) = Data('p_preTgt_RF');

     p_Sac_RF(index) = Data('p_Sac_RF');

      Sac_Good_RF(index) = Data('SacGood_RF');
     Sac_Bad_RF(index) = Data('SacBad_RF');   
    

      PSTH_Good_TgtOn_RF{index} = Data('PSTH_Good_TgtOn_RF_z');

     PSTH_Bad_TgtOn_RF{index} = Data('PSTH_Bad_TgtOn_RF_z');
     PSTH_Good_TgtOn_nonRF{index} = Data('PSTH_Good_TgtOn_nonRF_z');
     PSTH_Bad_TgtOn_nonRF{index} = Data('PSTH_Bad_TgtOn_nonRF_z');

     PSTH_SacGood_nonRF{index} = Data('PSTH_SacGood_nonRF_z');
     PSTH_SacBad_nonRF{index} = Data('PSTH_SacBad_nonRF_z');
     PSTH_SacGood_RF{index} = Data('PSTH_SacGood_RF_z');
     PSTH_SacBad_RF{index} = Data('PSTH_SacBad_RF_z');


     PSTH_Time_TgtOn_OneDR = Data('Time_TgtOn');
     PSTH_Time_SacOn_OneDR = Data('Time_SacOn');


     TgtGood_nonRF(index) = Data('TgtGood_nonRF');
     TgtBad_nonRF(index) = Data('TgtBad_nonRF');
     
     PreTgtGood_nonRF(index) = Data('PreTgtGood_nonRF');
     PreTgtBad_nonRF(index) = Data('PreTgtBad_nonRF');

     
     
     SacGood_nonRF(index) = Data('SacGood_nonRF');
     SacBad_nonRF(index) = Data('SacBad_nonRF');

     
     TgtAUC_nonRF(index) = Data('TgtAUC_nonRF');

     

     SacAUC_nonRF(index) = Data('SacAUC_nonRF');

     TgtAUC_RF(index) = Data('TgtAUC_RF');
     SacAUC_RF(index) = Data('SacAUC_RF');



     p_AUC_Tgt_RF(index) = Data('p_AUC_Tgt_RF');
     p_AUC_Sac_RF(index) = Data('p_AUC_Sac_RF');



     

     p_Tgt_nonRF(index) = Data('p_Tgt_nonRF');

     p_preTgt_nonRF(index) = Data('p_preTgt_nonRF');
     p_Sac_nonRF(index) = Data('p_Sac_nonRF');
     p_AUC_Tgt_nonRF(index) = Data('p_AUC_Tgt_nonRF');
     p_AUC_Sac_nonRF(index) =Data('p_AUC_Sac_nonRF');



     %Change into prefered psth
     Diff_TgtOn_curr = Tgt_Good_RF(index) - Tgt_Bad_RF(index);
     Diff_preTgtOn_curr = PreTgtGood_RF(index) - PreTgtBad_RF(index);

     if p_Tgt_RF(index) < 0.05
         if Diff_TgtOn_curr > 0
             PSTH_Prefer_TgtOn_RF{index} = PSTH_Good_TgtOn_RF{index};
             PSTH_nonPrefer_TgtOn_RF{index} = PSTH_Bad_TgtOn_RF{index};

             PSTH_Prefer_SacOn_RF{index} = PSTH_SacGood_RF{index};
             PSTH_nonPrefer_SacOn_RF{index} = PSTH_SacBad_RF{index};

             PSTH_Prefer_TgtOn_nonRF{index} = PSTH_Good_TgtOn_nonRF{index};
             PSTH_nonPrefer_TgtOn_nonRF{index} = PSTH_Bad_TgtOn_nonRF{index};

             PSTH_Prefer_SacOn_nonRF{index} = PSTH_SacGood_nonRF{index};
             PSTH_nonPrefer_SacOn_nonRF{index} = PSTH_SacBad_nonRF{index};

             




         else
             PSTH_Prefer_TgtOn_RF{index} = PSTH_Bad_TgtOn_RF{index};
             PSTH_nonPrefer_TgtOn_RF{index} = PSTH_Good_TgtOn_RF{index};

             PSTH_Prefer_SacOn_RF{index} = PSTH_SacBad_RF{index};
             PSTH_nonPrefer_SacOn_RF{index} = PSTH_SacGood_RF{index};

             PSTH_Prefer_TgtOn_nonRF{index} = PSTH_Bad_TgtOn_nonRF{index};
             PSTH_nonPrefer_TgtOn_nonRF{index} = PSTH_Good_TgtOn_nonRF{index};

             PSTH_Prefer_SacOn_nonRF{index} = PSTH_SacBad_nonRF{index};
             PSTH_nonPrefer_SacOn_nonRF{index} = PSTH_SacGood_nonRF{index};



         end
     elseif p_preTgt_RF(index) < 0.05
         if Diff_preTgtOn_curr > 0
             PSTH_Prefer_TgtOn_RF{index} = PSTH_Good_TgtOn_RF{index};
             PSTH_nonPrefer_TgtOn_RF{index} = PSTH_Bad_TgtOn_RF{index};

             PSTH_Prefer_SacOn_RF{index} = PSTH_SacGood_RF{index};
             PSTH_nonPrefer_SacOn_RF{index} = PSTH_SacBad_RF{index};

             PSTH_Prefer_TgtOn_nonRF{index} = PSTH_Good_TgtOn_nonRF{index};
             PSTH_nonPrefer_TgtOn_nonRF{index} = PSTH_Bad_TgtOn_nonRF{index};

             PSTH_Prefer_SacOn_nonRF{index} = PSTH_SacGood_nonRF{index};
             PSTH_nonPrefer_SacOn_nonRF{index} = PSTH_SacBad_nonRF{index};


         else
             PSTH_Prefer_TgtOn_RF{index} = PSTH_Bad_TgtOn_RF{index};
             PSTH_nonPrefer_TgtOn_RF{index} = PSTH_Good_TgtOn_RF{index};

             PSTH_Prefer_SacOn_RF{index} = PSTH_SacBad_RF{index};
             PSTH_nonPrefer_SacOn_RF{index} = PSTH_SacGood_RF{index};

             PSTH_Prefer_TgtOn_nonRF{index} = PSTH_Bad_TgtOn_nonRF{index};
             PSTH_nonPrefer_TgtOn_nonRF{index} = PSTH_Good_TgtOn_nonRF{index};

             PSTH_Prefer_SacOn_nonRF{index} = PSTH_SacBad_nonRF{index};
             PSTH_nonPrefer_SacOn_nonRF{index} = PSTH_SacGood_nonRF{index};



         end
     else
         if Diff_preTgtOn_curr > 0
            PSTH_Prefer_TgtOn_RF{index} = PSTH_Good_TgtOn_RF{index};
            PSTH_nonPrefer_TgtOn_RF{index} = PSTH_Bad_TgtOn_RF{index};

            PSTH_Prefer_SacOn_RF{index} = PSTH_SacGood_RF{index};
            PSTH_nonPrefer_SacOn_RF{index} = PSTH_SacBad_RF{index};

            PSTH_Prefer_TgtOn_nonRF{index} = PSTH_Good_TgtOn_nonRF{index};
            PSTH_nonPrefer_TgtOn_nonRF{index} = PSTH_Bad_TgtOn_nonRF{index};

            PSTH_Prefer_SacOn_nonRF{index} = PSTH_SacGood_nonRF{index};
            PSTH_nonPrefer_SacOn_nonRF{index} = PSTH_SacBad_nonRF{index};
         else
             PSTH_Prefer_TgtOn_RF{index} = PSTH_Bad_TgtOn_RF{index};
             PSTH_nonPrefer_TgtOn_RF{index} = PSTH_Good_TgtOn_RF{index};

             PSTH_Prefer_SacOn_RF{index} = PSTH_SacBad_RF{index};
             PSTH_nonPrefer_SacOn_RF{index} = PSTH_SacGood_RF{index};

             PSTH_Prefer_TgtOn_nonRF{index} = PSTH_Bad_TgtOn_nonRF{index};
             PSTH_nonPrefer_TgtOn_nonRF{index} = PSTH_Good_TgtOn_nonRF{index};

             PSTH_Prefer_SacOn_nonRF{index} = PSTH_SacBad_nonRF{index};
             PSTH_nonPrefer_SacOn_nonRF{index} = PSTH_SacGood_nonRF{index};
             
         end
     end

       Value_Vis_RF_OneDR(index) = sign(Tgt_Good_RF(index)-Tgt_Bad_RF(index));

    %    Value_Vis_nonRF_OneDR(index) = sign(Tgt_Good_nonRF(index)-Tgt_Bad_nonRF(index));
    Value_Vis_nonRF_OneDR(index) =NaN;

        Value_Sac_RF_OneDR(index) = sign(Sac_Good_RF(index)-Sac_Bad_RF(index));

     %   Value_Sac_nonRF_OneDR(index) = sign(Sac_Good_nonRF(index)-Sac_Bad_nonRF(index));
     Value_Sac_nonRF_OneDR(index) = NaN;

       

        
else
        OneDRIndex(index) = 0;

        Value_Vis_RF_OneDR(index) = NaN;

        Value_Vis_nonRF_OneDR(index) = NaN;

        Value_Sac_RF_OneDR(index) =NaN;

        Value_Sac_nonRF_OneDR(index) =NaN;




    
end


     
     index = index+1;   
    end %End of for file name
     
end


%Plot controller
PlotOptoEffect = 1;
PlotMemorySaccade = 0;
PlotDelaySaccade = 0;
PlotStableValue = 0;
PlotFlexibleValue= 0;
PlotOneDR = 0;

%Analysis controller 
AnalysisOS = 0;% For Video Free View
AnalysisMS = 0;% For memory saccade
AnalysisDS = 0;%For Delay saccade
AnalysisSV = 0;%For stable  value delay saccade
AnalysisOneDR = 0;%For  one direction reward saccade
AnalysisFV = 0; %For flexible value task

%Selection 


% OpticalStim_Mean(index) = mean(PSTH_Stim_Mean{index},'all');
   %  ControlStim_Mean(index) = mean(PSTH_Control_Mean{index},'all');
   %  FR_Mean_Video(index)

   A=[OpticalStim_Mean',ControlStim_Mean',FR_Mean_Video'];

%ManuelOut
CellIndex = 1:length(FR_Mean_Video);




%ManuelCheck=[550 469	478	477	473	479	425	408	353	351	320	239	240	217	213	201	178	163	160	139	141	68	43	48 115 132 137 157 181 200 221 496 516 529 538 581 623 97 131 159 315 325 329 343 350 400 427 435];
ManuelCheck=[496];
%ManuelCheck=CellIndex(OpticalStim_Mean<2);
ManuelOut = ~ismember(CellIndex,ManuelCheck);


Cri1 = OpticalStim_Mean>2;
%Cri1 = FR_Mean_Video>5;
Exclude = ~Cri1;
ExcludeNumber = sum(Exclude | ~ManuelOut);
disp(ExcludeNumber)



ScreenedIndex = CellIndex(Cri1 & ManuelOut);

ScreenCell = ismember(CellIndex,ScreenedIndex);

ScreenCellNo = CellIndex(ScreenCell);

Prop_Exclude = ExcludeNumber/length(Cri1);
disp(sprintf('Proportion of excluded neurons: %1.2f',Prop_Exclude));


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

PSTH_Tgt_RF = cell2mat(cellfun(@(x) x(1:140), PSTH_Tgt_RF(ScreenCell),'UniformOutput',false)');
PSTH_Sac_RF = cell2mat(PSTH_Sac_RF(ScreenCell)');

PSTH_Time_TgtOn = PSTH_Time_TgtOn(1:140);
%{
p_Cri20(index) = Data('p_cri50');
     p_Cri50(index) = Data('p_Cri50');
     Latency2(index) = Data('Latency2');

     PSTH_Stim_Norm2(index) = Data('psth_stim_normed2');

%}

%Other information
%try
%SM,SM_p,Latency,ChannelNO, CellType,FR_Mean_Video

All_Info = [SM',SM_p',Latency',ChannelFull', CellType',FR_Mean_Video',MemorySaccadeIndex',DelaySaccadeIndex',Prefer_TgtVector',p_TgtV_prefer', TuningStrength',AS',p_Cri20',p_Cri50',Latency2',AbsoluteIndex',p_whole_cri50',PMS_p'];

All_Info_Screened = All_Info(ScreenCell,:);

FilesNameScreened =FilesNameAll(ScreenCell);

SM = All_Info_Screened(:,1);
SM_p= All_Info_Screened(:,2);

Latency = All_Info_Screened(:,3);
ChannelScreened = All_Info_Screened(:,4);
CellType = All_Info_Screened(:,5);

FR_Mean_Video = All_Info_Screened(:,6);

MemorySaccadeIndex = All_Info_Screened(:,7);

DelaySaccadeIndex = All_Info_Screened(:,8);

Prefer_TgtVector = All_Info_Screened(:,9);

p_TgtV_prefer = All_Info_Screened(:,10);

TuningStrength =  All_Info_Screened(:,11);
AS =  All_Info_Screened(:,12);

p_Cri20 = All_Info_Screened(:,13);

p_Cri50 = All_Info_Screened(:,14);

Latency2 =  All_Info_Screened(:,15);

AbsoluteIndex = All_Info_Screened(:,16);

p_whole_cri50= All_Info_Screened(:,17);


PMS_p = All_Info_Screened(:,18);


p_Tgt_nonRF_stable=p_Tgt_nonRF_stable(ScreenCell)';




Prefer_SaccadeVector =  Prefer_SaccadeVector(ScreenCell,:);

p_SacV_prefer =  p_SacV_prefer(ScreenCell,:);
SacTuningStrength_All =  SacTuningStrength(ScreenCell,:);


VisTuning_Mean = VisTuning_Mean(ScreenCell,:);

VisTuning_Sem = VisTuning_Sem(ScreenCell,:);



SacTuning_Mean = SacTuning_Mean(ScreenCell,:);

SacTuning_Sem = SacTuning_Sem(ScreenCell,:);

Mean_FR_MS=Mean_FR_MS(ScreenCell)';

%{
   
ValueTaskInfo = [p_Tgt_RF_stable', p_Sac_RF_stable',TgtGood_RF_Stable',TgtBad_RF_Stable',TgtGood_nonRF_Stable',TgtBad_nonRF_Stable',...
    SacGood_RF_Stable',SacBad_RF_Stable',SacGood_nonRF_Stable',SacBad_nonRF_Stable',Value_Vis_RF',Value_Vis_nonRF',Value_Sac_RF', Value_Sac_nonRF',StableValueIndex'];

ValueTaskInfo_Screened = ValueTaskInfo(ScreenCell,:);

p_Tgt_RF_stable = ValueTaskInfo_Screened(:,1);
p_Sac_RF_stable = ValueTaskInfo_Screened(:,2);
TgtGood_RF_Stable = ValueTaskInfo_Screened(:,3);
TgtBad_RF_Stable = ValueTaskInfo_Screened(:,4);
TgtGood_nonRF_Stable = ValueTaskInfo_Screened(:,5);
TgtBad_nonRF_Stable = ValueTaskInfo_Screened(:,6);
SacGood_RF_Stable = ValueTaskInfo_Screened(:,7);
SacBad_RF_Stable = ValueTaskInfo_Screened(:,8);
SacGood_nonRF_Stable = ValueTaskInfo_Screened(:,9);
SacBad_nonRF_Stable = ValueTaskInfo_Screened(:,10);
Value_Vis_RF = ValueTaskInfo_Screened(:,11);
Value_Vis_nonRF = ValueTaskInfo_Screened(:,12);
Value_Sac_RF = ValueTaskInfo_Screened(:,13);
Value_Sac_nonRF = ValueTaskInfo_Screened(:,14);
StableValueIndex = ValueTaskInfo_Screened(:,15);

p_Sac_nonRF_stable=p_Sac_nonRF_stable(ScreenCell)';



 AUC_Stable_RF =  AUC_Stable_RF(ScreenCell);
 AUC_Stable_nonRF = AUC_Stable_nonRF(ScreenCell);

 pAUC_Stable_RF = pAUC_Stable_RF(ScreenCell);
 pAUC_Stable_nonRF = pAUC_Stable_nonRF(ScreenCell);

 


 PSTH_TgtOn_Prefer_RF_stable=PSTH_TgtOn_Prefer_RF_stable(ScreenCell) ;
 PSTH_TgtOn_nonPrefer_RF_stable = PSTH_TgtOn_nonPrefer_RF_stable(ScreenCell) ;
 PSTH_TgtOn_Prefer_nonRF_stable = PSTH_TgtOn_Prefer_nonRF_stable(ScreenCell);
 PSTH_TgtOn_nonPrefer_nonRF_stable = PSTH_TgtOn_nonPrefer_nonRF_stable(ScreenCell);

 PSTH_SacOn_Prefer_RF_stable = PSTH_SacOn_Prefer_RF_stable(ScreenCell);
 PSTH_SacOn_nonPrefer_RF_stable = PSTH_SacOn_nonPrefer_RF_stable(ScreenCell);
 PSTH_SacOn_Prefer_nonRF_stable= PSTH_SacOn_Prefer_nonRF_stable(ScreenCell);
 PSTH_SacOn_nonPrefer_nonRF_stable = PSTH_SacOn_nonPrefer_nonRF_stable(ScreenCell);

 


PSTH_TgtOn_Good_RF_stable =PSTH_TgtOn_Good_RF_stable(ScreenCell);
PSTH_TgtOn_Bad_RF_stable = PSTH_TgtOn_Bad_RF_stable(ScreenCell) ;
PSTH_TgtOn_Good_nonRF_stable= PSTH_TgtOn_Good_nonRF_stable(ScreenCell);
PSTH_TgtOn_Bad_nonRF_stable= PSTH_TgtOn_Bad_nonRF_stable(ScreenCell);

PSTH_SacOn_Good_RF_stable = PSTH_SacOn_Good_RF_stable(ScreenCell) ;
PSTH_SacOn_Bad_RF_stable = PSTH_SacOn_Bad_RF_stable(ScreenCell);
PSTH_SacOn_Good_nonRF_stable = PSTH_SacOn_Good_nonRF_stable(ScreenCell) ;
PSTH_SacOn_Bad_nonRF_stable = PSTH_SacOn_Bad_nonRF_stable(ScreenCell);





%% Flexible value 
ValueTaskInfoFlexible = [p_Tgt_RF_flexible', p_Sac_RF_flexible',TgtGood_RF_flexible',TgtBad_RF_flexible',TgtGood_nonRF_flexible',TgtBad_nonRF_flexible',...
    SacGood_RF_flexible',SacBad_RF_flexible',SacGood_nonRF_flexible',SacBad_nonRF_flexible',Value_Vis_RF_flexible',Value_Vis_nonRF_flexible',Value_Sac_RF_flexible', Value_Sac_nonRF_flexible',FlexibleValueIndex'];

ValueTaskInfo_Screened = ValueTaskInfo(ScreenCell,:);

p_Tgt_RF_flexible = ValueTaskInfo_Screened(:,1);
p_Sac_RF_flexible = ValueTaskInfo_Screened(:,2);
TgtGood_RF_flexible = ValueTaskInfo_Screened(:,3);
TgtBad_RF_flexible = ValueTaskInfo_Screened(:,4);
TgtGood_nonRF_flexible = ValueTaskInfo_Screened(:,5);
TgtBad_nonRF_flexible = ValueTaskInfo_Screened(:,6);
SacGood_RF_flexible = ValueTaskInfo_Screened(:,7);
SacBad_RF_flexible = ValueTaskInfo_Screened(:,8);
SacGood_nonRF_flexible = ValueTaskInfo_Screened(:,9);
SacBad_nonRF_flexible = ValueTaskInfo_Screened(:,10);
Value_Vis_RF_flexible = ValueTaskInfo_Screened(:,11);
Value_Vis_nonRF_flexible = ValueTaskInfo_Screened(:,12);
Value_Sac_RF_flexible = ValueTaskInfo_Screened(:,13);
Value_Sac_nonRF_flexible = ValueTaskInfo_Screened(:,14);
FlexibleValueIndex = ValueTaskInfo_Screened(:,15);



 PSTH_TgtOn_Prefer_RF_flexible=PSTH_TgtOn_Prefer_RF_flexible(ScreenCell) ;
 PSTH_TgtOn_nonPrefer_RF_flexible = PSTH_TgtOn_nonPrefer_RF_flexible(ScreenCell) ;
 PSTH_TgtOn_Prefer_nonRF_flexible = PSTH_TgtOn_Prefer_nonRF_flexible(ScreenCell);
 PSTH_TgtOn_nonPrefer_nonRF_flexible = PSTH_TgtOn_nonPrefer_nonRF_flexible(ScreenCell);

 PSTH_SacOn_Prefer_RF_flexible = PSTH_SacOn_Prefer_RF_flexible(ScreenCell);
 PSTH_SacOn_nonPrefer_RF_flexible = PSTH_SacOn_nonPrefer_RF_flexible(ScreenCell);
 PSTH_SacOn_Prefer_nonRF_flexible= PSTH_SacOn_Prefer_nonRF_flexible(ScreenCell);
 PSTH_SacOn_nonPrefer_nonRF_flexible = PSTH_SacOn_nonPrefer_nonRF_flexible(ScreenCell);

 


PSTH_TgtOn_Good_RF_flexible =PSTH_TgtOn_Good_RF_flexible(ScreenCell);
PSTH_TgtOn_Bad_RF_flexible = PSTH_TgtOn_Bad_RF_flexible(ScreenCell) ;
PSTH_TgtOn_Good_nonRF_flexible= PSTH_TgtOn_Good_nonRF_flexible(ScreenCell);
  PSTH_TgtOn_Bad_nonRF_flexible= PSTH_TgtOn_Bad_nonRF_flexible(ScreenCell);

        PSTH_SacOn_Good_RF_flexible = PSTH_SacOn_Good_RF_flexible(ScreenCell) ;
        PSTH_SacOn_Bad_RF_flexible = PSTH_SacOn_Bad_RF_flexible(ScreenCell);
        PSTH_SacOn_Good_nonRF_flexible = PSTH_SacOn_Good_nonRF_flexible(ScreenCell) ;
        PSTH_SacOn_Bad_nonRF_flexible = PSTH_SacOn_Bad_nonRF_flexible(ScreenCell);





%% One DR
OneDRTaskInfo = [p_Tgt_RF', p_Sac_RF', p_preTgt_RF',Tgt_Good_RF',Tgt_Bad_RF',PreTgtGood_RF',PreTgtBad_RF',OneDRIndex',Value_Vis_RF_OneDR',Value_Vis_nonRF_OneDR',Value_Sac_RF_OneDR',...
    Value_Sac_nonRF_OneDR'   ];

OneDRTaskScreened = OneDRTaskInfo(ScreenCell,:);

p_Tgt_RF = OneDRTaskScreened(:,1);
p_Sac_RF = OneDRTaskScreened(:,2);
p_preTgt_RF = OneDRTaskScreened(:,3);
Tgt_Good_RF = OneDRTaskScreened(:,4);
Tgt_Bad_RF = OneDRTaskScreened(:,5);
PreTgtGood_RF= OneDRTaskScreened(:,6);
PreTgtBad_RF= OneDRTaskScreened(:,7);

OneDRIndex = OneDRTaskScreened(:,8);

Value_Vis_RF_OneDR = OneDRTaskScreened(:,9);
Value_Vis_nonRF_OneDR = OneDRTaskScreened(:,10);
Value_Sac_RF_OneDR = OneDRTaskScreened(:,11);
Value_Sac_nonRF_OneDR = OneDRTaskScreened(:,12);


PSTH_Prefer_TgtOn_RF = PSTH_Prefer_TgtOn_RF(ScreenCell);
PSTH_nonPrefer_TgtOn_RF = PSTH_nonPrefer_TgtOn_RF(ScreenCell);

PSTH_Prefer_SacOn_RF= PSTH_Prefer_SacOn_RF(ScreenCell);
PSTH_nonPrefer_SacOn_RF= PSTH_nonPrefer_SacOn_RF(ScreenCell);

PSTH_Prefer_TgtOn_nonRF= PSTH_Prefer_TgtOn_nonRF(ScreenCell);
PSTH_nonPrefer_TgtOn_nonRF = PSTH_nonPrefer_TgtOn_nonRF(ScreenCell);

PSTH_Prefer_SacOn_nonRF = PSTH_Prefer_SacOn_nonRF (ScreenCell);
PSTH_nonPrefer_SacOn_nonRF= PSTH_nonPrefer_SacOn_nonRF(ScreenCell);



p_preTgt_nonRF = p_preTgt_nonRF(ScreenCell)';
p_Sac_nonRF = p_Sac_nonRF(ScreenCell)';
p_AUC_Tgt_nonRF = p_AUC_Tgt_nonRF(ScreenCell)';
p_AUC_Sac_nonRF =p_AUC_Sac_nonRF(ScreenCell)';
p_Tgt_nonRF = p_Tgt_nonRF(ScreenCell)';

TgtAUC_RF = TgtAUC_RF(ScreenCell)';
 SacAUC_RF = SacAUC_RF(ScreenCell)';

 TgtAUC_NonRF = TgtAUC_nonRF(ScreenCell)';
 SacAUC_NonRF = SacAUC_nonRF(ScreenCell)';

 p_AUC_Tgt_RF = p_AUC_Tgt_RF(ScreenCell)';
 p_AUC_Sac_RF = p_AUC_Sac_RF(ScreenCell)';




p_Response_interval_RF = p_Response_interval_RF(ScreenCell,:);
TgtAUC_interval_RF = TgtAUC_interval_RF(ScreenCell,:);
 p_Response_interval_nonRF = p_Response_interval_nonRF(ScreenCell,:);
TgtAUC_interval_nonRF = TgtAUC_interval_nonRF(ScreenCell,:);

FilesNameAll=FilesNameAll(ScreenCell)';

%{
catch
    disp('Sth wrong for value tasks');
end
%}

%{
 p_Tgt_RF_stable(index) = Data('p_Tgt_RF');


      %  PreTgtGood_RF(index) = Data('PreTgtGood_RF');
      %  PreTgtBad_RF(index) = Data('PreTgtBad_RF');
        p_Sac_RF_stable(index) = Data('p_Sac_RF');

  TgtGood_RF_Stable(index) = NaN;
        TgtBad_RF_Stable(index) = NaN;

         TgtGood_nonRF_Stable(index) = NaN;
        TgtBad_nonRF_Stable(index) = NaN;

        SacGood_RF_Stable(index) = NaN;
        SacBad_RF_Stable(index) = NaN;

         SacGood_nonRF_Stable(index) = NaN;
        SacBad_nonRF_Stable(index) = NaN;

        Value_Vis_RF(index) = NaN;

        Value_Vis_nonRF(index) = NaN;

        Value_Sac_RF(index) = NaN;

        Value_Sac_nonRF(index) = NaN;


%}

%PSTH_Tgt_RF_mat = cellfun(@(x) [x,NaN*ones(1,maxnum_ms-length(x))],PSTH_Tgt_RF,'UniformOutput',0);
%}
%%
%Separate cells into different modulation groups according to laser
%response


%Separate into 5 groups 


Fast_M = Latency<10;
Slow_M = Latency>=10;
Excitation = SM> 0;
Inhibition = SM<=0;

%Sig_m = SM_p < 0.05;
%Sig_m = ~isnan(Latency) & AS >=5;

%Sig_m = (p_whole_cri50<0.05 | SM_p<0.05 | PMS_p<0.05)&(~isnan(Latency));
Sig_m = SM_p<0.05&(~isnan(Latency));

%Sig_m = p_whole_cri50>-0.1;
%Latency
Non_Sig_m = ~Sig_m;

SepLevel = 0.5;
%{
Non_Sig_level2 = p_whole_cri50>SepLevel & SM_p>SepLevel;
Non_Sig_level1 = ~Non_Sig_level2  & Non_Sig_m;
%}
Non_Sig_level2 = SM_p>SepLevel;
Non_Sig_level1 = ~Non_Sig_level2  & Non_Sig_m;

%Non_Sig_m = Non_Sig_level2;


Ex_Fast_Sig = Fast_M & Excitation & Sig_m;

Ex_Slow_Sig = Slow_M & Excitation & Sig_m;


%In_Fast_Sig = Fast_M & Inhibition & Sig_m;

%In_Slow_Sig = Slow_M & Inhibition & Sig_m;

In_Sig = Inhibition & Sig_m;

%PSTH_Use = PSTH_norm;
%PSTH_Use = PSTH; 
%PSTH_Use = PSTH_z; 
%PSTH_Use = PSTH_norm2;

PSTH_Use = PSTH_norm;

ModuType = NaN*ones(1,length(Sig_m));
NSig_Level = NaN*ones(1,length(Sig_m));



%%
ModuType(Ex_Fast_Sig) =1;
ModuType(Ex_Slow_Sig) =2;

ModuType(In_Sig) =3;
%ModuType(In_Slow_Sig) =4;

ModuType(Non_Sig_m) =4;

NSig_Level(Non_Sig_level1&Excitation) =1;
NSig_Level(Non_Sig_level1&Inhibition) =2;

NSig_Level(Non_Sig_level2&Excitation) =3;
NSig_Level(Non_Sig_level2&Inhibition) =4;

NSig_Level=NSig_Level';


ModulString = {'Fast Excitation','Slow Excitation','Inhibition','No Effect'};

NonSigString = {'NonSig_Ex_l1','NonSig_In_l1','NonSig_Ex_l2','NonSig_In_l2'};

Latency_Ex = Latency(Ex_Fast_Sig|Ex_Slow_Sig);
Latency_In = Latency(In_Sig);





UniqueMol = unique(ModuType);

UniqueMol = UniqueMol(~isnan(UniqueMol));

ModuType=ModuType';%Transfer to column;

UniqueNonSig = unique(NSig_Level);
UniqueNonSig=UniqueNonSig(~isnan(UniqueNonSig));

%For bar plot
ExcitationNum = sum(Excitation & Sig_m);
InhibitionNum = sum(Inhibition & Sig_m);
NonSigNum = sum(~Sig_m);



%
%Find examples
%{
example_fast_saccade=FilesNameAll(Ex_Fast_Sig & CellType==3);
Channel_sel = ChannelScreened(Ex_Fast_Sig & CellType==3)';
A=Channel_sel ;
keyboard
%}
%%
%{
%Check = isnan(Latency) & Sig_m;

%A=PSTH_Use(Check,:);
A=PSTH_Use;
figure
%index = 1;
for index = 0:size(A,1)-1

    
    
    %if ModuType(i+1) == 2
    fig_index = floor(index+1/16)+1;
    subplot_index = mod(index+1,16)+1;
    figure(fig_index);
    subplot(4,4,subplot_index);

    

        plot(A(index+1,:));
      %  index=index+1;
   % end
   
end
%}
%{
figure
subplot(2,1,1)
edge = 0:2.5:150;
histogram(Latency(Excitation & Sig_m),edge );
subplot(2,1,2)
histogram(Latency(Inhibition & Sig_m),edge )
%}


%% Manuel Check for All cells 
%For each group
%{
Group1_name = AbsoluteIndex(ModuType == 1);

Group1_PSTH_Stim = PSTH_Use(ModuType == 1,:);

Group2_name = AbsoluteIndex(ModuType == 2);

Group2_PSTH_Stim = PSTH_Use(ModuType == 2,:);

Group3_name = AbsoluteIndex(ModuType == 3);

Group3_PSTH_Stim = PSTH_Use(ModuType == 3,:);

Group4_name = AbsoluteIndex(ModuType == 4);

Group4_PSTH_Stim = PSTH_Use(ModuType == 4,:);
%}

%index = 1;

%{
for i = 0:length(Group1_name )-1

    fig_index = floor(i/16)+1;
    subplot_index = mod(i,16)+1;
    figure(fig_index);
    subplot(4,4,subplot_index);
    

        plot(Group1_PSTH_Stim(i+1,:));
       % index=index+1;
       title(Group1_name(i+1))

end
%}
%{
for i = 0:length(Group2_name )-1

    fig_index = floor(i/16)+1;
    subplot_index = mod(i,16)+1;
    figure(fig_index);
    subplot(4,4,subplot_index);
    

        plot(Group2_PSTH_Stim(i+1,:));
       % index=index+1;
       title(Group2_name(i+1))

end

%}

%checkGroup(Group4_name,Group4_PSTH_Stim);

  
%AbsoluteIndex

%keyboard


for i = 1:length(UniqueMol)
    sel = ModuType == UniqueMol(i);

    ModuTypeCount(i) = sum(sel);

    PSTH_Stim_Group = PSTH_Use(sel,:);

    names = AbsoluteIndex(sel);

  

    PSTH_Stim_Mean_Group(i,:) = mean(PSTH_Stim_Group,1);

    PSTH_Stim_Sem_Group(i,:) = std(PSTH_Stim_Group,[],1)/sqrt(sum(sel));
    

    Latency_Mean_Group(i) = mean(Latency(sel));

    Latency_Sem_Group(i) = std(Latency(sel))/sqrt(sum(sel));

    MS_Mean_Group(i) = mean(abs(SM(sel)));
    MS_Sem_Group(i) = std(abs(SM(sel)))/sqrt(sum(sel));


end
ModuTypeCount

%% If ananlyze optostimulation effect 
if AnalysisOS

%Anova for latency test for all groups 
%[p,tbl,stats] = anova1(hogg);
[p_Anova_Latency,ANOVATAB,stats]= anova1(Latency(ModuType<4),ModuType(ModuType<4),'off');

disp(sprintf('Anova test for latency for groups:%d',p_Anova_Latency));

%Perform multiple comparison
results = multcompare(stats);

disp('Mutiple comparison between different groups for latency: ');
tbl = array2table(results,"VariableNames", ...
    ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"]);

%Anova for modulation strength
[p_Anova_Strength,ANOVATAB,stats]= anova1(abs(SM(ModuType<4)),ModuType(ModuType<4),'off');

disp(sprintf('Anova test for strength for groups:%d',p_Anova_Strength));

%Perform multiple comparison
results = multcompare(stats);

disp('Mutiple comparison between different groups for strength: ');
tbl = array2table(results,"VariableNames", ...
    ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"]);

close all;
end 
%}

%% For nonsig neurons 
for i = 1:length(UniqueNonSig)
    sel = NSig_Level == UniqueNonSig(i);
    NonSigCount(i) = sum(sel);


    PSTH_NonSigStim_Group = PSTH_Use(sel,:);

   

  

    PSTH_NonSigStim_Mean_Group(i,:) = mean(PSTH_NonSigStim_Group,1);

    PSTH_NonSigStim_Sem_Group(i,:) = std(PSTH_NonSigStim_Group,[],1)/sqrt(sum(sel));

   



end



%{
%Onebyone check
for i = 0:size(PSTH_Use,1)-1

    fig_index = floor(i/16)+1;
    subplot_index = mod(i,16)+1;


    figure(fig_index);

    subplot(4,4,subplot_index);
    if ModuType(i+1) == 4
    plot(PSTH_Use(i+1,:));
    title(sprintf('%d,Type:%d',ScreenCellNo(i+1),ModuType(i+1)));
    end

end
%}


%% For analysis of memory saccade

%Select neurons 

if AnalysisMS

    MemorySaccade_Sel = MemorySaccadeIndex ==1;

    CellTypeUnique = unique(CellType);
    CellTypeUnique = CellTypeUnique(~isnan(CellTypeUnique));

    %CellTypeUnique = CellTypeUnique(~isnan(CellTypeUnique) & CellTypeUnique<=3);

   % CellTypeUnique = ones(size(CellTypeUnique));%For a test

  % A=[Mean_FR_MS',MemorySaccadeIndex,CellType];

  Select_FR = Mean_FR_MS > 1;
  disp(sprintf('Memory Saccade exclude: %1.2f',sum(Select_FR==0)/length(Select_FR)));


 
  

  % UniqueMol = 1;%For a test;
 % ModuType = ones(size(ModuType));%For a test
    for i = 1:length(UniqueMol)
        for j = 1:length(CellTypeUnique)
            if UniqueMol(i) == 1 &  (CellTypeUnique(j)==3)
            sel = (ModuType == 1| ModuType == 2) & CellType == CellTypeUnique(j) & MemorySaccade_Sel & Select_FR;
          
            else
            sel = ModuType == UniqueMol(i) & CellType == CellTypeUnique(j) & MemorySaccade_Sel & Select_FR;
            end


            TypeCounts(i,j) = sum(sel);

            PSTH_Tgt_Mean{i,j} = nanmean(PSTH_Tgt_RF(sel,:),1);
            PSTH_Sac_Mean{i,j} = nanmean(PSTH_Sac_RF(sel,:),1);

            PSTH_Tgt_Sem{i,j} = nanstd(PSTH_Tgt_RF(sel,:),[],1)/sqrt(sum(sel));
            PSTH_Sac_Sem{i,j} = nanstd(PSTH_Sac_RF(sel,:),[],1)/sqrt(sum(sel));

           % Mean_FR_MS

           PSTH_Tgt_Each{i,j} = PSTH_Tgt_RF(sel,:);
            PSTH_Sac_Each{i,j} = PSTH_Sac_RF(sel,:);

           
        

        end



            

    end


%devider = repmat(sum( TypeCounts,2),1,size(TypeCounts,2));

devider = repmat(sum( TypeCounts,1),size(TypeCounts,1),1);
    
 TypeFreq =  TypeCounts./devider;


 
end 


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

%Manule Check for DS for each group 

%
%{
Group = 1;
TotalNum = size( VisTuning_Mean_Each{Group},1);
SaccadeVector= [45,    90,   135,   180,   225,   270,   315, 0];

for ii = 0:TotalNum-1
%
        fig_index = floor(ii/16)+1;
        subplot_index = mod(ii,16)+1;
        figure(fig_index);
        subplot(4,4,subplot_index);

        polarwitherrorbar([SaccadeVector/180*pi SaccadeVector(1)/180*pi],[VisTuning_Mean_Each{Group}(ii+1,:) VisTuning_Mean_Each{Group}(ii+1,1)],[VisTuning_Sem_Each{Group}(ii+1,:) VisTuning_Sem_Each{Group}(ii+1,1)],...
        '-b',3);
    hold on
    h_p = polarplot([VisPre_All{Group}(ii+1) VisPre_All{Group}(ii+1)]/180*pi,[0 VisTuningStrength{Group}(ii+1)],'r');

        set(h_p,'LineWidth',2.5);


        polarwitherrorbar([SaccadeVector/180*pi SaccadeVector(1)/180*pi],[SacTuning_Mean_Each{Group}(ii+1,:) SacTuning_Mean_Each{Group}(ii+1,1)],[SacTuning_Sem_Each{Group}(ii+1,:) SacTuning_Sem_Each{Group}(ii+1,1)],...
        '-g',3);
    hold on
    h_p = polarplot([SacPre_All{Group}(ii+1) SacPre_All{Group}(ii+1)]/180*pi,[0 SacTuningStrength_Sac{Group}(ii+1)],'m');

        set(h_p,'LineWidth',2.5);

        title(ScreenCellNo_DS_Group{Group}(ii+1));
        

%}
    %{
    p = VisPre_p_All{Group}(ii+1);
    psac = SacPre_p_All{Group}(ii+1);
subplot(1,2,1);
if p <0.05
   plotcolor = 'r';
   if p<0.05 && psac>0.05
       %Visual only tuning
       plotcolor = 'm';
   elseif p<0.05 && psac < 0.05
        plotcolor ='r';


   end
 h_p = polarplot([VisPre_All{Group}(ii+1) VisPre_All{Group}(ii+1)]/180*pi,[0 1],plotcolor);
 hold on
 
 set(h_p,'LineWidth',1.5);
    
set(gca,'ThetaDir' , 'counterclockwise');
set(gca,'ThetaZeroLocation','top')
set(gca,'FontSize',1,'FontWeight','Bold','LineWidth',3);
end

subplot(1,2,2);

 
 if psac < 0.05
    if psac<0.05 && p>0.05
       %Sac only tuning
       plotcolor = 'y';
   elseif p<0.05 && psac < 0.05
        plotcolor = 'g';


   end

  h_p = polarplot([SacPre_All{Group}(ii+1) SacPre_All{Group}(ii+1)]/180*pi,[0 1],plotcolor);
 hold on
    
set(gca,'ThetaDir' , 'counterclockwise');
set(gca,'ThetaZeroLocation','top')
set(gca,'FontSize',1,'FontWeight','Bold','LineWidth',3);
    
 set(h_p,'LineWidth',1.5);
    
end 

end %End of manule check
%}



%{
  VisTuning_Mean_Each{i} = VisTuning_Mean(sel,:);
         VisTuning_Sem_Each{i} = VisTuning_Sem(sel,:);
         VisPre_All{i} = Prefer_TgtVector(sel);
         VisPre_p_All{i} = p_TgtV_prefer (sel);




         SacTuning_Mean_Each{i} = SacTuning_Mean(sel,:);
         SacTuning_Sem_Each{i} = SacTuning_Sem(sel,:);
         SacPre_All{i} = Prefer_SaccadeVector(sel,1);
         SacPre_p_All{i} = p_SacV_prefer(sel,1);



%}


end %End of Analysis DS



%% Stable Value Task
if AnalysisSV

    StableValue_Sel = StableValueIndex ==1;

    for i = 1:length(UniqueMol)

         sel = ModuType == UniqueMol(i) & StableValue_Sel;

        % SigVal = (p_Tgt_RF_stable<0.05 | p_Sac_RF_stable<0.05) & sel ;
      %  SigVal = p_Tgt_RF_stable<0.05  & sel ;
%{
        p_Response_interval_RF(index,:)=Data('p_response_interval_RF');
TgtAUC_interval_RF(index,:) = Data('TgtAUC_interval_RF');

if isempty(Data('p_response_interval_nonRF'))
    p_Response_interval_nonRF(index,:) = NaN*ones(1,4);
    TgtAUC_interval_nonRF(index,:) = NaN*ones(1,4);

else
    p_Response_interval_nonRF(index,:) = Data('p_response_interval_nonRF');
    TgtAUC_interval_nonRF(index,:)=Data('TgtAUC_interval_nonRF');

end
%}
        %

       % p_Response_interval_RF_tmp = p_Response_interval_RF(sel,:);
        VisMat=p_Response_interval_RF<0.05;

        p_Vis = sum(VisMat,2)>0;

        VisMat_nonRF=p_Response_interval_nonRF<0.05;
        p_Vis_nonRF = sum(VisMat_nonRF,2)>0;

      %  SigVis = (p_Vis |p_Vis_nonRF) & sel ;

     %  SigVal = (p_Tgt_RF_stable<0.05)  & SigVis;
      SigVal = (p_Tgt_RF_stable<0.05|p_Sac_RF_stable<0.05)& sel ;
      NonSigVal = ~(p_Tgt_RF_stable<0.05|p_Sac_RF_stable<0.05)& sel ;

    %  SigVal = SigVis; %For visualization, will change back later

       % SigVal_nonRF = p_Tgt_nonRF_stable<0.05  & SigVis;
       SigVal_nonRF = (p_Tgt_nonRF_stable<0.05| p_Sac_nonRF_stable<0.05)& sel ;
       NonSigVal_nonRF = ~(p_Tgt_nonRF_stable<0.05| p_Sac_nonRF_stable<0.05)& sel ;


       SigValAll = SigVal|SigVal_nonRF;
       NonSigValAll = NonSigVal|NonSigVal_nonRF;
       

      
       


       % SigVis =SigVis&SigVal;


      %  VisMat=p_Response_interval_RF<0.05


        %  AUC_Stable_RF_Type{i}=AUC_Stable_RF( SigVis);
      %  AUC_Stable_nonRF_Type{i}=AUC_Stable_nonRF( SigVis);

      AUC_Stable_RF_Type{i}=AUC_Stable_RF( SigVal);
       AUC_Stable_nonRF_Type{i}=AUC_Stable_nonRF( SigVal_nonRF);

       GoodPrefer_RF = (AUC_Stable_RF>0.5)' & SigVal;
       GoodPrefer_nonRF = (AUC_Stable_nonRF>0.5)' & SigVal;

       BadPrefer_RF = (AUC_Stable_RF<0.5)' & SigVal;
       BadPrefer_nonRF = (AUC_Stable_nonRF<0.5)' & SigVal;


       

       

       GoodPrefer = GoodPrefer_RF|GoodPrefer_nonRF;
       BadPrefer = BadPrefer_RF|BadPrefer_nonRF;

       Value_RFOnly = SigVal&NonSigVal_nonRF;
       Value_nonRFOnly = NonSigVal & SigVal_nonRF;
       Value_Both = SigVal & SigVal_nonRF;

       Consistent = (GoodPrefer_RF & GoodPrefer_nonRF) | (BadPrefer_RF & BadPrefer_nonRF);
       InConsistent = (GoodPrefer_RF & BadPrefer_nonRF) | (BadPrefer_RF & GoodPrefer_nonRF);

       ConsistentCount(i)=sum(Consistent);
       InConsistentCount(i)=sum(InConsistent);


      
   %    Value_Vis_RF=AUC_Stable_RF>0.5;


         %Cell Type Counts

         ValueCellCount(i) = sum(SigValAll);
         NonValueCellCount(i) = sum(NonSigValAll);

      %   PropValueCell(i) =ValueCellCount(i) /sum(SigVis);
      PropValueCell(i) =ValueCellCount(i) /sum(sel);


        GoodPref_Vis_RF = Value_Vis_RF == 1 & SigVal;
         GoodPref_Vis_RF_Count(i) = sum(GoodPref_Vis_RF);
         GoodPref_Vis_RF_Prop(i) = sum(GoodPref_Vis_RF)/ValueCellCount(i);

        GoodPrefer_Count(i)=sum(GoodPrefer);
        BadPrefer_Count(i)=sum(BadPrefer);

        GoodPrefer_Prop(i) =GoodPrefer_Count(i)/(GoodPrefer_Count(i)+BadPrefer_Count(i));

        Value_RFOnlyCount(i) = sum(Value_RFOnly);
        Value_nonRFOnlyCount(i) = sum(Value_nonRFOnly);
        Value_BothCount(i) = sum(Value_Both);

        

        ValueInfo_Summary(i,:)=[ValueCellCount(i),NonValueCellCount(i),Value_RFOnlyCount(i),Value_nonRFOnlyCount(i),Value_BothCount(i),GoodPrefer_Count(i),BadPrefer_Count(i),ConsistentCount(i),InConsistentCount(i)];


        %GoodPref_Vis_RF = Value_Vis_RF == 1 & SigVal;

        %{
        GoodPref_RF_Count(i) = sum(GoodPrefer_RF);
        BadPref_RF_Count(i) = sum(BadPrefer_RF);
        GoodPref_RF_Prop(i) = GoodPref_RF_Count(i)/(GoodPref_RF_Count(i)+BadPref_RF_Count(i));



         
         GoodPref_nonRF_Count(i) = sum(GoodPrefer_nonRF);
         BadPref_nonRF_Count(i) = sum(BadPrefer_nonRF);
         GoodPref_nonRF_Prop(i) = sum(GoodPrefer_nonRF)/(sum(GoodPrefer_nonRF)+sum(BadPrefer_nonRF));
        %}

         

      


       


        %%PSTH for one by one check
%{
         p_Tgt_RF_stable_Mo{i}= p_Tgt_RF_stable(SigVis);
         p_Sac_RF_stable_Mo{i}= p_Sac_RF_stable(SigVis);

         p_Tgt_nonRF_stable_Mo{i} = p_Tgt_nonRF_stable(SigVis);

%}
       
        p_Tgt_RF_stable_Mo{i}= p_Tgt_RF_stable(SigVal);
         p_Sac_RF_stable_Mo{i}= p_Sac_RF_stable(SigVal);

         p_Tgt_nonRF_stable_Mo{i} = p_Tgt_nonRF_stable(SigVal_nonRF);

        %Value_Vis_RF
%{
            PSTH_TgtOn_Good_RF_stable_Mo{i} = cell2mat(PSTH_TgtOn_Good_RF_stable(SigVis)');
            PSTH_TgtOn_Bad_RF_stable_Mo{i} = cell2mat(PSTH_TgtOn_Bad_RF_stable(SigVis)');

            PSTH_TgtOn_Good_nonRF_stable_Mo{i} = cell2mat(PSTH_TgtOn_Good_nonRF_stable(SigVis)');
            PSTH_TgtOn_Bad_nonRF_stable_Mo{i} = cell2mat(PSTH_TgtOn_Bad_nonRF_stable(SigVis)');

             PSTH_SacOn_Good_RF_stable_Mo{i} = cell2mat(PSTH_SacOn_Good_RF_stable(SigVis)');
            PSTH_SacOn_Bad_RF_stable_Mo{i} = cell2mat(PSTH_SacOn_Bad_RF_stable(SigVis)');

            PSTH_SacOn_Good_nonRF_stable_Mo{i} = cell2mat(PSTH_SacOn_Good_nonRF_stable(SigVis)');
            PSTH_SacOn_Bad_nonRF_stable_Mo{i} = cell2mat(PSTH_SacOn_Bad_nonRF_stable(SigVis)');

%}

        PSTH_TgtOn_Good_RF_stable_Mo{i} = cell2mat(PSTH_TgtOn_Good_RF_stable(SigVal)');
            PSTH_TgtOn_Bad_RF_stable_Mo{i} = cell2mat(PSTH_TgtOn_Bad_RF_stable(SigVal)');

            PSTH_TgtOn_Good_nonRF_stable_Mo{i} = cell2mat(PSTH_TgtOn_Good_nonRF_stable(SigVal)');
            PSTH_TgtOn_Bad_nonRF_stable_Mo{i} = cell2mat(PSTH_TgtOn_Bad_nonRF_stable(SigVal)');

             PSTH_SacOn_Good_RF_stable_Mo{i} = cell2mat(PSTH_SacOn_Good_RF_stable(SigVal)');
            PSTH_SacOn_Bad_RF_stable_Mo{i} = cell2mat(PSTH_SacOn_Bad_RF_stable(SigVal)');

            PSTH_SacOn_Good_nonRF_stable_Mo{i} = cell2mat(PSTH_SacOn_Good_nonRF_stable(SigVal)');
            PSTH_SacOn_Bad_nonRF_stable_Mo{i} = cell2mat(PSTH_SacOn_Bad_nonRF_stable(SigVal)');   

          

         %PSTH
       

         PSTH_TgtOn_Prefer_RF_stable_sel = cell2mat(PSTH_TgtOn_Prefer_RF_stable(SigVal)');

         PSTH_TgtOn_Prefer_RF_stable_mean(i,:) = mean(PSTH_TgtOn_Prefer_RF_stable_sel,1);

         PSTH_TgtOn_Prefer_RF_stable_sem(i,:) = std(PSTH_TgtOn_Prefer_RF_stable_sel,1)/sqrt(sum(SigVal));


          PSTH_TgtOn_NonPrefer_RF_stable_sel = cell2mat(PSTH_TgtOn_nonPrefer_RF_stable(SigVal)');

         PSTH_TgtOn_NonPrefer_RF_stable_mean(i,:) = mean(PSTH_TgtOn_NonPrefer_RF_stable_sel ,1);

         PSTH_TgtOn_NonPrefer_RF_stable_sem(i,:) = std(PSTH_TgtOn_NonPrefer_RF_stable_sel ,1)/sqrt(sum(SigVal));


         PSTH_SacOn_Prefer_RF_stable_sel = cell2mat(PSTH_SacOn_Prefer_RF_stable(SigVal)');

         PSTH_SacOn_Prefer_RF_stable_mean(i,:) = mean(PSTH_SacOn_Prefer_RF_stable_sel,1);

         PSTH_SacOn_Prefer_RF_stable_sem(i,:) = std(PSTH_SacOn_Prefer_RF_stable_sel,1)/sqrt(sum(SigVal));


          PSTH_SacOn_NonPrefer_RF_stable_sel = cell2mat(PSTH_SacOn_nonPrefer_RF_stable(SigVal)');

         PSTH_SacOn_NonPrefer_RF_stable_mean(i,:) = mean(PSTH_SacOn_NonPrefer_RF_stable_sel ,1);

         PSTH_SacOn_NonPrefer_RF_stable_sem(i,:) = std(PSTH_SacOn_NonPrefer_RF_stable_sel ,1)/sqrt(sum(SigVal));


        %Non RF


          PSTH_TgtOn_Prefer_nonRF_stable_sel = cell2mat(PSTH_TgtOn_Prefer_nonRF_stable(SigVal_nonRF)');

         PSTH_TgtOn_Prefer_nonRF_stable_mean(i,:) = nanmean(PSTH_TgtOn_Prefer_nonRF_stable_sel,1);

         PSTH_TgtOn_Prefer_nonRF_stable_sem(i,:) = nanstd(PSTH_TgtOn_Prefer_nonRF_stable_sel,[],1)/sqrt(sum(SigVal_nonRF));


          PSTH_TgtOn_NonPrefer_nonRF_stable_sel = cell2mat(PSTH_TgtOn_nonPrefer_nonRF_stable(SigVal_nonRF)');

         PSTH_TgtOn_NonPrefer_nonRF_stable_mean(i,:) = nanmean(PSTH_TgtOn_NonPrefer_nonRF_stable_sel ,1);

         PSTH_TgtOn_NonPrefer_nonRF_stable_sem(i,:) = nanstd(PSTH_TgtOn_NonPrefer_nonRF_stable_sel ,[],1)/sqrt(sum(SigVal_nonRF));


         PSTH_SacOn_Prefer_nonRF_stable_sel = cell2mat(PSTH_SacOn_Prefer_nonRF_stable(SigVal_nonRF)');

         PSTH_SacOn_Prefer_nonRF_stable_mean(i,:) = nanmean(PSTH_SacOn_Prefer_nonRF_stable_sel,1);

         PSTH_SacOn_Prefer_nonRF_stable_sem(i,:) = nanstd(PSTH_SacOn_Prefer_nonRF_stable_sel,[],1)/sqrt(sum(SigVal_nonRF));


          PSTH_SacOn_NonPrefer_nonRF_stable_sel = cell2mat(PSTH_SacOn_nonPrefer_nonRF_stable(SigVal_nonRF)');

         PSTH_SacOn_NonPrefer_nonRF_stable_mean(i,:) = nanmean(PSTH_SacOn_NonPrefer_nonRF_stable_sel ,1);

         PSTH_SacOn_NonPrefer_nonRF_stable_sem(i,:) = nanstd(PSTH_SacOn_NonPrefer_nonRF_stable_sel ,[],1)/sqrt(sum(SigVal_nonRF));
        

%{
          PSTH_TgtOn_Prefer_RF_stable{index} =PSTH_TgtOn_Prefer_RF_stable{ScreenCell} ;
 PSTH_TgtOn_nonPrefer_RF_stable{index} = PSTH_TgtOn_nonPrefer_RF_stable{ScreenCell} ;
 PSTH_TgtOn_Prefer_nonRF_stable{index} = PSTH_TgtOn_Prefer_nonRF_stable{ScreenCell};
 PSTH_TgtOn_nonPrefer_nonRF_stable{index} = PSTH_TgtOn_nonPrefer_nonRF_stable{ScreenCell} ;

 PSTH_SacOn_Prefer_RF_stable{index} = PSTH_SacOn_Prefer_RF_stable{ScreenCell};
 PSTH_SacOn_nonPrefer_RF_stable{index} = PSTH_SacOn_nonPrefer_RF_stable{ScreenCell};
 PSTH_SacOn_Prefer_nonRF_stable{index} = PSTH_SacOn_Prefer_nonRF_stable{ScreenCell} ;
 PSTH_SacOn_nonPrefer_nonRF_stable{index} = PSTH_SacOn_nonPrefer_nonRF_stable{ScreenCell};



PSTH_TgtOn_Good_RF_stable{index} =PSTH_TgtOn_Good_RF_stable{ScreenCell};
PSTH_TgtOn_Bad_RF_stable{index} = PSTH_TgtOn_Bad_RF_stable{ScreenCell} ;
        PSTH_TgtOn_Good_nonRF_stable{index} = PSTH_TgtOn_Good_nonRF_stable{ScreenCell};
        PSTH_TgtOn_Bad_nonRF_stable{index} = PSTH_TgtOn_Bad_nonRF_stable{ScreenCell};

        PSTH_SacOn_Good_RF_stable{index} = PSTH_SacOn_Good_RF_stable{ScreenCell} ;
        PSTH_SacOn_Bad_RF_stable{index} = PSTH_SacOn_Bad_RF_stable{ScreenCell};
        PSTH_SacOn_Good_nonRF_stable{index} = PSTH_SacOn_Good_nonRF_stable{ScreenCell} ;
        PSTH_SacOn_Bad_nonRF_stable{index} = PSTH_SacOn_Bad_nonRF_stable{ScreenCell} ;

%}

    end


    
%% One by One check

%%Check the ROC
%{
 for i = 1:length(UniqueMol)
    % AUC_Stable_RF_Type{i}=AUC_Stable_RF(sel);
     %   AUC_Stable_nonRF_Type{i}=AUC_Stable_nonRF(sel);
    p_RF =  p_Tgt_RF_stable_Mo{i};
    p_nonRF = p_Tgt_nonRF_stable_Mo{i};

    auc_RF = AUC_Stable_RF_Type{i}';
    auc_nonRF = AUC_Stable_nonRF_Type{i}';

    sel_good_prefer = auc_RF>0.5 & p_RF<0.05;
     sel_bad_prefer = auc_RF<0.5& p_RF<0.05;

   %RF

   psth_good = PSTH_TgtOn_Good_RF_stable_Mo{i};
   psth_bad = PSTH_TgtOn_Bad_RF_stable_Mo{i};


   psth_total = [psth_good,psth_bad];
   [COEFF, SCORE, LATENT, TSQUARED, EXPLAINED, MU] = pca(psth_total');

   PC3 = COEFF(:,1:3);
   psth_good_trans = SCORE(1:size(psth_good,1),1:3);
   psth_bad_trans = SCORE(size(psth_good,1)+1:end,1:3);

   figure
%
   plot3(psth_good_trans(:,1),psth_good_trans(:,2),psth_good_trans(:,3),'-r');
   hold on
   plot3(psth_bad_trans(:,1),psth_bad_trans(:,2),psth_bad_trans(:,3),'-b');
%}
   

%}

  % keyboard

   %{

    psth_goodprefer_good_mean = mean(PSTH_TgtOn_Good_RF_stable_Mo{i}(sel_good_prefer,:),1);
    psth_goodprefer_bad_mean = mean(PSTH_TgtOn_Bad_RF_stable_Mo{i}(sel_good_prefer,:),1);

    psth_goodprefer_good_sem = std(PSTH_TgtOn_Good_RF_stable_Mo{i}(sel_good_prefer,:),1)/sqrt(sum(sel));
    psth_goodprefer_bad_sem = std(PSTH_TgtOn_Bad_RF_stable_Mo{i}(sel_good_prefer,:),1)/sqrt(sum(sel));


     psth_badprefer_good_mean = mean(PSTH_TgtOn_Good_RF_stable_Mo{i}(sel_bad_prefer,:),1);
    psth_badprefer_bad_mean = mean(PSTH_TgtOn_Bad_RF_stable_Mo{i}(sel_bad_prefer,:),1);

    psth_badprefer_good_sem = std(PSTH_TgtOn_Good_RF_stable_Mo{i}(sel_bad_prefer,:),1)/sqrt(sum(sel));
    psth_badprefer_bad_sem = std(PSTH_TgtOn_Bad_RF_stable_Mo{i}(sel_bad_prefer,:),1)/sqrt(sum(sel));


%Non-RF

     psth_goodprefer_good_nonRF_mean = nanmean(PSTH_TgtOn_Good_nonRF_stable_Mo{i}(sel_good_prefer,:),1);
    psth_goodprefer_bad_nonRF_mean = nanmean(PSTH_TgtOn_Bad_nonRF_stable_Mo{i}(sel_good_prefer,:),1);

    psth_goodprefer_good_nonRF_sem = nanstd(PSTH_TgtOn_Good_nonRF_stable_Mo{i}(sel_good_prefer,:),[],1)/sqrt(sum(sel));
    psth_goodprefer_bad_nonRF_sem = nanstd(PSTH_TgtOn_Bad_RF_stable_Mo{i}(sel_good_prefer,:),[],1)/sqrt(sum(sel));


     psth_badprefer_good_nonRF_mean = nanmean(PSTH_TgtOn_Good_nonRF_stable_Mo{i}(sel_bad_prefer,:),1);
    psth_badprefer_bad_nonRF_mean = nanmean(PSTH_TgtOn_Bad_nonRF_stable_Mo{i}(sel_bad_prefer,:),1);

    psth_badprefer_good_nonRF_sem = nanstd(PSTH_TgtOn_Good_nonRF_stable_Mo{i}(sel_bad_prefer,:),[],1)/sqrt(sum(sel));
    psth_badprefer_bad_nonRF_sem = nanstd(PSTH_TgtOn_Bad_nonRF_stable_Mo{i}(sel_bad_prefer,:),[],1)/sqrt(sum(sel));





     subplot(4,4,(i-1)*4+1);

      shadedErrorBar(PSTH_Time_TgtOn_stable ,psth_goodprefer_good_mean,psth_goodprefer_good_sem,'lineprops',{[1,0,0]},'transparent',1,'patchSaturation',0.3);
      hold on
       shadedErrorBar(PSTH_Time_TgtOn_stable ,psth_goodprefer_bad_mean,psth_goodprefer_bad_sem,'lineprops',{[0,0,1]},'transparent',1,'patchSaturation',0.3);

        set(findall(gca, 'Type', 'Line'),'LineWidth',3);

         subplot(4,4,(i-1)*4+2);

      shadedErrorBar(PSTH_Time_TgtOn_stable ,psth_badprefer_good_mean,psth_goodprefer_good_sem,'lineprops',{[1,0,0]},'transparent',1,'patchSaturation',0.3);
      hold on
       shadedErrorBar(PSTH_Time_TgtOn_stable ,psth_badprefer_bad_mean,psth_goodprefer_bad_sem,'lineprops',{[0,0,1]},'transparent',1,'patchSaturation',0.3);

        set(findall(gca, 'Type', 'Line'),'LineWidth',3);


         subplot(4,4,(i-1)*4+3);

      shadedErrorBar(PSTH_Time_TgtOn_stable ,psth_goodprefer_good_nonRF_mean,psth_goodprefer_good_nonRF_sem,'lineprops',{[1,0,0]},'transparent',1,'patchSaturation',0.3);
      hold on
       shadedErrorBar(PSTH_Time_TgtOn_stable ,psth_goodprefer_bad_nonRF_mean,psth_goodprefer_bad_nonRF_sem,'lineprops',{[0,0,1]},'transparent',1,'patchSaturation',0.3);

        set(findall(gca, 'Type', 'Line'),'LineWidth',3);



         subplot(4,4,(i-1)*4+4);

      shadedErrorBar(PSTH_Time_TgtOn_stable ,psth_badprefer_good_nonRF_mean,psth_badprefer_good_nonRF_sem,'lineprops',{[1,0,0]},'transparent',1,'patchSaturation',0.3);
      hold on
       shadedErrorBar(PSTH_Time_TgtOn_stable ,psth_badprefer_bad_nonRF_mean,psth_badprefer_bad_nonRF_sem,'lineprops',{[0,0,1]},'transparent',1,'patchSaturation',0.3);

        set(findall(gca, 'Type', 'Line'),'LineWidth',3);

   %}



     %{

     good_p = auc_RF(auc_RF>0.5 & p_RF<0.05);
     good_p_non = auc_nonRF(auc_RF>0.5 & p_RF<0.05);

     bad_p = auc_RF(auc_RF<0.5 & p_RF<0.05 );
     bad_p_non = auc_nonRF(auc_RF<0.5 & p_RF<0.05);

     
    plot(auc_nonRF,auc_RF,'ok','LineWidth',1);
    hold on
     plot(good_p_non,good_p,'or','LineWidth',3);
     hold on
     plot(bad_p_non,bad_p,'ob','LineWidth',3);
     

     





 end 
  %}
%{
keyboard;

Group = 2;

good = PSTH_TgtOn_Good_RF_stable_Mo{Group};
bad = PSTH_TgtOn_Bad_RF_stable_Mo{Group};

good_nonRF = PSTH_TgtOn_Good_nonRF_stable_Mo{Group};
bad_nonRF = PSTH_TgtOn_Bad_nonRF_stable_Mo{Group};

goodsac = PSTH_SacOn_Good_RF_stable_Mo{Group};
badsac = PSTH_SacOn_Bad_RF_stable_Mo{Group};

goodsac_nonRF = PSTH_TgtOn_Good_nonRF_stable_Mo{Group};
badsac_nonRF = PSTH_TgtOn_Bad_nonRF_stable_Mo{Group};

 p_Tgt = p_Tgt_RF_stable_Mo{Group};
p_Sac = p_Sac_RF_stable_Mo{Group} ;


name = size(good,1);


for ii = 0:name-1

        fig_index = floor(ii/16)+1;
        subplot_index = mod(ii,16)+1;
        figure(fig_index);
        subplot(4,4,subplot_index);
%
        plot(PSTH_Time_TgtOn_stable,good(ii+1,:),'-r');
        hold on
        plot(PSTH_Time_TgtOn_stable,bad(ii+1,:),'-b');
        plot(PSTH_Time_TgtOn_stable,good_nonRF(ii+1,:),'-m' );
        plot(PSTH_Time_TgtOn_stable,bad_nonRF(ii+1,:),'-g' );

      %{
        plot(PSTH_Time_SacOn_stable,goodsac(ii+1,:),'-r');
        hold on
        plot(PSTH_Time_SacOn_stable,badsac(ii+1,:),'-b');
        plot(PSTH_Time_SacOn_stable,goodsac_nonRF(ii+1,1:120),'-m' );
        plot(PSTH_Time_SacOn_stable,badsac_nonRF(ii+1,1:120),'-g' );
        %}
        if p_Tgt(ii+1)<0.05
            Vis = 1;
        else
            Vis = 0;
        end

        if p_Sac(ii+1)<0.05
            Sac = 1;
        else
            Sac = 0;
        end
        
        


       
       % index=index+1;
      title(sprintf("p_Tgt: %d,p_Sac: %d",Vis,Sac));

 end

%}

%{
PSTH_TgtOn_Good_RF_stable_Mo{i} = cell2mat(PSTH_TgtOn_Good_RF_stable(sel)');
            PSTH_TgtOn_Bad_RF_stable_Mo{i} = cell2mat(PSTH_TgtOn_Bad_RF_stable(sel)');

            PSTH_TgtOn_Good_nonRF_stable_Mo{i} = cell2mat(PSTH_TgtOn_Good_nonRF_stable(sel)');
            PSTH_TgtOn_Bad_nonRF_stable_Mo{i} = cell2mat(PSTH_TgtOn_Bad_nonRF_stable(sel)');

%}


end %End of stable value task


%% Flexible Value Task
if AnalysisFV 

    FlexibleValue_Sel = FlexibleValueIndex ==1;

    for i = 1:length(UniqueMol)

         sel = ModuType == UniqueMol(i) & StableValue_Sel;

         SigVal = (p_Tgt_RF_flexible<0.05 | p_Sac_RF_flexible<0.05) & sel ;

         %Cell Type Counts

         ValueCellCountFlexible(i) = sum(SigVal);
         PropValueCellFlexible(i) =ValueCellCountFlexible(i) /sum(sel);

         GoodPref_Vis_RF_flexible = Value_Vis_RF_flexible == 1 & SigVal;
         GoodPref_Vis_RF_Count_flexible(i) = sum(GoodPref_Vis_RF_flexible);
         GoodPref_Vis_RF_Prop_flexible(i) = sum(GoodPref_Vis_RF_flexible)/ValueCellCountFlexible(i);

         GoodPref_Vis_nonRF_flexible = Value_Vis_nonRF_flexible == 1 & SigVal;
         GoodPref_Vis_nonRF_Count_flexible(i) = sum(GoodPref_Vis_nonRF_flexible);
         GoodPref_Vis_nonRF_Prop_flexible(i) = sum(GoodPref_Vis_nonRF_flexible)/ValueCellCountFlexible(i);

         %PSTH
       

         PSTH_TgtOn_Prefer_RF_flexible_sel = cell2mat(PSTH_TgtOn_Prefer_RF_flexible(sel & p_Tgt_RF_flexible<0.05)');

         PSTH_TgtOn_Prefer_RF_flexible_mean(i,:) = mean(PSTH_TgtOn_Prefer_RF_flexible_sel,1);

         PSTH_TgtOn_Prefer_RF_flexible_sem(i,:) = std(PSTH_TgtOn_Prefer_RF_flexible_sel,1)/sqrt(sum(SigVal));


          PSTH_TgtOn_NonPrefer_RF_flexible_sel = cell2mat(PSTH_TgtOn_nonPrefer_RF_flexible(sel & p_Tgt_RF_flexible<0.05)');

         PSTH_TgtOn_NonPrefer_RF_flexible_mean(i,:) = mean(PSTH_TgtOn_NonPrefer_RF_flexible_sel ,1);

         PSTH_TgtOn_NonPrefer_RF_flexible_sem(i,:) = std(PSTH_TgtOn_NonPrefer_RF_flexible_sel ,1)/sqrt(sum(SigVal));


         PSTH_SacOn_Prefer_RF_flexible_sel = cell2mat(PSTH_SacOn_Prefer_RF_flexible(sel & p_Sac_RF_flexible<0.05)');

         PSTH_SacOn_Prefer_RF_flexible_mean(i,:) = mean(PSTH_SacOn_Prefer_RF_flexible_sel,1);

         PSTH_SacOn_Prefer_RF_flexible_sem(i,:) = std(PSTH_SacOn_Prefer_RF_flexible_sel,1)/sqrt(sum(SigVal));


          PSTH_SacOn_NonPrefer_RF_flexible_sel = cell2mat(PSTH_SacOn_nonPrefer_RF_flexible(sel & p_Sac_RF_flexible<0.05)');

         PSTH_SacOn_NonPrefer_RF_flexible_mean(i,:) = mean(PSTH_SacOn_NonPrefer_RF_flexible_sel ,1);

         PSTH_SacOn_NonPrefer_RF_flexible_sem(i,:) = std(PSTH_SacOn_NonPrefer_RF_flexible_sel ,1)/sqrt(sum(SigVal));

%{
          PSTH_TgtOn_Prefer_RF_flexible{index} =PSTH_TgtOn_Prefer_RF_flexible{ScreenCell} ;
 PSTH_TgtOn_nonPrefer_RF_flexible{index} = PSTH_TgtOn_nonPrefer_RF_flexible{ScreenCell} ;
 PSTH_TgtOn_Prefer_nonRF_flexible{index} = PSTH_TgtOn_Prefer_nonRF_flexible{ScreenCell};
 PSTH_TgtOn_nonPrefer_nonRF_flexible{index} = PSTH_TgtOn_nonPrefer_nonRF_flexible{ScreenCell} ;

 PSTH_SacOn_Prefer_RF_flexible{index} = PSTH_SacOn_Prefer_RF_flexible{ScreenCell};
 PSTH_SacOn_nonPrefer_RF_flexible{index} = PSTH_SacOn_nonPrefer_RF_flexible{ScreenCell};
 PSTH_SacOn_Prefer_nonRF_flexible{index} = PSTH_SacOn_Prefer_nonRF_flexible{ScreenCell} ;
 PSTH_SacOn_nonPrefer_nonRF_flexible{index} = PSTH_SacOn_nonPrefer_nonRF_flexible{ScreenCell};



PSTH_TgtOn_Good_RF_flexible{index} =PSTH_TgtOn_Good_RF_flexible{ScreenCell};
PSTH_TgtOn_Bad_RF_flexible{index} = PSTH_TgtOn_Bad_RF_flexible{ScreenCell} ;
        PSTH_TgtOn_Good_nonRF_flexible{index} = PSTH_TgtOn_Good_nonRF_flexible{ScreenCell};
        PSTH_TgtOn_Bad_nonRF_flexible{index} = PSTH_TgtOn_Bad_nonRF_flexible{ScreenCell};

        PSTH_SacOn_Good_RF_flexible{index} = PSTH_SacOn_Good_RF_flexible{ScreenCell} ;
        PSTH_SacOn_Bad_RF_flexible{index} = PSTH_SacOn_Bad_RF_flexible{ScreenCell};
        PSTH_SacOn_Good_nonRF_flexible{index} = PSTH_SacOn_Good_nonRF_flexible{ScreenCell} ;
        PSTH_SacOn_Bad_nonRF_flexible{index} = PSTH_SacOn_Bad_nonRF_flexible{ScreenCell} ;

%}

    end


    

    %p_Tgt_RF_stable
    








end %End of flexible value task


%% One DR
if AnalysisOneDR

    OneDRsel = OneDRIndex == 1;

     for i = 1:length(UniqueMol)

         sel = ModuType == UniqueMol(i) &  OneDRsel;

          SigVal = (p_Tgt_RF<0.05 | p_Sac_RF<0.05 | p_preTgt_RF < 0.05) & sel ;
          NonSigVal = ~(p_Tgt_RF<0.05 | p_Sac_RF<0.05 | p_preTgt_RF < 0.05)& sel ;

          SigVal_NonRF = (p_Tgt_nonRF<0.05 | p_Sac_nonRF<0.05 | p_preTgt_nonRF < 0.05) & sel ;
          NonSigVal_NonRF = ~(p_Tgt_nonRF<0.05 | p_Sac_nonRF<0.05 | p_preTgt_nonRF < 0.05)& sel ;

          SigVal_All = SigVal|SigVal_NonRF;
         % NonSigVall_All = ~SigVal_All;


         %% Cell Type Counts

         ValueCellCount_oneDR(i) = sum(SigVal_All);

         NonValueCellCount_oneDR(i) = sum(sel) - sum(SigVal_All);
         
         PropValueCell_OneDR(i) =ValueCellCount_oneDR(i) /sum(sel);
         %
         %% Value 
         Value_RF_Only = SigVal & NonSigVal_NonRF;
         Value_nonRF_Only = SigVal_NonRF & NonSigVal;

         Value_Both = SigVal & SigVal_NonRF;

         ValueRF_Count(i) = sum(Value_RF_Only);
         Value_nonRF_Count(i) = sum(Value_nonRF_Only);
         Value_Both_Count(i) = sum(Value_Both);

         GoodPrefer_RF = TgtAUC_RF>0.5 & SigVal;
         GoodPrefer_NonRF = TgtAUC_NonRF>0.5 & SigVal_NonRF;

         BadPrefer_RF = TgtAUC_RF<0.5 & SigVal;
         BadPrefer_NonRF = TgtAUC_NonRF<0.5 & SigVal_NonRF;

         ValuePrefer_Consistent = (GoodPrefer_RF & GoodPrefer_NonRF)|(BadPrefer_RF & BadPrefer_NonRF);
         ValuePrefer_InConsistent = (GoodPrefer_RF & BadPrefer_NonRF)|(BadPrefer_RF & GoodPrefer_NonRF);

         ValueConsistent_Count(i) = sum(ValuePrefer_Consistent);
         ValueInConsistent_Count(i) = sum(ValuePrefer_InConsistent);

         GoodPreferCount(i) = sum(GoodPrefer_RF|GoodPrefer_NonRF);
         BadPreferCount(i) = sum(BadPrefer_RF|BadPrefer_NonRF);

         ValueSummary_OneDR(i,:)=[ValueCellCount_oneDR(i),NonValueCellCount_oneDR(i),ValueRF_Count(i),Value_nonRF_Count(i),Value_Both_Count(i),ValueConsistent_Count(i),ValueInConsistent_Count(i),GoodPreferCount(i),BadPreferCount(i)];







         GoodPref_Vis_RF_OneDR = Value_Vis_RF_OneDR == 1 & SigVal;
         GoodPref_Vis_RF_OneDR(i) = sum(GoodPref_Vis_RF_OneDR);
         GoodPref_Vis_RF_Prop_OneDR(i) = sum(GoodPref_Vis_RF_OneDR)/ValueCellCount_oneDR(i) ;

         GoodPref_Vis_nonRF_OneDR = Value_Vis_nonRF_OneDR == 1 & SigVal;
         GoodPref_Vis_nonRF_Count_OneDR(i) = sum(GoodPref_Vis_nonRF_OneDR);
         GoodPref_Vis_nonRF_Prop_OneDR(i) = sum(GoodPref_Vis_nonRF_OneDR)/ValueCellCount_oneDR(i);


        % ValueRFOnly = 

         %}
          PSTH_TgtOn_Prefer_RF_sel = cell2mat(PSTH_Prefer_TgtOn_RF(sel & (p_Tgt_RF<0.05 |p_preTgt_RF<0.05))');

          Counts = sum(sel & (p_Tgt_RF<0.05 |p_preTgt_RF<0.05));

         PSTH_TgtOn_Prefer_RF_mean(i,:) = mean(PSTH_TgtOn_Prefer_RF_sel,1) ;

         PSTH_TgtOn_Prefer_RF_sem(i,:) = std(PSTH_TgtOn_Prefer_RF_sel,1)/sqrt(Counts);


          PSTH_TgtOn_NonPrefer_RF_sel = cell2mat(PSTH_nonPrefer_TgtOn_RF(sel & (p_Tgt_RF<0.05 |p_preTgt_RF<0.05))');

         PSTH_TgtOn_NonPrefer_RF_mean(i,:) = mean(PSTH_TgtOn_NonPrefer_RF_sel ,1);

         PSTH_TgtOn_NonPrefer_RF_sem(i,:) = std(PSTH_TgtOn_NonPrefer_RF_sel ,1)/sqrt(Counts);


         Counts_Sac = sum(sel & p_Sac_RF<0.05);
         PSTH_SacOn_Prefer_RF_sel = cell2mat(PSTH_Prefer_SacOn_RF(sel & p_Sac_RF<0.05)');

         PSTH_SacOn_Prefer_RF_mean(i,:) = mean(PSTH_SacOn_Prefer_RF_sel,1);

         PSTH_SacOn_Prefer_RF_sem(i,:) = std(PSTH_SacOn_Prefer_RF_sel,1)/sqrt(Counts_Sac );


          PSTH_SacOn_NonPrefer_RF_sel = cell2mat(PSTH_nonPrefer_SacOn_RF(sel & p_Sac_RF_stable<0.05)');
         

         % PSTH_SacOn_NonPrefer_RF_sel = cell2mat(PSTH_nonPrefer_SacOn_RF(sel & (p_Tgt_RF<0.05 |p_preTgt_RF<0.05))');


         PSTH_SacOn_NonPrefer_RF_mean(i,:) = mean(PSTH_SacOn_NonPrefer_RF_sel ,1);

         PSTH_SacOn_NonPrefer_RF_sem(i,:) = std(PSTH_SacOn_NonPrefer_RF_sel ,1)/sqrt(Counts_Sac );



        
       



     end



end %End of OneDR


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

    figtitlestr{FigureIndex}='PSTH_OptoModulation_Type';
    fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[50,100, 800,800],'Name',figtitlestr{FigureIndex});

    set(fig{FigureIndex}, 'PaperUnits', 'inches');
%formatFig(fig{FigureIndex}, [1.5 1.5], 'nature');
    formatFig(fig{FigureIndex}, [2 2], 'nature');
    AverageRegion = [0,100];
        area(AverageRegion,[1,1],...
        'FaceColor','k','FaceAlpha',.1,'EdgeAlpha',.1,'HandleVisibility','off');
        hold on

        AverageRegion = [0,100];
        area(AverageRegion,[-1,-1],...
        'FaceColor','k','FaceAlpha',.1,'EdgeAlpha',.1,'HandleVisibility','off');


    for i = 1:length(UniqueMol)
        %subplot(1,length(UniqueMol),i);
       % subplot(2,2,i);
      

        shadedErrorBar(PSTH_Opto_Time,PSTH_Stim_Mean_Group(i,:),PSTH_Stim_Sem_Group(i,:),'lineprops',{ColorOrder(i,:)},'transparent',1,'patchSaturation',0.3)
        set(findall(gca, 'Type', 'Line'),'LineWidth',0.5);
        
       
        
      %  title(ModulString(i));
        

      
        % legstring{i}=sprintf('N=%d', ModuTypeCount(i));

    end

   

     legend({'Fast Ex','Slow Ex','Inh','No effect'});
       set(gca,'LineWidth',0.5,'FontSize',5);%,'FontWeight','Bold');
        xlabel('Time from Opto onset(ms)');
        ylabel('Normalized Responses');

        box off;

     %% 

%{
    FigureStartNum = FigureStartNum+1;
    FigureIndex = FigureIndex+1;

    figtitlestr{FigureIndex}='Latency_Dist';
    fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[50,100, 600,400],'Name',figtitlestr{FigureIndex});

     set(fig{FigureIndex}, 'PaperUnits', 'inches');
    %formatFig(fig{FigureIndex}, [1.5 1.5], 'nature');
    formatFig(fig{FigureIndex}, [2 2], 'nature');


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
        %% Latency Distribution
         FigureStartNum = FigureStartNum+1;
    FigureIndex = FigureIndex+1;

        figtitlestr{FigureIndex}='Latency_Distribution';
   fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[50,100, 600,400],'Name',figtitlestr{FigureIndex});
   set(fig{FigureIndex}, 'PaperUnits', 'inches');
%formatFig(fig{FigureIndex}, [1.5 1.5], 'nature');
    formatFig(fig{FigureIndex}, [2 2], 'nature');

   Edges =[0:5:100];
   Edges_Fast =[0:5:10];

   subplot(2,1,1);
   histogram(Latency_Ex,Edges,'Facecolor',ColorOrder(2,:),'LineWidth',0.5);
   hold on
   histogram(Latency_Ex,Edges_Fast,'Facecolor',ColorOrder(1,:),'LineWidth',0.5);
 %  set(gca, 'YScale', 'log');
  % xlabel('Latency');
   ylabel('Number of neurons');
   set(gca,'LineWidth',0.5,'FontSize',5);%,'FontWeight','Bold');
   box off



   subplot(2,1,2);
   histogram(Latency_In,Edges,'Facecolor',ColorOrder(3,:),'LineWidth',0.5);
  % set(gca, 'YScale', 'log');
   xlabel('Latency');
   ylabel('Number of neurons');
   set(gca,'LineWidth',0.5,'FontSize',5);%,'FontWeight','Bold');
    box off


%% Pie plot 
    FigureStartNum = FigureStartNum+1;
    FigureIndex = FigureIndex+1;

   


        figtitlestr{FigureIndex}='Modulation';
   fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[50,100, 300,300],'Name',figtitlestr{FigureIndex});
% set(fig{FigureIndex}, 'PaperUnits', 'inches');
%formatFig(fig{FigureIndex}, [1.5 1.5], 'nature');
  %  formatFig(fig{FigureIndex}, [0.5 0.5], 'nature');

   h = piechart(ModuTypeCount,["FastExc","SlowExc","Inh","No-Effect"])%,'FontSize',5);%,["FastExc","SlowExc","Inh","No-Effect"]);
   %
   %{
  for i = 1:length(ModuTypeCount)
      modu_label{i} = sprintf('N = %d',ModuTypeCount(i));
  end
  h.Labels = modu_label;
   %}
   %h.Labels = '';
  h.ColorOrder = ColorOrder;
   h.StartAngle = -30;
   %
%{

    ExcitationNum = sum(Excitation & Sig_m);
InhibitionNum = sum(Inhibition & Sig_m);
NonSigNum = sum(~Sig_m);

%}
%% Modulation Strength
%{
MS_Mean_Group(i) = mean(abs(SM(sel)));
    MS_Sem_Group(i) = std(abs(SM(sel)))/sqrt(sum(sel));
%}
%{
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
%{
FigureStartNum = FigureStartNum+1;
    FigureIndex = FigureIndex+1;

     MS_TypeOrder = [0.8824    0.7451    0.4157;
                        0.2510    0.6902    0.6510;
                        0.6471    0.5020    0.9020;
                        0.5,0.5,0.5];

    figtitlestr{FigureIndex}='PSTH_MemorySac_Type';
    fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[50,100, 500,600],'Name',figtitlestr{FigureIndex});
     b = bar(TypeCounts,'stacked');

     for i = 1:length(b)
         subb = b(i);
         subb.FaceColor = 'Flat';
         subb.CData = [MS_TypeOrder(i,:)];
     end

     box off
%     xticks(UniqueMol(1:4));
    xticklabels({'Fast Ex','Slow Ex','In','No'});
    ylabel('Proportion');
    legend({'Vis-Only','Vis-Motor','Motor'});

     set(gca,'LineWidth',3,'FontSize',20,'FontWeight','Bold');
%}

 %% Population psth

 %for i = 1: length(UniqueMol)
    FigureStartNum = FigureStartNum+1;
    FigureIndex = FigureIndex+1;

    figtitlestr{FigureIndex}= sprintf('PSTH_MemorySac_Type:%d',UniqueMol(1));
    fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[50,100, 600,600],'Name',figtitlestr{FigureIndex});
     set(fig{FigureIndex}, 'PaperUnits', 'inches');
%formatFig(fig{FigureIndex}, [1.5 1.5], 'nature');
    formatFig(fig{FigureIndex}, [2 3.5], 'nature');

    typestring = {'Visual','Visual-Movement','Movement','Others'};

   for j = 1:length(CellTypeUnique)-1
        subplot(length(CellTypeUnique)-1,2,2*(j-1)+1);
      % subplot(1,2,1);

        %% Target onset 


        shadedErrorBar(PSTH_Time_TgtOn ,PSTH_Tgt_Mean{1,j},PSTH_Tgt_Sem{1,j},'lineprops',{[128,128,128]/255},'transparent',1,'patchSaturation',0.2)
        set(findall(gca, 'Type', 'Line'),'LineWidth',0.5);
   %      ylim([-0.5,3]);
       xlim([-250,300]);
       hold on
       
       ymax = ceil(max(max(PSTH_Tgt_Mean{1,j}+PSTH_Tgt_Sem{1,j}),max(PSTH_Sac_Mean{1,j}+PSTH_Sac_Sem{1,j})));
       ymin = floor(min(min(PSTH_Tgt_Mean{1,j}-PSTH_Tgt_Sem{1,j}),min(PSTH_Sac_Mean{1,j}-PSTH_Sac_Sem{1,j})));
       ylim([ymin,ymax]);
       plot([0,0],[ymin,ymax],'--k','LineWidth',0.5);
       if j == 2
            ylabel('Z-Score');
        end
        if j == 3
            xlabel('Time from event onset (ms)');
        end
        title(typestring(j));
        set(gca,'FontSize',5);
        

        %% Saccade onset
        subplot(length(CellTypeUnique)-1,2,2*(j-1)+2);

        
      % subplot(1,2,2);

        %Saccade on 

        shadedErrorBar(PSTH_Time_SacOn ,PSTH_Sac_Mean{1,j},PSTH_Sac_Sem{1,j},'lineprops',{[128,128,128]/255},'transparent',1,'patchSaturation',0.2)
        set(findall(gca, 'Type', 'Line'),'LineWidth',0.5);

         xlim([-450,200]);
        ylim([ymin,ymax]);
        hold on
        plot([0,0],[ymin,ymax],'--k','LineWidth',0.5);
        set(gca,'FontSize',5);
      %  ylim([-0.5,3]);


%PSTH_Tgt_Mean{i,j}


   end

% end
%Bar plot for fast activated neuron only
 FigureStartNum = FigureStartNum+1;
    FigureIndex = FigureIndex+1;

     figtitlestr{FigureIndex}='PSTH_MemorySac_ProjectionNeuron';
    fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[50,100, 500,600],'Name',figtitlestr{FigureIndex});
    set(fig{FigureIndex}, 'PaperUnits', 'inches');
%formatFig(fig{FigureIndex}, [1.5 1.5], 'nature');
    formatFig(fig{FigureIndex}, [2 3.5], 'nature');
    subplot(2,1,1);
     b = bar(TypeCounts(1,1:3),'FaceColor','k');

     
%         b.FaceColor = 'Flat';
        % b.CData = MS_TypeOrder(1,:);
     

     box off
     xticks(UniqueMol);
    xticklabels({'Vis','Vis-Mov','Mov'});
    ylabel('Number of neurons');
    
  %  legend({'Vis-Only','Vis-Motor','Motor'});

     set(gca,'LineWidth',0.5,'FontSize',5);%,'FontWeight','Bold');

subplot(2,1,2);
%% Frequency
      b = bar(TypeFreq(1,1:3),'FaceColor','k');

     
%         b.FaceColor = 'Flat';
        % b.CData = MS_TypeOrder(1,:);
     

     box off
     xticks(UniqueMol);
    xticklabels({'Vis','Vis-Mov','Mov'});
    ylabel('Proportion');
  %  legend({'Vis-Only','Vis-Motor','Motor'});

     set(gca,'LineWidth',0.5,'FontSize',5);%,'FontWeight','Bold');
 






 end %End of the if plot memory saccade

 %% 
 if PlotDelaySaccade

     FigureStartNum = FigureStartNum+1;
    FigureIndex = FigureIndex+1;


    figtitlestr{FigureIndex}='DelaySaccadeVectorDist';
    fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[50,100, 600,600],'Name',figtitlestr{FigureIndex});
   
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


 %%Stable Value    
  if PlotStableValue 


     % PSTH_TgtOn_Good_RF_stable{index}
%{
        FigureStartNum = FigureStartNum+1;
        FigureIndex = FigureIndex+1;


        figtitlestr{FigureIndex}='PSTH_Value_Raw';
fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[50,100, 500,600],'Name',figtitlestr{FigureIndex});
%}
        %%
 FigureStartNum = FigureStartNum+1;
        FigureIndex = FigureIndex+1;


        figtitlestr{FigureIndex}='ValueInformationSummary';
        fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[50,100, 500,600],'Name',figtitlestr{FigureIndex});

        for i = 1:length(UniqueMol)
            subplot(4,4,i);
            
            b = bar([1,2],ValueInfo_Summary(i,1:2));
            b.FaceColor='flat';
            b.CData = ColorOrder(i,:);
            xticklabels({'Value','NonValue'});

            subplot(4,4,4+i);
            b = bar([1,2,3],ValueInfo_Summary(i,3:5));
            b.FaceColor='flat';
            b.CData = ColorOrder(i,:);
            xticklabels({'ValueContraOnly','ValueIspiOnly','ValueBoth'});

            subplot(4,4,4*2+i);
            b = bar([1,2],ValueInfo_Summary(i,6:7));
            b.FaceColor='flat';
            b.CData = ColorOrder(i,:);
            xticklabels({'GoodPrefer','BadPrefer'});

            subplot(4,4,4*3+i);
            b = bar([1,2],ValueInfo_Summary(i,8:9));
            b.FaceColor='flat';
            b.CData = ColorOrder(i,:);
            xticklabels({'Consistent','InConsistent'});

             

           


        end
 % ValueInfo_Summary(i,:)=[ValueCellCount(i),NonValueCellCount(i),Value_RFOnlyCount(i),Value_nonRFOnlyCount(i),Value_BothCount(i),GoodPrefer_Count(i),BadPrefer_Count(i), ConsistentCount(i), InConsistentCount(i)];

     
 FigureStartNum = FigureStartNum+1;
        FigureIndex = FigureIndex+1;


        figtitlestr{FigureIndex}='ValueInformationSummary_EvsOthers';
        fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[50,100, 500,600],'Name',figtitlestr{FigureIndex});

        % Value Proportion
        FastEx = ValueInfo_Summary(1,1:2);
        Others = ValueInfo_Summary(2,1:2)+ValueInfo_Summary(3,1:2)+ValueInfo_Summary(4,1:2);
            subplot(4,1,1);
            
            %b = bar([1,2],FastEx);
           % b.FaceColor='flat';
           % b.CData = ColorOrder(1,:);

           pie(FastEx);
           colormap([1,0.8,0;0.5,0.5,0.5;0,0.8,1]);
           legend({'Value','Non-Value'})

           % xticklabels({'Value','NonValue'});
           % ylabel('Number of Neurons');
            box off

             set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');


          %  subplot(2,2,2);
            
           % b = bar([1,2],Others);
           % b.FaceColor='flat';
           % b.CData = ColorOrder(4,:);
         %   xticklabels({'Value','NonValue'});
         %   ylabel('Number of Neurons');
          %  box off

         %    set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');



        % RF vs NonRF
        FastEx = ValueInfo_Summary(1,3:5);
        Others = ValueInfo_Summary(2,3:5)+ValueInfo_Summary(3,3:5)+ValueInfo_Summary(4,3:5);
        subplot(4,1,2);
            % b = bar([1,2,3],FastEx);
           % b.FaceColor='flat';
           % b.CData = ColorOrder(1,:);
           pie(FastEx)
          % colormap([0.5,0,0;0,0.5,0])
           legend({'ValueContraOnly','ValueIspiOnly','ValueBoth'});
            %xticklabels({'ValueContraOnly','ValueIspiOnly','ValueBoth'});
           % ylabel('Number of Neurons');
          %  box off

             set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');

%{
           subplot(4,1,3);
            b = bar([1,2,3],Others);
            b.FaceColor='flat';
            b.CData = ColorOrder(4,:);
            xticklabels({'ValueContraOnly','ValueIspiOnly','ValueBoth'});
            ylabel('Number of Neurons');
            box off

             set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
%}
        
           %Good vs Bad

           FastEx = ValueInfo_Summary(1,6:7);
        Others = ValueInfo_Summary(2,6:7)+ValueInfo_Summary(3,6:7)+ValueInfo_Summary(4,6:7);

           subplot(4,1,3);
           % b = bar([1,2],FastEx);
           % b.FaceColor='flat';
           % b.CData = ColorOrder(1,:);
           % xticklabels({'GoodPrefer','BadPrefer'});
           % ylabel('Number of Neurons');
           % box off

          %   set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
          pie(FastEx);
          legend({'GoodPrefer','BadPrefer'});



%{
             subplot(4,1,6);
            b = bar([1,2],Others);
            b.FaceColor='flat';
            b.CData = ColorOrder(4,:);
            xticklabels({'GoodPrefer','BadPrefer'});

            ylabel('Number of Neurons');
            box off

             set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
%}
%Consistent vs InConsistent

           FastEx = ValueInfo_Summary(1,8:9);
        Others = ValueInfo_Summary(2,8:9)+ValueInfo_Summary(3,8:9)+ValueInfo_Summary(4,8:9);

           subplot(4,1,4);
           % b = bar([1,2],FastEx);
           pie(FastEx);
           legend({'Consistent','InConsistent'});
           %{
            b.FaceColor='flat';
            b.CData = ColorOrder(1,:);
            xticklabels({'Consistent','InConsistent'});
            ylabel('Number of Neurons');
            box off

             set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
           %}
%{
             subplot(4,2,8);
            b = bar([1,2],Others);
            b.FaceColor='flat';
            b.CData = ColorOrder(4,:);
            xticklabels({'Consistent','InConsistent'});
            ylabel('Number of Neurons');
            box off

             set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
%}
        %% 


        FigureStartNum = FigureStartNum+1;
        FigureIndex = FigureIndex+1;


        figtitlestr{FigureIndex}='PSTH_Contralateral_Value';
        fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[50,100, 500,600],'Name',figtitlestr{FigureIndex});
     

  
    for i = 1:length(UniqueMol)
    subplot(length(UniqueMol),4,(i-1)*4+1)
     shadedErrorBar(PSTH_Time_TgtOn_stable,PSTH_TgtOn_Prefer_RF_stable_mean(i,:),PSTH_TgtOn_Prefer_RF_stable_sem(i,:),'lineprops',{[255,194,10]/255},'transparent',1,'patchSaturation',0.3)
     set(findall(gca, 'Type', 'Line'),'LineWidth',3);


    hold on
    shadedErrorBar(PSTH_Time_TgtOn_stable,PSTH_TgtOn_NonPrefer_RF_stable_mean(i,:),PSTH_TgtOn_NonPrefer_RF_stable_sem(i,:),'lineprops',{[101,175,236]/255},'transparent',1,'patchSaturation',0.3);
   % ylim([-1,5]);

    set(findall(gca, 'Type', 'Line'),'LineWidth',3);
    xlabel('Time from Tgt Onset(ms)');
    ylabel('Z-Score');
    set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
    

    subplot(length(UniqueMol),4,(i-1)*4+2)
    shadedErrorBar(PSTH_Time_SacOn_stable,PSTH_SacOn_Prefer_RF_stable_mean(i,:),PSTH_SacOn_Prefer_RF_stable_sem(i,:),'lineprops',{[255,194,10]/255},'transparent',1,'patchSaturation',0.3);
    hold on
    shadedErrorBar(PSTH_Time_SacOn_stable,PSTH_SacOn_NonPrefer_RF_stable_mean(i,:),PSTH_SacOn_NonPrefer_RF_stable_sem(i,:),'lineprops',{[101,175,236]/255},'transparent',1,'patchSaturation',0.3);
   % ylim([-1,5]);
   set(findall(gca, 'Type', 'Line'),'LineWidth',3);

   xlabel('Time from Sac Onset(ms)');
   ylabel('Z-Score');
    set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');


     subplot(length(UniqueMol),4,(i-1)*4+3)
     shadedErrorBar(PSTH_Time_TgtOn_stable,PSTH_TgtOn_Prefer_nonRF_stable_mean(i,:),PSTH_TgtOn_Prefer_nonRF_stable_sem(i,:),'lineprops',{[255,194,10]/255},'transparent',1,'patchSaturation',0.3)
     set(findall(gca, 'Type', 'Line'),'LineWidth',3);


    hold on
    shadedErrorBar(PSTH_Time_TgtOn_stable,PSTH_TgtOn_NonPrefer_nonRF_stable_mean(i,:),PSTH_TgtOn_NonPrefer_nonRF_stable_sem(i,:),'lineprops',{[101,175,236]/255},'transparent',1,'patchSaturation',0.3);
  %  ylim([-1,5]);

    set(findall(gca, 'Type', 'Line'),'LineWidth',3);
    xlabel('Time from Tgt Onset(ms)');
    ylabel('Z-Score');
    set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
    

    subplot(length(UniqueMol),4,(i-1)*4+4)
    shadedErrorBar(PSTH_Time_SacOn_stable,PSTH_SacOn_Prefer_nonRF_stable_mean(i,:),PSTH_SacOn_Prefer_nonRF_stable_sem(i,:),'lineprops',{[255,194,10]/255},'transparent',1,'patchSaturation',0.3);
    hold on
    shadedErrorBar(PSTH_Time_SacOn_stable,PSTH_SacOn_NonPrefer_nonRF_stable_mean(i,:),PSTH_SacOn_NonPrefer_nonRF_stable_sem(i,:),'lineprops',{[101,175,236]/255},'transparent',1,'patchSaturation',0.3);
  %  ylim([-1,5]);
   set(findall(gca, 'Type', 'Line'),'LineWidth',3);

   xlabel('Time from Sac Onset(ms)');
   ylabel('Z-Score');
    set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');

    end 

%% 
     FigureStartNum = FigureStartNum+1;
        FigureIndex = FigureIndex+1;


        figtitlestr{FigureIndex}='ValuePreferenceSummary';
        fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[50,100, 500,600],'Name',figtitlestr{FigureIndex});

        subplot(2,1,1)
        b = bar(PropValueCell);
        b.FaceColor = 'flat';
        b.CData = ColorOrder;

        xticklabels({'ExFast','ExSlow','In','NoEffect'});
        ylabel('Proportion of value selective cells');
        box off;
        set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
        


        subplot(2,1,2)
        b=bar(GoodPref_Vis_RF_Prop);
        b.FaceColor = 'flat';
        b.CData = ColorOrder;
        xticklabels({'ExFast','ExSlow','In','NoEffect'});
        ylabel('Proportion of good preference cells');
        box off;
        hold on 
        plot([0,5],[0.5,0.5],'--k','LineWidth',3);
        set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');


%%   
       FigureStartNum = FigureStartNum+1;
        FigureIndex = FigureIndex+1;

          figtitlestr{FigureIndex}='AUC_Summary';
        fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[50,100, 500,600],'Name',figtitlestr{FigureIndex});

        for i = 1:length(UniqueMol)
            subplot(2,length(UniqueMol),i);
            auc_rf = AUC_Stable_RF_Type{i};
            auc_rf_p = p_Tgt_RF_stable_Mo{i};

            auc_nonrf = AUC_Stable_nonRF_Type{i};
            auc_nonrf_p = p_Tgt_nonRF_stable_Mo{i};

            auc_rf_sig = auc_rf(auc_rf_p<0.05);
            auc_rf_nonsig = auc_rf(~(auc_rf_p<0.05));

            edges = [0:0.1:1];
            values_rf_sig = histcounts(auc_rf_sig,edges);
            values_rf_nonsig = histcounts(auc_rf_nonsig,edges);


            b = bar(edges(1:end-1)+0.05,[values_rf_sig;values_rf_nonsig]','stack','LineWidth',3);
            b(1).FaceColor='r';
            b(2).FaceColor = 'w';
            set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
            legend({sprintf('Sig(N=%d)',length(auc_rf_sig)),sprintf('NonSig(N=%d)',length(auc_rf_nonsig))});
            title('AUC for RF');
            hold on;
            
            box off


            subplot(2,length(UniqueMol),length(UniqueMol)+i);

            auc_nonrf_sig = auc_nonrf(auc_nonrf_p<0.05);
            auc_nonrf_nonsig = auc_nonrf(~(auc_nonrf_p<0.05));

            edges = [0:0.1:1];
            values_nonrf_sig = histcounts(auc_nonrf_sig,edges);
            values_nonrf_nonsig = histcounts(auc_nonrf_nonsig,edges);

            b = bar(edges(1:end-1)+0.05,[values_nonrf_sig;values_nonrf_nonsig]','stack','LineWidth',3);
            b(1).FaceColor='r';
            b(2).FaceColor = 'w';
            set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
            legend({sprintf('Sig(N=%d)',length(auc_nonrf_sig)),sprintf('NonSig(N=%d)',length(auc_nonrf_nonsig))});
            title('AUC for nonRF');
            hold on;
            
            box off

        end


        FigureStartNum = FigureStartNum+1;
        FigureIndex = FigureIndex+1;

           figtitlestr{FigureIndex}='AUC_RFNONRF_CompareSummary';
        fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[50,100, 500,600],'Name',figtitlestr{FigureIndex});

        for i = 1:length(UniqueMol)
           
            auc_rf = AUC_Stable_RF_Type{i};
            auc_rf_p = p_Tgt_RF_stable_Mo{i};

            auc_nonrf = AUC_Stable_nonRF_Type{i};
            auc_nonrf_p = p_Tgt_nonRF_stable_Mo{i};

            auc_rf_sig = auc_rf(auc_rf_p<0.05|auc_nonrf_p<0.05);
            auc_rf_nonsig = auc_rf(~(auc_rf_p<0.05|auc_nonrf_p<0.05));

            auc_nonrf_sig = auc_nonrf(auc_rf_p<0.05|auc_nonrf_p<0.05);
            auc_nonrf_nonsig = auc_nonrf(~(auc_rf_p<0.05|auc_nonrf_p<0.05));

            subplot(2,2,i);
            plot([0.5,0.5],[0,1],'--k','LineWidth',3);
            plot([0,1],[0.5,0.5],'--k','LineWidth',3);

            plot(auc_nonrf_nonsig,auc_rf_nonsig,'ok');
            hold on
            plot(auc_nonrf_sig,auc_rf_sig,'.k','MarkerSize',25);

            xlim([0,1]);
            ylim([0,1]);
            hold on
            
            xlabel('AUC for Ipsilateral value responses');
            ylabel('AUC for Contralateral value responses');
            set(gca,'FontSize',15,'FontWeight','Bold');

           








        end





        %{
 p_Tgt_RF_stable_Mo{i}= p_Tgt_RF_stable(sel);
         p_Sac_RF_stable_Mo{i}= p_Sac_RF_stable(sel);

         p_Tgt_nonRF_stable_Mo{i} = p_Tgt_nonRF_stable(sel);


         AUC_Stable_RF_Type{i}=AUC_Stable_RF( SigVis);
        AUC_Stable_nonRF_Type{i}=AUC_Stable_nonRF( SigVis);
        %}
        

%{
%% Distribution of value ROC
 FigureStartNum = FigureStartNum+1;
 FigureIndex = FigureIndex+1;

 figtitlestr{FigureIndex}='ROC_Distribution_Stable_RF';
 fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[50,100, 500,600],'Name',figtitlestr{FigureIndex});
  for i = 1:length(UniqueMol)

      subplot(2,2,i);
      edges = [0:0.1:1];
      histogram(AUC_Stable_RF_Type{i},edges);
  end


FigureStartNum = FigureStartNum+1;
 FigureIndex = FigureIndex+1;

 figtitlestr{FigureIndex}='ROC_Distribution_Stable_nonRF';
 fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[50,100, 500,600],'Name',figtitlestr{FigureIndex});
  for i = 1:length(UniqueMol)

      subplot(2,2,i);
      edges = [0:0.1:1];
      histogram(AUC_Stable_nonRF_Type{i},edges);
  end

FigureStartNum = FigureStartNum+1;
 FigureIndex = FigureIndex+1;

 figtitlestr{FigureIndex}='RF_non_RF_Compare_Stable';
 fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[50,100, 500,600],'Name',figtitlestr{FigureIndex});
  for i = 1:length(UniqueMol)
    subplot(2,2,i);
      scatter(AUC_Stable_nonRF_Type{i},AUC_Stable_RF_Type{i});

      
     
  end

% AUC_Stable_RF_Type{i}=AUC_Stable_RF(sel);
% AUC_Stable_nonRF_Type{i}=AUC_Stable_nonRF(sel);
%}
  end %End of if plot stable value


   %%Flexible
  if PlotFlexibleValue 


        FigureStartNum = FigureStartNum+1;
        FigureIndex = FigureIndex+1;


        figtitlestr{FigureIndex}='PSTH_Contralateral_FlexibleValue';
        fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[50,100, 500,600],'Name',figtitlestr{FigureIndex});
     

  
    for i = 1:length(UniqueMol)
    subplot(length(UniqueMol),2,(i-1)*2+1)
     shadedErrorBar(PSTH_Time_TgtOn_flexible,PSTH_TgtOn_Prefer_RF_flexible_mean(i,:),PSTH_TgtOn_Prefer_RF_flexible_sem(i,:),'lineprops',{[255,194,10]/255},'transparent',1,'patchSaturation',0.3)
     set(findall(gca, 'Type', 'Line'),'LineWidth',3);


    hold on
    shadedErrorBar(PSTH_Time_TgtOn_flexible,PSTH_TgtOn_NonPrefer_RF_flexible_mean(i,:),PSTH_TgtOn_NonPrefer_RF_flexible_sem(i,:),'lineprops',{[101,175,236]/255},'transparent',1,'patchSaturation',0.3);
    ylim([-1,5]);

    set(findall(gca, 'Type', 'Line'),'LineWidth',3);
    xlabel('Time from Tgt Onset(ms)');
    ylabel('Z-Score');
    set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
    

    subplot(length(UniqueMol),2,i*2)
    shadedErrorBar(PSTH_Time_SacOn_flexible,PSTH_SacOn_Prefer_RF_flexible_mean(i,:),PSTH_SacOn_Prefer_RF_flexible_sem(i,:),'lineprops',{[255,194,10]/255},'transparent',1,'patchSaturation',0.3);
    hold on
    shadedErrorBar(PSTH_Time_SacOn_flexible,PSTH_SacOn_NonPrefer_RF_flexible_mean(i,:),PSTH_SacOn_NonPrefer_RF_flexible_sem(i,:),'lineprops',{[101,175,236]/255},'transparent',1,'patchSaturation',0.3);
    ylim([-1,5]);
   set(findall(gca, 'Type', 'Line'),'LineWidth',3);

   xlabel('Time from Sac Onset(ms)');
   ylabel('Z-Score');
    set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');

    end 


     FigureStartNum = FigureStartNum+1;
        FigureIndex = FigureIndex+1;


        figtitlestr{FigureIndex}='ValuePreferenceSummary';
        fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[50,100, 500,600],'Name',figtitlestr{FigureIndex});

        subplot(2,1,1)
        b = bar(PropValueCell);
        b.FaceColor = 'flat';
        b.CData = ColorOrder;

        xticklabels({'ExFast','ExSlow','In','NoEffect'});
        ylabel('Proportion of value selective cells');
        box off;
        set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
        


        subplot(2,1,2)
        b=bar(GoodPref_Vis_RF_Prop);
        b.FaceColor = 'flat';
        b.CData = ColorOrder;
        xticklabels({'ExFast','ExSlow','In','NoEffect'});
        ylabel('Proportion of good preference cells');
        box off;
        hold on 
        plot([0,5],[0.5,0.5],'--k','LineWidth',3);
        set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
        




  end %End of if plot flexible value



if PlotOneDR
   %  ValueSummary_OneDR(i,:)=[ValueCellCount_oneDR(i),NonValueCellCount_oneDR(i),ValueRF_Count(i),Value_nonRF_Count(i),Value_Both_Count(i),ValueConsistent_Count(i),ValueInConsistent_Count(i)];
%% Counts
 FigureStartNum = FigureStartNum+1;
 FigureIndex = FigureIndex+1;


 figtitlestr{FigureIndex}='ValueInfo_OneDR';
 fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[50,100, 500,600],'Name',figtitlestr{FigureIndex});

 for i = 1:length(UniqueMol)

            subplot(3,4,i);
             b = bar([1,2],ValueSummary_OneDR(i,1:2));
            b.FaceColor='flat';
            b.CData = ColorOrder(i,:);
          % pie(ValueSummary_OneDR(i,1:2))
            xticklabels({'Value','NonValue'});
          % colormap([1,0.8,0;1,1,1]);
           %legend({'Value','NonValue'});

            subplot(3,4,4+i);
            b = bar([1,2,3],ValueSummary_OneDR(i,3:5));
            b.FaceColor='flat';
            b.CData = ColorOrder(i,:);
            xticklabels({'ValueContraOnly','ValueIspiOnly','ValueBoth'});

            subplot(3,4,4*2+i);
            b = bar([1,2],ValueSummary_OneDR(i,6:7));
            b.FaceColor='flat';
            b.CData = ColorOrder(i,:);
            xticklabels({'Contra-Ipsi Consistent','Contra-Ipsi Inconsistent'});

 end

 %% Counts Few
 FigureStartNum = FigureStartNum+1;
 FigureIndex = FigureIndex+1;


 figtitlestr{FigureIndex}='ValueInfo_OneDR2Group';
        fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[50,100, 500,600],'Name',figtitlestr{FigureIndex});
        
% Value Proportion
        FastEx = ValueSummary_OneDR(1,1:2);
        Others = ValueSummary_OneDR(2,1:2)+ValueSummary_OneDR(3,1:2)+ValueSummary_OneDR(4,1:2);
            subplot(4,1,1);
            %{
            b = bar([1,2],FastEx);
            b.FaceColor='flat';
            b.CData = ColorOrder(1,:);
            xticklabels({'Value','NonValue'});
            ylabel('Number of Neurons');
            box off
            %}
            pie(FastEx);
            legend({'Value','NonValue'});
             set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');


          %  subplot(3,2,2);
            
            %{
            b = bar([1,2],Others);
            b.FaceColor='flat';
            b.CData = ColorOrder(4,:);
            xticklabels({'Value','NonValue'});
            
            

            ylabel('Number of Neurons');
            box off

             set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');

            %}

        % RF vs NonRF
        FastEx = ValueSummary_OneDR(1,3:5);
        Others = ValueSummary_OneDR(2,3:5)+ValueSummary_OneDR(3,3:5)+ValueSummary_OneDR(4,3:5);
        subplot(4,1,2);
        
        pie(FastEx);
            legend({'Value','NonValue'});
           %{ 
             b = bar([1,2,3],FastEx);
            b.FaceColor='flat';
            b.CData = ColorOrder(1,:);
            xticklabels({'ValueContraOnly','ValueIspiOnly','ValueBoth'});
            ylabel('Number of Neurons');
            box off

             set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
           %}
%{
           subplot(4,1,4);
            b = bar([1,2,3],Others);
            b.FaceColor='flat';
            b.CData = ColorOrder(4,:);
            xticklabels({'ValueContraOnly','ValueIspiOnly','ValueBoth'});
            ylabel('Number of Neurons');
            box off

             set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
%}
        
           %Good vs Bad

           FastEx = ValueSummary_OneDR(1,6:7);
        Others = ValueSummary_OneDR(2,6:7)+ValueSummary_OneDR(3,6:7)+ValueSummary_OneDR(4,6:7);

           subplot(4,1,3);
           pie(FastEx);
           legend({'Contra-Ipsi Consistent','Contra-Ipsi Inconsistent'});
           %{
            b = bar([1,2],FastEx);
            b.FaceColor='flat';
            b.CData = ColorOrder(1,:);
             xticklabels({'Contra-Ipsi Consistent','Contra-Ipsi Inconsistent'});
            ylabel('Number of Neurons');
            box off
           %}
             set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');

%{
             subplot(4,1,5);
             pi
            b = bar([1,2],Others);
            b.FaceColor='flat';
            b.CData = ColorOrder(4,:);
            xticklabels({'Contra-Ipsi Consistent','Contra-Ipsi Inconsistent'});

            ylabel('Number of Neurons');
            box off

             set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
%}
        FastEx = ValueSummary_OneDR(1,8:9);
        Others = ValueSummary_OneDR(2,8:9)+ValueSummary_OneDR(3,8:9)+ValueSummary_OneDR(4,8:9);

             subplot(4,1,4);
             pie( FastEx );
             legend({'Good','Bad'});



           %  BadPreferCount(i)



%% 
        FigureStartNum = FigureStartNum+1;
        FigureIndex = FigureIndex+1;


        figtitlestr{FigureIndex}='PSTH_Contralateral_OneDR';
        fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[50,100, 500,600],'Name',figtitlestr{FigureIndex});
     

  
    for i = 1:length(UniqueMol)
    subplot(length(UniqueMol),2,(i-1)*2+1)
     shadedErrorBar(PSTH_Time_TgtOn_OneDR,PSTH_TgtOn_Prefer_RF_mean(i,:),PSTH_TgtOn_Prefer_RF_sem(i,:),'lineprops',{[255,194,10]/255},'transparent',1,'patchSaturation',0.3)
     set(findall(gca, 'Type', 'Line'),'LineWidth',3);


    hold on
    shadedErrorBar(PSTH_Time_TgtOn_OneDR,PSTH_TgtOn_NonPrefer_RF_mean(i,:),PSTH_TgtOn_NonPrefer_RF_sem(i,:),'lineprops',{[101,175,236]/255},'transparent',1,'patchSaturation',0.3);
    

    set(findall(gca, 'Type', 'Line'),'LineWidth',3);
    xlabel('Time from Tgt Onset(ms)');
    ylabel('Z-Score');
    set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
    

    subplot(length(UniqueMol),2,i*2)
    shadedErrorBar(PSTH_Time_SacOn_OneDR,PSTH_SacOn_Prefer_RF_mean(i,:),PSTH_SacOn_Prefer_RF_sem(i,:),'lineprops',{[255,194,10]/255},'transparent',1,'patchSaturation',0.3);
    hold on
    shadedErrorBar(PSTH_Time_SacOn_OneDR,PSTH_SacOn_NonPrefer_RF_mean(i,:),PSTH_SacOn_NonPrefer_RF_sem(i,:),'lineprops',{[101,175,236]/255},'transparent',1,'patchSaturation',0.3);
    
   set(findall(gca, 'Type', 'Line'),'LineWidth',3);

   xlabel('Time from Sac Onset(ms)');
   ylabel('Z-Score');
    set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');

    end 

%
     FigureStartNum = FigureStartNum+1;
        FigureIndex = FigureIndex+1;


        figtitlestr{FigureIndex}='ValuePreferenceSummary';
        fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[50,100, 500,600],'Name',figtitlestr{FigureIndex});

        subplot(2,1,1)
        b = bar(PropValueCell_OneDR);
        b.FaceColor = 'flat';
        b.CData = ColorOrder;

        xticklabels({'ExFast','ExSlow','In','NoEffect'});
        ylabel('Proportion of value selective cells');
        box off;
        set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
        

%
        subplot(2,1,2)
        b=bar( GoodPref_Vis_RF_Prop_OneDR);
        b.FaceColor = 'flat';
        b.CData = ColorOrder;
        xticklabels({'ExFast','ExSlow','In','NoEffect'});
        ylabel('Proportion of good preference cells');
        box off;
        hold on 
        plot([0,5],[0.5,0.5],'--k','LineWidth',3);
        set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
        
%}


end %End of if plot oneDR




 


end %End of if figureplot

if OutputFlag == 1

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
OutputFigureName='SNr_VideoView_Neural_';


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

end %End of if outputFlag == 1

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

