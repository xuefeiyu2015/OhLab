function Batch_Opto_WaveFormClustering(Data,  StartFile, EndFile,ShowFigureFlag,OutputFlag);
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




    if isfield(OutputData,'Waveform')
      
       Waveform_tmp = OutputData.Waveform;
       wf = Waveform_tmp.DataStamp;

       

       if ~isempty(wf('AverageWaveForm')) & ~isempty(wf('TimePoint'))
        wf_ave{index} = wf('AverageWaveForm');
        
        wf_time{index}=wf('TimePoint');
        
        if length(wf('AverageWaveForm'))==71
        %   wf_index(index) = 1;
         wf_type(index) = 1;% presumed to be Plexon
        else
            wf_type(index) = 2;% presumed to be Blip
        end
       else
           wf_ave{index}=NaN*ones(1,180);
        wf_time{index}=NaN*ones(1,180);
        wf_lack{index}=FilesName{i}{j};
       % wf_index(index) = 0;
       wf_type(index)=0;
       end

    else
        wf_ave{index}=NaN*ones(1,180);
        wf_time{index}=NaN*ones(1,180);
        wf_lack{index}=FilesName{i}{j};
       % wf_index(index) = 0;
       wf_type(index)=0;


    end %End of if there is field waveform

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
PlotMemorySaccade = 1;
PlotDelaySaccade = 0;
PlotStableValue = 0;
PlotFlexibleValue= 0;
PlotOneDR = 0;

%Analysis controller 
AnalysisOS = 0;% For Video Free View
AnalysisMS = 1;% For memory saccade
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

PSTH_norm = cell2mat(PSTH_Stim_Norm(ScreenCell)'); %in use
PSTH_control_norm = cell2mat(PSTH_Control_Norm(ScreenCell)');


PSTH_norm2 =cell2mat(PSTH_Stim_Norm2(ScreenCell)');



PSTH_z = cell2mat(PSTH_Stim_Mean_z(ScreenCell)');


p_whole_cri50 = cell2mat(p_whole_cri50);

wf_type=wf_type(ScreenCell)';
wf_ave =wf_ave(ScreenCell)';
wf_time= wf_time(ScreenCell)';



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

%{
%Find examples
%figure
Inhi = PSTH_Use(Inhibition,:);

FilesNameAll_inh = FilesNameAll(Inhibition,:);
ChannelScreened_inh = ChannelScreened(Inhibition);

total_num = size(FilesNameAll_inh ,1);
%plot(Inhi')
%
figure_index = 0;
sub_index = 1;
for i = 1:total_num
    if mod(i,16)==0 || i == 1
        figure_index = figure_index+1;
        figure(figure_index);
        sub_index = 1;
    end
    subplot(4,4,sub_index);
    plot(PSTH_Opto_Time,Inhi(i,:),'-b');
    title(i);
    sub_index=sub_index+1;


end
%}
%{
keyboard
Example1 = [289,296,291,244,247,203,192,162,89];
export_path = FilesNameAll_inh(Example1);
A=export_path;
B=ChannelScreened_inh(Example1)';
keyboard
%}
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



wf_type_fast =  wf_type(Ex_Fast_Sig& wf_type==2);
wf_ave_fast = wf_ave(Ex_Fast_Sig & wf_type==2);
wf_time_fast = wf_time(Ex_Fast_Sig& wf_type==2);

wf_type_slow =  wf_type(Ex_Slow_Sig& wf_type==2);
wf_ave_slow = wf_ave(Ex_Slow_Sig & wf_type==2);
wf_time_slow = wf_time(Ex_Slow_Sig& wf_type==2);

wf_type_inh =  wf_type(In_Sig& wf_type==2);
wf_ave_inh = wf_ave(In_Sig& wf_type==2);
wf_time_inh = wf_time(In_Sig& wf_type==2);

wf_ave_fast_mat = cell2mat(cellfun(@(x) x(1:71),wf_ave_fast,'UniformOutput',0));
wf_ave_slow_mat = cell2mat(cellfun(@(x) x(1:71),wf_ave_slow,'UniformOutput',0));
wf_ave_inh_mat = cell2mat(cellfun(@(x) x(1:71),wf_ave_inh,'UniformOutput',0));


wf_ave_fast_mat = wf_ave_fast_mat./repmat(max(wf_ave_fast_mat,[],2),1,size(wf_ave_fast_mat,2));
wf_ave_slow_mat = wf_ave_slow_mat./repmat(max(wf_ave_slow_mat,[],2),1,size(wf_ave_slow_mat,2));
wf_ave_inh_mat = wf_ave_inh_mat./repmat(max(wf_ave_inh_mat,[],2),1,size(wf_ave_inh_mat,2));


wf_whole = [wf_ave_fast_mat;wf_ave_slow_mat;wf_ave_inh_mat];
[coeff, score, ~] = pca(wf_whole);
pca_data = score(:, 1:3);
colors = [1,0,0;1.0000,0.5490,0;0,0,1];
labels = zeros(size(wf_whole,1),1);
labels(1:size(wf_ave_fast_mat,1),1)=1;
labels(size(wf_ave_fast_mat,1)+1:size(wf_ave_fast_mat,1)+size(wf_ave_slow_mat,1),1)=2;
labels(size(wf_ave_fast_mat,1)+size(wf_ave_slow_mat,1)+1:size(wf_whole,1),1)=3;

figure
subplot(2,2,1)
plot(wf_ave_fast_mat','-r');
subplot(2,2,2)
plot(wf_ave_slow_mat','-y');
subplot(2,2,3)
plot(wf_ave_inh_mat','-b');
subplot(2,2,4)
hold on;
for i = 1:3
    scatter3(pca_data(labels == i, 1), pca_data(labels == i, 2), pca_data(labels == i, 3), ...
        50, colors(i, :), 'filled');
end
view(3);
xlabel('PC1');
ylabel('PC2');
zlabel('PC3');
title('Waveform Clustering in PCA Space');
legend('Fast', 'Slow', 'Inh');

keyboard

% wf_type(index) 


%
%Find examples
%{
example_fast_saccade=FilesNameAll(Ex_Fast_Sig & CellType==3);
Channel_sel = ChannelScreened(Ex_Fast_Sig & CellType==3)';
A=Channel_sel ;
keyboard
%}

end