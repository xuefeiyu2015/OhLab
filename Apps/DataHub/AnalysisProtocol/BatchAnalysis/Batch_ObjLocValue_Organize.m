function Batch_ObjLocValue_Organize(Data,  StartFile, EndFile,ShowFigureFlag,OutputFlag)
BatchFileName=Data.BatchFileName;
%Batch data file path
FilesName=Data.ResultFilePath(StartFile: EndFile);
RecordDate=Data.RecordDate(StartFile: EndFile);
FileNum=EndFile-StartFile+1;

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
    % tmp = load(FilesName{i}{j});
    % OutputData = tmp.OutputData;
     
     %OutputDataTemp=OutputData;

     FilesNameAll{index}=FilesName_Each{i};
     ChannelNumAll(index) = ChannelID(i, 1);

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


      %% %Stable Value
    % if isfield(OutputData,'StableValue')

          StableValueIndex(index) = 1;

          Data =OutputData.StableValue(1).DataStamp;

%{
         p_Tgt_RF_stable(index) = Data('p_Tgt_RF');
         if isempty(Data('p_Tgt_nonRF'))
              p_Tgt_nonRF_stable(index) = NaN;
         else
             
             p_Tgt_nonRF_stable(index) = Data('p_Tgt_nonRF');
         end

       %  Tgt_Good_RF_stable(index) = Data('Tgt_Good_RF');


      %  PreTgtGood_RF(index) = Data('PreTgtGood_RF');
      %  PreTgtBad_RF(index) = Data('PreTgtBad_RF');
        p_Sac_RF_stable(index) = Data('p_Sac_RF');
        if isempty(Data('p_Sac_nonRF'))
            p_Sac_nonRF_stable(index)=NaN;

        else
            p_Sac_nonRF_stable(index) = Data('p_Sac_nonRF');
        end
%}
        
%{
        PSTH_TgtOn_Prefer_RF_stable{index} = Data('PSTH_TgtOn_Prefer_RF_z');
        PSTH_TgtOn_nonPrefer_RF_stable{index} = Data('PSTH_TgtOn_nonPrefer_RF_z');
        PSTH_TgtOn_Prefer_nonRF_stable{index} = Data('PSTH_TgtOn_Prefer_nonRF_z');
        PSTH_TgtOn_nonPrefer_nonRF_stable{index} = Data('PSTH_TgtOn_nonPrefer_nonRF_z');

        PSTH_SacOn_Prefer_RF_stable{index} = Data('PSTH_SacOn_Prefer_RF_z');
        PSTH_SacOn_nonPrefer_RF_stable{index} = Data('PSTH_SacOn_nonPrefer_RF_z');
        PSTH_SacOn_Prefer_nonRF_stable{index} = Data('PSTH_SacOn_Prefer_nonRF_z');
        PSTH_SacOn_nonPrefer_nonRF_stable{index} = Data('PSTH_SacOn_nonPrefer_nonRF_z');

%}

        PSTH_TgtOn_Good_RF_stable{index} = Data('PSTH_TgtOn_Good_RF_z');
        PSTH_TgtOn_Bad_RF_stable{index} = Data('PSTH_TgtOn_Bad_RF_z');
        PSTH_TgtOn_Good_nonRF_stable{index} = Data('PSTH_TgtOn_Good_nonRF_z');
        PSTH_TgtOn_Bad_nonRF_stable{index} = Data('PSTH_TgtOn_Bad_nonRF_z');

        PSTH_SacOn_Good_RF_stable{index} = Data('PSTH_SacOn_Good_RF_z');
        PSTH_SacOn_Bad_RF_stable{index} = Data('PSTH_SacOn_Bad_RF_z');
        PSTH_SacOn_Good_nonRF_stable{index} = Data('PSTH_SacOn_Good_nonRF_z');
        PSTH_SacOn_Bad_nonRF_stable{index} = Data('PSTH_SacOn_Bad_nonRF_z');



        psth_rawtgt_good_RF = Data('PSTH_TgtOn_Good_RF');
        psth_rawtgt_bad_RF = Data('PSTH_TgtOn_Bad_RF');

        psth_rawsac_good_RF = Data('PSTH_SacOn_Good_RF');
        psth_rawsac_bad_RF = Data('PSTH_SacOn_Bad_RF');

        psth_rawtgt_good_NonRF = Data('PSTH_TgtOn_Good_nonRF');
        psth_rawtgt_bad_NonRF = Data('PSTH_TgtOn_Bad_nonRF');

        psth_rawsac_good_NonRF = Data('PSTH_SacOn_Good_nonRF');
        psth_rawsac_bad_NonRF = Data('PSTH_SacOn_Bad_nonRF');

        FR_Mean_stable(index) = mean([nanmean(psth_rawtgt_good_RF),nanmean(psth_rawtgt_bad_RF),nanmean(psth_rawsac_good_RF),nanmean(psth_rawsac_bad_RF),nanmean(psth_rawtgt_good_NonRF),nanmean(psth_rawtgt_bad_NonRF),nanmean(psth_rawsac_good_NonRF),nanmean(psth_rawsac_bad_NonRF)]);


        PSTH_Time_TgtOn_stable = Data('Time_TgtOn');
        PSTH_Time_SacOn_stable = Data('Time_SacOn');

       

        AUCPSTH_TgtOn_Value_RF_stable{index} = Data('AUCPSTH_TgtOn_Value_RF');
        
        AUCPSTH_TgtOn_Value_nonRF_stable{index} = Data('AUCPSTH_TgtOn_Value_nonRF');
        

        AUCPSTH_SacOn_Value_RF_stable{index} = Data('AUCPSTH_SacOn_Value_RF');
        
        AUCPSTH_SacOn_Values_nonRF_stable{index} = Data('AUCPSTH_SacOn_Value_nonRF');
        



        p_Vis_RF(index) = Data('p_Vis_RF');
        VisDiff_RF(index) = Data('VisDiff_RF');
        p_Vis_nonRF(index) = Data('p_Vis_nonRF');
        VisDiff_nonRF(index) = Data('VisDiff_nonRF');

        p_pureSac_RF(index) = Data('p_pureSac_RF');
        pureSacDiff_RF(index) = Data('pureSacDiff_RF');
        p_pureSac_nonRF(index) = Data('p_pureSac_nonRF');
        pureSacDiff_nonRF(index) = Data('pureSacDiff_nonRF');

        hPSTH_TgtOn_Value_RF{index}=Data('hPSTH_TgtOn_Value_RF');

      
       %peakAUC_Stable_RF(index) = Data('peakAUC_avg_RF');
       %p_peakAUC_Stable_RF(index) = Data('p_peakAUC_RF');

       %peakTime_AUC_Stable_RF(index) = Data('peakTime_AUC_RF');
       LatencyTimeAUC_Stable_RF(index) = Data('LatencyTimeAUC_RF');
       LatencyTimeAUC_Stable_ranksumRF(index) = Data('LatencyTimeAUCranksum_RF');

       TgtAUC_RF(index) = Data('TgtAUC_RF');
       p_AUC_Tgt_RF(index) = Data('p_AUC_Tgt_RF');
       p_AUC_Tgt_RF_ranksum(index) = Data('p_AUC_Tgt_RF_ranksum');
       TgtAUC_nonRF(index) = Data('TgtAUC_nonRF');
       p_AUC_Tgt_nonRF(index) = Data('p_AUC_Tgt_nonRF');
       p_AUC_Tgt_nonRF_ranksum(index) = Data('p_AUC_Tgt_nonRF_ranksum');
%
       psth_targeton_raw_stable{index} = Data('psth_targeton_raw');
       target_index_stable{index}      = Data('target_index');
       RF_index_stable{index}          = Data('RF_index');
       n_high_RF_stable(index)         = Data('n_high_RF');
       n_low_RF_stable(index)          = Data('n_low_RF');
       n_high_nonRF_stable(index)      = Data('n_high_nonRF');
       n_low_nonRF_stable(index)       = Data('n_low_nonRF');
%}



%end % Analysis all the data according to stable value

%Free Viewing Task

if isfield(OutputData,'FreeViewNeural')  & ~ismember(FilesNameAll(index),'Robin060923FRACSAC1')
    clear Data;
      Data =OutputData.FreeViewNeural.DataStamp;
      OptoIndex(index) = 1;
%

    PSTH_Stim_Norm{index} = Data('PSTH_Stim_Mean_norm');


     PSTH_Control_Norm{index} = Data('PSTH_Control_Mean_norm');
     PSTH_Stim_Mean_z{index} = Data('PSTH_Stim_Mean_z');

   %   PSTH_Stim_Mean_z{index} = Data('PSTH_Stim_Mean_z');%Use the z scored ones
   %  PSTH_Control_Mean_z{index} = Data('PSTH_Control_Mean_z');%Use the z scored ones

  %    PSTH_Stim_Norm{index} = Data('PSTH_Stim_Mean_z');
  %   PSTH_Control_Norm{index} = Data('PSTH_Control_Mean_z');

      PSTH_Stim_Mean{index} = Data('PSTH_Stim_Mean');
     PSTH_Control_Mean{index} = Data('PSTH_Control_Mean');

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
   
     PSTH_First{index} = Data('PSTH_Stim_Mean_First_z');
     PSTH_Second{index} = Data('PSTH_Stim_Mean_Second_z');

     p_whole_cri50{index} = Data('p_whole_cri50');



     


     AS_First(index) = Data('AS_First');
     AS_Second(index) = Data('AS_Second');

     PSTH_Opto_Time = Data('PSTH_Time');
else
     OptoIndex(index) = 0;
     PSTH_Stim_Norm{index} = NaN;
     PSTH_Control_Norm{index} = NaN;
     PSTH_Stim_Mean_z{index} = NaN;

   %   PSTH_Stim_Mean_z{index} = Data('PSTH_Stim_Mean_z');%Use the z scored ones
   %  PSTH_Control_Mean_z{index} = Data('PSTH_Control_Mean_z');%Use the z scored ones

  %    PSTH_Stim_Norm{index} = Data('PSTH_Stim_Mean_z');
  %   PSTH_Control_Norm{index} = Data('PSTH_Control_Mean_z');

      PSTH_Stim_Mean{index} = NaN;
     PSTH_Control_Mean{index} =NaN;

     FR_Mean_Video(index) = NaN;

     p_StimControl(index) = NaN;

     p_Cri20(index) = NaN;
     p_Cri50(index) = NaN;
     Latency2(index) = NaN;

     PSTH_Stim_Norm2{index} = NaN;


     
     AS(index) = NaN;
     Latency(index) = NaN;
%     p_AS(index) = Data('Ac_p');
     PMS(index) = NaN;
     PMS_p(index) = NaN;
     SM(index) = NaN;
     SM_p(index) = NaN;
   
     PSTH_First{index} = NaN;
     PSTH_Second{index} = NaN;

     p_whole_cri50{index} = NaN;




     


     AS_First(index) = NaN;
     AS_Second(index) = NaN;

  


      
end

%% For one DR task

if isfield(OutputData,'OneDR')

    OneDRIndex(index) = 1;
clear Data;
     Data =OutputData.OneDR.DataStamp;

    
    
%

      PSTH_Good_TgtOn_RF_oneDR{index} = Data('PSTH_Good_TgtOn_RF_z');
     PSTH_Bad_TgtOn_RF_oneDR{index} = Data('PSTH_Bad_TgtOn_RF_z');
     PSTH_Good_TgtOn_nonRF_oneDR{index} = Data('PSTH_Good_TgtOn_nonRF_z');
     PSTH_Bad_TgtOn_nonRF_oneDR{index} = Data('PSTH_Bad_TgtOn_nonRF_z');

     PSTH_SacGood_nonRF_oneDR{index} = Data('PSTH_SacGood_nonRF_z');
     PSTH_SacBad_nonRF_oneDR{index} = Data('PSTH_SacBad_nonRF_z');
     PSTH_SacGood_RF_oneDR{index} = Data('PSTH_SacGood_RF_z');
     PSTH_SacBad_RF_oneDR{index} = Data('PSTH_SacBad_RF_z');
%}
        psth_rawtgt_good_RF = Data('PSTH_Good_TgtOn_RF');
        psth_rawtgt_bad_RF = Data('PSTH_Bad_TgtOn_RF');

        psth_rawsac_good_RF = Data('PSTH_SacGood_RF');
        psth_rawsac_bad_RF = Data('PSTH_SacBad_RF');

        psth_rawtgt_good_NonRF = Data('PSTH_Good_TgtOn_nonRF');
        psth_rawtgt_bad_NonRF = Data('PSTH_Bad_TgtOn_nonRF');

        psth_rawsac_good_NonRF = Data('PSTH_SacGood_nonRF');
        psth_rawsac_bad_NonRF = Data('PSTH_SacBad_nonRF');

        FR_Mean_oneDR(index) = mean([nanmean(psth_rawtgt_good_RF),nanmean(psth_rawtgt_bad_RF),nanmean(psth_rawsac_good_RF),nanmean(psth_rawsac_bad_RF),nanmean(psth_rawtgt_good_NonRF),nanmean(psth_rawtgt_bad_NonRF),nanmean(psth_rawsac_good_NonRF),nanmean(psth_rawsac_bad_NonRF)]);




     PSTH_Time_TgtOn_OneDR = Data('Time_TgtOn');
     PSTH_Time_SacOn_OneDR = Data('Time_SacOn');

   
  
        AUCPSTH_TgtOn_Value_RF_oneDR{index} = Data('AUCPSTH_TgtOn_Value_RF');
        
        AUCPSTH_TgtOn_Value_nonRF_oneDR{index} = Data('AUCPSTH_TgtOn_Value_nonRF');
        

        AUCPSTH_SacOn_Value_RF_oneDR{index} = Data('AUCPSTH_SacOn_Value_RF');
        
       % AUCPSTH_SacOn_Values_nonRF_oneDR{index} = Data('AUCPSTH_SacOn_Value_nonRF');
        



        p_Vis_RF_oneDR(index) = Data('p_Vis_RF');
        VisDiff_RF_oneDR(index) = Data('VisDiff_RF');
        p_Vis_nonRF_oneDR(index) = Data('p_Vis_nonRF');
        VisDiff_nonRF_oneDR(index) = Data('VisDiff_nonRF');

        p_pureSac_RF_oneDR(index) = Data('p_pureSac_RF');
        pureSacDiff_RF_oneDR(index) = Data('pureSacDiff_RF');
        p_pureSac_nonRF_oneDR(index) = Data('p_pureSac_nonRF');
        pureSacDiff_nonRF_oneDR(index) = Data('pureSacDiff_nonRF');

        hPSTH_TgtOn_Value_RF_oneDR{index}=Data('hPSTH_TgtOn_Value_RF');

      
       %peakAUC_OneDR_RF(index) = Data('peakAUC_avg_RF');
       %p_peakAUC_OneDR_RF(index) = Data('p_peakAUC_RF');

       %peakTime_AUC_OneDR_RF(index) = Data('peakTime_AUC_RF');
       LatencyTimeAUC_OneDR_RF(index) = Data('LatencyTimeAUC_RF');
       LatencyTimeAUC_OneDR_ranksumRF(index) = Data('LatencyTimeAUC_ranksum_RF');

%
       psth_targeton_raw_oneDR{index} = Data('psth_targeton_raw');
       target_index_oneDR{index}      = Data('target_index');
       RF_index_oneDR{index}          = Data('RF_index');
       n_high_RF_oneDR(index)         = Data('n_high_RF');
       n_low_RF_oneDR(index)          = Data('n_low_RF');
       n_high_nonRF_oneDR(index)      = Data('n_high_nonRF');
       n_low_nonRF_oneDR(index)       = Data('n_low_nonRF');

       TgtAUC_oneDR_RF(index) = Data('TgtAUC_RF');
       p_AUC_Tgt_OneDR_RF(index) = Data('p_AUC_Tgt_RF');
       p_AUC_Tgt_OneDR_RF_ranksum(index) = Data('p_AUC_Tgt_RF_ranksum');
       TgtAUC_OneDR_nonRF(index) = Data('TgtAUC_nonRF');
       p_AUC_Tgt_OneDR_nonRF(index) = Data('p_AUC_Tgt_nonRF');
      p_AUC_Tgt_OneDR_nonRF_ranksum(index) = Data('p_AUC_Tgt_nonRF_ranksum');

%}


else
        OneDRIndex(index) = 0;

      

     PSTH_Good_TgtOn_RF_oneDR{index} =NaN;
     PSTH_Bad_TgtOn_RF_oneDR{index} = NaN;
     PSTH_Good_TgtOn_nonRF_oneDR{index} = NaN;
     PSTH_Bad_TgtOn_nonRF_oneDR{index} = NaN;

     PSTH_SacGood_nonRF_oneDR{index} = NaN;
     PSTH_SacBad_nonRF_oneDR{index} = NaN;
     PSTH_SacGood_RF_oneDR{index} = NaN;
     PSTH_SacBad_RF_oneDR{index} = NaN;
%}
        psth_rawtgt_good_RF = NaN;
        psth_rawtgt_bad_RF = NaN;

        psth_rawsac_good_RF = NaN;
        psth_rawsac_bad_RF = NaN;

        psth_rawtgt_good_NonRF = NaN;
        psth_rawtgt_bad_NonRF = NaN;

        psth_rawsac_good_NonRF = NaN;
        psth_rawsac_bad_NonRF = NaN;

        FR_Mean_oneDR(index) = NaN;



     

   
  
        AUCPSTH_TgtOn_Value_RF_oneDR{index} = NaN;
        
        AUCPSTH_TgtOn_Value_nonRF_oneDR{index} = NaN;
        

        AUCPSTH_SacOn_Value_RF_oneDR{index} = NaN;
        
       % AUCPSTH_SacOn_Values_nonRF_oneDR{index} = Data('AUCPSTH_SacOn_Value_nonRF');
        



        p_Vis_RF(index) = NaN;
        VisDiff_RF(index) = NaN;
        p_Vis_nonRF(index) = NaN;
        VisDiff_nonRF(index) = NaN;

        p_pureSac_RF(index) = NaN;
        pureSacDiff_RF(index) = NaN;
        p_pureSac_nonRF(index) = NaN;
        pureSacDiff_nonRF(index) =NaN;

        hPSTH_TgtOn_Value_RF{index}=NaN;

      
       %peakAUC_OneDR_RF(index) =NaN;
       %p_peakAUC_OneDR_RF(index) = NaN;

       %peakTime_AUC_OneDR_RF = NaN;
       LatencyTimeAUC_OneDR_RF(index) = NaN;
       LatencyTimeAUC_OneDR_ranksumRF(index) = NaN;
  


     

    
             

              p_Vis_RF_oneDR(index) = NaN;
        VisDiff_RF_oneDR(index) = NaN;
        p_Vis_nonRF_oneDR(index) = NaN;
        VisDiff_nonRF_oneDR(index) = NaN;

        p_pureSac_RF_oneDR(index) = NaN;
        pureSacDiff_RF_oneDR(index) = NaN;
        p_pureSac_nonRF_oneDR(index) = NaN;
        pureSacDiff_nonRF_oneDR(index) = NaN;

        hPSTH_TgtOn_Value_RF_oneDR{index}=NaN;
        hPSTH_TgtOn_Good_RF_oneDR{index}=NaN;
        hPSTH_TgtOn_Bad_RF_oneDR{index}=NaN;
        

        hPSTH_SacOn_Value_RF_oneDR{index}=NaN;
        hPSTH_SacOn_Good_RF_oneDR{index}=NaN;
        hPSTH_SacOn_Bad_RF_oneDR{index}=NaN;

        
        hPSTH_TgtOn_Value_nonRF_oneDR{index}=NaN;
        hPSTH_TgtOn_Good_nonRF_oneDR{index}=NaN;
        hPSTH_TgtOn_Bad_nonRF_oneDR{index}=NaN;

        hPSTH_SacOn_Value_nonRF_oneDR{index}=NaN;
        hPSTH_SacOn_Good_nonRF_oneDR{index}=NaN;
        hPSTH_SacOn_Bad_nonRF_oneDR{index}=NaN;
        


     
       
     

       Value_Vis_RF_OneDR(index) = NaN;

    %    Value_Vis_nonRF_OneDR(index) = sign(Tgt_Good_nonRF(index)-Tgt_Bad_nonRF(index));
    Value_Vis_nonRF_OneDR(index) =NaN;

        Value_Sac_RF_OneDR(index) = NaN;

     %   Value_Sac_nonRF_OneDR(index) = sign(Sac_Good_nonRF(index)-Sac_Bad_nonRF(index));
     Value_Sac_nonRF_OneDR(index) = NaN;


        %peakAUC_OneDR_RF(index) = NaN;
       %p_peakAUC_OneDR_RF(index) = NaN;

       %peakTime_AUC_OneDR_RF(index) = NaN;
       LatencyTimeAUC_OneDR_RF(index) = NaN;

       psth_targeton_raw_oneDR{index} = NaN;
       target_index_oneDR{index}      = NaN;
       RF_index_oneDR{index}          = NaN;
       n_high_RF_oneDR(index)         = NaN;
       n_low_RF_oneDR(index)          = NaN;
       n_high_nonRF_oneDR(index)      = NaN;
       n_low_nonRF_oneDR(index)       = NaN;


        TgtAUC_oneDR_RF(index) = NaN;
       p_AUC_Tgt_OneDR_RF(index) = NaN;
       p_AUC_Tgt_OneDR_RF_ranksum(index) = NaN;
       TgtAUC_OneDR_nonRF(index) = NaN;
       p_AUC_Tgt_OneDR_nonRF(index) = NaN;
       p_AUC_Tgt_OneDR_nonRF_ranksum(index) = NaN;




end

%% For passive viewing task 

if isfield(OutputData,'ObjectPassiveViewing')
    PV_Index(index)=1;
    clear Data;
    Data =OutputData.ObjectPassiveViewing.DataStamp;

    hPSTH_TgtOn_Value_RF_PV{index} =double(Data('hPSTH_TgtOn_Value_RF')); 
    hPSTH_TgtOn_Good_RF_PV{index} =double(Data('hPSTH_TgtOn_Good_RF')); 
    hPSTH_TgtOn_Bad_RF_PV{index} =double(Data('hPSTH_TgtOn_Bad_RF')); 
    Time_TgtOn_PV=Data('Time_TgtOn');

    %{
    PSTH_TgtOn_Good_RF_PV{index} =Data('PSTH_TgtOn_Good_RF_z'); 
    PSTH_TgtOn_Bad_RF_PV{index} =Data('PSTH_TgtOn_Bad_RF_z'); 
    %}
    PSTH_TgtOn_Good_RF_PV{index} =Data('PSTH_TgtOn_Good_RF'); 
    PSTH_TgtOn_Bad_RF_PV{index} =Data('PSTH_TgtOn_Bad_RF'); 

    AUCPSTH_TgtOn_Value_RF_PV{index} =Data('AUCPSTH_TgtOn_Value_RF'); 
   

    SelectTime_PV=[100,300];
    tmp = hPSTH_TgtOn_Value_RF_PV{index};

    hInterval = tmp(Time_TgtOn_PV>=SelectTime_PV(1) & Time_TgtOn_PV<=SelectTime_PV(2));
    if sum(hInterval)>5
        Value_PV_RF_Prop_Sig(index) = 1;
    else
        Value_PV_RF_Prop_Sig(index) = 0;

    end

  tmpGood = Data('PSTH_TgtOn_Good_RF_z');
  tmpBad = Data('PSTH_TgtOn_Bad_RF_z');

        Good_PV_RF = tmpGood(Time_TgtOn_PV>=100 & Time_TgtOn_PV<=300);
        Bad_PV_RF = tmpBad(Time_TgtOn_PV>=100 & Time_TgtOn_PV<=300);

        [h,p]=ttest2(Good_PV_RF,Bad_PV_RF);
        Value_PV_RF(index) = h;
        Value_PV_RF_Diff(index) = mean(Good_PV_RF)-mean(Bad_PV_RF);
        Value_PV_RF_ROC(index) = rocN(Good_PV_RF,Bad_PV_RF);

   
else
    PV_Index(index)=0;
    hPSTH_TgtOn_Value_RF_PV{index} = NaN*ones(1,160); 
    hPSTH_TgtOn_Good_RF_PV{index} = NaN*ones(1,160); 
    hPSTH_TgtOn_Bad_RF_PV{index} = NaN*ones(1,160); 

    PSTH_TgtOn_Good_RF_PV{index} = NaN*ones(1,160); 
    PSTH_TgtOn_Bad_RF_PV{index} = NaN*ones(1,160); 
    Value_PV_RF_Prop_Sig(index)=NaN;

    PSTH_TgtOn_Good_RF_PV{index} = NaN*ones(1,160); 
    PSTH_TgtOn_Bad_RF_PV{index} = NaN*ones(1,160); 

    AUCPSTH_TgtOn_Value_RF_PV{index} =NaN*ones(1,160); 


end



     
     index = index+1;   
    end %End of for file name
     
end


%%
%Separate cells into different modulation groups according to laser
%response


%Separate into 5 groups 


Fast_M = Latency<5;
Slow_M = Latency>=5;
Excitation = SM> 0;
Inhibition = SM<=0;

Sig_m = SM_p < 0.05&(~isnan(Latency));
%Sig_m = ~isnan(Latency) & AS >=5;

%Sig_m = (p_whole_cri50<0.05 | SM_p<0.05 | PMS_p<0.05)&(~isnan(Latency));
%Sig_m = p_whole_cri50>-0.1;
%Latency
Non_Sig_m = ~Sig_m;



%Non_Sig_m = Non_Sig_level2;


Ex_Fast_Sig = Fast_M & Excitation & Sig_m;

Ex_Slow_Sig = Slow_M & Excitation & Sig_m;



In_Sig = Inhibition & Sig_m;



ModuType = NaN*ones(1,length(Sig_m));



%%
ModuType(Ex_Fast_Sig) =1;
ModuType(Ex_Slow_Sig) =2;

ModuType(In_Sig) =3;
%ModuType(In_Slow_Sig) =4;

ModuType(Non_Sig_m) =4;




ModulString = {'Fast Excitation','Slow Excitation','Inhibition','No Effect'};

NonSigString = {'NonSig_Ex_l1','NonSig_In_l1','NonSig_Ex_l2','NonSig_In_l2'};

Latency_Ex = Latency(Ex_Fast_Sig|Ex_Slow_Sig);
Latency_In = Latency(In_Sig);



UniqueMol = unique(ModuType);

UniqueMol = UniqueMol(~isnan(UniqueMol));

ModuType=ModuType';%Transfer to column;



%Cell Selection
ResponseValid = (FR_Mean_oneDR>1 & FR_Mean_stable>1)';
VisualResponsive = p_Vis_RF<0.05 & VisDiff_RF>0;

Sel = OneDRIndex' == 1;
Sel_All = Sel & OptoIndex' ==1 & ResponseValid & VisualResponsive';
%Sel_All = OptoIndex' ==1;
Sel_Others = Sel_All & ModuType==4 ;
Sel_Fast = Sel_All & ModuType==1 ;

%%check psthAUC
 %AUCPSTH_TgtOn_Value_RF_stable{index} = Data('AUCPSTH_TgtOn_Value_RF');

ChannelFull_All = ChannelFull(Sel_All);
%% Population average
PSTH_TgtOn_Good_RF_Stable_All = cell2mat(PSTH_TgtOn_Good_RF_stable(Sel_All)');
PSTH_TgtOn_Bad_RF_Stable_All = cell2mat(PSTH_TgtOn_Bad_RF_stable(Sel_All)');

PSTH_TgtOn_Good_RF_oneDR_All = cell2mat(PSTH_Good_TgtOn_RF_oneDR(Sel_All)');
PSTH_TgtOn_Bad_RF_oneDR_All = cell2mat(PSTH_Bad_TgtOn_RF_oneDR(Sel_All)');



PSTH_TgtOn_Good_RF_OneDR_Others =cell2mat(PSTH_Good_TgtOn_RF_oneDR(Sel_Others)');
PSTH_TgtOn_Bad_RF_OneDR_Others=cell2mat(PSTH_Bad_TgtOn_RF_oneDR(Sel_Others)');

PSTH_TgtOn_Good_RF_OneDR_Fast =cell2mat(PSTH_Good_TgtOn_RF_oneDR(Sel_Fast)');
PSTH_TgtOn_Bad_RF_OneDR_Fast=cell2mat(PSTH_Bad_TgtOn_RF_oneDR(Sel_Fast)');




PSTH_TgtOn_Good_RF_Stable_Others = cell2mat(PSTH_TgtOn_Good_RF_stable(Sel_Others)');
PSTH_TgtOn_Bad_RF_Stable_Others = cell2mat(PSTH_TgtOn_Bad_RF_stable(Sel_Others)');

PSTH_TgtOn_Good_RF_Stable_Fast = cell2mat(PSTH_TgtOn_Good_RF_stable(Sel_Fast)');
PSTH_TgtOn_Bad_RF_Stable_Fast = cell2mat(PSTH_TgtOn_Bad_RF_stable(Sel_Fast)');





%Latency time
ValidBin = [-500,300];
LatencyTimeAUC_Stable_RF(LatencyTimeAUC_Stable_RF > ValidBin(2) | LatencyTimeAUC_Stable_RF < ValidBin(1))=NaN;
LatencyTimeAUC_Stable_ranksumRF(LatencyTimeAUC_Stable_ranksumRF > ValidBin(2) | LatencyTimeAUC_Stable_ranksumRF < ValidBin(1))=NaN;

%peakAUC_OneDR_RF_Fast =peakAUC_OneDR_RF(Sel_Fast);
%p_peakAUC_OneDR_RF_Fast = p_peakAUC_OneDR_RF(Sel_Fast);

%peakTime_AUC_OneDR_RF_Fast = peakTime_AUC_OneDR_RF(Sel_Fast);
LatencyTimeAUC_OneDR_RF_All = LatencyTimeAUC_OneDR_RF(Sel_All);
LatencyTimeAUC_OneDR_RF_Fast = LatencyTimeAUC_OneDR_RF(Sel_Fast);
LatencyTimeAUC_OneDR_RF_ranksumAll = LatencyTimeAUC_OneDR_ranksumRF(Sel_All);
LatencyTimeAUC_OneDR_RF_ranksumFast = LatencyTimeAUC_OneDR_ranksumRF(Sel_Fast);

%peakAUC_OneDR_RF_Others =peakAUC_OneDR_RF(Sel_Others);
%p_peakAUC_OneDR_RF_Others = p_peakAUC_OneDR_RF(Sel_Others);

%peakTime_AUC_OneDR_RF_Others = peakTime_AUC_OneDR_RF(Sel_Others);
LatencyTimeAUC_OneDR_RF_Others = LatencyTimeAUC_OneDR_RF(Sel_Others);
LatencyTimeAUC_OneDR_RF_ranksumOthers = LatencyTimeAUC_OneDR_ranksumRF(Sel_Others);
%peakAUC_Stable_RF_Fast =peakAUC_Stable_RF(Sel_Fast);
%p_peakAUC_Stable_RF_Fast = p_peakAUC_Stable_RF(Sel_Fast);

%peakTime_AUC_Stable_RF_Fast = peakTime_AUC_Stable_RF(Sel_Fast);
LatencyTimeAUC_Stable_RF_All = LatencyTimeAUC_Stable_RF(Sel_All);
LatencyTimeAUC_Stable_RF_Fast = LatencyTimeAUC_Stable_RF(Sel_Fast);
LatencyTimeAUC_Stable_RF_ranksumFast = LatencyTimeAUC_Stable_ranksumRF(Sel_Fast);
LatencyTimeAUC_Stable_RF_ranksumAll = LatencyTimeAUC_Stable_ranksumRF(Sel_All);


%peakAUC_Stable_RF_Others =peakAUC_Stable_RF(Sel_Others);
%p_peakAUC_Stable_RF_Others = p_peakAUC_Stable_RF(Sel_Others);

%peakTime_AUC_Stable_RF_Others = peakTime_AUC_Stable_RF(Sel_Others);
LatencyTimeAUC_Stable_RF_Others = LatencyTimeAUC_Stable_RF(Sel_Others);
LatencyTimeAUC_Stable_RF_ranksumOthers = LatencyTimeAUC_Stable_ranksumRF(Sel_Others);

%% Load the AUCtarget onset

TgtAUC_RF_All = TgtAUC_RF(Sel_All);
TgtAUC_RF_Fast = TgtAUC_RF(Sel_Fast);
TgtAUC_RF_Others = TgtAUC_RF(Sel_Others);

p_AUC_Tgt_RF_All = p_AUC_Tgt_RF(Sel_All);
p_AUC_Tgt_RF_Fast = p_AUC_Tgt_RF(Sel_Fast);
p_AUC_Tgt_RF_Others = p_AUC_Tgt_RF(Sel_Others);

p_AUC_Tgt_RF_ranksumAll = p_AUC_Tgt_RF_ranksum(Sel_All);
p_AUC_Tgt_RF_ranksumFast = p_AUC_Tgt_RF_ranksum(Sel_Fast);
p_AUC_Tgt_RF_ranksumOthers = p_AUC_Tgt_RF_ranksum(Sel_Others);


TgtAUC_RF_OneDRAll = TgtAUC_oneDR_RF(Sel_All);
TgtAUC_RF_OneDRFast = TgtAUC_oneDR_RF(Sel_Fast);
TgtAUC_RF_OneDR_Others = TgtAUC_oneDR_RF(Sel_Others);

p_AUC_Tgt_RF_OneDR_All = p_AUC_Tgt_OneDR_RF(Sel_All);
p_AUC_Tgt_RF_OneDR_Fast = p_AUC_Tgt_OneDR_RF(Sel_Fast);
p_AUC_Tgt_RF_OneDR_Others =p_AUC_Tgt_OneDR_RF(Sel_Others);

p_AUC_Tgt_RF_OneDR_ranksumAll = p_AUC_Tgt_OneDR_RF_ranksum(Sel_All);
p_AUC_Tgt_RF_OneDR_ranksumFast = p_AUC_Tgt_OneDR_RF_ranksum(Sel_Fast);
p_AUC_Tgt_RF_OneDR_ranksumOthers =p_AUC_Tgt_OneDR_RF_ranksum(Sel_Others);





%% Rectified AUC: fold AUC < 0.5 so all values represent distance from chance
rectifyAUC = @(x) x .* (x >= 0.5) + (1 - x) .* (x < 0.5);

% Rectify scalar TgtAUC
TgtAUC_RF_rect        = rectifyAUC(TgtAUC_RF);
TgtAUC_RF_rect_Fast   = TgtAUC_RF_rect(Sel_Fast);
TgtAUC_RF_rect_Others = TgtAUC_RF_rect(Sel_Others);
TgtAUC_RF_rect_All   = TgtAUC_RF_rect(Sel_All);

TgtAUC_oneDR_RF_rect        = rectifyAUC(TgtAUC_oneDR_RF);
TgtAUC_oneDR_RF_rect_Fast   = TgtAUC_oneDR_RF_rect(Sel_Fast);
TgtAUC_oneDR_RF_rect_Others = TgtAUC_oneDR_RF_rect(Sel_Others);
TgtAUC_oneDR_RF_rect_All = TgtAUC_oneDR_RF_rect(Sel_All);



% Rectify AUC PSTH: flip the entire time series based on scalar TgtAUC
% If TgtAUC < 0.5 for a neuron, apply 1-AUC to the whole sequence
AUCPSTH_TgtOn_Value_RF_stable_All = cell2mat(AUCPSTH_TgtOn_Value_RF_stable(Sel_All)');
AUCPSTH_TgtOn_Value_RF_stable_rect_All = AUCPSTH_TgtOn_Value_RF_stable_All ;
flip_stable = TgtAUC_RF_All < 0.5;
AUCPSTH_TgtOn_Value_RF_stable_rect_All(flip_stable',:) = 1-AUCPSTH_TgtOn_Value_RF_stable_rect_All(flip_stable',:);

AUCPSTH_TgtOn_Value_RF_oneDR_All = cell2mat(AUCPSTH_TgtOn_Value_RF_oneDR(Sel_All)');
AUCPSTH_TgtOn_Value_RF_oneDR_rect_All = AUCPSTH_TgtOn_Value_RF_oneDR_All;
flip_oneDR = TgtAUC_RF_OneDRAll < 0.5;
AUCPSTH_TgtOn_Value_RF_oneDR_rect_All(flip_oneDR',:) = 1-AUCPSTH_TgtOn_Value_RF_oneDR_rect_All(flip_oneDR',:);


%Define subtypes 

StableOnly = p_AUC_Tgt_RF_ranksumAll < 0.05 & p_AUC_Tgt_RF_OneDR_ranksumAll > 0.05 ;
OneDROnly = p_AUC_Tgt_RF_ranksumAll > 0.05 & p_AUC_Tgt_RF_OneDR_ranksumAll < 0.05 ;
Both = p_AUC_Tgt_RF_ranksumAll < 0.05 & p_AUC_Tgt_RF_OneDR_ranksumAll < 0.05 ;
None = p_AUC_Tgt_RF_ranksumAll > 0.05 & p_AUC_Tgt_RF_OneDR_ranksumAll > 0.05 ;


%Calculate the number 
StableOnly_num = sum(StableOnly);
OneDROnly_num = sum(OneDROnly);
Both_num = sum(Both);
%{
%Find out examples for illustration
% Each figure: 4x4 = 16 subplots, one per neuron. Title = file name.
% Rows: top = Stable task (Good/Bad RF), bottom = OneDR task (Good/Bad RF)
% Red = Good (high value), Blue = Bad (low value); solid = RF, dashed = nonRF

NeuronIndicesAll = find(Sel_All);
nNeuronsAll = length(NeuronIndicesAll);
nPerFig = 16;
figCountEx = 0;

for n = 1:nNeuronsAll
    idx = NeuronIndicesAll(n);
    plotPos = mod(n-1, nPerFig) + 1;

    if plotPos == 1
        figCountEx = figCountEx + 1;
        figure(500 + figCountEx);
        clf;
        set(gcf, 'Name', sprintf('Example PSTHs (%d/%d)', figCountEx, ceil(nNeuronsAll/nPerFig)));
    end

    subplot(4, 4, plotPos);
    hold on;
    % Stable task (obj): red = high value, blue = low value
    plot(PSTH_Time_TgtOn_stable, PSTH_TgtOn_Good_RF_stable{idx}, 'Color', 'r',            'LineWidth', 1);
    plot(PSTH_Time_TgtOn_stable, PSTH_TgtOn_Bad_RF_stable{idx},  'Color', 'b',            'LineWidth', 1);
    % OneDR task (loc): purple = high value, orange = low value
    plot(PSTH_Time_TgtOn_OneDR,  PSTH_Good_TgtOn_RF_oneDR{idx},  'Color', [0.5 0 0.8],   'LineWidth', 1);
    plot(PSTH_Time_TgtOn_OneDR,  PSTH_Bad_TgtOn_RF_oneDR{idx},   'Color', [1 0.5 0],     'LineWidth', 1);
    xline(0, 'LineStyle',':','color','k', 'LineWidth', 0.8);
    hold off;
    [~, fname] = fileparts(FilesNameAll{idx});
    if StableOnly(n),    typeStr = 'ObjOnly';
    elseif OneDROnly(n), typeStr = 'LocOnly';
    elseif Both(n),      typeStr = 'Both';
    else,                typeStr = 'None';
    end
    chanStr = sprintf('Ch%d', ChannelFull(idx));
    title({fname, [typeStr '  ' chanStr]}, 'Interpreter', 'none', 'FontSize', 6, 'FontWeight', 'normal');
    xlabel('Time (ms)', 'FontSize', 6);
    ylabel('z-score',   'FontSize', 6);
    set(gca, 'FontSize', 6);
end

keyboard
%}
%% Figures
FigureStartNum=150;
FigureIndex=1;


if ShowFigureFlag

if contains(BatchFileName, 'FEF'),     areaLabel_fig = 'FEF';
elseif contains(BatchFileName, 'SC'),  areaLabel_fig = 'SC';
else,                                  areaLabel_fig = 'Unknown'; 
end

%% Figure 1 Scatter plot of the raw AUC All between oneDR and stable
figtitlestr{FigureIndex}='Scatter of the raw AUC for all';

%fig = plotAUCScatter(x, y, px, py, xLabel, yLabel, titleStr, isRectified,FigureStartNum)
fig{FigureIndex} = plotAUCScatter(TgtAUC_RF_All, TgtAUC_RF_OneDRAll, ...
    p_AUC_Tgt_RF_ranksumAll, p_AUC_Tgt_RF_OneDR_ranksumAll, ...
    'Object AUC (stable)', 'Location AUC (oneDR)', ...
    figtitlestr{FigureIndex},  false,FigureStartNum);



%% Figure 2 Scatter plot of the reticfied AUC All between oneDR and stable
FigureStartNum=FigureStartNum+1;
FigureIndex=FigureIndex+1;

figtitlestr{FigureIndex}='Scatter of the rectified AUC for all';

fig{FigureIndex} = plotAUCScatter(TgtAUC_RF_rect_All, TgtAUC_oneDR_RF_rect_All, ...
    p_AUC_Tgt_RF_ranksumAll, p_AUC_Tgt_RF_OneDR_ranksumAll, ...
    'Object rectified AUC (stable)', 'Location rectified AUC (oneDR)', ...
    figtitlestr{FigureIndex},  true, FigureStartNum);

%% Figure 3 Veinn Graph
FigureStartNum=FigureStartNum+1;
FigureIndex=FigureIndex+1;

figtitlestr{FigureIndex}=sprintf(' Venn — %s', areaLabel_fig);

fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[200 200 420 400],'Name',figtitlestr{FigureIndex});


hold on; axis equal off;

col_red = [0.85 0.15 0.15];
col_pur = [0.55 0.10 0.80];
col_org = [0.95 0.55 0.05];
alpha   = 0.6;


[H, S] = venn([StableOnly_num, OneDROnly_num, Both_num], ...
    'FaceColor', {col_red , col_pur}, ...
    'FaceAlpha', {alpha, alpha}, ...
    'EdgeColor', 'k');

% Label each zone: exclusive stable, exclusive oneDR, overlap
zoneStrs = {num2str(StableOnly_num), num2str(OneDROnly_num), num2str(Both_num)};
for zi = 1:3
    text(S.ZoneCentroid(zi,1), S.ZoneCentroid(zi,2), zoneStrs{zi}, ...
        'HorizontalAlignment', 'center', 'FontSize', 14, 'FontWeight', 'bold');
end

legend([H(1), H(2)], {'Stable (object)', 'OneDR (location)'}, ...
    'Location', 'southoutside', 'Orientation', 'horizontal', 'FontSize', 10);
legend('boxoff');
title(sprintf('%s — AUC significance (Sel\\_All n=%d)', areaLabel_fig, sum(Sel_All)), ...
    'FontSize', 11);

%% Figures 4-6: Fast group scatter (raw + rectified) + Venn
FigureStartNum = FigureStartNum + 1;
FigureIndex = FigureIndex + 1;

StableOnly_Fast = p_AUC_Tgt_RF_ranksumFast < 0.05 & p_AUC_Tgt_RF_OneDR_ranksumFast > 0.05;
OneDROnly_Fast  = p_AUC_Tgt_RF_ranksumFast > 0.05 & p_AUC_Tgt_RF_OneDR_ranksumFast < 0.05;
Both_Fast       = p_AUC_Tgt_RF_ranksumFast < 0.05 & p_AUC_Tgt_RF_OneDR_ranksumFast < 0.05;

figtitlestr{FigureIndex} = sprintf('Scatter raw AUC — %s Fast', areaLabel_fig);
fig{FigureIndex} = plotAUCScatter(TgtAUC_RF_Fast, TgtAUC_RF_OneDRFast, ...
    p_AUC_Tgt_RF_ranksumFast, p_AUC_Tgt_RF_OneDR_ranksumFast, ...
    'Object AUC (stable)', 'Location AUC (oneDR)', ...
    figtitlestr{FigureIndex}, false, FigureStartNum);

FigureStartNum = FigureStartNum + 1;
FigureIndex = FigureIndex + 1;
figtitlestr{FigureIndex} = sprintf('Scatter rectified AUC — %s Fast', areaLabel_fig);
fig{FigureIndex} = plotAUCScatter(TgtAUC_RF_rect_Fast, TgtAUC_oneDR_RF_rect_Fast, ...
    p_AUC_Tgt_RF_ranksumFast, p_AUC_Tgt_RF_OneDR_ranksumFast, ...
    'Object rectified AUC (stable)', 'Location rectified AUC (oneDR)', ...
    figtitlestr{FigureIndex}, true, FigureStartNum);

FigureStartNum = FigureStartNum + 1;
FigureIndex = FigureIndex + 1;
figtitlestr{FigureIndex} = sprintf('Venn — %s Fast', areaLabel_fig);
fig{FigureIndex} = PrepareFigure(FigureStartNum, 'w', [200 200 420 400], 'Name', figtitlestr{FigureIndex});
hold on; axis equal off;
[H_f, S_f] = venn([sum(StableOnly_Fast), sum(OneDROnly_Fast), sum(Both_Fast)], ...
    'FaceColor', {col_red, col_pur}, ...
    'FaceAlpha', {alpha, alpha}, ...
    'EdgeColor', 'k');
zoneStrs_f = {num2str(sum(StableOnly_Fast)), num2str(sum(OneDROnly_Fast)), num2str(sum(Both_Fast))};
for zi = 1:3
    text(S_f.ZoneCentroid(zi,1), S_f.ZoneCentroid(zi,2), zoneStrs_f{zi}, ...
        'HorizontalAlignment', 'center', 'FontSize', 14, 'FontWeight', 'bold');
end
legend([H_f(1), H_f(2)], {'Stable (object)', 'OneDR (location)'}, ...
    'Location', 'southoutside', 'Orientation', 'horizontal', 'FontSize', 10);
legend('boxoff');
title(sprintf('%s — AUC significance (Fast n=%d)', areaLabel_fig, sum(Sel_Fast)), 'FontSize', 11);

%% Figures 7-9: Others group scatter (raw + rectified) + Venn
FigureStartNum = FigureStartNum + 1;
FigureIndex = FigureIndex + 1;

StableOnly_Others = p_AUC_Tgt_RF_ranksumOthers < 0.05 & p_AUC_Tgt_RF_OneDR_ranksumOthers > 0.05;
OneDROnly_Others  = p_AUC_Tgt_RF_ranksumOthers > 0.05 & p_AUC_Tgt_RF_OneDR_ranksumOthers < 0.05;
Both_Others       = p_AUC_Tgt_RF_ranksumOthers < 0.05 & p_AUC_Tgt_RF_OneDR_ranksumOthers < 0.05;

figtitlestr{FigureIndex} = sprintf('Scatter raw AUC — %s Others', areaLabel_fig);
fig{FigureIndex} = plotAUCScatter(TgtAUC_RF_Others, TgtAUC_RF_OneDR_Others, ...
    p_AUC_Tgt_RF_ranksumOthers, p_AUC_Tgt_RF_OneDR_ranksumOthers, ...
    'Object AUC (stable)', 'Location AUC (oneDR)', ...
    figtitlestr{FigureIndex}, false, FigureStartNum);

FigureStartNum = FigureStartNum + 1;
FigureIndex = FigureIndex + 1;
figtitlestr{FigureIndex} = sprintf('Scatter rectified AUC — %s Others', areaLabel_fig);
fig{FigureIndex} = plotAUCScatter(TgtAUC_RF_rect_Others, TgtAUC_oneDR_RF_rect_Others, ...
    p_AUC_Tgt_RF_ranksumOthers, p_AUC_Tgt_RF_OneDR_ranksumOthers, ...
    'Object rectified AUC (stable)', 'Location rectified AUC (oneDR)', ...
    figtitlestr{FigureIndex}, true, FigureStartNum);

FigureStartNum = FigureStartNum + 1;
FigureIndex = FigureIndex + 1;
figtitlestr{FigureIndex} = sprintf('Venn — %s Others', areaLabel_fig);
fig{FigureIndex} = PrepareFigure(FigureStartNum, 'w', [200 200 420 400], 'Name', figtitlestr{FigureIndex});
hold on; axis equal off;
[H_o, S_o] = venn([sum(StableOnly_Others), sum(OneDROnly_Others), sum(Both_Others)], ...
    'FaceColor', {col_red, col_pur}, ...
    'FaceAlpha', {alpha, alpha}, ...
    'EdgeColor', 'k');
zoneStrs_o = {num2str(sum(StableOnly_Others)), num2str(sum(OneDROnly_Others)), num2str(sum(Both_Others))};
for zi = 1:3
    text(S_o.ZoneCentroid(zi,1), S_o.ZoneCentroid(zi,2), zoneStrs_o{zi}, ...
        'HorizontalAlignment', 'center', 'FontSize', 14, 'FontWeight', 'bold');
end
legend([H_o(1), H_o(2)], {'Stable (object)', 'OneDR (location)'}, ...
    'Location', 'southoutside', 'Orientation', 'horizontal', 'FontSize', 10);
legend('boxoff');
title(sprintf('%s — AUC significance (Others n=%d)', areaLabel_fig, sum(Sel_Others)), 'FontSize', 11);



%Figure 4 population rectified AUC average for three subtypes
%% Figure: rectified AUC PSTH per significance group (stable + oneDR overlaid)
%Rectify according to the group type
AUCPSTH_TgtOn_Value_RF_stable_All_recgroup=AUCPSTH_TgtOn_Value_RF_stable_All;
AUCPSTH_TgtOn_Value_RF_oneDR_All_recgroup=AUCPSTH_TgtOn_Value_RF_oneDR_All;

%StableOnly

AUCPSTH_TgtOn_Value_RF_stable_All_recgroup(StableOnly & TgtAUC_RF_All<0.5,:) = 1 - AUCPSTH_TgtOn_Value_RF_stable_All_recgroup(StableOnly &TgtAUC_RF_All<0.5,:);
AUCPSTH_TgtOn_Value_RF_oneDR_All_recgroup(StableOnly & TgtAUC_RF_All<0.5,:) = 1 - AUCPSTH_TgtOn_Value_RF_oneDR_All_recgroup(StableOnly &TgtAUC_RF_All<0.5,:);

%OneDROnly
AUCPSTH_TgtOn_Value_RF_stable_All_recgroup(OneDROnly & TgtAUC_RF_OneDRAll<0.5,:) = 1 - AUCPSTH_TgtOn_Value_RF_stable_All_recgroup(OneDROnly  &TgtAUC_RF_OneDRAll<0.5,:);
AUCPSTH_TgtOn_Value_RF_oneDR_All_recgroup(OneDROnly & TgtAUC_RF_OneDRAll<0.5,:) = 1 - AUCPSTH_TgtOn_Value_RF_oneDR_All_recgroup(OneDROnly  &TgtAUC_RF_OneDRAll<0.5,:);

%Both
AUCPSTH_TgtOn_Value_RF_stable_All_recgroup(Both & TgtAUC_RF_All<0.5,:) = 1 - AUCPSTH_TgtOn_Value_RF_stable_All_recgroup(Both & TgtAUC_RF_All<0.5,:);
AUCPSTH_TgtOn_Value_RF_oneDR_All_recgroup(Both & TgtAUC_RF_OneDRAll<0.5,:) = 1 - AUCPSTH_TgtOn_Value_RF_oneDR_All_recgroup(Both & TgtAUC_RF_OneDRAll<0.5,:);


groupMasks  = {StableOnly, OneDROnly, Both, None};
groupLabels = {'Object-value only', 'Location only', 'Both significant', 'Neither'};

t_min_disp = -500;
t_max_disp = 300;

figure('Name', sprintf('Rectified AUC PSTH by sig group — %s', areaLabel_fig), ...
    'Position', [50 100 1400 350], 'Color', 'w');

for gi = 1:4
    gMask = groupMasks{gi};
    nG    = sum(gMask);
    subplot(1, 4, gi); hold on;

    if nG > 0
        % Stable (blue) — fill hidden from legend
        mat_s = AUCPSTH_TgtOn_Value_RF_stable_All_recgroup(gMask,:);
        avg_s = nanmean(mat_s, 1);
        sem_s = nanstd(mat_s, [], 1) / sqrt(nG);
        fill([PSTH_Time_TgtOn_stable, fliplr(PSTH_Time_TgtOn_stable)], ...
             [avg_s + sem_s, fliplr(avg_s - sem_s)], ...
             col_red, 'FaceAlpha', 0.25, 'EdgeColor', 'none', 'HandleVisibility', 'off');
        plot(PSTH_Time_TgtOn_stable, avg_s, '-', 'Color', col_red, 'LineWidth', 1.5, ...
             'DisplayName', 'Stable');

        % OneDR (orange) — fill hidden from legend
        mat_d = AUCPSTH_TgtOn_Value_RF_oneDR_All_recgroup(gMask,:);
        avg_d = nanmean(mat_d, 1);
        sem_d = nanstd(mat_d, [], 1) / sqrt(nG);
        fill([PSTH_Time_TgtOn_OneDR, fliplr(PSTH_Time_TgtOn_OneDR)], ...
             [avg_d + sem_d, fliplr(avg_d - sem_d)], ...
             col_pur, 'FaceAlpha', 0.25, 'EdgeColor', 'none', 'HandleVisibility', 'off');
        plot(PSTH_Time_TgtOn_OneDR, avg_d, '-', 'Color',col_pur, 'LineWidth', 1.5, ...
             'DisplayName', 'OneDR');
    end

    yline(0.5, 'LineStyle', '--', 'Color', 'k', 'LineWidth', 0.8);
    xline(0,   'LineStyle', '--', 'Color', 'k', 'LineWidth', 0.8);
    xlim([t_min_disp, t_max_disp]);
    xlabel('Time from target onset (ms)');
    ylabel('AUC (rectified)');
    title(sprintf('%s\nn=%d', groupLabels{gi}, nG));
    %legend('Location', 'northwest', 'FontSize', 8);
    %ylim([0.4,1]);
    box off;

    % Per-bin significance bars (ranksum vs 0.5 chance) after axes are set
    if nG > 0
        yL   = ylim;
        barH = 0.008 * diff(yL);
        barY_s = yL(2) + 0.02 * diff(yL);
        barY_d = yL(2) + 0.05 * diff(yL);
        plotSigBinsLocal(PSTH_Time_TgtOn_stable, mat_s, barY_s, barH, col_red);
        plotSigBinsLocal(PSTH_Time_TgtOn_OneDR,  mat_d, barY_d, barH, col_pur);
    end

    ylim([0.4,1]);

end
sgtitle(sprintf('%s — Rectified AUC PSTH by significance group (Sel\\_All)', areaLabel_fig));

%% Figure 5 
%% AUC correlation heatmaps: OneDR vs Stable (per-neuron, all time bins)
% Two heatmaps side-by-side:
%   (1) Raw (unrectified) AUC PSTH
%   (2) Rectified AUC PSTH — timeseries flipped per neuron if TgtAUC < 0.5

% --- Assemble raw AUCPSTH matrices [N x T] ---


T_stable = size(AUCPSTH_TgtOn_Value_RF_stable_All, 2);
T_oneDR  = size(AUCPSTH_TgtOn_Value_RF_oneDR_All,  2);
N_all    = size(TgtAUC_RF_All, 2);

fprintf('Computing raw AUC correlation matrix (%d x %d, n=%d neurons)...\n', ...
    T_stable, T_oneDR, N_all);

% --- (1) Correlation of raw AUC ---
auc_corr_raw = NaN(T_stable, T_oneDR);
pauc_corr_raw = NaN(T_stable, T_oneDR);
for ts = 1:T_stable
    v_s = AUCPSTH_TgtOn_Value_RF_stable_All(:, ts);
    for td = 1:T_oneDR
        v_d = AUCPSTH_TgtOn_Value_RF_oneDR_All(:, td);
        ok  = ~isnan(v_s) & ~isnan(v_d);
        if sum(ok) >= 3
            [r,p] = corrcoef(v_s(ok), v_d(ok));
            auc_corr_raw(ts, td) = r(1, 2);
            pauc_corr_raw(ts, td) = p(1, 2);
        end
    end
end

figure('Name', sprintf('AUC Correlation — Raw (%s)', areaLabel_fig), ...
       'Position', [100 100 580 480], 'Color', 'w');
imagesc(PSTH_Time_TgtOn_OneDR, PSTH_Time_TgtOn_stable, auc_corr_raw, [-1 1]);
axis xy square;
hold on;
xline(0, 'LineStyle','--','Color','w','LineWidth',1.2);
yline(0, 'LineStyle','--','Color','w','LineWidth',1.2);
contour(PSTH_Time_TgtOn_OneDR, PSTH_Time_TgtOn_stable, pauc_corr_raw, [0.05 0.05], 'r', 'LineWidth', 1.5);
cb = colorbar; cb.Label.String = 'Pearson r';
colormap(jet(256));
xlabel('Time — OneDR (ms)');
ylabel('Time — Stable (ms)');
title(sprintf('AUC correlation (raw): OneDR vs Stable, RF (%s)  [n=%d]', areaLabel_fig, N_all));
box off;

% --- (2) Correlation of group rectified AUC (flip timeseries if TgtAUC < 0.5) ---
fprintf('Computing rectified AUC correlation matrix...\n');

auc_corr_rect = NaN(T_stable, T_oneDR);
pauc_corr_rect = NaN(T_stable, T_oneDR);
for ts = 1:T_stable
    v_s = AUCPSTH_TgtOn_Value_RF_stable_All_recgroup(:, ts);
    for td = 1:T_oneDR
        v_d = AUCPSTH_TgtOn_Value_RF_oneDR_All_recgroup(:, td);
        ok  = ~isnan(v_s) & ~isnan(v_d);
        if sum(ok) >= 3
            [r,p] = corrcoef(v_s(ok), v_d(ok));
            auc_corr_rect(ts, td) = r(1, 2);
            pauc_corr_rect(ts, td) = p(1, 2);
        end
    end
end

figure('Name', sprintf('AUC Correlation — GroupRectified by TgtAUC (%s)', areaLabel_fig), ...
       'Position', [700 100 580 480], 'Color', 'w');
imagesc(PSTH_Time_TgtOn_OneDR, PSTH_Time_TgtOn_stable, auc_corr_rect, [-1 1]);
axis xy square;
hold on;
xline(0, 'LineStyle','--','Color','w','LineWidth',1.2);
yline(0, 'LineStyle','--','Color','w','LineWidth',1.2);
contour(PSTH_Time_TgtOn_OneDR, PSTH_Time_TgtOn_stable, pauc_corr_rect, [0.05 0.05], 'r', 'LineWidth', 1.5);
cb = colorbar; cb.Label.String = 'Pearson r';
colormap(jet(256));
xlabel('Time — OneDR (ms)');
ylabel('Time — Stable (ms)');
title(sprintf('AUC correlation (group rect by TgtAUC): OneDR vs Stable, RF (%s)  [n=%d]', areaLabel_fig, N_all));
box off;

% --- (3) Correlation of separate rectified AUC (flip timeseries if TgtAUC < 0.5) ---
fprintf('Computing sep rectified AUC correlation matrix...\n');

auc_corr_seprect = NaN(T_stable, T_oneDR);
pauc_corr_seprect = NaN(T_stable, T_oneDR);
for ts = 1:T_stable
    v_s = AUCPSTH_TgtOn_Value_RF_stable_rect_All(:, ts);
    for td = 1:T_oneDR
        v_d = AUCPSTH_TgtOn_Value_RF_oneDR_rect_All(:, td);
        ok  = ~isnan(v_s) & ~isnan(v_d);
        if sum(ok) >= 3
            [r,p] = corrcoef(v_s(ok), v_d(ok));
            auc_corr_seprect(ts, td) = r(1, 2);
            pauc_corr_seprect(ts, td) = p(1, 2);
        end
    end
end

figure('Name', sprintf('AUC Correlation — SepRectified by TgtAUC (%s)', areaLabel_fig), ...
       'Position', [700 100 580 480], 'Color', 'w');
imagesc(PSTH_Time_TgtOn_OneDR, PSTH_Time_TgtOn_stable, auc_corr_seprect, [-1 1]);
axis xy square;
hold on;
xline(0, 'LineStyle','--','Color','w','LineWidth',1.2);
yline(0, 'LineStyle','--','Color','w','LineWidth',1.2);
contour(PSTH_Time_TgtOn_OneDR, PSTH_Time_TgtOn_stable, pauc_corr_seprect, [0.05 0.05], 'r', 'LineWidth', 1.5);
cb = colorbar; cb.Label.String = 'Pearson r';
colormap(jet(256));
xlabel('Time — OneDR (ms)');
ylabel('Time — Stable (ms)');
title(sprintf('AUC correlation (sep rect by TgtAUC): OneDR vs Stable, RF (%s)  [n=%d]', areaLabel_fig, N_all));
box off;


end




%keyboard


%{

%% Population average: rectified AUC PSTH — setup


t_min_disp = max(PSTH_Time_TgtOn_stable(1),  PSTH_Time_TgtOn_OneDR(1));
t_max_disp = min(PSTH_Time_TgtOn_stable(end), PSTH_Time_TgtOn_OneDR(end));

% Compute All-neuron averages (Sel_All)

AUC_rect_stable_All_avg = nanmean(AUCPSTH_TgtOn_Value_RF_stable_rect, 1);
AUC_rect_stable_All_sem = nanstd(AUCPSTH_TgtOn_Value_RF_stable_rect, [], 1) / sqrt(sum(Sel_All));


AUC_rect_oneDR_All_avg = nanmean(AUCPSTH_TgtOn_Value_RF_oneDR_rect, 1);
AUC_rect_oneDR_All_sem = nanstd(AUCPSTH_TgtOn_Value_RF_oneDR_rect, [], 1) / sqrt(sum(Sel_All));

c_stable = [0.20 0.40 0.90];   % blue   — stable task
c_oneDR  = [0.90 0.40 0.10];   % orange — oneDR task
c_fast   = [0.85 0.15 0.15];   % red    — Fast group
c_others = [0.20 0.60 0.20];   % green  — Others group

%% Figure 1: All neurons — Stable vs OneDR
figure('Name', sprintf('AUC PSTH (All) — %s', areaLabel_fig), ...
    'Position', [100 100 600 420], 'Color', 'w');
hold on;
fill([PSTH_Time_TgtOn_stable, fliplr(PSTH_Time_TgtOn_stable)], ...
     [AUC_rect_stable_All_avg + AUC_rect_stable_All_sem, ...
      fliplr(AUC_rect_stable_All_avg - AUC_rect_stable_All_sem)], ...
     c_stable, 'FaceAlpha', 0.25, 'EdgeColor', 'none');
plot(PSTH_Time_TgtOn_stable, AUC_rect_stable_All_avg, '-', ...
     'Color', c_stable, 'LineWidth', 1.5, 'DisplayName', 'Stable');
fill([PSTH_Time_TgtOn_OneDR, fliplr(PSTH_Time_TgtOn_OneDR)], ...
     [AUC_rect_oneDR_All_avg + AUC_rect_oneDR_All_sem, ...
      fliplr(AUC_rect_oneDR_All_avg - AUC_rect_oneDR_All_sem)], ...
     c_oneDR, 'FaceAlpha', 0.25, 'EdgeColor', 'none');
plot(PSTH_Time_TgtOn_OneDR, AUC_rect_oneDR_All_avg, '-', ...
     'Color', c_oneDR, 'LineWidth', 1.5, 'DisplayName', 'OneDR');
yline(0.5,'LineStyle', '--', 'Color','k','LineWidth', 1);
xline(0, 'LineStyle',  '--','Color','k', 'LineWidth', 1);
xlabel('Time from target onset (ms)');
ylabel('AUC (rectified)');
title(sprintf('%s — All neurons  [n=%d]', areaLabel_fig, sum(Sel_All)));
legend('Location', 'northwest', 'FontSize', 9);
xlim([t_min_disp, t_max_disp]);

box off;

%% Figure 2: Fast vs Others — one subplot per task
figure('Name', sprintf('AUC PSTH (Fast vs Others) — %s', areaLabel_fig), ...
    'Position', [150 150 600 700], 'Color', 'w');

% --- Stable task: Fast vs Others ---
subplot(2,1,1); hold on;
fill([PSTH_Time_TgtOn_stable, fliplr(PSTH_Time_TgtOn_stable)], ...
     [AUC_rect_stable_Others_avg + AUC_rect_stable_Others_sem, ...
      fliplr(AUC_rect_stable_Others_avg - AUC_rect_stable_Others_sem)], ...
     c_others, 'FaceAlpha', 0.25, 'EdgeColor', 'none');
plot(PSTH_Time_TgtOn_stable, AUC_rect_stable_Others_avg, '-', ...
     'Color', c_others, 'LineWidth', 1.5, ...
     'DisplayName', sprintf('Others [n=%d]', sum(Sel_Others)));
fill([PSTH_Time_TgtOn_stable, fliplr(PSTH_Time_TgtOn_stable)], ...
     [AUC_rect_stable_Fast_avg + AUC_rect_stable_Fast_sem, ...
      fliplr(AUC_rect_stable_Fast_avg - AUC_rect_stable_Fast_sem)], ...
     c_fast, 'FaceAlpha', 0.25, 'EdgeColor', 'none');
plot(PSTH_Time_TgtOn_stable, AUC_rect_stable_Fast_avg, '-', ...
     'Color', c_fast, 'LineWidth', 1.5, ...
     'DisplayName', sprintf('Fast [n=%d]', sum(Sel_Fast)));
yline(0.5, 'LineStyle','--','Color','k', 'LineWidth', 1);
xline(0,   'LineStyle','--','Color','k',  'LineWidth', 1);
xlabel('Time from target onset (ms)');
ylabel('AUC (rectified)');
title(sprintf('%s — Stable task: Fast vs Others', areaLabel_fig));
legend('Location', 'northwest', 'FontSize', 9);
xlim([t_min_disp, t_max_disp]);
box off;

% --- OneDR task: Fast vs Others ---
subplot(2,1,2); hold on;
fill([PSTH_Time_TgtOn_OneDR, fliplr(PSTH_Time_TgtOn_OneDR)], ...
     [AUC_rect_oneDR_Others_avg + AUC_rect_oneDR_Others_sem, ...
      fliplr(AUC_rect_oneDR_Others_avg - AUC_rect_oneDR_Others_sem)], ...
     c_others, 'FaceAlpha', 0.25, 'EdgeColor', 'none');
plot(PSTH_Time_TgtOn_OneDR, AUC_rect_oneDR_Others_avg, '-', ...
     'Color', c_others, 'LineWidth', 1.5, ...
     'DisplayName', sprintf('Others [n=%d]', sum(Sel_Others)));
fill([PSTH_Time_TgtOn_OneDR, fliplr(PSTH_Time_TgtOn_OneDR)], ...
     [AUC_rect_oneDR_Fast_avg + AUC_rect_oneDR_Fast_sem, ...
      fliplr(AUC_rect_oneDR_Fast_avg - AUC_rect_oneDR_Fast_sem)], ...
     c_fast, 'FaceAlpha', 0.25, 'EdgeColor', 'none');
plot(PSTH_Time_TgtOn_OneDR, AUC_rect_oneDR_Fast_avg, '-', ...
     'Color', c_fast, 'LineWidth', 1.5, ...
     'DisplayName', sprintf('Fast [n=%d]', sum(Sel_Fast)));
yline(0.5, 'LineStyle','--','Color','k','LineWidth', 1);
xline(0,   'LineStyle','--','Color','k', 'LineWidth', 1);
xlabel('Time from target onset (ms)');
ylabel('AUC (rectified)');
title(sprintf('%s — OneDR task: Fast vs Others', areaLabel_fig));
legend('Location', 'northwest', 'FontSize', 9);
xlim([t_min_disp, t_max_disp]);
box off;

%% AUC PSTH average by significance & direction — Stable Value (RF)
% Fast group
PlotAUCPSTHBySignificance(AUCPSTH_TgtOn_Value_RF_stable(Sel_Fast), ...
    p_AUC_Tgt_RF(Sel_Fast), TgtAUC_RF(Sel_Fast), ...
    PSTH_Time_TgtOn_stable, 'Stable Value — AUC PSTH (RF) [Fast]');

% Others group
PlotAUCPSTHBySignificance(AUCPSTH_TgtOn_Value_RF_stable(Sel_Others), ...
    p_AUC_Tgt_RF(Sel_Others), TgtAUC_RF(Sel_Others), ...
    PSTH_Time_TgtOn_stable, 'Stable Value — AUC PSTH (RF) [Others]');

% Fast group
PlotAUCPSTHBySignificance(AUCPSTH_TgtOn_Value_RF_oneDR(Sel_Fast), ...
    p_AUC_Tgt_OneDR_RF(Sel_Fast), TgtAUC_oneDR_RF(Sel_Fast), ...
    PSTH_Time_TgtOn_OneDR, 'OneDR — AUC PSTH (RF) [Fast]');

% Others group
PlotAUCPSTHBySignificance(AUCPSTH_TgtOn_Value_RF_oneDR(Sel_Others), ...
    p_AUC_Tgt_OneDR_RF(Sel_Others), TgtAUC_oneDR_RF(Sel_Others), ...
    PSTH_Time_TgtOn_OneDR, 'OneDR — AUC PSTH (RF) [Others]');


%% Significance-based grouping: object value (stable) vs location (oneDR)
% Use ranksum p-values for higher power at small trial counts
alpha_sig = 0.05;

% Row-vector logical flags across all neurons
%sig_obj = (p_AUC_Tgt_RF_ranksum       < alpha_sig);   % significant in stable/object task
%sig_loc = (p_AUC_Tgt_OneDR_RF_ranksum < alpha_sig);   % significant in oneDR/location task

sig_obj = (p_AUC_Tgt_RF      < alpha_sig);   % significant in stable/object task
sig_loc = (p_AUC_Tgt_OneDR_RF < alpha_sig);   % significant in oneDR/location task

sel_row = Sel_All';   % row-vector to match orientation of sig_obj / sig_loc

g_obj_only = sig_obj & ~sig_loc & sel_row;   % object-value only
g_loc_only = ~sig_obj & sig_loc & sel_row;   % location only
g_both_sig = sig_obj &  sig_loc & sel_row;   % both significant
g_neither  = ~sig_obj & ~sig_loc & sel_row;  % neither

fprintf('--- AUC significance breakdown (Sel_All, n=%d) ---\n', sum(Sel_All));
fprintf('  Object-value only : %d\n', sum(g_obj_only));
fprintf('  Location only     : %d\n', sum(g_loc_only));
fprintf('  Both significant  : %d\n', sum(g_both_sig));
fprintf('  Neither           : %d\n', sum(g_neither));




%% Figure 1a — Scatter tgtAUC stable vs oneDR (raw, not rectified)
areaLabel_fig = 'Unknown';
if contains(BatchFileName, 'FEF'), areaLabel_fig = 'FEF'; end
if contains(BatchFileName, 'SC'),  areaLabel_fig = 'SC';  end

plotAUCScatter(TgtAUC_RF(Sel_All), TgtAUC_oneDR_RF(Sel_All), ...
    p_AUC_Tgt_RF(Sel_All), p_AUC_Tgt_OneDR_RF(Sel_All), ...
    'Object AUC (stable)', 'Location AUC (oneDR)', ...
    sprintf('tgtAUC stable vs oneDR — %s  [n=%d]', areaLabel_fig, sum(Sel_All)), false);

%% Figure 1b — Scatter tgtAUC stable vs oneDR (rectified)
plotAUCScatter(TgtAUC_RF_rect(Sel_All), TgtAUC_oneDR_RF_rect(Sel_All), ...
    p_AUC_Tgt_RF(Sel_All), p_AUC_Tgt_OneDR_RF(Sel_All), ...
    'Object AUC rect (stable)', 'Location AUC rect (oneDR)', ...
    sprintf('tgtAUC stable vs oneDR rectified — %s  [n=%d]', areaLabel_fig, sum(Sel_All)), true);

%keyboard
%}
%% Per-trial PSTH selection: sufficient trials in every RF condition (>=15)
TrialCount_OK = ...
    n_high_RF_stable >= 10 & ...
    n_low_RF_stable  >= 10 & ...
    n_high_RF_oneDR  >= 10 & ...
    n_low_RF_oneDR   >= 10;
TrialCount_OK = TrialCount_OK';   % column vector to match Sel_* vectors

Sel_Raw = Sel_All & TrialCount_OK;

% ModuType for selected neurons
ModuType_Sel = ModuType(Sel_Raw);

% psth_targeton_raw — cell arrays (variable n_trials per neuron)
psth_targeton_raw_stable_Sel = psth_targeton_raw_stable(Sel_Raw);
psth_targeton_raw_oneDR_Sel  = psth_targeton_raw_oneDR(Sel_Raw);

% target_index — cell arrays
target_index_stable_Sel = target_index_stable(Sel_Raw);
target_index_oneDR_Sel  = target_index_oneDR(Sel_Raw);

% RF_index — cell arrays
RF_index_stable_Sel = RF_index_stable(Sel_Raw);
RF_index_oneDR_Sel  = RF_index_oneDR(Sel_Raw);

% Split psth_targeton_raw into RF and nonRF trials using RF_index
% Each output is a cell array {1 x n_neurons}, each cell [n_trials x T]
splitPos = @(psth, idx)  psth( logical(idx), :);   % select rows where idx==1
splitNeg = @(psth, idx)  psth(~logical(idx), :);   % select rows where idx==0
selectEl = @(vec,  idx)  vec(  logical(idx));       % 1-D element select (any orientation)
selectElNeg = @(vec, idx) vec(~logical(idx));

% --- RF / nonRF split of psth ---
psth_targeton_raw_stable_RF    = cellfun(splitPos, psth_targeton_raw_stable_Sel, RF_index_stable_Sel, 'UniformOutput', false);
psth_targeton_raw_stable_nonRF = cellfun(splitNeg, psth_targeton_raw_stable_Sel, RF_index_stable_Sel, 'UniformOutput', false);

psth_targeton_raw_oneDR_RF    = cellfun(splitPos, psth_targeton_raw_oneDR_Sel, RF_index_oneDR_Sel, 'UniformOutput', false);
psth_targeton_raw_oneDR_nonRF = cellfun(splitNeg, psth_targeton_raw_oneDR_Sel, RF_index_oneDR_Sel, 'UniformOutput', false);

% --- RF / nonRF split of target_index (needed for high/low split below) ---
target_index_stable_RF    = cellfun(selectEl,    target_index_stable_Sel, RF_index_stable_Sel, 'UniformOutput', false);
target_index_stable_nonRF = cellfun(selectElNeg, target_index_stable_Sel, RF_index_stable_Sel, 'UniformOutput', false);

target_index_oneDR_RF    = cellfun(selectEl,    target_index_oneDR_Sel, RF_index_oneDR_Sel, 'UniformOutput', false);
target_index_oneDR_nonRF = cellfun(selectElNeg, target_index_oneDR_Sel, RF_index_oneDR_Sel, 'UniformOutput', false);

%% --- High / low value split (1=high, 0=low) within each RF/nonRF group ---
psth_targeton_raw_stable_RF_high    = cellfun(splitPos, psth_targeton_raw_stable_RF,    target_index_stable_RF,    'UniformOutput', false);
psth_targeton_raw_stable_RF_low     = cellfun(splitNeg, psth_targeton_raw_stable_RF,    target_index_stable_RF,    'UniformOutput', false);
psth_targeton_raw_stable_nonRF_high = cellfun(splitPos, psth_targeton_raw_stable_nonRF, target_index_stable_nonRF, 'UniformOutput', false);
psth_targeton_raw_stable_nonRF_low  = cellfun(splitNeg, psth_targeton_raw_stable_nonRF, target_index_stable_nonRF, 'UniformOutput', false);

psth_targeton_raw_oneDR_RF_high    = cellfun(splitPos, psth_targeton_raw_oneDR_RF,    target_index_oneDR_RF,    'UniformOutput', false);
psth_targeton_raw_oneDR_RF_low     = cellfun(splitNeg, psth_targeton_raw_oneDR_RF,    target_index_oneDR_RF,    'UniformOutput', false);
psth_targeton_raw_oneDR_nonRF_high = cellfun(splitPos, psth_targeton_raw_oneDR_nonRF, target_index_oneDR_nonRF, 'UniformOutput', false);
psth_targeton_raw_oneDR_nonRF_low  = cellfun(splitNeg, psth_targeton_raw_oneDR_nonRF, target_index_oneDR_nonRF, 'UniformOutput', false);

%% --- Z-score (RF): one mu and one sigma per neuron across all time bins and
%     all four trial conditions (stable high/low RF, oneDR high/low RF).
%     Applied only to the four RF high/low arrays below; nonRF arrays unchanged.
z_eps_rf = 1e-8;
N_rf_z = numel(psth_targeton_raw_stable_RF_high);
psth_RF_global_z_mu = nan(N_rf_z, 1);
psth_RF_global_z_sd = nan(N_rf_z, 1);
for nn = 1:N_rf_z
    pool_nn = [psth_targeton_raw_stable_RF_high{nn}(:); ...
               psth_targeton_raw_stable_RF_low{nn}(:); ...
               psth_targeton_raw_oneDR_RF_high{nn}(:); ...
               psth_targeton_raw_oneDR_RF_low{nn}(:)];
    if isempty(pool_nn)
        mu_nn = 0;
        sd_nn = 1;
    else
        mu_nn = mean(pool_nn);
        if numel(pool_nn) < 2
            sd_nn = 1;
        else
            sd_nn = std(pool_nn, 0);
            if ~(sd_nn > z_eps_rf)
                sd_nn = 1;
            else
                sd_nn = max(sd_nn, z_eps_rf);
            end
        end
    end
    psth_RF_global_z_mu(nn) = mu_nn;
    psth_RF_global_z_sd(nn) = sd_nn;

    psth_targeton_raw_stable_RF_high{nn} = (psth_targeton_raw_stable_RF_high{nn} - mu_nn) ./ sd_nn;
    psth_targeton_raw_stable_RF_low{nn}  = (psth_targeton_raw_stable_RF_low{nn}  - mu_nn) ./ sd_nn;
    psth_targeton_raw_oneDR_RF_high{nn}  = (psth_targeton_raw_oneDR_RF_high{nn}  - mu_nn) ./ sd_nn;
    psth_targeton_raw_oneDR_RF_low{nn}   = (psth_targeton_raw_oneDR_RF_low{nn}   - mu_nn) ./ sd_nn;
end
fprintf('[Batch_ObjLocValue_Organize] Z-scored RF PSTH (global mu/sd per neuron, 4 conditions pooled): N=%d\n', N_rf_z);

%
%% ── Build population PSTH matrices for TDR ───────────────────────────────────
% TDR_ValueContext_RF expects [nNeurons x nTimeBins] trial-averaged matrices
% on a single shared time axis. Crop stable and oneDR to their overlapping
% window first, then average across trials per neuron.
%{
% Overlapping time window (ms)
t_min_TDR = -500;%max(PSTH_Time_TgtOn_stable(1),  PSTH_Time_TgtOn_OneDR(1));
t_max_TDR = 300;%min(PSTH_Time_TgtOn_stable(end), PSTH_Time_TgtOn_OneDR(end));

idx_stable_TDR = PSTH_Time_TgtOn_stable >= t_min_TDR & PSTH_Time_TgtOn_stable <= t_max_TDR;
idx_oneDR_TDR  = PSTH_Time_TgtOn_OneDR  >= t_min_TDR & PSTH_Time_TgtOn_OneDR  <= t_max_TDR;

time_axis_TDR = PSTH_Time_TgtOn_stable(idx_stable_TDR);% / 1000;   % convert ms → seconds

% Average trials per neuron, crop to common window, stack → [nNeurons x T]
% cellArr(:) forces column orientation so cell2mat stacks rows, not columns.
%
% e.g. Others only: sel = ModuType_Sel == 4;
%      Fast only:   sel = ModuType_Sel == 1;
%      All neurons: sel = true(size(ModuType_Sel));

%sel =ModuType_Sel == 4;
sel = true(size(ModuType_Sel));

psth_targeton_raw_stable_RF_high = psth_targeton_raw_stable_RF_high(sel)';
psth_targeton_raw_stable_RF_low  = psth_targeton_raw_stable_RF_low(sel)';

psth_targeton_raw_oneDR_RF_high = psth_targeton_raw_oneDR_RF_high(sel)';
psth_targeton_raw_oneDR_RF_low  = psth_targeton_raw_oneDR_RF_low(sel)';


%% Preprocess (smooth)
[psth_targeton_raw_stable_RF_high, time_full] = PreprocessPSTH(psth_targeton_raw_stable_RF_high, PSTH_Time_TgtOn_stable, 'kernel', 'gaussian', 'win_ms', 20);
[psth_targeton_raw_stable_RF_low, ~]          = PreprocessPSTH(psth_targeton_raw_stable_RF_low',  PSTH_Time_TgtOn_stable,'kernel', 'gaussian', 'win_ms', 20);




%% Preprocess oneDR (smooth + select same time window [-500, 300] ms)
[psth_targeton_raw_oneDR_RF_high, time_full_oneDR] = PreprocessPSTH(psth_targeton_raw_oneDR_RF_high', PSTH_Time_TgtOn_OneDR, 'kernel', 'gaussian', 'win_ms', 20);
[psth_targeton_raw_oneDR_RF_low, ~]               = PreprocessPSTH(psth_targeton_raw_oneDR_RF_low,   PSTH_Time_TgtOn_OneDR, 'kernel', 'gaussian', 'win_ms', 20);




avgTrial = @(cellArr, tidx) cell2mat( ...
    cellfun(@(m) mean(m(:, tidx), 1), cellArr(:), 'UniformOutput', false) );

psth_high_value_object_RF   = avgTrial(psth_targeton_raw_stable_RF_high, idx_stable_TDR);
psth_low_value_object_RF    = avgTrial(psth_targeton_raw_stable_RF_low,  idx_stable_TDR);
psth_high_value_location_RF = avgTrial(psth_targeton_raw_oneDR_RF_high,  idx_oneDR_TDR);
psth_low_value_location_RF  = avgTrial(psth_targeton_raw_oneDR_RF_low,   idx_oneDR_TDR);
%}


%{

ComputeValueAxisOrthogonality( ...
    psth_high_value_object_RF, psth_low_value_object_RF, ...
    psth_high_value_location_RF, psth_low_value_location_RF, ...
    time_full_oneDR, varargin)
%}
%{
psth_high_value_object_RF = PSTH_TgtOn_Good_RF_Stable_All(:,idx_stable_TDR);
psth_low_value_object_RF = PSTH_TgtOn_Bad_RF_Stable_All(:,idx_stable_TDR);

psth_high_value_location_RF = PSTH_TgtOn_Good_RF_oneDR_All(:,idx_stable_TDR);
psth_low_value_location_RF = PSTH_TgtOn_Bad_RF_oneDR_All(:, idx_stable_TDR);
%}
%{
psth_high_value_object_RF = PSTH_TgtOn_Good_RF_Stable_Fast(:,idx_stable_TDR);
psth_low_value_object_RF = PSTH_TgtOn_Bad_RF_Stable_Fast(:,idx_stable_TDR);

psth_high_value_location_RF = PSTH_TgtOn_Good_RF_OneDR_Fast(:,idx_stable_TDR);
psth_low_value_location_RF = PSTH_TgtOn_Bad_RF_OneDR_Fast(:, idx_stable_TDR);
%}
%{
psth_high_value_object_RF = PSTH_TgtOn_Good_RF_Stable_Others(:,idx_stable_TDR);
psth_low_value_object_RF = PSTH_TgtOn_Bad_RF_Stable_Others(:,idx_stable_TDR);

psth_high_value_location_RF = PSTH_TgtOn_Good_RF_OneDR_Others(:,idx_stable_TDR);
psth_low_value_location_RF = PSTH_TgtOn_Bad_RF_OneDR_Others(:, idx_stable_TDR);
%}
%{
%% ── Call TDR ─────────────────────────────────────────────────────────────────
TDR_RF = TDR_ValueContext_RF( ...
    psth_high_value_object_RF,   psth_low_value_object_RF, ...
    psth_high_value_location_RF, psth_low_value_location_RF, ...
    time_axis_TDR, false);

%}

%{
Results = ComputeValueAxisOrthogonality( ...
    psth_high_value_object_RF, psth_low_value_object_RF, ...
    psth_high_value_location_RF, psth_low_value_location_RF, ...
    time_axis_TDR);
%}
if OutputFlag




%% Save 4 RF PSH matrices
% Determine brain area from batch file name
if contains(BatchFileName, 'FEF')
    AreaLabel = 'FEF';
elseif contains(BatchFileName, 'SC')
    AreaLabel = 'SC';
else
    AreaLabel = 'Unknown';
end

SaveRoot = '/Users/xuefeiyu/Documents/XuefeiFile/WorkRelated/Program_Matlab_Local/OHLab/DataAnalysis/Results/ObjectLocationPSTH';
SaveDir  = fullfile(SaveRoot, AreaLabel);
if ~exist(SaveDir, 'dir')
    mkdir(SaveDir);
end

% All neurons (with ModuType for downstream group selection)
SaveFile = fullfile(SaveDir, sprintf('PopPSTH_RF_%s.mat', AreaLabel));
save(SaveFile, ...
    'psth_targeton_raw_stable_RF_high',    'psth_targeton_raw_stable_RF_low', ...
    'psth_targeton_raw_stable_nonRF_high', 'psth_targeton_raw_stable_nonRF_low', ...
    'psth_targeton_raw_oneDR_RF_high',     'psth_targeton_raw_oneDR_RF_low', ...
    'psth_targeton_raw_oneDR_nonRF_high',  'psth_targeton_raw_oneDR_nonRF_low', ...
    'psth_RF_global_z_mu', 'psth_RF_global_z_sd', ...
    'ModuType_Sel', ...
    'PSTH_Time_TgtOn_stable', 'PSTH_Time_TgtOn_OneDR', '-v7.3');
fprintf('Saved PopPSTH RF -> %s  (%d neurons)\n', SaveFile, sum(Sel_Raw));
%{
% Save inputs for standalone value decoder script (Batch_ObjLocValue_Decoder.m)
SaveFile_Decoder = fullfile(SaveDir, sprintf('ValueDecoder_Input_%s.mat', AreaLabel));
save(SaveFile_Decoder, ...
    'psth_high_value_object_RF', 'psth_low_value_object_RF', ...
    'psth_high_value_location_RF', 'psth_low_value_location_RF', ...
    'time_axis_TDR', '-v7.3');
fprintf('Saved ValueDecoder inputs -> %s\n', SaveFile_Decoder);
%}
end








end


%% Local helper functions
function fig = plotAUCScatter(x, y, px, py, xLabel, yLabel, titleStr, isRectified,FigureStartNum)
% Plot scatter of tgtAUC values with marginal histograms.
% Marker scheme (all filled circles, p < 0.05):
%   red    : sig stable only  (px<0.05, py>=0.05)
%   purple : sig oneDR only   (px>=0.05, py<0.05)
%   orange : sig in both      (px<0.05, py<0.05)
%   open gray circle: neither significant
% Top histogram:   all (gray outline) + sig stable only + both (red/orange).
% Right histogram: all (gray outline) + sig oneDR only + both (purple/orange).
%
% INPUT:
%   x, y        - AUC values [1 x N] or [N x 1]
%   px, py      - p-values [1 x N] or [N x 1]
%   xLabel      - x-axis label
%   yLabel      - y-axis label
%   titleStr    - figure title
%   isRectified - logical; if true axes span [0.5, 1.0]; if false [0, 1]

x = x(:)'; y = y(:)'; px = px(:)'; py = py(:)';

sig_x = px < 0.05;
sig_y = py < 0.05;

grp_xonly = sig_x & ~sig_y;   % red
grp_yonly = ~sig_x & sig_y;   % purple
grp_both  = sig_x & sig_y;    % orange
grp_none  = ~sig_x & ~sig_y;  % open gray

if isRectified
    axLim     = [0.5, 1.0];
    histEdges = 0.5:0.02:1;
else
    axLim     = [0, 1];
    histEdges = 0:0.05:1;
end

mSz     = 50;
col_red = [0.85 0.15 0.15];
col_pur = [0.55 0.10 0.80];
col_org = [0.95 0.55 0.05];
alpha   = 0.6;

pos_sc  = [0.15  0.12  0.52  0.52];
pos_top = [0.15  0.67  0.52  0.20];
pos_rgt = [0.70  0.12  0.20  0.52];

fig=PrepareFigure(FigureStartNum,'w',[50,100, 800,800],'Name',titleStr);

set(fig, 'PaperUnits', 'inches');

% --- Scatter ---
ax_sc = axes('Position', pos_sc); hold on;

scatter(x(grp_none),  y(grp_none),  mSz, 'o', ...
    'MarkerEdgeColor', [0.6 0.6 0.6], 'MarkerFaceColor', 'none', 'LineWidth', 0.8);

h_red = scatter(x(grp_xonly), y(grp_xonly), mSz, 'o', 'filled', ...
    'MarkerEdgeColor', [0.6 0.6 0.6], 'MarkerFaceColor', col_red,'LineWidth', 0.8);
h_red.MarkerFaceAlpha = alpha;

h_pur = scatter(x(grp_yonly), y(grp_yonly), mSz, 'o', 'filled', ...
    'MarkerEdgeColor', [0.6 0.6 0.6], 'MarkerFaceColor', col_pur,'LineWidth', 0.8);
h_pur.MarkerFaceAlpha = alpha;

h_org = scatter(x(grp_both), y(grp_both), mSz, 'o', 'filled', ...
    'MarkerEdgeColor', [0.6 0.6 0.6], 'MarkerFaceColor', col_org,'LineWidth', 0.8);
h_org.MarkerFaceAlpha = alpha;

if ~isRectified
    plot(axLim, axLim, '--', 'Color', [0.6 0.6 0.6], 'LineWidth', 0.8);
end

xline(0.5, 'LineStyle', ':', 'Color', [0.6 0.6 0.6], 'LineWidth', 0.8);
yline(0.5, 'LineStyle', ':', 'Color', [0.6 0.6 0.6], 'LineWidth', 0.8);
xlim(axLim); ylim(axLim);
xlabel(xLabel); ylabel(yLabel);
title(titleStr);
set(ax_sc, 'TickDir', 'out', 'Box', 'off', 'PlotBoxAspectRatioMode', 'auto', 'DataAspectRatioMode', 'auto');

% --- Top histogram (stable AUC) ---
ax_top = axes('Position', pos_top); hold on;
histogram(x, histEdges, 'FaceColor', 'none', ...
    'EdgeColor', [0.7 0.7 0.7], 'LineWidth', 0.8);
h_tr = histogram(x(sig_x), histEdges, 'FaceColor', col_red, 'EdgeColor', 'none');
h_tr.FaceAlpha = alpha;
set(ax_top, 'XLim', axLim, 'XTickLabel', {}, 'TickDir', 'out', 'Box', 'off');
ylabel('Count');

% --- Right histogram (oneDR AUC, horizontal) ---
ax_rgt = axes('Position', pos_rgt); hold on;
histogram(y, histEdges, 'Orientation', 'horizontal', ...
    'FaceColor', 'none', 'EdgeColor', [0.7 0.7 0.7], 'LineWidth', 0.8);
h_pr = histogram(y(sig_y), histEdges, 'Orientation', 'horizontal', ...
    'FaceColor', col_pur, 'EdgeColor', 'none');
h_pr.FaceAlpha = alpha;
set(ax_rgt, 'YLim', axLim, 'YTickLabel', {}, 'TickDir', 'out', 'Box', 'off');
xlabel('Count');

% Link data limits for pan/zoom sync
linkaxes([ax_sc, ax_top], 'x');
linkaxes([ax_sc, ax_rgt], 'y');

% Apply formatting, then align histogram data areas to scatter data area.
% Each axes reserves different margins for tick/axis labels (TightInset).
% We compute the scatter's actual data area, then back-calculate histogram
% Positions so their data areas have the same left/width (top) and
% bottom/height (right).
formatFig(fig, [5 5], 'nature');
set(ax_sc,  'Position', pos_sc);
set(ax_top, 'Position', pos_top);
set(ax_rgt, 'Position', pos_rgt);
xlim(ax_sc,  axLim); ylim(ax_sc,  axLim);
xlim(ax_top, axLim);
ylim(ax_rgt, axLim);
drawnow;

ti_sc  = ax_sc.TightInset;   % [left bottom right top] label margins
data_l = pos_sc(1) + ti_sc(1);
data_b = pos_sc(2) + ti_sc(2);
data_w = pos_sc(3) - ti_sc(1) - ti_sc(3);
data_h = pos_sc(4) - ti_sc(2) - ti_sc(4);

ti_top = ax_top.TightInset;
set(ax_top, 'Position', [data_l - ti_top(1), pos_top(2), ...
                          data_w + ti_top(1) + ti_top(3), pos_top(4)]);

ti_rgt = ax_rgt.TightInset;
set(ax_rgt, 'Position', [pos_rgt(1), data_b - ti_rgt(2), ...
                          pos_rgt(3), data_h + ti_rgt(2) + ti_rgt(4)]);


end



function val = safe_nanmean(x, idx)
% Compute nanmean over indexed bins, returning NaN if x is too short
    if length(x) >= length(idx)
        val = nanmean(x(idx));
    else
        val = NaN;
    end
end

function plotSigBinsLocal(tAxis, mat, yPos, barH, clr)
% Plot a horizontal significance bar above a trace.
% Significant bins: ranksum of column values vs 0.5 (chance), p < 0.05.
%   tAxis — 1-by-T time vector
%   mat   — N-by-T matrix of AUC values
%   yPos  — y centre of bar
%   barH  — half-height of bar
%   clr   — [r g b] color
    nT = size(mat, 2);
    sigBins = false(1, nT);
    for t = 1:nT
        v = mat(:, t);
        v = v(~isnan(v));
        if length(v) >= 3
           sigBins(t) = signrank(v - 0.5) < 0.05;
           %sigBins(t) = ranksum(v,0.5) < 0.05;
        end
    end
    if ~any(sigBins), return; end
    dt     = median(diff(tAxis));
    sigIdx = find(sigBins);
    breaks = find(diff(sigIdx) > 1);
    starts = [sigIdx(1), sigIdx(breaks + 1)];
    ends   = [sigIdx(breaks), sigIdx(end)];
    for k = 1:length(starts)
        x1 = tAxis(starts(k)) - dt/2;
        x2 = tAxis(ends(k))   + dt/2;
        fill([x1 x2 x2 x1], [yPos-barH yPos-barH yPos+barH yPos+barH], ...
            clr, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    end
end
