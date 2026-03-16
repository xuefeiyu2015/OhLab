function Batch_StableOneDR_Opto(Data,  StartFile, EndFile,ShowFigureFlag,OutputFlag)
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






 


      %% %Stable Value
    % if isfield(OutputData,'StableValue')

          StableValueIndex(index) = 1;

          Data =OutputData.StableValue(1).DataStamp;


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

        
%
        PSTH_TgtOn_Prefer_RF_stable{index} = Data('PSTH_TgtOn_Prefer_RF_z');
        PSTH_TgtOn_nonPrefer_RF_stable{index} = Data('PSTH_TgtOn_nonPrefer_RF_z');
        PSTH_TgtOn_Prefer_nonRF_stable{index} = Data('PSTH_TgtOn_Prefer_nonRF_z');
        PSTH_TgtOn_nonPrefer_nonRF_stable{index} = Data('PSTH_TgtOn_nonPrefer_nonRF_z');

        PSTH_SacOn_Prefer_RF_stable{index} = Data('PSTH_SacOn_Prefer_RF_z');
        PSTH_SacOn_nonPrefer_RF_stable{index} = Data('PSTH_SacOn_nonPrefer_RF_z');
        PSTH_SacOn_Prefer_nonRF_stable{index} = Data('PSTH_SacOn_Prefer_nonRF_z');
        PSTH_SacOn_nonPrefer_nonRF_stable{index} = Data('PSTH_SacOn_nonPrefer_nonRF_z');

%

        PSTH_TgtOn_Good_RF_stable{index} = Data('PSTH_TgtOn_Good_RF_z');
        PSTH_TgtOn_Bad_RF_stable{index} = Data('PSTH_TgtOn_Bad_RF_z');
        PSTH_TgtOn_Good_nonRF_stable{index} = Data('PSTH_TgtOn_Good_nonRF_z');
        PSTH_TgtOn_Bad_nonRF_stable{index} = Data('PSTH_TgtOn_Bad_nonRF_z');

        PSTH_SacOn_Good_RF_stable{index} = Data('PSTH_SacOn_Good_RF_z');
        PSTH_SacOn_Bad_RF_stable{index} = Data('PSTH_SacOn_Bad_RF_z');
        PSTH_SacOn_Good_nonRF_stable{index} = Data('PSTH_SacOn_Good_nonRF_z');
        PSTH_SacOn_Bad_nonRF_stable{index} = Data('PSTH_SacOn_Bad_nonRF_z');


%}
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

       


        %{
        PSTH_TgtOn_Good_RF_stable{index} = Data('PSTH_TgtOn_Good_RF');
        PSTH_TgtOn_Bad_RF_stable{index} = Data('PSTH_TgtOn_Bad_RF');
        PSTH_TgtOn_Good_nonRF_stable{index} = Data('PSTH_TgtOn_Good_nonRF');
        PSTH_TgtOn_Bad_nonRF_stable{index} = Data('PSTH_TgtOn_Bad_nonRF');

        PSTH_SacOn_Good_RF_stable{index} = Data('PSTH_SacOn_Good_RF');
        PSTH_SacOn_Bad_RF_stable{index} = Data('PSTH_SacOn_Bad_RF');
        PSTH_SacOn_Good_nonRF_stable{index} = Data('PSTH_SacOn_Good_nonRF');
        PSTH_SacOn_Bad_nonRF_stable{index} = Data('PSTH_SacOn_Bad_nonRF');
        %}

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

        

        tmp = hPSTH_TgtOn_Value_RF{index};
        continuousBin = 5;
        Latency_Value_Stable(index)=NaN;
        tmp_time = NaN;

        for b = 1:length(hPSTH_TgtOn_Value_RF{index})
            if sum(tmp(b:min(b+continuousBin-1,length(tmp))))==continuousBin
                if isnan(Latency_Value_Stable(index)) 
                    if  PSTH_Time_TgtOn_stable(b)<0
                        tmp_time = PSTH_Time_TgtOn_stable(b);
                    else
                        Latency_Value_Stable(index) = PSTH_Time_TgtOn_stable(b);
                    end
                end
            end
        end

        if isnan(Latency_Value_Stable(index))
            Latency_Value_Stable(index)=tmp_time;
        end


 %Average responses from 100 to 300

        Good_stable_RF = psth_rawtgt_good_RF(PSTH_Time_TgtOn_stable>=100 & PSTH_Time_TgtOn_stable<=300);
        Bad_stable_RF = psth_rawtgt_bad_RF(PSTH_Time_TgtOn_stable>=100 & PSTH_Time_TgtOn_stable<=300);

        [h,p]=ttest2(Good_stable_RF,Bad_stable_RF);
        Value_Stable_RF(index) = h;
        Value_Stable_RF_Diff(index) = mean(Good_stable_RF)-mean(Bad_stable_RF);
        Value_Stable_RF_ROC(index) = rocN(Good_stable_RF,Bad_stable_RF);

        




        hPSTH_TgtOn_Good_RF{index}=Data('hPSTH_TgtOn_Good_RF');
        hPSTH_TgtOn_Bad_RF{index}=Data('hPSTH_TgtOn_Bad_RF');
        

        hPSTH_SacOn_Value_RF{index}=Data('hPSTH_SacOn_Value_RF');
        hPSTH_SacOn_Good_RF{index}=Data('hPSTH_SacOn_Good_RF');
        hPSTH_SacOn_Bad_RF{index}=Data('hPSTH_SacOn_Bad_RF');

        
        hPSTH_TgtOn_Value_nonRF{index}=Data('hPSTH_TgtOn_Value_nonRF');
        hPSTH_TgtOn_Good_nonRF{index}=Data('hPSTH_TgtOn_Good_nonRF');
        hPSTH_TgtOn_Bad_nonRF{index}=Data('hPSTH_TgtOn_Bad_nonRF');

        hPSTH_SacOn_Value_nonRF{index}=Data('hPSTH_SacOn_Value_nonRF');
        hPSTH_SacOn_Good_nonRF{index}=Data('hPSTH_SacOn_Good_nonRF');
        hPSTH_SacOn_Bad_nonRF{index}=Data('hPSTH_SacOn_Bad_nonRF');
        

%Select out neurons which shows significant responses
%Select out neurons which shows significant value difference from 100 to 300
SelectTime_Stable=[100,300];
 tmp = hPSTH_TgtOn_Value_RF{index};

hInterval = tmp(PSTH_Time_TgtOn_stable>=SelectTime_Stable(1) & PSTH_Time_TgtOn_stable<=SelectTime_Stable(2));
if sum(hInterval)>5
    Value_Stable_RF_Prop_Sig(index) = 1;
else
    Value_Stable_RF_Prop_Sig(index) = 0;

end

        

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


%{

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

     Tgt_Good_RF_oneDR(index) = Data('TgtGoodRF');
     Tgt_Bad_RF_oneDR(index) = Data('TgtBad_RF');   
     p_Tgt_RF_oneDR(index) = Data('p_Tgt_RF');






     PreTgtGood_RF_oneDR(index) = Data('PreTgtGood_RF');
     PreTgtBad_RF_oneDR(index) = Data('PreTgtBad_RF');
     p_preTgt_RF_oneDR(index) = Data('p_preTgt_RF');

     p_Sac_RF_oneDR(index) = Data('p_Sac_RF');

      Sac_Good_RF_oneDR(index) = Data('SacGood_RF');
     Sac_Bad_RF_oneDR(index) = Data('SacBad_RF');   
    
%å

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



%{
      PSTH_Good_TgtOn_RF_oneDR{index} = Data('PSTH_Good_TgtOn_RF');
     PSTH_Bad_TgtOn_RF_oneDR{index} = Data('PSTH_Bad_TgtOn_RF');
     PSTH_Good_TgtOn_nonRF_oneDR{index} = Data('PSTH_Good_TgtOn_nonRF');
     PSTH_Bad_TgtOn_nonRF_oneDR{index} = Data('PSTH_Bad_TgtOn_nonRF');

     PSTH_SacGood_nonRF_oneDR{index} = Data('PSTH_SacGood_nonRF');
     PSTH_SacBad_nonRF_oneDR{index} = Data('PSTH_SacBad_nonRF');
     PSTH_SacGood_RF_oneDR{index} = Data('PSTH_SacGood_RF');
     PSTH_SacBad_RF_oneDR{index} = Data('PSTH_SacBad_RF');
%}
     PSTH_Time_TgtOn_OneDR = Data('Time_TgtOn');
     PSTH_Time_SacOn_OneDR = Data('Time_SacOn');

     %Average responses from 150 to 350

        Good_location_RF = psth_rawtgt_good_RF(PSTH_Time_TgtOn_OneDR>=-100 & PSTH_Time_TgtOn_OneDR<=100);
        Bad_location_RF = psth_rawtgt_bad_RF(PSTH_Time_TgtOn_stable>=-100 & PSTH_Time_TgtOn_stable<=100);

        [h,p]=ttest2(Good_location_RF,Bad_location_RF);
        Value_Location_RF(index) = h;
        Value_Location_RF_Diff(index) = mean(Good_location_RF)-mean(Bad_location_RF);
        Value_Location_RF_ROC(index) = rocN(Good_location_RF,Bad_location_RF);

        AUCPSTH_TgtOn_Value_RF_oneDR{index} = Data('AUCPSTH_TgtOn_Value_RF');
        
        AUCPSTH_TgtOn_Value_nonRF_oneDR{index} = Data('AUCPSTH_TgtOn_Value_nonRF');
        

        AUCPSTH_SacOn_Value_RF_oneDR{index} = Data('AUCPSTH_SacOn_Value_RF');
        
        %AUCPSTH_SacOn_Values_nonRF_oneDR{index} = Data('AUCPSTH_SacOn_Value_nonRF');
        


     TgtGood_nonRF_oneDR(index) = Data('TgtGood_nonRF');
     TgtBad_nonRF_oneDR(index) = Data('TgtBad_nonRF');
     
     PreTgtGood_nonRF_oneDR(index) = Data('PreTgtGood_nonRF');
     PreTgtBad_nonRF_oneDR(index) = Data('PreTgtBad_nonRF');

     
     
     SacGood_nonRF_oneDR(index) = Data('SacGood_nonRF');
     SacBad_nonRF_oneDR(index) = Data('SacBad_nonRF');

     
     TgtAUC_nonRF_oneDR(index) = Data('TgtAUC_nonRF');

     

     SacAUC_nonRF_oneDR(index) = Data('SacAUC_nonRF');

     TgtAUC_RF_oneDR(index) = Data('TgtAUC_RF');
     SacAUC_RF_oneDR(index) = Data('SacAUC_RF');



     p_AUC_Tgt_RF_oneDR(index) = Data('p_AUC_Tgt_RF');
     p_AUC_Sac_RF_oneDR(index) = Data('p_AUC_Sac_RF');



     

     p_Tgt_nonRF_oneDR(index) = Data('p_Tgt_nonRF');

     p_preTgt_nonRF_oneDR(index) = Data('p_preTgt_nonRF');
     p_Sac_nonRF_oneDR(index) = Data('p_Sac_nonRF');
     p_AUC_Tgt_nonRF_oneDR(index) = Data('p_AUC_Tgt_nonRF');
     p_AUC_Sac_nonRF_oneDR(index) =Data('p_AUC_Sac_nonRF');

      p_Vis_RF_oneDR(index) = Data('p_Vis_RF');
        VisDiff_RF_oneDR(index) = Data('VisDiff_RF');
        p_Vis_nonRF_oneDR(index) = Data('p_Vis_nonRF');
        VisDiff_nonRF_oneDR(index) = Data('VisDiff_nonRF');

        p_pureSac_RF_oneDR(index) = Data('p_pureSac_RF');
        pureSacDiff_RF_oneDR(index) = Data('pureSacDiff_RF');
        p_pureSac_nonRF_oneDR(index) = Data('p_pureSac_nonRF');
        pureSacDiff_nonRF_oneDR(index) = Data('pureSacDiff_nonRF');

        hPSTH_TgtOn_Value_RF_oneDR{index}=Data('hPSTH_TgtOn_Value_RF');
        hPSTH_TgtOn_Good_RF_oneDR{index}=Data('hPSTH_TgtOn_Good_RF');
        hPSTH_TgtOn_Bad_RF_oneDR{index}=Data('hPSTH_TgtOn_Bad_RF');
        

        hPSTH_SacOn_Value_RF_oneDR{index}=Data('hPSTH_SacOn_Value_RF');
        hPSTH_SacOn_Good_RF_oneDR{index}=Data('hPSTH_SacOn_Good_RF');
        hPSTH_SacOn_Bad_RF_oneDR{index}=Data('hPSTH_SacOn_Bad_RF');

        
        hPSTH_TgtOn_Value_nonRF_oneDR{index}=Data('hPSTH_TgtOn_Value_nonRF');
        hPSTH_TgtOn_Good_nonRF_oneDR{index}=Data('hPSTH_TgtOn_Good_nonRF');
        hPSTH_TgtOn_Bad_nonRF_oneDR{index}=Data('hPSTH_TgtOn_Bad_nonRF');

        hPSTH_SacOn_Value_nonRF_oneDR{index}=Data('hPSTH_SacOn_Value_nonRF');
        hPSTH_SacOn_Good_nonRF_oneDR{index}=Data('hPSTH_SacOn_Good_nonRF');
        hPSTH_SacOn_Bad_nonRF_oneDR{index}=Data('hPSTH_SacOn_Bad_nonRF');

        tmp = hPSTH_TgtOn_Value_RF_oneDR{index};
        continuousBin = 5;
        Latency_Value_OneDR(index)=NaN;

        for b = 1:length(hPSTH_TgtOn_Value_RF_oneDR{index})
            if sum(tmp(b:min(b+continuousBin-1,length(tmp))))==continuousBin
                if isnan(Latency_Value_OneDR(index))
                    Latency_Value_OneDR(index) = PSTH_Time_TgtOn_OneDR(b);
                end
            end
        end
        

        SelectTime_Location=[-100,100];
        tmp = hPSTH_TgtOn_Value_RF_oneDR{index};

        

        hInterval = tmp(PSTH_Time_TgtOn_OneDR>=SelectTime_Location(1) & PSTH_Time_TgtOn_OneDR<=SelectTime_Location(2));
if sum(hInterval)>5
    Value_OneDR_RF_Prop_Sig(index) = 1;
else

    Value_OneDR_RF_Prop_Sig(index) = 0;
end

tmp_good = PSTH_Good_TgtOn_RF_oneDR{index};
tmp_bad = PSTH_Bad_TgtOn_RF_oneDR{index};
PSTH_TgtOn_Good_RF_oneDR_sel = tmp_good(PSTH_Time_TgtOn_OneDR>=SelectTime_Location(1) & PSTH_Time_TgtOn_OneDR<=SelectTime_Location(2));
PSTH_TgtOn_Bad_RF_oneDR_sel = tmp_bad(PSTH_Time_TgtOn_OneDR>=SelectTime_Location(1) & PSTH_Time_TgtOn_OneDR<=SelectTime_Location(2));

Mean_Value_Diff(index) = mean(PSTH_TgtOn_Good_RF_oneDR_sel)-mean(PSTH_TgtOn_Bad_RF_oneDR_sel);

        


     %Change into prefered psth
     Diff_TgtOn_curr_oneDR = Tgt_Good_RF_oneDR(index) - Tgt_Bad_RF_oneDR(index);
     Diff_preTgtOn_curr_oneDR = PreTgtGood_RF_oneDR(index) - PreTgtBad_RF_oneDR(index);

     if p_Tgt_RF_oneDR(index) < 0.05
         if Diff_TgtOn_curr_oneDR > 0
             PSTH_Prefer_TgtOn_RF_oneDR{index} = PSTH_Good_TgtOn_RF_oneDR{index};
             PSTH_nonPrefer_TgtOn_RF_oneDR{index} = PSTH_Bad_TgtOn_RF_oneDR{index};

             PSTH_Prefer_SacOn_RF_oneDR{index} = PSTH_SacGood_RF_oneDR{index};
             PSTH_nonPrefer_SacOn_RF_oneDR{index} = PSTH_SacBad_RF_oneDR{index};

             PSTH_Prefer_TgtOn_nonRF_oneDR{index} = PSTH_Good_TgtOn_nonRF_oneDR{index};
             PSTH_nonPrefer_TgtOn_nonRF_oneDR{index} = PSTH_Bad_TgtOn_nonRF_oneDR{index};

             PSTH_Prefer_SacOn_nonRF_oneDR{index} = PSTH_SacGood_nonRF_oneDR{index};
             PSTH_nonPrefer_SacOn_nonRF_oneDR{index} = PSTH_SacBad_nonRF_oneDR{index};

             




         else
             PSTH_Prefer_TgtOn_RF_oneDR{index} = PSTH_Bad_TgtOn_RF_oneDR{index};
             PSTH_nonPrefer_TgtOn_RF_oneDR{index} = PSTH_Good_TgtOn_RF_oneDR{index};

             PSTH_Prefer_SacOn_RF_oneDR{index} = PSTH_SacBad_RF_oneDR{index};
             PSTH_nonPrefer_SacOn_RF_oneDR{index} = PSTH_SacGood_RF_oneDR{index};

             PSTH_Prefer_TgtOn_nonRF_oneDR{index} = PSTH_Bad_TgtOn_nonRF_oneDR{index};
             PSTH_nonPrefer_TgtOn_nonRF_oneDR{index} = PSTH_Good_TgtOn_nonRF_oneDR{index};

             PSTH_Prefer_SacOn_nonRF_oneDR{index} = PSTH_SacBad_nonRF_oneDR{index};
             PSTH_nonPrefer_SacOn_nonRF_oneDR{index} = PSTH_SacGood_nonRF_oneDR{index};



         end
     elseif p_preTgt_RF_oneDR(index) < 0.05
         if Diff_preTgtOn_curr_oneDR > 0
             PSTH_Prefer_TgtOn_RF_oneDR{index} = PSTH_Good_TgtOn_RF_oneDR{index};
             PSTH_nonPrefer_TgtOn_RF_oneDR{index} = PSTH_Bad_TgtOn_RF_oneDR{index};

             PSTH_Prefer_SacOn_RF_oneDR{index} = PSTH_SacGood_RF_oneDR{index};
             PSTH_nonPrefer_SacOn_RF_oneDR{index} = PSTH_SacBad_RF_oneDR{index};

             PSTH_Prefer_TgtOn_nonRF_oneDR{index} = PSTH_Good_TgtOn_nonRF_oneDR{index};
             PSTH_nonPrefer_TgtOn_nonRF_oneDR{index} = PSTH_Bad_TgtOn_nonRF_oneDR{index};

             PSTH_Prefer_SacOn_nonRF_oneDR{index} = PSTH_SacGood_nonRF_oneDR{index};
             PSTH_nonPrefer_SacOn_nonRF_oneDR{index} = PSTH_SacBad_nonRF_oneDR{index};


         else
             PSTH_Prefer_TgtOn_RF_oneDR{index} = PSTH_Bad_TgtOn_RF_oneDR{index};
             PSTH_nonPrefer_TgtOn_RF_oneDR{index} = PSTH_Good_TgtOn_RF_oneDR{index};

             PSTH_Prefer_SacOn_RF_oneDR{index} = PSTH_SacBad_RF_oneDR{index};
             PSTH_nonPrefer_SacOn_RF_oneDR{index} = PSTH_SacGood_RF_oneDR{index};

             PSTH_Prefer_TgtOn_nonRF_oneDR{index} = PSTH_Bad_TgtOn_nonRF_oneDR{index};
             PSTH_nonPrefer_TgtOn_nonRF_oneDR{index} = PSTH_Good_TgtOn_nonRF_oneDR{index};

             PSTH_Prefer_SacOn_nonRF_oneDR{index} = PSTH_SacBad_nonRF_oneDR{index};
             PSTH_nonPrefer_SacOn_nonRF_oneDR{index} = PSTH_SacGood_nonRF_oneDR{index};



         end
     else
         if Diff_preTgtOn_curr_oneDR > 0
            PSTH_Prefer_TgtOn_RF_oneDR{index} = PSTH_Good_TgtOn_RF_oneDR{index};
            PSTH_nonPrefer_TgtOn_RF_oneDR{index} = PSTH_Bad_TgtOn_RF_oneDR{index};

            PSTH_Prefer_SacOn_RF_oneDR{index} = PSTH_SacGood_RF_oneDR{index};
            PSTH_nonPrefer_SacOn_RF_oneDR{index} = PSTH_SacBad_RF_oneDR{index};

            PSTH_Prefer_TgtOn_nonRF_oneDR{index} = PSTH_Good_TgtOn_nonRF_oneDR{index};
            PSTH_nonPrefer_TgtOn_nonRF_oneDR{index} = PSTH_Bad_TgtOn_nonRF_oneDR{index};

            PSTH_Prefer_SacOn_nonRF_oneDR{index} = PSTH_SacGood_nonRF_oneDR{index};
            PSTH_nonPrefer_SacOn_nonRF_oneDR{index} = PSTH_SacBad_nonRF_oneDR{index};
         else
             PSTH_Prefer_TgtOn_RF_oneDR{index} = PSTH_Bad_TgtOn_RF_oneDR{index};
             PSTH_nonPrefer_TgtOn_RF_oneDR{index} = PSTH_Good_TgtOn_RF_oneDR{index};

             PSTH_Prefer_SacOn_RF_oneDR{index} = PSTH_SacBad_RF_oneDR{index};
             PSTH_nonPrefer_SacOn_RF_oneDR{index} = PSTH_SacGood_RF_oneDR{index};

             PSTH_Prefer_TgtOn_nonRF_oneDR{index} = PSTH_Bad_TgtOn_nonRF_oneDR{index};
             PSTH_nonPrefer_TgtOn_nonRF_oneDR{index} = PSTH_Good_TgtOn_nonRF_oneDR{index};

             PSTH_Prefer_SacOn_nonRF_oneDR{index} = PSTH_SacBad_nonRF_oneDR{index};
             PSTH_nonPrefer_SacOn_nonRF_oneDR{index} = PSTH_SacGood_nonRF_oneDR{index};
             
         end
     end

       Value_Vis_RF_OneDR(index) = sign(Tgt_Good_RF_oneDR(index)-Tgt_Bad_RF_oneDR(index));

        %Value_Vis_nonRF_OneDR(index) = sign(Tgt_Good_nonRF_oneDR(index)-Tgt_Bad_nonRF_oneDR(index));
    Value_Vis_nonRF_OneDR(index) =NaN;

        Value_Sac_RF_OneDR(index) = sign(Sac_Good_RF_oneDR(index)-Sac_Bad_RF_oneDR(index));

     %   Value_Sac_nonRF_OneDR(index) = sign(Sac_Good_nonRF(index)-Sac_Bad_nonRF(index));
     Value_Sac_nonRF_OneDR(index) = NaN;

       

        
else
        OneDRIndex(index) = 0;

        Value_Vis_RF_OneDR(index) = NaN;

        Value_Vis_nonRF_OneDR(index) = NaN;

        Value_Sac_RF_OneDR(index) =NaN;

        Value_Sac_nonRF_OneDR(index) =NaN;


     Tgt_Good_RF_oneDR(index) = NaN;
     Tgt_Bad_RF_oneDR(index) = NaN;   
     p_Tgt_RF_oneDR(index) =NaN;






     PreTgtGood_RF_oneDR(index) = NaN;
     PreTgtBad_RF_oneDR(index) =NaN;
     p_preTgt_RF_oneDR(index) = NaN;

     p_Sac_RF_oneDR(index) = NaN;

      Sac_Good_RF_oneDR(index) = NaN;
     Sac_Bad_RF_oneDR(index) =NaN;   
    


      PSTH_Good_TgtOn_RF_oneDR{index} =NaN*ones(1,160);
     PSTH_Bad_TgtOn_RF_oneDR{index} = NaN*ones(1,160);
     PSTH_Good_TgtOn_nonRF_oneDR{index} = NaN*ones(1,160);
     PSTH_Bad_TgtOn_nonRF_oneDR{index} = NaN*ones(1,160);

     AUCPSTH_TgtOn_Value_RF_oneDR{index} =NaN*ones(1,160);
        
        AUCPSTH_TgtOn_Value_nonRF_oneDR{index} = NaN*ones(1,160);
        

        AUCPSTH_SacOn_Value_RF_oneDR{index} =NaN*ones(1,160);
        
        AUCPSTH_SacOn_Values_nonRF_oneDR{index} = NaN*ones(1,160);
        

     PSTH_SacGood_nonRF_oneDR{index} = NaN;
     PSTH_SacBad_nonRF_oneDR{index} = NaN;
     PSTH_SacGood_RF_oneDR{index} = NaN;
     PSTH_SacBad_RF_oneDR{index} = NaN;


    


     TgtGood_nonRF_oneDR(index) = NaN;
     TgtBad_nonRF_oneDR(index) = NaN;
     
     PreTgtGood_nonRF_oneDR(index) = NaN;
     PreTgtBad_nonRF_oneDR(index) = NaN;

     
     
     SacGood_nonRF_oneDR(index) = NaN;
     SacBad_nonRF_oneDR(index) = NaN;

     
     TgtAUC_nonRF_oneDR(index) = NaN;

     

     SacAUC_nonRF_oneDR(index) = NaN;

     TgtAUC_RF_oneDR(index) = NaN;
     SacAUC_RF_oneDR(index) = NaN;



     p_AUC_Tgt_RF_oneDR(index) = NaN;
     p_AUC_Sac_RF_oneDR(index) = NaN;



     

     p_Tgt_nonRF_oneDR(index) = NaN;

     p_preTgt_nonRF_oneDR(index) = NaN;
     p_Sac_nonRF_oneDR(index) = NaN;
     p_AUC_Tgt_nonRF_oneDR(index) = NaN;
     p_AUC_Sac_nonRF_oneDR(index) = NaN;



     

    
             PSTH_Prefer_TgtOn_RF_oneDR{index} = NaN;
             PSTH_nonPrefer_TgtOn_RF_oneDR{index} = NaN;

             PSTH_Prefer_SacOn_RF_oneDR{index} = NaN;
             PSTH_nonPrefer_SacOn_RF_oneDR{index} = NaN;

             PSTH_Prefer_TgtOn_nonRF_oneDR{index} = NaN;
             PSTH_nonPrefer_TgtOn_nonRF_oneDR{index} = NaN;

             PSTH_Prefer_SacOn_nonRF_oneDR{index} = NaN;
             PSTH_nonPrefer_SacOn_nonRF_oneDR{index} = NaN;

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





%Selection 
%{
%ManuelOut
CellIndex = 1:length(FR_Mean_Video);
%ManuelCheck=[550 469	478	477	473	479	425	408	353	351	320	239	240	217	213	201	178	163	160	139	141	68	43	48 115 132 137 157 181 200 221 496 516 529 538 581 623 97 131 159 315 325 329 343 350 400 427 435];
ManuelCheck=[];
ManuelOut = ~ismember(CellIndex,ManuelCheck);


Cri1 = FR_Mean_Video>1;



ScreenedIndex = CellIndex(Cri1 & ManuelOut);

ScreenCell = ismember(CellIndex,ScreenedIndex);

ScreenCellNo = CellIndex(ScreenCell);
%}

%% For VideoFreeView

%{

FilesNameScreened =FilesNameAll(ScreenCell);

SM = SM(ScreenCell);
SM_p= SM_p(ScreenCell);

Latency = Latency(ScreenCell);
%}
%ChannelScreened = ChannelScreened(ScreenCell);
%CellType = CellType(ScreenCell);

%FR_Mean_Video = FR_Mean_Video(ScreenCell);

%MemorySaccadeIndex = MemorySaccadeIndex(ScreenCell);

%DelaySaccadeIndex = DelaySaccadeIndex(ScreenCell);

%Prefer_TgtVector = Prefer_TgtVector(ScreenCell);

%p_TgtV_prefer = p_TgtV_prefer(ScreenCell);

%TuningStrength =  TuningStrength(ScreenCell);
%{
AS =  AS(ScreenCell)';

p_Cri20 = p_Cri20(ScreenCell);

p_Cri50 = p_Cri50(ScreenCell);

Latency2 =  Latency2(ScreenCell);

AbsoluteIndex = AbsoluteIndex(ScreenCell);

p_whole_cri50= p_whole_cri50(ScreenCell);


PMS_p = PMS_p(ScreenCell);

PSTH = cell2mat(PSTH_Stim_Mean(ScreenCell)'); 
PSTH_control = cell2mat(PSTH_Control_Mean(ScreenCell)');

PSTH_norm = cell2mat(PSTH_Stim_Norm(ScreenCell)'); 
PSTH_control_norm = cell2mat(PSTH_Control_Norm(ScreenCell)');


PSTH_norm2 =cell2mat(PSTH_Stim_Norm2(ScreenCell)');



PSTH_z = cell2mat(PSTH_Stim_Mean_z(ScreenCell)');


p_whole_cri50 = cell2mat(p_whole_cri50);
%}
%maxnum_ms = max(cellfun(@numel,PSTH_Tgt_RF));

%% For Stable Value

 %{  
ValueTaskInfo = [p_Tgt_RF_stable', p_Sac_RF_stable',TgtGood_RF_Stable',TgtBad_RF_Stable',TgtGood_nonRF_Stable',TgtBad_nonRF_Stable',...
    SacGood_RF_Stable',SacBad_RF_Stable',SacGood_nonRF_Stable',SacBad_nonRF_Stable',Value_Vis_RF',Value_Vis_nonRF',Value_Sac_RF', Value_Sac_nonRF',StableValueIndex'];

ValueTaskInfo_Screened = ValueTaskInfo(ScreenCell,:);

%p_Tgt_RF_stable = ValueTaskInfo_Screened(:,1);
%p_Sac_RF_stable = ValueTaskInfo_Screened(:,2);
TgtGood_RF_Stable = ValueTaskInfo_Screened(:,3);
TgtBad_RF_Stable = ValueTaskInfo_Screened(:,4);
TgtGood_nonRF_Stable = ValueTaskInfo_Screened(:,5);
TgtBad_nonRF_Stable = ValueTaskInfo_Screened(:,6);
%SacGood_RF_Stable = ValueTaskInfo_Screened(:,7);
%SacBad_RF_Stable = ValueTaskInfo_Screened(:,8);
%SacGood_nonRF_Stable = ValueTaskInfo_Screened(:,9);
%SacBad_nonRF_Stable = ValueTaskInfo_Screened(:,10);
%Value_Vis_RF = ValueTaskInfo_Screened(:,11);
%Value_Vis_nonRF = ValueTaskInfo_Screened(:,12);
%Value_Sac_RF = ValueTaskInfo_Screened(:,13);
%Value_Sac_nonRF = ValueTaskInfo_Screened(:,14);
StableValueIndex = ValueTaskInfo_Screened(:,15);

%p_Sac_nonRF_stable=p_Sac_nonRF_stable(ScreenCell)';



 %AUC_Stable_RF =  AUC_Stable_RF(ScreenCell);
 %AUC_Stable_nonRF = AUC_Stable_nonRF(ScreenCell);

 %pAUC_Stable_RF = pAUC_Stable_RF(ScreenCell);
 %pAUC_Stable_nonRF = pAUC_Stable_nonRF(ScreenCell);

 


 %PSTH_TgtOn_Prefer_RF_stable=PSTH_TgtOn_Prefer_RF_stable(ScreenCell) ;
 %PSTH_TgtOn_nonPrefer_RF_stable = PSTH_TgtOn_nonPrefer_RF_stable(ScreenCell) ;
 %PSTH_TgtOn_Prefer_nonRF_stable = PSTH_TgtOn_Prefer_nonRF_stable(ScreenCell);
 %PSTH_TgtOn_nonPrefer_nonRF_stable = PSTH_TgtOn_nonPrefer_nonRF_stable(ScreenCell);

 %PSTH_SacOn_Prefer_RF_stable = PSTH_SacOn_Prefer_RF_stable(ScreenCell);
 %PSTH_SacOn_nonPrefer_RF_stable = PSTH_SacOn_nonPrefer_RF_stable(ScreenCell);
 %PSTH_SacOn_Prefer_nonRF_stable= PSTH_SacOn_Prefer_nonRF_stable(ScreenCell);
 %PSTH_SacOn_nonPrefer_nonRF_stable = PSTH_SacOn_nonPrefer_nonRF_stable(ScreenCell);

 


PSTH_TgtOn_Good_RF_stable =PSTH_TgtOn_Good_RF_stable(ScreenCell);
PSTH_TgtOn_Bad_RF_stable = PSTH_TgtOn_Bad_RF_stable(ScreenCell) ;
PSTH_TgtOn_Good_nonRF_stable= PSTH_TgtOn_Good_nonRF_stable(ScreenCell);
PSTH_TgtOn_Bad_nonRF_stable= PSTH_TgtOn_Bad_nonRF_stable(ScreenCell);

%PSTH_SacOn_Good_RF_stable = PSTH_SacOn_Good_RF_stable(ScreenCell) ;
%PSTH_SacOn_Bad_RF_stable = PSTH_SacOn_Bad_RF_stable(ScreenCell);
%PSTH_SacOn_Good_nonRF_stable = PSTH_SacOn_Good_nonRF_stable(ScreenCell) ;
%PSTH_SacOn_Bad_nonRF_stable = PSTH_SacOn_Bad_nonRF_stable(ScreenCell);



 %}






%{
%% One DR
OneDRTaskInfo = [p_Tgt_RF', p_Sac_RF', p_preTgt_RF',Tgt_Good_RF',Tgt_Bad_RF',PreTgtGood_RF',PreTgtBad_RF',OneDRIndex',Value_Vis_RF_OneDR',Value_Vis_nonRF_OneDR',Value_Sac_RF_OneDR',...
    Value_Sac_nonRF_OneDR'   ];

OneDRTaskScreened = OneDRTaskInfo(ScreenCell,:);

%p_Tgt_RF = OneDRTaskScreened(:,1);
%p_Sac_RF = OneDRTaskScreened(:,2);
%p_preTgt_RF = OneDRTaskScreened(:,3);
Tgt_Good_RF = OneDRTaskScreened(:,4);
Tgt_Bad_RF = OneDRTaskScreened(:,5);
%PreTgtGood_RF= OneDRTaskScreened(:,6);
%PreTgtBad_RF= OneDRTaskScreened(:,7);



OneDRIndex = OneDRTaskScreened(:,8);

TgtGood_nonRF = TgtGood_nonRF(ScreenCell)';
TgtBad_nonRF =TgtBad_nonRF(ScreenCell)';

%Value_Vis_RF_OneDR = OneDRTaskScreened(:,9);
%Value_Vis_nonRF_OneDR = OneDRTaskScreened(:,10);
%Value_Sac_RF_OneDR = OneDRTaskScreened(:,11);
%Value_Sac_nonRF_OneDR = OneDRTaskScreened(:,12);


%PSTH_Prefer_TgtOn_RF = PSTH_Good_TgtOn_RF(ScreenCell);
%PSTH_nonPrefer_TgtOn_RF = PSTH_Bad_TgtOn_RF(ScreenCell);

%PSTH_Prefer_SacOn_RF= PSTH_SacGood_RF(ScreenCell);
%PSTH_nonPrefer_SacOn_RF= PSTH_SacBad_RF(ScreenCell);

%PSTH_Prefer_TgtOn_nonRF= PSTH_Good_TgtOn_nonRF(ScreenCell);
%PSTH_nonPrefer_TgtOn_nonRF = PSTH_Bad_TgtOn_nonRF(ScreenCell);

%PSTH_Prefer_SacOn_nonRF = PSTH_SacGood_nonRF(ScreenCell);
%PSTH_nonPrefer_SacOn_nonRF= PSTH_SacBad_nonRF(ScreenCell);

PSTH_Good_TgtOn_RF = PSTH_Good_TgtOn_RF(ScreenCell);
PSTH_Bad_TgtOn_RF = PSTH_Bad_TgtOn_RF(ScreenCell);

PSTH_Good_TgtOn_nonRF= PSTH_Good_TgtOn_nonRF(ScreenCell);
PSTH_Bad_TgtOn_nonRF = PSTH_Bad_TgtOn_nonRF(ScreenCell);








%p_preTgt_nonRF = p_preTgt_nonRF(ScreenCell)';
%p_Sac_nonRF = p_Sac_nonRF(ScreenCell)';
%p_AUC_Tgt_nonRF = p_AUC_Tgt_nonRF(ScreenCell)';
%p_AUC_Sac_nonRF =p_AUC_Sac_nonRF(ScreenCell)';
%p_Tgt_nonRF = p_Tgt_nonRF(ScreenCell)';

%TgtAUC_RF = TgtAUC_RF(ScreenCell)';
 %SacAUC_RF = SacAUC_RF(ScreenCell)';

 %TgtAUC_NonRF = TgtAUC_nonRF(ScreenCell)';
 %SacAUC_NonRF = SacAUC_nonRF(ScreenCell)';

 %p_AUC_Tgt_RF = p_AUC_Tgt_RF(ScreenCell)';
 %p_AUC_Sac_RF = p_AUC_Sac_RF(ScreenCell)';




%p_Response_interval_RF = p_Response_interval_RF(ScreenCell,:);
%TgtAUC_interval_RF = TgtAUC_interval_RF(ScreenCell,:);
% p_Response_interval_nonRF = p_Response_interval_nonRF(ScreenCell,:);
%TgtAUC_interval_nonRF = TgtAUC_interval_nonRF(ScreenCell,:);

%}
%%
%Separate cells into different modulation groups according to laser
%response


%Separate into 5 groups 


Fast_M = Latency<10;
Slow_M = Latency>=10;
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


%In_Fast_Sig = Fast_M & Inhibition & Sig_m;

%In_Slow_Sig = Slow_M & Inhibition & Sig_m;

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


%{
%For bar plot
ExcitationNum = sum(Excitation & Sig_m);
InhibitionNum = sum(Inhibition & Sig_m);
NonSigNum = sum(~Sig_m);
%}
%Cell Selection
ResponseValid = (FR_Mean_oneDR>1 & FR_Mean_stable>1)';

Sel = OneDRIndex' == 1;
Sel_All = Sel & OptoIndex' ==1 & ResponseValid;
%Sel_All = OptoIndex' ==1;
Sel_Others = Sel_All & ModuType==4 ;
Sel_Fast = Sel_All & ModuType==1 ;


PSTH_TgtOn_Good_RF_Stable = cell2mat(PSTH_TgtOn_Good_RF_stable');
PSTH_TgtOn_Bad_RF_Stable = cell2mat(PSTH_TgtOn_Bad_RF_stable');

PSTH_TgtOn_Good_RF_Stable=(PSTH_TgtOn_Good_RF_Stable')';
PSTH_TgtOn_Bad_RF_Stable=(PSTH_TgtOn_Bad_RF_Stable')';


PSTH_TgtOn_Diff_RF_Stable = PSTH_TgtOn_Good_RF_Stable -PSTH_TgtOn_Bad_RF_Stable;

PSTH_SacOn_Good_RF_Stable = cell2mat(PSTH_SacOn_Good_RF_stable');
PSTH_SacOn_Bad_RF_Stable = cell2mat(PSTH_SacOn_Bad_RF_stable');

PSTH_SacOn_Diff_RF_Stable = PSTH_SacOn_Good_RF_Stable -PSTH_SacOn_Bad_RF_Stable;

[value,sorted_index]=sort(Latency_Value_Stable,'Ascend');
%[value,sorted_index]=sort(Latency_Value_OneDR,'Ascend');
Latency_Sel_Fast_Stable = Latency_Value_Stable(Sel_Fast);
Latency_Sel_Fast_OneDR = Latency_Value_OneDR(Sel_Fast);

Latency_Sel_Others_Stable = Latency_Value_Stable(Sel_Others);
Latency_Sel_Others_OneDR = Latency_Value_OneDR(Sel_Others);

Latency_Sel_All_Stable = Latency_Value_Stable(Sel_All);
Latency_Sel_All_OneDR = Latency_Value_OneDR(Sel_All);
%{
NumberValue_Fast_Stable = sum(~isnan(Latency_Sel_Fast_Stable));
PropValue_Fast_Stable = NumberValue_Fast_Stable/length(Latency_Sel_Fast_Stable);

NumberValue_Fast_OneDR = sum(~isnan(Latency_Sel_Fast_OneDR));
PropValue_Fast_OneDR = NumberValue_Fast_OneDR/length(Latency_Sel_Fast_OneDR);

NumberValue_Other_Stable = sum(~isnan(Latency_Sel_Others_Stable));
PropValue_Other_Stable = NumberValue_Other_Stable/length(Latency_Sel_Others_Stable);

NumberValue_Other_OneDR = sum(~isnan(Latency_Sel_Others_OneDR));
PropValue_Other_OneDR = NumberValue_Other_OneDR/length(Latency_Sel_Others_OneDR);
%}

%Distribution of time of value separation for whole FEF population
bin_edges = -400:50:400;
counts_Stable_All = histcounts(Latency_Sel_All_Stable, bin_edges);
counts_OneDR_All = histcounts(Latency_Sel_All_OneDR, bin_edges);



%keyboard

%{
PSTH_TgtOn_Good_RF_Stable_early = PSTH_TgtOn_Good_RF_Stable(Latency_Value_Stable<0,:);
PSTH_TgtOn_Bad_RF_Stable_early = PSTH_TgtOn_Bad_RF_Stable(Latency_Value_Stable<0,:);

number = size(PSTH_TgtOn_Good_RF_Stable_early,1);
%}
%{
figure
for i = 1:number
    if mod(i,16)==1
        figure
        index = 1;
    end
    subplot(4,4,index)
    plot(PSTH_Time_TgtOn_stable,PSTH_TgtOn_Good_RF_Stable_early(i,:),'-r');
    hold on
    plot(PSTH_Time_TgtOn_stable,PSTH_TgtOn_Bad_RF_Stable_early(i,:),'-b');
    index = index+1;
    

end
%}
%{
figure
subplot(2,1,1)
histogram(Latency_Sel_Fast_Stable)
subplot(2,1,2)
histogram(Latency_Sel_Fast_OneDR)
%}
PSTH_Good_TgtOn_RF_oneDR = (cell2mat(PSTH_Good_TgtOn_RF_oneDR')')';
PSTH_Bad_TgtOn_RF_oneDR = (cell2mat(PSTH_Bad_TgtOn_RF_oneDR')')';

PSTH_TgtOn_Diff_RF_oneDR = PSTH_Good_TgtOn_RF_oneDR - PSTH_Bad_TgtOn_RF_oneDR;


%% Proportions

%% Fast group
hPSTH_TgtOn_Value_RF_sel_Fast = cell2mat(hPSTH_TgtOn_Value_RF(Sel_Fast)');%(Sel_All);

Prop_TgtOn_Value_RF_Fast = sum(hPSTH_TgtOn_Value_RF_sel_Fast,1)/size(hPSTH_TgtOn_Value_RF_sel_Fast,1);

hPSTH_SacOn_Value_RF_sel_Fast = cell2mat(hPSTH_SacOn_Value_RF(Sel_Fast)');%(Sel_All);

Prop_SacOn_Value_RF_Fast = sum(hPSTH_SacOn_Value_RF_sel_Fast,1)/size(hPSTH_SacOn_Value_RF_sel_Fast,1);


hPSTH_TgtOn_Value_nonRF{2}=logical(zeros(1,160));
hPSTH_TgtOn_Value_nonRF{3}=logical(zeros(1,160));

hPSTH_SacOn_Value_nonRF{2}=logical(zeros(1,120));
hPSTH_SacOn_Value_nonRF{3}=logical(zeros(1,120));


hPSTH_TgtOn_Value_nonRF_sel_Fast = cell2mat(hPSTH_TgtOn_Value_nonRF(Sel_Fast)');%(Sel_All);

Prop_TgtOn_Value_nonRF_Fast = sum(hPSTH_TgtOn_Value_nonRF_sel_Fast,1)/size(hPSTH_TgtOn_Value_nonRF_sel_Fast,1);

hPSTH_SacOn_Value_nonRF_sel_Fast = cell2mat(hPSTH_SacOn_Value_nonRF(Sel_Fast)');%(Sel_All);

Prop_SacOn_Value_nonRF_Fast = sum(hPSTH_SacOn_Value_nonRF_sel_Fast,1)/size(hPSTH_SacOn_Value_nonRF_sel_Fast,1);

%% Other group

hPSTH_TgtOn_Value_RF_sel_Others = cell2mat(hPSTH_TgtOn_Value_RF(Sel_Others)');%(Sel_All);

Prop_TgtOn_Value_RF_Others = sum(hPSTH_TgtOn_Value_RF_sel_Others,1)/size(hPSTH_TgtOn_Value_RF_sel_Others,1);

hPSTH_SacOn_Value_RF_sel_Others = cell2mat(hPSTH_SacOn_Value_RF(Sel_Others)');%(Sel_All);

Prop_SacOn_Value_RF_Others = sum(hPSTH_SacOn_Value_RF_sel_Others,1)/size(hPSTH_SacOn_Value_RF_sel_Others,1);



hPSTH_TgtOn_Value_nonRF_sel_Others = cell2mat(hPSTH_TgtOn_Value_nonRF(Sel_Others)');%(Sel_All);

Prop_TgtOn_Value_nonRF_Others = sum(hPSTH_TgtOn_Value_nonRF_sel_Others,1)/size(hPSTH_TgtOn_Value_nonRF_sel_Others,1);

hPSTH_SacOn_Value_nonRF_sel_Others = cell2mat(hPSTH_SacOn_Value_nonRF(Sel_Others)');%(Sel_All);

Prop_SacOn_Value_nonRF_Others = sum(hPSTH_SacOn_Value_nonRF_sel_Others,1)/size(hPSTH_SacOn_Value_nonRF_sel_Others,1);


%Whole group 
hPSTH_TgtOn_Value_RF_sel_All = cell2mat(hPSTH_TgtOn_Value_RF(Sel_All)');%(Sel_All);

Prop_TgtOn_Value_RF_All = sum(hPSTH_TgtOn_Value_RF_sel_All,1)/size(hPSTH_TgtOn_Value_RF_sel_All,1);

hPSTH_SacOn_Value_RF_sel_All = cell2mat(hPSTH_SacOn_Value_RF(Sel_All)');%(Sel_All);

Prop_SacOn_Value_RF_All = sum(hPSTH_SacOn_Value_RF_sel_All,1)/size(hPSTH_SacOn_Value_RF_sel_All,1);



hPSTH_TgtOn_Value_nonRF_sel_All = cell2mat(hPSTH_TgtOn_Value_nonRF(Sel_All)');%(Sel_All);

Prop_TgtOn_Value_nonRF_All = sum(hPSTH_TgtOn_Value_nonRF_sel_All,1)/size(hPSTH_TgtOn_Value_nonRF_sel_All,1);

hPSTH_SacOn_Value_nonRF_sel_All = cell2mat(hPSTH_SacOn_Value_nonRF(Sel_All)');%(Sel_All);

Prop_SacOn_Value_nonRF_All = sum(hPSTH_SacOn_Value_nonRF_sel_All,1)/size(hPSTH_SacOn_Value_nonRF_sel_All,1);






%Good vs Bad
hPSTH_TgtOn_Good_RF_sel = cell2mat(hPSTH_TgtOn_Good_RF(Sel_Fast)');

Prop_TgtOn_Good_RF = sum(hPSTH_TgtOn_Good_RF_sel,1)/size(hPSTH_TgtOn_Good_RF_sel,1);

hPSTH_SacOn_Good_RF_sel = cell2mat(hPSTH_SacOn_Good_RF(Sel_Fast)');%(Sel_All);

Prop_SacOn_Good_RF = sum(hPSTH_SacOn_Good_RF_sel,1)/size(hPSTH_SacOn_Good_RF_sel,1);

hPSTH_TgtOn_Bad_RF_sel = cell2mat(hPSTH_TgtOn_Bad_RF(Sel_Fast)');

Prop_TgtOn_Bad_RF = sum(hPSTH_TgtOn_Bad_RF_sel,1)/size(hPSTH_TgtOn_Bad_RF_sel,1);

hPSTH_SacOn_Bad_RF_sel = cell2mat(hPSTH_SacOn_Bad_RF(Sel_Fast)');%(Sel_All);

Prop_SacOn_Bad_RF = sum(hPSTH_SacOn_Bad_RF_sel,1)/size(hPSTH_SacOn_Bad_RF_sel,1);


hPSTH_TgtOn_Good_nonRF{2}=logical(zeros(1,160));
hPSTH_TgtOn_Good_nonRF{3}=logical(zeros(1,160));

hPSTH_SacOn_Good_nonRF{2}=logical(zeros(1,120));
hPSTH_SacOn_Good_nonRF{3}=logical(zeros(1,120));

hPSTH_TgtOn_Bad_nonRF{2}=logical(zeros(1,160));
hPSTH_TgtOn_Bad_nonRF{3}=logical(zeros(1,160));

hPSTH_SacOn_Bad_nonRF{2}=logical(zeros(1,120));
hPSTH_SacOn_Bad_nonRF{3}=logical(zeros(1,120));

hPSTH_TgtOn_Good_nonRF_sel = cell2mat(hPSTH_TgtOn_Good_nonRF(Sel_Fast)');

Prop_TgtOn_Good_nonRF = sum(hPSTH_TgtOn_Good_nonRF_sel,1)/size(hPSTH_TgtOn_Good_nonRF_sel,1);

hPSTH_SacOn_Good_nonRF_sel = cell2mat(hPSTH_SacOn_Good_nonRF(Sel_Fast)');%(Sel_All);

Prop_SacOn_Good_nonRF = sum(hPSTH_SacOn_Good_nonRF_sel,1)/size(hPSTH_SacOn_Good_nonRF_sel,1);

hPSTH_TgtOn_Bad_nonRF_sel = cell2mat(hPSTH_TgtOn_Bad_nonRF(Sel_Fast)');

Prop_TgtOn_Bad_nonRF = sum(hPSTH_TgtOn_Bad_nonRF_sel,1)/size(hPSTH_TgtOn_Bad_nonRF_sel,1);

hPSTH_SacOn_Bad_nonRF_sel = cell2mat(hPSTH_SacOn_Bad_nonRF(Sel_Fast)');%(Sel_All);

Prop_SacOn_Bad_nonRF = sum(hPSTH_SacOn_Bad_nonRF_sel,1)/size(hPSTH_SacOn_Bad_nonRF_sel,1);

%OneDR

%% Fast group
hPSTH_TgtOn_Value_RF_sel_oneDR_Fast = cell2mat(hPSTH_TgtOn_Value_RF_oneDR(Sel_Fast)');%(Sel_All);

Prop_TgtOn_Value_RF_oneDR_Fast = sum(hPSTH_TgtOn_Value_RF_sel_oneDR_Fast,1)/size(hPSTH_TgtOn_Value_RF_sel_oneDR_Fast,1);

hPSTH_SacOn_Value_RF_sel_oneDR_Fast = cell2mat(hPSTH_SacOn_Value_RF_oneDR(Sel_Fast)');%(Sel_All);

Prop_SacOn_Value_RF_oneDR_Fast = sum(hPSTH_SacOn_Value_RF_sel_oneDR_Fast,1)/size(hPSTH_SacOn_Value_RF_sel_oneDR_Fast,1);



hPSTH_TgtOn_Value_nonRF_sel_oneDR_Fast = cell2mat(hPSTH_TgtOn_Value_nonRF_oneDR(Sel_Fast)');%(Sel_All);

Prop_TgtOn_Value_nonRF_oneDR_Fast = sum(hPSTH_TgtOn_Value_nonRF_sel_oneDR_Fast,1)/size(hPSTH_TgtOn_Value_nonRF_sel_oneDR_Fast,1);

hPSTH_SacOn_Value_nonRF_sel_oneDR_Fast = cell2mat(hPSTH_SacOn_Value_nonRF_oneDR(Sel_Fast)');%(Sel_All);

Prop_SacOn_Value_nonRF_oneDR_Fast = sum(hPSTH_SacOn_Value_nonRF_sel_oneDR_Fast,1)/size(hPSTH_SacOn_Value_nonRF_sel_oneDR_Fast,1);

%% Others group
hPSTH_TgtOn_Value_RF_sel_oneDR_Others = cell2mat(hPSTH_TgtOn_Value_RF_oneDR(Sel_Others)');%(Sel_All);

Prop_TgtOn_Value_RF_oneDR_Others = sum(hPSTH_TgtOn_Value_RF_sel_oneDR_Others,1)/size(hPSTH_TgtOn_Value_RF_sel_oneDR_Others,1);
hPSTH_SacOn_Value_RF_sel_oneDR_Others = cell2mat(hPSTH_SacOn_Value_RF_oneDR(Sel_Others)');%(Sel_All);

Prop_SacOn_Value_RF_oneDR_Others = sum(hPSTH_SacOn_Value_RF_sel_oneDR_Others,1)/size(hPSTH_SacOn_Value_RF_sel_oneDR_Others,1);



hPSTH_TgtOn_Value_nonRF_sel_oneDR_Other = cell2mat(hPSTH_TgtOn_Value_nonRF_oneDR(Sel_Others)');%(Sel_All);

Prop_TgtOn_Value_nonRF_oneDR_Other = sum(hPSTH_TgtOn_Value_nonRF_sel_oneDR_Other,1)/size(hPSTH_TgtOn_Value_nonRF_sel_oneDR_Other,1);

hPSTH_SacOn_Value_nonRF_sel_oneDR_Other = cell2mat(hPSTH_SacOn_Value_nonRF_oneDR(Sel_Others)');%(Sel_All);

Prop_SacOn_Value_nonRF_oneDR_Other = sum(hPSTH_SacOn_Value_nonRF_sel_oneDR_Other,1)/size(hPSTH_SacOn_Value_nonRF_sel_oneDR_Other,1);

% Statistics 




% All group
hPSTH_TgtOn_Value_RF_sel_oneDR_All = cell2mat(hPSTH_TgtOn_Value_RF_oneDR(Sel_All)');%(Sel_All);

Prop_TgtOn_Value_RF_oneDR_All = sum(hPSTH_TgtOn_Value_RF_sel_oneDR_All,1)/size(hPSTH_TgtOn_Value_RF_sel_oneDR_All,1);
hPSTH_SacOn_Value_RF_sel_oneDR_All = cell2mat(hPSTH_SacOn_Value_RF_oneDR(Sel_All)');%(Sel_All);

Prop_SacOn_Value_RF_oneDR_All = sum(hPSTH_SacOn_Value_RF_sel_oneDR_All,1)/size(hPSTH_SacOn_Value_RF_sel_oneDR_All,1);



hPSTH_TgtOn_Value_nonRF_sel_oneDR_All = cell2mat(hPSTH_TgtOn_Value_nonRF_oneDR(Sel_All)');%(Sel_All);

Prop_TgtOn_Value_nonRF_oneDR_All = sum(hPSTH_TgtOn_Value_nonRF_sel_oneDR_All,1)/size(hPSTH_TgtOn_Value_nonRF_sel_oneDR_All,1);

hPSTH_SacOn_Value_nonRF_sel_oneDR_All = cell2mat(hPSTH_SacOn_Value_nonRF_oneDR(Sel_All)');%(Sel_All);

Prop_SacOn_Value_nonRF_oneDR_All = sum(hPSTH_SacOn_Value_nonRF_sel_oneDR_All,1)/size(hPSTH_SacOn_Value_nonRF_sel_oneDR_All,1);


%% 


%{
%Good vs Bad
hPSTH_TgtOn_Good_RF_sel_oneDR = cell2mat(hPSTH_TgtOn_Good_RF_oneDR(Sel_Fast)');

Prop_TgtOn_Good_RF_oneDR = sum(hPSTH_TgtOn_Good_RF_sel_oneDR,1)/size(hPSTH_TgtOn_Good_RF_sel_oneDR,1);

hPSTH_SacOn_Good_RF_sel_oneDR = cell2mat(hPSTH_SacOn_Good_RF_oneDR(Sel_Fast)');%(Sel_All);

Prop_SacOn_Good_RF_oneDR = sum(hPSTH_SacOn_Good_RF_sel_oneDR,1)/size(hPSTH_SacOn_Good_RF_sel_oneDR,1);

hPSTH_TgtOn_Bad_RF_sel_oneDR = cell2mat(hPSTH_TgtOn_Bad_RF_oneDR(Sel_Fast)');

Prop_TgtOn_Bad_RF_oneDR = sum(hPSTH_TgtOn_Bad_RF_sel_oneDR,1)/size(hPSTH_TgtOn_Bad_RF_sel_oneDR,1);

hPSTH_SacOn_Bad_RF_sel_oneDR = cell2mat(hPSTH_SacOn_Bad_RF_oneDR(Sel_Fast)');%(Sel_All);

Prop_SacOn_Bad_RF_oneDR = sum(hPSTH_SacOn_Bad_RF_sel_oneDR,1)/size(hPSTH_SacOn_Bad_RF_sel_oneDR,1);


hPSTH_TgtOn_Good_nonRF_sel_oneDR = cell2mat(hPSTH_TgtOn_Good_nonRF_oneDR(Sel_Fast)');

Prop_TgtOn_Good_nonRF_oneDR = sum(hPSTH_TgtOn_Good_nonRF_sel_oneDR,1)/size(hPSTH_TgtOn_Good_nonRF_sel_oneDR,1);

hPSTH_SacOn_Good_nonRF_sel_oneDR = cell2mat(hPSTH_SacOn_Good_nonRF_oneDR(Sel_Fast)');%(Sel_All);

Prop_SacOn_Good_nonRF_oneDR = sum(hPSTH_SacOn_Good_nonRF_sel_oneDR,1)/size(hPSTH_SacOn_Good_nonRF_sel_oneDR,1);

hPSTH_TgtOn_Bad_nonRF_sel_oneDR = cell2mat(hPSTH_TgtOn_Bad_nonRF_oneDR(Sel_Fast)');

Prop_TgtOn_Bad_nonRF_oneDR = sum(hPSTH_TgtOn_Bad_nonRF_sel_oneDR,1)/size(hPSTH_TgtOn_Bad_nonRF_sel_oneDR,1);

hPSTH_SacOn_Bad_nonRF_sel_oneDR = cell2mat(hPSTH_SacOn_Bad_nonRF_oneDR(Sel_Fast)');%(Sel_All);

Prop_SacOn_Bad_nonRF_oneDR = sum(hPSTH_SacOn_Bad_nonRF_sel_oneDR,1)/size(hPSTH_SacOn_Bad_nonRF_sel_oneDR,1);
%}





%{

%%For Fast
%Stable Value Only
hPSTH_TgtOn_Value_RF_sel_Fast_Or = hPSTH_TgtOn_Value_RF_sel_Fast | ~hPSTH_TgtOn_Value_RF_sel_oneDR_Fast;
hPSTH_TgtOn_Value_RF_sel_Fast_stableOnly = hPSTH_TgtOn_Value_RF_sel_Fast & ~hPSTH_TgtOn_Value_RF_sel_oneDR_Fast;
Prop_TgtOn_Value_RF_Fast_stableOnly = sum(hPSTH_TgtOn_Value_RF_sel_Fast_stableOnly,1)/size(hPSTH_TgtOn_Value_RF_sel_Fast_stableOnly,1);
%Prop_TgtOn_Value_RF_Fast_stableOnly = sum(hPSTH_TgtOn_Value_RF_sel_Fast_stableOnly,1);%./sum(hPSTH_TgtOn_Value_RF_sel_Fast_Or,1);

%OneDR Only

hPSTH_TgtOn_Value_RF_sel_Fast_OneDROnly = ~hPSTH_TgtOn_Value_RF_sel_Fast & hPSTH_TgtOn_Value_RF_sel_oneDR_Fast;
Prop_TgtOn_Value_RF_Fast_OneDROnly = sum(hPSTH_TgtOn_Value_RF_sel_Fast_OneDROnly,1)/size(hPSTH_TgtOn_Value_RF_sel_Fast_OneDROnly,1);
%Prop_TgtOn_Value_RF_Fast_OneDROnly = sum(hPSTH_TgtOn_Value_RF_sel_Fast_OneDROnly,1);%./sum(hPSTH_TgtOn_Value_RF_sel_Fast_Or,1);

%Stable Value and OneDR Both
hPSTH_TgtOn_Value_RF_sel_Fast_Both = hPSTH_TgtOn_Value_RF_sel_Fast & hPSTH_TgtOn_Value_RF_sel_oneDR_Fast;
Prop_TgtOn_Value_RF_Fast_Both = sum(hPSTH_TgtOn_Value_RF_sel_Fast_Both,1)/size(hPSTH_TgtOn_Value_RF_sel_Fast_Both,1);
%Prop_TgtOn_Value_RF_Fast_Both = sum(hPSTH_TgtOn_Value_RF_sel_Fast_Both,1);%./sum(hPSTH_TgtOn_Value_RF_sel_Fast_Or,1);

%Stable Value Only
hPSTH_SacOn_Value_RF_sel_Fast_Or = hPSTH_SacOn_Value_RF_sel_Fast | ~hPSTH_SacOn_Value_RF_sel_oneDR_Fast;
hPSTH_SacOn_Value_RF_sel_Fast_stableOnly = hPSTH_SacOn_Value_RF_sel_Fast & ~hPSTH_SacOn_Value_RF_sel_oneDR_Fast;
Prop_SacOn_Value_RF_Fast_stableOnly = sum(hPSTH_SacOn_Value_RF_sel_Fast_stableOnly,1)/size(hPSTH_SacOn_Value_RF_sel_Fast_stableOnly,1);
%Prop_TgtOn_Value_RF_Fast_stableOnly = sum(hPSTH_TgtOn_Value_RF_sel_Fast_stableOnly,1);%./sum(hPSTH_TgtOn_Value_RF_sel_Fast_Or,1);

%OneDR Only

hPSTH_SacOn_Value_RF_sel_Fast_OneDROnly = ~hPSTH_SacOn_Value_RF_sel_Fast & hPSTH_SacOn_Value_RF_sel_oneDR_Fast;
Prop_SacOn_Value_RF_Fast_OneDROnly = sum(hPSTH_SacOn_Value_RF_sel_Fast_OneDROnly,1)/size(hPSTH_SacOn_Value_RF_sel_Fast_OneDROnly,1);
%Prop_TgtOn_Value_RF_Fast_OneDROnly = sum(hPSTH_TgtOn_Value_RF_sel_Fast_OneDROnly,1);%./sum(hPSTH_TgtOn_Value_RF_sel_Fast_Or,1);

%Stable Value and OneDR Both
hPSTH_SacOn_Value_RF_sel_Fast_Both = hPSTH_SacOn_Value_RF_sel_Fast & hPSTH_SacOn_Value_RF_sel_oneDR_Fast;
Prop_SacOn_Value_RF_Fast_Both = sum(hPSTH_SacOn_Value_RF_sel_Fast_Both,1)/size(hPSTH_SacOn_Value_RF_sel_Fast_Both,1);
%Prop_TgtOn_Value_RF_Fast_Both = sum(hPSTH_TgtOn_Value_RF_sel_Fast_Both,1);%./sum(hPSTH_TgtOn_Value_RF_sel_Fast_Or,1);




% For others
%Stable Value Only
hPSTH_TgtOn_Value_RF_sel_Others_Or = hPSTH_TgtOn_Value_RF_sel_Others | ~hPSTH_TgtOn_Value_RF_sel_oneDR_Others;
hPSTH_TgtOn_Value_RF_sel_Others_stableOnly = hPSTH_TgtOn_Value_RF_sel_Others & ~hPSTH_TgtOn_Value_RF_sel_oneDR_Others;
Prop_TgtOn_Value_RF_Others_stableOnly = sum(hPSTH_TgtOn_Value_RF_sel_Others_stableOnly,1)/size(hPSTH_TgtOn_Value_RF_sel_Others_stableOnly,1);
%Prop_TgtOn_Value_RF_Fast_stableOnly = sum(hPSTH_TgtOn_Value_RF_sel_Fast_stableOnly,1);%./sum(hPSTH_TgtOn_Value_RF_sel_Fast_Or,1);

%OneDR Only

hPSTH_TgtOn_Value_RF_sel_Others_OneDROnly = ~hPSTH_TgtOn_Value_RF_sel_Others & hPSTH_TgtOn_Value_RF_sel_oneDR_Others;
Prop_TgtOn_Value_RF_Others_OneDROnly = sum(hPSTH_TgtOn_Value_RF_sel_Others_OneDROnly,1)/size(hPSTH_TgtOn_Value_RF_sel_Others_OneDROnly,1);
%Prop_TgtOn_Value_RF_Fast_OneDROnly = sum(hPSTH_TgtOn_Value_RF_sel_Fast_OneDROnly,1);%./sum(hPSTH_TgtOn_Value_RF_sel_Fast_Or,1);

%Stable Value and OneDR Both
hPSTH_TgtOn_Value_RF_sel_Others_Both = hPSTH_TgtOn_Value_RF_sel_Others & hPSTH_TgtOn_Value_RF_sel_oneDR_Others;
Prop_TgtOn_Value_RF_Others_Both = sum(hPSTH_TgtOn_Value_RF_sel_Others_Both,1)/size(hPSTH_TgtOn_Value_RF_sel_Others_Both,1);
%Prop_TgtOn_Value_RF_Fast_Both = sum(hPSTH_TgtOn_Value_RF_sel_Fast_Both,1);%./sum(hPSTH_TgtOn_Value_RF_sel_Fast_Or,1);

%% Saccade on
%Stable Value Only
hPSTH_SacOn_Value_RF_sel_Others_Or = hPSTH_SacOn_Value_RF_sel_Others | ~hPSTH_SacOn_Value_RF_sel_oneDR_Others;
hPSTH_SacOn_Value_RF_sel_Others_stableOnly = hPSTH_SacOn_Value_RF_sel_Others & ~hPSTH_SacOn_Value_RF_sel_oneDR_Others;
Prop_SacOn_Value_RF_Others_stableOnly = sum(hPSTH_SacOn_Value_RF_sel_Others_stableOnly,1)/size(hPSTH_SacOn_Value_RF_sel_Others_stableOnly,1);
%Prop_TgtOn_Value_RF_Fast_stableOnly = sum(hPSTH_TgtOn_Value_RF_sel_Fast_stableOnly,1);%./sum(hPSTH_TgtOn_Value_RF_sel_Fast_Or,1);

%OneDR Only

hPSTH_SacOn_Value_RF_sel_Others_OneDROnly = ~hPSTH_SacOn_Value_RF_sel_Others & hPSTH_SacOn_Value_RF_sel_oneDR_Others;
Prop_SacOn_Value_RF_Others_OneDROnly = sum(hPSTH_SacOn_Value_RF_sel_Others_OneDROnly,1)/size(hPSTH_SacOn_Value_RF_sel_Others_OneDROnly,1);
%Prop_TgtOn_Value_RF_Fast_OneDROnly = sum(hPSTH_TgtOn_Value_RF_sel_Fast_OneDROnly,1);%./sum(hPSTH_TgtOn_Value_RF_sel_Fast_Or,1);

%Stable Value and OneDR Both
hPSTH_SacOn_Value_RF_sel_Others_Both = hPSTH_SacOn_Value_RF_sel_Others & hPSTH_SacOn_Value_RF_sel_oneDR_Others;
Prop_SacOn_Value_RF_Others_Both = sum(hPSTH_SacOn_Value_RF_sel_Others_Both,1)/size(hPSTH_SacOn_Value_RF_sel_Others_Both,1);
%Prop_TgtOn_Value_RF_Fast_Both = sum(hPSTH_TgtOn_Value_RF_sel_Fast_Both,1);%./sum(hPSTH_TgtOn_Value_RF_sel_Fast_Or,1);


%Select out value selective cells:
%Value_Stable_RF_sel_All = Value_Stable_RF==1 & Sel_All;

Value_Stable_RF_ROC_All=Value_Location_RF_ROC(Sel_All);
Value_Location_RF_ROC_All=Value_Stable_RF_ROC(Sel_All);

Value_Stable_RF_ROC_All_Sig = Value_Stable_RF_ROC(Value_Stable_RF'==1&Sel_All);
Value_Location_RF_ROC_All_StableSig = Value_Location_RF_ROC(Value_Stable_RF'==1&Sel_All);

Value_Location_RF_ROC_All_Sig = Value_Location_RF_ROC(Value_Location_RF'==1&Sel_All);
Value_Stable_RF_ROC_All_LocationSig = Value_Stable_RF_ROC(Value_Location_RF'==1&Sel_All);

Value_Location_RF_ROC_All_BothSig = Value_Location_RF_ROC(Value_Location_RF'==1&Value_Stable_RF'==1&Sel_All);
Value_Stable_RF_ROC_All_BothSig = Value_Stable_RF_ROC(Value_Location_RF'==1&Value_Stable_RF'==1&Sel_All);


%{
figure
plot(Value_Stable_RF_ROC_All,Value_Location_RF_ROC_All,'ok');
hold on
plot(Value_Stable_RF_ROC_All_Sig,Value_Location_RF_ROC_All_StableSig,'or','MarkerFaceColor','r');
plot(Value_Stable_RF_ROC_All_LocationSig,Value_Location_RF_ROC_All_Sig,'ob','MarkerFaceColor','b');
plot(Value_Stable_RF_ROC_All_BothSig,Value_Location_RF_ROC_All_BothSig,'og','MarkerFaceColor','g');

xlabel('Object Value');
ylabel('Location Value');
%}
%}

%Statistics 
NumberValue_Fast_Stable = sum(Value_Stable_RF_Prop_Sig(Sel_Fast));
TotalNumberValue_Fast_Stable = length(Value_Stable_RF_Prop_Sig(Sel_Fast));

NumberValue_Other_Stable = sum(Value_Stable_RF_Prop_Sig(Sel_Others));
TotalNumberValue_Other_Stable = length(Value_Stable_RF_Prop_Sig(Sel_Others));

 
Prop_Value_Fast_Stable = NumberValue_Fast_Stable/TotalNumberValue_Fast_Stable;
Prop_Value_Other_Stable = NumberValue_Other_Stable/TotalNumberValue_Other_Stable;

[chi_squared, df, p_value_stableprop, is_significant] = chisq_prop_test(NumberValue_Fast_Stable,TotalNumberValue_Fast_Stable-NumberValue_Fast_Stable,NumberValue_Other_Stable,TotalNumberValue_Other_Stable-NumberValue_Other_Stable);

disp('Significance for proportion difference between fast and others:')
p_value_stableprop


NumberValue_Fast_OneDR = sum(Value_OneDR_RF_Prop_Sig(Sel_Fast));
TotalNumberValue_Fast_OneDR = length(Value_OneDR_RF_Prop_Sig(Sel_Fast));

NumberValue_Other_OneDR = sum(Value_OneDR_RF_Prop_Sig(Sel_Others));
TotalNumberValue_Other_OneDR = length(Value_OneDR_RF_Prop_Sig(Sel_Others));

Prop_Value_Fast_OneDR = NumberValue_Fast_OneDR/TotalNumberValue_Fast_OneDR;
Prop_Value_Other_OneDR = NumberValue_Other_OneDR/TotalNumberValue_Other_OneDR;

[chi_squared, df, p_value_oneDRprop, is_significant] = chisq_prop_test(NumberValue_Fast_OneDR,TotalNumberValue_Fast_OneDR-NumberValue_Fast_OneDR,NumberValue_Other_OneDR,TotalNumberValue_Other_OneDR-NumberValue_Other_OneDR);

disp('Significance for proportion difference between fast and others:')
p_value_oneDRprop


NumberValue_Fast_PV = sum(Value_PV_RF_Prop_Sig(Sel_Fast&PV_Index'==1));
TotalNumberValue_Fast_PV = length(Value_PV_RF_Prop_Sig(Sel_Fast & PV_Index'==1));
NumberValue_Other_PV = sum(Value_PV_RF_Prop_Sig(Sel_Others & PV_Index'==1));
TotalNumberValue_Other_PV = length(Value_PV_RF_Prop_Sig(Sel_Others & PV_Index'==1));

Prop_Value_Fast_PV = NumberValue_Fast_PV/TotalNumberValue_Fast_PV;
Prop_Value_Other_PV = NumberValue_Other_PV/TotalNumberValue_Other_PV;

[chi_squared, df, p_value_PVprop, is_significant] = chisq_prop_test(NumberValue_Fast_PV,TotalNumberValue_Fast_PV-NumberValue_Fast_PV,NumberValue_Other_PV,TotalNumberValue_Other_PV-NumberValue_Other_PV);

disp('Significance for proportion difference between fast and others:')
p_value_PVprop



Prop_Sig_Location_RF = sum(Value_OneDR_RF_Prop_Sig(Sel_Fast))/length(Value_OneDR_RF_Prop_Sig(Sel_Fast));
Prop_Sig_Stable_RF = sum(Value_Stable_RF_Prop_Sig(Sel_Fast))/length(Value_Stable_RF_Prop_Sig(Sel_Fast));

Value_OneDR_RF_Prop_SigOnly = Value_OneDR_RF_Prop_Sig & ~Value_Stable_RF_Prop_Sig;
Num_Sig_LocationOnly_RF = sum(Value_OneDR_RF_Prop_SigOnly(Sel_Fast));

Value_Stable_RF_Prop_SigOnly = Value_Stable_RF_Prop_Sig & ~Value_OneDR_RF_Prop_Sig;
Num_Sig_StableOnly_RF = sum(Value_Stable_RF_Prop_SigOnly(Sel_Fast));

Value_Both_RF_Prop_SigOnly = Value_Stable_RF_Prop_Sig & Value_OneDR_RF_Prop_Sig;
Num_Sig_Both_RF = sum(Value_Both_RF_Prop_SigOnly(Sel_Fast));

Prop_Sig_LocationOnly_RF = Num_Sig_LocationOnly_RF/(Num_Sig_LocationOnly_RF+Num_Sig_StableOnly_RF+Num_Sig_Both_RF);
Prop_Sig_StableOnly_RF = Num_Sig_StableOnly_RF/(Num_Sig_LocationOnly_RF+Num_Sig_StableOnly_RF+Num_Sig_Both_RF);
Prop_Sig_Both_RF = Num_Sig_Both_RF/(Num_Sig_LocationOnly_RF+Num_Sig_StableOnly_RF+Num_Sig_Both_RF);

%For others
Prop_Sig_Location_RF_Others = sum(Value_OneDR_RF_Prop_Sig(Sel_Others))/length(Value_OneDR_RF_Prop_Sig(Sel_Others));
Prop_Sig_Stable_RF_Others = sum(Value_Stable_RF_Prop_Sig(Sel_Others))/length(Value_Stable_RF_Prop_Sig(Sel_Others));


Num_Sig_LocationOnly_RF_Others = sum(Value_OneDR_RF_Prop_SigOnly(Sel_Others));
Num_Sig_StableOnly_RF_Others = sum(Value_Stable_RF_Prop_SigOnly(Sel_Others));
Num_Sig_Both_RF_Others = sum(Value_Both_RF_Prop_SigOnly(Sel_Others));

Prop_Sig_LocationOnly_RF_Others = Num_Sig_LocationOnly_RF_Others/(Num_Sig_LocationOnly_RF_Others+Num_Sig_StableOnly_RF_Others+Num_Sig_Both_RF_Others);
Prop_Sig_StableOnly_RF_Others = Num_Sig_StableOnly_RF_Others/(Num_Sig_LocationOnly_RF_Others+Num_Sig_StableOnly_RF_Others+Num_Sig_Both_RF_Others);
Prop_Sig_Both_RF_Others = Num_Sig_Both_RF_Others/(Num_Sig_LocationOnly_RF_Others+Num_Sig_StableOnly_RF_Others+Num_Sig_Both_RF_Others);


%Population PSTH
PSTH_TgtOn_Good_RF_Stable_Fast_Value = PSTH_TgtOn_Good_RF_Stable(Value_Stable_RF_Prop_Sig' & Sel_Fast,:);
PSTH_TgtOn_Bad_RF_Stable_Fast_Value = PSTH_TgtOn_Bad_RF_Stable(Value_Stable_RF_Prop_Sig' & Sel_Fast,:);

%PSTH_TgtOn_Diff_RF_Stable_Fast_Value=PSTH_TgtOn_Good_RF_Stable_Fast_Value-PSTH_TgtOn_Bad_RF_Stable_Fast_Value;



Value_Stable_RF_Diff_Fast_Value=Value_Stable_RF_Diff(Value_Stable_RF_Prop_Sig' & Sel_Fast);
GoodPrefer_Stable = Value_Stable_RF_Diff_Fast_Value>0;
BadPrefer_Stable = ~GoodPrefer_Stable;


PSTH_TgtOn_Good_RF_Stable_Fast_Value_GoodPrefer = mean(PSTH_TgtOn_Good_RF_Stable_Fast_Value(GoodPrefer_Stable,:));
PSTH_TgtOn_Bad_RF_Stable_Fast_Value_GoodPrefer = mean(PSTH_TgtOn_Bad_RF_Stable_Fast_Value(GoodPrefer_Stable,:));

PSTH_TgtOn_Good_RF_Stable_Fast_Value_GoodPrefer_Sem = std(PSTH_TgtOn_Good_RF_Stable_Fast_Value(GoodPrefer_Stable,:))/(4*sqrt(sum(GoodPrefer_Stable)));
PSTH_TgtOn_Bad_RF_Stable_Fast_Value_GoodPrefer_Sem = std(PSTH_TgtOn_Bad_RF_Stable_Fast_Value(GoodPrefer_Stable,:))/(4*sqrt(sum(GoodPrefer_Stable)));

PSTH_TgtOn_Good_RF_Stable_Fast_Value_BadPrefer = mean(PSTH_TgtOn_Good_RF_Stable_Fast_Value(BadPrefer_Stable,:));
PSTH_TgtOn_Bad_RF_Stable_Fast_Value_BadPrefer = mean(PSTH_TgtOn_Bad_RF_Stable_Fast_Value(BadPrefer_Stable,:));

PSTH_TgtOn_Good_RF_Stable_Fast_Value_BadPrefer_Sem = std(PSTH_TgtOn_Good_RF_Stable_Fast_Value(BadPrefer_Stable,:))/(4*sqrt(sum(BadPrefer_Stable)));
PSTH_TgtOn_Bad_RF_Stable_Fast_Value_BadPrefer_Sem = std(PSTH_TgtOn_Bad_RF_Stable_Fast_Value(BadPrefer_Stable,:))/(4*sqrt(sum(BadPrefer_Stable)));

%% Stable AUC
AUCPSTH_TgtOn_Value_RF_stable = cell2mat(AUCPSTH_TgtOn_Value_RF_stable');
AUCPSTH_TgtOn_Value_RF_stable_Fast=AUCPSTH_TgtOn_Value_RF_stable(Value_Stable_RF_Prop_Sig' & Sel_Fast,:);
AUCPSTH_TgtOn_Value_RF_stable_Other=AUCPSTH_TgtOn_Value_RF_stable(Value_Stable_RF_Prop_Sig' & Sel_Others,:);

% Fast
AUCPSTH_TgtOn_RF_Stable_Fast_Value_GoodPrefer = mean(AUCPSTH_TgtOn_Value_RF_stable_Fast(GoodPrefer_Stable,:));

AUCPSTH_TgtOn_RF_Stable_Fast_Value_GoodPrefer_Sem = std(AUCPSTH_TgtOn_Value_RF_stable_Fast(GoodPrefer_Stable,:))/(sqrt(sum(GoodPrefer_Stable)));

AUCPSTH_TgtOn_RF_Stable_Fast_Value_BadPrefer = mean(AUCPSTH_TgtOn_Value_RF_stable_Fast(BadPrefer_Stable,:));

AUCPSTH_TgtOn_RF_Stable_Fast_Value_BadPrefer_Sem = std(AUCPSTH_TgtOn_Value_RF_stable_Fast(BadPrefer_Stable,:))/(sqrt(sum(BadPrefer_Stable)));


Value_Stable_RF_Diff_Others_Value=Value_Stable_RF_Diff(Value_Stable_RF_Prop_Sig' & Sel_Others);
GoodPrefer_Stable_Others = Value_Stable_RF_Diff_Others_Value>0;
BadPrefer_Stable_Others = ~GoodPrefer_Stable_Others;

% Others
AUCPSTH_TgtOn_RF_Stable_Others_Value_GoodPrefer = mean(AUCPSTH_TgtOn_Value_RF_stable_Other(GoodPrefer_Stable_Others,:));

AUCPSTH_TgtOn_RF_Stable_Others_Value_GoodPrefer_Sem = std(AUCPSTH_TgtOn_Value_RF_stable_Other(GoodPrefer_Stable_Others,:))/(sqrt(sum(GoodPrefer_Stable_Others)));

AUCPSTH_TgtOn_RF_Stable_Others_Value_BadPrefer = mean(AUCPSTH_TgtOn_Value_RF_stable_Other(BadPrefer_Stable_Others,:));

AUCPSTH_TgtOn_RF_Stable_Others_Value_BadPrefer_Sem = std(AUCPSTH_TgtOn_Value_RF_stable_Other(BadPrefer_Stable_Others,:))/(sqrt(sum(BadPrefer_Stable_Others)));

%% Passive Viewing
hPSTH_TgtOn_Value_RF_PV=cell2mat(hPSTH_TgtOn_Value_RF_PV');

hPSTH_TgtOn_Value_RF_PV_Fast=hPSTH_TgtOn_Value_RF_PV(Sel_Fast&PV_Index'==1,:);

Prop_Sig_PV_Fast=sum(hPSTH_TgtOn_Value_RF_PV_Fast)/size(hPSTH_TgtOn_Value_RF_PV_Fast,1);



hPSTH_TgtOn_Value_RF_PV_Others=hPSTH_TgtOn_Value_RF_PV(Sel_Others&PV_Index'==1,:);

Prop_Sig_PV_Others=sum(hPSTH_TgtOn_Value_RF_PV_Others)/size(hPSTH_TgtOn_Value_RF_PV_Others,1);

PSTH_TgtOn_Good_RF_PV=(cell2mat(PSTH_TgtOn_Good_RF_PV')')';
PSTH_TgtOn_Bad_RF_PV=(cell2mat(PSTH_TgtOn_Bad_RF_PV')')';

Sel_Sig_PV_Fast = Value_PV_RF_Prop_Sig'==1 &PV_Index'==1&Sel_Fast==1;
Sel_Sig_PV_Others = Value_PV_RF_Prop_Sig'==1 &PV_Index'==1&Sel_Others==1;
Good_sel_PV = Value_PV_RF_Diff>0;
Bad_sel_PV = ~Good_sel_PV;

PSTH_TgtOn_Good_RF_PV_Fast_Sig_Good = PSTH_TgtOn_Good_RF_PV(Good_sel_PV'&Sel_Sig_PV_Fast,:);
PSTH_TgtOn_Bad_RF_PV_Fast_Sig_Good = PSTH_TgtOn_Bad_RF_PV(Good_sel_PV'&Sel_Sig_PV_Fast,:);

PSTH_TgtOn_Good_RF_PV_Fast_Sig_Bad = PSTH_TgtOn_Good_RF_PV(Bad_sel_PV'&Sel_Sig_PV_Fast,:);
PSTH_TgtOn_Bad_RF_PV_Fast_Sig_Bad = PSTH_TgtOn_Bad_RF_PV(Bad_sel_PV'&Sel_Sig_PV_Fast,:);



PSTH_TgtOn_Good_RF_PV_Fast_Sig_Good_Mean=mean(PSTH_TgtOn_Good_RF_PV_Fast_Sig_Good,1);
PSTH_TgtOn_Bad_RF_PV_Fast_Sig_Good_Mean=mean(PSTH_TgtOn_Bad_RF_PV_Fast_Sig_Good,1);

PSTH_TgtOn_Good_RF_PV_Fast_Sig_Good_Sem=std(PSTH_TgtOn_Good_RF_PV_Fast_Sig_Good,[],1)/(4*sqrt(size(PSTH_TgtOn_Good_RF_PV_Fast_Sig_Good,1)));
PSTH_TgtOn_Bad_RF_PV_Fast_Sig_Good_Sem=std(PSTH_TgtOn_Bad_RF_PV_Fast_Sig_Good,[],1)/(4*sqrt(size(PSTH_TgtOn_Good_RF_PV_Fast_Sig_Good,1)));

PSTH_TgtOn_Good_RF_PV_Fast_Sig_Bad_Mean=mean(PSTH_TgtOn_Good_RF_PV_Fast_Sig_Bad,1);
PSTH_TgtOn_Bad_RF_PV_Fast_Sig_Bad_Mean=mean(PSTH_TgtOn_Bad_RF_PV_Fast_Sig_Bad,1);

PSTH_TgtOn_Good_RF_PV_Fast_Sig_Bad_Sem=std(PSTH_TgtOn_Good_RF_PV_Fast_Sig_Bad,[],1)/(4*sqrt(size(PSTH_TgtOn_Good_RF_PV_Fast_Sig_Bad,1)));
PSTH_TgtOn_Bad_RF_PV_Fast_Sig_Bad_Sem=std(PSTH_TgtOn_Bad_RF_PV_Fast_Sig_Bad,[],1)/(4*sqrt(size(PSTH_TgtOn_Good_RF_PV_Fast_Sig_Bad,1)));

%OneDR
Mean_Value_Diff_OneDR=Mean_Value_Diff;
Sel_Sig_OneDR_Fast = Value_OneDR_RF_Prop_Sig'==1 &OneDRIndex'==1 & Sel_Fast==1;
Sel_Sig_OneDR_Others = Value_OneDR_RF_Prop_Sig'==1 &OneDRIndex'==1 & Sel_Others==1;

GoodSel_OneDR = Mean_Value_Diff_OneDR>0;
BadSel_OneDR = ~GoodSel_OneDR;

PSTH_TgtOn_Good_RF_OneDR_Fast_Sig_Good = PSTH_Good_TgtOn_RF_oneDR(GoodSel_OneDR'&Sel_Sig_OneDR_Fast,:);
PSTH_TgtOn_Bad_RF_OneDR_Fast_Sig_Good = PSTH_Bad_TgtOn_RF_oneDR(GoodSel_OneDR'&Sel_Sig_OneDR_Fast,:);

PSTH_TgtOn_Good_RF_OneDR_Fast_Sig_Bad = PSTH_Good_TgtOn_RF_oneDR(BadSel_OneDR'&Sel_Sig_OneDR_Fast,:);
PSTH_TgtOn_Bad_RF_OneDR_Fast_Sig_Bad = PSTH_Bad_TgtOn_RF_oneDR(BadSel_OneDR'&Sel_Sig_OneDR_Fast,:);



PSTH_TgtOn_Good_RF_OneDR_Fast_Sig_Good_Mean=mean(PSTH_TgtOn_Good_RF_OneDR_Fast_Sig_Good,1);
PSTH_TgtOn_Bad_RF_OneDR_Fast_Sig_Good_Mean=mean(PSTH_TgtOn_Bad_RF_OneDR_Fast_Sig_Good,1);

PSTH_TgtOn_Good_RF_OneDR_Fast_Sig_Good_Sem=std(PSTH_TgtOn_Good_RF_OneDR_Fast_Sig_Good,[],1)/(4*sqrt((size(PSTH_TgtOn_Good_RF_OneDR_Fast_Sig_Good,1))));
PSTH_TgtOn_Bad_RF_OneDR_Fast_Sig_Good_Sem=std(PSTH_TgtOn_Bad_RF_OneDR_Fast_Sig_Good,[],1)/(4*sqrt((size(PSTH_TgtOn_Good_RF_OneDR_Fast_Sig_Good,1))));

PSTH_TgtOn_Good_RF_OneDR_Fast_Sig_Bad_Mean=mean(PSTH_TgtOn_Good_RF_OneDR_Fast_Sig_Bad,1);
PSTH_TgtOn_Bad_RF_OneDR_Fast_Sig_Bad_Mean=mean(PSTH_TgtOn_Bad_RF_OneDR_Fast_Sig_Bad,1);

PSTH_TgtOn_Good_RF_OneDR_Fast_Sig_Bad_Sem=std(PSTH_TgtOn_Good_RF_OneDR_Fast_Sig_Bad,[],1)/(4*(sqrt(size(PSTH_TgtOn_Good_RF_OneDR_Fast_Sig_Bad,1))));
PSTH_TgtOn_Bad_RF_OneDR_Fast_Sig_Bad_Sem=std(PSTH_TgtOn_Bad_RF_OneDR_Fast_Sig_Bad,[],1)/(4*sqrt(size(PSTH_TgtOn_Good_RF_OneDR_Fast_Sig_Bad,1)));


AUCPSTH_TgtOn_RF_oneDR=cell2mat(AUCPSTH_TgtOn_Value_RF_oneDR');

AUCPSTH_TgtOn_RF_OneDR_Fast_Sig_Good = AUCPSTH_TgtOn_RF_oneDR(GoodSel_OneDR'&Sel_Sig_OneDR_Fast,:);
AUCPSTH_TgtOn_RF_OneDR_Fast_Sig_Bad = AUCPSTH_TgtOn_RF_oneDR(BadSel_OneDR'&Sel_Sig_OneDR_Fast,:);

AUCPSTH_TgtOn_RF_OneDR_Others_Sig_Good = AUCPSTH_TgtOn_RF_oneDR(GoodSel_OneDR'&Sel_Sig_OneDR_Fast,:);
AUCPSTH_TgtOn_RF_OneDR_Others_Sig_Bad = AUCPSTH_TgtOn_RF_oneDR(BadSel_OneDR'&Sel_Sig_OneDR_Fast,:);

AUCPSTH_TgtOn_RF_OneDR_Fast_Sig_Good_mean = mean(AUCPSTH_TgtOn_RF_oneDR(GoodSel_OneDR'&Sel_Sig_OneDR_Fast,:),1);
AUCPSTH_TgtOn_RF_OneDR_Fast_Sig_Bad_mean = mean(AUCPSTH_TgtOn_RF_oneDR(BadSel_OneDR'&Sel_Sig_OneDR_Fast,:),1);

AUCPSTH_TgtOn_RF_OneDR_Fast_Sig_Good_sem = std(AUCPSTH_TgtOn_RF_oneDR(GoodSel_OneDR'&Sel_Sig_OneDR_Fast,:),1)/sqrt(sum(GoodSel_OneDR'&Sel_Sig_OneDR_Fast));
AUCPSTH_TgtOn_RF_OneDR_Fast_Sig_Bad_sem = std(AUCPSTH_TgtOn_RF_oneDR(BadSel_OneDR'&Sel_Sig_OneDR_Fast,:),1)/sqrt(sum(GoodSel_OneDR'&Sel_Sig_OneDR_Fast));


AUCPSTH_TgtOn_RF_OneDR_Others_Sig_Good_mean = mean(AUCPSTH_TgtOn_RF_oneDR(GoodSel_OneDR'&Sel_Sig_OneDR_Others,:),1);
AUCPSTH_TgtOn_RF_OneDR_Others_Sig_Bad_mean = mean(AUCPSTH_TgtOn_RF_oneDR(BadSel_OneDR'&Sel_Sig_OneDR_Others,:),1);

AUCPSTH_TgtOn_RF_OneDR_Others_Sig_Good_sem = std(AUCPSTH_TgtOn_RF_oneDR(GoodSel_OneDR'&Sel_Sig_OneDR_Others,:),1)/sqrt(sum(GoodSel_OneDR'&Sel_Sig_OneDR_Others));
AUCPSTH_TgtOn_RF_OneDR_Others_Sig_Bad_sem = std(AUCPSTH_TgtOn_RF_oneDR(BadSel_OneDR'&Sel_Sig_OneDR_Others,:),1)/sqrt(sum(GoodSel_OneDR'&Sel_Sig_OneDR_Others));


% Passive Viewing Task

AUCPSTH_TgtOn_RF_PV=cell2mat(AUCPSTH_TgtOn_Value_RF_PV');


AUCPSTH_TgtOn_RF_PV_Fast_Sig_Good = AUCPSTH_TgtOn_RF_PV(Good_sel_PV'&Sel_Sig_PV_Fast,:);
AUCPSTH_TgtOn_RF_PV_Fast_Sig_Bad = AUCPSTH_TgtOn_RF_PV(Bad_sel_PV'&Sel_Sig_PV_Fast,:);

AUCPSTH_TgtOn_RF_PV_Others_Sig_Good = AUCPSTH_TgtOn_RF_PV(Good_sel_PV'&Sel_Sig_PV_Others,:);
AUCPSTH_TgtOn_RF_PV_Others_Sig_Bad = AUCPSTH_TgtOn_RF_PV(Bad_sel_PV'&Sel_Sig_PV_Others,:);

AUCPSTH_TgtOn_RF_PV_Fast_Sig_Good_mean=mean(AUCPSTH_TgtOn_RF_PV_Fast_Sig_Good,1);
AUCPSTH_TgtOn_RF_PV_Fast_Sig_Good_sem=std(AUCPSTH_TgtOn_RF_PV_Fast_Sig_Good,[],1)/sqrt(sum(Good_sel_PV'&Sel_Sig_PV_Fast));

AUCPSTH_TgtOn_RF_PV_Fast_Sig_Bad_mean=mean(AUCPSTH_TgtOn_RF_PV_Fast_Sig_Bad,1);
AUCPSTH_TgtOn_RF_PV_Fast_Sig_Bad_sem=std(AUCPSTH_TgtOn_RF_PV_Fast_Sig_Bad,[],1)/sqrt(sum(Bad_sel_PV'&Sel_Sig_PV_Fast));

AUCPSTH_TgtOn_RF_PV_Others_Sig_Good_mean=mean(AUCPSTH_TgtOn_RF_PV_Others_Sig_Good,1);
AUCPSTH_TgtOn_RF_PV_Others_Sig_Good_sem=std(AUCPSTH_TgtOn_RF_PV_Others_Sig_Good,[],1)/sqrt(sum(Good_sel_PV'&Sel_Sig_PV_Others));

AUCPSTH_TgtOn_RF_PV_Others_Sig_Bad_mean=mean(AUCPSTH_TgtOn_RF_PV_Others_Sig_Bad,1);
AUCPSTH_TgtOn_RF_PV_Others_Sig_Bad_sem=std(AUCPSTH_TgtOn_RF_PV_Others_Sig_Bad,[],1)/sqrt(sum(Bad_sel_PV'&Sel_Sig_PV_Others));

%{
figure
subplot(1,2,1)
plot(PSTH_TgtOn_Good_RF_OneDR_Fast_Sig_Good_Mean,'-r');
hold on
plot(PSTH_TgtOn_Bad_RF_OneDR_Fast_Sig_Good_Mean,'-b');

subplot(1,2,2)
plot(PSTH_TgtOn_Good_RF_OneDR_Fast_Sig_Bad_Mean,'-r');
hold on
plot(PSTH_TgtOn_Bad_RF_OneDR_Fast_Sig_Bad_Mean,'-b');
%}


%keyboard


%keyboard
%PSTH_TgtOn_Good_RF_Stable_Value = 
%Value_Stable_RF_Diff

%PSTH_Time_TgtOn_stable

%{
figure
subplot(1,2,1);
imagesc(PSTH_Time_TgtOn_stable,1:length(sorted_index(Sel_Fast)),PSTH_TgtOn_Diff_RF_Stable(sorted_index(Sel_Fast),:));
colormap jet;
caxis([-5 5]);
subplot(1,2,2);
imagesc(PSTH_Time_TgtOn_OneDR,1:length(sorted_index(Sel_Fast)),PSTH_TgtOn_Diff_RF_oneDR(sorted_index(Sel_Fast),:));
colormap jet;
caxis([-5 5]);
keyboard
%}

%{
PSTH=[PSTH_TgtOn_Good_RF_Stable(Sel_Fast,:),PSTH_TgtOn_Bad_RF_Stable(Sel_Fast,:)];
PSTH = zscore(PSTH, 0, 2); % Normalize each row
[coeff, score, latent] = pca(PSTH');
%
% Visualize the explained variance
figure;
explained_variance = latent / sum(latent) * 100; % Percentage of variance explained
bar(explained_variance);
xlabel('Principal Component');
ylabel('Explained Variance (%)');
title('Variance Explained by Principal Components');

% Visualize the first two principal components
figure;
plot(PSTH_Time_TgtOn_stable,score(1:160, 1), '-r');
hold on
plot(PSTH_Time_TgtOn_stable,score(161:320, 1), '-b');
xlabel('Time');
ylabel('PC1');
title('PCA: First  Principal Components');

figure;
plot(PSTH_Time_TgtOn_stable,score(1:160,2), '-r');
hold on
plot(PSTH_Time_TgtOn_stable,score(161:320,2), '-b');
xlabel('Time');
ylabel('PC2');
title('PCA: Second Principal Components');

figure;
plot(PSTH_Time_TgtOn_stable,score(1:160,3), '-r');
hold on
plot(PSTH_Time_TgtOn_stable,score(161:320,3), '-b');
xlabel('Time');
ylabel('PC3');
title('PCA: Thrid Principal Components');
figure;
plot(PSTH_Time_TgtOn_stable,score(1:160,4), '-r');
hold on
plot(PSTH_Time_TgtOn_stable,score(161:320,4), '-b');
xlabel('Time');
ylabel('PC4');
title('PCA: Forth Principal Components');
%}

%{
% Visualize the dynamics in the first three PCs
figure;
plot3(score(1:160, 1), score(1:160, 2), score(1:160, 3), '-r');
hold on
plot3(score(161:320, 1), score(161:320, 2), score(161:320, 3), '-b');
xlabel('PC1');
ylabel('PC2');
zlabel('PC3');
title('PCA: First Three Principal Components');
grid on;

%}
%keyboard
%{

PSTH=[PSTH_Good_TgtOn_RF_oneDR(Sel_Fast,:),PSTH_Bad_TgtOn_RF_oneDR(Sel_Fast,:)];
labels = [ones(1, 160), 2 * ones(1, 160)];
%PSTH=[PSTH_Good_TgtOn_RF_oneDR,PSTH_Bad_TgtOn_RF_oneDR];
PSTH = zscore(PSTH, 0, 2); % Normalize each row
[coeff, score, latent] = pca(PSTH');
%}
%{
%{
% Select the first three PCs
num_PCs = 3; % Number of PCs to consider
coeff_first_three = coeff(:, 1:num_PCs); % Coefficients for the first 3 PCs

% Compute neuron contributions to each PC
variance_explained = latent(1:num_PCs)' ./ sum(latent); % Variance explained by each PC
neuron_contributions = sum((coeff_first_three.^2) .* variance_explained, 2);
neuron_contributions = neuron_contributions / sum(neuron_contributions);
%}
% Visualize the explained variance
figure;
explained_variance = latent / sum(latent) * 100; % Percentage of variance explained
cumulative_variance = cumsum(explained_variance);
plot(cumulative_variance,'-ob');
xlabel('Principal Component');
ylabel('Explained Variance (%)');
title('Variance Explained by Principal Components');

% Visualize the first two principal components
figure;
plot(PSTH_Time_TgtOn_stable,score(1:160, 1), '-r');
hold on
plot(PSTH_Time_TgtOn_stable,score(161:320, 1), '-b');
xlabel('Time');
ylabel('PC1');
title('PCA: First  Principal Components');

figure;
plot(PSTH_Time_TgtOn_stable,score(1:160,2), '-r');
hold on
plot(PSTH_Time_TgtOn_stable,score(161:320,2), '-b');
xlabel('Time');
ylabel('PC2');
title('PCA: Second Principal Components');

figure;
plot(PSTH_Time_TgtOn_stable,score(1:160,3), '-r');
hold on
plot(PSTH_Time_TgtOn_stable,score(161:320,3), '-b');
xlabel('Time');
ylabel('PC2');
title('PCA: Thrid Principal Components');

figure;
plot(PSTH_Time_TgtOn_stable,score(1:160,4), '-r');
hold on
plot(PSTH_Time_TgtOn_stable,score(161:320,4), '-b');
xlabel('Time');
ylabel('PC2');
title('PCA: Forth Principal Components');

% Visualize the dynamics in the first three PCs
figure;
plot3(score(1:160, 1), score(1:160, 2), score(1:160, 3), '-r');
hold on
plot3(score(161:320, 1), score(161:320, 2), score(161:320, 3), '-b');

scatter3(score(1, 1), score(1, 2), score(1, 3), 150, 'red', 'filled');
scatter3(score(161, 1), score(161, 2), score(161, 3), 150, 'blue', 'filled');



xlabel('PC1');
ylabel('PC2');
zlabel('PC3');
title('PCA: First Three Principal Components');
grid on;



keyboard
%}



%{
Sel_Vis_Contra = Sel & p_Tgt_RF_stable<0.05;
Sel_Vis_Ispi = Sel & p_Tgt_nonRF_stable<0.05;

Sel_Vis_Opto_Contra = (Sel_All & p_Tgt_RF_stable<0.05)'; % Directly use this population, minor number difference 
Sel_Vis_Opto_Ispi = (Sel_All & p_Tgt_nonRF_stable<0.05)';% Directly use this population, minor number difference 
%}
%Visual & Saccade responses significance
%{
 p_Vis_RF
 VisDiff_RF
 p_Vis_nonRF
 VisDiff_nonRF

 p_pureSac_RF
 pureSacDiff_RF
 p_pureSac_nonRF
 pureSacDiff_nonRF
%}



%{
disp('Number of data with contralateral visual responses for only value analysis:');
disp(sum(Sel_Vis_Contra));
disp('Number of data with contralateral visual responses for value and opto analysis:');
disp(sum(Sel_Vis_Opto_Contra));

disp('Number of data with ipsilateral visual responses for only value analysis:');
disp(sum(Sel_Vis_Ispi));
disp('Number of data with ipsilateral visual responses for value and opto analysis:');
disp(sum(Sel_Vis_Opto_Ispi));

Sel_Vis_Contral_Fast = Sel_Vis_Opto_Contra & ModuType==1;
Sel_Vis_Ipsi_Fast = Sel_Vis_Opto_Ispi & ModuType==1;

disp('Number of data fast responsive contral:');
disp(sum(Sel_Vis_Contral_Fast));
disp('Number of data fast responsive ipsi:');
disp(sum(Sel_Vis_Ipsi_Fast));
%}





%{
 hPSTH_TgtOn_Value_RF(index,:)=Data('hPSTH_TgtOn_Value_RF');
 hPSTH_TgtOn_Good_RF(index,:)=Data('hPSTH_TgtOn_Good_RF');
 hPSTH_TgtOn_Bad_RF(index,:)=Data('hPSTH_TgtOn_Bad_RF');
        

        hPSTH_SacOn_Value_RF(index,:)=Data('hPSTH_SacOn_Value_RF');
        hPSTH_SacOn_Good_RF(index,:)=Data('hPSTH_SacOn_Good_RF');
        hPSTH_SacOn_Bad_RF(index,:)=Data('hPSTH_SacOn_Bad_RF');

        hPSTH_TgtOn_Value_nonRF(index,:)=Data('hPSTH_TgtOn_Value_nonRF');
        hPSTH_TgtOn_Good_nonRF(index,:)=Data('hPSTH_TgtOn_Good_nonRF');
        hPSTH_TgtOn_Bad_nonRF(index,:)=Data('hPSTH_TgtOn_Bad_nonRF');

        hPSTH_SacOn_Value_nonRF(index,:)=Data('hPSTH_SacOn_Value_nonRF');
        hPSTH_SacOn_Good_nonRF(index,:)=Data('hPSTH_SacOn_Good_nonRF');
        hPSTH_SacOn_Bad_nonRF(index,:)=Data('hPSTH_SacOn_Bad_nonRF');

%}







%% Figures
%% First,focus on contralateral stimulus
FigureStartNum = 200;
FigureIndex=1;
figtitlestr{FigureIndex}='Proportion_Of_SignificantUnits_FEF_Stable_Contra';
fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[50,100, 1000,600],'Name',figtitlestr{FigureIndex});


%figure
subplot(1,2,1)
plot(PSTH_Time_TgtOn_stable,Prop_TgtOn_Value_RF_Fast,'-r','LineWidth',3);
hold on
plot(PSTH_Time_TgtOn_stable,Prop_TgtOn_Value_RF_Others,'-k','LineWidth',3);
%plot(PSTH_Time_TgtOn_stable,Prop_TgtOn_Value_nonRF,'-b','LineWidth',3);
xlabel('Time from target onset (ms)');
ylabel('Proportion of significant neurons');
set(gca,'FontSize',15,'FontWeight','Bold','LineWidth',3);
ylim([0,0.4]);
box off;
legend({'Fast','Others'})
subplot(1,2,2)
plot(PSTH_Time_SacOn_stable,Prop_SacOn_Value_RF_Fast,'-r','LineWidth',3);
hold on
plot(PSTH_Time_SacOn_stable,Prop_SacOn_Value_RF_Others,'-k','LineWidth',3);
%plot(PSTH_Time_SacOn_stable,Prop_SacOn_Value_nonRF,'-b','LineWidth',3);
xlabel('Time from saccade onset(ms)');
%ylabel('Proportion of significant neurons');
set(gca,'FontSize',15,'FontWeight','Bold','LineWidth',3);
set(gca, 'YColor', 'none');
ylim([0,0.4]);
box off





%keyboard
%{
%figure
subplot(2,3,2)
plot(PSTH_Time_TgtOn_stable,Prop_TgtOn_Good_RF,'-r','LineWidth',3);
hold on
plot(PSTH_Time_TgtOn_stable,Prop_TgtOn_Bad_RF,'-b','LineWidth',3);
xlabel('Time from target onset (ms)');
legend({'Good','Bad'});
title('Contra');

subplot(2,3,5)
plot(PSTH_Time_SacOn_stable,Prop_SacOn_Good_RF,'-r','LineWidth',3);
hold on
plot(PSTH_Time_SacOn_stable,Prop_SacOn_Bad_RF,'-b','LineWidth',3);
xlabel('Time from saccade onset (ms)');


%figure
subplot(2,3,3)
plot(PSTH_Time_TgtOn_stable,Prop_TgtOn_Good_nonRF,'-r','LineWidth',3);
hold on
plot(PSTH_Time_TgtOn_stable,Prop_TgtOn_Bad_nonRF,'-b','LineWidth',3);
legend({'Good','Bad'});
xlabel('Time from target onset (ms)');
title('Ipsi');
subplot(2,3,6)
plot(PSTH_Time_SacOn_stable,Prop_SacOn_Good_nonRF,'-r','LineWidth',3);
hold on
plot(PSTH_Time_SacOn_stable,Prop_SacOn_Bad_nonRF,'-b','LineWidth',3);
xlabel('Time from saccade onset (ms)');
%}
FigureStartNum = FigureStartNum+1;
FigureIndex = FigureIndex+1;
figtitlestr{FigureIndex}='Proportion_Of_SignificantUnits_FEF_Stable_FastOthers';
fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[50,100, 1000,600],'Name',figtitlestr{FigureIndex});

b = bar([1,2],[Prop_Value_Fast_Stable,Prop_Value_Other_Stable]);
b.FaceColor='flat';
b.CData=[1,0,0;0,0,0];
ylabel('Proportion of value neurons');
xticklabels({'Fast','Others'});
title('Object Value');
set(gca,'LineWidth',3,'FontSize',15,'FontWeight','bold');

box off

FigureStartNum = FigureStartNum+1;
FigureIndex = FigureIndex+1;
figtitlestr{FigureIndex}='Proportion_Of_SignificantUnits_FEF_OneDR_FastOthers';
fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[50,100, 1000,600],'Name',figtitlestr{FigureIndex});

b = bar([1,2],[Prop_Value_Fast_OneDR,Prop_Value_Other_OneDR]);
b.FaceColor='flat';
b.CData=[1,0,0;0,0,0];
ylabel('Proportion of value neurons');
xticklabels({'Fast','Others'});
title('Location Value');
set(gca,'LineWidth',3,'FontSize',15,'FontWeight','bold');

box off

FigureStartNum = FigureStartNum+1;
FigureIndex = FigureIndex+1;
figtitlestr{FigureIndex}='Proportion_Of_SignificantUnits_FEF_PV_FastOthers';
fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[50,100, 1000,600],'Name',figtitlestr{FigureIndex});

b = bar([1,2],[Prop_Value_Fast_PV,Prop_Value_Other_PV]);
b.FaceColor='flat';
b.CData=[1,0,0;0,0,0];
ylabel('Proportion of value neurons');
xticklabels({'Fast','Others'});
title('PV');
set(gca,'LineWidth',3,'FontSize',15,'FontWeight','bold');

box off




FigureStartNum = FigureStartNum+1;
FigureIndex = FigureIndex+1;
figtitlestr{FigureIndex}='Proportion_Of_SignificantUnits_FEF_Stable_Contra_All';
fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[50,100, 1000,600],'Name',figtitlestr{FigureIndex});

%figure
subplot(1,2,1)
plot(PSTH_Time_TgtOn_stable,Prop_TgtOn_Value_RF_All,'-k','Color','#785EF0','LineWidth',3);

%plot(PSTH_Time_TgtOn_OneDR,Prop_TgtOn_Value_nonRF_oneDR,'-b');
set(gca,'FontSize',15,'FontWeight','Bold','LineWidth',3);
xlabel('Time from target onset (ms)');
ylabel('Proportion of significant neurons');
title('Object value');

ylim([0,0.4]);
box off;
subplot(1,2,2)
plot(PSTH_Time_SacOn_OneDR,Prop_SacOn_Value_RF_oneDR_All,'-k','Color','#785EF0','LineWidth',3);

set(gca,'FontSize',15,'FontWeight','Bold','LineWidth',3);
%hold on
%plot(PSTH_Time_SacOn_OneDR,Prop_SacOn_Value_nonRF_oneDR,'-b');
xlabel('Time from saccade onset (ms)');

set(gca, 'YColor', 'none');
box off
ylim([0,0.4]);





%%
FigureStartNum = FigureStartNum+1;
FigureIndex = FigureIndex+1;
figtitlestr{FigureIndex}='FEF_CellTypePieChart';
fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[50,100, 600,600],'Name',figtitlestr{FigureIndex});

labels = {'Location Only RF', 'Stable Only RF', 'Both RF'};
subplot(1,2,1);
h=piechart([Prop_Sig_LocationOnly_RF,Prop_Sig_StableOnly_RF,Prop_Sig_Both_RF],labels);
colors = {'#785EF0'; 
          '#785EF0';
          '#E69F00'
          }; 
for k = 1:length(h)
    if mod(k, 2) == 0
        h(k).FaceColor = colors(ceil(k/2), :);
    end
end

subplot(1,2,2);
h=piechart([Prop_Sig_LocationOnly_RF_Others,Prop_Sig_StableOnly_RF_Others,Prop_Sig_Both_RF_Others],labels);
colors = {'#785EF0'; 
          '#785EF0';
          '#E69F00'
          }; 
for k = 1:length(h)
    if mod(k, 2) == 0
        h(k).FaceColor = colors(ceil(k/2), :);
    end
end

%{
%% All
FigureStartNum = FigureStartNum+1;
FigureIndex = FigureIndex+1;
figtitlestr{FigureIndex}='Proportion_Of_SignificantUnits_FEF_OneDR_Contra_All';
fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[50,100, 1000,600],'Name',figtitlestr{FigureIndex});

%figure
subplot(1,2,1)
plot(PSTH_Time_TgtOn_OneDR,Prop_TgtOn_Value_RF_oneDR_All,'-k','Color','#009E73','LineWidth',3);

%plot(PSTH_Time_TgtOn_OneDR,Prop_TgtOn_Value_nonRF_oneDR,'-b');
set(gca,'FontSize',15,'FontWeight','Bold','LineWidth',3);
xlabel('Time from target onset (ms)');
ylabel('Proportion of significant neurons');

ylim([0,0.4]);
box off;
subplot(1,2,2)
plot(PSTH_Time_SacOn_OneDR,Prop_SacOn_Value_RF_oneDR_All,'-k','Color','#009E73','LineWidth',3);

set(gca,'FontSize',15,'FontWeight','Bold','LineWidth',3);
%hold on
%plot(PSTH_Time_SacOn_OneDR,Prop_SacOn_Value_nonRF_oneDR,'-b');
xlabel('Time from saccade onset (ms)');
title('Location value');
set(gca, 'YColor', 'none');
box off
ylim([0,0.4]);
%}
%{
%% Distribution of value separation starting time 
FigureStartNum = FigureStartNum+1;
FigureIndex = FigureIndex+1;
figtitlestr{FigureIndex}='Proportion_Of_SignificantUnits_FEF_OneDR_Contra_All';
fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[50,100, 400,600],'Name',figtitlestr{FigureIndex});

subplot(2,1,1)
bar(bin_edges(1:end-1),counts_Stable_All,'FaceColor','#785EF0');
box off
ylabel('Number of units');
set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');


subplot(2,1,2)
bar(bin_edges(1:end-1),counts_OneDR_All,'FaceColor','#009E73');
box off
xlabel('Starting time of value separation');
set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
%}

%% Population Average of PSTH_Stable
FigureStartNum = FigureStartNum+1;
FigureIndex = FigureIndex+1;
figtitlestr{FigureIndex}='PSTH_Fast_Stable_RF';
fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[50,100, 400,600],'Name',figtitlestr{FigureIndex});

subplot(1,2,1)
shadedErrorBar(PSTH_Time_TgtOn_stable,PSTH_TgtOn_Good_RF_Stable_Fast_Value_GoodPrefer,PSTH_TgtOn_Good_RF_Stable_Fast_Value_GoodPrefer_Sem,'lineProps',{[1    0    0]},'transparent',1,'patchSaturation',0.3);
hold on
shadedErrorBar(PSTH_Time_TgtOn_stable,PSTH_TgtOn_Bad_RF_Stable_Fast_Value_GoodPrefer,PSTH_TgtOn_Bad_RF_Stable_Fast_Value_GoodPrefer_Sem,'lineProps',{[0    0    1]},'transparent',1,'patchSaturation',0.3);
set(findall(gca, 'Type', 'Line'),'LineWidth',3);
xlabel('Time from target onset (ms)');
 ylabel('z-score');
 title('High value prefer');
 set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');

subplot(1,2,2)
shadedErrorBar(PSTH_Time_TgtOn_stable,PSTH_TgtOn_Good_RF_Stable_Fast_Value_BadPrefer,PSTH_TgtOn_Good_RF_Stable_Fast_Value_BadPrefer_Sem,'lineProps',{[1    0    0]},'transparent',1,'patchSaturation',0.3);
hold on
shadedErrorBar(PSTH_Time_TgtOn_stable,PSTH_TgtOn_Bad_RF_Stable_Fast_Value_BadPrefer,PSTH_TgtOn_Bad_RF_Stable_Fast_Value_BadPrefer_Sem,'lineProps',{[0    0    1]},'transparent',1,'patchSaturation',0.3);
 set(findall(gca, 'Type', 'Line'),'LineWidth',3);
 xlabel('Time from target onset (ms)');
 title('Low value prefer');
 set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');



%% Passive view
FigureStartNum = FigureStartNum+1;
FigureIndex = FigureIndex+1;
figtitlestr{FigureIndex}='PSTH_PV_RF';
fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[50,100, 400,600],'Name',figtitlestr{FigureIndex});
plot(Time_TgtOn_PV,Prop_Sig_PV_Fast,'-r','LineWidth',3);
hold on
plot(Time_TgtOn_PV,Prop_Sig_PV_Others,'-k','LineWidth',3);
legend({'Fast','Others'});
xlabel('Time from target onset(ms)');
ylabel('Proportion of units');
set(gca,'FontSize',15,'FontWeight','Bold','LineWidth',3);


box off;

%Population average psth
FigureStartNum = FigureStartNum+1;
FigureIndex = FigureIndex+1;
figtitlestr{FigureIndex}='PSTH_PV_RF';
fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[50,100, 400,600],'Name',figtitlestr{FigureIndex});


subplot(1,2,1)
shadedErrorBar(Time_TgtOn_PV,PSTH_TgtOn_Good_RF_PV_Fast_Sig_Good_Mean,PSTH_TgtOn_Good_RF_PV_Fast_Sig_Good_Sem,'lineProps',{[1    0    0]},'transparent',1,'patchSaturation',0.3)
hold on
shadedErrorBar(Time_TgtOn_PV,PSTH_TgtOn_Bad_RF_PV_Fast_Sig_Good_Mean,PSTH_TgtOn_Bad_RF_PV_Fast_Sig_Good_Sem,'lineProps',{[0    0    1]},'transparent',1,'patchSaturation',0.3)
set(findall(gca, 'Type', 'Line'),'LineWidth',3);
xlabel('Time from target onset (ms)');
 ylabel('z-score');
 title('High value prefer');
 set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');

subplot(1,2,2)
shadedErrorBar(Time_TgtOn_PV,PSTH_TgtOn_Good_RF_PV_Fast_Sig_Bad_Mean,PSTH_TgtOn_Good_RF_PV_Fast_Sig_Bad_Sem,'lineProps',{[1    0    0]},'transparent',1,'patchSaturation',0.3)
hold on
shadedErrorBar(Time_TgtOn_PV,PSTH_TgtOn_Bad_RF_PV_Fast_Sig_Bad_Mean,PSTH_TgtOn_Bad_RF_PV_Fast_Sig_Bad_Sem,'lineProps',{[0    0    1]},'transparent',1,'patchSaturation',0.3)
set(findall(gca, 'Type', 'Line'),'LineWidth',3);
xlabel('Time from target onset (ms)');
 ylabel('z-score');
 title('Low value prefer');
 set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');


%Population average psth
FigureStartNum = FigureStartNum+1;
FigureIndex = FigureIndex+1;
figtitlestr{FigureIndex}='PSTH_OneDR_RF';
fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[50,100, 400,600],'Name',figtitlestr{FigureIndex});


subplot(1,2,1)
shadedErrorBar(PSTH_Time_TgtOn_OneDR,PSTH_TgtOn_Good_RF_OneDR_Fast_Sig_Good_Mean,PSTH_TgtOn_Good_RF_OneDR_Fast_Sig_Good_Sem,'lineProps',{[1    0    0]},'transparent',1,'patchSaturation',0.3)
hold on
shadedErrorBar(PSTH_Time_TgtOn_OneDR,PSTH_TgtOn_Bad_RF_OneDR_Fast_Sig_Good_Mean,PSTH_TgtOn_Bad_RF_OneDR_Fast_Sig_Good_Sem,'lineProps',{[0    0    1]},'transparent',1,'patchSaturation',0.3)
set(findall(gca, 'Type', 'Line'),'LineWidth',3);
xlabel('Time from target onset (ms)');
 ylabel('z-score');
 title('High value prefer');
 set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');

subplot(1,2,2)
shadedErrorBar(PSTH_Time_TgtOn_OneDR,PSTH_TgtOn_Good_RF_OneDR_Fast_Sig_Bad_Mean,PSTH_TgtOn_Good_RF_OneDR_Fast_Sig_Bad_Sem,'lineProps',{[1    0    0]},'transparent',1,'patchSaturation',0.3)
hold on
shadedErrorBar(PSTH_Time_TgtOn_OneDR,PSTH_TgtOn_Bad_RF_OneDR_Fast_Sig_Bad_Mean,PSTH_TgtOn_Bad_RF_OneDR_Fast_Sig_Bad_Sem,'lineProps',{[0    0    1]},'transparent',1,'patchSaturation',0.3)
set(findall(gca, 'Type', 'Line'),'LineWidth',3);
xlabel('Time from target onset (ms)');
 ylabel('z-score');
 title('Low value prefer');
 set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');


%% AUC as time
% Stable 
FigureStartNum = FigureStartNum+1;
FigureIndex = FigureIndex+1;
figtitlestr{FigureIndex}='AUC_PSTH_stable';
fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[50,100, 400,600],'Name',figtitlestr{FigureIndex});


subplot(1,2,1)
shadedErrorBar(PSTH_Time_TgtOn_stable,AUCPSTH_TgtOn_RF_Stable_Fast_Value_GoodPrefer,AUCPSTH_TgtOn_RF_Stable_Fast_Value_GoodPrefer_Sem,'lineProps',{[1    0    0]},'transparent',1,'patchSaturation',0.3)
hold on
shadedErrorBar(PSTH_Time_TgtOn_stable,AUCPSTH_TgtOn_RF_Stable_Others_Value_GoodPrefer,AUCPSTH_TgtOn_RF_Stable_Others_Value_GoodPrefer_Sem,'lineProps',{[0    0    0]},'transparent',1,'patchSaturation',0.3)

set(findall(gca, 'Type', 'Line'),'LineWidth',3);
legend({'Fast','Others'});
xlabel('Time from target onset (ms)');
 ylabel('AUC');
 title('High value prefer');
 set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');

subplot(1,2,2)
shadedErrorBar(PSTH_Time_TgtOn_stable,AUCPSTH_TgtOn_RF_Stable_Fast_Value_BadPrefer,AUCPSTH_TgtOn_RF_Stable_Fast_Value_BadPrefer_Sem,'lineProps',{[0    0    1]},'transparent',1,'patchSaturation',0.3)
hold on
shadedErrorBar(PSTH_Time_TgtOn_stable,AUCPSTH_TgtOn_RF_Stable_Others_Value_BadPrefer,AUCPSTH_TgtOn_RF_Stable_Others_Value_BadPrefer_Sem,'lineProps',{[0    0    0]},'transparent',1,'patchSaturation',0.3)

set(findall(gca, 'Type', 'Line'),'LineWidth',3);
legend({'Fast','Others'});
xlabel('Time from target onset (ms)');
 ylabel('AUC');
 title('Low value prefer');
 set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');

%OneDR

FigureStartNum = FigureStartNum+1;
FigureIndex = FigureIndex+1;
figtitlestr{FigureIndex}='AUC_PSTH_oneDR';
fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[50,100, 400,600],'Name',figtitlestr{FigureIndex});


subplot(1,2,1)
shadedErrorBar(PSTH_Time_TgtOn_OneDR,AUCPSTH_TgtOn_RF_OneDR_Fast_Sig_Good_mean,AUCPSTH_TgtOn_RF_OneDR_Fast_Sig_Good_sem,'lineProps',{[1    0    0]},'transparent',1,'patchSaturation',0.3)
hold on
shadedErrorBar(PSTH_Time_TgtOn_OneDR,AUCPSTH_TgtOn_RF_OneDR_Others_Sig_Good_mean,AUCPSTH_TgtOn_RF_OneDR_Others_Sig_Good_sem,'lineProps',{[0    0    0]},'transparent',1,'patchSaturation',0.3)

set(findall(gca, 'Type', 'Line'),'LineWidth',3);
legend({'Fast','Others'});
xlabel('Time from target onset (ms)');
 ylabel('AUC');
 title('High value prefer');
 set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');

subplot(1,2,2)
shadedErrorBar(PSTH_Time_TgtOn_OneDR,AUCPSTH_TgtOn_RF_OneDR_Fast_Sig_Bad_mean,AUCPSTH_TgtOn_RF_OneDR_Fast_Sig_Bad_sem,'lineProps',{[1    0    0]},'transparent',1,'patchSaturation',0.3)
hold on
shadedErrorBar(PSTH_Time_TgtOn_OneDR,AUCPSTH_TgtOn_RF_OneDR_Others_Sig_Bad_mean,AUCPSTH_TgtOn_RF_OneDR_Others_Sig_Bad_sem,'lineProps',{[0    0    0]},'transparent',1,'patchSaturation',0.3)

set(findall(gca, 'Type', 'Line'),'LineWidth',3);
legend({'Fast','Others'});
xlabel('Time from target onset (ms)');
 ylabel('AUC');
 title('Low value prefer');
 set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');

% Passive Viewing 
FigureStartNum = FigureStartNum+1;
FigureIndex = FigureIndex+1;
figtitlestr{FigureIndex}='AUC_PSTH_PV';
fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[50,100, 400,600],'Name',figtitlestr{FigureIndex});


subplot(1,2,1)
shadedErrorBar(Time_TgtOn_PV,AUCPSTH_TgtOn_RF_PV_Fast_Sig_Good_mean,AUCPSTH_TgtOn_RF_PV_Fast_Sig_Good_sem,'lineProps',{[1    0    0]},'transparent',1,'patchSaturation',0.3)
hold on
shadedErrorBar(Time_TgtOn_PV,AUCPSTH_TgtOn_RF_PV_Others_Sig_Good_mean,AUCPSTH_TgtOn_RF_PV_Others_Sig_Good_sem,'lineProps',{[0    0    0]},'transparent',1,'patchSaturation',0.3)

set(findall(gca, 'Type', 'Line'),'LineWidth',3);
legend({'Fast','Others'});
xlabel('Time from target onset (ms)');
 ylabel('AUC');
 title('High value prefer');
 set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');

subplot(1,2,2)
shadedErrorBar(Time_TgtOn_PV,AUCPSTH_TgtOn_RF_PV_Fast_Sig_Bad_mean,AUCPSTH_TgtOn_RF_PV_Fast_Sig_Bad_sem,'lineProps',{[1    0    0]},'transparent',1,'patchSaturation',0.3)
hold on
shadedErrorBar(Time_TgtOn_PV,AUCPSTH_TgtOn_RF_PV_Others_Sig_Bad_mean,AUCPSTH_TgtOn_RF_PV_Others_Sig_Bad_sem,'lineProps',{[0    0    0]},'transparent',1,'patchSaturation',0.3)


set(findall(gca, 'Type', 'Line'),'LineWidth',3);
legend({'Fast','Others'});
xlabel('Time from target onset (ms)');
 ylabel('AUC');
 title('Low value prefer');
 set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');


keyboard


%{
%% Proportion contralateral
FigureStartNum = FigureStartNum+1;
FigureIndex = FigureIndex+1;
figtitlestr{FigureIndex}='Proportion_Of_SignificantUnits_FEF_TwoTasks_Contra';
fig{FigureIndex}=PrepareFigure(FigureStartNum,'w',[50,100, 1200,600],'Name',figtitlestr{FigureIndex});

subplot(1,4,1)
plot(PSTH_Time_TgtOn_OneDR,Prop_TgtOn_Value_RF_Fast_stableOnly,'-r','LineWidth',3)
hold on
plot(PSTH_Time_TgtOn_OneDR,Prop_TgtOn_Value_RF_Fast_OneDROnly,'-b','LineWidth',3)
plot(PSTH_Time_TgtOn_OneDR,Prop_TgtOn_Value_RF_Fast_Both,'-g','LineWidth',3);
legend({'Object Value Only','Location Value Only','Both'});
set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
xlabel('Time from target onset');
ylabel('Proportion of significant units');
title('For Fast group');
ylim([0,0.3]);
box off;

subplot(1,4,2)
plot(PSTH_Time_SacOn_OneDR,Prop_SacOn_Value_RF_Fast_stableOnly,'-r','LineWidth',3)
hold on
plot(PSTH_Time_SacOn_OneDR,Prop_SacOn_Value_RF_Fast_OneDROnly,'-b','LineWidth',3)
plot(PSTH_Time_SacOn_OneDR,Prop_SacOn_Value_RF_Fast_Both,'-g','LineWidth',3);
%legend({'Object Value Only','Location Value Only','Both'});
set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
xlabel('Time from saccade onset');
%ylabel('Proportion of significant units');
set(gca, 'YColor', 'none');
ylim([0,0.3]);
box off;

subplot(1,4,3)
plot(PSTH_Time_TgtOn_OneDR,Prop_TgtOn_Value_RF_Others_stableOnly,'-r','LineWidth',3)
hold on
plot(PSTH_Time_TgtOn_OneDR,Prop_TgtOn_Value_RF_Others_OneDROnly,'-b','LineWidth',3)
plot(PSTH_Time_TgtOn_OneDR,Prop_TgtOn_Value_RF_Others_Both,'-g','LineWidth',3);
set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
xlabel('Time from target onset');
ylabel('Proportion of significant units');
title('For Others group');
ylim([0,0.3]);
box off;

subplot(1,4,4)
plot(PSTH_Time_SacOn_OneDR,Prop_SacOn_Value_RF_Others_stableOnly,'-r','LineWidth',3)
hold on
plot(PSTH_Time_SacOn_OneDR,Prop_SacOn_Value_RF_Others_OneDROnly,'-b','LineWidth',3)
plot(PSTH_Time_SacOn_OneDR,Prop_SacOn_Value_RF_Others_Both,'-g','LineWidth',3);
set(gca,'LineWidth',3,'FontSize',15,'FontWeight','Bold');
set(gca, 'YColor', 'none');
xlabel('Time from saccade onset');
ylim([0,0.3]);
box off;

%{
subplot(2,3,2)
plot(PSTH_Time_TgtOn_OneDR,Prop_TgtOn_Good_RF_oneDR,'-r');
hold on
plot(PSTH_Time_TgtOn_OneDR,Prop_TgtOn_Bad_RF_oneDR,'-b');
xlabel('Time from target onset (ms)');
title('Contra');
legend({'Good','Bad'})
subplot(2,3,5)
plot(PSTH_Time_SacOn_OneDR,Prop_SacOn_Good_RF_oneDR,'-r');
hold on
plot(PSTH_Time_SacOn_OneDR,Prop_SacOn_Bad_RF_oneDR,'-b');

xlabel('Time from saccade onset (ms)');



subplot(2,3,3)
plot(PSTH_Time_TgtOn_OneDR,Prop_TgtOn_Good_nonRF_oneDR,'-r');
hold on
plot(PSTH_Time_TgtOn_OneDR,Prop_TgtOn_Bad_nonRF_oneDR,'-b');
xlabel('Time from target onset (ms)');
legend({'Good','Bad'});
title('Ipsi');
subplot(2,3,6)
plot(PSTH_Time_SacOn_OneDR,Prop_SacOn_Good_nonRF_oneDR,'-r');
hold on
plot(PSTH_Time_SacOn_OneDR,Prop_SacOn_Bad_nonRF_oneDR,'-b')

xlabel('Time from saccade onset (ms)');
%}
%}
keyboard





end