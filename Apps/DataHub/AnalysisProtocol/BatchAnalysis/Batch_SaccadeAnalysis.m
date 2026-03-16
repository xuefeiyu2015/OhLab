function Batch_SaccadeAnalysis(Data,  StartFile, EndFile,ShowFigureFlag,OutputFlag);
%Batch analysis for saccades in ITI, delay saccade & memory saccade; by
%Xuefei Yu 12272022
ECode;
BatchFileName=Data.BatchFileName;
%Batch data file path
FilesPath=Data.ResultFilePath(StartFile: EndFile);
FilesName=Data.FileName;
RecordDate=Data.RecordDate(StartFile: EndFile);
FileNum=EndFile-StartFile+1;
TotalNeuronNum=sum(cellfun(@numel,FilesPath));
ChannelID=Data.ChannelNumber;


%Keynames
 keysname_Saccade ={'PSTH_FR_SacOn_ITI_Time','PSTH_FR_SacOn_ITI','PSTH_FR_SacOff_ITI_Time','PSTH_FR_SacOff_ITI',...
        'SacVec','SacVectTuning_PreDurPost_mean','SacVectTuning_PreDurPost_sem','SacVectTuning_p','Pref_Vect','Pref_Vect_Length',...
        'SacOnLocMean','SacOnLocSem'}; 
 keysname_DS={
    'PSTH_FR_TgtOn_Time','PSTH_FR_TgtOn',...
    'PSTH_FR_SacOn_Time','PSTH_FR_SacOn',...
    'PSTH_FR_SacOff_ITI_Time','PSTH_FR_SacOff_ITI',...
    'SacVec','SacVectTuning_PreDurPost_mean','SacVectTuning_PreDurPost_sem','SacVectTuning_p','Pref_Vect','Pref_Vect_Length',...
    'TgtVec','TgtVectTuning_mean','TgtVectTuning_sem','TgtVectTuning_p','Pref_TgtVect','Pref_TgtVect_Length',...
    'SacOnLocMean','SacOnLocSem'};

 keysname_MS={'UniqueAngle','PSTH_Time_TgtOn','TgtOn_Loc_Mean','TgtOn_Loc_Sem',...
    'PSTH_Time_SacOn','SacOn_Loc_Mean','SacOn_Loc_Sem'};



%Load file
ChannelFull=[];
FileIndex=1;
for i=1:length(FilesPath)
  clear OutputData;
  ChannelFull=[ChannelFull,ChannelID(i,1):ChannelID(i,2)];
  
  for j=1:length(FilesPath{i})
      currFileName=FilesPath{i}{j};
      load(currFileName);
     %OutputDataTemp=OutputData; 
      ExploreSac{FileIndex} = {};
      SponSac{FileIndex} = {};
      MemSac{FileIndex} = {};
      DSSac{FileIndex} = {};

      %For each plot
      PSTH_timeOn_ES=[];
      PSTH_SacOn_ES =[];

      PSTH_timeOff_ES=[];
      PSTH_SacOff_ES=[];




      PSTH_TimeOn_SS=[];
      PSTH_TimeOn_DS=[];
      PSTH_TimeOn_MS=[];



     if isfield(OutputData,'SaccadeAnalysis')
         %For video free view and delay saccade
        
        
         for m = 1:length(OutputData.SaccadeAnalysis)
             Task_tmp=OutputData.SaccadeAnalysis(m).TaskCode;
          if Task_tmp == VIDEO_FREEVIEW_TASK
              tmp = OutputData.SaccadeAnalysis(m).DataStamp;
              ExploreSac{FileIndex} = values(tmp, keysname_Saccade);

              PSTH_timeOn_ES =ExploreSac{FileIndex}{1};
              PSTH_SacOn_ES = ExploreSac{FileIndex}{2};

           %   PSTH_timeOff_ES =ExploreSac{FileIndex}{3};
           %   PSTH_SacOff_ES = ExploreSac{FileIndex}{4};



             
              

              


          elseif Task_tmp == SACSTIM_TASK || Task_tmp == SCC_TASK
              tmp = OutputData.SaccadeAnalysis(m).DataStamp;
              SponSac{FileIndex}= values(tmp, keysname_Saccade);

              PSTH_timeOn_SS =SponSac{FileIndex}{1};
              PSTH_SacOn_SS = SponSac{FileIndex}{2};

             % PSTH_timeOff_SS =SponSac{FileIndex}{3};
             % PSTH_SacOff_SS = SponSac{FileIndex}{4};



          end

        end

     end

     if isfield(OutputData,'DelaySaccadeTuning')
         %For delay saccade
         tmp = OutputData.DelaySaccadeTuning.DataStamp;
         DSSac{FileIndex}= values(tmp, keysname_DS);

     end

     if isfield(OutputData,'MemorySaccade')
         %For delay saccade
         tmp = OutputData.MemorySaccade.DataStamp;
         MemSac{FileIndex}= values(tmp, keysname_MS);

     end

     
      FileIndex=FileIndex+1;
      
  end


  %Compare the overall saccade psth;




end
















%Assume that keys are all the same for each file, use the first file to
%abstract keys for all files;
KeysLib=keys(FourBlockSceneCompare{1});
AllValues=cellfun(@(x) values(x),FourBlockSceneCompare,'uniform',0);
%Reorganize the value according to each key
for i=1:length(KeysLib)
  %  currKey=KeysLib{i};
   KeyValues{i}=cellfun(@(x) x{i},AllValues,'uniform',0);
 
end

SSISceneKey={'Time_PSTH_SceneOn','Time_PSTH_SceneOff','SSI_Group','SSI_BlockIndex','SceneTypeIndex','SSI_GroupOff',' p_SSI_GroupOff',' p_SSI_GroupOff',... 
    };

%Time 
SSISceneOnTime=KeyValues{ismember(KeysLib,SSISceneKey{1})};
 
 SSISceneOn_TimeMat=fillinnancol(SSISceneOnTime);
 SSISceneOn_TimeMean=nanmean(SSISceneOn_TimeMat,1);
 
 SSISceneOffTime=KeyValues{ismember(KeysLib,SSISceneKey{2})};
 
 SSISceneOff_TimeMat=fillinnancol(SSISceneOffTime);
 SSISceneOff_TimeMean=nanmean(SSISceneOff_TimeMat,1);
 
%SSI
SceneGroup={[1000:1007],[1010:1017],[1054:1057,1060:1063]};
SceneGroupStr={'FlexibleScenes','StableScenes','DuelScenes'};

SSIOnData=KeyValues{ismember(KeysLib,SSISceneKey{3})};
SSIOnBlockData=KeyValues{ismember(KeysLib,SSISceneKey{4})};
SSISceneTypeData=KeyValues{ismember(KeysLib,SSISceneKey{5})};
%SSI_SceneOnTime=cellfun(@(x) cellfun(@(y) y(:,[2:3]),x,'uniform',0),KeyValues('Time_PSTH_SceneOn'),'uniform',0);%Normalize
% PSTHSceneData_Time=cellfun(@(x) cellfun(@(y) y(:,1),x,'uniform',0),KeyValues(Sel_tmp),'uniform',0);


end