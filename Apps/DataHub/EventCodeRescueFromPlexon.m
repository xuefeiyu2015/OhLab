%Rescue eventcode from Plexon file
ECode;
%by Xuefei Yu, 05242021
fileNamePack='Robin052321FORAGEXPLOITL1';
%Get File Information
num = regexp(fileNamePack, '\d+', 'match');
RecordDateOriginal=cell2mat(num(1));
RecordDate=datetime(RecordDateOriginal,'InputFormat','MMddyy');
FileNum=cell2mat(num (2));
%To get the monkey name
RemainFiles=erase(fileNamePack,RecordDateOriginal);
Delimiter=find(isstrprop(RemainFiles,'upper')==1);
MonkeyName=RemainFiles(Delimiter(1):Delimiter(2)-1);

TaskCode=329;

PlexonFileName=strcat(MonkeyName,RecordDateOriginal,'N',num2str(FileNum),'T',num2str(TaskCode),'.mat');


%Setup output directory
Workingdirectory=pwd;
%Set up output path according to the system
OperationSystem = computer;%Get the operation system information

if strcmp(OperationSystem(1:3),"PCW")  
    %PC
     MarkerFolder='DataHub';
     SysDelimeter='\';
elseif strcmp(OperationSystem(1:3),"MAC")  
    %Mac
     MarkerFolder='DataAnalysis';
     SysDelimeter='/';
end
 
PlexonDataFolder='PlexonData';
PlexonParsedDataFolder='ParsedData';

Flag=strfind(Workingdirectory,MarkerFolder);
PlexonPath=strcat(Workingdirectory(1:Flag+length(MarkerFolder)),PlexonDataFolder,SysDelimeter,PlexonParsedDataFolder,SysDelimeter);

ExtraSpikePath=strcat(PlexonPath,RecordDateOriginal,SysDelimeter);

    
   
    cd(ExtraSpikePath);
   
        
       load(PlexonFileName);
       %Measure the basic information
       TrialNum_Plexon=su_task.TrialNum;
       SpikeTime_Extra_Pre=su_task.SpikeTiming;       
       SpikeWaveForm_Med=su_task.MedWaveform;
       Event_Plexon_ori=su_task.Event;
       EventTiming_Plexon=su_task.EventTiming;
       
      
      
       TrialIndex_Plexon=su_task.TrialIndex;
       
       %Transfer event code into matrix
       MaxNumOfEvent=max(cellfun(@(x) numel(x),Event_Plexon));

      Event_Plexon=cellfun(@(x) x',Event_Plexon_ori,'uniform',0);

      EventDataChannel=cell2mat(cellfun(@(x) [x,NaN*ones(1,MaxNumOfEvent-numel(x))],Event_Plexon,'uniform',0));
      
      %Good trials
      LastEvent=cellfun(@(x) x(end),Event_Plexon);
      GoodTrialsFlag=LastEvent==CORRECTCD;
      %Reproduce scene index
      SceneOnMarker=SCENE_ON;
      SceneOnLoc=FindOutMarkerLoc(EventDataChannel,SceneOnMarker);
      SceneIndexLoc=SceneOnLoc+1;
      SceneIndex=FindOutMarker(EventDataChannel,SceneIndexLoc);
      SceneIndex=SceneIndex-800;
      
      %Reproduce object index
      ObjMarker=FIXOFF;
      FixOffMarker=FindOutMarkerLoc(EventDataChannel,ObjMarker);
      ObjLoc=FixOffMarker+1;
      ObjectType=FindOutMarker(EventDataChannel,ObjLoc);
      
      %Remove those which are not objecttype code
      ObjectTypeMarker=[GOODOBJ  ,BADOBJ];
      
      ObjectType=RemoveMarker(ObjectType,ObjectTypeMarker);
      
      ObjIndexLoc=FixOffMarker+2;
      ObjectIndex=FindOutMarker(EventDataChannel,ObjIndexLoc);
      ObjectIndexMarker=[800:899];
      ObjectIndex=RemoveMarker(ObjectIndex,ObjectIndexMarker);
      ObjectIndex=ObjectIndex-800;
    
      TrialIndex=EventDataChannel(:,4)-1000;
      BlockType=NaN*ones(size(TrialIndex,1),1);
      BlockType(ismember(SceneIndex,[37:40]))=11;%s'FixLoc,RandSceneSeq';
      BlockType(ismember(SceneIndex,[41:44]))=1;%'RandLoc,RandSceneSeq'
      DuelSceneTrialIndicator=find(ismember(SceneIndex,[45:48])==1);
      DuelSceneTrial=TrialIndex(DuelSceneTrialIndicator);
      
      BlockType(ismember(SceneIndex,[45:48])'&TrialIndex<=185)=11;%'FixLoc,RandSceneSeq';
      BlockType(ismember(SceneIndex,[45:48])'&TrialIndex>185)=1;%'RandLoc,RandSceneSeq'
      
      
      %Spikes
      SpikeTime=su_task.SpikeTiming; 
      
      %Spike1
     s
      
       
       keyboard
       
      % A=TrialIndex_Plexon;
      % B=TrialIndexBlip;
       %{
       %Load lfp
       lfpName=strcat(PlexonFileName(1:end-4),'lfp');
       load(lfpName);   
       LFP=lfp_task.LFP;
       LFPTime_Pre=lfp_task.LFPTime;
       LFPSamplingRate=lfp_task.LFPSamplingRate;
      % SpikeWaveFormTime_Extra=su_task.WaveformTime;
      % SR=su_task.SpikeWaveFormSamplingRate;
      %
      
       %First check the correspondance of the trial number
       if TrialNum_Plexon==TotalTrials
           %Align Event Code
           try
           Plexon_StartMarker=PLX_ON;
           AlignTimeBlip=cellfun(@(x,y) y(find(x==Plexon_StartMarker,1,'first')),EventCodeForAll,EventTimeForAll);%in ms
           AlignTimePlexon=cellfun(@(x,y) y(find(x==Plexon_StartMarker,1,'first'))*1000,Event_Plexon,EventTiming_Plexon);%transform s into ms
           catch
               disp('Align code switch to ENABLE');
                Plexon_StartMarker=ENABLECD;
                AlignTimeBlip=cellfun(@(x,y) y(find(x==Plexon_StartMarker,1,'first')),EventCodeForAll,EventTimeForAll);%in ms
                AlignTimePlexon=cellfun(@(x,y) y(find(x==Plexon_StartMarker,1,'first'))*1000,Event_Plexon,EventTiming_Plexon);%transform s into ms
               
           end
           
          
          %Change to the common time
         RefTimeDiff=AlignTimeBlip-AlignTimePlexon;
          Index=num2cell([1:TotalTrials]');
          
          SpikeTime_Extra=cellfun(@(x) cellfun(@(y,z) y'*1000+RefTimeDiff(z),x, Index,'uniform',0),SpikeTime_Extra_Pre,'uniform',0);%in ms
          %Reorder the spike time according to the recording channel number;
       %   ClusterID_ori=[fillClustertmp;cell2mat(su_task.clusterID)'];
       %  ChannelID_ori=[fillChanneltmp;cell2mat(su_task.channelID)'];
       
       
            ClusterID_ori=su_task.clusterID;
            ChannelID_ori=su_task.channelID;
     
        
       
       
         ClusterID=[fillClustertmp;ClusterID_ori];
         ChannelID=[fillChanneltmp;ChannelID_ori];
          
          LFPTime=cellfun(@(x,y) x*1000+RefTimeDiff(y),LFPTime_Pre,Index,'uniform',0);
          LFPValidIndex=1;
        
          %Event Marker Match Check
         EventTime_AfterAlign_Plexon=cellfun(@(x,y) x*1000+RefTimeDiff(y),EventTiming_Plexon,Index,'uniform',0);%in ms
        
          PlexonMatch=cellfun(@EventMatch, Event_Plexon,EventCodeForAll,'uniform',0);
          BlipMatch=cellfun(@EventMatch,EventCodeForAll,Event_Plexon,'uniform',0);
          
          
          %Event check
          PlexonEventCommon_Check=cellfun(@(x,y) y(logical(x))',PlexonMatch,Event_Plexon,'uniform',0);
          BlipEventCommon_Check=cellfun(@(x,y) y(logical(x))',BlipMatch,EventCodeForAll,'uniform',0);
          
          
          PlexonEventTime_Check=cellfun(@(x,y) y(logical(x))',PlexonMatch,EventTime_AfterAlign_Plexon,'uniform',0);
          BlipEventTime_Check=cellfun(@(x,y) y(logical(x)),BlipMatch,EventTimeForAll,'uniform',0);
          EventTimeDifferenceCheck=cellfun(@(x,y) abs(x-y),PlexonEventTime_Check,BlipEventTime_Check,'uniform',0);
        %  A=cellfun(@max,EventTimeDifferenceCheck);
        %  B=TrialIndex_Plexon;
        %   C=[A,B];
           
          LargestTimeDescrepency=max(cellfun(@max,EventTimeDifferenceCheck));%in ms
          DescrepencyTolerance=5;%in ms, the maximum temporal descrepency between Blip events and the Plexon events
          if LargestTimeDescrepency>DescrepencyTolerance
              disp(sprintf('Error!!%d ms descrepency between Blip and Plexon,above %d ms please check where is wrong!!',LargestTimeDescrepency,DescrepencyTolerance));
     
              return;
          end
          
          disp(sprintf('Largest temporal descrepency is %d ms',LargestTimeDescrepency));
          
          %Combine Plexon spikes with original spikes
          for i=1:size(SpikeTime_Extra,2)
              
          SpikeTime(:,PlexonStartingChannel+i-1)=SpikeTime_Extra{i};
          end
          %{
          %Packing the waveform sample to accomondate the old styles using
          %the median waveform
        %  WaveFormSampleNum=1000;
          NumOfFake=2;
          %SpikeWaveFormSample_Wform=cell2mat(cellfun(@(x) x(1:min(size(x,1),WaveFormSampleNum),:),SpikeWaveForm_Extra,'uniform',0)');
          
          SpikeWaveFormSample_Wform=SpikeWaveForm_Med;
          SpikeWaveFormSample_Num=size(SpikeWaveForm_Med,1);
         % SpikeWaveFormSample_Num=cell2mat(cellfun(@(x) min(size(x,1),WaveFormSampleNum),SpikeWaveForm_Extra,'uniform',0)');
         SpikeWaveFormSample_Num_Array=1:sum(SpikeWaveFormSample_Num)+NumOfFake;
       %   SpikeWaveFormSample_Num_Array=1:sum(SpikeWaveFormSample_Num)+NumOfFake;
          SpikeWaveFormSample_Type=-1*ones(sum(SpikeWaveFormSample_Num)+NumOfFake,1); %Add some fake data about the unsorted spikes
          for i=1:length(SpikeWaveFormSample_Num)
              SpikeWaveFormSample_Type((i-1)*SpikeWaveFormSample_Num(max(i-1,1))+1:i*SpikeWaveFormSample_Num(i))=-1+i;
    
          end
          %Add some fake data about the unsorted spikes
         SpikeWaveFormSample_Wform=[SpikeWaveFormSample_Wform;NaN*ones(NumOfFake,size(SpikeWaveFormSample_Wform,2))];
         SampleSR=SR*ones(sum(SpikeWaveFormSample_Num)+NumOfFake,1);
         
         SpikeWaveFormSample.Wform=SpikeWaveFormSample_Wform;
         SpikeWaveFormSample.SR=SampleSR;
         SpikeWaveFormSample.Num=SpikeWaveFormSample_Num_Array;
         SpikeWaveFormSample.Type=SpikeWaveFormSample_Type;
         %}
         
       else
           %Try to rescue with trialid
           %Match trial index from Blip and Plexon
          %TrialIndex_Plexon
           disp('TrialNum Dismatch!!!');
           msgbox('TrialNum Dismatch!!Please check the Plexon parsing manually!','PlexonParsingError','error');
       end
       
       
      
        
   %}
    
function loc=FindOutMarkerLoc(EventCode,Marker);
        for i=1:size(EventCode,1)
            valid=ismember(EventCode(i,:),Marker);
            if sum(valid)
                loc{i}=find(valid==1);
            else
                loc{i}=NaN;
            end
            
            
        end
        %Change into a matrix
        maxnum=max(cellfun(@numel,loc));
        if maxnum>1
        loc=cell2mat(cellfun(@(x) [x,NaN*ones(1,maxnum-numel(x))]',loc,'uniform',0))';
        else
            loc=cell2mat(loc)';
        end
      
end
    
function marker=FindOutMarker(EventCode,loc);
        for i=1:size(EventCode,1)
            if sum(~isnan(loc(i,:)))>0
                marker{i}=EventCode(i,loc(i,~isnan(loc(i,:))));
                
            else
                marker{i}=NaN;
            end
            
        end
        %Change into a matrix
        maxnum=max(cellfun(@(x) size(x,2),marker));
        if maxnum>1
        marker=cell2mat(cellfun(@(x) [x,NaN*ones(1,maxnum-size(x,2))]',marker,'uniform',0))';
        else
            marker=cell2mat(marker);
        end
      
end
function marker=RemoveMarker(Marker_ori,removemarker)

marker=Marker_ori;
for i=1:size(marker,1)
    marker(i,~ismember(marker(i,~isnan(marker(i,:))),removemarker))=NaN;
end


end
 
 
 



