function OptoWaveformPlot(Data, StartTrial, EndTrial,ShowFigureFlag,OutputFlag);
% Load basic task information
ECode;

TaskType=Data.TaskType;
DataPath= Data.Path;
FileName=Data.FileName;
TaskCode=Data.TaskCode;

Selection=StartTrial:EndTrial;%Select the trials according to the pannel

useReISI = 0;
PlotExample = 1;
      
SegEvent = [ON_STIM];
SegEventEnd = [OFF_STIM];
SegEventStr={'StimOn'};

PlotStartOffset=[100];%in ms
PlotStopOffset=[100];%in ms

%Sequence information from the library
ParaLib=Data.ParaLib;
if ismember('stimOn',keys(ParaLib))
   stimOnFlag = ParaLib('stimOn');
   TrialStimOn = unique(stimOnFlag(stimOnFlag(:,2) == 1,1));
   stimOn = [Selection',0*ones(length(Selection),1)];
   stimOn(ismember(stimOn(:,1),TrialStimOn),2) = 1;
  
   stimOn= stimOn(:,2);
%   stimOn=stimOn(Selection,:);
   
end
if ismember('stimPower',keys(ParaLib))
   stimPower_ori = ParaLib('stimPower');
   stimPower_ori=stimPower_ori(:,2);
   
  
   
end

if ismember('preStimDur,StimDur,PostStimDur',keys(ParaLib))
   stimdur_info = ParaLib('preStimDur,StimDur,PostStimDur');
   stimDur=stimdur_info(:,3);
   
   stimDur=stimDur(Selection,:);
   
end


%Manuel Modification of the power by mistake
if ismember(FileName,{'Robin061423VIDEOVIEW1.rst','Robin062223VIDEOVIEW1.rst','Robin071323VIDEOVIEW1.rst','Robin072023VIDEOVIEW1.rst','Robin072723VIDEOVIEW1.rst'})
    stimPower = RescuePower(FileName,stimPower_ori);
    
else
    stimPower = stimPower_ori;
end

 stimPower=stimPower(Selection,:);

%Load event data
EventChannel=Data.EventChannel(Selection,:);
EventTimeChannel=Data.EventTimeChannel(Selection,:);




%Load the raw waveform
FileName=Data.FileName;
%Setup output directory
Workingdirectory=pwd;


%Set up output path according to the system
OperationSystem = computer;%Get the operation system information

if strcmp(OperationSystem(1:3),"PCW")  
    %PC
     MarkerFolder='DataHub';
elseif strcmp(OperationSystem(1:3),"MAC")  
    %Mac
     MarkerFolder='DataAnalysis';
 end
Flag=strfind(Workingdirectory,MarkerFolder);
BasicDirectory=Workingdirectory(1:Flag+length(MarkerFolder));


%Load the raw waveform 
%OutputFolerName='ForagingTask';
OutPath=strcat(BasicDirectory,'Results',BasicDirectory(end));%,OutputFolerName);
if ~exist(OutPath)
    disp('Please export the waveform file first');
    return

    %mkdir(OutPath)
end
cd(OutPath);


%Get File Information
num = regexp(FileName, '\d+', 'match');
RecordDateOriginal=cell2mat(num(1));

NeuronNum=cell2mat(num(2));



%To get the monkey name
RemainFiles=erase(FileName,RecordDateOriginal);
Delimiter=find(isstrprop(RemainFiles,'upper')==1);
MonkeyName=RemainFiles(Delimiter(1):Delimiter(2)-1);


%% Load the spike channel
SpikeChannelWhole=Data.SpikeChannel;
SpikeTimeWhole=Data.SpikeTimeData(:,:,Selection);
ChannelIDWhole=Data.ChannelID;
ClusterIDWhole=Data.ClusterID;

MultiChannel=0;
if length(SpikeChannelWhole)>1
    MultiChannel=1;
    
end

for spk=1:length(SpikeChannelWhole)
    SpikeChannel=SpikeChannelWhole(spk);
    
        SpikeTime=squeeze(SpikeTimeWhole(spk,:,:))';

        ChannelID=ChannelIDWhole(spk);
        ClusterID=ClusterIDWhole(spk);


        WaveFileName=sprintf('%s_%s_N%s_C%s.mat',MonkeyName,string(RecordDateOriginal),string(NeuronNum),string(SpikeChannel));
        if exist(WaveFileName)

            load(WaveFileName);
        else
            disp('No waveform file. Export waveform first!');
            return
    
        end

        Waveform = OutputData;
        try
            WaveformRaw = Waveform.Waveform;

        catch
             disp('No waveform field. Export waveform to the outputfile first!');
             

        end

       WaveformData = WaveformRaw.DataStamp;

       WaveformTime = WaveformData('TimePoint');
       Waveform = WaveformData('WaveForm');
       


       %{'TimePoint'}    {'WaveForm'}
       keyboard




end


%Data



keyboard 
end