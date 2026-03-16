function BatchFileOut=CheckBatchFile(BatchFilePath,BatchFileName,RawPath,PackPath,RepackFlag,WaveFormFlag);
%Function to check whether all files in the batch file has been packed 
%Xuefei Yu 06032020
%Load bacth file

BatchFilePath_Whole=[BatchFilePath,BatchFileName];
BatchFile_temp=importdata(BatchFilePath_Whole);

ChannelNumber=[];
FileName=[];
FileType=[];
RecordDateFile=[];
ResultFilePath=[];
WaveformFileFlag=[];
WaveformFile=[];
WaveformPath=[];
TrialInterval=[];
TaskType=[];
TaskCode=[];
%Check whether it's a nested file
for i=1:length(BatchFile_temp)
    %Check whether it's another txt file
    currFile=BatchFile_temp(i);
    
   % if contains(currFile,'txt')
   try 
        BatchFile=importdata([BatchFilePath currFile{1}]);%If another txt file
        FileTypeStr{i}=currFile{1};
   catch
       BatchFile= BatchFile_temp;
        FileTypeStr{i}=BatchFileName;
   end
    
     ChannelNumber=[ChannelNumber;BatchFile.data(:,1:2)];
     TrialInterval=[TrialInterval;BatchFile.data(:,3:4)];
     %Remove spaces 
     temp=BatchFile.textdata;
     temp=temp(~cellfun('isempty',temp));
     FileName=[FileName;temp];
     FileType=[FileType;i*ones(length(ChannelNumber),1)];
     
    
end



%Output path
Flag=strfind(BatchFilePath,'Results');
ResultDirectory=BatchFilePath(1:Flag+length('Results'));


 Progressbar=waitbar(0,'Checking Batchfile....');   
%Check whether all files are packed and reporduce the raw file path
for i=1:length(FileName)
   
    
    FileCurr=FileName{i};
    %Get the date string
    num = regexp(FileCurr, '\d+', 'match');
    RecordDate=cell2mat(num(1));
    %Get File Information
     NeuronNum=cell2mat(num(2));
     SpikeChannel=ChannelNumber(i,:);


   %To get the monkey name
   RemainFiles=erase(FileCurr,RecordDate);
   Delimiter=find(isstrprop(RemainFiles,'upper')==1);
   MonkeyName=RemainFiles(Delimiter(1):Delimiter(2)-1);
   
   %File Path
  
    RawFilePath{i}=strcat(RawPath,RecordDate,RawPath(end),[FileCurr,'.rst']);
    RawFilePathWithoutFile{i}=strcat(RawPath,RecordDate,RawPath(end));
    PackedFilePath{i}=strcat(PackPath,[FileCurr,'.mat']);
    if ~contains(FileCurr,'Robin022924VIDEOVIEW1')
       
    if ~exist(PackedFilePath{i}) || RepackFlag
        disp('Re-packing data... Please wait!')
       %Repack data
              %Present a waitbar to indicate the progress
               ShowStr=sprintf('Packing file: %s ; %d remaining...',FileCurr,length(FileName)-i);
                waitbar(i/length(FileName), Progressbar,ShowStr); 
                if i==1
                   [GoodData,AllData]=LoadAndExportData(RawFilePath{i},PackPath,[FileCurr,'.mat']); 
                    TaskType=GoodData.TaskType;
                    TaskCode=GoodData.TaskCode;
                else
                    LoadAndExportData(RawFilePath{i},PackPath,[FileCurr,'.mat']); 
                    
                end
    else
        if i==1
            cd(PackPath);
           load([FileName{1},'.mat']);
           TaskType=GoodData.TaskType;
           TaskCode=GoodData.TaskCode;

        end
   
    end
    else
        cd(PackPath);
           load([FileName{end},'.mat']);
           TaskType=GoodData.TaskType;
           TaskCode=GoodData.TaskCode;

    end
    
    %Save output path
    RecordDateFile{i}=datestr(datenum(RecordDate,'mmddyy'));
    numindex=1;
    
    for j=SpikeChannel(1):SpikeChannel(2)
        OutputFileName=sprintf('%s_%s_N%s_C%s.mat',MonkeyName,string(RecordDate),string(NeuronNum),string(j));
    
        ResultFilePath{i}{numindex}=strcat(ResultDirectory,OutputFileName);
        numindex=numindex+1;
    end
    
    %Search for waveform file
     if WaveFormFlag
 
     %Automatically search for the waveform file
      RawFolder=[RawPath,RecordDate,RawPath(end)];
      WaveformFileAuto=[RecordDate 'SPK' NeuronNum '.wv'];
      %Also from the plaxon data
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
      Flag=strfind(RawPath,MarkerFolder);
      PlexonWaveformPath=strcat(RawPath(1:Flag+length(MarkerFolder)),PlexonDataFolder,SysDelimeter,PlexonParsedDataFolder,SysDelimeter,RecordDate,SysDelimeter);
      PlexonFileName=strcat(MonkeyName,RecordDate,'N',num2str(NeuronNum),'T',num2str(TaskCode),'wv','.mat');
      if isfile([RawFolder,WaveformFileAuto])
          WaveformFile{i}=WaveformFileAuto;
          WaveformFileFlag(i)=1;
           WaveformPath{i}=[RawFolder,WaveformFileAuto];
      
      elseif isfile([PlexonWaveformPath,PlexonFileName])
           WaveformFile{i}=PlexonFileName;
          WaveformFileFlag(i)=1;
           WaveformPath{i}=[PlexonWaveformPath,PlexonFileName];
        %  WaveformOrder{i}=
      else
          WaveformFile{i}=[];
          WaveformFileFlag(i)=0;
           WaveformPath{i}=[];

      end
                    
     end
    
  
    
end

 
close(Progressbar);
%msgbox('Batch file packing checking completed!');

%Get the file task type from the first(sample) file

%load([FileName{1},'.mat']);

 cd(PackPath);
BatchFileOut.PackedFilePath=PackedFilePath;
BatchFileOut.RawFilePath=RawFilePath;
BatchFileOut.ChannelNumber=ChannelNumber;
BatchFileOut.FileName=FileName;
BatchFileOut.TaskType=TaskType;
BatchFileOut.TaskCode=TaskCode;
BatchFileOut.FileNumber=length(FileName);
BatchFileOut.RawFilePathWithoutFile=RawFilePathWithoutFile;
BatchFileOut.RecordDate=RecordDateFile;
BatchFileOut.ResultFilePath=ResultFilePath;
BatchFileOut.TrialInterval=TrialInterval;

BatchFileOut.FileType=FileType;
BatchFileOut.FileTypeStr=FileTypeStr;

BatchFileOut.WaveFormFile=WaveformFile;
BatchFileOut.WaveFormFileFlag=WaveformFileFlag;
BatchFileOut.WaveFormPath=WaveformPath;

 



end