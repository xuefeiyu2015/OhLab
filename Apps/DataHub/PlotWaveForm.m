%Plot Waveform in a new window YXF 11082019
function PlotWaveForm(WaveformData,WaveformName,SpikeChannel,ShowFigureFlag,OutputFlag)

WaveformName=WaveformName(2:end-3);
FileName=WaveformData.RawFileName;
%Order=WaveformData.Order;
%Reload Waveform Data
WaveformPath=WaveformData.Path;
DataType=WaveformPath(end-3:end);
if contains(DataType,'mat')
    %Load from Plexon
    WaveformData=WaveformPlexon(WaveformPath);
else
    %load from Blip original file
     WaveformData=LoadWaveform(WaveformPath);
    
end
figtitlestr{1}=WaveformName;
 fig{1}=PrepareFigure(1,'w',[50,500, 400,400],'Name',figtitlestr{1});

 
%Plot the waveform
            SpikeType_All=length(unique(WaveformData.Type))-1;
            if SpikeChannel>3 & SpikeChannel<=3+SpikeType_All+1
                    SpikeType=SpikeChannel-3;
                    OriginalSpikeChannel = SpikeChannel;
                elseif SpikeChannel<=3 
                    SpikeType=min(SpikeChannel,3);
                    OriginalSpikeChannel = SpikeChannel;
                    
                elseif SpikeChannel==3+size(WaveformData.Wform,1)+1
                    %Only online spikes
                     SpikeType=1:min(SpikeChannel,3);
                     OriginalSpikeChannel = 1:SpikeType;
                elseif SpikeChannel==3+size(WaveformData.Wform,1)+2
                    %All spikes;
                    SpikeType=1:SpikeType_All;
                    OriginalSpikeChannel = 1:SpikeType_All;
            end
           %{     
            if isempty(SpikeChannel)
                SpikeType=1:SpikeType_All;
            else
               SpikeType=SpikeChannel; 
            end
            %}
           
            for i=1:length(SpikeType)
            %T0=[0: unique( WaveformData.SR( WaveformData.Type==i-1)):(unique( WaveformData.DataPointNum(WaveformData.Type==i-1))-1)*unique(WaveformData.SR(WaveformData.Type==i-1))];
            T0=0:1/unique( WaveformData.SR( WaveformData.Type==SpikeType(i)-1)):1/unique( WaveformData.SR( WaveformData.Type==SpikeType(i)-1))*size(WaveformData.Wform(WaveformData.Type==SpikeType(i)-1,:),2);
            T0=T0(2:end)*10^3;
            %{
            RawTime=WaveformData.WTime(WaveformData.Type==i-1,:);
           AlignTime=WaveformData.WTime2Threshold(WaveformData.Type==i-1,:);
           WaveForm=WaveformData.Wform(WaveformData.Type==i-1,:);
           AlignTimeToOnset=AlignTime-RawTime;
           %}
            if ShowFigureFlag
              %  figtitlestr{i}=WaveformName;
              %  fig{i}=PrepareFigure(1+i-1,'w',[50,500, 400,400],'Name',figtitlestr{i});
          %    subplot(ceil(length(SpikeType)/ceil(sqrt(length(SpikeType)))),ceil(sqrt(length(SpikeType))),i);

            
                 % plot(T0,WaveformData.Wform(WaveformData.Type==i-1,:)','color',[0.3+(i-1)/2,0.8-(i-1)/2,0.1])
                 Waveform_ori=WaveformData.Wform(WaveformData.Type==SpikeType(i)-1,:);
                 %Remove extremem values 
                Index_tmp=sum(abs(Waveform_ori)>500,2);
                 %Remove waveform
                 Waveform_screened=Waveform_ori;
              %   Waveform_screened(Index_tmp>0,:)=[];
                 
                plot(T0,WaveformData.Wform(WaveformData.Type==SpikeType(i)-1,:)','color','#B7BEC9')
              % plot(T0,Waveform_screened(1,:)','color','#B7BEC9')
                TitleString=sprintf("Neuron %d Wave Form",SpikeType(i));
                %title(TitleString);
                xlim([0,3]);
                hold on
                if isfield(WaveformData,'MedWave')
                     plot(T0,WaveformData.MedWave(SpikeType(i),:),'-','LineWidth',3,'Color','k');
                else
                plot(T0,median(WaveformData.Wform(WaveformData.Type==SpikeType(i)-1,:),1,'omitnan'),'-','LineWidth',3,'Color','k');
                end
                box off
                axis off
            end
            %For output
            TimePointOutput{SpikeType(i)}=T0;
            AverageWaveform{SpikeType(i)}=nanmean(WaveformData.Wform(WaveformData.Type==SpikeType(i)-1,:));
            RawWaveform{SpikeType(i)}=WaveformData.Wform(WaveformData.Type==SpikeType(i)-1,:);
            
            
            end
            
 %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Output to
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%files%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

%%%%%%%%%%%Export data%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%OutputFolerName='ForagingTask';
OutPath=strcat(BasicDirectory,'Results',BasicDirectory(end));%,OutputFolerName);
if ~exist(OutPath)
    mkdir(OutPath)
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

for i=1:length(OriginalSpikeChannel)
    ChannelNum=OriginalSpikeChannel(i);

%Output file name
OutputFileName=sprintf('%s_%s_N%s_C%s.mat',MonkeyName,string(RecordDateOriginal),string(NeuronNum),string(ChannelNum));
OutputFigureName=sprintf('%s_%s_N%s_C%s',MonkeyName,string(RecordDateOriginal),string(NeuronNum),string(ChannelNum));




%Load the old file if exist
if exist(OutputFileName)

    load(OutputFileName);
    
end


%%%%%%%%%%%%%%%%%%%%OutputData For current analysis
OutputData_New=[];


%%%%%%%%%Save overall scene response characteristics%%%%%%%%%%%%%%%%%%%%%%%
%Organize the data
NameStr={'TimePoint','AverageWaveForm','WaveForm'};

DataLib={TimePointOutput{SpikeType(i)},AverageWaveform{SpikeType(i)}, RawWaveform{SpikeType(i)}};


DataStamp=containers.Map(NameStr,DataLib);

OutputData_New.Waveform.DataStamp=DataStamp;


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
    
    for j=1:TotalFigure
        OutputFigureName=strcat(OutputFigureName,figtitlestr{j},'.jpg');
        if ~isempty(fig{j})
        saveas(fig{j},OutputFigureName);
        end
        
    end

    disp('Figures have been exported to the neuron folder');

end
OutputData.Waveform=OutputData_New.Waveform; 





if OutputFlag
    cd(OutPath);
    save(OutputFileName,'OutputData');
    disp('Data Exported');
else
    disp('No data export');
end 

end
            
end