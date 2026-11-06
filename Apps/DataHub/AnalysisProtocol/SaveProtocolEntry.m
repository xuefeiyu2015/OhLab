function OutputFiles = SaveProtocolEntry(Data, ProtocolMarker, Entry)
%SAVEPROTOCOLENTRY Write OutputData.(ProtocolMarker) = Entry into the
%per-neuron result files Results/[Monkey]_[Date]_N#_C#.mat, one file per
%channel in Data.SpikeChannel (batch mode passes a channel range, e.g.
%4:46). Other fields already in each file are preserved.
%
% Results folder: lab convention, found from the current folder
% ('DataAnalysis' on Mac, 'DataHub' on PC) + 'Results'.
%
% Inputs:
%   Data           : DataHub Data struct (FileName, SpikeChannel)
%   ProtocolMarker : field name in OutputData, e.g. 'OneDR_Performance'
%   Entry          : struct to store
% Output:
%   OutputFiles    : {nChannel x 1} full paths written
%
% Xuefei Yu

%% Output directory
Workingdirectory = pwd;
OperationSystem = computer;
if strcmp(OperationSystem(1:3), "PCW")
    MarkerFolder = 'DataHub';
else
    MarkerFolder = 'DataAnalysis';
end
Flag = strfind(Workingdirectory, MarkerFolder);
BasicDirectory = Workingdirectory(1:Flag(1)+length(MarkerFolder));
OutPath = strcat(BasicDirectory, 'Results', BasicDirectory(end));
if ~exist(OutPath, 'dir')
    mkdir(OutPath);
end

%% File names: [Monkey]_[Date]_N#_C#.mat, one per channel
FileName = Data.FileName(1:end-4);
num = regexp(FileName, '\d+', 'match');
RecordDateOriginal = num{1};
NeuronNum = num{2};
if isfield(Data, 'SpikeChannel') && ~isempty(Data.SpikeChannel)
    ChannelList = Data.SpikeChannel(:)';
else
    ChannelList = 1;
end
RemainFiles = erase(FileName, RecordDateOriginal);
Delimiter = find(isstrprop(RemainFiles, 'upper'));
MonkeyName = RemainFiles(Delimiter(1):Delimiter(2)-1);
OutputFiles = arrayfun(@(c) fullfile(OutPath, sprintf('%s_%s_N%s_C%d.mat', ...
    MonkeyName, RecordDateOriginal, NeuronNum, c)), ChannelList, 'UniformOutput', false)';

%% Merge into each channel's file
for k = 1:numel(OutputFiles)
    OutputData = struct();
    if exist(OutputFiles{k}, 'file')
        Loaded = load(OutputFiles{k}, 'OutputData');
        if isfield(Loaded, 'OutputData')
            OutputData = Loaded.OutputData;
        end
    end
    OutputData.(ProtocolMarker) = Entry;
    save(OutputFiles{k}, 'OutputData');
end
end
