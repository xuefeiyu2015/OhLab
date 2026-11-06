% RunBatchObjLoc.m
% Runner script for Batch_ObjLocValue_Organize — FEF and SC separately.
% Builds the Data struct directly from batch txt files to avoid
% CheckBatchFile's waitbar (which requires a display).
%
% Run headless from terminal:
%   /Applications/MATLAB_R2025b.app/bin/matlab -nodisplay -nosplash -batch \
%     "cd('/Users/xuefeiyu/Documents/XuefeiFile/WorkRelated/Program_Matlab_Local/OHLab/DataAnalysis/DataHub_YXF/Apps/DataHub'); RunBatchObjLoc"
%
% Xuefei Yu — 2026-03-23

%% Paths
DataHubPath = '/Users/xuefeiyu/Documents/XuefeiFile/WorkRelated/Program_Matlab_Local/OHLab/DataAnalysis/DataHub_YXF/Apps/DataHub';
addpath(genpath(DataHubPath));

BatchFileDir = '/Users/xuefeiyu/Documents/XuefeiFile/WorkRelated/Program_Matlab_Local/OHLab/DataAnalysis/Results/BatchFile/';
ResultDir    = '/Users/xuefeiyu/Documents/XuefeiFile/WorkRelated/Program_Matlab_Local/OHLab/DataAnalysis/Results/';

%% Batch files and function to run (pre-set these variables to override)
if ~exist('BatchFiles', 'var'),     BatchFiles = {'FEF_BothRA_StableValue.txt', 'SC_BothRA_StableValue.txt'}; end
if ~exist('BatchFunc', 'var'),      BatchFunc = @Batch_ObjLocValue_Organize; end   % or @Batch_ObjLocValue_SingleNeuron
if ~exist('ShowFigureFlag', 'var'), ShowFigureFlag = 0; end
if ~exist('OutputFlag', 'var'),     OutputFlag = 1; end

for b = 1:length(BatchFiles)
    BatchFileName = BatchFiles{b};
    disp('==========================================================');
    disp(['Running batch: ', BatchFileName]);
    disp('==========================================================');

    %% Parse batch txt file
    % Format: FileName  startChan  endChan  startTrial  endTrial (tab-separated)
    BatchData     = importdata(fullfile(BatchFileDir, BatchFileName));
    FileNames     = BatchData.textdata;
    FileNames     = FileNames(~cellfun('isempty', FileNames));
    ChannelNumbers = BatchData.data(:, 1:2);
    TrialIntervals = BatchData.data(:, 3:4);

    %% Build ResultFilePath, RecordDate
    ResultFilePath = cell(length(FileNames), 1);
    RecordDate     = cell(length(FileNames), 1);

    for i = 1:length(FileNames)
        FileCurr = FileNames{i};

        % Parse date string and neuron number from filename
        num           = regexp(FileCurr, '\d+', 'match');
        RecordDateStr = cell2mat(num(1));
        NeuronNum     = cell2mat(num(2));

        % Parse monkey name (text between first and second uppercase letter)
        RemainFiles = erase(FileCurr, RecordDateStr);
        Delimiter   = find(isstrprop(RemainFiles, 'upper') == 1);
        MonkeyName  = RemainFiles(Delimiter(1):Delimiter(2)-1);

        RecordDate{i} = datestr(datenum(RecordDateStr, 'mmddyy'));

        % Build per-channel result file paths
        SpikeChannel = ChannelNumbers(i, :);
        numindex = 1;
        for j = SpikeChannel(1):SpikeChannel(2)
            OutputFileName = sprintf('%s_%s_N%s_C%d.mat', ...
                MonkeyName, RecordDateStr, NeuronNum, j);
            ResultFilePath{i}{numindex} = fullfile(ResultDir, OutputFileName);
            numindex = numindex + 1;
        end
    end

    %% Build Data struct
    Data                = struct();
    Data.BatchFileName  = strrep(BatchFileName, '.txt', '');
    Data.ResultFilePath = ResultFilePath;
    Data.RecordDate     = RecordDate;
    Data.FileName       = FileNames;
    Data.ChannelNumber  = ChannelNumbers;
    Data.TrialInterval  = TrialIntervals;
    if exist('SingleNeuronParams', 'var')      % optional overrides for Batch_ObjLocValue_SingleNeuron
        Data.SingleNeuronParams = SingleNeuronParams;
    end

    %% Run batch (all files, no figures, save output)
    StartFile       = 1;
    EndFile         = length(FileNames);
    BatchFunc(Data, StartFile, EndFile, ShowFigureFlag, OutputFlag);

    disp(['Completed: ', BatchFileName]);
end

disp('All batches done.');
