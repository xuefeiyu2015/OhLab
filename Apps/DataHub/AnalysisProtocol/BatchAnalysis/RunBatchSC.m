% RunBatchSC.m — Run Batch_ObjLocValue_Organize for SC only.
% Called headless after FEF completed successfully.
% Xuefei Yu — 2026-03-23

addpath(genpath('/Users/xuefeiyu/Documents/XuefeiFile/WorkRelated/Program_Matlab_Local/OHLab/DataAnalysis/DataHub_YXF/Apps/DataHub'));

BatchFileDir = '/Users/xuefeiyu/Documents/XuefeiFile/WorkRelated/Program_Matlab_Local/OHLab/DataAnalysis/Results/BatchFile/';
ResultDir    = '/Users/xuefeiyu/Documents/XuefeiFile/WorkRelated/Program_Matlab_Local/OHLab/DataAnalysis/Results/';
BatchFileName = 'SC_BothRA_StableValue.txt';

disp(['Running batch: ', BatchFileName]);

BatchData      = importdata(fullfile(BatchFileDir, BatchFileName));
FileNames      = BatchData.textdata;
FileNames      = FileNames(~cellfun('isempty', FileNames));
ChannelNumbers = BatchData.data(:, 1:2);
TrialIntervals = BatchData.data(:, 3:4);

ResultFilePath = cell(length(FileNames), 1);
RecordDate     = cell(length(FileNames), 1);

for i = 1:length(FileNames)
    FileCurr      = FileNames{i};
    num           = regexp(FileCurr, '\d+', 'match');
    RecordDateStr = cell2mat(num(1));
    NeuronNum     = cell2mat(num(2));
    RemainFiles   = erase(FileCurr, RecordDateStr);
    Delimiter     = find(isstrprop(RemainFiles, 'upper') == 1);
    MonkeyName    = RemainFiles(Delimiter(1):Delimiter(2)-1);
    RecordDate{i} = datestr(datenum(RecordDateStr, 'mmddyy'));
    SpikeChannel  = ChannelNumbers(i, :);
    numindex = 1;
    for j = SpikeChannel(1):SpikeChannel(2)
        OutputFileName = sprintf('%s_%s_N%s_C%d.mat', MonkeyName, RecordDateStr, NeuronNum, j);
        ResultFilePath{i}{numindex} = fullfile(ResultDir, OutputFileName);
        numindex = numindex + 1;
    end
end

Data                = struct();
Data.BatchFileName  = strrep(BatchFileName, '.txt', '');
Data.ResultFilePath = ResultFilePath;
Data.RecordDate     = RecordDate;
Data.FileName       = FileNames;
Data.ChannelNumber  = ChannelNumbers;
Data.TrialInterval  = TrialIntervals;

Batch_ObjLocValue_Organize(Data, 1, length(FileNames), 0, 1);

disp('SC batch completed.');
