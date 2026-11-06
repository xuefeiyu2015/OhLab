function Data = BuildObjLocBatchData(BatchFileDir, ResultDir, BatchFileName)
% BuildObjLocBatchData  Build the DataHub batch struct from a batch .txt file
% (no GUI / waitbar). Same logic as RunBatchObjLoc.
%
% Inputs
%   BatchFileDir  — folder with batch .txt files
%   ResultDir     — folder with per-neuron result .mat files
%   BatchFileName — e.g. 'FEF_BothRA_StableValue.txt'
% Output: struct Data with BatchFileName, ResultFilePath, RecordDate,
%         FileName, ChannelNumber, TrialInterval
%
% Xuefei Yu, 2026-11-05

B = importdata(fullfile(BatchFileDir, BatchFileName));
FileNames = B.textdata(~cellfun('isempty', B.textdata));
nF = numel(FileNames);
ResultFilePath = cell(nF,1);
RecordDate     = cell(nF,1);
for i = 1:nF
    f   = FileNames{i};
    num = regexp(f, '\d+', 'match');
    rem = erase(f, num{1});
    d   = find(isstrprop(rem, 'upper'));
    mk  = rem(d(1):d(2)-1);
    RecordDate{i} = datestr(datenum(num{1}, 'mmddyy'));
    chans = B.data(i,1):B.data(i,2);
    ResultFilePath{i} = arrayfun(@(j) fullfile(ResultDir, sprintf('%s_%s_N%s_C%d.mat', mk, num{1}, num{2}, j)), ...
                                 chans, 'UniformOutput', false);
end
Data.BatchFileName  = strrep(BatchFileName, '.txt', '');
Data.ResultFilePath = ResultFilePath;
Data.RecordDate     = RecordDate;
Data.FileName       = FileNames;
Data.ChannelNumber  = B.data(:,1:2);
Data.TrialInterval  = B.data(:,3:4);
end
