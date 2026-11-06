% CheckMissingKeys.m
% Report which keys are missing from Robin060722FRACSAC8 FreeViewNeural DataStamp
% compared to a neuron that loaded successfully.

ResultDir = '/Users/xuefeiyu/Documents/XuefeiFile/WorkRelated/Program_Matlab_Local/OHLab/DataAnalysis/Results/';

% The problematic file
ProblemFile = fullfile(ResultDir, 'Robin_060722_N8_C2.mat');

% Find a reference file with a complete FreeViewNeural DataStamp
% Use the first SC neuron that loaded OK (Robin042722FRACSAC1 or similar)
% We'll search for any SC result file that has FreeViewNeural
BatchFileDir = '/Users/xuefeiyu/Documents/XuefeiFile/WorkRelated/Program_Matlab_Local/OHLab/DataAnalysis/Results/BatchFile/';
BatchData = importdata(fullfile(BatchFileDir, 'SC_BothRA_StableValue.txt'));
FileNames = BatchData.textdata;
FileNames = FileNames(~cellfun('isempty', FileNames));
ChannelNumbers = BatchData.data(:, 1:2);

% Find a reference file with complete FreeViewNeural
RefFile = '';
for i = 1:length(FileNames)
    FileCurr = FileNames{i};
    num = regexp(FileCurr, '\d+', 'match');
    RecordDateStr = cell2mat(num(1));
    NeuronNum = cell2mat(num(2));
    RemainFiles = erase(FileCurr, RecordDateStr);
    Delimiter = find(isstrprop(RemainFiles, 'upper') == 1);
    MonkeyName = RemainFiles(Delimiter(1):Delimiter(2)-1);
    j = ChannelNumbers(i, 1);
    fname = fullfile(ResultDir, sprintf('%s_%s_N%s_C%d.mat', MonkeyName, RecordDateStr, NeuronNum, j));
    if ~exist(fname, 'file'), continue; end
    tmp = load(fname, 'OutputData');
    if isfield(tmp, 'OutputData') && isfield(tmp.OutputData, 'FreeViewNeural')
        ds = tmp.OutputData.FreeViewNeural.DataStamp;
        if isKey(ds, 'PSTH_Stim_Mean_norm')
            RefFile = fname;
            break;
        end
    end
end

if isempty(RefFile)
    disp('Could not find a reference file with complete FreeViewNeural DataStamp.');
    return;
end

fprintf('Reference file: %s\n', RefFile);
fprintf('Problem file:   %s\n\n', ProblemFile);

% Load both DataStamps
refData  = load(RefFile,     'OutputData');
probData = load(ProblemFile, 'OutputData');

refKeys  = keys(refData.OutputData.FreeViewNeural.DataStamp);
probKeys = keys(probData.OutputData.FreeViewNeural.DataStamp);

missing = setdiff(refKeys, probKeys);
extra   = setdiff(probKeys, refKeys);

fprintf('Keys in reference but MISSING from problem file (%d):\n', length(missing));
for k = 1:length(missing)
    fprintf('  - %s\n', missing{k});
end

if ~isempty(extra)
    fprintf('\nKeys in problem file but not in reference (%d):\n', length(extra));
    for k = 1:length(extra)
        fprintf('  + %s\n', extra{k});
    end
end
