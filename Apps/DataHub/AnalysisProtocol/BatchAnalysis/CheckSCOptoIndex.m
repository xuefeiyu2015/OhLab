% CheckSCOptoIndex.m
% Report which SC neurons are missing optogenetics (FreeViewNeural) data.
% Xuefei Yu — 2026-03-23

BatchFileDir = '/Users/xuefeiyu/Documents/XuefeiFile/WorkRelated/Program_Matlab_Local/OHLab/DataAnalysis/Results/BatchFile/';
ResultDir    = '/Users/xuefeiyu/Documents/XuefeiFile/WorkRelated/Program_Matlab_Local/OHLab/DataAnalysis/Results/';
BatchFileName = 'SC_BothRA_StableValue.txt';

BatchData      = importdata(fullfile(BatchFileDir, BatchFileName));
FileNames      = BatchData.textdata;
FileNames      = FileNames(~cellfun('isempty', FileNames));
ChannelNumbers = BatchData.data(:, 1:2);

fprintf('\n%-35s  %-6s  %s\n', 'Neuron File', 'Chan', 'Opto Status');
fprintf('%s\n', repmat('-', 1, 70));

noOpto      = {};
incompleteOpto = {};

for i = 1:length(FileNames)
    FileCurr      = FileNames{i};
    num           = regexp(FileCurr, '\d+', 'match');
    RecordDateStr = cell2mat(num(1));
    NeuronNum     = cell2mat(num(2));
    RemainFiles   = erase(FileCurr, RecordDateStr);
    Delimiter     = find(isstrprop(RemainFiles, 'upper') == 1);
    MonkeyName    = RemainFiles(Delimiter(1):Delimiter(2)-1);
    SpikeChannel  = ChannelNumbers(i, :);

    for j = SpikeChannel(1):SpikeChannel(2)
        OutputFileName = sprintf('%s_%s_N%s_C%d.mat', MonkeyName, RecordDateStr, NeuronNum, j);
        FilePath = fullfile(ResultDir, OutputFileName);

        if ~exist(FilePath, 'file')
            fprintf('%-35s  C%-5d  FILE MISSING\n', FileCurr, j);
            noOpto{end+1} = OutputFileName;
            continue;
        end

        tmp = load(FilePath, 'OutputData');
        if ~isfield(tmp, 'OutputData')
            fprintf('%-35s  C%-5d  NO OutputData\n', FileCurr, j);
            noOpto{end+1} = OutputFileName;
            continue;
        end

        od = tmp.OutputData;

        if ~isfield(od, 'FreeViewNeural')
            fprintf('%-35s  C%-5d  NO FreeViewNeural\n', FileCurr, j);
            noOpto{end+1} = OutputFileName;
        else
            % Check key completeness
            ds = od.FreeViewNeural.DataStamp;
            if ~isKey(ds, 'PSTH_Stim_Mean_norm')
                fprintf('%-35s  C%-5d  INCOMPLETE DataStamp\n', FileCurr, j);
                incompleteOpto{end+1} = OutputFileName;
            end
        end
    end
end

fprintf('\n--- Summary ---\n');
fprintf('Missing FreeViewNeural entirely: %d neurons\n', length(noOpto));
fprintf('FreeViewNeural present but incomplete DataStamp: %d neurons\n', length(incompleteOpto));
