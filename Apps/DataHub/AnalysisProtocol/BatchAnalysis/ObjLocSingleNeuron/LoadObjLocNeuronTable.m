function T = LoadObjLocNeuronTable(FilesName, FileNameEach, ChannelID, ExcludeOptoFiles, LoadTrials)
% LoadObjLocNeuronTable  Load per-neuron object (StableValue) and location
% (OneDR) value fields, plus opto (FreeViewNeural) fields, from result files.
%
% Inputs
%   FilesName        — cell{nFile} of cell{nChan} result .mat paths (Data.ResultFilePath)
%   FileNameEach     — cell{nFile} session names (Data.FileName)
%   ChannelID        — nFile-by-2 [startChan endChan] (Data.ChannelNumber)
%   ExcludeOptoFiles — cellstr of session names whose opto data are ignored
%   LoadTrials       — (optional, default false) also load RF single-trial rates
%
% Output: struct T, one entry per neuron (N), fields:
%   FileName{N}, Channel(N), HasStable(N), HasOneDR(N), HasOpto(N)
%   Obj.* / Loc.*  (task = StableValue / OneDR):
%       TgtAUC(N), pRanksum(N), pVis(N), VisDiff(N), FRMean(N),
%       nHighRF(N), nLowRF(N)  (trial counts, RF high/low value),
%       AUCPSTH(N,T), Time(1,T)
%       If LoadTrials: Trials(N) struct with .Rate (nRFtrial-by-T raw FR,
%       unsmoothed 10 ms bins, from 'psth_targeton_raw') and .IsHigh (logical);
%       TrialTime(1,T) (= Time_TgtOn)
%   Opto.Latency(N), Opto.SM(N), Opto.SM_p(N)
%   A missing task leaves NaN in that task's fields only.
%
% Xuefei Yu, 2026-11-05

if nargin < 5, LoadTrials = false; end
nN = sum(cellfun(@numel, FilesName));
T.FileName  = cell(nN,1);
T.Channel   = nan(nN,1);
T.HasStable = false(nN,1);
T.HasOneDR  = false(nN,1);
T.HasOpto   = false(nN,1);
T.Obj  = InitTask(nN, LoadTrials);
T.Loc  = InitTask(nN, LoadTrials);
T.Opto = struct('Latency', nan(nN,1), 'SM', nan(nN,1), 'SM_p', nan(nN,1));
aucObj = cell(nN,1);
aucLoc = cell(nN,1);

rawKeysObj = {'PSTH_TgtOn_Good_RF','PSTH_TgtOn_Bad_RF','PSTH_SacOn_Good_RF','PSTH_SacOn_Bad_RF', ...
              'PSTH_TgtOn_Good_nonRF','PSTH_TgtOn_Bad_nonRF','PSTH_SacOn_Good_nonRF','PSTH_SacOn_Bad_nonRF'};
rawKeysLoc = {'PSTH_Good_TgtOn_RF','PSTH_Bad_TgtOn_RF','PSTH_SacGood_RF','PSTH_SacBad_RF', ...
              'PSTH_Good_TgtOn_nonRF','PSTH_Bad_TgtOn_nonRF','PSTH_SacGood_nonRF','PSTH_SacBad_nonRF'};

n = 0;
for i = 1:numel(FilesName)
    chans = ChannelID(i,1):ChannelID(i,2);
    for j = 1:numel(FilesName{i})
        n = n + 1;
        T.FileName{n} = FileNameEach{i};
        T.Channel(n)  = chans(min(j, numel(chans)));

        S = load(FilesName{i}{j}, 'OutputData');
        OD = S.OutputData;

        if isfield(OD, 'StableValue')
            T.HasStable(n) = true;
            [T.Obj, aucObj{n}] = FillTask(T.Obj, n, OD.StableValue(1).DataStamp, rawKeysObj, LoadTrials);
        end
        if isfield(OD, 'OneDR')
            T.HasOneDR(n) = true;
            [T.Loc, aucLoc{n}] = FillTask(T.Loc, n, OD.OneDR.DataStamp, rawKeysLoc, LoadTrials);
        end
        if isfield(OD, 'FreeViewNeural') && ~ismember(FileNameEach{i}, ExcludeOptoFiles)
            M = OD.FreeViewNeural.DataStamp;
            T.HasOpto(n)      = true;
            T.Opto.Latency(n) = M('ResponseLatency');
            T.Opto.SM(n)      = M('StimModulation');
            T.Opto.SM_p(n)    = M('StimModulation_p');
        end
    end
end

T.Obj.AUCPSTH = StackRows(aucObj, numel(T.Obj.Time), 'StableValue');
T.Loc.AUCPSTH = StackRows(aucLoc, numel(T.Loc.Time), 'OneDR');
end

% ---------------------------------------------------------------------------
function S = InitTask(nN, loadTrials)
S = struct('TgtAUC', nan(nN,1), 'pRanksum', nan(nN,1), 'pVis', nan(nN,1), ...
           'VisDiff', nan(nN,1), 'FRMean', nan(nN,1), 'nHighRF', nan(nN,1), 'nLowRF', nan(nN,1), ...
           'AUCPSTH', [], 'Time', []);
if loadTrials
    S.Trials    = repmat(struct('Rate', [], 'IsHigh', false(0,1)), nN, 1);
    S.TrialTime = [];
end
end

function [S, aucRow] = FillTask(S, n, M, rawKeys, loadTrials)
% Copy one neuron's scalar fields from a DataStamp map into task struct S.
S.TgtAUC(n)   = M('TgtAUC_RF');
S.pRanksum(n) = M('p_AUC_Tgt_RF_ranksum');
S.pVis(n)     = M('p_Vis_RF');
S.VisDiff(n)  = M('VisDiff_RF');
S.FRMean(n)   = mean(cellfun(@(k) mean(M(k), 'omitnan'), rawKeys));
S.nHighRF(n)  = M('n_high_RF');
S.nLowRF(n)   = M('n_low_RF');
aucRow        = double(M('AUCPSTH_TgtOn_Value_RF'));
t             = M('Time_TgtOn');
if isempty(S.Time)
    S.Time = t(:)';
end
if loadTrials
    P  = double(M('psth_targeton_raw'));
    y  = M('target_index'); rf = M('RF_index');
    if size(P,1) ~= numel(y) || size(P,2) ~= numel(t)
        error('LoadObjLocNeuronTable:TrialMismatch', ...
            'Neuron %d: psth_targeton_raw is %dx%d, expected %dx%d.', n, size(P,1), size(P,2), numel(y), numel(t));
    end
    inRF = rf(:) == 1;
    S.Trials(n).Rate   = P(inRF,:);
    S.Trials(n).IsHigh = y(inRF) == 1;
    S.Trials(n).IsHigh = S.Trials(n).IsHigh(:);
    if isempty(S.TrialTime), S.TrialTime = t(:)'; end
end
end

function A = StackRows(rows, nT, taskName)
% Stack per-neuron AUC rows into N-by-T; missing neurons -> NaN row.
A = nan(numel(rows), nT);
for k = 1:numel(rows)
    if isempty(rows{k}), continue; end
    if numel(rows{k}) ~= nT
        error('LoadObjLocNeuronTable:TimeMismatch', ...
            '%s AUCPSTH for neuron %d has %d bins, expected %d.', taskName, k, numel(rows{k}), nT);
    end
    A(k,:) = rows{k}(:)';
end
end
