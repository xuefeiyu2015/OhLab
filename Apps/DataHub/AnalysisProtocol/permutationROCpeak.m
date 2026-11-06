function p = permutationROCpeak(TrialIndex,psth,time,peakTime,peakAUC,window,NROC,Npermute)
% Permutation test for peak AUC around a specified time window.
% TrialIndex : vector of labels (1 vs 0)
% psth       : trials x time firing rates
% time       : 1 x T vector of time points
% peakTime   : time (ms) of observed peak AUC
% peakAUC    : observed mean AUC within the window
% window     : half-window size (ms) for averaging (e.g., 50 => [-50,+50])
% NROC       : number of thresholds for rocN (default 100)
% Npermute   : number of permutations (default 1000)

if nargin < 8 || isempty(Npermute), Npermute = 1000; end
if nargin < 7 || isempty(NROC),     NROC     = 100;  end

TrialIndex = TrialIndex(:);
time = time(:)';

winMask = time >= (peakTime - window) & time <= (peakTime + window);
if ~any(winMask) || isnan(peakAUC)
    p = NaN;
    return;
end

permPeakAUC = nan(Npermute,1);

for k = 1:Npermute
    % shuffle labels
    permLabel = TrialIndex(randperm(length(TrialIndex)));

    % time-resolved AUC for this shuffle
    auc_t = nan(1, size(psth,2));
    for t = 1:size(psth,2)
        x = psth(permLabel==1, t);
        y = psth(permLabel==0, t);
        auc_t(t) = rocN(x, y, NROC, 0);
    end

    % average AUC within the window around the (true) peak time
    permPeakAUC(k) = mean(auc_t(winMask), 'omitnan');
end

% Two-tailed p based on distance from 0.5
obsDiff  = abs(peakAUC - 0.5);
permDiff = abs(permPeakAUC - 0.5);
p = (sum(permDiff >= obsDiff) + 1) / (Npermute + 1);

end
