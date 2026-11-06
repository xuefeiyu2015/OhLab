function sortIdx = ComputeAUCSortOrder(aucScalar, sortMode)
% Purpose : row order for AUC heatmaps, based on a per-neuron scalar AUC.
% Inputs  : aucScalar — N-vector of AUC values (e.g. object-task TgtAUC_RF)
%           sortMode  — 'signed' (descending AUC: high-prefer top, low-prefer
%                       bottom) or 'abs' (descending |AUC-0.5|)
% Outputs : sortIdx   — N-by-1 index vector; NaN AUCs placed last
    aucScalar = aucScalar(:);
    switch lower(sortMode)
        case 'signed', key = aucScalar;
        case 'abs',    key = abs(aucScalar - 0.5);
        otherwise, error('ComputeAUCSortOrder: unknown sortMode "%s"', sortMode);
    end
    [~, sortIdx] = sort(key, 'descend', 'MissingPlacement', 'last');
end
