function [X_out, time_out] = PreprocessPSTH(X_in, time_in, varargin)
% PreprocessPSTH  Smooth and downsample PSTH data for decoding.
%
% No padding is applied — the caller should pass a wider time window than
% needed and trim afterward (e.g. pass [-600,400] then select [-500,300]).

% Inputs:
%   X_in     {N x 1} cell array — each cell is [trials x nTimeBins_in] for one neuron
%   time_in  [1 x nTimeBins_in]  time axis in ms
%
% Optional name-value:
%   'win_ms'     50      Boxcar: window width in ms; Gaussian: sigma in ms (default 50)
%   'kernel'     'boxcar'  Kernel type: 'boxcar' or 'gaussian' (default 'boxcar')
%   'ds_factor'    1     downsample factor (default 1)
%
% Outputs:
%   X_out    {N x 1} cell array — each cell is [trials x nTimeBins_out]
%   time_out [1 x nTimeBins_out]
%
% Example:
%   [X_obj_ds, time_ds] = PreprocessPSTH(psth_targeton_raw_stable_Others_RF_high, PSTH_Time_TgtOn_stable);
%   [X_obj_ds, time_ds] = PreprocessPSTH(X_in, time_in, 'kernel', 'gaussian', 'win_ms', 30);

%% Parse options
p = inputParser;
addParameter(p, 'win_ms',    50, @isnumeric);
addParameter(p, 'kernel', 'boxcar', @(x) ismember(x, {'boxcar','gaussian'}));
addParameter(p, 'ds_factor',  1, @isnumeric);
parse(p, varargin{:});
win_ms      = p.Results.win_ms;
kernel_type = p.Results.kernel;
ds_factor   = p.Results.ds_factor;

%% Validate
assert(iscell(X_in), 'X_in must be a cell array {N x 1}');
nTimeBins = numel(time_in);
assert(size(X_in{1}, 2) == nTimeBins, ...
    'Each cell in X_in must have size(cell,2) == numel(time_in)');
nNeurons = numel(X_in);

% Infer step size from time axis
dt_ms = round((time_in(2) - time_in(1)));  % ms

%% Build symmetric kernel (centered)
win_bins = round(win_ms / dt_ms);
if mod(win_bins, 2) == 0
    win_bins = win_bins + 1;  % ensure odd length for symmetry
end

switch kernel_type
    case 'boxcar'
        kernel = ones(1, win_bins) / win_bins;
    case 'gaussian'
        % win_ms is used directly as sigma; window spans ±3*sigma
        sd_bins = win_ms / dt_ms;
        half = ceil(3 * sd_bins);
        win_bins = 2 * half + 1;  % override win_bins for gaussian
        x = -half : half;
        kernel = exp(-0.5 * (x / sd_bins).^2);
        kernel = kernel / sum(kernel);  % normalize
end

fprintf('[PreprocessPSTH] %s=%dms (%d bins) | ds=%dx\n', ...
    kernel_type, win_ms, win_bins, ds_factor);

%% Smooth and downsample — loop over neurons
time_out = time_in(1 : ds_factor : end);
nTimeBins_out = numel(time_out);
X_out = cell(nNeurons, 1);

for n = 1:nNeurons
    slice = single(X_in{n});             % [nTrials x nTimeBins]

    % Convolve with 'same' — no padding, edges handled naturally
    % Since the input window is wider than the analysis window,
    % edge artifacts will be trimmed by the caller.
    
        slice_conv = conv2(slice, kernel, 'same');
    

    % Downsample
    X_out{n} = slice_conv(:, 1 : ds_factor : end);  % [nTrials x nTimeBins_out]
    
end

fprintf('[PreprocessPSTH] Output: %d neurons x %d bins (%.0f to %.0f ms, %dms step)\n', ...
    nNeurons, nTimeBins_out, ...
    time_out(1), time_out(end), ...
    round((time_out(2)-time_out(1))));

end
