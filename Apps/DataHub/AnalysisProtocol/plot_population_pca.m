function [scores_struct, explained_variance] = plot_population_pca( ...
    psth_object_low, psth_object_high, psth_location_low, psth_location_high, ...
    n_components, zscore_flag)
% Plot population PCA trajectories for four conditions.
% Inputs:
%   psth_object_low    : [n_neurons x n_time]
%   psth_object_high   : [n_neurons x n_time]
%   psth_location_low  : [n_neurons x n_time]
%   psth_location_high : [n_neurons x n_time]
%   n_components       : number of PCs to retain (default 3)
%   zscore_flag        : if true, z-score each neuron across all conditions (default true)
%
% Outputs:
%   scores_struct      : struct with fields object_low/high, location_low/high (n_time x n_components)
%   explained_variance : variance explained ratio for retained components
%
% Note: rows = neurons, columns = time points.

if nargin < 5 || isempty(n_components), n_components = 3; end
if nargin < 6 || isempty(zscore_flag), zscore_flag = true; end

% basic checks
conds = {psth_object_low, psth_object_high, psth_location_low, psth_location_high};
nt = cellfun(@(x) size(x,2), conds);
if numel(unique(nt)) ~= 1
    error('All PSTH inputs must have the same number of time points.');
end

% concatenate conditions along time for joint PCA
X = [psth_object_low, psth_object_high, psth_location_low, psth_location_high]; % [neurons x (4*time)]

% z-score per neuron if requested
if zscore_flag
    mu = mean(X, 2, 'omitnan');
    sd = std(X, 0, 2, 'omitnan');
    sd(sd == 0) = 1;
    X = (X - mu) ./ sd;
end

% observations = time points, features = neurons
[coeff, score, latent] = pca(X', 'Centered', false); %#ok<ASGLU>

% variance explained
latent = latent(:);
explained_variance = latent ./ sum(latent);

% keep requested components
n_components = min(n_components, size(score, 2));
score = score(:, 1:n_components);
explained_variance = explained_variance(1:n_components);

% split scores back to conditions (each block nt(1) long)
tlen = nt(1);
idx1 = 1:tlen;
idx2 = tlen + (1:tlen);
idx3 = 2*tlen + (1:tlen);
idx4 = 3*tlen + (1:tlen);

scores_struct.object_low    = score(idx1, :);
scores_struct.object_high   = score(idx2, :);
scores_struct.location_low  = score(idx3, :);
scores_struct.location_high = score(idx4, :);

% plotting
colors = struct('object_low', [0 0 1], ...       % blue
                'object_high',[1 0 0], ...       % red
                'location_low',[0 0.6 0], ...    % green
                'location_high',[0.6 0 0.6]);    % purple

figure;
if n_components >= 3
    hold on;
    plot3(scores_struct.object_low(:,1), scores_struct.object_low(:,2), scores_struct.object_low(:,3), 'Color', colors.object_low, 'LineWidth', 2);
    plot3(scores_struct.object_high(:,1), scores_struct.object_high(:,2), scores_struct.object_high(:,3), 'Color', colors.object_high, 'LineWidth', 2);
    plot3(scores_struct.location_low(:,1), scores_struct.location_low(:,2), scores_struct.location_low(:,3), 'Color', colors.location_low, 'LineWidth', 2);
    plot3(scores_struct.location_high(:,1), scores_struct.location_high(:,2), scores_struct.location_high(:,3), 'Color', colors.location_high, 'LineWidth', 2);
    xlabel('PC1'); ylabel('PC2'); zlabel('PC3');
    grid on; axis vis3d;
else
    hold on;
    plot(scores_struct.object_low(:,1), scores_struct.object_low(:,2), 'Color', colors.object_low, 'LineWidth', 2);
    plot(scores_struct.object_high(:,1), scores_struct.object_high(:,2), 'Color', colors.object_high, 'LineWidth', 2);
    plot(scores_struct.location_low(:,1), scores_struct.location_low(:,2), 'Color', colors.location_low, 'LineWidth', 2);
    plot(scores_struct.location_high(:,1), scores_struct.location_high(:,2), 'Color', colors.location_high, 'LineWidth', 2);
    xlabel('PC1'); ylabel('PC2'); grid on;
end

legend({'Object Low','Object High','Location Low','Location High'}, 'Location', 'best');
title('Population PCA Trajectory');

end
