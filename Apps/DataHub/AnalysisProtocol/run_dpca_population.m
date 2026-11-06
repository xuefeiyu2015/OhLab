function [scores_struct, explained, axisOverlap] = run_dpca_population(obj_low, obj_high, loc_low, loc_high, nComp)
% run_dpca_population  dPCA with correct condition ordering and separate trajectories.

if nargin < 5 || isempty(nComp), nComp = 10; end

%% 1) ensure toolbox
if ~isfolder('dpca-master')
    system('git clone https://github.com/elifesciences-publications/dpca-master.git');
end
addpath(genpath('dpca-master'));

%% 2) input check
conds = {obj_low,obj_high,loc_low,loc_high};
nt = cellfun(@(x) size(x,2), conds);
if numel(unique(nt)) ~= 1
    error('All inputs must have the same number of time points.');
end
N = size(obj_low,1); T = nt(1);

%% 3) tensor: Neuron x Source x Value x Time
X = zeros(N,2,2,T);
X(:,1,1,:) = obj_low;
X(:,1,2,:) = obj_high;
X(:,2,1,:) = loc_low;
X(:,2,2,:) = loc_high;

% z-score across all conditions/time per neuron
Xflat = reshape(X,N,[]);
mu = mean(Xflat,2,'omitnan');
sd = std(Xflat,0,2,'omitnan'); sd(sd==0)=1;
Xflat = (Xflat - mu)./sd;
X_centered = reshape(Xflat,size(X));

%% reorder to (N x T x S x V) so time is slowest in flattening blocks
X_tv = permute(X_centered,[1 4 2 3]); % N x T x S x V
X_lin = reshape(X_tv, N, T*4);        % order: for each time, S then V

%% 4) dPCA
combinedParams = {{3,[1 3],[2 3],[1 2 3]}, {1,[1 2]}, {2}}; % time, source, value
[W,V,whichMarg] = dpca(reshape(X_tv,[N 2 2 T]), nComp, 'combinedParams', combinedParams); %#ok<ASGLU>

%% explained variance (robust)
expl = dpca_explainedVariance(reshape(X_tv,[N 2 2 T]), W, V, 'combinedParams', combinedParams);
if isstruct(expl)
    if isfield(expl,'explainedVar') && isfield(expl.explainedVar,'total')
        explained = expl.explainedVar.total(:);
    elseif isfield(expl,'totalVarExplained')
        explained = expl.totalVarExplained(:);
    else
        explained = expl(:);
    end
else
    explained = expl(:);
end
explained = explained(1:min(nComp,numel(explained)));

%% 5) overlap (value axis)
valIdx = find(whichMarg==3,1,'first'); if isempty(valIdx), valIdx=1; end
axis_obj = mean(obj_high - obj_low,2); axis_obj = axis_obj./(norm(axis_obj)+eps);
axis_loc = mean(loc_high - loc_low,2); axis_loc = axis_loc./(norm(axis_loc)+eps);
cosSim = dot(axis_obj, axis_loc); cosSim = max(min(cosSim,1),-1);
axisOverlap.cosine = cosSim;
axisOverlap.angle_deg = acosd(cosSim);
axisOverlap.valueDecoder = W(:,valIdx);

%% 6) project and reshape back by condition (no zigzag)
proj = (W(:,1:min(3,nComp))' * X_lin)'; % (T*4) x pcs
proj = reshape(proj, T, 2, 2, []);       % time x source x value x pc (S fastest after time)
proj = permute(proj,[4 2 3 1]);          % pc x source x value x time

scores_struct.object_low    = squeeze(proj(:,1,1,:)).';
scores_struct.object_high   = squeeze(proj(:,1,2,:)).';
scores_struct.location_low  = squeeze(proj(:,2,1,:)).';
scores_struct.location_high = squeeze(proj(:,2,2,:)).';

%% 7) plot separate trajectories
colors = {[0 0 1],[1 0 0],[0 0.6 0],[0.6 0 0.6]};
if size(proj,1) >=3
    figure; hold on;
    plot3(scores_struct.object_low(:,1), scores_struct.object_low(:,2), scores_struct.object_low(:,3), 'Color', colors{1}, 'LineWidth', 2);
    plot3(scores_struct.object_high(:,1), scores_struct.object_high(:,2), scores_struct.object_high(:,3), 'Color', colors{2}, 'LineWidth', 2);
    plot3(scores_struct.location_low(:,1), scores_struct.location_low(:,2), scores_struct.location_low(:,3), 'Color', colors{3}, 'LineWidth', 2);
    plot3(scores_struct.location_high(:,1), scores_struct.location_high(:,2), scores_struct.location_high(:,3), 'Color', colors{4}, 'LineWidth', 2);
    xlabel('dPC1'); ylabel('dPC2'); zlabel('dPC3'); grid on; axis vis3d;
else
    figure; hold on;
    plot(scores_struct.object_low(:,1), scores_struct.object_low(:,2), 'Color', colors{1}, 'LineWidth', 2);
    plot(scores_struct.object_high(:,1), scores_struct.object_high(:,2), 'Color', colors{2}, 'LineWidth', 2);
    plot(scores_struct.location_low(:,1), scores_struct.location_low(:,2), 'Color', colors{3}, 'LineWidth', 2);
    plot(scores_struct.location_high(:,1), scores_struct.location_high(:,2), 'Color', colors{4}, 'LineWidth', 2);
    xlabel('dPC1'); ylabel('dPC2'); grid on;
end
legend({'Obj Low','Obj High','Loc Low','Loc High'}, 'Location','best');
title('dPCA population trajectories');

end
