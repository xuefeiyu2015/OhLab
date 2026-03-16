% Load fisheriris dataset
load fisheriris

% Define the number of dimensions to reduce to
k = 2;

% Define the target class labels
target_labels = {'setosa', 'versicolor'};

% Find the indices of the target class labels in the dataset
target_indices = ismember(species, target_labels);

% Define the colors for the target class labels
colors = lines(numel(target_labels));

% Convert logical indices to integer indices
target_indices = find(target_indices);

% Map integer indices to colors
target_colors = colors(target_indices, :);

% Perform targeted dimension reduction using t-SNE
X_target = tsne(meas(target_indices, :), 'NumDimensions', k);

% Plot the reduced data points
scatter(X_target(:, 1), X_target(:, 2), 15, target_colors, 'filled');
