function fold_idx = AssignFolds(n, k)
% AssignFolds  Randomly assign n trials to k folds (balanced, stratified).
%
% Distributes trials as evenly as possible across k folds using a
% round-robin assignment on a random permutation, so fold sizes differ
% by at most 1.
%
% Inputs:
%   n  - number of trials
%   k  - number of folds
% Output:
%   fold_idx [n x 1] - fold label (1..k) for each trial

fold_idx = zeros(n, 1);
perm = randperm(n);
for i = 1:n
    fold_idx(perm(i)) = mod(i - 1, k) + 1;
end
end
