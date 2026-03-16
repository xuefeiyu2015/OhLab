function z = nanzscore(v)
%
%  z = nanzscore(v)
% 
% Fixes the mean/nanmean bug in zscore. Thanks, Mathworks!
%
% last modified 2014-may-04
% dbtm

z = [v-nanmean(v)]./nanstd(v);
