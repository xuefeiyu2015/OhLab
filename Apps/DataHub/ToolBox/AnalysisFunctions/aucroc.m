function [auc dprime] = aucroc(noise,sig)
%
% [auc dprime] = aucroc(noise,sig)
%
% Computes area under the curve of Receiver Operator Characteristic (ROC).
%
% AUC       Area Under Cuve (between 0 and 1; 0.5 means no discrimination).
% d-prime   Discriminibility for an ideal observer.
%
% created 2018-jul-31
% dbtm

if ndims(noise)>2 | ndims(sig)>2 
    error('Input arguments must be 2D vectors.');
end
if size(noise,1)>1
    noise=noise';
end
if size(sig,1)>1
    sig=sig';
end
steps = unique(sort([sig noise]));
steps = [min(steps)-1 steps max(steps)+1];

for t=1:length(steps)
    hits(t) = length(find(sig>=steps(t)))/length(sig);
    
    fa(t) = length(find(noise>=steps(t)))/length(noise);
end
hits = fliplr(hits);
fa = fliplr(fa);
%auc = sum(hits/length(sig))/length(hits);
%dprime = [];
%threshold = steps(find(max(hits./fa)));

dh = [hits(2:end)+hits(1:end-1)]/2;
df = fa(2:end)-fa(1:end-1);
auc = sum(dh.*df);

zh = icdf('norm',hits,0,1);
zf = icdf('norm',fa,0,1);
dprime = [mean(sig)-mean(noise)]/sqrt(0.5*[std(sig)^2 + std(noise)^2]);

