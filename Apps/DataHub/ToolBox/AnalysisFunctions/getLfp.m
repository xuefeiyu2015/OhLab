function [lfp sem] = getLfp(dat,channel,win)

% [LFP SEM] = getLfp(DAT,CHANNEL,WIN)
%
% Returns a voltage trace LFP averaged across all trials for the field
% potential CHANNEL listed in DAT.  WIN specifies the time window around
% time zero.  WIN must be a 2 element vector = [PRE POST].
%
% last modified 2010-sept-15
% dbtm

t = dat.t;
a = dat.lfp(:,:,channel);

pre = find(t==win(1));
post = find(t==win(2));

time = t(pre:post); 
analog = a(:,pre:post);

lfp = nanmean(analog);

sd = nanstd(analog);
n = sum(~isnan(analog));
sem = sd./sqrt(n);


