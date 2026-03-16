function [spike_count fr trial_length] = countSpikes(dat,channel,event,win);

% [SPIKE_COUNT FR TRIAL_LENGTH = countSpikes(dat,channel,event,win);
%
% Returns the trial-by-trial spike-count occurring within a specified
% integration window.  DAT Must be a hanuman data structure.  CHANNEL
% indicates the spike-channel.  EVENT must be a key (column) indicating an
% event timestamp.  WIN must be a two-element vector specifying the
% beginning and the end of the integration window.
%
% SPIKE_COUNT is a vector with one element for each trial in DAT, and
% contains the raw spike-count.  FR is the firing rate in Hz.  TRIAL_LENGTH
% is the integration window in milliseconds.  TRIAL_LENGTH will be a
% shorter time period than the input argument WIN on trials where the trail
% ends before the integration window closes (or trials that begin after the
% integration window opens.
%
% last modified 2006-oct-19
% dbtm

% compute raw spike-count

eval(unpackHeader(dat));

%dat = alignTimebase(dat,event);
spikes = dat.s(:,:,channel);
H = histc(spikes',win(1):win(2))';
spike_count = sum(H')';

% compute trial length
%[Ntrials junk] = size(dat.c);
%win_start = repmat(win(1),Ntrials,1);
%win_stop = repmat(win(2),Ntrials,1);
%start_trial = max([win_start dat.c(:,T_START)]')';
%stop_trial = min([win_stop dat.c(:,T_STOP)]')';
%trial_length = stop_trial - start_trial;

% firing rate in Hz
trial_length = win(2)-win(1);
fr = 1000*spike_count./trial_length;
