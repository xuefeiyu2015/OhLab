function [psth sem]  = getPsth(dat,channel,win,varargin);

% VER0.4 NOTE: NOW ONLY PROVIDES PSTH OUTPUT.
%
%function [psth sem stddev ntt]  = getPsth(dat,channel,win,varargin);
%
% [PSTH] = getPsth(DAT,CHANNEL,WINDOW);
%
% Hanuman version 0.4 function.
%
% Returns a peri-stimulus time histogram (PSTH) and corresponding
% time-base (TIME).  DAT must be a hanuman database with a spike data field
% (DAT.S).  CHANNEL is a scalar that specifies which spike channel is to be
% processed (i.e., DAT.S{CHANNEL} is selected).  WINDOW must be a two
% element vector specifying the beginning and end of the time-window.  The
% time-base in DAT should be aligned before this function is called (use
% alignTimebase.m)
%
% Optional input parameter BINWIDTH can be used to specify... wait for
% it... the binwidth.  Default value is 10 msec):
%   PSTH = ras2psth(DAT,CHANNEL,WINDOW,BINWIDTH);  
%
% Extra output parameters:
%   [PSTH SEM STDDEV NTT] = ras2psth(DAT,CHANNEL,WINDOW);


% NTRIALS returns the number of trials for each element of TIME.  If some
% trials terminate before the integration window closes, NTRIALS will get
% smaller as more trials drop out.  STDEV and SEM are the standard
% deviation and standard error of the mean around PSTHas a function of
% time. 
%
% 2013-dec-21 Cleaned up handling of SEM output argument.
% 
% last modified 2013-dec-21
% dbtm

pre = win(1);
post = win(2);

if ~isempty(varargin)
    binwidth = varargin{1};
else
    binwidth = 10;
end

% pad time-window.
time = pre-binwidth:binwidth:post+binwidth;

% % get the number of trials at each timestep
% ntt = get_n_trials_time(dat,win,binwidth);
% blanks = find(ntt==0);
% ntt(blanks) = NaN;
% 
% construct matrix of 1's and 0's
spikes = dat.s(:,:,channel);
H = histc(spikes',time)';

[Ntrials junk] = size(H);
if Ntrials>1
    sumspikes = sum(H);
elseif Ntrials==1
    sumspikes = H;
elseif Ntrials==0
    sumspikes = nan(1,length(time));
end   

% crop time-window
sumspikes([1 end]) = [];

%psth = sumspikes./ntt;
psth = sumspikes./Ntrials;
psth = psth*(1000/binwidth);

foo = [];
stddev = nanstd(H)*(1000/binwidth);
stddev([1 end]) = [];
sem = stddev./sqrt(Ntrials);
