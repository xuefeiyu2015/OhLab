function [sdf sem stddev ntt]  = getSdf(dat,channel,win,varargin);


%
% [SDF TIME] = getSdf(DAT,CHANNEL,WINDOW);
%
% Returns a smoothed spike density function (SDF) and corresponding
% time-base (TIME).  DAT must be a hanuman database with a spike data field
% (DAT.S).  CHANNEL is a scalar that specifies which spike channel is to be
% processed (i.e., DAT.S{CHANNEL} is selected).  WINDOW must be a two
% element vector specifying the beginning and end of the time-window.  The
% time-base in DAT should be aligned before this function is called (use
% alignTimebase.m)
%
% Optional input parameter KERNEL can be used to specify sigma of gaussian
% kernel (default value is 10 msec):
%   SDF = getSdf(DAT,CHANNEL,WINDOW,KERNEL);
%
% Extra output parameters:
%   [SDF SEM stddev ntt] = getSdf(DAT,CHANNEL,WINDOW);
% NTRIALS returns the number of trials for each element of TIME.  If some
% trials terminate before the integration window closes, NTRIALS will get
% smaller as more trials drop out.  STDEV and SEM are the standard
% deviation and standard error of the mean around SDF as a function of
% time. 
%
% 2013-dec-21 Cleaned up handling of SEM output argument.
%
% last modified 2012-dec-21
% dbtm

pre = win(1);
post = win(2);

if ~isempty(varargin)
    sigma = varargin{1};
else
    sigma = 10;
end

% if isempty(dat)
%     warning('Input argument DAT is empty!!!  Returning NaN outputs.');
%     foo = nan(1,post-pre);
%     sdf = foo;
%     sem = foo;
%     stddev = foo;
%     ntt = foo;
%     return
% end

% get psth
binwidth = 1;
[psth err] = getPsth(dat,channel,win,1);

% pad psth w/mirror-image of itself (to reduce edge-artifacts).
mirror = fliplr(psth);
padded = [mirror psth mirror];

% convolve psth with kernel
kwidth = -sigma*3:3*sigma;
%kernel = normpdf(kwidth,0,sigma);
kernel = openGaussian(kwidth,0,sigma);
s = conv(padded,kernel);
s(1:length(mirror)) = [];
s(length(s)-length(mirror)+1:length(s)) = [];
s(1:floor(length(kernel)/2)) = [];
s(length(psth)+1:length(s)) = [];

sdf = s;


% APPLY EQUIVALENT TREATMENT TO SEM.
% pad psth w/mirror-image of itself (to reduce edge-artifacts).
mirror = fliplr(err);
padded = [mirror err mirror];

% convolve err with kernel
kwidth = -sigma*3:3*sigma;
%kernel = normpdf(kwidth,0,sigma);
kernel = openGaussian(kwidth,0,sigma);
s = conv(padded,kernel);
s(1:length(mirror)) = [];
s(length(s)-length(mirror)+1:length(s)) = [];
s(1:floor(length(kernel)/2)) = [];
s(length(err)+1:length(s)) = [];

sem = s;



% get the number of trials at each timestep
%ntt = get_n_trials_time(dat,win,binwidth);
%blanks = find(ntt==0);
%ntt(blanks) = NaN;
% 
% ntt = size(dat.s,1);
% 
% try
%     stddev = nanstd(sdf);
%     sem = stddev./sqrt(ntt);
% catch
%     ;
%     %    disp('mathworks sucks');
% end