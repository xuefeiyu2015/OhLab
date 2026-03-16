function [peaks starttimes stoptimes] = findSaccades(v,varargin);
%
% [peaks start stop] = findSaccades(v,threshold,noise);
%
% Detects saccades in eye velocity vector V. Saccades are defined
% as events with peak velocity > 200 deg/sec. Start and stop are
% the times at wich eye velocity exceeds "noise" level of 50
% deg/sec.
%
% Returns the indices corresponding to saccade PEAKS, START times,
% and STOP times. Thus peak saccade amplitude is:
% V(PEAKS)
%
% For a vector T that is the time base corresponding to V:
% saccade duration = T(STOP)-T(START)
%
% Ignores eye movements occurring at the beginning and end of V for
% which start and stop times are not defined.
% 
% last modified 2011-mar-24
% dbtm

if isempty(varargin)
    thresh = 200;
    noise = 50;
elseif length(varargin)==2
    thresh = varargin{1};
    noise = varargin{2};
else
    error;help findSaccades
end
    
peaks = intersect(find(v>thresh),localmax(v));
peaks = peaks';
for p = 1:length(peaks)
    try
        starttimes(p) = max(find(v(1:peaks(p))<noise));
    catch
        starttimes(p) = NaN;
    end
    try
        stoptimes(p) = min(find(v(peaks(p):end)<noise))+peaks(p)-1;
    catch
        stoptimes(p) = NaN;
    end
end

if ~isempty(peaks)
    duds = union(find(isnan(starttimes)),find(isnan(stoptimes)));
    peaks(duds) = [];
    starttimes(duds) = [];
    stoptimes(duds) = [];
else
    starttimes = [];
    stoptimes = [];
end


