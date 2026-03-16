function ts = getTs(codes,times,c,varargin);
%
% TS = getTs(CODES,TIMES,C);
%
% Utility function used in converting relatively raw Cortex or Plexon
% datafiles into Hanuman format database.  Input argument CODES must be a
% vector of all codes dropped on a single trial, TIMES must be the
% corresponding timestamps,and C must be the code for the timestamp you're
% searching for.  Returns TS, timestamp of code C.  Returns NaN if no
% element of C occurs in CODES.
%
% Example of usage:  To retrieve the time of stimulus onset:
%   data = ctx2mat(fname);
%   times = data{t}{2};     
%   codes = data{t}{3};     
%   stim_on = 13;
%   cond = getTs(codes,times,stim_on);
%
% Optional input arguments: Use these if you expect more than one instance
% of code C per trial:
%
% TS_TRAIN = getTs(CODES,TIMES,C,'train');
% TS = getTs(CODES,TIMES,C,'first');
%
% 'train' option:
% If you expect code C to be dropped multiple times in the course of a
% single trial (for instance, if C is a spike, or if you use codes to mark
% the refresh rate), the function retrieves a vector TS if you pass the
% character string 'train' as a fourth optional input argument.  Examples:
%   spikes = getTs(codes,times,1,'train');
%   movie_frame_advances = getTs(codes,times,77,'train');
%
% 'first' option:
% If you expect more than one instance of C per trial but you're only
% interested in the first one, use this syntax:
%   reward_ts = getTs(codes,times,15,'first');
%
% see also: getCode.m 
%
% last modified 2007-mar-07
% dbtm 

if length(varargin)>0
    if strcmp(varargin{1},'train')
        train_flag=1;
        first_flag=0;
    elseif strcmp(varargin{1},'first')
        first_flag=1;
        train_flag=0;
    else
        error('unknown input argument.  Try "train" or "first" instead.')
    end
else
    train_flag=0;
    first_flag=0;
end

i = find(ismember(codes,c));
if isempty(i)
    ts = NaN;
elseif length(i)==1
    ts = times(i);
elseif length(i)>1
    if train_flag==1
        ts = times(i);              
    else
        ts = times(i(1));
        if first_flag==0            
            warning(['One code ' num2str(c) 'expected, more than one found.  ' ...
            'Returning only the first code.']); 
        end
    end
end

