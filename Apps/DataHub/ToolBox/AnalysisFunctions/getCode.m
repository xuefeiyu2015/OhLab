function n = getCode(codes,cc);
%
% N = getCode(CODES,CC);
%
% Utility function used in converting relatively raw Cortex or Plexon
% datafiles into Hanuman format database.  Input argument CODES must be a
% vector of codes from a single trial, and CC must be a vector specifying
% the range of expected values of the code you're searching for.  Returns
% N, the member of CC that is found in CODES.  Note that N is the event
% code itself, not the timestamp of the event code.  Returns NaN if no
% element of CC occurs in CODES.  Issues a warning if more than one element
% occurs, or if a single element occurs more than once.
%
% Example of usage:  Suppose you're using a Cortex program with 8
% conditions, and you want to know what condition was used on trial t.
% Assuming your timing file drops the condition number as an event code in
% the datafile, you can do it like this:
%
%   data = ctx2mat(fname);
%   codes = data{t}{3};     
%   possible_conditions = [1:8];
%   cond = getCode(codes, possible_conditions);
%
% see also: getTs.m 
%
% last modified 2007-mar-07
% dbtm 

i = find(ismember(codes,cc));
if isempty(i)
    n = NaN;
elseif length(i)==1
    n = codes(i);
elseif length(i)>1
    warning(['One code ' num2str(cc) 'expected, more than one found.  ' ...
        'Returning only the first code.']); 
    n = codes(i(1));
end
