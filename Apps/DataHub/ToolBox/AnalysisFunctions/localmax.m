function m = localmax(x,varargin)
%
% m = localmax(x)
%
% Returns a vector M containing the indices of all local maxima in
% the vector x.
% 
% m = localmax(x,'>=')
% Including the character string '>=' as a second input argument
% finds all points in v that are greater than or equal to the
% neighboring points, excluding the lowest points in v.
%
% last modified 2011-mar-24
% dbtm

current = x(2:end-1);
pre = x(1:end-2);
post = x(3:end);
if isempty(varargin)
    m = find(current>pre & current>post)+1;
elseif isequal(varargin{1},'>=')
    m = find(current>=pre & current>=post & current>min(x))+1;    
else
    help(localmax);
    error;
end

