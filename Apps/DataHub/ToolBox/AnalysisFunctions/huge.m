function huge(varargin)
%
% 
% 
% last modified 2012-dec-18
% dbtm

huge = [1 1 2375 1179];

if ~isempty(varargin)
    han=varargin{1};
else
    han = gcf;
end
set(han,'Position',huge);