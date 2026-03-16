function big(varargin)
%
% USAGE:
% figure(999);big;
% This function called without input arguments resizes the current
% figure. The original default size is:
% [LEFT BOTTOM WIDTH HEIGHT]
% [ 40    40    1200   750 ]
%
% You can change the default resizing behavior by including input arguments:
%
% big('reset')
% big('r')
% If the character string 'reset' or 'r' is passed as an input
% argument, the default size is reset to match the current figure.
%
% big(FIG);
% If a scalar is passed as input argument, the default size is set
% to match figure(FIG).
%
% big([LEFT BOTTOM WIDTH HEIGHT])
% A 4-element vector can also be passed to reset the default size explicitly.
% 
% last modified 2010-dec-17
% dbtm

persistent pos;

if isempty(pos)
    pos = [40 40 1200 750];
end

if nargin>0
    if ischar(varargin{1});
        % get current figure
        pos = get(gcf,'Position');
    elseif isnumeric(varargin{1})
        if length(varargin{1})==1
            % get specified figure
            pos = get(varargin{1},'Position');
        elseif length(varargin{1})==4
            pos = [varargin{1}];
        end
    end
elseif nargin==0
    set(gcf,'Position',pos);
end
