function h = xline_special(y,varargin);
%
% h = xline_special(y);
%
% Plots a horizontall line on the current figure at point Y. Returns
% the handle h for that line.
%
% (Renamed from xline to xline_special so it no longer shadows MATLAB's
%  builtin xline, which draws a VERTICAL line and has a different signature.)
%
% h = xline_special(y,'color')
% Plots a line specified by the character string 'color' (default
% is black.
%
% h = xline_special(y,'PropertyName',PropertyValue)
% is equivalent to this:
% h = xline_special(y);
% set(h,'PropertyName',PropertyValue)
%
% 2012-apr-30: Revised to automatically transpose column matrices.
%
% last modified 2012-may-16
% dbtm

if nargin==2
    c = varargin{1};
    if (~ischar(c)|(length(c)~=1))
       help xline_special;
       error;
    end    
else
    c = 'k';    
end

if size(y,1)~=1
    y=y';
end

yy = [y; y];

x = get(gca,'XLim');
xx = repmat(x',[1 size(yy,2)]);
h = plot(xx,yy,c);

if nargin>=3
    npairs = nargin-1;
    if mod(npairs,2)~=0
        error('Number of properties and values don''t match!');
    end
    while(~isempty(varargin))
        property = varargin{1};
        value = varargin{2};
        if ~ischar(property)
            help xline_special;
            error;
        end
    set(h,property,value);
    varargin(1:2) = [];
    end
end









