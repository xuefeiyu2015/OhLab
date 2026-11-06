function h = yline(x,varargin);
%
% h = yline(x);
%
% Plots a vertical line on the current figure at point X. Returns
% the handle h for that line.
% 
% h = yline(x,'color')
% Plots a line specified by the character string 'color' (default
% is black.
%
% h = yline(x,'PropertyName',PropertyValue)
% is equivalent to this:
% h = yline(x);
% set(h,'PropertyName',PropertyValue)
%
% 2012-apr-30: Revised to automatically transpose column matrices.
%
% 
%
% last modified 2012-apr-30
% dbtm

if nargin==2
    c = varargin{1};
    if (~ischar(c)|(length(c)~=1))
       help yline;
       error;
    end    
else
    c = 'k';    
end

if size(x,1)~=1
    x=x';
end

xx = [x; x];

y = get(gca,'YLim');
yy = repmat(y',[1 size(xx,2)]);
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
            help yline;
            error;
        end
    set(h,property,value);
    varargin(1:2) = [];
    end
end









