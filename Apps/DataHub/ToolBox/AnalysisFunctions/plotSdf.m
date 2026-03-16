function h = plotSdf(sdf,win,varargin);

% H = plotSdf(SDF,[PRE POST],[STYLE]);
%
% Default STYLE = 'k-'
%
% New function for plotting sdfs.
% Will add new features as I think of them and time permits.
%
% 2007-jul-06 Bug repair: 
% The function now handles the case where max(sdf)==0
%
% last modified 2007-jul-06
% dbtm

pre = win(1);
post = win(2);

if length(varargin)>=1
    style = varargin{1};
else
    style = 'k-';
end
if length(varargin)>=2
    ymax = varargin{2};
else
    if max(sdf)>0
        ymax = max(sdf)*1.25;
    else
        ymax = 1;
    end
end


samples = length(sdf);
window = post-pre;
binwidth = window/(samples-1);

time = [pre:binwidth:post];
h = plot(time,sdf,style,'LineWidth',2);    
axis([pre post 0 ymax])
hold on;
    
