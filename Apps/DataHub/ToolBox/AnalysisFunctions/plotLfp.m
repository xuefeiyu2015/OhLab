function h = plotLfp(lfp,win,varargin);

% H = plotLfp(LFP,[PRE POST],[YMIN YMAX],[STYLE]);
%
% Default STYLE = 'k-'
%
% Function for plotting LFPs.
%
% last modified 2006-sept-9
% dbtm

%disp('smoothed lfp');
%lfp = gaussSmooth(lfp);

pre = win(1);
post = win(2);

if length(varargin)>=1
    ymax = max(varargin{1});
    ymin = min(varargin{1});
else
    ymax = max(lfp)*1.25;
    ymin = min(lfp)*1.25;
end
if length(varargin)>=2
    style = varargin{2};
else
    style = 'k-';
end

samples = length(lfp);
window = post-pre;
binwidth = window/(samples-1);

time = [pre:binwidth:post];
h = plot(time,lfp,style,'LineWidth',2);    
axis([pre post ymin ymax])
hold on;
    
