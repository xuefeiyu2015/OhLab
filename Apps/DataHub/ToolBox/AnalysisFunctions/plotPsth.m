function h = plotPsth(psth,win,varargin)

% h = plotPsth(psth,[pre post])
%
% To specify the Y-axis range, use optional input argument:
% h = plotPsth(psth,[pre post], Ymax)
%
% New function for plotting psths.
% Will add new features as I think of them and time permits.
%
% last modified 041209
% dbtm

if ~isempty(varargin)
    ymax = varargin{1};
else
    ymax = max(psth)*1.25;
end

pre = win(1);
post = win(2);

% determine binwidth
samples = length(psth);
window = post+(-pre);
%window = post-pre;
binwidth = window/(samples-1);

% plot
time = [pre:binwidth:post];
h = bar(time,psth,1);    
set(h,'FaceColor','k')
axis([pre post 0 ymax])
hold on;
line([0 0],[0 ymax],'Color','k');
