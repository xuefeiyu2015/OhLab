function h = plotRasters(spikes,win);

% H = plotRasters(DAT,[PRE POST])
%
% Utility function for making rasterplots.
%
% DAT must be a NxM matrix in which each row N is a separate trial.
% Time-stamps marking the occurrence of each spike are contained in columns
% SPK_BEGIN to SPK_END in DAT (most of these elements will be NaNs).
%
% WIN specifies the time window over which spikes are to be counted (in
% milliseconds).
% 
% last modified 2006-sept-9
% dbtm

%global SPK_BEGIN  SPK_END

SPIKEHEIGHT = 1;

pre = win(1);
post = win(2);

%spikes = dat(:,SPK_BEGIN:SPK_END);
[sr sc] = size(spikes);
y = -[1:sr]';
Y = repmat(y,1,sc);

N = sr*sc;
SpikeVec(1:N) = spikes(1:N);
TrialVec(1:N) = Y(1:N);
blanks = find(isnan(SpikeVec));
SpikeVec(blanks) = [];
TrialVec(blanks) = [];

s = SpikeVec;
y0 = TrialVec;
y1 = TrialVec+SPIKEHEIGHT;
h = line([s;s],[y0;y1],'Color','k');
axis([pre post -sr-1 0]);
% axis off;
hold on;
line([0 0],[-sr-1 0],'Color','k');
