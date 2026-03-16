function h = rasBitmap(dat,spk,win);

% h = rasBitmap(dat,spk,win);
%
% Utility function for making a raster plot as a bitmap.
%
% Intened to replace the vector-based utility function plotRasters.
%
% last modified 2013-jan-02
% dbtm



%load /data/mcmahond/fishRugs/mas/sieglinde20110525b.mat
%load /data/mcmahond/dogsButterflies/mas/sdbsig003a.mat

unpack;
%spk = 1
%win = [-100 500]
%dat = selectTrials(dat,STIM,3);
spikes = dat.s(:,:,spk);

SPIKEHEIGHT = 1;

pre = win(1);
post = win(2);
[sr sc] = size(spikes);
y = [1:sr]';
Y = repmat(y,1,sc);

mask = find(spikes<win(1) | spikes>win(2));
spikes(mask) = NaN;

N = sr*sc;
SpikeVec(1:N) = spikes(1:N)';
TrialVec(1:N) = Y(1:N)';
crop = find(isnan(SpikeVec)); 
TrialVec(crop) = [];
SpikeVec(crop) = [];
SpikeVec = SpikeVec-win(1);


ras = zeros(sr,length([win(1):win(2)]));
index = two2one(SpikeVec,TrialVec,ras);
ras(index) = 1;

%figure;
imagesc(abs(ras-1));colormap(gray)


