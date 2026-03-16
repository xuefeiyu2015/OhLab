function smoothed = gaussSmooth(x,sigma)
%
% smoothed = gaussSmooth(x,sigma)
%
% Applies smoothing to X via convolution with gaussian kernel of width
% SIGMA.
%
% last modified 2013-dec-20
% dbtm

mirror = fliplr(x);
padded = [mirror x mirror];
kwidth = -sigma*3:3*sigma;
kernel = normpdf(kwidth,0,sigma);
s = conv(padded,kernel);
s(1:length(mirror)) = [];
s(length(s)-length(mirror)+1:length(s)) = [];
s(1:floor(length(kernel)/2)) = [];
s(length(x)+1:length(s)) = [];
smoothed = s;