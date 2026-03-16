function smoothx = smoothfun(x,sigma);

%smoothx = smoothfun(x,sigma);
%
% Smooths function X by convolution with Gaussian kernel of width SIGMA.
% SIGMA is in units of points in X.
%
% created 2018-jul-18
% dbtm

% pad x w/mirror-image of itself (to reduce edge-artifacts).
mirror = fliplr(x);
padded = [mirror x mirror];

% convolve x with kernel
kwidth = -sigma*3:3*sigma;
%kernel = normpdf(kwidth,0,sigma);
kernel = openGaussian(kwidth,0,sigma);
s = conv(padded,kernel);
s(1:length(mirror)) = [];
s(length(s)-length(mirror)+1:length(s)) = [];
s(1:floor(length(kernel)/2)) = [];
s(length(x)+1:length(s)) = [];

smoothx = s;

