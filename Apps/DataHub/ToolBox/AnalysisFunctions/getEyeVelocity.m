function [v vt] = getEyeVelocity(x,y,t);
%
% [v vt] = getEyeVelocity(x,y,t);
%
% Computes the first derivative of eye position over time dR/dt.
%
% INPUT ARGUMENTS
% x and y       horizontal and vertical eye position vectors
% t             time base for x and y
%
% OUTPUT ARGUMENTS
% V  : eye velocity
% VT : time base of V
%
% Note that V will have one fewer elements than X and Y.
%
% last modified 2018-aug-28:
% got rid of superfluous sample period as an input argument.
% dbtm

x=gaussSmooth(x,3);
y=gaussSmooth(y,3);

SamplePeriod = t(2)-t(1);
deltax = x(2:end)-x(1:end-1);
deltay = y(2:end)-y(1:end-1);
r = sqrt([deltax.^2] + [deltay.^2]);
v = (1000/SamplePeriod)*r; % convert to degrees per second.
vt = t(1:end-1)+(0.5*SamplePeriod);
