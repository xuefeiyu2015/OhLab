function rgb = redblue(varargin)
%
% RGB = redblue
%
% Function for making a red-blue color map in which magnitude is
% indicated by color saturation as follows:
% -1 pure blue
%  0 pure white
% +1 pure red
% 
% default = 21 rows. To change this (e.g., 101 rows) use optional
% input argument:
% RGB = redblue(NROWS);
% 
% last modified 2014-apr-28
% dbtm

if isempty(varargin)
    nrows = 21;
else
    nrows = varargin{1};
end

step = 1/floor(nrows/2)

if mod(nrows,2)==0
    g = [1:-step:step/2]';
else
    g = [1:-step:0]';
end
b = g;
r = ones(length(g),1);
red = [r g b];

blue = fliplr(red);
blue = flipud(blue);
if mod(nrows,2)==1
    blue(end,:) = [];
end

rgb = [blue;red];

