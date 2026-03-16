function h = plotSdfErrorShading(t,sdf,sem,cc);
%
% h = plotSdfErrorShading(t,sdf,sem,cc);
%
% Utility function for plotting errorbars around SDF (or anything else) as
% a shaded region of specified color. draws patch object on the current
% figure.
%
% INPUT ARGUMENTS:
% T values along the x-axis
% SDF the mean of the signal
% SEM the errorbars around the signal.
% cc the color of the shaded patch.
%
% OUTPUT ARGUMENT:
% H handle of patch object drawn onto current figure.
% 
% last modified 2013-dec-21
% dbtm
%

tt = [t fliplr(t)];
yy = [sdf-sem fliplr(sdf+sem)];
h = patch(tt,yy,cc);
set(h,'EdgeColor',cc);

% % r = ones(11,1);
% % b = [1:-0.1:0]';
% % g=b;
% % pink = [r g b];
% % figure;
% % imagesc([1:10]');
% % colormap(pink);
% % 
% % blue = pink(:,[3 2 1]);
% % figure;
% % imagesc([1:10]');
% % colormap(blue);
% 
% pink = [1 0.8 0.8];
% blue = [0.8 0.8 1];
% 
% x = [0:0.01:10];
% y = sin(x);
% y2 = cos(x);
% figure(222);clf;
% hold on;
% 
% xx = [x fliplr(x)];
% yy = [y-0.1 fliplr(y)+0.1];
% h = patch(xx,yy,pink);
% set(h,'EdgeColor',pink)
% 
% plot(x,y,'r','LineWidth',2);
% 
% yy2 = [y2-0.1 fliplr(y2)+0.1];
% 
% h = patch(xx,yy2,blue);
% set(h,'EdgeColor',blue)
% 
% plot(x,y2,'b','LineWidth',2);
% 
