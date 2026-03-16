%Demon of 2-D kernal density
%Xuefei Yu
% Generate some normally distributed data
x = randn(50,1);
y = randn(50,1);
% Estimate a continuous pdf from the discrete data
[pdfx xi]= ksdensity(x);
[pdfy yi]= ksdensity(y);
% Create 2-d grid of coordinates and function values, suitable for 3-d plotting

[xxi,yyi]     = meshgrid(xi,yi);
[pdfxx,pdfyy] = meshgrid(pdfx,pdfy);
% Calculate combined pdf, under assumption of independence
pdfxy = pdfxx.*pdfyy; 
% Plot the results
figure
subplot(2,2,1)
mesh(xxi,yyi,pdfxy)
subplot(2,2,2)
contourf(xxi,yyi,pdfxy);
set(gca,'XLim',[min(xi) max(xi)])
set(gca,'YLim',[min(yi) max(yi)])
