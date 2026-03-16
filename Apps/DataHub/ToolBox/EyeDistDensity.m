function [xxi,yyi,pdfxy]=EyeDistDensity(x,y);
%To measure the distribution of eye position on the screen
%Xuefei Yu 03232020
%Reshape data
x_reshape=[];
y_reshape=[];
for i=1:size(x,1)
    x_reshape=[x_reshape,x{i}];
    y_reshape=[y_reshape,y{i}];
    
    
    
end


% Estimate a continuous pdf from the discrete data
[N,Xbins,Ybins]=hist2d( x_reshape,y_reshape,50,'probability') ;

[xxi,yyi]= meshgrid(Xbins,Ybins);

%contourf(xxi,yyi,N);


[pdfx xi]= ksdensity(x_reshape);

[pdfy yi]= ksdensity(y_reshape);

% Create 2-d grid of coordinates and function values, suitable for 3-d plotting
[xxi,yyi]= meshgrid(xi,yi);
[pdfxx,pdfyy]= meshgrid(pdfx,pdfy);

% Calculate combined pdf, under assumption of independence
pdfxy =  pdfxx.*pdfyy;

%{
figure
subplot(2,2,1)
mesh(xxi,yyi,pdfxy)
subplot(2,2,2)
contourf(xxi,yyi,pdfxy);
set(gca,'XLim',[min(xi) max(xi)])
set(gca,'YLim',[min(yi) max(yi)])
%}


end