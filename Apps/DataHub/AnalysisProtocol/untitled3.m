figure(1)
xd=rand(1,50);
yd=rand(1,50);

n=49;
xi=linspace(min(xd(:)),max(xd(:)),n);
yi=linspace(min(yd(:)),max(yd(:)),n);

xr=interp1(xi,1:numel(xi),xd,'nearest');
yr=interp1(yi,1:numel(yi),yd,'nearest');

z=accumarray([xr,yr],1,[n,n]);

figure(2)
surf(z)