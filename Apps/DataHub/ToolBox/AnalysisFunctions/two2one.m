function index = two2one(x,y,m)
%
%  INDEX = two2one(ROWS,COLUMNS,MATRIX)
% 
% Converts 2D coordinates specified by paired vectors ROWS and
% COLUMNS referenencing MATRIX to 1D vector INDEX.
%
% last modified 2013-jan-02
% dbtm 

%x = [1:5];
%y = [2:2:10];
%m = zeros(10,5);

if length(x)~=length(y)
    error('x and y must be the same length');
end

if size(m,2)>10^6
    error('m has more than 10^6 columns!!!')
end

if size(x,1)==1
    x=x';
end
if size(y,1)==1
    y=y';
end

ind = [1:size(m,1)*size(m,2)]';
r = repmat([1:size(m,1)]' , [size(m,2) 1]);
cstep = 1/size(m,1);
c = ceil(cstep:cstep:size(m,2))';
lut = [ind r c];
lut(:,4) = [10^6]*r+c;
yx = [10^6]*y+x;

index = find(ismember(lut(:,4),yx));