function v = sqdiag(r,varargin);
%
% % V = sqdiag(R);
% Returns a vector V containing the diagonal of a square
% matrix R. Useful for returning the diagonal of a correlation
% matrix, for instance, in which the diagonal represents odd
% vs. even split halves.
%
% V = sqdiag(R,FLAG);
% Optional input argument for specifying off-diagonal components.
% 
% Possible values of FLAG (only the first letter matters):
%    'diagonal' Returns the X=Y diagonal (default)
%    'top'   Returns only the X>Y values
%    'bottom' Returns only the X<Y values
%    'offdiagonal' Returns both 'top' and 'bottom'
%    'all' Returns all cells 
% 
% last modified 2014-apr-28
% dbtm

if isempty(varargin)
    flag = 'd';
else
    flag = lower(varargin{1}(1));
    if ~ismember(flag,'dtboa');
        help sqdiag;
        error;
    end
end

if size(r,1)~=size(r,2)
    error('R must be a square matrix');
end

N = size(r,1);
% remove the diagonal and redundant symmetrical parts of rr
Nmat = repmat([1:N],[N 1]);

diag = r(find(Nmat==Nmat'))';
top = r(find(Nmat>Nmat'))';
bottom = r(find(Nmat<Nmat'))';
offdiag = r(find(Nmat~=Nmat'))';
all = r(1:size(r,1)^2);

switch flag
  case 'd'
    v = diag;
  case 't'
    v = top;
  case 'b'
    v = bottom;
  case 'o'
    v = offdiag;
  case 'a'
    v = all;
end



