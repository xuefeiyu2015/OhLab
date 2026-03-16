function dat = sortTrials(dat,key,varargin);
%
% dat = sortTrials(dat,key);
%
% Resorts the trails in data a structure.
% INPUT ARGUMENTS:
% DAT is a hanuman 0.5 data structure
% KEY specified the column to sort by.
% 
% By default sorting is done in descending order.  To sort in ascending
% order, Use third optional input argument -1:
%
% dat = sortTrials(dat,key,-1);
%
% To restore DAT to the default order, use:
% DAT = sortTrials(DAT,TR);
%
% 2012-apr-16: added blp handling field.
% 
% last modified 2012-apr-16
% dbtm


if ~isempty(varargin)
    if ~ismember(varargin{1},[-1 1]);
        help sortTrials;
        error;
    end
    sign = varargin{1};
else
    sign = +1;
end    

Ntrials = size(dat.c,1);
mat = [[1:Ntrials]' dat.c(:,key)];
mat = sortrows(mat,sign*2);
ord = mat(:,1);

dat.c = dat.c(ord,:);

if isfield(dat,'s')
    dat.s = dat.s(ord,:,:);
end

if isfield(dat,'lfp')
    dat.lfp = dat.lfp(ord,:,:);
end

if isfield(dat,'blp')
    dat.blp = dat.blp(ord,:,:);
end

