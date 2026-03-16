function dat = selectTrials(dat,varargin)

% dat = selectTrials(dat,key1,value1, ... keyN,valueN)
%
% Function for selecting trials from a hanuman database.
% Examples of usage:
%   LagSixteenSymm = selectTrials(dat,LAG,16,TYPE,SYMM);
% 
% 2011-mar-22 Added dat.eye field.
% 2012-apr-16 Added dat.blp field.
% 2013-apr-06 Added dat.mua field.
% 2013-jul-30 Added dat.wf field.
% 

% last modified 2013-jul-30
% dbtm


if mod(length(varargin),2)~=0 % if not iseven
    error('KEY and VALUE input arguments must be paired.');
end

Npairs = length(varargin)/2;
for n=1:Npairs
    key = varargin{2*(n-1)+1};
    value = varargin{2*(n-1)+2};
    dat = binary_select_trials(dat,key,value);
end

function dat = binary_select_trials(dat,key,value)

unpack;
daysOld = unique(dat.c(:,DATE))';
trials = find(ismember(dat.c(:,key),value));
dat.c = dat.c(trials,:);
daysNew  = unique(dat.c(:,DATE))';
%if ~isequal(daysOld,daysNew)
%    
%end
dat.h.session = unique(dat.c(:,SESSION));

if isfield(dat,'s')
    dat.s = dat.s(trials,:,:);    
end

if isfield(dat,'lfp')
    dat.lfp = dat.lfp(trials,:,:);    
end

if isfield(dat,'eye')
    dat.eye = dat.eye(trials,:,:);
end

if isfield(dat,'blp')
    
    dat.blp = dat.blp(trials,:,:);
end

if isfield(dat,'mua')
    dat.mua = dat.mua(trials,:,:);
end

if isfield(dat,'Zs')
    for n=1:length(dat.Zs)
        dat.Zs{n} = dat.Zs{n}(trials,:,:);
    end
end
    
if isfield(dat,'Zf')
    for n=1:length(dat.Zf)
        dat.Zf{n} = dat.Zf{n}(trials,:,:);
    end
end

if isfield(dat,'sc')
    dat.sc = dat.sc(trials,:,:);
end

% if isfield(dat,'wf')
%     if ~isequal(daysOld,daysNew)
%         gone = find(~ismember(daysOld,daysNew));
%         dat.wf.v(:,:,gone) = [];
%         Nsig = length(dat.wf.ts);
%         for s=1:length(Nsig)
%             dat.wf.ts{s}(gone) = [];
%         end
%     end
% end
