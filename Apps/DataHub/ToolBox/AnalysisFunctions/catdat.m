function dat = catdat(varargin)
% 
% DAT = catdat(dat1,dat2,dat3,...,datN);
% 
% Function for concatinating two or more hanuman data structures.  
% This utility is under development - the header field contains the
% fnames and dates that went into the  combined field, but is
% otherwise a copy of the first arugment's header field, dat1.h.  
%
% To print the files in the header in a more readible format:
%
% sprintf(dat.h.fname)
% 
% Output copied to dat.h.catFnames
% 
% 2013-jul-30 Revised handling of dat.wf and dat.h.session fields.
%
% last modified 2013-jul-30
% dbtm

if length(varargin)<1
    help catdat
    error
elseif length(varargin)==1
    dat = varargin{1};
elseif length(varargin)==2
    dat = binary_catdat(varargin{1},varargin{2});
elseif length(varargin)>2
    dat1 = binary_catdat(varargin{1},varargin{2});
    dat2 = varargin{3};
    dat = binary_catdat(dat1,dat2);
    therest = varargin;
    therest(1:3) = [];
    while ~isempty(therest)
        dat = binary_catdat(dat,therest{1});
        therest(1) = [];                
    end    
end
unpack;
dat.h.session = unique(dat.c(:,SESSION))';

function dat = binary_catdat(dat1,dat2);

if isfield(dat1,'h')
    header = dat1.h;
    header.fname = [dat1.h.fname '\n' dat2.h.fname];
    header.date = [dat1.h.date ; dat2.h.date];
    header.catFnames = sprintf(header.fname);
    dat.h = header;
end

if isfield(dat1,'c')
    CC = [dat1.c ; dat2.c];    
    dat.c = CC;
end

if isfield(dat1,'s')
    y1 = size(dat1.s,2);
    z1 = size(dat1.s,3);    
    y2 = size(dat2.s,2);
    z2 = size(dat2.s,3);    
    if ~isequal(dat1.h.snames,dat2.h.snames)
        error('Mismatch between header.snames field between dat1 and dat2!!!');
    end
    if isequal([y1 z1],[y2 z2])
        SS = [dat1.s ; dat2.s];    
    else 
        error('Not all dat.s fields of comparable dimension!');
    end
    dat.s = SS;
end

if isfield(dat1,'wf')
    Nsig1 = size(dat1.h.snames,1);
    Nsig2 = size(dat2.h.snames,1);
    if (Nsig1>1)|(Nsig2>1)
        warning(['Attempting to catdat WF fields with more than one ' ...
                 'spike!!! I''m not sure if that will work. -dbtm, ' ...
                 '2013-jul-31)']);
    end
    if ~isequal(Nsig1,Nsig2)
        error('Mismatch between header.snames field between dat1 and dat2!!!');
    end
    if ~isequal(length(dat1.wf.t),length(dat2.wf.t))
        error('Mismatch between wf.t timestamps between dat1 and dat2!!!');
    end
    days1 = floor(datenum(dat1.h.date))';
    days1 = unique(days1);
    days2 = floor(datenum(dat2.h.date))';
    days2 = unique(days2);    
    if isequal(days1,days2)
        WF = dat1.wf;
    else
        WF = dat1.wf;
        % get rid of redundant entries in wf.v and wf.ts
        dupes = find(ismember(days2,days1));
        if ~isempty(dupes)
            days2(dupes) = [];
            dat2.wf.v(:,:,dupes) = [];
            for s = 1:Nsig1
                WF.ts{s} = dat1.wf.ts{s};
                dat2.wf.ts{s}(dupes) = [];
            end
        end
        for s = 1:Nsig1        
            for d = 1:length(days2)
                WF.ts{s}{d+length(days1)} = dat2.wf.ts{s}{d};
            end
        end
    end
    WF.v = cat(3,dat1.wf.v,dat2.wf.v);
    % sort wf.v in ascending order of days
    daysAll = [days1 days2];
    [b reord] = unique(daysAll);
    daysAll = daysAll(reord);
    WF.v = WF.v(:,:,reord);
    % do likewise unto wf.ts        
    for s = 1:Nsig1
        for d = 1:length(daysAll)
            WF.ts{s}(:) = WF.ts{s}(reord);
        end
    end
    dat.wf = WF;
end
    
%end

if isfield(dat1,'lfp')
    y1 = size(dat1.lfp,2);
    z1 = size(dat1.lfp,3);    
    y2 = size(dat2.lfp,2);
    z2 = size(dat2.lfp,3);    
    if isequal([y1 z1],[y2 z2])
        LFP = [dat1.lfp ; dat2.lfp];    
    else 
        error('Not all dat.lfp fields of comparable dimension!');
    end
    dat.lfp = LFP;
    try
        isequal(dat1.t,dat2.t);
        dat.t = dat1.t;
    catch
        error('Not all dat.t fields match! Check dat.h.Taligned');
    end
end

if isfield(dat1,'blp')
    y1 = size(dat1.blp,2);
    z1 = size(dat1.blp,3);    
    y2 = size(dat2.blp,2);
    z2 = size(dat2.blp,3);    
    if isequal([y1 z1],[y2 z2])
        BLP = [dat1.blp ; dat2.blp];    
    else 
        error('Not all dat.blp fields of comparable dimension!');
    end
    dat.blp = BLP;
    try
        isequal(dat1.t,dat2.t);
        dat.t = dat1.t;
    catch
        error('Not all dat.t fields match! Check dat.h.Taligned');
    end
end

