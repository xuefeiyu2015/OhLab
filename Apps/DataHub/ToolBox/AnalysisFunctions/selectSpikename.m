function dat = selectSpikename(dat,sig);
%
% dat = selectSpikename(dat,sig);
% 
% Removes all spikes but one from DAT. The remaining spike is
% specified by text string SIG, which must be a member of
% dat.h.snames.
%
% 2012-dec-16 Revised to also select dat.wf field (if present).
%
% 2013-jul-29 dat.wf.v now matrix not structure; using common
% timebase for dat.wf.t
% 
% last modified 2013-jul-31
% dbtm

spk = strmatch(sig,dat.h.snames);

if isempty(spk)
    warning(['Spike ' sig ' not found in ' dat.h.fname]);    
    dummy = nan(size(dat.s,1),size(dat.s,2));
    dat.s = dummy;
else
    dat.s = dat.s(:,:,spk);
end

if isfield(dat,'wf')
    if isempty(spk)
        dummyV = nan(1,length(dat.wf.v(1,:)));
        dummyT = nan(1,length(dat.wf.t(1,:)));        
        dat.wf.v = dummyV;
        dat.wf.t = dummyT;
        dat.wf.ts = [];
        Ndays = size(dat.wf.v,3);
        for d=1:Ndays
            dat.wf.ts{1}{d} = [];
        end
    else
        dat.wf.t = dat.wf.t;
        dat.wf.v = dat.wf.v(spk,:,:);
        dat.wf.ts= dat.wf.ts(spk);
    end
end

dat.h.snames = sig;

if isfield(dat.h,'chan')
    chan = str2num(sig(4:6));
    dat.h.chan = chan;
end
