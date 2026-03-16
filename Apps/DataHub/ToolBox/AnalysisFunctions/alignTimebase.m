function dat = alignTimebase(dat,AlignEvent)
%
% DAT = alignTimebase(DAT,AlignEvent);
%
% Hanuman version 0.4 function for aligning a event timestamps, spike
% times, and analog time series relative to a particular event (such as
% stimulus onset or saccade).   
%
% INPUT ARGUMENTS 
% AlignEvent is a scalar variable is a key that specifies the column in
% dat.c that contains the timestamp for a particular event.  
%
% DAT is a hanuman data structure.  DAT must contain at least a .c field
% (which contains conditions and event timestamps), otherwise an error is
% returned.  Likewise an error is returned is any of the rows (trials) in
% dat.c contain NaN values in the column (key) referenced by the input
% argument AlignEvent. Spikes and LFPs are aligned if DAT contains fields
% .s and .lfp respectively.  
%
% EXAMPLE OF USAGE
% filename = 'poopy1.nex'
% dat = build_hanuman_database(filename);
% set_poopy_global;
% visual_response = alignTimebase(dat,STIMULUS_ONSET);
% saccadic_response = alignTimebase(dat,SACCADE_ONSET);
%
% 2011-mar-14: added feature for locking timebase. Function doesn't
% proceed if the field dat.h.timebase_locked exists.
% 
% 2011-mar-22: added handling for dat.eye field.
% 
% last modified 2011-mar-22
% dbtm
%

TIMESTAMPS = dat.h.timestamps;

if isfield(dat.h,'timebase_locked');
    error('Field dat.h.timebase_locked detected: Timebase cannot be shifted.');
end


if ~isfield(dat,'c')
    error('Database structure must contain dat.c field.');
end

Talign = dat.c(:,AlignEvent);
missing = find(isnan(Talign));
if length(missing>0)
    error('Timestamp to be aligned on contains NaNs');
end

% Shift event timestamps
T_align_events = repmat(Talign,[1 length(TIMESTAMPS)]);
dat.c(:,TIMESTAMPS) = dat.c(:,TIMESTAMPS) - T_align_events;


% Shift spike timestamps
if isfield(dat,'s')
    [r c s] = size(dat.s);
    T_align_spikes = repmat(Talign,[1 c s]);
    dat.s = dat.s - T_align_spikes;
end

% Shift lfp series
if isfield(dat,'lfp')
    if ~isempty(dat.lfp)
        Tpast = dat.c(:,dat.h.Taligned);
        Tfuture = dat.c(:,AlignEvent);
        deltaT = Tpast - Tfuture;  
        [r c s] = size(dat.lfp);
        for t = 1:length(deltaT)
            lfp = dat.lfp(t,:,:);
            if deltaT(t)<0
                shift = abs(deltaT(t));
                pad = nan(1,shift,s);
                crop = 1:shift;
                lfp(:,crop,:) = [];
                lfp = [lfp pad];
                dat.lfp(t,:,:) = lfp;
            elseif deltaT(t)>0
                shift = abs(deltaT(t));
                pad = nan(1,shift,s);
                crop = [length(lfp)-shift+1 : length(lfp)];
                lfp(:,crop,:) = [];
                lfp = [pad lfp];
                dat.lfp(t,:,:) = lfp;
            end
        end
        %       dat.h.Taligned = AlignEvent;
    end
end

% % Shift eye-trace series
if isfield(dat,'eye')
    if ~isempty(dat.eye)
        Tpast = dat.c(:,dat.h.Taligned);
        Tfuture = dat.c(:,AlignEvent);
        deltaT = Tpast - Tfuture;  
        [r c s] = size(dat.eye);
        for t = 1:length(deltaT)
            eye = dat.eye(t,:,:);
            if deltaT(t)<0
                shift = abs(deltaT(t));
                pad = nan(1,shift,s);
                crop = 1:shift;
                eye(:,crop,:) = [];
                eye = [eye pad];
                dat.eye(t,:,:) = eye;
            elseif deltaT(t)>0
                shift = abs(deltaT(t));
                pad = nan(1,shift,s);
                crop = [length(eye)-shift+1 : length(eye)];
                eye(:,crop,:) = [];
                eye = [pad eye];
                dat.eye(t,:,:) = eye;
            end
        end
        %        dat.h.Taligned = AlignEvent;
    end
end

dat.h.Taligned = AlignEvent;



