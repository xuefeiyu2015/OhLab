function amp = getSaccadeAmplitude(x,y,t,tstart,tstop);
%
%amp =  getSaccadeAmplitude(x,y,t,tstart,tstop);
%
% Returns the amplitude in degrees for all saccade times specified
% by tstart and tstop. x,y,and t are vectors of analog horizontal
% position, vertical position, and corresponding time base. Units
% in x and y must already be converted to degrees.
%
% tstart and tstop can be obtained using the function
% findSaccades.m  
%
% last modified 2011-mar-24
% dbtm


%blinkthresh = 0.5*[1:1000]+5;

if isempty(tstart)
    amp = [];
    %    blink = [];
else
    duration = round(tstop-tstart);
    for s = 1:length(tstart)
        onset = max(find(t<tstart(s)));
        offset = min(find(t>tstop(s)));
        deltax = x(offset)-x(onset);
        deltay = y(offset)-y(onset);
        amp(s) = sqrt(deltax^2 + deltay^2);
        
        %       % decide if the saccade is a blink
        %cutoff = blinkthresh(duration);
        %if amp(s)>cutoff
        %    blink(s) = 1;
        %else
        %    blink(s) = 0;
        %end
    end
end