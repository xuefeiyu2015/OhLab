function output = notchfilter(raw);
%
% output = notchfilter(lfp);
%
% Applies a band-stop filter (59-61 Hz) to input matrix LFP.
%
% Utility function for getting rid of line noise in LFP
% traces. Apply to the continuous raw data before snipping it up
% into trial-length segments. 
%
% LFP must be an F x T matrix with F rows for each channel and T
% columns for each time sample.
% 
% last modified 2012-may-02
% dbtm

Nchan = size(raw,2);
Nsam = size(raw,1);

Srate = 1000;
n = 2^13;
%x=1:n;
%freq = Srate*(0:(n/2))/n;


lhs = [1:n:Nsam]';
rhs = lhs-1+n;         
rhs(end) = Nsam;
%clear nmax

output = nan(Nsam,Nchan);

for c=1:Nchan
    for i=1:length(lhs)
    %    for i=1:1    
        win = lhs(i):rhs(i);
        [B A] = cheby1(5,0.5,[59 61]/[1000/2],'stop');
        output(win,c) = filter(B,A,raw(win,c));    
        %        lfp = [lfp; output];
    end
end
