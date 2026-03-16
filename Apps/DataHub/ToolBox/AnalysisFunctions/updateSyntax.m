% script for updating syntax of dada code to standard format. Performs
% global replace to convert old syntax to new syntax for all utility
% functions in standard library, for instance:
%
% OLD             NEW
% selectTrials   selectTrails
% getSdf         getSdf
% countSpikes    countSpikes
%
% operates on all *.m files in current working directory.
%
% 2018-oct-3
% dbtm

n=1;
% L(n).old = 'alignTimebase'; L(n).new = 'alignTimebase';n=n+1;
% L(n).old = 'countSpikes'; L(n).new = 'countSpikes';n=n+1;
% L(n).old = 'findSaccades'; L(n).new = 'findSaccades';n=n+1;
% L(n).old = 'getCode'; L(n).new = 'getCode';n=n+1;
% L(n).old = 'getEyeVelocity'; L(n).new = 'getEyeVelocity';n=n+1;
% L(n).old = 'getLfp'; L(n).new = 'getLfp';n=n+1;
% L(n).old = 'getPsth'; L(n).new = 'getPsth';n=n+1;
% L(n).old = 'getSaccadeAmplitude'; L(n).new = 'getSaccadeAmplitude';n=n+1;
% L(n).old = 'getSdf'; L(n).new = 'getSdf';n=n+1;
% L(n).old = 'getTs'; L(n).new = 'getTs';n=n+1;
% L(n).old = 'plotLfp'; L(n).new = 'plotLfp';n=n+1;
% L(n).old = 'plotRasters'; L(n).new = 'plotRasters';n=n+1;
% L(n).old = 'plotSdf'; L(n).new = 'plotSdf';n=n+1;
% L(n).old = 'selectSpikename'; L(n).new = 'selectSpikename';n=n+1;
% L(n).old = 'sortTrials'; L(n).new = 'sortTrials';n=n+1;
% L(n).old = 'selectTrials'; L(n).new = 'selectTrials';n=n+1;
% L(n).old = 'unpackHeader'; L(n).new = 'unpackHeader';n=n+1;

%L(n).old = ''; L(n).new = '';n=n+1;

%L(n).old = ''for ''; L(n).new = 'for dopa';n=n+1;
L(n).old = ''; L(n).new = '';n=n+1;
%L(n).old = 'Fractals'; L(n).new = 'dopaFractals';n=n+1;
%L(n).old = 'pp.dada'; L(n).new = 'pp.dada';n=n+1;


for n=1:length(L)
    old = L(n).old;
    new = L(n).new;
    perlcmd = ['!perl -pi -w -e ''s/' old '/' new '/g;'' *.m'];
    disp(perlcmd);
    eval(perlcmd)
end
!perl -pi -w -e 's/gaussSmooth/gaussSmooth/g;' *.m


