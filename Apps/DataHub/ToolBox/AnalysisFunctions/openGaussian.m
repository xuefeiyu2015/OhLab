function y = openGaussian(x,mu,sigma)
%
% Open source alternative to the mathworks normpdf
%
% last modified 2012-aug-9
% dbtm

try
    y = normpdf(x,mu,sigma);
catch
    %   disp('bypassing stats toolbox. Eat me, mathworks!');
    y = 1/[sigma*sqrt(2*pi)]*exp([-[x-mu].^2]/[2*sigma.^2]);
end

