function p = normcdf(x,mu,sigma)
%NORMCDF Find the cdf value of a Gaussian distributed number.
%
%   This function is called by datalikelihood.m
%
%   This function should not be called manually by the user.
%
%   Date: 24/7 - 2007, ver. 0.5
%   HMM geolocation toolbox, IMM and DIFRES

p = .5*erfc((mu-x)./(sigma*sqrt(2)));
