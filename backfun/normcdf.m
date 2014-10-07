function p = normcdf(x,mu,sigma)
%NORMCDF Find the cdf value of a Gaussian distributed number.
%
%   This function is called by datalikelihood.m
%
%   This function should not be called manually by the user.
%
%   Date: 24/7 - 2007, ver. 0.5
%   HMM geolocation toolbox, IMM and DIFRES
% gwc 
% there may be opportunity to improve this function.  Currently
% this assigns a likelihood of 1 if the fish is above bottom
% anywhere in the water column.  
% However, it is most likely that the fish is near-bottom
% an extreme implementation would be the following:
%p = erfc(((mu-x).^2)./(sigma*sigma*sqrt(2)));
% actually normpdf (built in matlab function) will do the trick


p = .5*erfc((mu-x)./(sigma*sqrt(2)));
