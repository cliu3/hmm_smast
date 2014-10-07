function [loglikval,phi] = likelihood(s,db,td,LIK)
%LIKELIHOOD Evaluate the likelihood function of the parameters.
%
%   This function is called by hmmgeolocate.m
%
%   This function should not be called manually by the user.
%
%   Date: 12/12 - 2007, ver. 0.51
%   HMM geolocation toolbox, DTU Informatics and DTU Aqua

if size(s) ==1
    S = [s s];
else
    S = s;
end

[phi,normaliser] = hmmfilter(S,db,td,LIK);
loglikval = sum(-log(normaliser));