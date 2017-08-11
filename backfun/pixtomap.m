function [lat long] = pixtomap(M,px,py)
%PIXTOMAP Uses a mapmatrix to convert from indices to lat/long.
%
%   This function is called by samptrack.m and mptrack.m
%
%   This function should not be called manually by the user.
%
%   Date: 24/7 - 2007, ver. 0.5
%   HMM geolocation toolbox, IMM and DIFRES

lgt = length(px);
latlong =[ones(lgt,1) py ones(lgt,1) px]*M;
lat  = latlong(:,1);
long = latlong(:,2);