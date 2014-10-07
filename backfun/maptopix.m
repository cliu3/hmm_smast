function [px py] = maptopix(M,lat,lon)
%MAPTOPIX Uses a mapmatrix to convert from lat/long to indices.
%
%   This function is called by samptrack.m and mptrack.m
%
%   This function should not be called manually by the user.
%
%   Date: 24/7 - 2007, ver. 0.5
%   HMM geolocation toolbox, IMM and DIFRES

px = (lon - M(3,2))/M(4,2);
py = (lat - M(1,1))/M(2,1);