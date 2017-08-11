function longfac = deglong(lat)
%DEGLONG returns the length of one degree of longitude at the latitude
%        given as input to the function.
%
%   This function is called by readdb.m, samptrack.m and mptrack.m
%
%   This function should not be called manually by the user.
%
%   Date: 12/12 - 2007, ver. 0.51
%   HMM geolocation toolbox, IMM and DIFRES

lat = lat*pi/180;
a = 6378.137;
b = 6356.7523;
longfac = pi/180*cos(lat).*sqrt( (a^4*cos(lat).^2 + b^4*sin(lat).^2)./(a^2*cos(lat).^2 + b^2*sin(lat).^2) );