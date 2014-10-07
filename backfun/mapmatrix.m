function M = mapmatrix(y,x,dy,dx)
%MAPMATRIX Create a matrix that converts indices to lat/long.
%
%   This function is called by samptrack.m and mptrack.m
%
%   This function should not be called manually by the user.
%
%   Date: 24/7 - 2007, ver. 0.5
%   HMM geolocation toolbox, IMM and DIFRES


% y is pos (1,1) in latitude array
% x is pos (1,1) in longitude array
% dy is increment in latitude direction
% dx is increment in longitude direction

M      = zeros(4,2);
M(1,1) = y-dy;
M(2,1) = dy;
M(3,2) = x-dx;
M(4,2) = dx;