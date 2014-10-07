function [pdf,normConst]=normalise(pdf)
%NORMALISE Normalises a distribution to sum to 1.
%
%   This function is omnipresent throughout the toolbox.
%
%   This function should not be called manually by the user.
%
%   Date: 24/7 - 2007, ver. 0.5
%   HMM geolocation toolbox, IMM and DIFRES


normConst = sum(pdf(:));%+eps^20;
pdf = pdf/normConst;