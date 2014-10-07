function paddtrack(track,varargin)
%PADDTRACK  Add a track to a plot created by pplottrack.
%   PADDTRACK(TRACK,OPTIONS)
%
%   - TRACK a track struct created eg. with mptrack.
%
%     Optional arguments
%
%   - OPTIONS specify options as for the basic plot function
%   see help plot for further details.
%
%   The function requires the M_Map package including the 
%   high resolution coast line.
%   See http://www.eos.ubc.ca/~rich/map.html
%
%  EXAMPLE   
%   PADDTRACK(samptrack,'-xk','linewidth',4)
%
%   Date: 13/7 - 2007, ver. 0.5
%   HMM geolocation toolbox, IMM and DIFRES

hold on
m_plot(track.long,track.lat,varargin{:})
hold off