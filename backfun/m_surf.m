function [cs,h]=m_surf(long,lat,data,varargin);
%M_SURF Draws a surface on a map
%   M_SURF(LONG,LAT,DATA,...) draws a surface on a map. Behaviour
%   is the same as for M_CONTOUR except that LONG and LAT vectors or
%   matrices must be specified.
%
%   [CS,H]=M_SURF(...) returns a contour matrix C and a vector
%   H of handles to LINE or PATCH objects for use by CLABEL.
%
%   This function is called by plottingfancy.m
%
%   This function should not be called manually by the user.
%
%   Date: 12/12 - 2007, ver. 0.51
%   HMM geolocation toolbox, IMM and DIFRES

global MAP_PROJECTION

% Have to have initialized a map first

if isempty(MAP_PROJECTION),
  disp('No Map Projection initialized - call M_PROJ first!');
  return;
end;

if min(size(long))==1 & min(size(lat))==1,
 [long,lat]=meshgrid(long,lat);
end;

[X,Y]=m_ll2xy(long,lat,'clip','on');

i=isnan(X);      % For these we set the *data* to NaN...
data(i)=NaN;

                 % And then recompute positions without clipping. THis
                 % is necessary otherwise contouring fails (X/Y with NaN
                 % is a no-no. 
if any(i(:)), [X,Y]=m_ll2xy(long,lat,'clip','off'); end;  

if any(~i(:)),
 [h]=surf(X,Y,data,varargin{:});shading interp; view(2)
 set(h,'tag','m_surf');
else
  cs=[];h=[];
end;

if nargout==0,
 clear cs h
end;

