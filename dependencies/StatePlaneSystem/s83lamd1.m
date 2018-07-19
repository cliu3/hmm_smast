function [NORTH,EAST,CONV,KP] = s83lamd1(FI,LAM,CM,EO,NB,RAD,ER,RF,F,ESQ,E,SINFO,RB,K,KO,NO,G)

% function [NORTH,EAST,CONV,KP] = s83lamd1(FI,LAM,CM,EO,NB,RAD,ER,RF,F,ESQ,E,SINFO,RB,K,KO,NO,G)
% convert geodesic to state plane co-ordinates
% original BY E. CARLSON / subroutine
% 15.03.12 M.P.

% LAMBERT CONFORMAL CONIC PROJECTION, 2 STANDARD PARALLELS
% CONVERSION OF GEODETIC COORDINATES TO GRID COORDINATES

% ************************ SYMBOLS AND DEFINITIONS *********************
%       Latitude positive north, longitude positive west.  All angles
%         are in radian measure.
%       FI, LAM are latitude and longitude respectively.
%       NORTH, EAST are northing and easting coordinates respectively.
%       NORTH EQUALS Y PLANE AND EAST EQUALS THE X PLANE.
%       CONV is convergence.
%       KP is point scale factor.
%       ER is equatorial radius of the ellipsoid (= major semiaxis).
%       ESQ is the square of first eccentricity of the ellipsoid.
%       E is first eccentricity.
%       CM is the central meridian of the projection zone.
%       EO is false easting value at the central meridian.
%       NB is false northing for the southernmost parallel of the
%           projection, usually zero.
%       SINFO = sin(FO), where FO is the central parallel.  This is a
%         precomputed value.
%       RB is mapping radius at the southernmost latitude. This is a
%         precomputed value.
%       K is mapping radius at the equator.  This is a precomputed
%         value.
%
% ***********************************************************************

      SINLAT=sin(FI);
      COSLAT=cos(FI);
      CONV=(CM-LAM)*SINFO;

      Q=(log((1+SINLAT)/(1-SINLAT))-E*log((1+E*SINLAT)/(1-E*SINLAT)))/2.0;
      RPT=K/exp(SINFO*Q);
      NORTH=NB+RB-RPT*cos(CONV);
      EAST=EO+RPT*sin(CONV);
      WP=sqrt(1.0-ESQ*SINLAT^2);
      KP=WP*SINFO*RPT/(ER*COSLAT);

