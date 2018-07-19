function [NORTH,EAST,CONV,KP] = s83tmgrid(FI,LAM,SF,OR,CM,FE,FN,RAD,ER,RF,F,ESQ,EPS,R,A,B,C,U,V,W,SO)

% function [NORTH,EAST,CONV,KP] = s83tmgrid(FI,LAM,SF,OR,CM,FE,FN,RAD,ER,RF,F,ESQ,EPS,R,A,B,C,U,V,W,SO)
% convert geodesic to state plane co-ordinates
% original BY E. CARLSON / subroutine
% 15.03.12 M.P.

% TRANSVERSE MERCATOR PROJECTION
% CONVERSION OF GEODETIC COORDINATES TO GRID COORDINATES

% *****************  SYMBOLS AND DEFINITIONS *************************
%   Latitude positive north, longitude positive west.  All angles are
%     in radian measure.
%   N, E are northing and easting coordinates respectively.
%   LAT, LON are latitude and longitude respectively.
%   CONV is convergence.
%   KP is point scale factor.
%   ER is equatorial radius of the ellipsoid (= major semiaxis).
%   ESQ is the square of first eccentricity of the ellipsoid.
%   EPS is the square of second eccentricity of the ellipsoid.
%   CM is the central meridian of the projection zone.
%   FE is false easting value at the central meridian.
%   FN is "false northing" at the southernmost latitude, usually zero.
%   SF is scale factor at the central meridian.
%   SO is meridional distance (multiplied by the scale factor) from
%     the equator to the southernmost parallel of latitude for the zone.
%   R is the radius of the rectifying sphere (used for computing
%     meridional distance from latitude and vice versa).
%   A, B, C, U, V, W are other precomputed constants for determination
%     of meridional distance from latitude and vice versa.
%
%   The formula used in this subroutine gives geodetic accuracy within
%   zones of 7 degrees in east-west extent.  Within State transverse
%   Mercator projection zones, several minor terms of the equations
%   may be omitted (see a separate NGS publication).  If programmed
%   in full, the subroutine can be used for computations in surveys
%   extending over two zones.
%
% *********************************************************************

      OM=FI + A*sin(2.0*FI) + B*sin(4.0*FI) + C*sin(6.0*FI);
      S=R*OM*SF;
      SINFI=sin(FI);
      COSFI=cos(FI);
      TN=SINFI/COSFI;
      TS=TN^2;
      ETS=EPS*COSFI^2;
      L=(LAM-CM)*COSFI;
      LS=L*L;
      RN=SF*ER/sqrt(1.0-ESQ*SINFI^2);

      A2=RN*TN/2.0;
      A4=(5.0-TS+ETS*(9.0+4.0*ETS))/12.0;
      A6=(61.0+TS*(TS-58.0)+ETS*(270.0-330.0*TS))/360.0;
      A1=-RN;
      A3=(1.0-TS+ETS)/6.0;
      A5=(5.0+TS*(TS-18.0)+ETS*(14.0-58.0*TS))/120.0;
      A7=(61.0-479.0*TS+179.0*TS^2-TS^3)/5040.0;
      NORTH=S-SO + A2*LS*(1.0+LS*(A4+A6*LS)) +FN;
      EAST=FE + A1*L*(1.0+ LS*(A3+LS*(A5+A7*LS)));

% *** CONVERGENCE
      C1=-TN;
      C3=(1.0+3.0*ETS+2.0*ETS^2)/3.0;
      C5=(2.0-TS)/15.0;
      CONV=C1*L*(1.0+LS*(C3+C5*LS));

% *** POINT SCALE FACTOR
      F2=(1.0+ETS)/2.0;
      F4=(5.0-4.0*TS+ETS*( 9.0-24.0*TS))/12.0;
      KP=SF*(1.0+F2*LS*(1.0+F4*LS));

