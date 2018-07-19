function [EPS,R,A,B,C,U,V,W,SO] = s83tconst(SF,OR,RAD,ER,RF,F,ESQ)

% function [EPS,R,A,B,C,U,V,W,SO] = s83tconst(SF,OR,RAD,ER,RF,F,ESQ)
% convert geodesic to state plane co-ordinates
% original BY E. CARLSON / subroutine
% 15.03.12 M.P.

% TRANSVERSE MERCATOR PROJECTION
% PRECOMPUTATION OF CONSTANTS

% ******************** SYMBOLS AND DEFINITIONS  **********************
%   ER is equatorial radius of the ellipsoid (= major semiaxis).
%   RF is reciprocal of flattening of the ellipsoid.
%   SF is scale factor of the central meridian.
%   OR is southernmost parallel of latitude (in radians) for which
%     the northing coordinate is zero at the central meridian.
%   R, A, B, C, U, V, W are ellipsoid constants used for computing
%     meridional distance from latitude and vice versa.
%   SO is meridional distance (multiplied by the scale factor) from
%     the equator to the southernmost parallel of latitude.
% ********************************************************************

      %F=1.0/RF;
      %ESQ=(F+F-F^2);
      EPS=ESQ/(1.-ESQ);
      PR=(1.-F)*ER;
      EN=(ER-PR)/(ER+PR);
      A=-1.500*EN + (9./16.)*EN^3;
      B= 0.937500*EN^2 - (15./32.)*EN^4;
      C=-(35./48.)*EN^3;
      U=1.500*EN - (27./32.)*EN^3;
      V=1.312500*EN^2 - (55./32.)*EN^4;
      W=(151./96.)*EN^3;
      R=ER*(1.-EN)*(1.-EN^2)*(1.+2.2500*EN^2+(225./64.)*EN^4);
      OMO=OR + A*sin(2.*OR) + B*sin(4.*OR) + C*sin(6.*OR);
      SO=SF*R*OMO;
