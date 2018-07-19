function [E, N] = LatLongToStatePlane(lat, lon, ICODE)
% function [E, N, C, K, S, P] = LatLongToStatePlane(lat, lon, ICODE)

% function [E, N, C, K, S, P] = LatLongToStatePlane(lat, lon, ICODE)
% convert geodesic to state plane co-ordinates
% original BY E. CARLSON / main program
%
% input:
% - lat: geogr. Latitude (N+/S-)
% - lon: geogr. Longitude (W-/E+)
% - ICODE: numeric Code for State Plane
% output:
% - E: Easting
% - N: Northing
% - C: Convergence
% - K: Scale Factor
% - S: State Plane No.
% - P: Projection Type
%
% 15.03.12 M.P.

% Bogenmasz
RHO = pi/180.0;
RAD = 1/RHO;

% Ellipsoid Constants
ER=6378137.00;
RF=298.257222101;
F=1.00/RF;
ESQ=(F+F-F*F);

% load STATE PLANE COORDINATE TABLES
[IZC,AP,SPCC,UTMC,ZN] = s83tblspc(0);

% conversion to STATE PLANE COORDINATES
[NORTH,EAST,CONV,KP,IZ,mess] = s83gppc83(lat,lon,ICODE,IZC,SPCC,UTMC,AP,ZN,RAD,ER,RF,F,ESQ);

% Rueckgabe
E = EAST;
N = NORTH;
C = CONV;
K = KP;
S = ZN(IZ).s;
P = AP(IZ).s;
% mess
