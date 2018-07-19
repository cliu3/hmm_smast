function [lat, lon] = StatePlaneToLatLong(E, N, ICODE)
% function [lat, lon, C, K, S, P] = StatePlaneToLatLong(E, N, ICODE)

% function [lat, lon, C, K, S, P] = StatePlaneToLatLong(E, N, ICODE)
% convert state plane to geodesic co-ordinates
% original BY E. CARLSON / main program
%
% input:
% - E: Easting
% - N: Northing
% - ICODE: numeric Code for State Plane
% output:
% - lat: geogr. Latitude (N+/S-)
% - lon: geogr. Longitude (W-/E+)
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
[lat,lon,CONV,KP,IZ,mess] = s83pcgp83(N,E,ICODE,IZC,SPCC,UTMC,AP,ZN,RAD,ER,RF,F,ESQ);

% Rueckgabe
lat = lat;
lon = lon;
C = CONV;
K = KP;
S = ZN(IZ).s;
P = AP(IZ).s;
% mess
