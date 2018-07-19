function [NORTH,EAST,CONV,KP,IZ,mess] = s83gppc83(lat,lon,ICODE,IZC,SPCC,UTMC,AP,ZN,RAD,ER,RF,F,ESQ)

% function [NORTH,EAST,CONV,KP,IZ,mess] = s83gppc83(lat,lon,ICODE,IZC,SPCC,UTMC,AP,ZN,RAD,ER,RF,F,ESQ)
% convert geodesic to state plane co-ordinates
% original BY E. CARLSON / subroutine
% 15.03.12 M.P.

% Ellipsoid Constants
%ER=6378137.00;
%RF=298.257222101;
%F=1.00/RF;
%ESQ=(F+F-F*F);

% intgp.for

% DIRECTION OF LONGITUDE - E OR W
% ***/***/***
%if (sgn(lon) = 1)
%   EWFLAG = 1;
%else
%   EWFLAG = 0;
%end

%ISEC = SLAT * 1.0E5 + 0.5;
%JSEC = SLON * 1.0E5 + 0.5;

%FI=(LD+(LM+SLAT/60.00)/60.00)/RAD;
%LAM=(LOD+(LOM+SLON/60.00)/60.00)/RAD;

FI=lat/RAD;
LAM=-lon/RAD;

% conversion
[NORTH,EAST,CONV,KP,IZ,mess] = s83drgppc(FI,LAM,ICODE,IZC,SPCC,UTMC,AP,ZN,RAD,ER,RF,F,ESQ);

CONV=CONV*RAD;

