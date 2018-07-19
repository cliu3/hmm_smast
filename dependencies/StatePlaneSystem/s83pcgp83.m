function [lat,lon,CONV,KP,IZ,mess] = s83pcgp83(NORTH,EAST,ICODE,IZC,SPCC,UTMC,AP,ZN,RAD,ER,RF,F,ESQ)

% function [lat,lon,CONV,KP,IZ,mess] = s83pcgp83(NORTH,EAST,ICODE,IZC,SPCC,UTMC,AP,ZN,RAD,ER,RF,F,ESQ)
% convert state plane to geodesic co-ordinates
% original BY E. CARLSON / subroutine
% 15.03.12 M.P.

% Ellipsoid Constants
%ER=6378137.00;
%RF=298.257222101;
%F=1.00/RF;
%ESQ=(F+F-F*F);

% intpc.for

% conversion
[FI,LAM,CONV,KP,IZ,mess] = s83drpcgp(NORTH,EAST,ICODE,IZC,SPCC,UTMC,AP,ZN,RAD,ER,RF,F,ESQ);

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

lat=FI*RAD;
lon=-LAM*RAD;
CONV=CONV*RAD;

