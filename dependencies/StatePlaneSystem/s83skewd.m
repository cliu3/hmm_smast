function [NORTH,EAST,CONV,KP] = s83skewd(FI,LAM,GAMC,FN,FE,RAD,ER,RF,F,ESQ,A,B,C,D,SGO,CGO,SGC,CGC,XI,LONO,F0,F2,F4,F6)

% function [NORTH,EAST,CONV,KP] = s83skewd(FI,LAM,GAMC,FN,FE,RAD,ER,RF,F,ESQ,A,B,C,D,SGO,CGO,SGC,CGC,XI,LONO,F0,F2,F4,F6)
% convert geodesic to state plane co-ordinates
% original BY E. CARLSON / subroutine
% 15.03.12 M.P.

% OBLIQUE MERCATOR PROJECTION
% CONVERSION OF GEODETIC COORDINATES TO GRID COORDINATES

      E=sqrt(ESQ);
      SINB=sin(FI);
      COSB=cos(FI);
      DL=(LAM-LONO)*B;
      SINDL=sin(DL);
      COSDL=cos(DL);
      Q=(log((1+SINB)/(1-SINB)) - E*log((1+E*SINB)/(1-E*SINB)))/2.0;
      R=sinh(B*Q+C);
      S=cosh(B*Q+C);
      U=D*atan((CGO*R-SGO*SINDL)/COSDL);
      V=D*log((S-SGO*R-CGO*SINDL)/(S+SGO*R+CGO*SINDL))/2.0;
      NORTH=U*CGC-V*SGC+FN;
      EAST=U*SGC+V*CGC+FE;
      CONV=atan((SGO-CGO*SINDL*R)/(CGO*COSDL*S))-GAMC;
      KP=XI*sqrt(1-ESQ*SINB^2)*cos(U/D)/COSB/COSDL;
