function [FI,LAM,CONV,KP] = s83skewr(NORTH,EAST,GAMC,FN,FE,RAD,ER,RF,F,ESQ,A,B,C,D,SGO,CGO,SGC,CGC,XI,LONO,F0,F2,F4,F6)

% function [FI,LAM,CONV,KP] = s83skewr(NORTH,EAST,GAMC,FN,FE,RAD,ER,RF,F,ESQ,A,B,C,D,SGO,CGO,SGC,CGC,XI,LONO,F0,F2,F4,F6)
% convert geodesic to state plane co-ordinates
% original BY E. CARLSON / subroutine
% 15.03.12 M.P.

% OBLIQUE MERCATOR PROJECTION
% CONVERSION OF GRID COORDS TO GEODETIC COORDS

      U=SGC*(EAST-FE)+CGC*(NORTH-FN);
      V=CGC*(EAST-FE)-SGC*(NORTH-FN);
      R=sinh(V/D);
      S=cosh(V/D);
      SINE=sin(U/D);
      Q=(log((S-SGO*R+CGO*SINE)/(S+SGO*R-CGO*SINE))/2.00-C)/B;
      EX=exp(Q);
      XR=atan((EX-1.00)/(EX+1.00))*2.00;
      CS=cos(XR);

      LAT=XR+(F0+F2*CS*CS+F4*CS^4+F6*CS^6)*CS*sin(XR);
      LON=LONO-atan((SGO*SINE+CGO*R)/cos(U/D))/B;
      
% ***
      
      FI = LAT;
      LAM = LON;

      E=sqrt(ESQ);
      SINB=sin(FI);
      COSB=cos(FI);
      DL=(LAM-LONO)*B;
      SINDL=sin(DL);
      COSDL=cos(DL);
      Q=(log((1+SINB)/(1-SINB)) - E*log((1+E*SINB)/(1-E*SINB)))/2.0;
      R=SINH(B*Q+C);
      S=COSH(B*Q+C);
      U=D*atan((CGO*R-SGO*SINDL)/COSDL);
      V=D*log((S-SGO*R-CGO*SINDL)/(S+SGO*R+CGO*SINDL))/2.0;
      CONV=atan((SGO-CGO*SINDL*R)/(CGO*COSDL*S))-GAMC;
      KP=XI*sqrt(1-ESQ*SINB^2)*cos(U/D)/COSB/COSDL;

