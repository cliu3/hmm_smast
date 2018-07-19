function [A,B,C,D,SGO,CGO,SGC,CGC,XI,LONO,F0,F2,F4,F6] = s83oconst(LONC,LATC,GAMC,KC,RAD,ER,RF,F,ESQ)

% function [A,B,C,D,SGO,CGO,SGC,CGC,XI,LONO,F0,F2,F4,F6] = s83oconst(LONC,LATC,GAMC,KC,RAD,ER,RF,F,ESQ)
% convert geodesic to state plane co-ordinates
% original BY E. CARLSON / subroutine
% 15.03.12 M.P.

% OBLIQUE MERCATOR PROJECTION
% COMPUTATIONS OF CONSTANTS

      E=sqrt(ESQ);
      EPS=ESQ/(1.00-ESQ);
      E2=ESQ;
      E4=E2*E2;
      E6=E2^3;
      E8=E2^4;

      C2=E2/2.00+5.00*E4/24.00+E6/12.00+13.00*E8/360.00;
      C4=7.00*E4/48.00+29.00*E6/240.00+811.00*E8/11520.00;
      C6=7.00*E6/120.00+81.00*E8/1120.00;
      C8=4279.00*E8/161280.00;

      F0=2.00*C2-4.00*C4+6.00*C6-8.00*C8;
      F2=8.00*C4-32.00*C6+80.00*C8;
      F4=32.00*C6-192.00*C8;
      F6=128.00*C8;

      SINB=sin(LATC);
      COSB=cos(LATC);
      B=sqrt(1.00+EPS*COSB^4);
      W=sqrt(1.00-ESQ*SINB*SINB);
      A=B*ER*sqrt(1.00-ESQ)/(W*W);
      QC=s83qqinline(SINB,E)
      C=s83coshinline(B*sqrt(1.00-ESQ)/W/COSB)-B*QC
      test=acosh(B*sqrt(1.00-ESQ)/W/COSB)-B*QC
      D=A*KC/B;

      SGC=sin(GAMC);
      CGC=cos(GAMC);
      SGO=SGC*COSB*ER/(A*W);
      CGO=sqrt(1.00-SGO*SGO);
      LONO=LONC+asin(SGO*sinh(B*QC+C)/CGO)/B
      EF=-SGO;
      G=CGO;
      H=EF/G;
      XI=A*KC/ER;
