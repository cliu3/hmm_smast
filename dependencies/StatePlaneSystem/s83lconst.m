function [E,SINFO,RB,K,KO,NO,G] = s83lconst(NB,FIS,FIN,FIB,RAD,ER,RF,F,ESQ)

% function [E,SINFO,RB,K,KO,NO,G] = s83lconst(NB,FIS,FIN,FIB,RAD,ER,RF,F,ESQ)
% convert geodesic to state plane co-ordinates
% original BY E. CARLSON / subroutine
% 15.03.12 M.P.

% LAMBERT CONFORMAL CONIC PROJECTION, 2 STANDARD PARALLELS
% PRECOMPUTATION OF CONSTANTS

% ************************ SYMBOLS AND DEFINITIONS *********************
%       Latitude positive north, in radian measure.
%       ER is equatorial radius of the ellipsoid (= major semiaxis).
%       RF is reciprocal of flattening of the ellipsoid.
%       FIS, FIN, FIB are respecitvely the latitudes of the south
%         standard parallel, the north standard parallel, and the
%         southernmost parallel.
%       ESQ is the square of first eccentricity of the ellipsoid.
%       E is first eccentricity.
%       SINFO = sin(FO), where FO is the central parallel.
%       RB is mapping radius at the southernmost latitude.
%       K is mapping radius at the equator.
%       NB is false northing for the southernmost parallel.
%       KO is scale factor at the central parallel.
%       NO is northing of intersection of central meridian and parallel.
%       G is a constant for computing chord-to-arc corrections.
% ***********************************************************************

      %F=1./RF;
      %ESQ=F+F-F^2;
      E=sqrt(ESQ);
      SINFS=sin(FIS);
      COSFS=cos(FIS);
      SINFN=sin(FIN);
      COSFN=cos(FIN);
      SINFB=sin(FIB);

      QS=s83qinline(E,SINFS);
      QN=s83qinline(E,SINFN);
      QB=s83qinline(E,SINFB);
      W1=sqrt(1.-ESQ*SINFS^2);
      W2=sqrt(1.-ESQ*SINFN^2);
      SINFO=log(W2*COSFS/(W1*COSFN))/(QN-QS);
      K=ER*COSFS*exp(QS*SINFO)/(W1*SINFO);
      RB=K/exp(QB*SINFO);
      QO=s83qinline(E,SINFO);
      RO=K/exp(QO*SINFO);
      COSFO=sqrt(1.-SINFO^2);
      KO=sqrt(1.-ESQ*SINFO^2)*(SINFO/COSFO)*RO/ER;
      NO=RB+NB-RO;
      G=(1-ESQ*SINFO^2)^2/(2*(ER*KO)^2)*(1-ESQ);

