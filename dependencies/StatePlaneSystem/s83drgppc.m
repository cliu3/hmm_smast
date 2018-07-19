function [NORTH,EAST,CONV,KP,IZ,mess] = s83drgppc(FI,LAM,ICODE,IZC,SPCC,UTMC,AP,ZN,RAD,ER,RF,F,ESQ)

% function [NORTH,EAST,CONV,KP,IZ,mess] = s83drgppc(FI,LAM,ICODE,IZC,SPCC,UTMC,AP,ZN,RAD,ER,RF,F,ESQ)
% convert geodesic to state plane co-ordinates
% original BY E. CARLSON / subroutine
% 15.03.12 M.P.

% DIRECTION OF LONGITUDE - E OR W
% ***/***/***
%if (sign(LAM) == 1)
%   EWFLAG = 1;
%else
%   EWFLAG = 0;
%end

%if (EWFLAG == 1)
%   LAM = (360.00/RAD) - LAM;
%end

% DO 10
%for J=1:3
   if (ICODE == 0)
		mess = 'CODE NULL';
      NORTH = -999.0;
      EAST = -999.0;
      CONV = -999.0;
      KP = -999.0;
      IZ = 0;
      return
   end
   
   IZ=0;
   % DO 20 
   for I=1:135
      if (IZC(I) == ICODE)
         IZ=I;
      end
   % 20     CONTINUE
   end
   
   if (IZ == 0)
      %WRITE(6,30) ICODE(J)
      %30     FORMAT(' IMPROPER STATE ZONE CODE-',I4)
      %GO TO 10
      mess = 'IMPROPER STATE ZONE CODE';
      NORTH = -999.0;
      EAST = -999.0;
      CONV = -999.0;
      KP = -999.0;
      IZ = 0;
      return
   elseif (AP(IZ).s == 'N')
      %WRITE(6,40)ICODE(J)
      %40     FORMAT(' THE ZONE CONSTANTS ARE NOT ','YET AVAILABLE FOR -',I4)
      %GO TO 10
      mess = ['THE ZONE CONSTANTS ARE NOT YET AVAILABLE FOR' num2str(I4)];
      NORTH = -999.0;
      EAST = -999.0;
      CONV = -999.0;
      KP = -999.0;
      IZ = 0;
      return
   elseif (AP(IZ).s == 'L')
      % PERFORM LAMBERT CONIC CONVERSION
      
      %  GET ALL THE ZONE CONSTANCES ****
      CM=SPCC(IZ,1)/RAD;
      EO=SPCC(IZ,2);
      NB=SPCC(IZ,3);
      FIS=SPCC(IZ,4)/RAD;
      FIN=SPCC(IZ,5)/RAD;
      FIB=SPCC(IZ,6)/RAD;
      
      %  FIND ZONE NAME  ********
      ZONE=ZN(IZ).s;
      
      %      COMPUTE ALL CONSTANCES FOR PROJECTION
      %CALL LCONST(ER,RF,FIS,FIN,FIB,ESQ,E,SINFO,RB,K,KO,NO,G,NB);
      [E,SINFO,RB,K,KO,NO,G] = s83lconst(NB,FIS,FIN,FIB,RAD,ER,RF,F,ESQ);
      
      %      CONVERT LAT AND LONG TO PCS
      %CALL LAMD1(FI,LAM,NORTH,EAST,CONV,KP,ER,ESQ,E,CM,EO,NB,SINFO,RB,K);
      [NORTH,EAST,CONV,KP] = s83lamd1(FI,LAM,CM,EO,NB,RAD,ER,RF,F,ESQ,E,SINFO,RB,K,KO,NO,G);

      %      PRINT OUTPUT
      %s83formgp(CARDR,NORTH,EAST,CONV,KP,ZONE,FILFLAG,J);
      
   elseif (AP(IZ).s == 'T')
      % PERFORM TRANSVERSE MERCATOR
      CM=SPCC(IZ,1)/RAD;
      FE=SPCC(IZ,2);
      OR=SPCC(IZ,3)/RAD;
      SF=1.D0-1.D0/SPCC(IZ,4);
      FN=SPCC(IZ,5);
      
      % FIND ZONE NAME  ********
      ZONE=ZN(IZ).s;
      
      if ((ZONE == 'HI 5') | (ZONE == 'GU  '))
         SF= 1.000;
      end
      
      %      COMPUTE  ALL CONSTANCES FOR PROJECTION
      % CALL TCONST(ER,RF,SF,OR,ESQ,EPS,R,A,B,C,U,V,W,SO,CM,FE,FN);
      [EPS,R,A,B,C,U,V,W,SO] = s83tconst(SF,OR,RAD,ER,RF,F,ESQ);
      
      %      CONVERT LAT AND LONG TO PCS
      %CALL TMGRID(FI,LAM,NORTH,EAST,CONV,KP,ER,ESQ,EPS,CM,FE,FN,SF,SO,R,A,B,C,U,V,W);
      [NORTH,EAST,CONV,KP] = s83tmgrid(FI,LAM,SF,OR,CM,FE,FN,RAD,ER,RF,F,ESQ,EPS,R,A,B,C,U,V,W,SO);
      
      %      PRINT OUTPUT
      %s83formgp(CARDR,NORTH,EAST,CONV,KP,ZONE,FILFLAG,J);
      
   elseif (AP(IZ).s == 'O')
      % PERFORM OBLIQUE MERCATOR
      LONC=SPCC(IZ,1)/RAD;
      FE=SPCC(IZ,2);
      FN=SPCC(IZ,3);
      GAMC=SPCC(IZ,4);
      LATC=SPCC(IZ,5)/RAD;
      KC=1.D0-1.D0/SPCC(IZ,6);
      
      % FIND ZONE NAME  ********
      ZONE=ZN(IZ).s;
      
      %      COMPUTE ALL CONSTANCES FOR PROJECTION
      %CALL OCONST(ER,RF,A,B,C,D,SGO,CGO,GAMC,SGC,CGC,XI,KC,LONO,F0,F2,F4,F6,LATC,LONC,ESQ);
      [A,B,C,D,SGO,CGO,SGC,CGC,XI,LONO,F0,F2,F4,F6] = s83oconst(LONC,LATC,GAMC,KC,RAD,ER,RF,F,ESQ);
      
      %      CONVERT LAT AND LONG TO PCS
      %CALL SKEWD(FI,LAM,U,V,NORTH,EAST,CONV,KP,B,C,D,SGO,CGO,GAMC,CGC,SGC,XI,E,ESQ,LONO,FN,FE);
      [NORTH,EAST,CONV,KP] = s83skewd(FI,LAM,GAMC,FN,FE,RAD,ER,RF,F,ESQ,A,B,C,D,SGO,CGO,SGC,CGC,XI,LONO,F0,F2,F4,F6);

      %      PRINT OUTPUT
      %s83formgp(CARDR,NORTH,EAST,CONV,KP,ZONE,FILFLAG,J);
      
   end
      
%  10  CONTINUE
%end

mess = 'done';
