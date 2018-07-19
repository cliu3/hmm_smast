function [Q] = s83qinline(E,S)

% function [E,SINFO,RB,K,KO,NO,G] = s83lconst(NB,FIS,FIN,FIB,RAD,ER,RF,F,ESQ)
% convert geodesic to state plane co-ordinates
% original BY E. CARLSON / subroutine
% 15.03.12 M.P.

Q =(log((1+S)/(1-S))-E*log((1+E*S)/(1-E*S)))/2.0;
