function [T] = s83tinline(X,Y)

% function function [T] = s83tinline(X,Y)
% convert geodesic to state plane co-ordinates
% original BY E. CARLSON / subroutine
% 15.03.12 M.P.

T =X+Y/60.00;
