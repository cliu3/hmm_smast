function [COSHI] = s83coshinline(X)

% function [COSHI] = s83coshinline(X)
% convert geodesic to state plane co-ordinates
% original BY E. CARLSON / subroutine
% 15.03.12 M.P.

COSHI=log(X+sqrt(X*X-1));
