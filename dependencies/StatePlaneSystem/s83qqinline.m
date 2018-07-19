function [QQ] = s83qqinline(X,E)

% function [QQ] = s83qqinline(X,E)
% convert geodesic to state plane co-ordinates
% original BY E. CARLSON / subroutine
% 15.03.12 M.P.

QQ=(log((1.00+X)/(1.00-X))-E*log((1.00+E*X)/(1.00-E*X)))/2.00;
