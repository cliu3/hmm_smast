function y = gausspdf(x,mu,invcholSigma,const)
%GAUSSPDF  Find the pdf value in a Gaussian distribution.
%   [PDFVAL] = GAUSSPDF(X,MU,INVCHOLSIGMA,CONST)
%
%   - X independent variable.
%   - MU mean of the Gaussian distribution.
%   - INVCHOLSIGMA inverse of the cholesky factorised covariance matrix.
%   - CONST the multiplicative constant in the Gaussian distribution.
%
%   Output
%
%   - PDFVAL value of the pdf at X.
%
%   Background function
%
%   Date: 12/7 - 2007, ver. 0.4
%   HMM geolocation toolbox, IMM and DIFRES

if length(x) ~= length(mu), error('x and mu must be of the same length'),end
if length(x) ~= size(invcholSigma,1), error('x and incholSigma must be of the same length'), end

sqrtxSigma = (x-mu)*invcholSigma;
xSigma = sum(sqrtxSigma.^2);
y = const .* exp(-0.5*xSigma);