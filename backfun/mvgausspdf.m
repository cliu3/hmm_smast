function p = mvgausspdf(x,mu,sigma)
%MVGAUSSPDF  Find the pdf value in a Gaussian distribution.
%   [PDFVAL] = MVGAUSSPDF(X,MU,SIGMA)
%
%   - X independent variable.
%   - MU mean of the Gaussian distribution.
%   - SIGMA covariance matrix.
%
%   Output
%
%   - PDFVAL value of the pdf at X.
%
%   Background function
%
%   Date: 31/7 - 2007, ver. 0.5
%   HMM geolocation toolbox, IMM and DIFRES

d = length(x);

[s,err] = chol(sigma);
if err ~= 0
    error('Covariance matrix was not positive definite. Bad parameters for the datalikelihood may be defined.');
end
x = (x - mu) / s;

p = exp(-0.5*sum(x.^2) - sum(log(diag(s))) - d*log(2*pi)/2);