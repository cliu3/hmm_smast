function [rmse rsquare ampli out Yhat1 mwh alpha beta] = lssinfit(ons,cost,sint,ts,lim)
%LSSINFIT Fit a sinewave to input data by LS.
%
%   This function is called by tidebehavextr.m
%
%   This function should not be called manually by the user.
%
%   Date: 24/7 - 2007, ver. 0.5
%   HMM geolocation toolbox, IMM and DIFRES


%p=12.42; %period in hours
%w=2*pi/(p/24); % Angular frequency

out=0;
X=[ons cost sint];
Y=ts; Y2=Y; [n m]=size(X); 
% n is number of observations
% m is number of paramters

% Solve normal equations
theta = (X'*X)\X'*Y;
Yhat1=X*theta; % predictions
res=Yhat1-Y; % residuals

rsquare = 1 - sum(res.^2)./sum((Y-mean(Y)).^2);
rmse = sqrt(sum(res.^2)/(n-m));
ampli = sqrt(theta(2)^2 + theta(3)^2);
lengthres = length(res);
df = n-m-1;
S = sum(res.^2)/(df);
mwh = theta(1);
alpha = theta(2);
beta = theta(3);
if nargin == 5
    resvar=zeros(n,1);
    for i=1:n
        %tmp = res; tmp(i)=[]; % remove residual before finding variance
        %resvar(i) = sum(tmp.^2)/(df);
        resvar(i) = S - res(i)^2/df;
    end

    H=X/(X'*X)*X'; % hat matrix
    studres = res./sqrt(resvar.*(1-diag(H))); % studentized residuals

    % Remove outliers
    index=find(abs(studres)>lim);
    outFound=length(index);
    if outFound < 7 & outFound > 0
        out = outFound;
        Y(index)=[]; X(index,:)=[];
        theta = (X'*X)\X'*Y;
        Yhat=X*theta;
        res=Yhat-Y;
        rsquare = 1-sum(res.^2)./sum((Y-mean(Y)).^2);
        rmse = sqrt(mean(res.^2));
    end

    ampli = sqrt(theta(2)^2 + theta(3)^2);
end