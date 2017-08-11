function p=makeplotstandard(phitest)
%MAKEPLOTSTANDARD Convert the input distribution to a 2D cdf.
%
%   This function is called by pplotdistr.m
%
%   This function should not be called manually by the user.
%
%   Date: 24/7 - 2007, ver. 0.5
%   HMM geolocation toolbox, IMM and DIFRES


% Make plot standard i.e. make 2D cdf for phi
[row,col]=size(phitest);% Row and column numbers
p=phitest(:);           % convert from matrix to column
[q,i]=sort(p);          % sort p ascending
r=cumsum(q);            % make cumulated sum
p(i)=r;                 % store r in the correct index positions
p=reshape(p,row,col);