function [smooth] = smoothing(s,phi,db,td,pred)
%SMOOTHING Perform the backward sweep of the filtering.
%
%   This function is called by hmmgeolocate.m
%
%   This function should not be called manually by the user.
%
%   Date: 7/8 - 2009, ver. 0.52
%   HMM geolocation toolbox, DTU Informatics and DTU Aqua

% Initialise
[row,col] = size(db.depth);
icalc     = length(td.d24);
smooth    = zeros(row,col,icalc); % Density function
smooth(:,:,icalc) = phi(:,:,icalc); % Last est is a smoothed est as well

par1.covmat = 2*s(1)*eye(2);
par2.covmat = 2*s(2)*eye(2);
kern1 = makekern2(par1);
kern2 = makekern2(par2);


for i=icalc:-1:2
    if     td.behav(i-1) == 1
        ratio = smooth(:,:,i)./(pred(:,:,i)+eps^20);
        smooth(:,:,i-1) = phi(:,:,i-1) .* conv2(ratio,kern1,'same');
    elseif td.behav(i-1) == 2
        ratio = smooth(:,:,i)./(pred(:,:,i)+eps^20);
        smooth(:,:,i-1) = phi(:,:,i-1) .* conv2(ratio,kern2,'same');
    end
    [smooth(:,:,i-1),NO_USE] = normalise(smooth(:,:,i-1));
end