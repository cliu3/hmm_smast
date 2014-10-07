function tr=proboftrack(tr,result,LIK)
%PROBOFTRACK  Calculate the probability of a track given the observations.
%             The result is found in the tr.avgP field in the output.
%   [TR] = PROBOFTRACK(TR,RESULT,LIK)
%
%   - TR Track to calculate probability of.
%   - RESULT a result struct created with hmmgeoloc.
%   - LIK an output from the datalikelihood function.
%
%  EXAMPLE   
%   newmpt = PROBOFTRACK(mpt,result,LIK);
%   Here the newmpt.avgP contains the estimated probability.
%
%   Date: 12/12 - 2007, ver. 0.51
%   HMM geolocation toolbox, IMM and DIFRES

[row col icalc] = size(result.smooth);
no = size(tr.lat,2);
tr.P = ones(no,icalc);
tr.P_trans = ones(no,icalc-1);
tr.L = zeros(no,icalc);
names = fieldnames(LIK);
names = names(~strcmp(names,'type'));
names = names(~strcmp(names,'mode'));
numnames = length(names);

% Combine all data likelihood to one array in Ltotal %
Ltotal = ones(row,col,icalc-1);
for j = 1:numnames
    Ltotal = Ltotal .* LIK.(names{j});
end

% Define transition probabilities (convolution kernels)
% s = result.D*result.D2s;
% unc    = sqrt(2*s(1)); 
% ks = ceil(unc*10+1); ks = ks + mod(ks,2) + 1;
% ks1  = max([15 ks]); %ksize = 21;
% kern1  = gausskern(ks1,unc);
% ks1 = ceil(ks1/2);
% unc    = sqrt(2*s(2));
% ks = ceil(unc*10+1); ks = ks + mod(ks,2) + 1;
% ks2  = max([15 ks]); %ksize = 21;
% kern2  = gausskern(ks2,unc);
% ks2 = ceil(ks2/2);

s = result.D*result.D2s;
par1.covmat = 2*s(1)*eye(2);
par2.covmat = 2*s(2)*eye(2);
kern1 = makekern2(par1);
kern2 = makekern2(par2);
ks1 = ceil(max(size(kern1))/2);
ks2 = ceil(max(size(kern2))/2);

pred = zeros(row,col,icalc);
for i = 2:icalc
    pred(:,:,i) = normalise(result.pred(:,:,i));
end

% Cycle through tracks
for k = 1:no
    % Cycle through days for specific track
    for i = 1:icalc-1
        x = tr.long_pix_clean(i,k);
        xp = tr.long_pix_clean(i+1,k);
        y = tr.lat_pix_clean(i,k);
        yp = tr.lat_pix_clean(i+1,k);
        % tr.P(k,i) = result.smooth(y,x,i);
        %dlong = abs(tr.long_pix_clean(i+1,k) - tr.long_pix_clean(i,k));
        %dlat  = abs(tr.lat_pix_clean(i+1,k) - tr.lat_pix_clean(i,k));
        dlong = abs(xp - x);
        dlat  = abs(yp - y);
        switch result.behav(i)
            case 1
                ks = ks1; kern = kern1;
            case 2
                ks = ks2; kern = kern2;
        end
        tr.P_trans(k,i) = kern(ks+dlat,ks+dlong);
        tr.P(k,i) = result.smooth(y,x,i);
        tr.P2(k,i) = result.phi(y,x,i) / pred(yp,xp,i+1);
        tr.L(k,i+1) = log(Ltotal(tr.lat_pix_clean(i+1,k),tr.long_pix_clean(i+1,k),i)*kern(ks+dlat,ks+dlong));
    end
    i = icalc;
    %tr.P(k,i) = result.smooth(tr.lat_pix_clean(i,k),tr.long_pix_clean(i,k),i);
    x = tr.long_pix_clean(i,k);
    y = tr.lat_pix_clean(i,k);
    tr.P(k,i) = result.smooth(y,x,i);
    tr.P2(k,i) = result.phi(y,x,i);
    tr.avgP(k) = mean(tr.P(k,:));
    tr.avgP_trans(k) = mean(tr.P_trans(k,:));
    tr.avgL(k) = mean(tr.L(k,:));
    tr.logtagprob(k) = sum(log(tr.P2(k,:))) + sum(log(tr.P_trans(k,:)));
end

