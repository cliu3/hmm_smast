function [samptracks] = samptrack(result,LIK,no)
%SAMPTRACK  Sample random tracks from a geolocation result. 
%   [SAMPTRACKS] = SAMPTRACK(RESULT,LIK,NO)
%
%   - RESULT output from the hmmgeolocate function.
%   - LIK output from the datalikelihood function.
%
%     Optional arguments
%
%   - NO the number of tracks to be sampled.
%   default is 1
%
%     Output
%
%   - SAMPTRACKS a struct containing the coordinates for the sampled track(s).
%
%  EXAMPLE   
%   [samptracks] = SAMPTRACK(result,LIK,2);
%
%   Date: 22/10 - 2008, ver. 0.58
%   HMM geolocation toolbox, DTU Informatics and DTU Aqua

if nargin < 3, no = 1; end
if nargin < 2, error('too few inputs! type help samptrack for help'), end

if ~isfield(result,'DBname')
    result.DBname = 'tidaldb.mat';
end
disp(['Loading DB:' result.DBname])
load(result.DBname)

[row,col,icalc] = size(result.phi);
s = result.D * result.D2s;
if length(s) == 1, s = [s s]; end
samptracks.land = result.land;
samptracks.maplat  = result.maplat;
samptracks.maplong = result.maplong;
samptracks.lat_pix_clean   = zeros(icalc,no);
samptracks.long_pix_clean  = zeros(icalc,no);
samptracks.lat_pix         = zeros(icalc,no);
samptracks.long_pix        = zeros(icalc,no);
samptracks.lat_clean       = zeros(icalc,no);
samptracks.long_clean      = zeros(icalc,no);
samptracks.lat             = zeros(icalc,no);
samptracks.long            = zeros(icalc,no);
samptracks.time            = result.time;
samptracks.P               = zeros(no,icalc);
samptracks.avgP            = zeros(1,no);
samptracks.L               = zeros(no,icalc);
samptracks.avgL            = zeros(1,no);
samptracks.steps           = zeros(icalc-1,no);
samptracks.length          = zeros(1,no);

dlong = (result.maplong(1,col)-result.maplong(1,1))/(col-1);
dlat  = (result.maplat(row,1)-result.maplat(1,1))/(row-1);
R = mapmatrix(result.maplat(1,1),result.maplong(1,1),dlat, dlong);

% Setup kernels
par1.covmat = 2*s(1)*eye(2);
par2.covmat = 2*s(2)*eye(2);
kern1 = makekern2(par1);
kern2 = makekern2(par2);
ks1 = ceil(max(size(kern1))/2);
ks2 = ceil(max(size(kern2))/2);

% Loop over number of tracks to be generated
for k = 1:no
    % Find terminal (recapture) position
    distr = result.phi(:,:,icalc);
    distr = normalise(distr);
    cdf   = cumsum(distr(:));    %Make cumulated sum = cdf
    index = sum(cdf<rand)+1;     %Find index of sample, +1 because Matlab index start at 1
    [samptracks.lat(icalc,k),samptracks.long(icalc,k)] = ind2sub([row,col],index);

    for i = icalc:-1:2
        % Sample entire track
        if     result.behav(i-1) == 1
            ks = ks1; kern = kern1;
        elseif result.behav(i-1) == 2
            ks = ks2; kern = kern2;
        end
        distr = zeros(ks,ks);
        % Instead of convoluting a dirac delta do something much more
        % complicated, but also faster
        % Construct the correct kernel and copy it to the distribution
        kminlat  = 1+max([ceil(ks/2)-samptracks.lat(i,k) 0]);
        kmaxlat  = min([ks-(samptracks.lat(i,k)+floor(ks/2)-row) ks]);
        kminlong = 1+max([ceil(ks/2)-samptracks.long(i,k) 0]);
        kmaxlong = min([ks-(samptracks.long(i,k)+floor(ks/2)-col) ks]);
        
        mminlat  = max([samptracks.lat(i,k)-floor(ks/2) 1]);
        mmaxlat  = min([samptracks.lat(i,k)+floor(ks/2) row]);
        mminlong = max([samptracks.long(i,k)-floor(ks/2) 1]);
        mmaxlong = min([samptracks.long(i,k)+floor(ks/2) col]);

        distr = result.phi(mminlat:mmaxlat,mminlong:mmaxlong,i-1) .* kern(kminlat:kmaxlat,kminlong:kmaxlong);
        [drow dcol] = size(distr);

        distr = normalise(distr);
        cdf   = cumsum(distr(:));    %Make cumulated sum = cdf
        index = sum(cdf<rand)+1;     %Find index of sample, +1 because Matlab index start at 1
        [lat long] = ind2sub([drow,dcol],index);
        samptracks.lat(i-1,k)  = lat + mminlat -1;
        samptracks.long(i-1,k) = long + mminlong -1;
    end
    % Defining initial position
    %samptracks.xy = [xsamp ysamp];
    samptracks.long_pix_clean(:,k) = samptracks.long(:,k); samptracks.long_pix(:,k) = samptracks.long(:,k);
    samptracks.lat_pix_clean(:,k)  = samptracks.lat(:,k);  samptracks.lat_pix(:,k)  = samptracks.lat(:,k);
    % Make path more clear by adding a random number btw -.5 and .5
    samptracks.long_pix(2:end,k) = samptracks.long_pix_clean(2:end,k) + rand(icalc-1,1)-0.5;
    samptracks.lat_pix(2:end,k)  = samptracks.lat_pix_clean(2:end,k)  + rand(icalc-1,1)-0.5;
    [samptracks.lat_clean(:,k) samptracks.long_clean(:,k)] = pixtomap(R,samptracks.long_pix_clean(:,k),samptracks.lat_pix_clean(:,k));
    [samptracks.lat(:,k) samptracks.long(:,k)] = pixtomap(R,samptracks.long_pix(:,k),samptracks.lat_pix(:,k));
    % Find probability of track %
    for i = 1:icalc
        samptracks.P(k,i) = result.smooth(samptracks.lat_pix_clean(i,k),samptracks.long_pix_clean(i,k),i);
    end
    samptracks.avgP(k) = mean(samptracks.P(k,:));
    % Find step sizes %
    %samptracks.steps2(:,k) = (db.h * sqrt(diff(samptracks.long_pix_clean(:,k)).^2 + diff(samptracks.lat_pix_clean(:,k)).^2));
    %samptracks.length(k)  = sum(samptracks.steps2(:,k));
    
    lats = samptracks.lat_clean(1:end-1,k)+diff(samptracks.lat_clean(:,k));
    samptracks.steps(:,k) = sqrt( (diff(samptracks.lat_clean(:,k)).*deglong(0)).^2 ...
    + (diff(samptracks.long_clean(:,k)).*deglong(lats)).^2);
    samptracks.length(k)  = sum(samptracks.steps(:,k));
    if no < 11, disp(sprintf('Done track %i',k)), end
end

samptracks = proboftrack(samptracks,result,LIK);