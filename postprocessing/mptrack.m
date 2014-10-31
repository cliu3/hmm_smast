function [mpt] = mptrack(tagno)
%MPTRACK  Find the Most Probable Track.
%   [MPT] = MPTRACK(TAGNO)
%
%   - TAGNO indentifier as string for the tag.
%
%     Output
%
%   - MPT a struct containing the coordinates for the mpt.
%
%   the function assumes the following files are available
%
%     tagdataTAGNO.mat
%     datalikelihoodTAGNO.mat
%     resultTAGNO.mat
%     (tidaldb.mat)
%
%   See the reference manual for reference.
%
%  EXAMPLE   
%   mpt = MPTRACK('2255');
%
%   Date: 4/12 - 2008, ver. 0.58
%   HMM geolocation toolbox, DTU Informatics and DTU Aqua

disp(sprintf('\n\n=== Finding MPT for tag #%s ===',tagno))
warning('off')
filename = ['tagdata' tagno '.mat'];
disp(sprintf('\n\nLoading %s...',filename))
load(filename)
filename = ['datalikelihood' tagno '.mat'];
disp(sprintf('Loading %s...',filename))
load(filename)
if exist('L'), LIK = L; clear L; end
filename = ['result' tagno '.mat'];
disp(sprintf('Loading %s...\n',filename))
load(filename)

if ~isfield(td,'DBname')
    td.DBname = 'tidaldb.mat';
end
disp(['Loading DB:' td.DBname])
load(td.DBname)

disp(sprintf('Number of days: %i\n',length(td.d24)))

if length(result.D) == 1
    result.D(2) = result.D(1);
end
s = result.D*result.D2s;
[row,col,icalc]=size(result.phi);
names = fieldnames(LIK);
names = names(~strcmp(names,'type'));
names = names(~strcmp(names,'mode'));
numnames = length(names);

par1.covmat = 2*s(1)*eye(2);
par2.covmat = 2*s(2)*eye(2);


[~,x_rec]=min(abs(db.long(1,:)-td.catch_long));
[~,y_rec]=min(abs(db.lat(:,1)-td.catch_lat));

% Setup output struct
mpt.land    = db.land;
mpt.maplat  = db.lat;
mpt.maplong = db.long;
mpt.tagno   = td.tagno;
mpt.lat_pix_clean   = zeros(icalc,1);
mpt.long_pix_clean  = zeros(icalc,1);
mpt.lat_pix         = zeros(icalc,1);
mpt.long_pix        = zeros(icalc,1);
mpt.lat_clean       = zeros(icalc,1);
mpt.long_clean      = zeros(icalc,1);
mpt.lat             = zeros(icalc,1);
mpt.long            = zeros(icalc,1);


dlong = (result.maplong(1,col)-result.maplong(1,1))/(col-1);
dlat  = (result.maplat(row,1)-result.maplat(1,1))/(row-1);
R = mapmatrix(result.maplat(1,1),result.maplong(1,1),dlat, dlong);

% Setup state metric
M = log(result.phi(:,:,1)); % initialise
% Find relevant position (ones with finite probability)
subject = (M~=-inf);

zro = zeros(row,col);
theend = icalc;

% Define track array, containing the mpt leading to each pos.
Tprevx = zro; Tprevy = zro;
Tprevx(td.y0,td.x0) = td.x0; Tprevy(td.y0,td.x0) = td.y0;

Ltotal = ones(row,col,icalc-1);
for j = 1:numnames
    Ltotal = Ltotal .* LIK.(names{j});
end
clear LIK, clear result

disp('Starting iterations...')
disp(sprintf('Day   1 -  9...')); tic;
% Viterbi algorithm
for j = 2:theend
    if ~mod(j,10), disp(sprintf('\b done! time = %4.3f',toc)); 
                   disp(sprintf('Day %3.0i -%3.0i...',j,j+9)); tic; end
    Mtemp = log(zro); % Mtemp starts with -inf
    Ttempx = -1+zro; Ttempy = Ttempx;
    Tx = zeros(row,col,j); Ty = Tx;
    for x = 1:col
        for y = 1:row
            if subject(y,x)
                
                par1.u=db.h.*[(x_rec-x) (y_rec-y)]./(2*(1+theend-j));
                par2.u=par1.u;
                
                kern1 = makekern2(par1);
                kern2 = makekern2(par2);
                ks1 = min(size(kern1));
                ks2 = min(size(kern2));
                
                switch td.behav(j-1)
                    case 1
                        ks = ks1; kern = kern1;
                    case 2
                        ks = ks2; kern = kern2;
                end
                
                kminlat  = 1 + max([ceil(ks/2)-y 0]);
                kmaxlat  =     min([ks-(y+floor(ks/2)-row) ks]);
                kminlong = 1 + max([ceil(ks/2)-x 0]);
                kmaxlong =     min([ks-(x+floor(ks/2)-col) ks]);
                klat = kminlat:kmaxlat; klong = kminlong:kmaxlong;
                
                mminlat  = max([y-floor(ks/2) 1]);
                mmaxlat  = min([y+floor(ks/2) row]);
                mminlong = max([x-floor(ks/2) 1]);
                mmaxlong = min([x+floor(ks/2) col]);
                mlat = mminlat:mmaxlat; mlong = mminlong:mmaxlong;

                % Branch matrix for current position
                B = log(Ltotal(mlat,mlong,j-1) .* kern(klat,klong));
                                
                % Total probability of the possible tracks from (x,y)
                Msub = B + M(y,x); % sum of current state and branch metric
                
                % Get relevant area from large array
                Mupdate  = Mtemp(mlat,mlong);
                Txupdate = Ttempx(mlat,mlong);
                Tyupdate = Ttempy(mlat,mlong);
                
                % Find states to be updated
                update = Mupdate<Msub;
                
                % Update in small arrays
                Mupdate(update)  = Msub(update);
                Txupdate(update) = x;
                Tyupdate(update) = y;
                
                % Transfer small arrays to large arrays
                Mtemp(mlat,mlong)  = Mupdate;
                Ttempx(mlat,mlong) = Txupdate;
                Ttempy(mlat,mlong) = Tyupdate;
            end
        end
    end
    Mtemp(db.land) = -inf;
    % Swap tracks
    subject = (Mtemp~=-inf);
    subject(db.land) = 0;
    for x = 1:col
        for y = 1:row
            if subject(y,x)
                Tx(y,x,1:j-1) = Tprevx(Ttempy(y,x),Ttempx(y,x),:);
                Ty(y,x,1:j-1) = Tprevy(Ttempy(y,x),Ttempx(y,x),:);
                Tx(y,x,j) = x; Ty(y,x,j) = y;
            end
        end
    end
    % Update the state metrics
    M = Mtemp;
    %[val ind] = max(M(:));
    %[ym xm]  = ind2sub([row col],ind);
    %LON = shiftdim(Tx(ym,xm,:),2);
    %LAT = shiftdim(Ty(ym,xm,:),2);
    %imagesc(mpt.land), hold on, plot(LON,LAT,'w'), axis ij, drawnow
    % Store current tracks
    Tprevx = Tx; Tprevy = Ty;
    %M(M==-inf)=0;
    %if ~mod(j,10), disp(sprintf('\b done! time = %4.3f',toc)); end
end

%%
M(db.land) = -inf;
[val ind] = max(M(:));
[ym xm]  = ind2sub([row col],ind);
mpt.long_pix_clean = shiftdim(Tx(ym,xm,:),2); mpt.long_pix = mpt.long_pix_clean;
mpt.lat_pix_clean  = shiftdim(Ty(ym,xm,:),2); mpt.lat_pix = mpt.lat_pix_clean;

mpt.long_pix(2:end) = mpt.long_pix_clean(2:end) + rand(theend-1,1)-0.5;
mpt.lat_pix(2:end)  = mpt.lat_pix_clean(2:end)  + rand(theend-1,1)-0.5;
[mpt.lat mpt.long] = pixtomap(R,mpt.long_pix(:),mpt.lat_pix(:));
[mpt.lat_clean mpt.long_clean] = pixtomap(R,mpt.long_pix_clean(:),mpt.lat_pix_clean(:));
mpt.time = td.time_plot(td.d24);

clear Tx, clear Ty, clear Ltotal
load(['result' tagno '.mat']);
load(['datalikelihood' tagno '.mat']);
mpt = proboftrack(mpt,result,LIK);

%mpt.steps  = (db.h * sqrt(diff(mpt.long_pix_clean).^2 + diff(mpt.lat_pix_clean).^2));
%mpt.length = sum(mpt.steps);
lats = mpt.lat_clean(1:end-1)+diff(mpt.lat_clean);
mpt.steps = sqrt( (diff(mpt.lat_clean).*deglong(0)).^2 ...
    + (diff(mpt.long_clean).*deglong(lats)).^2);
mpt.length  = sum(mpt.steps);

% interpolate depth onto the track
mpt.depth = interp2(db.long,db.lat,db.depth,mpt.long_clean,mpt.lat_clean);

%% Creating *.mat file
filename = sprintf('mpt%s',tagno);
disp(sprintf('Saving -> %s.mat <- in\n%s',filename,cd))
save(filename,'mpt');

%% Writing to text file
filename = sprintf('mpt%s.txt',tagno);
disp(sprintf('Saving -> %s <- in\n%s',filename,cd))
fid = fopen(filename,'wt');
fprintf(fid,'%s\n','UTCdate     long       lat');
for i=1:numel(mpt.long);
  fprintf(fid,'%s %12.8f %12.8f\n',datestr(mpt.time(i),'mm/dd/yy HH:MM:SS'), mpt.long(i),mpt.lat(i)); 
end;
fclose(fid);

disp(sprintf('\nDONE mptracking!\n'))
warning('on');
