function datalikelihood(tagno,type,iter,recap)
%DATALIKELIHOOD  Compute datalikelihood from preprocessed tag data
%   DATALIKELIHOOD(TAGNO,TYPE,ITER,RECAP)
%
%   - TAGNO identify the preprocessed data file from the tidebehavextr
%   function to search for in the current directory, eg. TAGNO = '2255'
%   loads tagdata2255.mat.
%
%     Optional arguments
%
%   - TYPE when set to 'fast' relaxes some variance parameters
%   to decrease computation time (database roughnesses). 'full' 
%   includes all variance parameters.
%   default is 'fast'.
%   - ITER when set to 'off' omits iteration output to the screen
%   default is 'on'.
%   - RECAP when set to 'no' the information from recapture position
%   is omitted.
%   default is to use the recapture position.
%
%   DEPENDENCIES - the function needs access to the following files
%
%     tagdataTAGNO.mat
%     tidaldb.mat
%
%  EXAMPLE
%   DATALIKELIHOOD('2255','fast','off')
%
%   Date: 21/10 - 2008, ver. 0.52
%   HMM geolocation toolbox, DTU Informatics and DTU Aqua

if nargin < 4, recap = 'yes'; end
if nargin < 3, iter = 'on'; end
if nargin < 2, type = 'fast'; end
filename = ['tagdata' tagno '.mat'];
disp(sprintf('\n\nLoading %s...',filename))
load(filename), db=1;
if ~isfield(td,'DBname')
    td.DBname = 'tidaldb.mat';
end
disp(['Loading DB:' td.DBname])
load(td.DBname),
dbdir = which(td.DBname); 
td.dbdir = dbdir;
save([td.dbdir(1:end-length(td.DBname)) td.DBname(1:end-4) '_BAK.mat'],'db');
LDB = length(td.DBname);
if (db.lat(1,1) -db.lat(end,end))  < 0, db = flipdb(db,'lat'); save([td.dbdir(1:end-LDB) td.DBname],'db'); end
if (db.long(1,1)-db.long(end,end)) > 0, db = flipdb(db,'long');save([td.dbdir(1:end-LDB) td.DBname],'db'); end
% dbdir = which(td.DBname); save([dbdir(1:end-length(td.DBname)) td.DBname '_BAK.mat'],'db');
% td.dbdir = dbdir;
% LDB = length(td.DBname);
% if (db.lat(1,1) -db.lat(end,end))  < 0, db = flipdb(db,'lat'); save([td.dbdir(1:end-LDB) td.DBname],'db'); end
% if (db.long(1,1)-db.long(end,end)) > 0, db = flipdb(db,'long');save([td.dbdir(1:end-LDB) td.DBname],'db'); end
% %dbdir = which('tidaldb.mat'); save([dbdir(1:end-11) 'tidaldb_BAK.mat'],'db');
% %if (db.lat(1,1) -db.lat(end,end))  < 0, db = flipdb(db,'lat'); save([dbdir(1:end-11) 'tidaldb.mat'],'db'); end
% %if (db.long(1,1)-db.long(end,end)) > 0, db = flipdb(db,'long');save([dbdir(1:end-11) 'tidaldb.mat'],'db'); end
%load('temptidaldb.mat'), disp('USING TEMPTIDALDB')
disp(sprintf('\n=== Compute datalikelihood for tag #%s ===',td.tagno))
LIK.mode = [td.behavrsq' 1-td.behavrsq'];
disp(sprintf('Computation type: %s',type))
disp(sprintf('Display iterations: %s',iter))

days            = 1:length(td.d24)-1;
[row,col,modes] = size(db.amp);
modes           = 1:modes;
LIK.type        = type;

amplitude = zeros(row,col,modes(end));
argument  = zeros(row,col,modes(end));
freq      = zeros(row,col,modes(end));
for mode=modes
    amplitude(:,:,mode) = td.f(mode) .* db.amp(:,:,mode);
    argument(:,:,mode)  = td.G(mode) - db.phase(:,:,mode);
    freq(:,:,mode)      = ones(row,col) .* db.freq(mode);
end

LIK.tide = zeros(row,col,days(end));

% Load datalikelihood parameters
datalikparam

if strcmp(type,'fast')
    invcovs = zeros(td.tideFL,td.tideFL,days(end));
    sigma_tid = invcovs;
    consts  = ones(days(end));
    disp('Setting up covariance matrices for positions...')
    for k = days
        %whos s_E lambda epsilon s_eta_tid s_e
        sigma_tid(:,:,k) = s_E + epsilon(k)^2 * lambda.^c + s_eta_tid; 
        %sigma_tid(:,:,k) = s_E + 0.4^2 * lambda.^c + s_eta_tid; 
        [cholSigma,pdI] = chol(sigma_tid(:,:,k));
        if pdI == 0
            invcovs(:,:,k) = inv(cholSigma);
            invsqrtdetsigma = inv(prod(diag(cholSigma)));
            consts(k) = sqrt((2*pi)^(-td.tideFL)) * invsqrtdetsigma;
        end
    end
end

%%
disp('Computing observational likelihood matrix...')
L2=zeros(row,col);
for k = days
    tic
    if td.tide(k) % Tidal data found
        intv = td.tideBestfit(k):td.tideBestfit(k)+td.tideFL-1;
        t  = td.time(intv);
        ts = td.depth(intv);
        for i = 1:row
            for j = 1:col
                if ~db.land(i,j)
                    for mode=modes
                        temp(mode,:) = amplitude(i,j,mode) * cos( freq(i,j,mode)*t + argument(i,j,mode) );
                    end
                    % Predicted time series for this grid cell
                    mu = -sum(temp,1)+db.depth(i,j);
                    if strcmp(type,'fast') % "fast" computation
                        LIK.tide(i,j,k) = gausspdf(ts,mu,invcovs(:,:,k),consts(k));
                    else
                        sigma = s_E + epsilon(k)^2 * lambda.^c + s_e(i,j)*cospattern + s_eta(i,j);
                        LIK.tide(i,j,k) = mvnpdf(ts,mu,sigma);
                    end
                end
            end
        end
    else % Tidal data not found
        [mindepth indx] = min(td.depth(td.d24(k):td.d24(k+1)-1));
        tidal = amplitude .* cos( freq.*td.time(td.d24(k)-1+indx) + argument );
        % Compute tidal contribution
        mu_depth = -sum(tidal,3) + db.depth;
        LIK.tide(:,:,k)  = normcdf(mindepth*ones(row,col),mu_depth,sqrt(s_eta))...
                     ./(eps+normcdf(zeros(row,col),mu_depth,sqrt(s_eta)));
    end
    endtime=toc;
    if ~strcmp(iter,'off'),disp(sprintf('Done day %i of %i in %3.2f sec',k,days(end),endtime)), end
end
LIK.tide(isnan(LIK.tide)) = 0;

if sum(isnan(LIK.tide(:))) ~= 0, warning('NaN found in LIK.tide!'), end

row
col
td
%% Add recapture position
if ~strcmp(recap,'no')
    distr  = zeros(row,col); distr(td.y1,td.x1) = 1;
    unc    = td.catch_unc/db.h; 
%     ks = ceil(unc*10+1); ks = ks + mod(ks,2) + 1;
%     ksize  = max([15 ks]);
%     kern   = gausskern(ksize,unc);
    par.covmat = unc^2 * eye(2);
    kern   = makekern2(par);
    Lcatch = convn(distr,kern,'same'); %imagesc(Lcatch);
    LIK.tide(:,:,end) = LIK.tide(:,:,end) .* Lcatch;
end

%% Creating *.mat file
filename = sprintf('datalikelihood%s',td.tagno);
disp(sprintf('Saving -> %s.mat <- in\n%s',filename,cd))
save(filename,'LIK')
disp(sprintf('\nDONE with datalikelihood! \n\nNow run --> hmmgeolocate \n\nto create/update the geolocation!\n'))
