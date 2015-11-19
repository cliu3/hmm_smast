function hmmgeolocate(tagno,mode,viewres,Duser,ext)
%HMMGEOLOCATE  Obtain geolocation by filtering preprocessed data
%   HMMGEOLOCATE(TAGNO,MODE,VIEWRES,DUSER,EXT)
%
%   - TAGNO indentifier as string for the tag to geolocate.
%
%     Optional arguments
%
%   - MODE number of behaviour modes to use (1 or 2).
%   default is 2.
%   - VIEWRES plots the marginal distributions consecutively
%   in an animation when the geolocation has finished
%   successfully (by using the viewdistr function).
%   default is 'on'.
%   - DUSER user defined diffusivity vector e.g DUSER = [10 100].
%   if omitted the diffusivity is estimated with maximum likelihood.
%   (- EXT self estimate behaviour - IS NOT OPERATIONAL YET!)
%
%   DEPENDENCIES - the function needs access to the following files
%
%     tagdataTAGNO.mat
%     datalikelihoodTAGNO.mat
%     tidaldb.mat
%     cmap.mat
%
%   and creates as output the file resultTAGNO.mat in the current folder.
%
%  EXAMPLES
%   HMMGEOLOCATE('2255',2,'on')
%   HMMGEOLOCATE('1432',[],'on',[10 100])
%
%   Date: 22/10 - 2008, ver. 0.55
%   HMM geolocation toolbox, DTU Informatics and DTU Aqua

disp(sprintf('\n\n=== Geolocating tag #%s ===\n',tagno))
filename = ['tagdata' tagno '.mat'];
disp(sprintf('\n\nLoading %s...',filename))
load(filename)
filename = ['datalikelihood' tagno '.mat'];
disp(sprintf('Loading %s...\n',filename))
load(filename)
if exist('L','var'), LIK = L; clear L; end
if ~isfield(td,'DBname')
    td.DBname = 'tidaldb.mat';
end
disp(['Loading DB:' td.DBname])
load(td.DBname),
%load('tidaldb.mat')
%load temptidaldb, load likttemp2255_10; L.temp = Ltemp;
load cmap

if nargin < 3 | isempty(viewres),  viewres = 'on'; end
if nargin < 2 | isempty(mode),  mode = 2; end
if length(unique(td.behav)) == 1, mode = 1; end
if nargin < 4 | isempty(Duser)
    if mode == 1
        Duser = [60 60];
    else
        Duser = [10  100];
    end
    mp = 0.5; % Start guess for p in (1-p)*mode1 + p*mode2, behaviour switching
    estimate = 1;
else
    estimate = 0;
    if isstruct(Duser)
        mp = Duser.mp;
        Duser = Duser.Duser;
    else
        % It is assumed that no mode probability (mp) is defined
        if length(Duser) == 1
            disp(sprintf('Using user defined diffusivity, one mode:\nD = %8.4f',Duser))
            Duser = [Duser Duser];
        elseif length(Duser) == 2
            disp(sprintf('Using user defined diffusivity, two modes:\nD = [%8.4f, %8.4f]',Duser(1),Duser(2)))
        end
    end
end
if nargin < 5, ext = 0; end

ext = 0; %%% EXTENDED BEHVAIOUR ESTIMATION NOT YET OPERATIONAL!

disp(sprintf('==Tag #%s==',tagno));

% Number of days to calculate forward
[row,col]=size(db.depth);
icalc = length(td.d24);
if mode == 1
    td.behav = ones(1,icalc); disp('Using ONE behavioural mode!')
elseif mode == 2
    disp('Using TWO behavioural modes!')
end
result.behav   = td.behav; % Store behaviour vector for eg. track sampling
result.land    = db.land;  % Store land for plotting of track
result.maplong = db.long; result.maplat = db.lat;
result.time    = td.time_plot(td.d24);
result.tagno   = tagno;
result.DBname  = td.DBname;
result.dbdir   = td.dbdir;

%% Set up parameters %%
% k = time step = 1 day
k=1;

while 1
    %% Initial guess on D %%
    if db.h <= 0, error('db.h is less than or equal to zero, ie. the database i degenerate!'), end
    if ext
        result.D = Duser; % Unit: km^2/d?gn % 2255 18/2
        
    else
        result.D = Duser; % Unit: km^2/d?gn % 2255 18/2
        D2s  = k/db.h^2; result.D2s = D2s;
        sEst = result.D*D2s;
        s = kstest(sEst,25); % Calc size of conv kernel, stop if large
        result.D = s/D2s;
    end
    % [s I] = max(result.D)*D2s;
    % unc    = sqrt(2*s);
    % ks = ceil(unc*10+1); ks = ks + mod(ks,2) + 1;
    % if ks > 100,
    %     fac = (25/ks)^2;
    %     result.D = [result.D(I) result.D(I)]*fac;
    %     sEst = result.D*D2s;
    %     warning('Diffusivity too large, overwriting! new D= %f, %f',result.D(1),result.D(2))
    % end
    % s = max(result.D(unique(td.behav))*D2s);
    unc    = sqrt(2*s);
    ks = ceil(unc*10+1);
    ks = ks + mod(ks,2) + 1;
    
    %% Plot likelihood function
    nevals = 5; % Number of evaluation points on likelihood curve
    %LB = [2 2]*D2s; UB = [225 225]*D2s;
    LB = [2 2]*D2s; UB = [300 300]*D2s;
    
    %LB = [.2443 .2443]/10000*D2s; UB = [32.8 32.8]/10000*D2s;
    %LB = kstest(LB,7); UB = kstest(UB,81);
    
    %plotlikelihoodfunction
    %return
    
    %% Find MLE using builtin fmincon (optimisation toolbox)
    if estimate == 1
        disp('Estimating D...'), tic
        ds = 0.00001;
        
        if mode == 1
            [sEst,loglik] = fminbnd(@(s) likelihood(s,db,td,LIK),LB(1),UB(1), ...
                optimset('TolX',1e-4,'MaxFunEvals',30,'Display','iter'));
            sEst = [sEst sEst];
            time_fmins = toc; result.D  = sEst/(D2s);
            hess=(likelihood(sEst(1)+ds,db,td,LIK)+likelihood(sEst(1)-ds,db,td,LIK)-2*loglik)/ds^2;
            result.MLvar = [1/(hess*D2s^2) 1/(hess*D2s^2)];
            disp(sprintf('\nMLE found!\nDhat = %1.4f,\t stdev = %f\ntime spent: %1.4f sec\n',result.D(1), sqrt(result.MLvar(1)),time_fmins))
        else
            guess = result.D*D2s;
            [sEst,loglik] = fminsearchbnd(@(s) likelihood(s,db,td,LIK),guess,LB,UB, ...
                optimset('TolX',1e-4,'MaxFunEvals',30,'Display','iter'));
            time_fmins = toc; result.D  = sEst/(D2s);
            hess1=(likelihood(sEst+[ds 0],db,td,LIK)+likelihood(sEst+[-ds 0],db,td,LIK)-2*loglik)/ds^2;
            hess2=(likelihood(sEst+[0 ds],db,td,LIK)+likelihood(sEst+[0 -ds],db,td,LIK)-2*loglik)/ds^2;
            result.MLvar(1) = 1/(hess1*D2s^2);
            result.MLvar(2) = 1/(hess2*D2s^2);
            disp(sprintf('\nMLE found!\nDhat1 = %1.4f,\t stdev1 = %f\nDhat2 = %1.4f,\t stdev2 = %f\ntime spent: %1.4f sec\n',result.D(1), sqrt(result.MLvar(1)),result.D(2), sqrt(result.MLvar(2)),time_fmins))
        end
        result.loglikval = loglik;
    end
    
    if result.D(1)>result.D(2)
        fprintf('\nD(1)>D(2), MLE failed... \nUsing default D\n');
        result.D=guess./D2s;
    end
    %% Prediction %%
    disp(sprintf('Using D: [%f %f]',result.D(1),result.D(2))),
    if ext == 0
        disp('Predicting...'),
        %tic, [result.phi,normaliser] = hmmfilter(sEst,db,td,LIK); tt=toc;
        tic, [result.phi,normaliser,result.pred,isDtoosmall] = hmmfilter(sEst,db,td,LIK); tt=toc;
    else
        disp('Predicting (extended version)...')
        par.s = sEst;
        par.mp = mp;
        tic, [result.phi,normaliser,result.pred] = hmmfiltermode(par,td,LIK); tt=toc;
    end
    
    if isDtoosmall==0
        break
    else
        Duser = Duser.*2;
        disp('D is too small, doubled D and trying again...')
    end
end
disp(sprintf('\b done in %3.2f sec!',tt))

%viewdistr(squeeze(result.phi.p(:,:,1,:)))

%% Smoothing %%
disp('Smoothing...')
%tic, [result.smooth] = smoothing(sEst,result.phi,db,td); tt=toc;
tic, [result.smooth] = smoothing(sEst,result.phi,db,td,result.pred); tt=toc;
disp(sprintf('\b done in %3.2f sec!',tt))

result.phi_plot = zeros(row,col,icalc); result.smooth_plot = result.phi_plot;
for i=1:icalc
    post = result.smooth(:,:,i);
    post(db.land) = -0.1*max(post(:));
    result.smooth_plot(:,:,i) = post;
    post = result.phi(:,:,i);
    post(db.land) = -0.1*max(post(:));
    result.phi_plot(:,:,i) = post;
end

%% Utilisation Distribution %%
result.UD = sum(result.smooth,3);
result.UD = normalise(result.UD);
post = result.UD;
post(db.land) = -0.1*max(post(:));
result.UD_plot = post;

%% Creating *.mat file
filename = sprintf('result%s',td.tagno);
disp(sprintf('Saving -> %s.mat <- in\n%s',filename,cd))
%try
%    save(filename,'result')
%catch EM
    save(filename,'result','-v7.3');  %for large bathymetric databases / long tags struct is too large for v7 datatype
%end
disp(sprintf('\nDONE geolocating!\n'))

%% View smoothed distribution
%if ~strcmp(viewres,'off'), close all, viewdistr(result.smooth_plot); end
if ~strcmp(viewres,'off'), close all,
    figure, set(gcf,'position',[50 150 900 400])
    subplot(121), viewdistr(result.UD_plot), title('Utilisation Distribution')
    %subplot(121), viewdistr(result.UD_plot,[],[],[],[],'fancylock',result.land), title('Utilisation Distribution')
    subplot(122), viewdistr(result.smooth_plot);
    %subplot(122), viewdistr(result.smooth,[],[],[],[],'fancylock',result.land);
end

function sEst=kstest(s,siz)
sEst = s;
[s I] = max(s);
unc = sqrt(2*s);
ks = ceil(unc*10+1); ks = ks + mod(ks,2) + 1;
if ks > 100,
    fac = (siz/ks)^2;
    sEst = [s s]*fac
    warning('Diffusivity too large, overwriting!')
end


