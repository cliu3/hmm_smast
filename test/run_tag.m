%clear all;
%close all;
addpath(genpath('../'));
%addpath(genpath('../../preprocess/'));
%addpath('/opt/matlab/googleearth');

global fvcom_tidaldb % path to fvcom tidal database
fvcom_tidaldb = 'data/fvcomdb_gom3_v2.mat';
global bottom_temperature  % path to fvcom bottom temperature
bottom_temperature   = 'data/gom3_btemp_davged_MayJun_2010.nc';

%ptags = [7, 8];
ptags = 7;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% tag-specific paremeters  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global std_temp_offset tag_depth_range tag_depth_accu tag_temp_accu low_fit mod_fit
std_temp_offset=2.0; %higher value is more inclusive
tag_depth_range = 250; % in meters
tag_depth_accu = 0.008; % fraction of depth renge
tag_temp_accu = 0.1; % in degree C
low_fit = 13; %tidal fit duration (hours) for low activity
mod_fit = 5; %tidal fit duration (hours) for moderate activity
D_low = 1; % random walk diffusioncoefficient for low activity, in km^2/day
D_high = 10; % random walk diffusioncoefficient for high activity, in km^2/day

tag_num_range = ptags;

global tideLV
% tideLV  = [RMSE upper bound, R^2 lower bound, AMPLITUDE lower, AMPLITUDE upper]
tideLV  = [0.42 0.85 0.2 2.0];

% main loop over tags
for tag_num=tag_num_range
    clear tag;
    clear mpt;
    
    
    
    recap = 'yes';  %yes/no, use/do not use the recapture info
    %if the uncertainty < 0 , this will automatically be disabled
    
    z_off_bottom = 20.0; %max off-bottom extent in meters
    %set to -99.0 to use original Petersen likelihood estimator
    
    %do_parts = 6;
    
    do_parts(1) = 1; %1 generate a new tidaldb, =0 use tidaldb.mat
    do_parts(2) = 1; %2 strip
    do_parts(3) = 1; %3 behavior
    do_parts(4) = 1; %4 likelihood
    do_parts(5) = 1; %5 cliu likelihood & tidal constraint
    do_parts(6) = 1; %6 geolocate
    do_parts(7) = 1; %7 most probable track
    do_parts(8) = 0; %8 make a movie
    
    fast_likelihood = 1; %=1 use fast scheme, =0 use more accurate scheme
    
    tagname = [num2str(tag_num) '_raw'];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% generate the tidal database                                    %%%
    %%% only run the command if the raw tidal database (harmonics)     %%%
    %%% have been changed or you wish to change the lon/lat bounds     %%%
    %%% or resolution
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if(do_parts(1)==1);
        %  gen_tidaldb(-71,-70,42,43,.02,{'M2','S2','N2','K1','O1'}); %for fixed tags
        %gen_tidaldb(-71,-66,40,45,.015,{'M2','N2','S2','O1','K1','K2','P1','Q1'}); %all GOM
        %gen_tidaldb(-71,-66,40,45,.05,{'M2','N2','S2','O1','K1','K2','P1','Q1'}); %all GOM
        %readdb
        %finddbvars
        gen_tidaldb_draft(-71,-66,40,45,.05);
        
    end;
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Read in the raw data file from an SMAST-format tag             %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    path_to_tags = 'preprocessing/' ;
    fprintf('loading %s\n',[path_to_tags tagname]);
    if(exist([path_to_tags tagname '.mat']))
        
        load([path_to_tags tagname]);
        
        tagid = [num2str(tag_num) '_' tag.tag_id];
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Subsample the tag data and shift to hmm time                   %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if(do_parts(2)==1)
            smast_datastrip(tag)
        end;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Extract tidal and behaviour information from rawXXXX.mat       %%%
        %%% tagdataxxxx.mat is created                                     %%%
        % tidebehavextr(tagno,tideFL,tideLV,behavFL,behavLV,DBname)
        %
        %   inputs:
        %      tidFL is the length of the fitting interval in hours (def=10)
        %      tideLV - goodness of fit thresholds for detecting tide
        %            [RMSE, R^2, AMPLITUDE]:  def = [0.42 0.85 0.6]
        %              fit is found if criteria are met:
        %                  rmse_computed < RMSE
        %                  r^2_computed  > R^2
        %                  amp_computed  > AMPLITUDE
        %      behavFL lenght of fitting interval for behaviour - 16
        %      behavLV - goodness of fit thresholds for detecting behavior
        %            [RMSE, R^2, AMPLITUDE]:  def = [0.42 0.85 0.6]
        %      dbname - alternate tidal base: def = tidaldb.mat
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if(do_parts(3)==1)
            tidebehavextr(tagid,mod_fit,tideLV,low_fit,tideLV);  %default
        end;
        %tidebehavextr(tagid);
        % make strict criteria so no tide is found
        %tidebehavextr(tagid,10,[0.1 0.99 0.2],16,[0.1 0.99 0.2]);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Calculate the datalikelihood using the tidal database          %%%
        %%% datalikelihoodxxxx.mat is created                              %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if(do_parts(4)==1)
            % manually override a recap='on' if there is no recap [lon/lat]
            if (tag.recap_uncertainty_km<0)
                recap='no';
                fprintf('============setting recap to NO\n');
            end;
            
            if(fast_likelihood)
                datalikelihood(tagid,'fast', 'on',recap,z_off_bottom);
            else
                datalikelihood(tagid,'full', 'on',recap,z_off_bottom);
            end;
            
        end;
        
        if(do_parts(5)==1)
            likelihood_cliu(tag_num,path_to_tags,tagname)
            tidal_rmse_cliu(tag_num,path_to_tags,tagname)
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Geolocate the tag                                              %%%
        %%% resultxxxx.mat is created                                      %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if(do_parts(6)==1)
            %hmmgeolocate(tagid,2,'on',[],true)
            %hmmgeolocate(tagid,2,'on',[10. 100.]);
            hmmgeolocate(tagid,2,'on',[D_low, D_high]);
        end;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Find the most probable track of the fish                       %%%
        %%% mptxxxx.mat is created                                         %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if(do_parts(7)==1)
            mpt=mptrack(tagid);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%% Plot the most probable track                                   %%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            close all
            plottrack(mpt);
            %write_ge_track(mpt);
            write_ge_track_UD(tagid)
        end;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Create an avi file that shows how the probability distribution %%%
        %%% evolves in time                                                %%%
        %%% geolocationxxxx.avi is created                                 %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if(do_parts(8)==1)
            %write_ge_track_UD(tagid)
            pause on;close all;
            avimaker(tagid)
        end
    else
        error('tag file does not exist, stopping...');
    end;
end;
%return

