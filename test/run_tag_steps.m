%% Initialization 

clear variables;
close all;
addpath(genpath('../'));


global fvcom_tidaldb % path to fvcom tidal database
fvcom_tidaldb = 'data/fvcomdb_gom3_v2.mat';
global bottom_temperature  % path to fvcom bottom temperature
bottom_temperature   = 'data/gom3_btemp_davged_MayJun_2010.nc';

%ptags = [7, 8];
ptags = 7;

% tag_num_range = ptags;
tag_num = ptags;

global std_temp_offset tag_depth_range tag_depth_accu tag_temp_accu low_fit mod_fit
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% tag-specific paremeters  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
std_temp_offset=2.0; %higher value is more inclusive
tag_depth_range = 250; % in meters
tag_depth_accu = 0.008; % fraction of depth renge
tag_temp_accu = 0.1; % in degree C

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     other paremeters     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
low_fit = 13; %tidal fit duration (hours) for low activity
mod_fit = 5; %tidal fit duration (hours) for moderate activity
D_low = 1; % random walk diffusion coefficient for low activity, in km^2/day
D_high = 10; % random walk diffusion coefficient for high activity, in km^2/day



global tideLV
% tideLV  = [RMSE upper bound, R^2 lower bound, AMPLITUDE lower, AMPLITUDE upper]
tideLV  = [0.42 0.85 0.2 2.0];

% main loop over tags
% for tag_num=tag_num_range


    clear tag;
    clear mpt;
    
    
    
    recap = 'yes';  %yes/no, use/do not use the recapture info
    %if the uncertainty < 0 , this will automatically be disabled
    
    z_off_bottom = 20.0; %max off-bottom extent in meters
    %set to -99.0 to use original Petersen likelihood estimator
    
    

    
    fast_likelihood = 1; %=1 use fast scheme, =0 use more accurate scheme
    
    tagname = [num2str(tag_num) '_raw'];

    %% Part 1
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% generate the tidal database                                    %%%
    %%% only run the command if the raw tidal database (harmonics)     %%%
    %%% have been changed or you wish to change the lon/lat bounds     %%%
    %%% or resolution
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    


        gen_tidaldb_draft(-71,-66,40,45,.05);
        

    
    %% Part 2
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Read in the raw data file from an SMAST-format tag             %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    path_to_tags = 'preprocessing/' ;
    fprintf('loading %s\n',[path_to_tags tagname]);
    if(~exist([path_to_tags tagname '.mat'], 'file'))
        error('tag file does not exist, stopping...');
    end
        load([path_to_tags tagname]);
        
        tagid = [num2str(tag_num) '_' tag.tag_id];
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Subsample the tag data and shift to hmm time                   %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            smast_datastrip(tag)
    %% Part 3
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
        
            tidebehavextr(tagid,mod_fit,tideLV,low_fit,tideLV);  %default
        
        %tidebehavextr(tagid);
        % make strict criteria so no tide is found
        %tidebehavextr(tagid,10,[0.1 0.99 0.2],16,[0.1 0.99 0.2]);
        
        %% Part 4
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Calculate the datalikelihood using the tidal database          %%%
        %%% datalikelihoodxxxx.mat is created                              %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
            % manually override a recap='on' if there is no recap [lon/lat]
            if (tag.recap_uncertainty_km<0)
                recap='no';
                fprintf('============setting recap to NO\n');
            end
            
            if(fast_likelihood)
                datalikelihood(tagid,'fast', 'on',recap,z_off_bottom);
            else
                datalikelihood(tagid,'full', 'on',recap,z_off_bottom);
            end
            
        
        %% Part 5
        if(do_parts(5)==1)
            likelihood_cliu(tag_num,path_to_tags,tagname)
            tidal_rmse_cliu(tag_num,path_to_tags,tagname)
        end
        
        %% Part 6
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Geolocate the tag                                              %%%
        %%% resultxxxx.mat is created                                      %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        

            hmmgeolocate(tagid,2,'on',[D_low, D_high]);
        
        %% Part 7
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Find the most probable track of the fish                       %%%
        %%% mptxxxx.mat is created                                         %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
            mpt=mptrack(tagid);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%% Plot the most probable track                                   %%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            close all
            plottrack(mpt);
            %write_ge_track(mpt);
            write_ge_track_UD(tagid)
        
        %% Part 8
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Create an avi file that shows how the probability distribution %%%
        %%% evolves in time                                                %%%
        %%% geolocationxxxx.avi is created                                 %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
            %write_ge_track_UD(tagid)
            pause on;close all;
            avimaker(tagid)
        

% end


