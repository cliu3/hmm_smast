clear all, close all, clc

addpath(genpath('../../preprocess/'));
addpath(genpath('../'));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% generate the tidal database                                    %%%
%%% only run the command if the raw tidal database (harmonics)     %%%
%%% have been changed or you wish to change the lon/lat bounds     %%%
%%% or resolution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

gen_tidaldb(-71,-69,42,44,.05,{'M2','N2','S2','O1','K1','K2','P1','Q1'})

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Read the database text files                                   %%%
%%% rawtidaldb.mat is created                                      %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

readdb

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Read the rawtidaldb.mat and compute database variances         %%%
%%% tidaldb.mat is created                                         %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
finddbvars

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Read in the raw data file from an SMAST-format tag             %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%path_to_data='C:\Documents and Settings\localadmin\My Documents\MATLAB\program for the header\data\';
path_to_data = '~/Dropbox/Geolocation/projects/cod_zemeckis/fixed_sensor_data/';

% datafile = 'D10137est.txt';
% datafile = 'D10079est.txt';
% datafile = 'D08386est.txt';
%datafile = 'D03008est.txt';
% datafile = '9999.txt'; %idealized bathymetry track
datafile='S12068.mat';  %Dougs mooring data
load([path_to_data,datafile],'tag')
tagdata=tag;
clear tag;
%tagdata=readtag(path_to_data,datafile);
tagid=tagdata.tag_id;
tagdata.flag_recpos=1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Subsample the tag data and shift to hmm time                   %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
smast_datastrip(tagdata)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Extract tidal and behaviour information from rawXXXX.mat       %%%
%%% tagdataxxxx.mat is created                                     %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tidebehavextr(tagid); 
% make strict criteria so no tide is found
%tidebehavextr(tagid,10,[0.1 0.99 0.2],16,[0.1 0.99 0.2]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Calculate the datalikelihood using the tidal database          %%%
%%% datalikelihoodxxxx.mat is created                              %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
recap='no';
if isfield(tagdata,'flag_recpos')
    if (tagdata.flag_recpos==1)
        recap='on';
    end;
end;
%recap = 'no';  %option for user to manually override and ignore recap 
datalikelihood(tagid,'fast', 'on',recap);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Geolocate the tag                                              %%%
%%% resultxxxx.mat is created                                      %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%hmmgeolocate(tagid,2,'on');
hmmgeolocate(tagid,2,'on',[1 1]);

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
%return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Create an avi file that shows how the probability distribution %%%
%%% evolves in time                                                %%%
%%% geolocationxxxx.avi is created                                 %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pause on
avimaker(tagid)
