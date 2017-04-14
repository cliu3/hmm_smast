function smast_datastrip(tagdata,dt)
%close all; clear all;
%load ../tag_data/S10440;
%tagdata = tag;
% function smast_datastrip
%
% read tag data structure generated by readtag (see readtag.m)
%
% DESCRIPTION:
%  Read a tag data structure created by readtag.m
%  Sample the tag data at a discrete time step
%  Shift the time to yearday
%  Plot some stuff
%
% INPUT
%  tagdata:  path to an SMAST-format tag text file
%  dt:  time step in minutes.  This argument is optional.  If it is not present
%       hmm will use the max timestep from the SMAST-format tag file
%
% OUTPUT:
%    fvcomdb.mat
%
% EXAMPLE USAGE
%    smast_datastrip(tagdata,dt)
%
% Author(s):  
%    Martin Pedersen
%
% Revision History
%    8/2010: Klavdija Jenko (UMASSD):  modified routine to read in tag structure
%                                      generated by readtag
% Revision History
%    5/2012: Geoff Cowles (UMASSD):  modified to read Matlab-based SMAST tags
%==============================================================================

% set the time step (use argument dt if present)
if nargin < 2, dt = tagdata.max_intvl_seconds/60; end; 
%dt = round(dt);
td.deltat = dt;

%-------------------------------------------------------------------
% Set tag timeseries of time and depth, tag number, and time step
%-------------------------------------------------------------------
%td.dt = int(dt);
if isfield(tagdata,'fish_id')
    td.tagno = [num2str(tagdata.fish_id) '_'  tagdata.tag_id];
else
    td.tagno = tagdata.tag_id;
end
td.unit = 'm';
td.time_org = tagdata.dnum;
td.depth_org = -tagdata.depth; %depth of tag record is positive up
td.temp_org =tagdata.temp ;
td.time = []; td.depth = []; td.temp = [];
disp(sprintf('\n\n=== Commencing datastrip for tag #%s ===',td.tagno))

%-------------------------------------------------------------------
% Set release/recapture/recapture uncertainty information
%-------------------------------------------------------------------
td.rel_long   = tagdata.release_lon;
td.rel_lat    = tagdata.release_lat;
td.catch_long = tagdata.recapture_lon;
td.catch_lat  = tagdata.recapture_lat;
td.catch_unc  = tagdata.recap_uncertainty_km;  
% determine recapture availability
if (tagdata.recap_uncertainty_km > 0 && floor(tagdata.dnum(end)) >= floor(tagdata.recapture_dnum) )
    td.recap = 'yes';
else
    td.recap = 'no';
end

% Initialise time and depth records
td.depth = td.depth_org;
td.time = td.time_org;
plot(td.time_org-td.time_org(1))
figure
plot(td.time_org-td.time_org(1),td.depth)
figure
plot(diff(td.time_org)*60*24) 
title('time interval from raw tag time (minutes)')
figure
datestr(td.time_org(1))
datestr(td.time_org(end))

% % gwc - note
% % this subsampling is a mess, converting to integers and comparing will still suffer 
% % from floating point depending on the format of the raw tag time
% % we should consider either resplining the raw data to a common time (why would that be bad)
% % or transferring to the common time with an exact threshold (within +/- a few seconds is fine)
% % fixed - CL 5/26/2016 
% 
% % Subsample in 'dt' min sample intervals
% disp('Subsampling...')
% SF = 1/(24*60/dt); % Sample frequency in days
% td.time = td.time(1):SF:td.time(end);
% save junk td
% % integer vectors containing sample freq time and original time
% t  = round(td.time_org*1e7);
% tn = round(td.time*1e7);
% % allocate space for depth and temperature on subsample time vector
% td.depth = zeros(1,length(tn));
% td.temp = zeros(1,length(tn));
% if isempty(td.temp_org), td.temp_org = zeros(1,length(td.depth_org)); end  %no temp data in tag, make fake temp data
% j = 1; p=0;
% 
% %if SF < max(diff(td.time_org)), disp(sprintf('There might be missing data or selected sample rate, %3.1fmin, is faster than data i.e. data is removed!',dt)), pause(0.01), end
% % gwc, deal with precision problem
% if (diff(td.time_org) - SF > 1./(24*3600)  ), disp(sprintf('There might be missing data or selected sample rate, %3.1fmin, is faster than data i.e. data is removed!',dt)), pause(0.01), end
% % subsampling does not interpolate.  If time from the raw matches exactly (using the integer) subsample time, then it will use the data
% % from the tag at that time.  If it finds no time stamp at that time, it will give it no data and give an eror
% % thus, the subsample time (dt) must be an integer multiple of the largest time interval from the raw tag
% mis = 0;
% for i=1:length(tn)
%     k=0; f=0;
%     while k<dt*(6+1) && f~=1 && i <= numel(tn)
% %        if tn(i) == t(j)
%         if (abs(tn(i) - t(j)) < 1000) %time stamps are within a second should be OK
%             td.depth(i) = td.depth_org(j);
%             td.temp(i) = td.temp_org(j);
%             f=1;
%         elseif tn(i) < t(j) % if data missing
%             tn(i)       = []; % remove expected data
%             td.time(i)  = []; % remove expected data
%             td.depth(i) = []; % remove expected data
%             td.temp(i)  = []; % remove expected data
%             j=j-1;
%             mis = mis + 1;
%             if mis == 1000, disp('Program execution might be slow due to missing data, or badly selected sample rate!'),pause(0.01), end
%         end
%         k=k+1; 
%         j=j+1;
%     end
%     if j>length(td.time_org), break, end
% end
% disp(sprintf('...subsampling removed %i datapoints.',length(td.time_org)-length(td.time)))


% Make td.time_plot and shift td.time
year = str2num(datestr(td.time(1),10)); 
td.time_plot = td.time;
td.time = td.time - datenum(year,1,1,0,0,0);

% Creating *.mat file
filename = sprintf('raw%s',td.tagno);
disp(sprintf('Saving -> %s.mat <- in\n%s',filename,cd))
save(filename,'td')
disp(sprintf('\nDONE with datastrip! \n\nNow run --> tidebehavextr \n\nto create/update tidal and behaviour extraction!\n'))

%% Plot result %%
plot(td.time_org,td.depth_org,'-'), datetick('x','QQ-YY')
xlabel('Date'), ylabel(sprintf('Unit: %s','Depth [m]'))
title(sprintf('Unprocessed depth record for tag #%s',td.tagno))
figure, plot(td.time_plot,td.depth,'-'), datetick('x','QQ-YY')
xlabel('Date'), ylabel('Depth, m')
title(sprintf('Subsampled and trimmed depth record for tag #%s',td.tagno))

