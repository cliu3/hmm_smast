function [tag] = process_lotek(tag)
% Read raw tag data and generate a Matlab file containing all information 
%
% function [] = process_tag(raw_tagfile,metadata)  
%
% DESCRIPTION:
%    Read raw tag data and dump to matlab file for tag: tagid 
%
% INPUT 
%   raw_tagile  = textfile containing raw data  
%   metadata    = matlab file containing tag metadata
%
% OUTPUT:
%   matlab file containing tag data + metadata in standardized format 
%   where time series are standardized to GMT.
%
% EXAMPLE USAGE
%    process_lotek_tag('S10440.DAT','S10440_meta'); 
%
% Author(s):  
%    Geoff Cowles (University of Massachusetts Dartmouth)
%
% Revision history
%   
%==============================================================================

% used for testing
%clear all; close all;
%raw_tagfile = 'S10440.DAT';
%metadata = 'S10440_meta';

% set parameters
depth_cutoff = 10.;  %time series doesn't start until depth exceeds this value
                     %time series ends last time depth is less than this value
                     %this is used to trim data from the tag where the fish 
                     %is not in the water

% make sure the tagfile exists
if(~exist(tag.datafile));
  fprintf('tag file does not exist %s\n',tag.datafile);
  error('stopping')
end;

% read the data from the tag
nheader = 5; 

[ttime,temp,depth] = textread(tag.datafile,'%f %f %f','headerlines',nheader);

% set arguments
tag.process_date = datestr(now);

% convert from psi to depth
depth = 1/(1.4504) * depth;

% convert time to Matlab/GMT
ntimes = numel(ttime);

if(strcmp(tag.tzone,'UTC'));
  time_shift_hrs = 0.0;
elseif(strcmp(tag.tzone,'EDT'));
  time_shift_hrs = 4.0;
elseif(strcmp(tag.tzone,'EST'));
  time_shift_hrs = 5.0;
else
  fprintf('tag time zone is %s\n',tag.tzone);
  error('not setup to shift from that time zone');
end;

%smast-preprocessed tags, note time is shifted from excel time (days since 1900,1,0,0,0,0)
%into the datenum time (days since 0000,1,0,0,0,0).  
%note, we also subtract 1 day because Excel (Windows version, not Mac version) 
%incorrectly believes that 1900 was a leap year
%which it was not.  Datenum correctly calculates all the leap years.
dnum   = (ttime -1.0 + datenum(1900,1,0,0,0,0))  + time_shift_hrs/24.;
 

fprintf('tag turned on at %s UTC\n',datestr(dnum(1)));
fprintf('\n');


    
% trim the tag to have data from only when fish is in the water
pts = find(depth > depth_cutoff);
pts = min(pts):max(pts);
tag.dnum_raw = dnum;
tag.temp_raw = temp;
tag.depth_raw = depth;

tag.dnum = dnum(pts);
tag.temp = temp(pts);
tag.depth = depth(pts);
tag.days_at_large = tag.dnum(end)-tag.dnum(1);


fprintf('fish in water at %s UTC\n',datestr(tag.dnum(pts(1))));
% check the time interval 
tag.min_intvl_seconds = (min(diff(tag.dnum))*3600*24);
tag.max_intvl_seconds = (max(diff(tag.dnum))*3600*24);

%figure
%npts = ceil((24*3600)/tag.min_intvl_seconds);
%plot((tag.dnum(1:npts)-floor(tag.dnum(1)))*24-4,tag.depth(1:npts));
%grid on
%error('stop');

