function [tag] = process_staroddi(tag)
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
depth_cutoff = 12.;  %time series doesn't start until depth exceeds this value
                     %time series ends last time depth is less than this value
                     %this is used to trim data from the tag where the fish 
                     %is not in the water
if tag.fish_id==71
    depth_cutoff = 20.;
end
% make sure the tagfile exists
if(~exist(tag.datafile));
  fprintf('tag file does not exist %s\n',tag.datafile);
  error('stopping')
end;

% read the data from the tag
nheader = 14; 
if(tag.fish_id > 46); nheader = 15; end;

%[cnt,time1,time2,tchar,depth] = textread(tag.datafile,'%d %s %s %s %f','headerlines',nheader);
[cnt,time1,time2,tchar,dchar] = textread(tag.datafile,'%d %s %s %s %s','headerlines',nheader);

% set arguments
tag.process_date = datestr(now);
  

% convert time to Matlab/GMT
ntimes = numel(cnt);

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
dnum   = datenum([char(time1) char(time2)],'dd.mm.yyHH:MM:SS') + time_shift_hrs/24.;

% process the depth data
depth = -999*ones(ntimes,1);
for i=1:ntimes
  if(~strcmp(char(dchar(i)),'____'))
    depth(i) = str2num(char(dchar(i)));
  end;
end;

% process the temperature data
temp = -999*ones(ntimes,1);
for i=1:ntimes
  if(~strcmp(char(tchar(i)),'____'))
    temp(i) = str2num(char(tchar(i)));
  end;
end;

tag.temp_raw = temp;
% spline temp values to all times
okpts = find(temp > -999);
temp = interp1(dnum(okpts),temp(okpts),dnum,'spline');
    
% trim the tag to have data from only when fish is in the water
pts = find(depth(dnum<tag.recapture_dnum+1) > depth_cutoff);
pts = min(pts):max(pts);
pts = setdiff(pts, find(depth==-999));
tag.dnum_raw = dnum;
tag.depth_raw = depth;

tag.dnum = dnum(pts);
tag.temp = temp(pts);
tag.temp_raw = tag.temp_raw(pts);
tag.depth = depth(pts);
tag.days_at_large = tag.dnum(end)-tag.dnum(1);

% check the time interval 
tag.min_intvl_seconds = (min(diff(tag.dnum))*3600*24);
tag.max_intvl_seconds = (max(diff(tag.dnum))*3600*24);

% plot a quick figure to check against a tide chart 
%figure
%npts = ceil((24*3600)/tag.min_intvl_seconds);
%plot((tag.dnum(1:npts)-floor(tag.dnum(1)))*24-4,tag.depth(1:npts));
%grid on
%error('stop');



