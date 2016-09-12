clear all;close all;
addpath(genpath('../../dependencies/t_tide/')); 
% process all the tags 

% raw tag directory
tagdir = 'raw_tag_data';

% remove all the existing tag files
% system('rm *.mat');

% load the project summary info from Doug's spreadsheet 
%fname = 'Inventory of SCCZ DST Recaptures.xlsx';
fname = 'Inventory_of_tags.xlsx';
[numeric,txt,raw]= xlsread(fname);
[ntags,ncols] = size(numeric);


ptags=[7, 8];

tagset = ptags;

% loop over tags
for i=tagset  
  clear tag;

  tag.datafile  = [tagdir '/' strtrim(char(txt(i+1,14)))];
  if(strcmp(tag.datafile,[tagdir '/NONE'])); 
    continue
  end;

  % set project metadata
  common_meta;

  % set tag-specific metadata 
  tag.fish_id   = numeric(i,1);
  tag.tag_id    = char(txt(i+1,2));
  tag.type      = char(txt(i+1,3));
  tag.length    = numeric(i,10); 
  tag.sex       = char(txt(i+1,11));
  tag.maturity  = numeric(i,12);
  tag.release_dnum = (numeric(i,8) + datenum(1900,1,0,0,0,0));
  tag.recapture_lon = numeric(i,19); 
  tag.recapture_lat = numeric(i,18); 
  tag.recapture_dnum = (numeric(i,15) + datenum(1900,1,0,0,0,0));
  tag.recap_uncertainty_km = numeric(i,21); 


  % for stationaey tags release location = recap location
  if ismember(i,63:89)
      tag.release_lon = tag.recapture_lon;
      tag.release_lat = tag.recapture_lat;
  end
  
  %set time-dependent fields
  if(strcmp(tag.type,'STAR-ODDI'));
    tag = process_staroddi(tag);
  elseif(strcmp(tag.type,'LOTEK'));
    %fprintf('something up with LOTEK tags, skipping \n');
    %continue
    tag = process_lotek(tag);
  else
    fprintf('type is %s\n',tag.type);
    error('not setup to read the type');
  end;

  %check the tidal harmonics over first few days
  %fname = [num2str(tag.fish_id) '_' tag.tag_id];
  fname = [num2str(tag.fish_id) '_raw'];
  save(fname,'tag');

  check_tag(tag);
end;
