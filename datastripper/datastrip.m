function YTdatastrip(tagdata,dt)

%DATASTRIP  Strip data file and save in *.mat file.
%   DATASTRIP(TAGNO,FILETS,FILERELCATCH,UNIT,TYP)
%
%   - TAGNO as string eg. TAGNO = '2255' => raw2255.mat.
%   - FILETS pointer to file with time and depth data
%   as string eg. FILETS = '2255PRES.CSV'.
%   - FILERELCATCH pointer to file with release and recapture data
%   as string eg. FILERELCATCH = 'relcatch2255.CSV'.
%
%   -- Units --
%   UNIT = 'm'    - meter,                  1 m = 1 m.
%   UNIT = 'dbar' - decabar,                1 m = 1.0194 dbar.
%   UNIT = 'psi'  - pounds per square inch, 1 m = 0.7028 psi.
%
%   -- Data types --
%   TYP = '1' eg. 2001/03/30 00:01:00,9.698552   (#2255)
%   TYP = '2' eg. 24/03/99 00:01,-0.110853       (#1432)
%   TYP = '3' eg. 06/10/1999 17:21,1.835654      (#2324)
%   TYP = '4' eg. 21.11.2003 00:01,-0.146        (#6448, #872, #897, #6433)
%   TYP = '5' eg. 11.03.2005 00:01,-1.15,6.152   (#1186)
%   TYP = '6' eg. 04/12/05 00:01:00 16.326       (#1950)
%   TYP = '7' eg. 16.03.2003  16:10:00,57.2263   (#3446)
%   TYP = '8' eg. 06/02/2003  13:50	43.57 7.594  (NNS_65)
%   TYP = '9' eg. 38113.68449074  20.11   -1.84  (Excel date since 1900, T(C) PSI) - YT
%
%  EXAMPLE   
%   DATASTRIP('2255','2255PRES.CSV','2255relcatch.CSV','dbar','1')
%   results in the file raw2255.mat being stored in the current dir.
%
%   Date: 6/12 - 2007, ver. 0.51
%   HMM geolocation toolbox, DTU Informatics and DTU Aqua
%
%

if nargin < 2, dt = tagdata.timestep_max; end

td.dt = dt; 
td.tagno = tagdata.tag_ID; 
td.unit = 'm'; 
td.time_org = tagdata.time; 
td.depth_org = tagdata.depth; 
td.temp_org =tagdata.temp ; 
td.time = []; td.depth = []; td.temp = [];
disp(sprintf('\n\n=== Commencing datastrip for tag #%s ===',td.tagno))
%disp(sprintf('TYPE = %s, UNIT = %s\n',typ,unit))

% Open release and recapture file
% dat = fopen(filerelcatch,'r');
% if dat==-1, error('No file named %s found!',filerelcatch), end
% disp(sprintf('Accessing %s and searching for data...',filerelcatch))
% 
% lin = fgetl(dat);
% b = fscanf(dat,'%f,%f,%f,%f,%f',inf);
td.rel_long   = tagdata.rel_lon;
td.rel_lat    = tagdata.rel_lat;
td.catch_long = tagdata.rec_lon;
td.catch_lat  = tagdata.rec_lat;
td.catch_unc  = tagdata.tow_length; %recapture location uncertainty in km?  
td.flag_recpos=tagdata.flag_recpos;
% Open depth data file
% dat = fopen(filets,'r');
% if dat==-1, error('No file named %s found!',filets), end
% disp(sprintf('Accessing %s and searching for data...',filets))
% 
% % Make definitions based on type
% switch typ
%     case '1' % 2001/03/30 00:01:00,9.698552 (#2255)
%         form = '%d/%d/%d %d:%d:%d,%f\n';
%         stringlength = 7;
%     case '2' % 24/03/99 00:01,-0.110853 (#1432)
%         form = '%d/%d/%d %d:%d,%f\n';
%         stringlength = 6;
%     case '3' % 06/10/1999 17:21,1.835654 (#2324)
%         form = '%d/%d/%d %d:%d,%f\n';
%         stringlength = 6;
%     case '4' % 21.11.2003 00:01,-0.146 (#872)
%         form = '%d.%d.%d %d:%d,%f\n';
%         stringlength = 6;
%     case '5' % 11.03.2005 00:01,-1.15,6.152 (#1186)
%         form = '%d.%d.%d %d:%d,%f,%f\n';
%         stringlength = 7;
%     case '6' % 04/12/05 00:01:00	16.326 (#1950)
%         form = '%d/%d/%d %d:%d:%d %f\n';
%         stringlength = 7;
%     case '7' % 16.03.2003  16:10:00,57.22638 (#3446)
%         form = '%d.%d.%d %d:%d:%d,%f\n';
%         stringlength = 7;
%     case '8' % 06/02/2003	13:50	43.57	7.594 (NNS_65)
%         form = '%d/%d/%d %d:%d %f %f\n'
%         stringlength = 7;
%     case '9' % 38113.68449074  20.11   -1.84 
%         form = '%f %f %f\n';
%         stringlength = 3;
%     otherwise
%         error('No type stated!')
% end
% 
% % Find data and save in a
% disp('Finding and arranging data in file...')
% a = []; lin = 1; m=0;
% while lin ~= -1
%     n=0;
%     a_temp = fscanf(dat,form,inf);
%     while length(a_temp) < stringlength*2
%         lin=fgetl(dat);
%         a_temp = fscanf(dat,form,inf);
%         n=n+1;
%         if n>2000, break, end
%         if lin == -1, break, end
%     end
%     %a_temp(end-20:end)
%     a=[a; a_temp(1:end-stringlength)];  %gwc, why cut the last data point
%     m=m+1;
%     if m>2000, break, end
%     if m>2, warning('Data file may have missing data points!'), end
% end
% fclose(dat);
% 
% 
% % Reshape to matrix
% a=reshape(a,stringlength,length(a)/stringlength);
% 
% 
% % Create vectors of time, depth and temperature if available
% switch typ
%     case '1' % (#2255)
%         td.time_org  = datenum(a(1,:), a(2,:),  a(3,:), a(4,:), a(5,:), a(6,:));
%         td.depth_org = a(7,:);
%     case '2' % (#1432)
%         td.time_org  = datenum(a(3,:)+1900, a(2,:),  a(1,:), a(4,:), a(5,:), 0);
%         td.depth_org = a(6,:);
%     case '3' % (#2324)
%         td.time_org  = datenum(a(3,:), a(2,:),  a(1,:), a(4,:), a(5,:), 0);
%         td.depth_org = a(6,:);
%     case '4' % (#872)
%         td.time_org  = datenum(a(3,:), a(2,:),  a(1,:), a(4,:), a(5,:), 0);
%         td.depth_org = a(6,:);
%     case '5' % (#1186)
%         td.time_org  = datenum(a(3,:), a(2,:),  a(1,:), a(4,:), a(5,:), 0);
%         td.depth_org = a(6,:);
%         td.temp_org =   a(7,:);
%     case '6' % (#1950)
%         td.time_org  = datenum(a(3,:)+2000, a(2,:), a(1,:), a(4,:), a(5,:), a(6,:));
%         td.depth_org = a(7,:);
%     case '7' % (#3446)
%         td.time_org  = datenum(a(3,:), a(2,:), a(1,:), a(4,:), a(5,:), a(6,:));
%         td.depth_org = a(7,:);
%     case '8' 
%         td.time_org  = datenum(a(3,:), a(2,:),  a(1,:), a(4,:), a(5,:), 0);
%         td.depth_org = a(6,:);
%         td.temp_org =   a(7,:);
%     case '9' 
%         %smast-preprocessed tags, note time is shifted from excel time (days since 1900,1,0,0,0,0)
%         %into the datenum time (days since 0000,1,0,0,0,0).  
%         %note, we also subtract 1 day because Excel (Windows version, not Mac version) 
%         %incorrectly believes that 1900 was a leap year
%         %which it was not.  Datenum correctly calculates all the leap years.
%         td.time_org  = a(1,:) - 1.  + datenum(1900,1,0,0,0,0);  %shift from excel days to datenum 
%         td.temp_org  = a(2,:);
%         td.depth_org = a(3,:);
%         td.time_org  = td.time_org + 5/24.; %shift to GMT from EST
%     otherwise
%         error('No type stated!')
% end

% Assure that time starts at a 1 minute
% inittime = datestr(td.time_org(1));
% if length(inittime) == 11, numremove = 1; else
% numremove = mod(11-str2num(inittime(17)),10); end
% if numremove ~= 0
%     td.time_org(1:numremove)  = [];
%     td.depth_org(1:numremove) = [];
% end

% Initialise time and depth records
td.depth = td.depth_org;
td.time = td.time_org;
plot(td.time_org-td.time_org(1))
figure
plot(td.time_org-td.time_org(1),td.depth)
figure
plot(diff(td.time_org)*60*24) 
title('time interval from raw tag time')
figure
datestr(td.time_org(1))
datestr(td.time_org(end))

% gwc - note
% this subsampling is a mess, converting to integers and comparing will still suffer 
% from floating point depending on the format of the raw tag time
% we should consider either resplining the raw data to a common time (why would that be bad)
% or transferring to the common time with an exact threshold (within +/- a few seconds is fine)
  

% Subsample in 'dt' min sample intervals
disp('Subsampling...')
SF = 1/(24*60/dt); % Sample frequency in days
td.time = td.time(1):SF:td.time(end);
save junk td
% integer vectors containing sample freq time and original time
t  = round(td.time_org*1e7);
tn = round(td.time*1e7);
% allocate space for depth and temperature on subsample time vector
td.depth = zeros(1,length(tn));
td.temp = zeros(1,length(tn));
if isempty(td.temp_org), td.temp_org = zeros(1,length(td.depth_org)); end  %no temp data in tag, make fake temp data
j = 1; p=0;

%if SF < max(diff(td.time_org)), disp(sprintf('There might be missing data or selected sample rate, %3.1fmin, is faster than data i.e. data is removed!',dt)), pause(0.01), end
% gwc, deal with precision problem
if (diff(td.time_org) - SF > 1./(24*3600)  ), disp(sprintf('There might be missing data or selected sample rate, %3.1fmin, is faster than data i.e. data is removed!',dt)), pause(0.01), end
% subsampling does not interpolate.  If time from the raw matches exactly (using the integer) subsample time, then it will use the data
% from the tag at that time.  If it finds no time stamp at that time, it will give it no data and give an eror
% thus, the subsample time (dt) must be an integer multiple of the largest time interval from the raw tag
mis = 0;
for i=1:length(tn)
    k=0; f=0;
    while k<dt*(6+1) && f~=1
%        if tn(i) == t(j)
        if (abs(tn(i) - t(j)) < 100) %time stamps are within a second should be OK
            td.depth(i) = td.depth_org(j);
            td.temp(i) = td.temp_org(j);
            f=1;
        elseif tn(i) < t(j) % if data missing
            tn(i)       = []; % remove expected data
            td.time(i)  = []; % remove expected data
            td.depth(i) = []; % remove expected data
            td.temp(i)  = []; % remove expected data
            j=j-1;
            mis = mis + 1;
            if mis == 1000, disp('Program execution might be slow due to missing data, or badly selected sample rate!'),pause(0.01), end
        end
        k=k+1; 
        j=j+1;
    end
    if j>length(td.time_org), break, end
end
disp(sprintf('...subsampling removed %i datapoints.',length(td.time_org)-length(td.time)))

% Calibrate tag
% disp('Calibrating and converting units...')
% switch td.unit
%     case 'm',    td.depth = -1 * td.depth;
%     case 'dbar', td.depth = -1/(0.981) * td.depth;
%     case 'psi',  td.depth = -1/(1.4504) * td.depth;
%     otherwise, error('No unit defined!')
% end

%td.depth = td.depth-td.depth(1); % Make sea level equal to zero, KJ
%td.depth = td.depth-max(td.depth); % Make sea level equal to zero

% Find time of release and recapture
% disp('Finding time of release and recapture and trimming...')
% lengthbefore = length(td.time);
% diffdepth = find(abs(diff(td.depth)) > 0.5); % The first and last time the depth changes by more than 0.5 m are the times of release and recapture respectively
% release   = min(diffdepth);
% recapture = max(diffdepth)+1;
% atliberty = release:recapture;
% td.time   = td.time(atliberty);
% td.depth  = td.depth(atliberty);
% td.temp   = td.temp(atliberty);
% trimmed = lengthbefore-length(td.time);
% disp(sprintf('...trimmed %i datapoints (%1.1f days of data).',trimmed,trimmed*SF)) % SF is sample frequency

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

