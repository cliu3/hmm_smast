function tidebehavextr(tagno,tideFL,tideLV,behavFL,behavLV,DBname)
%TIDEBEHAVEXTR  Extract tidal and behaviour information from a raw data file.
%   TIDEBEHAVEXTR(TAGNO,TIDEFL,TIDELV,BEHAVFL,BEHAVLV,DBNAME)
%
%   - TAGNO identify the raw data file from the datastrip function to  
%   search for in the current directory, eg. TAGNO = '2255' loads raw2255.mat.
%
%     Optional arguments
%
%   - TIDEFL is the length of the fitting interval in hours for the
%   tidal extraction, eg. TIDEFL = 9.
%   default value is 10.
%   - TIDELV is the limit values that determine whether an interval
%   contains tidal information, TIDELV = [RMSE RSQUARE AMPLITUDE]
%   RMSE is the root mean square error in metres, eg. RMSE = 0.3
%   RSQUARE is the coefficient of determination, eg. RSQUARE = 0.9
%   AMPLITUDE is the amplitude of the fit in metres, eg. RSQUARE = 0.1
%   default values are TIDELV = [0.42 0.85 0.6].
%   - BEHAVFL is the length of the fitting interval in hours for the
%   behaviour extraction, eg. BEHAVFL = 15
%   default value is 16.
%   - BEHAVLV is the limit values that determine whether an interval
%   contains tidal information, BEHAVLV = [RMSE RSQUARE AMPLITUDE]
%   default values are BEHAVLV = [0.42 0.85 0.6].
%   - DBNAME is the file name of an alternative tidal data base
%   default is 'tidaldb.mat'
%
%   DEPENDENCIES - the function needs access to the following files
%
%     rawTAGNO.mat
%     tidaldb.mat
%     lssinfit.m
%
%  EXAMPLE
%   TIDEBEHAVEXTR('2255',10,[0.42 0.85 0.6],16,[0.42 0.85 0.6])
%   results in the file tagdata2255.mat being stored in the current dir
%
%   Date: 28/11 - 2008, ver. 0.57
%   HMM geolocation toolbox, DTU Informatics and DTU Aqua

warning('on')
% Define default values
if nargin < 2 || isempty(tideFL),  tideFL  = 10; end
if nargin < 3 || isempty(tideLV),  tideLV  = [0.42 0.85 0.6]; end
if nargin < 4 || isempty(behavFL), behavFL = 16; end
if nargin < 5 || isempty(behavLV), behavLV = [0.42 0.85 0.6]; end
if nargin < 6 || isempty(DBname), DBname = 'tidaldb.mat'; end
% Error messages
if tideFL <= 0, error('Bad input for tideFL!'), end
if (tideLV(1) <= 0 || (tideLV(2)<=0 || tideLV(2)>=1) || tideLV(3)<=0 ), error('Bad input in tideLV!'), end
if behavFL <= 0, error('Bad input for behavFL!'), end
if (behavLV(1) <= 0 || (behavLV(2)<=0 || behavLV(2)>=1) || behavLV(3)<=0 ), error('Bad input in behavLV!'), end

filename = ['raw' tagno '.mat'];
disp(sprintf('\n\nLoading %s...',filename))
load(filename), 
td.DBname = DBname; db=1;
disp(['Loading DB:' td.DBname])
load(td.DBname)
dbdir = which(td.DBname); 
td.dbdir = dbdir;
save([td.dbdir(1:end-length(td.DBname)) td.DBname '_BAK.mat'],'db');
LDB = length(td.DBname);
if (db.lat(1,1) -db.lat(end,end))  < 0, db = flipdb(db,'lat'); save([td.dbdir(1:end-LDB) td.DBname],'db'); end
if (db.long(1,1)-db.long(end,end)) > 0, db = flipdb(db,'long');save([td.dbdir(1:end-LDB) td.DBname],'db'); end

% LDB = length(DBname);
% if (db.lat(1,1) -db.lat(end,end))  < 0, db = flipdb(db,'lat'); save([dbdir(1:end-LDB) DBname],'db'); end
% if (db.long(1,1)-db.long(end,end)) > 0, db = flipdb(db,'long');save([dbdir(1:end-LDB) DBname],'db'); end
disp(sprintf('\n=== Processing raw data of tag #%s ===',td.tagno))

if isfield(td,'deltat'), td.dt = td.deltat; disp('Creating field dt from detalt.'), end
if ~isfield(td,'dt'), td.dt = 10; disp('No dt field found in td struct, assuming dt=10 as default'), end

% Change dims of td.time, td.depth and td.temp
[A B] = size(td.time);  if A > 1, td.time =td.time';  end
[A B] = size(td.depth); if A > 1, td.depth=td.depth'; end
[A B] = size(td.temp);  if A > 1, td.temp =td.temp';  end

% Find indices of 24 hours
SR = 24*60/td.dt;
first24 = ceil(SR*(ceil(td.time(1))+1/1440-td.time(1)));
if first24 == 1
    td.d24 = [first24:SR:length(td.time) length(td.time)];
else
    td.d24 = [1 first24:SR:length(td.time) length(td.time)];
end
lengthtime = length(td.time); 
p=12.42; %period in hours
w=2*pi/(p/24); % Angular frequency
sint = sin(w*td.time)';
cost = cos(w*td.time)';
ts = td.depth';

td.diffs = ones(length(td.d24)-1,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%  TIDAL INFORMATION  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Extract summary statistics of sine fit
disp(sprintf('Extracting tidal information...\nIn %1.2f hour intervals using limit values:\n%1.2f > rmse, %1.2f < rsquare, %1.2f < amplitude\n\n',tideFL,tideLV(1),tideLV(2),tideLV(3)))
td.tideLV = tideLV; 
td.tideFL = round(60/td.dt*tideFL); % td.tideFL is fitlength in sample points
ons = ones(td.tideFL,1);
loop = 1:(lengthtime-td.tideFL);
rmse = zeros(1,(lengthtime-td.tideFL)); rsquare = rmse; ampli = rmse; out = rmse; tj=0;
tic
for i=loop
    intv = i:td.tideFL+i-1;
    [rmse(i) rsquare(i) ampli(i) out(i)]=lssinfit(ons,cost(intv), sint(intv),ts(intv));
    if ~mod(i,floor(lengthtime/100)), disp(sprintf('\b.')), end
    if ~mod(i,floor(lengthtime/9.99)), disp(sprintf('\b%1.0f%% (%1.2f sec)\n',100*i/loop(end),toc)),tic, end
end
disp(sprintf('\b100%% (%1.2f sec)\n\n',toc))

%% Find intervals with tidal information according to criteria
crit = [(rmse<td.tideLV(1) & rsquare>td.tideLV(2) & ampli>td.tideLV(3)) zeros(1,td.tideFL)];
i=1;
td.tideFound = zeros(1,lengthtime);
while i < length(crit)+1
    if crit(i) == 1
        td.tideFound(i:i+td.tideFL-1) = 1;
        i=i+td.tideFL-1;
    end
    i=i+1;
end

%% Determine the best fits
td.tideUsed = zeros(1,lengthtime); td.tide = zeros(1,length(td.d24)-1); td.tideBestfit = td.tide; td.rmse = td.tide+1;
lengthrmse = length(rmse);
for i=1:length(td.d24)-1
    intv = td.d24(i):min([td.d24(i+1)-1 lengthrmse]);
    if sum(crit(intv)) > 0
        [td.rmse(i) minind] = min(rmse(intv));
        td.tideBestfit(i) = minind + td.d24(i) -1;
        td.tideUsed(td.tideBestfit(i):td.tideBestfit(i)+td.tideFL-1) = 1;
        td.tide(i) = 1;
    end  
    if tj==1, break, end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%  BEHAVIOUR  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Extract summary statistics of sine fit
disp(sprintf('Classifying behaviour...\nIn %1.2f hour intervals using limit values:\n%1.2f > rmse, %1.2f < rsquare, %1.2f < amplitude\n\n',behavFL,behavLV(1),behavLV(2),behavLV(3)))
td.behavLV = behavLV; 
td.behavFL = round(60/td.dt*behavFL); % td.behavFL is fitlength in sample points
ons = ones(td.behavFL,1);
loop = 1:(lengthtime-td.behavFL);
rmse = zeros(1,(lengthtime-td.behavFL)); rsquare = rmse; ampli = rmse;  out = rmse; tj=0;
outlim = 3; % data with stud. residuals above this level are outliers
tic
for i=loop
    intv = i:td.behavFL+i-1;
    [rmse(i) rsquare(i) ampli(i) out(i)]=lssinfit(ons,cost(intv),sint(intv),ts(intv),outlim);
    if ~mod(i,floor(lengthtime/100)), disp(sprintf('\b.')), end
    if ~mod(i,floor(lengthtime/9.99)), disp(sprintf('\b%1.0f%% (%1.2f sec)\n',100*i/loop(end),toc)),tic, end
end

disp(sprintf('\b100%% (%1.2f sec)\n\n',toc))

% Make behaviour vector
crit = [(rmse<td.behavLV(1) & rsquare>td.behavLV(2) & ampli>td.behavLV(3)) zeros(1,td.behavFL)];
td.behavFound = ones(1,lengthtime)+1; 
td.behav = ones(1,length(td.d24)-1)+1; tj = 0;
td.behavrsq = zeros(1,length(td.behav)) + 0.5; % Uniformly distributed between the two modes
lengthrmse = length(rmse);

for i=1:length(td.d24)-1
    intv = td.d24(i):min([td.d24(i+1)-1 lengthrmse]);
    if sum(crit(intv)) > 0
        [minval minind] = min(rmse(intv));
        globind = minind + (td.d24(i):td.d24(i)+td.behavFL-1);
        td.behavFound(globind) = 1;
        td.behav(i) = 1;
        td.behavrsq(i) = rsquare(intv(minind));
    end  
    if tj==1, break, end
end

%% Determine values of f and G (computed by nodal.exe)
year = str2num(datestr(td.time_plot(1),10));
switch year
    case 1991
        td.f = [0.9821 1.0000 0.9821 1.1606 1.1063 1.0660 0.9646];
        td.G = [0.813 0.500 336.997 215.717 340.254 17.556 1.627]*pi/180;
    case 1992
        td.f = [0.9939 1.0000 0.9939 1.0658 1.0539 1.0336 0.9878];
        td.G = [101.810 0.500 349.273 217.542 79.507 18.603 203.620]*pi/180;
    case 1993
        td.f = [1.0064 1.0000 1.0064 0.9679 0.9923 0.9955 1.0129];
        td.G = [178.173 0.500 323.843 219.500 154.226 19.808 356.345]*pi/180;
    case 1994
        td.f = [1.0183 1.0000 1.0183 0.8798 0.9277 0.9557 1.0370];
        td.G = [278.692 0.500 335.642 217.446 255.567 18.999 197.385]*pi/180;
    case 1995
        td.f = [1.0282 1.0000 1.0282 0.8106 0.8685 0.9196 1.0573];
        td.G = [18.997 0.500 347.224 213.391 358.457 17.066 37.994]*pi/180;
    case 1996
        td.f = [1.0350 1.0000 1.0350 0.7659 0.8250 0.8932 1.0713];
        td.G = [119.140 0.500 358.645 207.676 102.886 14.119 238.281]*pi/180;
    case 1997
        td.f = [1.0378 1.0000 1.0378 0.7479 0.8061 0.8818 1.0771];
        td.G = [194.809 0.500 332.527 202.906 183.006 11.496 29.619]*pi/180;
    case 1998
        td.f = [1.0364 1.000 1.0364 0.7573 0.8160 0.8878 1.0741];
        td.G = [294.848 0.500 343.843 196.037 288.636 7.801 229.696]*pi/180;
    case 1999 %(#1432)
        td.f = [1.0308 1.0000 1.0308 0.7936 0.8524 0.9098 1.0625];
        td.G = [34.955 0.500 355.228 189.917 33.476 4.596 69.910]*pi/180;
    case 2000
        td.f = [1.0218 1.0000 1.0218 0.8553 0.9077 0.9435 1.0440];
        td.G = [135.201 0.500 6.752 185.279 136.880 2.316 270.402]*pi/180;
    case 2001 %(#2255)
        td.f = [1.0103 1.0000 1.0103 0.9385 0.9719 0.9829 1.0207];
        td.G = [211.263 0.500 341.026 184.538 213.322 2.136 62.527]*pi/180;
    case 2002
        td.f = [0.9978 1.0000 0.9978 1.0343 1.0351 1.0219 0.9957];
        td.G = [311.942 0.500 352.983 183.888 313.770 2.047 263.883]*pi/180;
    case 2003 %(#872)
        td.f = [0.9857 1.0000 0.9857 1.1314 1.0908 1.0564 0.9716];
        td.G = [52.858 0.500 5.176 185.16 53.25 2.853 105.717]*pi/180;
    case 2004 %(#6448)
        td.f = [0.9753 1.0000 0.9753 1.2172 1.1351 1.0837 0.9512];
        td.G = [153.998 0.500 17.592 188.045 152.095 4.341 307.995]*pi/180;
    case 2005 %(#1186)
        td.f = [0.9677 1.0000 0.9677 1.2806 1.1656 1.1024 0.9364];
        td.G = [230.942 0.500 352.752 194.118 225.149 7.301 101.885]*pi/180;
    case 2006
        td.f = [0.9637 1.000 0.9637 1.3134 1.1810 1.1118 0.9288];
        td.G = [332.377 0.500 5.464 198.946 323.309 9.558 304.754]*pi/180;
    case 2007
        td.f = [0.9639 1.0000 0.9639 1.3122 1.1804 1.1114 0.9291];
        td.G = [73.845 0.500 18.208 204.011 61.375 11.908 147.691]*pi/180;
    case 2008
        td.f = [0.9681 1.0000 0.9681 1.2771 1.1640 1.1015 0.9372];
        td.G = [175.279 0.500 30.921 208.804 159.552 14.152 350.558]*pi/180;
    case 2009
        td.f = [0.9759 1.0000 0.9759 1.2119 1.1324 1.0821 0.9524];
        td.G = [252.205 0.500 6.058 214.809 232.622 17.085 144.411]*pi/180;
    case 2010
        td.f = [0.9865 1.000 0.9865 1.1250 1.0873 1.0543 0.9732];
        td.G = [353.329 0.500 18.457 217.595 331.501 18.532 346.657]*pi/180;
end

%% Transform release and recapture positions to pixel coords
[row,col] = size(db.depth);
dlong = (db.long(1,col)-db.long(1,1))/(col-1);
dlat  = (db.lat(row,1)-db.lat(1,1))/(row-1);
R = mapmatrix(db.lat(1,1),db.long(1,1),dlat, dlong);
[td.x1 td.y1] = maptopix(R,td.catch_lat,td.catch_long);
[td.x0 td.y0] = maptopix(R,td.rel_lat,td.rel_long);
td.x1 = round(td.x1); td.y1 = round(td.y1); td.x0 = round(td.x0); td.y0 = round(td.y0);
% Check if release or recapture are on land
if db.land(td.y0,td.x0)
    warning('Release position is on land!')
    newx = td.x0 + [-1 0 1]; newx(newx < 1) = [];
    newy = td.y0 + [-1 0 1]; newy(newy < 1) = [];
    for x = newx, for y = newy
            if ~db.land(y,x), td.x0 = x; td.y0 = y; end
    end, end
end
if db.land(td.y1,td.x1)
    warning('Recapture position is on land!')
    newx = td.x1 + [-1 0 1]; newx(newx < 1) = [];
    newy = td.y1 + [-1 0 1]; newy(newy < 1) = [];
    for x = newx, for y = newy
            if ~db.land(y,x), td.x1 = x; td.y1 = y; end
    end, end
end

td

%% Creating *.mat file
filename = sprintf('tagdata%s',td.tagno);
disp(sprintf('Saving -> %s.mat <- in\n%s',filename,cd))
save(filename,'td')
disp(sprintf('\nDONE with tidal and behaviour extraction! \n\nNow run --> datalikelihood \n\nto create/update the datalikelihood matrix!\n'))

%% Plotting 
%t = td.time_plot;
t = td.time;
close all
parts = [1 round(length(t)/2) length(t)];
figure, hold on
for i = 1:2
    tt  = t(parts(i):parts(i+1));
    tsF = td.tideFound(parts(i):parts(i+1));
    tsU = td.tideUsed(parts(i):parts(i+1));
    [f_tsf f_sf] = stairs(tt,tsF);
    f_tsf = [f_tsf(1);f_tsf;f_tsf(end)]; f_sf = [0;f_sf;0];
    pl=patch(f_tsf,f_sf*min(td.depth),'g');
    set(pl,'EdgeColor',[0.6 1 0.4])
    set(pl,'FaceColor',[0.5 1 0.3]) % ligth green
    [f_tsf f_sf] = stairs(tt,tsU);
    f_tsf = [f_tsf(1);f_tsf;f_tsf(end)]; f_sf = [0;f_sf;0];
    pl=patch(f_tsf,f_sf*min(td.depth),'g');
    set(pl,'EdgeColor','none')
    set(pl,'FaceColor',[0 0.8 0]) % ligth green
end
plot(t,td.depth,'b'), 
xlabel('Julian day'),
%datetick('x','keeplimits'), xlabel('Date'),
axis tight, , title('Result of tidal classification')
 ylabel('Depth, m')
hold off
set(gcf,'position',[50 100 850 220])

figure, hold on
for i = 1:2
    tt  = t(parts(i):parts(i+1));
    tsF = td.behavFound(parts(i):parts(i+1));
    [f_tsf f_sf] = stairs(tt,abs(tsF-2));
    f_tsf = [f_tsf(1);f_tsf;f_tsf(end)]; f_sf = [0;f_sf;0];
    pl=patch(f_tsf,f_sf*min(td.depth),'g');
    set(pl,'EdgeColor',[0.6 1 0.4])
    set(pl,'FaceColor',[0.5 1 0.3]) % ligth green
end
plot(t,td.depth,'b'), 
xlabel('Julian day'),
%datetick('x','keeplimits'), xlabel('Date'),
axis tight,  title('Result of behaviour classification')
 ylabel('Depth, m')
hold off
set(gcf,'position',[50 340 850 220])

% plot release and recapture and show database area
figure, cmap = [1 1 1;0. 0.7 0.]; % white water, green land
surf(db.long,db.lat,db.land-1);
colormap(cmap), shading flat, hold on, view(2), axis tight, grid on
rel = plot(td.rel_long,td.rel_lat,'v','markersize',10,'markerfacecolor','g','markeredgecolor','k');
rec = plot(td.catch_long,td.catch_lat,'^','markersize',10,'markerfacecolor','r','markeredgecolor','k');
legend([rel rec],'Release position','Recapture position','location','best');
xlabel('Longitude, deg'), ylabel('Latitude, deg'), title('Release and recapture')
hold off
set(gcf,'position',[150 100 500 350])