function readdb
%READDB Read the plain database text files and store in a .mat file.
%   READDB() 
%
%   DEPENDENCIES - the function needs access to the following files
%
%     dbinfo.dat
%     latlongdep.dat
%     const*.dat (where * indicates constituent number)
%
%   See the reference manual for further information on this function.
%
%   Date: 12/12 - 2007, ver. 0.51
%   HMM geolocation toolbox, IMM and DIFRES

disp(sprintf('\n\n=== Reading the database text files ==='))

%% Read database informations %%
if ~exist('dbinfo.dat','file'), error('dbinfo.dat does not exist or is in the wrong dir'), end
fid = fopen('dbinfo.dat','r');
fgetl(fid);
landindicator = fscanf(fid,'%f',1); % Get land indicator
fgetl(fid); fgetl(fid);
rowcol = fscanf(fid,'%f',2);        % Read rows and columns
fgetl(fid); fgetl(fid);
noconst = fscanf(fid,'%i',1);       % Read number of constituents
fclose(fid);
row = rowcol(1); col = rowcol(2);

%% Read lat, long and depth %%
if ~exist('latlongdep.dat','file'), error('latlongdep.dat does not exist or is in the wrong dir'), end
fid = fopen('latlongdep.dat','r');
fgetl(fid);
db.lat = reshape(fscanf(fid,'%f',prod(rowcol)),col,row)';   % Read latitude
fgetl(fid); fgetl(fid);
db.long = reshape(fscanf(fid,'%f',prod(rowcol)),col,row)';  % Read longitude
fgetl(fid); fgetl(fid);
db.depth = -reshape(fscanf(fid,'%f',prod(rowcol)),col,row)'; % Read bathymetry
fclose(fid);
db.land = db.depth == -landindicator;
db.depth(db.land) = 0;
% hmin = (db.long(1,2)-db.long(1,1))*(111.320 + 0.373*sin(db.lat(1,1)*pi/180)^2)*cos(db.lat(1,1)*pi/180);
% hmax = (db.long(1,2)-db.long(1,1))*(111.320 + 0.373*sin(db.lat(size(db.depth,1),1)*pi/180)^2)*cos(db.lat(size(db.depth,1),1)*pi/180);
hmin = (db.long(1,2)-db.long(1,1))*deglong(db.lat(1,1));
hmax = (db.long(1,2)-db.long(1,1))*deglong(db.lat(end,1));

db.h = mean([hmin hmax]);

%% Read amplitude and phase for constituents %%
db.amp = zeros(row,col,noconst); db.phase = db.amp;
for i=1:noconst
    filename = ['const' num2str(i) '.dat'];
    if ~exist(filename,'file'), error([filename ' does not exist or is in the wrong dir']), end
    fid = fopen(filename,'r');
    fgetl(fid);
    db.name(i) = {fscanf(fid,'%s',1)};              % Read constituent name
    fgetl(fid); fgetl(fid);
    db.freq(i) = fscanf(fid,'%f',1) *pi/180*24;     % Read frequency and convert from deg/hour to rad/day
    fgetl(fid); fgetl(fid);
    db.amp(:,:,i) = reshape(fscanf(fid,'%f',prod(rowcol)),col,row)'*0.01;   % Read amplitude and convert from cm to m
    fgetl(fid); fgetl(fid);
    db.phase(:,:,i) = reshape(fscanf(fid,'%f',prod(rowcol)),col,row)'*pi/180;  % Read phase and convert from deg to rad
    fclose(fid);
end
db.amp(db.amp == landindicator*0.01) = 0;
db.phase(find(db.phase == landindicator*pi/180)) = 0;


%% Store in mat file %%
save rawtidaldb db

disp(sprintf('\nDONE! \n\nNow run --> finddbvars \n\nto calculate database variances!\n'))