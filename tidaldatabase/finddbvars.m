function finddbvars
%FINDDBVARS Reads the rawtidaldb and calculates database variances.
%   FINDDBVARS() 
%
%   DEPENDENCIES - the function needs access to the following files
%
%     rawtidaldb.mat
%
%   See the reference manual for further information on this function.
%
%   Date: 31/7 - 2007, ver. 0.5
%   HMM geolocation toolbox, IMM and DIFRES

disp(sprintf('\n\n=== Commencing calculation of database variances ==='))
load rawtidaldb

[row,col]=size(db.depth);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculate bathymetry roughness (variance) %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Z=repmat(NaN,row+2,col+2);
Z(2:end-1,2:end-1) = db.depth;

db.bathro = zeros(row,col);
disp('Computing the bathymetry variance ...')
for x=2:col+1
    for y=2:row+1
        dybder=[Z(y-1,x-1) Z(y+1,x+1) Z(y-1,x) Z(y+1,x) ...
                Z(y-1,x+1) Z(y+1,x-1) Z(y,x-1) Z(y,x+1)];
        % variance in a uniform distribution
        maxdyb = max(dybder)-min(dybder);
        db.bathro(y-1,x-1) = 1/12*maxdyb^2; 
    end
end
db.bathro(db.land) = 0;




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculate tidal roughness (variance)      %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Create time vector %%
startdate = datenum([2001 1 1 00 00 00]);
t= 0:1/144:3;

%% Store f and G from nodal.exe (fortran program)
% values of f and G that set t=0 at 1/1-01 00:01
f = [1.0103 1.0000 1.0103 0.9385 0.9719 0.9829 1.0207];
G = [211.263 0.500 341.026 184.538 213.322 2.136 62.527]*pi/180;
    
%% Tidal predictions %%
lng=0; plt=0;
pred = zeros(length(t),row,col);
disp('Predicting the tide at all positions...')
for plat=1:row
    for plong=1:col
        if ~db.land(plat,plong)
            for mode=1:length(db.freq)
                temp(mode,:)  = f(mode)*db.amp(plat+plt,plong+lng,mode) * ...
                                cos( db.freq(mode)*t + G(mode) - db.phase(plat+plt,plong+lng,mode));
            end
            pred(:,plat,plong) = db.depth(plat+plt,plong+lng) - sum(temp,1);
        end
    end
end

db.tidalro = zeros(row,col);
%% Computing tidal variance %%
lng=82; 
plt=68;
delta = [-1 -1 -1 0  0  1 1 1;
         -1  0  1 -1 1 -1 0 1];
disp('Computing the tidal variance...')
for x=2:col-1
    for y=2:row-1
        vars=[];
        for k=1:8
            xx=x+delta(2,k);
            yy=y+delta(1,k);
            if ~db.land(y,x) & ~db.land(y+delta(1,k),x+delta(2,k))
                vars = [vars var(pred(:,y,x) - pred(:,yy,xx))];
                % variance in a uniform distribution
                db.tidalro(y,x) = max(vars);
            end
        end
    end
end
db.tidalro(db.land) = 0;


db2=db;
%% Store in mat file %%
save tidaldb db

disp(sprintf('\nDONE! \n\nNow run --> datastrip \n\nto extract the raw tag data!\n'))