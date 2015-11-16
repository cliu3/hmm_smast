function gen_tidaldb(lonmin,lonmax,latmin,latmax,delta)
%lonmin=-71;lonmax=-66;latmin=39;latmax=44;delta=.1;

global fvcom_tidaldb
load(fvcom_tidaldb)

% lon lat meshgrid
[lon,lat]=meshgrid(lonmin:delta:lonmax,latmax:-delta:latmin);


%-------------------------------------------------------------
% find cell with nearest cell center
%-------------------------------------------------------------
nelems = fvcom.nelems;
radlist = zeros(nelems,1);
cnt = 0;
[x,y] = my_project(lon,lat,'forward'); 
[ny, nx] = size(x);
incell = zeros(ny,nx);
for j=1:nx
for i=1:ny
  xpos=x(i,j);
  ypos=y(i,j);
  incell(i,j) = 0;
  radlist = sqrt((fvcom.xc-xpos).^2 + (fvcom.yc-ypos).^2);
  ii = 0;
  found = 0;
  while(ii <= 10 & found==0) 
    [minval,minloc]   = min(radlist);
    xtri    = fvcom.x(fvcom.tri(minloc,1:3));
    ytri    = fvcom.y(fvcom.tri(minloc,1:3));
    if(isintriangle(xtri,ytri,xpos,ypos));  
      incell(i,j)    = minloc;
      cnt = cnt + 1;
      found = 1;
    end;  
    radlist(minloc) = 1e6;
    ii = ii + 1;
  end;
end;
fprintf('j %d\n',j);
end;

fprintf('%d of %d points in the domain\n',cnt,nx*ny);

land = zeros(ny,nx);
land = (incell==0);

db.lat = lat;
db.long = lon;
db.land = land;
db.depth = zeros(ny,nx); %fake

hmin = (db.long(1,2)-db.long(1,1))*deglong(db.lat(1,1));
hmax = (db.long(1,2)-db.long(1,1))*deglong(db.lat(end,1));

db.h = mean([hmin hmax]);

save tidaldb db