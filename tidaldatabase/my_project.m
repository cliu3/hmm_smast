function [out_east,out_north] = my_project(in_east,in_north,direction) 

% Sample user-defined projection and inverse projection of (lon,lat) to (x,y) 
% Copy to my_project (not a member of the toolbox) and modify to suite you
%
% function [out_east,out_north] = my_project(in_east,in_north,direction) 
%
% DESCRIPTION:
%    Define projections between geographical and Euclidean coordinates 
%
% INPUT: 
%   in_east   = 1D vector containing longitude (forward) x (reverse)
%   in_north  = 1D vector containing latitude  (forward) y (reverse)
%   direction = ['forward' ;  'inverse']
%           
% OUTPUT:
%   (lon,lat) or (x,y) depending on choice of forward or reverse projection
%
% EXAMPLE USAGE
%    [lon,lat] = my_project(x,y,'reverse') 
%
% Author(s):  
%    Geoff Cowles (University of Massachusetts Dartmouth)
%
% Revision history
%   
%==============================================================================

%subname = 'my_project';
%fprintf('\n')
%fprintf(['begin : ' subname '\n'])

%------------------------------------------------------------------------------
% Parse input arguments
%------------------------------------------------------------------------------

ProjectDirection = 'forward';

if(direction == 'forward')
	ProjectDirection = 'forward';
        lon = in_east;
        lat = in_north;
else
	ProjectDirection = 'inverse';
        x = in_east;
        y = in_north;
end;



%------------------------------------------------------------------------------
% Perform the projection:  USER DEFINED 
% Example:  project/inverse project to state plane 1802
%------------------------------------------------------------------------------
if ispc
    if(ProjectDirection == 'forward')
        prj4_params = '-f  "%.12f" +proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs  +to +proj=tmerc +lat_0=42d50 +lon_0=-70d10 +k=0.9999666666666667 +x_0=900000 +y_0=0 +ellps=GRS80 +units=m +no_defs';
        [x,y] = cs2cs(lon, lat, prj4_params);
    else
        prj4_params = '-f  "%.12f" +proj=tmerc +lat_0=42d50 +lon_0=-70d10 +k=0.9999666666666667 +x_0=900000 +y_0=0 +ellps=GRS80 +units=m +no_defs +to +proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs';
        [lon,lat] = cs2cs(x, y, prj4_params);
    end;
    
else
    
    if(ProjectDirection == 'forward')
%         [x,y] = sp_proj('1802','forward',lon,lat,'m');
        [x,y] = arrayfun(@(a,b) LatLongToStatePlane(a,b,1802), lat, lon);
    else
%         [lon,lat] = sp_proj('1802','inverse',x,y,'m');
        [lat, lon] = arrayfun(@(a,b) StatePlaneToLatLong(a,b,1802), x, y);
    end;
    
end
%------------------------------------------------------------------------------
% Skagit, UTM, Zone 10 (see http://www.dmap.co.uk/utmworld.htm)
%------------------------------------------------------------------------------
%m_proj('UTM','longitude',[-123,-120],'latitude',[47,49],'zone',10,'hemisphere','north','ellipsoid','wgs84')
%m_proj get
%[x,y] = m_ll2xy(-122.530820 , 48.363114);
%fprintf('x %f y %f\n',x,y-1e7);
%fprintf('should be 534752, 5356766.\n')
%deltay = 1e7;
%
%if(ProjectDirection == 'forward')
%%	fprintf('Projecting from (lon,lat) to (x,y)\n');
%	[x,y]=m_ll2xy(lon,lat); 
%	y = y - deltay; %why?
%else
%%	fprintf('Inverse Projecting from (x,y) to (lon,lat)\n')
%	[lon,lat]=m_xy2ll(x,y+deltay); 
%end;
%
 
% set the output
if(ProjectDirection == 'forward')
  out_east = x;
  out_north = y;
else
  out_east = lon;
  out_north = lat;
end;

