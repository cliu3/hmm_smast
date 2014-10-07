function plotud(result,range,typ)
%PLOTUD  Plot a utilisation distribution.
%   PLOTUD(RESULT,RANGE,TYPE)
%
%   - RESULT output from the hmmgeolocate function.
%
%     Optional arguments
%
%   - RANGE defines the range of DAYS to be plotted eg. 1:10. 
%   default is plotting of all days.
%   - TYPE either 'plain', 'fancy' or 'log'.
%   default is 'plain'.
%
%   DEPENDENCIES - the function needs access to the following files
%
%     cmap.mat
%    (cmapfancy.mat)
%
%  EXAMPLES
%   PLOTUD(result,1:80,'log')
%
%   Date: 14/12 - 2007, ver. 0.56
%   HMM geolocation toolbox, IMM and DIFRES

if nargin < 3 || isempty(typ)
    typ = 'plain';
end
if nargin < 2 || isempty(range),
    range = 1:size(result.smooth,3);
end

if strcmp(typ,'fancy') || strcmp(typ,'fancylock'), 
    load cmapfancy
else 
    load cmap, 
end

UD = sum(result.smooth(:,:,range),3);
UD(result.land) = 0;
UD = normalise(UD);

switch typ
    case 'fancy'
        UD(result.land) = 0.5*max(UD(:));
        surf(result.maplong,result.maplat,UD); colormap(cmapfancy)
        axis tight, view(2), shading flat
        if length(range) > 1, title(sprintf('UD. Fancy. Days: %i - %i',range(1),range(end))), end
    case 'log'
        UD = log(UD); UD(result.land)=0;
        UD(result.land) = min(UD(:)) -abs(0.1*(max(UD(:))-min(UD(:))));
        surf(result.maplong,result.maplat,UD);
        axis tight, view(2), shading flat
        if length(range) > 1, title(sprintf('UD. Log. Days: %i - %i',range(1),range(end))), end
    case 'plain'
        UD(result.land) = -0.1*max(UD(:));
        surf(result.maplong,result.maplat,UD); colormap(cmap)
        axis tight, view(2), shading flat
        if length(range) > 1, title(sprintf('UD. Days: %i - %i',range(1),range(end))), end
end