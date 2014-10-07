function pplotdistr(result,field,typ,day,lonrange,latrange,pos)
%PPLOTDISTR  Plot a distribution to be used in a paper.
%   PPLOTDISTR(RESULT,FIELD,TYP,DAY,LONRANGE,LATRANGE,POS)
%
%   - RESULT a result struct created with hmmgeoloc.
%   - FIELD the field in the result struct to plot.
%   - TYP resolution of the coast line. Options:
%   'c'=crude, 'l'=low, 'i'=intermediate, 'h'=high, 'f'=full
%   - DAY the day number to be plotted.
%
%     Optional arguments
%
%   - LONRANGE longitude range of plot area.
%   - LATRANGE latitude range of plot area.
%   - POS marks the position given by POS = [lat long].
%
%   The function requires the M_Map package including the 
%   high resolution coast line.
%   See http://www.eos.ubc.ca/~rich/map.html
%
%  EXAMPLE   
%   PPLOTDISTR(result,'smooth','h',100,[-10 8],[48 60],[56 6]);
%
%   Date: 12/12 - 2007, ver. 0.52
%   HMM geolocation toolbox, DTU Informatics and Aqua

if nargin < 6 || (isempty(latrange) || isempty(lonrange))
    lonrange = [-10 8];
    latrange = [48 60];
    disp('Using default long and latitude ranges in pplotdistr.')
end

proj = 'Gall-Peters'; %Rectangular
if ~exist('m_proj.m','file'), close all, error('The M_map package seems not to be installed properly!'), end
m_proj(proj,'lon',lonrange,'lat',latrange);
m_grid('box','fancy','tickdir','in','linestyle','none'); hold on
distr = result.(field)(:,:,day);
%distr = makeplotstandard(distr);
distr = makeplotstandard(normalise(distr));
[y,i] = max(distr(:));
lon = result.maplong(:); lon = lon(i);
lat = result.maplat(:); lat = lat(i);
if ~strcmp(field,'Ldepth')
    m_contourf(result.maplong,result.maplat,distr,[.05 .5])
    load distrpaperplot, colormap(distrpaperplot)
else
    m_contourf(result.maplong,result.maplat,distr,[eps eps])
    load Ldistr, colormap(Ldistr)
end

m_plot(lon,lat,'.','markeredgecolor','k','markerfacecolor','k','markersize',10)
if exist('pos'), m_plot(pos(2),pos(1),'*','markeredgecolor','k','markerfacecolor','k','markersize',10), end
switch typ
    case 'c', m_gshhs_c('patch',[.5 .5 .5]);
    case 'l', m_gshhs_l('patch',[.5 .5 .5]);
    case 'i', m_gshhs_i('patch',[.5 .5 .5]);
    case 'h', m_gshhs_h('patch',[.5 .5 .5]);
    case 'f', m_gshhs_f('patch',[.5 .5 .5]);
end

%xlabel('Longitude'), 
%ylabel('Latitude')
hold off