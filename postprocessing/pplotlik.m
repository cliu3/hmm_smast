function pplotlik(LIK,result,typ,day,lonrange,latrange)
%PPLOTLIK  Plot a data likelihood to be used in a paper.
%   PPLOTLIK(LIK,RESULT,TYP,DAY,LONRANGE,LATRANGE)
%
%   - LIK data likelihood to be plotted.
%   - RESULT a result struct created with hmmgeoloc.
%   - TYP resolution of the coast line. Options:
%   'c'=crude, 'l'=low, 'i'=intermediate, 'h'=high, 'f'=full
%   - DAY the day number to be plotted.
%
%     Optional arguments
%
%   - LONRANGE longitude range of plot area.
%   - LATRANGE latitude range of plot area.
%
%   The function requires the M_Map package including the 
%   high resolution coast line.
%   See http://www.eos.ubc.ca/~rich/map.html
%
%  EXAMPLE   
%   PPLOTLIK(LIK,result,'h',100,[-10 8],[48 60]);
%
%   Date: 12/12 - 2007, ver. 0.51
%   HMM geolocation toolbox, IMM and DIFRES

if nargin < 6
    lonrange = [-10 8];
    latrange = [48 60];
    disp('Using default long and latitude ranges in pplotlik.')
end

load distrpaperplot
proj = 'Gall-Peters'; %Rectangular
%close all
m_proj(proj,'lon',lonrange,'lat',latrange);
switch typ
    case 'c', m_gshhs_c('patch',[.5 .5 .5]);
    case 'l', m_gshhs_l('patch',[.5 .5 .5]);
    case 'i', m_gshhs_i('patch',[.5 .5 .5]);
    case 'h', m_gshhs_h('patch',[.5 .5 .5]);
    case 'f', m_gshhs_f('patch',[.5 .5 .5]);
end

m_grid('box','fancy','tickdir','in','linestyle','none');
hold on
names = fieldnames(LIK);
names = names(~strcmp(names,'type'));
numnames = length(names);
[row col days] = size(result.smooth);
% Combine all data likelihood to one array in Ltotal %
Ltotal = ones(row,col,days-1);
for j = 1:numnames
    Ltotal = Ltotal .* LIK.(names{j});
end
distr = Ltotal(:,:,day);
distr = makeplotstandard(normalise(distr));
[y,i] = max(distr(:));
lon = result.maplong(:); lon = lon(i);
lat = result.maplat(:); lat = lat(i);
m_contourf(result.maplong,result.maplat,distr,[0.05 0.5])
load distrpaperplot, colormap(distrpaperplot)
switch typ
    case 'c', m_gshhs_c('patch',[.5 .5 .5]);
    case 'l', m_gshhs_l('patch',[.5 .5 .5]);
    case 'i', m_gshhs_i('patch',[.5 .5 .5]);
    case 'h', m_gshhs_h('patch',[.5 .5 .5]);
    case 'f', m_gshhs_f('patch',[.5 .5 .5]);
end
xlabel('Longitude'), ylabel('Latitude')
hold off