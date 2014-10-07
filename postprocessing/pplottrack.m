function pplottrack(track,typ,lonrange,latrange,dr,opttrack,db,landcolor)
%PPLOTTRACK  Plot a track to be used in a paper
%   PPLOTTRACK(TRACK,TYP,LONRANGE,LATRANGE,DAYRANGE,OPTTRACK,DB,LANDCOLOUR)
%
%   - TRACK a track struct created eg. with mptrack.
%
%     Optional arguments
%
%   - TYP resolution of the coast line. Options:
%   'c'=crude, 'l'=low, 'i'=intermediate, 'h'=high, 'f'=full.
%   default is 'l'.
%   - LONRANGE longitude range of plot area.
%   default is [-10 8].
%   - LATRANGE latitude range of plot area.
%   default is [48 60].
%   - DAYRANGE the range of days to plot eg. 1:50.
%   default is all days within the time at liberty.
%   - OPTTRACK add another track to the plot.
%   - DB specify the bathymetry database.
%   default is the North Sea.
%   - LANDCOLOUR define the colour of the land areas.
%   default is grey i.e. [0.5 0.5 0.5] in rgb.
%
%   The function requires the M_Map package including the 
%   high resolution coast line.
%   See http://www.eos.ubc.ca/~rich/map.html
%
%  EXAMPLE   
%   PPLOTTRACK(mpt,'h',[-10 8],[48 60]);
%   PPLOTTRACK(tr,'i',[13.5 16.5],[54 56.5],1:50,[],db,'g')
%
%   Date: 22/10 - 2008, ver. 0.52
%   HMM geolocation toolbox, DTU Informatics and DTU Aqua

if nargin < 2 || isempty(typ), typ = 'l'; end
if nargin < 3 || isempty(lonrange)
    lonrange = [min(track.long(:))-1 max(track.long(:))+1];
end
if nargin < 4 || isempty(latrange)
    latrange = [min(track.lat(:))-1 max(track.lat(:))+1];
end
if nargin >= 5 && ~isempty(dr),
    track.long = track.long(dr); track.lat = track.lat(dr);
    if nargin >= 6, opttrack.long = opttrack.long(dr); opttrack.lat = opttrack.lat(dr); end
end
if nargin < 6, opttrack = []; end
if nargin < 7, db = []; end
if nargin < 8, landcolor = [.5 .5 .5]; end

%clf
proj = 'Gall-Peters'; %Rectangular
m_proj(proj,'lon',lonrange,'lat',latrange);
if isstruct(db), m_contourf(db.long,db.lat,db.depth,'linestyle','none'), hold on, colormap gray, end
switch typ
    case 'c', m_gshhs_c('patch',landcolor);
    case 'l', m_gshhs_l('patch',landcolor);
    case 'i', m_gshhs_i('patch',landcolor);
    case 'h', m_gshhs_h('patch',landcolor);
    case 'f', m_gshhs_f('patch',landcolor);
end
hold on
if isstruct(opttrack)
    m_plot(opttrack.long,opttrack.lat,'color',[.7 .7 .7])
    m_plot(opttrack.long(1),opttrack.lat(1),'o','markersize',3,'markeredgecolor','k','markerfacecolor','k')
    m_plot(opttrack.long(end),opttrack.lat(end),'o','markersize',3,'markeredgecolor','k','markerfacecolor',[.7 .7 .7])
end
m_plot(track.long,track.lat,'color','k')
m_plot(track.long(1),track.lat(1),'o','markeredgecolor','k','markerfacecolor','g','markersize',10)
%m_plot(track.long(end),track.lat(end),'o','markeredgecolor','k','markerfacecolor',[.7 .7 .7])
m_plot(track.long(end),track.lat(end),'o','markeredgecolor','k','markerfacecolor','r','markersize',10)
m_grid('box','fancy','tickdir','in','linewidth',10,'linestyle','none');
xlabel('Longitude'), 
ylabel('Latitude')
hold off