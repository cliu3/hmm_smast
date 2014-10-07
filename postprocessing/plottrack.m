function [handle] = plottrack(track,linsty,pltyp,opttrack,clean)
%PLOTTRACK  Plot a track (eg. mpt or sampled track)
%   HANDLE = PLOTTRACK(TRACK,LINSTY,PLTYP,OPTTRACK)
%
%   - TRACK a struct containing an output from eg. samptrack.
%
%     Optional arguments
%
%   - LINSTY line style for the tracks. Type help plot for more info
%   - PLTYP specifies the plot type, either '1d' or '2d'
%   - OPTTRACK input another track for comparison
%   - CLEAN do not show track with random variation within each grid cell
%   default is a non-clean track (clean = 0)
%
%  EXAMPLE   
%   PLOTTRACK(track,'x-','1d',track2,1)
%   where the variables track and track2 is created eg. by samptrack
%   or mptrack.
%
%   Date: 22/10 - 2008, ver. 0.52
%   HMM geolocation toolbox, DTU Informatics and DTU Aqua

if nargin < 2 || isempty(linsty), linsty = '-b'; end
if nargin < 3 || isempty(pltyp), pltyp = '2d'; end
if nargin < 4 || isempty(opttrack), opttrack = []; end
if nargin < 5, clean = 0; end

if clean ~= 0, track.long = track.long_clean; track.lat = track.lat_clean; end

cmap = [1 1 1;0. 0.7 0.]; % white water, green land
lw = 1;
if strcmp(pltyp,'2d') || strcmp(pltyp,'2D')
    switch size(track.long,2) 
        case 0
            error('Wrong track input!')
        otherwise
            a = surf(track.maplong,track.maplat,track.land-1);
            colormap(cmap), shading flat, hold on
            view(2), axis tight, grid on
            if ~isempty(opttrack), opt = plot(opttrack.long,opttrack.lat,'r--','linewidth',1.5); if size(opttrack.long,2)>1, lw = 2; end, end
            st = plot(track.long,track.lat,linsty,'linewidth',lw);
            rel = plot(track.long(1),track.lat(1),'v','markersize',10,'markerfacecolor','g','markeredgecolor','k');
            rec = plot(track.long(end,:),track.lat(end,:),'^','markersize',10,'markerfacecolor','r','markeredgecolor','k');
            hold off
            if ~isempty(opttrack), legend([rel rec st(end) opt(end)],'Release position','Recapture position','Track','Optional track','location','best')
            else legend([rel rec st(end)],'Release position','Estimated recapture position','Track','location','best'), end
            xlabel('Longitude, deg'), ylabel('Latitude, deg'), title('plottrack')
            handle = a;
    end
elseif strcmp(pltyp,'1d') || strcmp(pltyp,'1D')
    switch size(track.long,2)
        case 0
            error('Wrong track input')
        otherwise
            a=subplot(211); hold on
            if ~isempty(opttrack), opt = plot(opttrack.long,'r--','linewidth',1.5); if size(opttrack.long,2)>1, lw = 2; end,end
            st = plot(1:length(track.long),track.long,linsty,'linewidth',lw); hold off
            xlabel('Time, days at liberty'), ylabel('Longitude, deg'), axis tight
            if ~isempty(opttrack), legend([st(end) opt(end)],'Track','Optional track','location','best')
            else legend(st(end),'Track','location','best'), end
            b=subplot(212); hold on
            if ~isempty(opttrack), opt = plot(opttrack.time,opttrack.lat,'r--','linewidth',1.5); end
            st = plot(track.time,track.lat,linsty,'linewidth',lw); datetick('x'); hold off
            xlabel('Time, date'), ylabel('Latitude, deg'), axis tight
            if ~isempty(opttrack), legend([st(end) opt(end)],'Track','Optional track','location','best')
            else legend(st(end),'Track','location','best'), end
            handle = [a b];
    end
end