function avimaker(tagno,o)
%AVIMAKER  Create an avi-file in based on a geolocation result.
%   AVIMAKER(TAGNO,OPTIONS)
%
%   - TAGNO indentifier as string for the tag to create avi of.
%
%     Optional arguments
%
%   - OPTIONS struct containing options in the fields (not all fields
%             need to be specified):
%
%   - RATE an integer that indentifies the sample rate of the avi
%   e.g. RATE = 3, stores every third frame.
%   default is 1.
%   - FPS number of frames per second. 
%   default is 5.
%   - RANGE defines the range of DAYS to be plotted. 
%   default is plotting of all days.
%   - MODE play animation 'backward' or 'forward' 
%   default is 'forward'.
%   - NO number of times to display the animation.
%   default is 1.
%   - COMP a string indicating the compressor to use. 
%   see help avifile for options (compression).
%   default is 'Cinepak'.
%   - TRACK plot the MPT on top of the animation. Input as string saying 
%   either 'on' or 'off'. Requires access to an mptTAGNO.mat file.
%   This file is created by the mptrack function.
%   default is 'off'.
%   - MOVNAME as string for a custom name for the avi-file.
%   default is eg. geolocation2255.avi (when TAGNO = '2255').
%   - ZOOM define the area of the domain to be plotted.
%   in the form [minlat maxlat minlong max long
%   e.g. ZOOM = [52 55 0 4]
%   default is ZOOM = [], i.e. no zoom.
%   - TYPE choose 'plain', 'fancy' or 'bw'
%   default is 'plain'.
%   - LOCK 'on' locks the colorscale, works only with TYPE='fancy'.
%   default is 'off'.
%
%   DEPENDENCIES - the function needs access to the following files
%
%     tagdataTAGNO.mat
%     resultTAGNO.mat
%     tidaldb.mat
%     cmap.mat
%     plotting.m
%     m_surf.m
%    (mptTAGNO.mat if TRACK = 'on')
%    (cmapfancy.mat, Ldistr.mat is TYP='fancy' or 'bw' respectively)
%    (M_Map package including the high-res coastline if TYP='fancy or 'bw')
%
%  EXAMPLE
%   AVIMAKER('2255')
%
%   options.rate = 3;
%   options.fps = 10;
%   options.range = 10:30;
%   options.mode = 'forward';
%   options.no = 2;
%   options.track = 'on';
%   options.lock = 'on';
%   options.type = 'fancy'
%   AVIMAKER('2255',options)
%
%   Date: 18/12 - 2007, ver. 0.53
%   HMM geolocation toolbox, IMM and DIFRES

clear mex, close all
filename = ['tagdata' tagno '.mat'];
load(filename)
filename = ['result' tagno '.mat'];
load(filename)
load('tidaldb.mat'),
[row,col,numbstor] = size(result.smooth);

if nargin < 2,
    o.rate = 1;
    o.fps = 5;
    o.range = 1:numbstor;
    o.mode = 'forward';
    o.no = 1;
    o.comp = 'cinepak';
    o.track = 'off';
    o.movname = ['geolocation' tagno];
    o.zoom = [1 col 1 row]; zm = o.zoom;
    o.type = 'plain';
    o.lock = 'off';
else
    if ~isfield(o,'rate'),    o.rate = 1; end
    if ~isfield(o,'fps'),     o.fps  = 5; end
    if ~isfield(o,'range'),   o.range  = 1:numbstor; end
    if ~isfield(o,'mode'),    o.mode = 'forward'; end
    if ~isfield(o,'no'),      o.no   = 1; end
    if ~isfield(o,'comp'),    o.comp = 'cinepak'; end
    if ~isfield(o,'track'),   o.track = 'off'; end
    if ~isfield(o,'movname'), o.movname = ['geolocation' tagno]; end
    if ~isfield(o,'zm'),      o.zoom = [1 row 1 col]; zm = o.zoom; else
        zm = o.zoom;
        dlong = (result.maplong(1,col)-result.maplong(1,1))/(col-1);
        dlat  = (result.maplat(row,1)-result.maplat(1,1))/(row-1);
        R = mapmatrix(result.maplat(1,1),result.maplong(1,1),dlat, dlong);
        zminp = zm; clear zm
        [zminp(3:4) zminp(1:2)] = maptopix(R,zminp(1:2),zminp(3:4));
        zm(1) = max([floor(zminp(3)) 1]);
        zm(2) = min([ceil(zminp(4)) col]);
        zm(3) = max([floor(zminp(2)) 1]);
        zm(4) = min([ceil(zminp(1)) row]);
    end
end
if ~isfield(o,'type'),     o.type = 'plain'; end
if ~isfield(o,'lock'),    o.lock = 'off'; end
if isunix, o.comp = 'none'; disp(sprintf('Compression is set to %s. (Because you are running UNIX)\nYou may experience problems with memory as no compression can be used.',o.comp)), end

if sum(strcmp(o.type,{'plain','bw','fancy'})) == 0, o.type = 'plain'; warning('Bad value for type, using type = "plain"'), end

if strcmp(o.mode,'backward')
    o.range = o.range(numbstor:-1:1);
end
if strcmp(o.track,'on')
    load(['mpt' tagno '.mat']);
else
    mpt = [];
end


fig = figure; set(fig,'position',[50 100 900 600])
mov = avifile(o.movname,'fps',o.fps,'quality',100,'compression',o.comp);
set(fig,'NextPlot','replacechildren');

for i = 1:o.no
    disp(sprintf('Sequence %i of %i...',i,o.no))
    for day = o.range(1:o.rate:end)
        clf
        switch o.type
            case 'plain'
                load cmap
                plotting(day,result,td,cmap,o.rate,mpt,zm)
            case 'bw'
                load Ldistr
                plottingbw(day,result,td,Ldistr,o.rate,mpt,zm)
            case 'fancy'
                if ~exist('m_proj.m','file'), close all, error('Cannot make a "fancy" plot because the M_map package seems not to be installed properly!'), end
                load cmapfancy
                plottingfancy(day,result,td,cmapfancy,o.rate,mpt,zm,o.lock)
        end
        F = getframe(gcf);
        mov = addframe(mov,F);
        disp(sprintf('Storing %i of %i',day,numbstor))
    end
end
mov=close(mov);
disp(sprintf('Stored animation in -> %s.mat <- in\n%s',o.movname,cd))