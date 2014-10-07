function viewdistr(distr,fps,range,mode,no,typ,land)
%VIEWDISTR  Plot a probability distribution in time.
%   VIEWDISTR(DISTR,FPS,RANGE,MODE,NO,TYPE,LAND)
%
%   - DISTR an array of size ROWxCOLxDAYS found as a subvar
%   in a result structure output from the hmmgeolocate function.
%
%     Optional arguments
%
%   - FPS number of frames per second. 
%   default is max.
%   - RANGE defines the range of DAYS to be plotted. 
%   default is plotting of all days.
%   - MODE play animation 'backward' or 'forward' 
%   default is 'forward'.
%   - NO number of times to show the animation.
%   default is 1.
%   - TYPE if 'fancy' or 'fancylock' the LAND var needs to be specified.
%   default is 'plain'.
%   - LAND a matrix containing land indicators.
%   default is empty. (required when using "fancy" or "fancylock")
%
%   DEPENDENCIES - the function needs access to the following files
%
%     cmap.mat
%    (cmapfancy.mat)
%
%  EXAMPLES
%   VIEWDISTR(result.smooth_plot,[],[],'backward')
%   this plots the entire distribution in a forward sweep.
%
%   VIEWDISTR(result.phi_plot,10,50:3:200,'backward',2)
%   this plots the distribution twice at 10 frames per second  
%   in a backwards sweep starting at day 200 end jumping 3 at a 
%   time until reaching day 50.
%
%   Date: 14/12 - 2007, ver. 0.56
%   HMM geolocation toolbox, DTU Informatics and DTU Aqua

if nargin < 7 || isempty(land)
    land = [];
    if nargin == 6 && (strcmp(typ,'fancy') || strcmp(typ,'fancylock'))
       disp('It is recommended to input a land array eg. db.land or result.land when using "fancy" or "fancylock".') 
    end
end
if nargin < 6 || isempty(typ)
    typ = 'plain';
end
if nargin < 5 || isempty(no) || no < 1, 
    no = 1;
end
if nargin < 4 || isempty(mode), 
    mode = 'forward';
end
if nargin < 3 || isempty(range),
    range = 1:size(distr,3);
end
if nargin < 2 || isempty(fps), 
    fps = 0;
end
if strcmp(typ,'fancy') || strcmp(typ,'fancylock'), 
    load cmapfancy
    cax(1) = 0; cax(2) = max(distr(:));
else 
    load cmap, 
end

% sweep direction
if strcmp(mode,'backward')
    range = range(length(range):-1:1);
end

% fps
if fps <= 0, 
    pause('off');
else
    pause('on');
end

delay = 1/(fps*1.1+eps);

switch typ
    case 'fancy'
        for j = 1:no
            for i = range
                post = distr(:,:,i);
                post(land) = 0.5*max(post(:));
                imagesc(post); 
                if length(range) > 1, title(sprintf('viewdistr, day %i',i)), end
                colormap(cmapfancy), drawnow, 
                pause(delay)
            end
        end
    case 'fancylock'
        for j = 1:no
            for i = range
                post = distr(:,:,i);
                maxp = max(post(:));
                if maxp > cax(2)/32, 
                    post(land) = .5*cax(2);
                    imagesc(post); 
                    caxis(cax);
                else
                    post(land) = .5*cax(2)/32;
                    imagesc(post);
                    caxis([0 cax(2)/32]), 
                end
                colorbar
                if length(range) > 1, title(sprintf('viewdistr, day %i',i)), end
                colormap(cmapfancy), drawnow, 
                pause(delay)
            end
        end
    case 'log'
        for j = 1:no
            for i = range
                post = log(distr(:,:,i));
                imagesc(post); 
                if length(range) > 1, title(sprintf('viewdistr, day %i',i)), end
                drawnow, 
                pause(delay)
            end
        end
    case 'plain'
        for j = 1:no
            for i = range
                post = distr(:,:,i);
                imagesc(post); 
                if length(range) > 1, title(sprintf('viewdistr, day %i',i)), end
                colormap(cmap), drawnow, 
                pause(delay)
            end
        end
end