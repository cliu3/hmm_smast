function p = probofvisit(tagno,area,dr,plotflag)
%PROBOFVISIT  Estimate the probability that the fish visited some area.
%   [p] = PROBOFVISIT(TAGNO,AREA,DAYRANGE)
%
%   - TAGNO indentifier as string for the tag.
%   - AREA specified by the four corners in the rectangle.
%   It is possible to use multiple rectangles to define area (see ex.).
%
%      Optional arguments
%
%   - DAYRANGE Specify the range of day to base the calculation on.
%   default is all days.
%   - PLOTFLAG Set to 1 for a plot of the area.
%
%   the function assumes the following files are available
%
%     resultTAGNO.mat
%     samptrack.m
%
%  EXAMPLE   
%   p = PROBOFVISIT('2255',[1 2 53.5 54.5; 1 3 54.5 55], 12:46);
%
%   Date: 7/7 - 2008, ver. 0.53
%   HMM geolocation toolbox, IMM and DIFRES

filename = ['result' tagno '.mat'];
disp(sprintf('\n\nLoading %s...',filename))
load(filename)
if nargin < 3, dr = 1:size(result.phi,3); end
if nargin < 4, plotflag = 0; end
disp(sprintf('\n=== Calculating the probability of tag #%s ===',tagno))
disp(sprintf('having visited the area:'))
for k = 1:size(area,1)
    disp(sprintf('%3.1f - %3.1f long, %3.1f - %3.1f lat',area(k,1), area(k,2), area(k,3), area(k,4)))
end
disp(sprintf('\nWithin the time from\n %s to %s...\n',datestr(result.time(dr(1))), datestr(result.time(dr(end)))))
filename = ['datalikelihood' tagno '.mat'];
disp(sprintf('Loading %s...\n',filename))
load(filename)
if exist('L','var'), LIK = L; clear L; end

n = 1000;
disp('Estimating probability')
tr = samptrack(result,LIK,n);
tr.long = tr.long(dr,:);
tr.lat  = tr.lat(dr,:);

longind = zeros(size(tr.long));
latind = longind;
for k = 1:size(area,1)
    longind = longind + (tr.long > area(k,1) & tr.long < area(k,2));
    latind  = latind  + (tr.lat > area(k,3)  & tr.lat < area(k,4));
end
longind = longind > 0;
latind  = latind > 0;

visits = sum(sum(longind.*latind,1) > 0);
p = visits/n;

if plotflag == 1,
    proj = 'Gall-Peters'; %Rectangular
    lonrange = [-10 8];
    latrange = [48 60];
    if ~exist('m_proj.m','file'), 
        surf(result.maplong,result.maplat,result.land-1); cmap = [1 1 1;0. 0.7 0.];
        colormap(cmap), 
        shading flat, hold on
        axis tight, view(2), shading flat
        hold on
        for k = 1:size(area,1)
            plot([area(k,1) area(k,2) area(k,2) area(k,1) area(k,1)],[area(k,3) area(k,3) area(k,4) area(k,4) area(k,3)],'color','k','linewidth',2)
        end
        xlabel('Longitude'), ylabel('Latitude')
        title(sprintf('\nTime span\n %s to %s',datestr(result.time(dr(1))), datestr(result.time(dr(end)))))
        hold off
    else
        m_proj(proj,'lon',lonrange,'lat',latrange);
        m_gshhs_l('patch',[.5 .5 .5]); 
        m_grid('box','fancy','tickdir','in','linewidth',10,'linestyle','none');
        hold on
        for k = 1:size(area,1)
            m_plot([area(k,1) area(k,2) area(k,2) area(k,1) area(k,1)],[area(k,3) area(k,3) area(k,4) area(k,4) area(k,3)],'color','k')
        end
        xlabel('Longitude'), ylabel('Latitude')
        title(sprintf('\nTime span\n %s to %s',datestr(result.time(dr(1))), datestr(result.time(dr(end)))))
        hold off
    end
end