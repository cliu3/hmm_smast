function [] = check_harmonics(tag); 
% Check harmonics from a fixed bottom tag and compare with FVCOM database at the same location
%
% function [] = check_harmonics(tag)  
%
% DESCRIPTION:
%   Check harmonics from a fixed bottom tag and compare with FVCOM database at the same location
%
% INPUT 
%   tagfile = tagid (will be used to open tagid.mat) 
%
% OUTPUT:
%   screen dump of harmonics comparison  
%
% EXAMPLE USAGE
%    load S12689.mat;
%    check_harmonics(tag); 
%
% Author(s):  
%    Geoff Cowles (University of Massachusetts Dartmouth)
%
% Revision history
%   
%==============================================================================

% parameters
dd = 1.;  %distance in degrees to set zoombox around release/recapture area

% get screen info for figure positioning
set(0,'Units','pixels') ;
scnsize = get(0,'ScreenSize');


deltat_minutes = tag.min_intvl_seconds/60.;
ranger = 20:20+ceil((60*24*3)/deltat_minutes); 
zeta = tag.depth(ranger); 
time = tag.dnum(ranger); 

% tag interval time
fprintf('\n\n');
fprintf('==== reporting information for tag %s =====\n',tag.tag_id);
fprintf('fish id %d\n',tag.fish_id);
fprintf('release time : %s\n',datestr(tag.dnum(1)));
fprintf('recap   time : %s\n',datestr(tag.dnum(end)));
fprintf('days at large: %d\n',round(tag.dnum(end)-tag.dnum(1)));
fprintf('minimum time interval in seconds: %d\n',tag.min_intvl_seconds);
fprintf('maximum time interval in seconds: %d\n',tag.max_intvl_seconds);

fig1 = figure;
subplot(3,1,1)
plot(tag.dnum-tag.dnum(1),-tag.depth);
ylabel('depth (m)')

subplot(3,1,2)
plot(tag.dnum-tag.dnum(1),tag.temp);
ylabel('temp(C)')

subplot(3,1,3)
%plot((time-floor(time(1)))*24-4,zeta); 
plot(time-time(1),zeta,'bx-');
ylabel('depth (m)')
xlabel('days to process tide signal');

% now check the harmonics
[amp_tag,pha_tag] = wrap_ttide_harmonics(zeta,deltat_minutes,time(1)) ;

% interpolate from FVCOM database at release lon/lat
% find nearest node in FVCOM domain 
load ../data/fvcomdb_gom3_v2.mat;
rad = ((fvcom.lon - tag.release_lon).^2 + (fvcom.lat - tag.release_lat).^2);
[dist,node] = min(rad);
fprintf('nearest nodes is %d at distance %f km \n',node,dist*75.);

% report harmonics comparison
fprintf('   component    fvcom_amp(m)   tag_amp(m)   fvcom_pha(degG)  tag_pha(deg G)\n');

comps = {'M2','N2','S2','O1','K1','K2','P1','Q1'};
ncomps = numel(comps);
for n=1:3;  %ncomps
    icomp = 0;
    for i=1:fvcom.ncomps
      if(strcmp(char(comps(n)),fvcom.comps(i))==1);
        icomp = i;
      end;
    end;
    if(icomp==0);
      fprintf('component %s does not seem to exist in the fvcom database\n',char(comps(n)));
      error('stopping ...\n');
    end;
    amp_mdl(n) = fvcom.amp(node,icomp)*.01; %convert from cm to m
    pha_mdl(n) = fvcom.phase(node,icomp);
    if(amp_tag(n) > 0.)
      fprintf('%18s %18.2f %18.2f %20.0f %20.0f\n',char(comps(n)),amp_mdl(n),amp_tag(n),pha_mdl(n),pha_tag(n));
    end;
end; %comp loop

% set figure position
position = get(fig1,'Position');
outerpos = get(fig1,'OuterPosition'); borders = outerpos-position;
edge = -borders(1)/2;
pos1 = [edge,...
        scnsize(4) * (3/3),...
        scnsize(3)/2 - edge,...
        scnsize(4)/3];
set(fig1,'OuterPosition',pos1)

