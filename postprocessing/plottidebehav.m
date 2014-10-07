function plottidebehav(tagno)
%PLOTTIDEBEHAV  Plot the tidal and behaviour classification.
%   PLOTTIDEBEHAV(TAGNO)
%
%   - TAGNO indentifier as string for the tag.
%
%   the function assumes the following files are available
%
%     tagdataTAGNO.mat
%
%  EXAMPLE   
%   PLOTTIDEBEHAV('2255');
%
%   Date: 12/12 - 2007, ver. 0.51
%   HMM geolocation toolbox, IMM and DIFRES

load(['tagdata' tagno]);

t = td.time_plot;
close all

% Plot tidal classification
[f_tsf f_sf]=stairs(t(1:length(td.tideUsed)),td.tideUsed);
f_tsf = [f_tsf(1);f_tsf;f_tsf(end)]; f_sf = [0;f_sf;0];
p2=fill(f_tsf,f_sf*0.5*min(td.depth),[0 .8 0]); 
set(p2,'EdgeColor',[0 .8 0])
hold on

% Plot behaviour classification
[f_tsf f_sf]=stairs(t(1:length(td.behavFound)),abs(td.behavFound-2));
f_tsf = [f_tsf(1) f_tsf' f_tsf(end)]; f_sf = [0 f_sf' 0];
pl=fill(f_tsf,f_sf*0.5*min(td.depth)+(0.5*min(td.depth)),'g');
set(pl,'EdgeColor',[0.6 1 0.4],'FaceColor',[0.6 1 0.4])

% Plot depth record
plot(t,td.depth,'k'), 
datetick('x','keeplimits'), title('Result of tidal and behaviour classification')
xlabel('Date'), ylabel('Depth, m')

hold off

legend('tidal','behaviour','depth record','location','best')