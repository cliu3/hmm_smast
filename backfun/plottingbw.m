function plottingbw(day,result,td,Ldistr,rate,mpt,zm)
%PLOTTINGBW Plot the result of a geolocation.
%
%   This function is used by avimaker.m
%
%   This function should not be called manually by the user.
%
%   Date: 23/8 - 2007, ver. 0.5
%   HMM geolocation toolbox, IMM and DIFRES


if day == 1, bday=1; else bday = day-1; end
if td.behav(bday)==2, c1='none'; cc1='k'; c2=[0 0.6 0]; cc2='w'; s1=10; s2=14; 
else c1=[0 0.6 0]; cc1='w'; c2='none'; cc2='k'; s1=14; s2=10; end
days = size(result.smooth,3);
post = result.smooth(:,:,day);
post = makeplotstandard(normalise(post));
[y,i] = max(post(:));
lon = result.maplong(:); lon = lon(i)
lat = result.maplat(:); lat = lat(i)

% Map axes
proj = 'Gall-Peters'; %Rectangular
m_proj(proj,'lon',result.maplong(1,zm(3:4)),'lat',result.maplat(zm(2:-1:1),1)');
m_contourf(result.maplong,result.maplat,post,[.05 .5]), colormap(Ldistr),hold on
%m_plot(lon,lat)%,'markeredgecolor','k','markerfacecolor','k','markersize',10)
if ~isempty(mpt), m_plot(mpt.long(1:day),mpt.lat(1:day),'k'); m_plot(mpt.long(day),mpt.lat(day),'*k'); end
rel=m_plot(td.rel_long,td.rel_lat,'^','markersize',10,'markerfacecolor','g','markeredgecolor','k');
rec=m_plot(td.catch_long,td.catch_lat,'v','markersize',10,'markerfacecolor','r','markeredgecolor','k');
m_gshhs_i('patch',[.5 .5 .5]),
m_grid('box','fancy','tickdir','in','linestyle','none');
set(gca,'position',[0.2 0.07 0.8 0.73]);
legend([rel rec],'Release pos.','Reported recapture pos.','location','southeast')
xlabel('Longitude'), ylabel('Latitude')

% Time series axes
axes('position',[0.06 0.84 0.94 0.16])
[f_tsf f_sf]=stairs(td.time_plot,td.tideFound);
pl=fill(f_tsf,f_sf*min(td.depth),'g'); axis tight
set(pl,'EdgeColor','none')
set(pl,'FaceColor',[0.6 1 0.4]) % gr0n
datetick('x','mmm','keeplimits'), ylabel('Depth, m')
hold on
plot(td.time_plot,td.depth), axis([min(td.time_plot) max(td.time_plot) min(td.depth) max(td.depth)])
plot(td.time_plot(td.d24(day))*[1 1],[min(td.depth) max(td.depth)],'r','linewidth',2)
hold off

% Text axes
axes('position',[0 0 0.2 0.8]), axis([0 1 0 1]), axis off
text(0.05,0.95,sprintf('HMMgeolocate'),'Color','k','FontSize',16,'edgecolor','k','backgroundcolor','w')
text(0.05,0.88,sprintf('- Tag -'),'Color','k','FontSize',12)
text(0.05,0.83,sprintf('%s',td.tagno),'Color','k','FontSize',18)
text(0.05,0.78,sprintf('Current date'),'Color','k','FontSize',12)
text(0.05,0.73,sprintf('%s',datestr(td.time_plot(td.d24(day)),1)),'Color','k','FontSize',16)
text(0.05,0.68,sprintf('Day #: %i',day),'Color','k','FontSize',12)
if length(result.D) == 1
    text(0.05,0.60,sprintf('- Diffusivity -'),'Color','w','FontSize',12,'backgroundcolor',[0 0.6 0])
    text(0.05,0.55,sprintf('%3.2f km^2/day',result.D(1)),'Color','w','FontSize',14,'backgroundcolor',[0 0.6 0])
else
    text(0.05,0.60,sprintf('- D, low activity -'),'Color',cc1,'FontSize',12,'backgroundcolor',c1)
    text(0.05,0.55,sprintf('%3.2f km^2/day',result.D(1)),'Color',cc1,'FontSize',s1,'backgroundcolor',c1)
    text(0.05,0.45,sprintf('- D, high activity -'),'Color',cc2,'FontSize',12,'backgroundcolor',c2)
    text(0.05,0.40,sprintf('%3.2f km^2/day',result.D(2)),'Color',cc2,'FontSize',s2,'backgroundcolor',c2)
end
text(0.05,0.16,sprintf('Sample rate: %i',rate),'Color','k','FontSize',10) 
text(0.05,0.12,sprintf('version 0.6'),'Color','k','FontSize',10) 
text(0.05,0.08,sprintf('Created:  %s',date),'Color','k','FontSize',10) 
text(0.05,0.03,'IMM & DIFRES','Color','k','FontSize',13) 