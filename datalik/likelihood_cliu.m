function likelihood_cliu(fish_no,path_to_tags,tagname)
% Construction of likelihood function after (Le Bris et al, 2013 eq (2))
% using daily max depth and depth where tidal signal is detected, with
% coorespoding temperature.


addpath('../../hmm_smast/backfun/')
tag_name=[num2str(fish_no),'_raw'];
load([path_to_tags tagname]);
tagno=[num2str(fish_no),'_',tag.tag_id];
global tideLV
%tideLV  = [0.22 0.85 0.2 2.0];

mean_ampli=[];std_ampli=[];mean_phase=[];std_phase=[];
Twindow = 5;  %time window = 13 h
% for Twindow=10:0.5:20
nwindow = floor(Twindow*3600/tag.min_intvl_seconds); % window size in data point numbers

%p: M2 period in hours
p = 12.420601;
w=2*pi/(p/24); % Angular frequency

%tag.dnum=tag.dnum+4/24; %convert time to UTC


%load tidaldb.mat

ntimes = numel(tag.dnum);

int_dnum = floor(tag.dnum);
dbeg = int_dnum(1);
dend = int_dnum(end);
days = dbeg:dend;
ndays = numel(days);

%%% !!!!!!!!!!!!!!!!!!!!!!!!!!
%%% !!!!!!!!!!!comment this
%ndays=30;
%%% !!!!!!!!!!!!!!!!!!!!!!!!!!

sint = sin(w*tag.dnum);
cost = cos(w*tag.dnum);


figure(1);plot(tag.dnum,tag.depth,'b-')
hold on
%loop over day
td_detected=nan(size(tag.dnum));
td_used=td_detected;
day_tidal_depth=nan(size(days));
day_tidal_depth_temp=day_tidal_depth;
day_max_depth=nan(size(days));
for i=1:ndays;
    fprintf(['day: ' num2str(i) ' of ' num2str(ndays) ' \n'])
    days_idx=find(int_dnum == days(i));
    rmse=[];rsquare=[];ampli=[];
    if (days_idx(1)+nwindow > ntimes)
        [day_max_depth(i),day_max_dep_ind]=max(tag.depth(days_idx));
        day_max_depth_temp(i)=tag.temp(days_idx(day_max_dep_ind));
        break
    end
    [day_max_depth(i),day_max_dep_ind]=max(tag.depth(days_idx));
    day_max_depth_temp(i)=tag.temp(days_idx(day_max_dep_ind));
    %move window for each data point
    for j=1:numel(days_idx)
        if (days_idx(j)+nwindow > ntimes)
            break
        end
        intv=days_idx(j):min(ntimes,days_idx(j)+nwindow-1);
        [rmse(j) rsquare(j) ampli(j) jnk pred{j} mwh(j) alpha beta]=lssinfit(ones(numel(intv),1), cost(intv), sint(intv),tag.depth(intv));
        %phase(j) = deg2rad(191.25) - acos(alpha/ampli(j)); %phase is in Greenwich phase, 191.25 is phase lag of matlab datenum=0
        
        crit = (rmse(j)<tideLV(1) & rsquare(j)>tideLV(2) & ampli(j)>tideLV(3) & ampli(j)<tideLV(4));
        if crit==1
            td_detected(intv)=1;
        end
        
        %         figure(1);plot(tag.dnum(intv),pred{j},'r-')
        %         xlim([min(tag.dnum(intv)) max(tag.dnum(intv))])
    end
    
    % Find intervals with tidal information according to criteria
    crit = (rmse<tideLV(1) & rsquare>tideLV(2) & ampli>tideLV(3) & ampli<tideLV(4));
    %crit=ones(numel(rmse));
    %find best fit for each day and reconstruct corresponding fvcom signal
    if (sum(crit)>0)
        idx=find(rmse==min(rmse(crit)));
        idx=idx(1);
        intv=days_idx(idx):min(ntimes,days_idx(idx)+nwindow-1);
        td_used(intv)=1;
        day_tidal_depth(i)=mean(tag.depth(intv));
        day_tidal_depth_temp(i)=mean(tag.temp(intv));
    end
    
    
    
end

% figure(1);plot(tag.dnum.*td_detected,tag.depth.*td_detected,'g','LineWidth',3);
% figure(1);plot(tag.dnum.*td_used,tag.depth.*td_used,'y','LineWidth',3);
% figure(1);plot(tag.dnum.*td_used,tag.temp.*td_used,'r','LineWidth',1);
% xlim([min(tag.dnum),max(tag.dnum)])
% figure(2);hold on;
% plot(days,day_max_depth);
% plot(days,day_tidal_depth,'k-');

%% =============
% likelihood
% ==============


% warning off
close all;
% addpath('../hmm_smast/backfun/')
%load 12_raw.mat
addpath('../')



%load ~/Dropbox/Geolocation/projects/cod_zemeckis/tag_data/vemco.mat


% ==== Load  FVCOM  ====
global fvcom_tidaldb
load(fvcom_tidaldb)

% search within radius
%search_rad=200000; %m
search_rad=-1;

[xt,yt]=my_project(tag.release_lon,tag.release_lat,'forward');
if (search_rad>0)
    node_idx=find( sqrt((xt-fvcom.x).^2 + (yt-fvcom.y).^2)<=search_rad );
else
    node_idx=1:numel(fvcom.x);
end


% load bottom temperature
global bottom_temperature
% time
time_mjd = double(ncread(bottom_temperature,'time'));
ntimes = numel(time_mjd);
time_mdl = floor(time_mjd + datenum(1858,11,17,0,0,0));
h = ncread(bottom_temperature,'h');
nverts = numel(h);
% bottom temperature
fprintf('loading temperature data ... ');
t = ncread(bottom_temperature,'temp',[1 1],[nverts ntimes]);
fprintf('done loading temperature data\n\n');

% ================
% determine edges
nEdges = fvcom.nelems*3;
edge = zeros(nEdges,2);
icnt = 1;
for i=1:fvcom.nelems
    edge(icnt  ,1:2) = fvcom.tri(i,1:2);
    edge(icnt+1,1:2) = fvcom.tri(i,2:3);
    edge(icnt+2,1:2) = fvcom.tri(i,[3,1]);
    icnt = icnt + 3;
end;

% determine nodes surrounding nodes (no specific order)
ntsn = zeros(fvcom.nverts,1);
nbsn = nan(fvcom.nverts,12);

for i=1:nEdges
    i1 = edge(i,1);
    i2 = edge(i,2);
    [lmin,loc] = min(abs(nbsn(i1,:)-i2));
    if(lmin ~= 0);
        ntsn(i1) = ntsn(i1)+1;
        nbsn(i1,ntsn(i1)) = i2;
    end;
    [lmin,loc] = min(abs(nbsn(i2,:)-i1));
    if(lmin ~= 0);
        ntsn(i2) = ntsn(i2)+1;
        nbsn(i2,ntsn(i2)) = i1;
    end;
end;

% ==============


% compute depth std for neighboring nodes
std_dep=nan(size(fvcom.dep));
for nd=1:numel(node_idx)
    % progress output
    if (mod(nd,500)==0)
        fprintf('node: %d/%d\n',nd,fvcom.nverts)
    end
    
    nnode_list=nbsn(node_idx(nd),:);
    nnode_list=nnode_list(isfinite(nnode_list));
    std_dep(node_idx(nd))=std(fvcom.dep(nnode_list )-fvcom.dep(node_idx(nd)));
    
    
end



%% loop over days, calculate daily likelihood distribution
global std_temp_offset tag_depth_range tag_depth_accu tag_temp_accu
if isempty(std_temp_offset)
    std_temp_offset=2.0; %higher value is more inclusive
end
if isempty(tag_depth_range)
    tag_depth_range = 250; % in meter
end
if isempty(tag_depth_accu)
    tag_depth_accu = 0.008; % fraction of depth renge
end
if isempty(tag_temp_accu)
    tag_temp_accu = 0.1; % in degree C
end

% fig3=figure('units','normalized','position',[.05 .05 .6 .9]);
% plot_axis = [8e5,11e5,-2e5,2e5];
% pause on
% plot_axis = [8.2e5,9.5e5,-0.9e5,0.5e5];

% recap location
[xr,yr]=my_project(tag.recapture_lon,tag.recapture_lat,'forward');
dist_r = ( (fvcom.x-xr).^2+(fvcom.y-yr).^2 ).^0.5;

ObsLh=nan(ndays,numel(node_idx));

tide = zeros(1,ndays); % tide: activity level classification
% 2 - low activity
% 1 - moderate activity
% 0 - high activity

for i=1:ndays;
    %for i=1:12


    
    
    if isfinite(day_tidal_depth(i))
        tide(i)=1;
        %ObsLh_dep_tidal = normcdf((day_tidal_depth(i)+250*0.008)*ones(size(fvcom.dep)),fvcom.dep,std_dep)-...
         ObsLh_dep_tidal = normcdf((day_tidal_depth(i)+tag_depth_range*tag_depth_accu)*ones(size(fvcom.dep)),fvcom.dep,std_dep)-...
            normcdf((day_tidal_depth(i)-tag_depth_range*tag_depth_accu)*ones(size(fvcom.dep)),fvcom.dep,std_dep);
        ObsLh_dep_tidal = ObsLh_dep_tidal ./ max(ObsLh_dep_tidal);
        
        ObsLh_dep_total=ObsLh_dep_tidal;
    else
        tide(i)=0;
        ObsLh_dep = normcdf( -day_max_depth(i)*ones(size(fvcom.dep)), -fvcom.dep,std_dep) ./ ...
            normcdf(zeros(size(fvcom.dep)),-fvcom.dep,std_dep);
%         %
%         figure(1);clf;
%         patch('Vertices',[fvcom.x,fvcom.y],'Faces',fvcom.tri,'Cdata',ObsLh_dep,'edgecolor','none','facecolor','interp');
%         axis equal;%axis(plot_axis);
%         colorbar()
%         H = text(.82e6,1.7e5,['day: ' num2str(i) ' of ' num2str(ndays) ' ']);
%         set(H,'FontSize',16,'Color','k');
%         %pause
        %
        ObsLh_dep_total=ObsLh_dep;
        
    end

    
    % compute temp std for neighboring nodes
    std_temp=nan(size(fvcom.dep));
    
    fprintf('computing temp std for day %d\n',i)
    [~,iframe] = min(abs(int_dnum(i)-time_mdl));
    for nd=1:numel(node_idx)
        
        
        nnode_list=nbsn(node_idx(nd),:);
        nnode_list=nnode_list(isfinite(nnode_list));
        std_temp(node_idx(nd))=std(t(nnode_list,iframe )-t(node_idx(nd),iframe));
        
        
    end
    std_temp=std_temp+std_temp_offset;
    
    
    if isfinite(day_tidal_depth(i))
        ObsLh_temp_tidal = normcdf((day_tidal_depth_temp(i)+tag_temp_accu)*ones(size(fvcom.dep)),t(:,iframe),std_temp)-...
            normcdf((day_tidal_depth_temp(i)-tag_temp_accu)*ones(size(fvcom.dep)),t(:,iframe),std_temp);
        ObsLh_temp_tidal = ObsLh_temp_tidal ./ max(ObsLh_temp_tidal);
        ObsLh_temp_total=ObsLh_temp_tidal;
    else
        ObsLh_temp = normcdf((day_max_depth_temp(i)+tag_temp_accu)*ones(size(fvcom.dep)),t(:,iframe),std_temp)-...
            normcdf((day_max_depth_temp(i)-tag_temp_accu)*ones(size(fvcom.dep)),t(:,iframe),std_temp);
        ObsLh_temp = ObsLh_temp./max(ObsLh_temp);
        ObsLh_temp_total=ObsLh_temp;
        
        
    end

    
    %%%%%%%%%%%%%%%%%%%%%%
    % recapture location attraction likelihood
    %%%%%%%%%%%%%%%%%%%%%%
    
    t_remain=ndays-i;
    sigma = max( 1000*tag.recap_uncertainty_km, 0.5*25000*t_remain);
    AttLh = normpdf(dist_r,0,sigma); %25000: typical cod swimming speed (30 cm/s)
    AttLh = AttLh./max(AttLh);
    
    
    
    ObsLh(i,:)=ObsLh_dep_total.*ObsLh_temp_total.*AttLh;
    
    % release location treatment
    if i==1
        [xl,yl]=my_project(tag.release_lon,tag.release_lat,'forward');
        dist_rl = ( (fvcom.x-xl).^2+(fvcom.y-yl).^2 ).^0.5;
        [~,rel_idx] = find(dist_rl==min(dist_rl));
        ObsLh(i,rel_idx) = 1;
    end
    %
%         figure(2);clf;
%         patch('Vertices',[fvcom.x,fvcom.y],'Faces',fvcom.tri,'Cdata',ObsLh(i,:),'edgecolor','none','facecolor','interp');
%         axis equal;%axis(plot_axis);
%         colorbar()
%         H = text(.82e6,1.7e5,['day: ' num2str(i) ' of ' num2str(ndays) ' ']);
%         set(H,'FontSize',16,'Color','k');
%         pause
    %
    %     figure(fig3);clf;
    %     patch('Vertices',[fvcom.x,fvcom.y],'Faces',fvcom.tri,'Cdata',ObsLh,'edgecolor','none','facecolor','interp');
    %     figure(fig3);axis equal;axis(plot_axis);
    %     hold on;
    %     plot(xt,yt,'gp','MarkerSize',10,'MarkerFaceColor','none','MarkerEdgeColor','g')
    %
    %     colorbar()
    %     H = text(.82e6,1.7e5,['day: ' num2str(i) ' of ' num2str(ndays) ' ']);
    %     %H = text(.825e6,4e4,['day: ' num2str(i) ' of ' num2str(ndays) ' ']);
    %     set(H,'FontSize',16,'Color','k');
    %     if tide==1
    %         H = text(.82e6,1.6e5,'Tide: Yes ');
    %         set(H,'FontSize',16,'Color','g');
    %     else
    %         H = text(.82e6,1.6e5,'Tide: No ');
    %         set(H,'FontSize',16,'Color','r');
    %     end
    %
    %     vst=find(i==vemco.dnum-dbeg+1 & vemco.FISHID==fish_no);
    %     figure(fig3);plot(vemco.x(vst),vemco.y(vst),'w+')
    %     H = text(.82e6,1.5e5,vemco.STATION(vst));
    %     set(H,'FontSize',16,'Color','k');
    
    %     %export_fig(['22_out/22_',num2str(i,'%02d'),'.png']);
    %     saveas(fig3,[tag_name,'_out/',tag_name,'_',num2str(i,'%02d'),'.png'],'png');
    %     hold off
end

filename = sprintf('ObsLh%s',tagno);
disp(sprintf('Saving -> %s.mat <- \n',filename))
save(filename,'ObsLh','tide')


end

