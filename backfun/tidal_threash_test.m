function tidal_threash_test(fish_no)
%calculate rmse between tag and fvcom, create tidal signal constraint on
%likelihood funcion

tag_name=[num2str(fish_no),'_raw'];
load([tag_name,'.mat'])
tagno=[num2str(fish_no),'_',tag.tag_id];

global tideLV
%tideLV  = [0.42 0.85 0.2 2.0];
%tideLV  = [0.42 0.92 0.2 2.0];
%tideLV  = [0.42 0.95 0.2 2.0];

mean_ampli=[];std_ampli=[];mean_phase=[];std_phase=[];
Twindow = 13;  %time window = 13 h
% for Twindow=10:0.5:20
nwindow = floor(Twindow*3600/tag.min_intvl_seconds); % window size in data point numbers

%p: M2 period in hours
p = 12.420601;
w=2*pi/(p/24); % Angular frequency

%tag.dnum=tag.dnum+4/24; %convert time to UTC

plot_axis = [8e5,11e5,-2e5,2e5];


ntimes = numel(tag.dnum);

int_dnum = floor(tag.dnum);
dbeg = int_dnum(1);
dend = int_dnum(end);
days = dbeg:dend;
ndays = numel(days);

sint = sin(w*tag.dnum);
cost = cos(w*tag.dnum);


% ==== Load & interp FVCOM  ====
global fvcom_tidaldb
load(fvcom_tidaldb)

%search_rad=2000; %m
search_rad=-1; % minus value: use all nodes


%Define tidal constituents and values
inpcon = fvcom.comps;



filename=['ObsLh_',num2str(fish_no),'.mat'];
disp(sprintf('Loading %s...\n',filename))
%load(filename)



% fprintf('reconstructing FVCOM tidal signal ... \n');
% 
% %find FVCOM harmonic constants at grid points within radius from tag location
% 
% [xt,yt]=my_project(tag.release_lon,tag.release_lat,'forward');
% [xr,yr]=my_project(tag.recapture_lon,tag.recapture_lat,'forward');
% if (search_rad>0)
%     node_idx=find( sqrt((xt-fvcom.x).^2 + (yt-fvcom.y).^2)<=search_rad );
% else
%     node_idx=1:numel(fvcom.x);
% end
% [~,node_tag]=min( sqrt((xt-fvcom.x).^2 + (yt-fvcom.y).^2) );
% pha1=zeros(numel(node_idx),numel(fvcom.comps));
% amp1=zeros(numel(node_idx),numel(fvcom.comps));
% 
% 
% 
% %Load names,freq from ttide database and merge it to variable tidecon
% %following ttide format
% ttstuff = load('t_constituents.mat');
% ncon     = length(inpcon);
% % Allocate names and freq to later load from ttide database
% names   = cell(ncon,1);
% freq    = zeros(ncon,1);
% tidecon1 = zeros(ncon,4);
% 
% for nd=1:numel(node_idx)
%     if (mod(nd,500)==0)
%         fprintf('node: %d/%d\n',nd,numel(node_idx))
%     end
%     pha1(nd,:)=fvcom.phase(node_idx(nd),:);
%     amp1(nd,:)=fvcom.amp(node_idx(nd),:)*0.01; %cm to m
%     for ic = 1:length(inpcon)
%         names(ic) = inpcon(ic);
%         idf = strcmp(ttstuff.const.name,inpcon(ic));
%         freq(ic,:)  = ttstuff.const.freq(idf,:);
%         tidecon1(ic,:) = [amp1(nd,ic) 0.0 pha1(nd,ic) 0.0];
%     end
%     
%     % create timeseries with the defined tidal harmonics info with ttide
%     eta1{nd}= t_predic(tag.dnum,names,freq,tidecon1);
%     
% end



%% fitting
figure(100);plot(tag.dnum,-tag.depth);hold on;
set(gca,'xtick',[min(floor(tag.dnum)):max(floor(tag.dnum))])
datetick('x','mmm dd','keepticks')
thresh=0.2461;
rmse_con=ones(ndays,numel(fvcom.x));
%figh=figure('units','normalized','position',[.05 .05 .6 .9]);
%loop over day
rmse_tag=[];
day_ampli=nan(ndays);
for i=1:ndays;
    fprintf(['day: ' num2str(i) ' of ' num2str(ndays) ' \n'])
    days_idx=find(int_dnum == days(i));
    rmse=[];rsquare=[];ampli=[];
    if (days_idx(1)+nwindow > ntimes)
        break
    end
    [day_depth(i),day_max_dep_ind]=max(tag.depth(days_idx));
    day_temp(i)=tag.temp(days_idx(day_max_dep_ind));
    %move window for each data point
    for j=1:numel(days_idx)
        if (days_idx(j)+nwindow > ntimes)
            break
        end
        intv=days_idx(j):min(ntimes,days_idx(j)+nwindow-1);
        [rmse(j) rsquare(j) ampli(j) jnk jnk mwh(j) alpha beta]=lssinfit(ones(numel(intv),1), cost(intv), sint(intv),tag.depth(intv));
        %phase(j) = deg2rad(191.25) - acos(alpha/ampli(j)); %phase is in Greenwich phase, 191.25 is phase lag of matlab datenum=0
        
        %             figure(1);plot(tag.dnum(intv),tag.depth(intv),'bx-')
        %             hold on
        %             figure(1);plot(tag.dnum(intv),pred,'r-')
        %             xlim([min(tag.dnum(intv)) max(tag.dnum(intv))])
    end
    
    % Find intervals with tidal information according to criteria
    crit = (rmse<tideLV(1) & rsquare>tideLV(2) & ampli>tideLV(3) & ampli<tideLV(4));
    %crit=ones(numel(rmse));
    %find best fit for each day and reconstruct corresponding fvcom signal
    if (sum(crit)>0)
        tide(i)=2;
        
        idx=find(rmse==min(rmse(crit)));
        idx=idx(1);
        intv=days_idx(idx):min(ntimes,days_idx(idx)+nwindow-1);
        time=tag.dnum(intv);
        day_ampli(i) = ampli(idx);
        
        eta_tag{i}=tag.depth(intv)-mean(tag.depth(intv));
        %eta_tag{i}=tag.depth(intv);
        
        % nonlinear sine fit
        f=fit(time,eta_tag{i},'sin1');
        eta_tag_fit{i}=f(time);
        
        figure(100);plot(time,-(eta_tag_fit{i}+mean(tag.depth(intv))),'g','LineWidth',3);hold on;
        %plot(time,eta_tag_fit{i},'g');
        
        
%         eta_tagnode=eta1{find(node_idx==node_tag)}(intv)-mean(eta1{find(node_idx==node_tag)}(intv));
%         %plot(time,eta_tagnode,'r');
%         
%         % calculate rmse_tag
%         
%         rmse_tag(i)=rms(eta_tagnode-eta_tag_fit{i});
%         
%         % calculate rmse map
%         rmse_eta=nan(size(fvcom.x));
%         
%         for nd=1:numel(node_idx)
%             
%             %calculate tidal range % time of high tide
%             eta1_window=eta1{nd}(intv)-mean(eta1{nd}(intv));
%             rmse_eta(node_idx(nd))= sqrt(mean( (eta_tag_fit{i}-eta1_window).^2));
%             
%             
%             % figure(1);plot(time,eta1{nd},'r');hold on
%         end
%         
%         rmse_con(i,rmse_eta<=thresh)=1;
%         rmse_con(i,rmse_eta>thresh)=0;
%         
%         %         figure(figh)
%         %         clf
%         %         patch('Vertices',[fvcom.x,fvcom.y],'Faces',fvcom.tri,'Cdata',rmse_eta,'edgecolor','none','facecolor','interp');
%         %         axis equal;%axis(plot_axis);
%         %         caxis([0 0.3]);
%         %         [a,b]=min(rmse_eta);
%         %         hold on
%         %         %plot(fvcom.x(b),fvcom.y(b),'ro')
%         %         %plot(xt,yt,'ko')
%         %         colorbar()
%         %
%         %         caxis([0,thresh]);
%         
%         
%         H = text(.82e6,1.7e5,['day: ' num2str(i) ' of ' num2str(ndays) ' ']);
%         set(H,'FontSize',16,'Color','k');
        %
        %         %pause(1)
        %
        %         figure(100);plot(time,eta1{find(node_idx==b)}(intv),'r');
        
        
    end
    
end
figure()
plot(day_ampli,'x-')
title('daily amplitude')
fprintf(['days with tidal fits: ' num2str(numel(find(tide==2))) ' of ' num2str(ndays) ' \n'])

% 
% ObsLh=ObsLh.*rmse_con';
% tide(sum(ObsLh)==0)=0;
% 
% disp(sprintf('Saving -> %s.mat <- \n',filename))
% save(filename,'ObsLh','tide','days')
end


