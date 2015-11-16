% cliu
% export daily likelihood plots
function plot_likelihood(fish_no,plot_mpt)
if ~exist('plot_mpt'), plot_mpt = 0; end

filename=['ObsLh',num2str(fish_no),'.mat'];
load(filename)
filename=['mpt',num2str(fish_no),'.mat'];
load(filename)
days = mpt.time;


global fvcom_tidaldb
load(fvcom_tidaldb)

dir_name=[num2str(fish_no), '_out'];
if ~exist(dir_name, 'dir')
    mkdir(dir_name);
end

ndays=size(ObsLh,1);

if plot_mpt==1
    [mpt_x,mpt_y] = my_project(mpt.long,mpt.lat,'forward');
end

% make plots
for d=1:ndays
    
    H2 = figure(2);clf;hold on
    set(H2,'Position', [100, 100, 1024, 768]);
    patch('Vertices',[fvcom.x,fvcom.y],'Faces',fvcom.tri,'Cdata',ObsLh(d,:),'edgecolor','none','facecolor','interp');
    H = text(6.1959e5,2.0322e5,['day: ' num2str(d) ' of ' num2str(ndays),'  ', datestr(days(d),'mmm dd yyyy')]);
    if plot_mpt==1
        plot(mpt_x,mpt_y,'w-')
    end
    export_fig([dir_name,'/likelihood_',num2str(d,'%04d'),'.png']);
    
    
end

end