function interp_regular(fish_no,path_to_tags)
% interp_regular interpolate likelihood distribution from ObsLh to the
% regular HMM computational grid.
% This function is invoked when tidal exclusion is not necessary and
% disabled.
%
% - Chang Liu
%
global fvcom_tidaldb
load(fvcom_tidaldb)
load tidaldb.mat

tag_name=[num2str(fish_no),'_raw'];
load([path_to_tags tag_name]);
tagno=[num2str(fish_no),'_',tag.tag_id];
% load ObsLh
filename = ['ObsLh' tagno '.mat'];
disp(sprintf('Loading %s...\n',filename))
load(filename)
% interpolate onto regular grid
filename = ['datalikelihood' tagno '.mat'];
fprintf('Interpolating likelihood onto regular grid ... \n');
disp(sprintf('Loading %s...\n',filename))
load(filename)
%[fvcom_lon,fvcom_lat]=my_project(fvcom.x,fvcom.y,'inverse');

ndays = size(ObsLh,1);
for i=1:ndays
    F = TriScatteredInterp(fvcom.lon,fvcom.lat, ObsLh(i,:)');
    %F = TriScatteredInterp(fvcom_lon,fvcom_lat, ObsLh(i,:)');
    TempLh=F(db.long,db.lat);
    TempLh(db.land)=0;
    TempLh(isnan(TempLh))=0;
    LIK.tide(:,:,i)=TempLh;
end
filename = sprintf('datalikelihood%s',tagno);
disp(sprintf('Saving -> %s.mat <- \n',filename))
save(filename,'LIK')
end