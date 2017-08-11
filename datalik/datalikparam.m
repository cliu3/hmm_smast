%DATALIKPARAM SCRIPT
%   Script used by datalikelihhood for defining the
%   variance parameters to used in the computation.
%
%   Date: 7/12 - 2007, ver. 0.61
%   HMM geolocation toolbox, IMM and DIFRES



%% White noise stdev %%
E = 0.2;

%% Parameters for AR(1) term %%
% Time varying variance
epsilon = td.rmse;                        
% Time constant variance
%epsilon = 0.4 * ones(1,length(td.rmse)); 
% "Forgetting coefficient"
lambda  = 0.05^(1/40);                    

%% Variables for the cos term %%
ptime = 360/(db.freq(1)/(60/td.dt*24)*180/pi);
[a b] = meshgrid(1:td.tideFL); c=abs(a-b);
cospattern = cos(2*pi/ptime * c);


%% Defining the variance parameters %%
% Cos (tidalroughness)
if isfield(db,'tidalro')
    s_e =        db.tidalro.^2+eps^20;
elseif isfield(db,'tidro')
    s_e =        db.tidro.^2+eps^20;
else
    error('Tidal roughness is missing in tidaldb.mat!')
end
% Roughness
if isfield(db,'bathro')
    s_eta =      db.bathro+eps^20;  
elseif isfield(db,'rough')
    s_eta =      db.rough+eps^20;  
else
    error('Bahtymetry roughness is missing in tidaldb.mat!')
end
s_eta_tid =  10^2;
% White noise
s_E =        E^2 * eye(td.tideFL);               