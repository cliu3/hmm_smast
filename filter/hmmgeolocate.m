function hmmgeolocate(tagno,mode,viewres,Duser,ext,GEN)
%HMMGEOLOCATE  Obtain geolocation by filtering preprocessed data 
%   HMMGEOLOCATE(TAGNO,MODE,VIEWRES,DUSER,EXT)
%
%   - TAGNO indentifier as string for the tag to geolocate.
%
%     Optional arguments
%
%   - MODE number of behaviour modes to use (1 or 2).
%   default is 2.
%   - VIEWRES plots the marginal distributions consecutively
%   in an animation when the geolocation has finished 
%   successfully (by using the viewdistr function).
%   default is 'on'.
%   - DUSER user defined diffusivity vector e.g DUSER = [10 100].
%   if omitted the diffusivity is estimated with maximum likelihood.
%   (- EXT self estimate behaviour - IS NOT OPERATIONAL YET!)
%
%   DEPENDENCIES - the function needs access to the following files
%
%     tagdataTAGNO.mat
%     datalikelihoodTAGNO.mat
%     tidaldb.mat
%     cmap.mat
%
%   and creates as output the file resultTAGNO.mat in the current folder.
%
%  EXAMPLES
%   HMMGEOLOCATE('2255',2,'on')
%   HMMGEOLOCATE('1432',[],'on',[10 100])
%
%   Date: 22/10 - 2008, ver. 0.55
%   HMM geolocation toolbox, DTU Informatics and DTU Aqua

if nargin < 2 | isempty(mode),  mode = 2; end

if nargin < 3 | isempty(viewres),  viewres = 'on'; end

if nargin < 5 | isempty(ext), ext = false; end


if ext
    if nargin < 4 | isempty(Duser)
        if mode == 1
            Duser.Duser = [60 60];
        else
            Duser.Duser = [10  100];
        end
        Duser.mp = 0.5; % Start guess for p in (1-p)*mode1 + p*mode2, behaviour switching
        Duser.estimate = 1;
    else
        Duser.estimate = 0;
        % It is assumed that no mode probability (mp) is defined
        if length(Duser.Duser) == 1
            disp(sprintf('Using user defined diffusivity, one mode:\nD = %8.4f',Duser.Duser))
            D = Duser.Duser;
            Duser.Duser = [D D];
        elseif length(Duser) == 2
            disp(sprintf('Using user defined diffusivity, two modes:\nD = [%8.4f, %8.4f]',Duser.Duser(1),Duser.Duser(2)))
        end
    end
    if nargin < 6 | isempty(GEN), GEN = false; end
    hmmgeolocate_mode(tagno,mode,viewres,Duser,GEN)
else
    % No behaviour switching and no generator use (28/11-08)
    if nargin < 4 | isempty(Duser)
        if mode == 1
            Duser = [60 60];
        else
            Duser = [10  100];
        end
    else
        if length(Duser) == 1
            disp(sprintf('Using user defined diffusivity, one mode:\nD = %8.4f',Duser))
            Duser = [Duser Duser];
        elseif length(Duser) == 2
            disp(sprintf('Using user defined diffusivity, two modes:\nD = [%8.4f, %8.4f]',Duser(1),Duser(2)))
        end
    end
    % run hmmgeolocate with the defined parameters
    hmmgeolocate1(tagno,mode,viewres,Duser)
end