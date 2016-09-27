function [x,y]=cs2cs(u, v, prj4_params)
%CS2CS - cartographic coordinate system filter
%  This is a wrapper-function for cs2cs.
%
% Syntax:  [x,y] = cs2cs(u, v, prj4_params)
% 
% Inputs: 
%    u - vector with horizontal input coordinates
%    v - vector with vertical input coordinates
%    prj4_params - parameters of cs2cs
% 
% Outputs:
%    x - vector with horicontal output coordinates
%    y - vector with vertical output coordinates
%

% Author: Erwin Nindl
% Email: nine-cs2cs-2013@wirdorange.org
% Website: https://github.com/nine/matlab-cs2cs
% April 2013; Last revision: 12-May-2004

  % call cs2cs from matlab
  % inspired by
  % http://marinescience.wiki.otago.ac.nz/Cs2csFromMatlab
  %------------------------------------------------------------------------
  
  % input validation
  %------------------------------------------------------------------------  
  if nargin<3
    error('please provide all params');
  end
  %if ~isvector(u)
      [nx, ny]=size(u);
  %end
  u = u(:);
  v = v(:);
  if length(u)~=length(v)
    error('input coordinates must have the same length');
  end
  
  if ~isempty(u)
    % cs2cs: operating system dependent stuff
    %----------------------------------------------------------------------
    if ispc() % MS windows
      f = filesep();
      [pathstr,~,~] =  fileparts(mfilename('fullpath'));
      proj_path     = [pathstr f 'util' f 'proj' f 'bin'];
      proj_lib_path = [pathstr f 'util' f 'proj' f 'nad'];
      clear pathstr;
      clear f;
      setenv('PATH', [getenv('PATH') ';' proj_path]);
      setenv('PROJ_LIB', proj_lib_path);
    else
      if isunix() || ismac()
        if unix('which cs2cs')
          error('binary of cs2cs not found in path');
        end
      else % unknown OS
        error('operating system not supported');
      end
    end
  
    % call cs2cs
    %----------------------------------------------------------------------
    tmp_file1 = tempname(); % infile
    tmp_file2 = tempname(); % outfile
    dlmwrite(tmp_file1, [u, v], 'delimiter', ' ', 'precision', '%.6f');

    [status,result] = system(['cs2cs ' prj4_params ' < ' tmp_file1 ' > ' tmp_file2]);
    if status~=0
      error(['Calling cs2cs: ' result]);
    end

    M = dlmread(tmp_file2);
    x = M(:,1);
    y = M(:,2);
    x = reshape(x, nx, ny);
    y = reshape(y, nx, ny);

    % cleanup
    %----------------------------------------------------------------------
    delete(tmp_file1);
    delete(tmp_file2);
  else
    x = [];
    y = [];
  end

end %eof
