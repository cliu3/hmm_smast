function [kern kern3] = makekern2(par)
% MAKEKERN2
% [kern kern3] = makekern2(par)
% Inputs: par
% par must have the field covmat which should be a 2 by 2 covariance matrix
if ~isfield(par,'u'), par.u = [0 0]; end


% Setup initial kernel (kernel base) for x direction
ux = abs(par.u(1));
kbl_x = ceil(ux)*2+1;
kb_x = zeros(1,kbl_x);
inp = [1-ux+floor(ux) ux-floor(ux)];
if inp(2) == 0
    kb_x(kbl_x) = 1;
else
    kb_x(kbl_x-1:kbl_x) = inp;
end
if sign(par.u(1)) == -1, kb_x = fliplr(kb_x); end
kb_var_x = sum(kb_x.*((-(kbl_x-1)/2:(kbl_x-1)/2) - par.u(1)).^2);
rv(1) = par.covmat(1,1) - kb_var_x; % Remaining variance
%rv(1) = 0.25 - kb_var_x; % Remaining variance
%%%%% if rv(1) < 0, error('Diffusion is too small compared to advection in x-dir! (makekern2)'), end
% % mk = [0.25*rv(1) 1-0.5*rv(1) 0.25*rv(1)];
% % kb_x = convn(convn(kb_x,mk),mk);
% % kb_var_x = sum(kb_x.*((-(kbl_x+3)/2:(kbl_x+3)/2) - par.u(1)).^2);

% Setup initial kernel (kernel base) for y direction
uy = abs(par.u(2));
kbl_y = ceil(uy)*2+1;
kb_y = zeros(1,kbl_y);
inp = [1-uy+floor(uy) uy-floor(uy)];
if inp(2) == 0
    kb_y(kbl_y) = 1;
else
    kb_y(kbl_y-1:kbl_y) = inp;
end
if sign(par.u(2)) == -1, kb_y = fliplr(kb_y); end
kb_var_y = sum(kb_y.*((-(kbl_y-1)/2:(kbl_y-1)/2) - par.u(2)).^2);
rv(2) = par.covmat(2,2) - kb_var_y; % Remaining variance
%rv(2) = 0.25 - kb_var_y; % Remaining variance
%%%%% if rv(2) < 0, error('Diffusion is too small compared to advection in y-dir! (makekern2)'), end
% Make sure all schemes have 0.25 variance to begin with since som values
% for advection forces the scheme to a minimum variance of 0.25 (when u =
% 0.5 or 1.5 or 2.5 or....)
% % mk = [0.25*rv(2) 1-0.5*rv(2) 0.25*rv(2)];
% % kb_y = convn(convn(kb_y,mk),mk);
% % kb_var_y = sum(kb_y.*((-(kbl_y+3)/2:(kbl_y+3)/2) - par.u(2)).^2);

N = [1 1];
kb{1} = kb_x;
kb{2} = kb_y;
for j = 1:2 % Cycle the two directions
    %disp('============')
    if rv(j)<0
        D = 0.5*par.covmat(j,j);
        %D = 0.5*(par.covmat(j,j)-0.25); % The 0.25 is to cope with numerical diffusion from advection (see above)
    else
    D = 0.5*rv(j); % The 0.25 is to cope with numerical diffusion from advection (see above)
    end
    n = floor(D/0.5);
    mink = 1;
    endk = 1;
    while (endk > 0.0001 | mink < 0 | (D/n)> 0.4) % max probability in end point of 1D kernel
        if n == 100
            warning('Kernel is large! (n>100)')
        end
        if n == 200
            error('Kernel too large! (n>200)')
        end
        n = n+1;
        kern2 = [D/n 1-2*D/n D/n];
        kern3 = kb{j};
        for i=1:n
            kern3 = conv(kern3,kern2);
        end
        endk = max([kern3(1) kern3(end)]);
        mink = min(kern3);
        %D/n
        %kern3
    end
    % Convolve to take care of advection (increase kernel size)
%     n = n+ua;
%     kern2 = [D/n 1-2*D/n D/n];
%     kern3 = 1;
%     for uadd = 1:n
%         kern3 = conv(kern3,kern2);
%     end
    % put a limit on kernel size
    if length(kern3) > 100
        inds = kern3 < 1e-7;
        if sum(inds) > 0,
            l = sum(inds(1:n));
            r = sum(inds(n+2:end));
            lr = min([l r]);
            kern3 = kern3([lr+1:length(kern3)-lr]);
        end
    end
    % Store kernel
    K{j} = kern3;
    N(j) = n;
%     mink
%     endk
%     D/n
% kern3
end
% N
% K{1}
% K{2}
% % % The largest 1D-kernel determines the size of the 2D-kernel
% % % These calsulations are to make the kernel square in size
% % [L I] = max([length(K{1}) length(K{2})]);
% % if I == 2, 
% %     Is = 1; 
% % elseif I == 1, 
% %     Is = 2; 
% % end
% % 
% % n = N(I);
% % D = 0.5*rv(Is);
% % kern2 = [D/n 1-2*D/n D/n];
% % kern3 = 1;
% % for i=1:n
% %     kern3 = conv(kern3,kern2);
% % end
% % 
% % K{Is} = kern3;
% K{1}
% K{2}
kern = K{2}'*K{1};

if sum(kern(:)) < 0.99, warning('kernel does not sum to one, it sums to %f (makekern2)',sum(kern(:))),end
if min(kern(:)) < 0, warning('kernel has elements that are below zero!'); end
%KS = size(kern);
%if KS(1)~=KS(2), error('kernel is not square! (makekern2)'), end




