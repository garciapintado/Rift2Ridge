function [Z] = function_sealand(Z,no_nodes,F_loc,MM_loc,De,nexp,dt,kappa, ...
    Dx,sealevel,ksea,kappa_s,kdecay,slim,nip,nnodel)
% [Z] = FUNCTION_SEALAND(Z,NO_NODES,F_LOC,MM_LOC,DE,NEXP,DT,DX,SEALEVEL,
% KSEA,KDEKAY) calculates a new topography Z solving with FEM the 
% sediment transport partial differential equation. This
% model can solve all together linear hill-slope diffusion, fluvial
% sediment transport and submarine sediment transport. Z in the output is
% the topography after the time step, Z in the input is the topography
% before the time step, no_nodes is the number of nodes, F_loc is the
% source and Neuman's boundary conditions added together, MM_loc is a part
% of the Stiffness matrix containing the integration of the multiplication
% of the shape functions evaluated at the integration points, DE is the
% fluvial transport coefficient(c x alpha^nexp), NEXP is the exponent of
% the water flux for fluvial transport, DT is the time step, DX is a vector
% of element sizes, SEALEVEL is the sea level, KSEA is
% the submarine diffusion coefficient and KDECAY is the parameter that
% conditions the decay of KSEA with the water depth. FEM is done here with
% linear shape functions and 2 integration points. This function cannot be
% used at its current state with different shape functions and integration
% points.

% Author: John Armitage, Institut de Physique du Globe de Paris
%
% Edited: Miguel Andres-Martinez, University of Bremen
%         andresma@uni-bremen.de

%==========================================================================
% INITIALIZATION
%==========================================================================
Z0 = Z;
% Calculate number of elements
no_el       = no_nodes-1;
% Nodes per element
no_nodes_el = 2;

% Coordinates X
X = [0; cumsum(Dx)];
XX = X;

% Connectivity matrix
nodes   = [(1:no_nodes-1)' (2:no_nodes)'];
% Initialize diffusivity coefficient (different for each element) 
% ((d1+d2)/2) where d1 is the diffusivity in local node 1 and d2 in local
% node 2
D       = zeros(no_el,1);

% Integration points
if nip==1
    IP_x{1} = 0;
    w{1} = 2;
elseif nip==2
    IP_x{1} = -sqrt(1/3);
    IP_x{2} = sqrt(1/3);
    w{1} = 1;
    w{2} = 1;
elseif nip==3
    IP_x{1} = -sqrt(3/5);
    IP_x{2} = 0;
    IP_x{3} = sqrt(3/5);
    w{1} = 5/9;
    w{2} = 8/9;
    w{3} = 5/9;
end

% Shape functions
if nnodel==2
    for ip = 1:nip
        N{ip} = [1/2*(1-IP_x{ip}) 1/2*(1+IP_x{ip})];
        dNdx{ip} = [-1/2 1/2];
    end
elseif nnodel==3
    for ip = 1:nip
        N{ip} = [1/2*IP_x{ip}.*(IP_x{ip}-1) 1-IP_x{ip}.^2 1/2*IP_x{ip}.*(IP_x{ip}+1)];
        dNdx{ip} = [IP_x{ip}-1 -2*IP_x{ip} IP_x{ip}+1];
    end
    X = [X; (X(1:end-1)+X(2:end))/2];
    nodes = [nodes (no_nodes+(1:no_el))'];
end

%==========================================================================
% CALCULATE THE SLOPE DIRECTION AND WATER FLUX
%==========================================================================
[x,sortx] = sort(X);
z = interp1(XX,Z,x);
dx = diff(x);
% Water flux
Dx_sign = -dx.*(sign(diff(z)));
% diff(Z) explained
% -----------------
%
%                       . Z4 = 5
%                      /          Z       = topography defined in the nodes
%           .Z1 = 1   /           diff(z) = defined in each element
%          / \       . Z3 = 0
%         /   \     /
%              \   /
%               \.Z2 = -1
%
% diff(Z) = [-2     1   5]
%
% sign(diff(Z)) returns a vector of ones with the sign of diff(Z):
%   sign(diff(Z)) = [-1 1 1]

% Find top of the hills and depocenter nodes
indx = [find(diff(diff(z)>0)~=0)+1; length(z)];
% figure(2); clf
% plot(X,Z,'o-')
% hold on
% plot(X(indx),Z(indx),'x')

X_d = zeros(size(x));
ind0 = 1;
for n = 1:length(indx)
    inde            = indx(n);
    dx_sign         = Dx_sign(ind0:inde-1);
    length_river    = sum(dx_sign);
    if all(dx_sign<0)
        x_dn = cumsum([0; dx_sign])-length_river;
    elseif all(dx_sign>0)
        x_dn = cumsum([0; dx_sign]);
    else
        error('Bad calculation of distance to the drainage divide')
    end
    x_dn([1 end]) = 0;
    X_d(ind0:inde) = x_dn;
    ind0 = inde;
end

X_d(sortx) = X_d;
z(sortx) = z;
x(sortx) = x;

Dn = zeros(size(X));
% Fluvial transport diffusivity defined in the nodes
Dn = kappa + De*(X_d.^nexp);
% Submarine transport diffusivity defined in the nodes
% Dn(z<sealevel) = kappa_s+ksea*exp(-kdecay*abs(sealevel-z(z<sealevel))); %+1/365.25/24/3600;
Dn_sea = kappa_s+ksea*exp(-kdecay*abs(sealevel-z));

Dn(z<sealevel) = Dn_sea(z<sealevel);
% % Averaging of the diffusion coefficient
% if avrg_diff
%     Dn = [Dn(1:2); (Dn(1:end-4)+Dn(2:end-3)+Dn(3:end-2)+Dn(4:end-1)+Dn(5:end))/3; Dn(end-1:end)];
% end
D_l = ones(1,no_el);
for n = 1
% Diffusivity in the integration points
D = (Dn(1:end-1)+Dn(2:end))/2;

% Slope limit in the ocean for turbidites
zz = [Z(1:end-1)'; Z(2:end)'];
D_lim = zeros(1,size(zz,2));
for ip = 1:nip
    D_lim = D_lim + 0.5*dNdx{ip}*zz*2./Dx';
end
D_l = D_l.^(D_lim'.*slim.^(-1));
el_sea = (sum(zz)/2)<sealevel;
%D(el_sea) = D(el_sea).*D_lim(el_sea)./(D(el_sea)+D_lim(el_sea));
D(el_sea) = 1e-20/3600/24/365.25;

% Setup system matrix
KM          = zeros(no_nodes,no_nodes);
MM          = zeros(no_nodes,no_nodes);
F           = zeros(no_nodes,1);

% % Plot diffusivity and scaled topography with sealevel
% figure(3); clf
% line([min(x) max(x)],[sealevel sealevel],'LineStyle','--')
% hold on
% line(x,z,'Color',[0.85 0.325 0.098],'Marker','.','MarkerSize',10)
% ax1 = gca; % current axes
% ax1_pos = ax1.Position;
% ax2 = axes('Position',ax1_pos,'XAxisLocation','top','YAxisLocation', ...
%     'right','Color','none');
% linkaxes([ax1 ax2],'x')
% line(x,Dn,'Parent',ax2,'Color','k','Marker','.','MarkerSize',10)
% grid on
% drawnow

for iel=1:no_el
    if nip>0
        KM_loc = zeros(nnodel,nnodel);
        MM_loc = zeros(nnodel,nnodel);
        for ip = 1:nip
            KM_loc = KM_loc + dNdx{ip}'*N{ip}*Dn(nodes(iel,:))*dNdx{ip}*2/Dx(iel)*w{ip};
            MM_loc = MM_loc + N{ip}'*N{ip}/2*w{ip};
        end
    else
        KM_loc       = D(iel)*[1/Dx(iel) -1/Dx(iel); -1/Dx(iel) 1/Dx(iel)];
    end
    
    for i=1:no_nodes_el
        ii                              = nodes(iel,i);
        for j=1:no_nodes_el
            jj                          = nodes(iel,j);
            KM(ii,jj)                   = KM(ii,jj) + KM_loc(i,j);
            MM(ii,jj)                   = MM(ii,jj) + Dx(iel)*MM_loc(i,j); %%%NEW%%%
        end
        F(ii)                           = F(ii)     + F_loc(i);
    end
end

F_tot                                   =       F   + 1/dt*MM*Z0; %%%NEW%%%
KK_TOT                                  = 1/dt*MM   +      KM; %%%NEW%%%

% zero gradient

U                     = chol(KK_TOT);
L                     = U';

Z                     = U\(L\F_tot);

% Dirichlet

% Free                 = 1:no_nodes;
% Free([1 no_nodes])   = [];
% Z                    = zeros(no_nodes,1);
% Z([1 no_nodes])      = [0 0];
% F_tot                = F_tot - KK_TOT*Z;
% 
% U                    = chol(KK_TOT(Free,Free));
% L                    = U';
% 
% Z(Free)              = U\(L\F_tot(Free));

% line(x,Z,'Parent',ax1,'Color','g','Marker','.','MarkerSize',10)
% drawnow
end
end
