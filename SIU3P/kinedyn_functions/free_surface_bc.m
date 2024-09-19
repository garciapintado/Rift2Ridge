function [KKfs_zz,KKfs_xz,RHSfs_z] = free_surface_bc(MESH,PHYSICS,NUMSCALE,DensEl,alpha,beta,dt)

%==============================================================
% BOUNDARY CONDITIONS FOR A FREE SURFACE
%==============================================================

%______________________________________________________________
%                 Top_nodes x
%     x-----x-----x-----x-----x-----x-----x-----x-----x Surface
%                  \ els_top / \ els_top /
%                   \       /   \       /
%                    o     o     o     o
%                     \   /       \   /
%                      \ /         \ /
%                       o           o
%
%______________________________________________________________

is_triangle = ismember(MESH.nnodel,[3 6 7]);
if is_triangle
    if MESH.nnodel==3
        nnod_edge = 2;
    else
        nnod_edge = 3;
    end
else
    if MESH.nnodel==4
        nnod_edge = 2;
    else
        nnod_edge = 3;
    end
end
nods_top = find(ismember(MESH.PointID,MESH.PointID_top));
  % nodes along top boundary incl. corners
els_top  = find(sum(ismember(MESH.EL2NOD,nods_top),1)==nnod_edge);
  % elements with "nnod_edge" nodes at top boundary
nseg     = length(els_top);

% Density at top boundary
if size(DensEl,2)>1
    DensEl = mean(DensEl,2);
end

% Gravitational force at integration point
del_Dens = DensEl(els_top);%-PHYSICS.fs_dens_top;
Fg_el    = PHYSICS.g * NUMSCALE.Bscale * del_Dens;

% Lithostatic pressure gradient
zmax_top   = max(MESH.GCOORD(2,nods_top));
dz_top     = (zmax_top-MESH.GCOORD(2,nods_top));
Pref_top   = dz_top * NUMSCALE.Dens0 * PHYSICS.g;
P_dens_top = dz_top * PHYSICS.fs_dens_top * PHYSICS.g;
P_top      = zeros(MESH.nnod,1);
P_top(nods_top) = NUMSCALE.Bscale * (P_dens_top-Pref_top);

% Setup for integration along top boundary
nip_fs = 3; % Number of integration points.
Ip_fs  = [-sqrt(3/5), 0, sqrt(3/5)]; % 1D integration points
Ipw_fs = [5/9; 8/9; 5/9];            % 1D integration weights
[N_fs] = shp_line_int(Ip_fs, nnod_edge); % 1D shape functions *SUBFUNCTION*

% Allocate storage for matrices and rhs vector
FSi_zz = zeros(nnod_edge*nnod_edge,nseg,'uint32');
  % Matrix of indices i for the free surface elements in the
  % vertical component.
FSj_zz = zeros(nnod_edge*nnod_edge,nseg,'uint32');
  % Matrix of indices j for the free surface elements in the
  % vertical component.
FSv_zz = zeros(nnod_edge*nnod_edge,nseg);
  % Matrix of values for the correction of the matrix KK for
  % the free surface in the vertical component.
if beta>0
    FSi_xz = zeros(nnod_edge*nnod_edge,nseg,'uint32');
        % Matrix of indices i for the free surface elements in the
        % x component.
    FSj_xz = zeros(nnod_edge*nnod_edge,nseg,'uint32');
        % Matrix of indices j for the free surface elements in the
        % x component. i.e. (3*3, elements on top) 
    FSv_xz = zeros(nnod_edge*nnod_edge,nseg);
        % Matrix of values for the correction of the matrix A for
        % the free surface in the x component.
end
RHS_i = zeros(nnod_edge,nseg);
RHS_v = zeros(nnod_edge,nseg);

% Loop over segments along top boundary
for iseg = 1:nseg
    el            = els_top(iseg);
    nods_top_el   = double(intersect(MESH.EL2NOD(:,el)',nods_top));
    [Xsort, perm] = sort(MESH.GCOORD(1,nods_top_el));
    nods_top_el   = nods_top_el(perm);
    dx_el         = max(Xsort)-min(Xsort);

    % Form 2x2 or 3x3 matrix for integrating along the top boundary
    NN_integ_e = zeros(nnod_edge);
    for ip_fs = 1:nip_fs % Integration loop            
        NN_integ_e = NN_integ_e + N_fs{ip_fs}*N_fs{ip_fs}'*...
            Ipw_fs(ip_fs)*dx_el/2; % Ni*Nj*w*dx
        % /2 because master 1D line element has length 2
    end
    
    % Calculate the force terms resulting from surface topogarphy
    FSrhs_e = NN_integ_e * P_top(nods_top_el);
    
    % Calculate matrices for the correction of K for each element.
    FSv_zz_e = NN_integ_e*Fg_el(iseg)*alpha*dt;
        % Afs for y = inc_density*g*alpha*dt*integral(Ni*Nj*dx)
    
    % Calculate matrices for the correction of r.h.s for each element.
    if beta>0
        h_fs     = MESH.GCOORD(2,nods_top_el); % Heights of the nodes.
        dhdx     = (h_fs(end)-h_fs(1))/dx_el; % Slope at each node.
        FSv_xz_e = NN_integ_e*Fg_el(iseg)*alpha*beta*-(dhdx)*dt;
    end
    
    % Put element FS_zz into storage
    idof           = 2*ones(nnod_edge,1) * nods_top_el;
    jdof           = idof';
    FSi_zz(:,iseg) = idof(:);
    FSj_zz(:,iseg) = jdof(:);
    FSv_zz(:,iseg) = FSv_zz_e(:);
    
    % Put element FSrhs into storage
    RHS_i(:,iseg)  = 2*nods_top_el';
    RHS_v(:,iseg)  = FSrhs_e(:);
    
    % Put element FS_xz into storage (idof is the same as in FS_zz)
    if beta>0
        FSi_xz(:,iseg) = idof(:);
        FSj_xz(:,iseg) = jdof(:)-1; % minus 1 because x-component
        FSv_xz(:,iseg) = FSv_xz_e(:);
    end
end

% % Remove upper triangular part before sparse command because stiffness
% % matrix only contains lower + diagonal part
% switch nnod_edge
%     case 2
%         ind_rm           = 3;
%         FSv_zz(ind_rm,:) = [];
%         FSi_zz(ind_rm,:) = [];
%         FSj_zz(ind_rm,:) = [];
%     case 3
%         %
%         %  / k1  k4  k7 \   Remove           / k1         \
%         % |  k2  k5  k8  |  non-diagonal => |  k2  k5      |
%         %  \ k3  k6  k9 /   terms            \ k3  k6  k9 /
%         ind_rm           = [4 7 8];
%         FSv_zz(ind_rm,:) = [];
%         FSi_zz(ind_rm,:) = [];
%         FSj_zz(ind_rm,:) = [];
%     otherwise
%         error(' Only 2 or 3 nodes per element can be on the boundary.');
% end

nUdof   = 2*max(MESH.EL2NOD(:));
KKfs_zz = sparse3(FSi_zz, FSj_zz, FSv_zz, nUdof, nUdof); 
  % sparsification of FSv_zz.

if beta>0
    KKfs_xz = sparse3(FSi_xz, FSj_xz, FSv_xz, nUdof, nUdof); 
      % sparsification of FSv_xz
else
    KKfs_xz = [];
end

% Assembly of topography force vector
nnod    = size(MESH.GCOORD,2);
RHSfs_z = accumarray(RHS_i(:),RHS_v(:),[2*nnod 1]);
RHSfs_z(:)=0;

end % END OF FUNCTION free_surface_bc

% #########################################################################
%                               SUBFUNCTIONS
% #########################################################################

function  [N] = shp_line_int(ipx,nnodel)

nip  = length(ipx);
N    = cell(nip,1);
for i=1:nip
    eta = ipx(i);
    switch nnodel
        case 2
            SHP = [0.5*(1 - eta);
                   0.5*(1 + eta)];
        case 3
            SHP = [eta*(eta-1)/2;
                   1-eta^2;
                   eta*(eta+1)/2];
    end
    N{i} = SHP;
end

end % END OF SUBFUNCTION shp_line_int