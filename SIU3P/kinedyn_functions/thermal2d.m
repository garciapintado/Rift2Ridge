function T = thermal2d(T, VAR, MESH, SETTINGS, PHYSICS, NUMSCALE, TBC, dt)
% Usage: T = thermal2d(T, VAR, MESH, SETTINGS, PHYSICS, NUMSCALE, TBC, dt)
% 
% Purpose: Finite element thermal solver, cartesian coordinates
%
% Input:
% 1.  T        : [colvector] : temperature field before heat conduction
% 2.  VAR      : [structure] : major variable fields (Dens, Cp, dQdt, etc)
% 3.  MESH     : [structure] : FE mesh parameters
% 4.  SETTINGS : [structure] : model parameters
% 5.  PHYSICS  : [structure] : physical properties
% 6.  NUMSCALE : [structure] : numerical scaling parameters
% 7.  TBC      : [structure] : boundary conditions on diffusion problem
% 8.  dt       : [scalar]    : time over which is diffused
%
% Output:
%   T        : [colvector] : temperature field after heat conduction
%
% Part of M2TRI - 2D FINITE ELEMENT MANTLE CONVECTION CODE
%
% Based on MILAMIN: MATLAB-based FEM solver for large problems, Version 1.0
% Copyright (C) 2007, M. Dabrowski, M. Krotkiewski, D.W. Schmid
% University of Oslo, PHYSICS of Geological Processes
% http://milamin.org
% See License file for terms of use.
%
% Email contact: jhasenclever@geomar.de
% For numerical methods see online Ph.D. thesis of J.Hasenclever
% (http://www.sub.uni-hamburg.de/opus/volltexte/2010/4873/)
%
% JH Aug 2015
%

%==========================================================================
% MODEL INFO
%==========================================================================

if isfield(SETTINGS,'nelblk')
  nelblk = SETTINGS.nelblk;
else
  nelblk = 1200;  
end
if isfield(SETTINGS,'nnodel')
  nnodel = SETTINGS.nnodel;
else
  nnodel = 6;  % default kinedyn
end

nvert    = 3; % triangle has three vertex nodes
% [EL2NOD,PhaseID] = trimesh_p2_to_p1(MESH.EL2NOD(1:6,:),MESH.PhaseID, SETTINGS.MODEL);

EL2NOD   = MESH.EL2NOD(1:nnodel,:);

%PhaseID  = MESH.PhaseID;                                 % JGP note: commented as not used 

nel      = size(EL2NOD,2);
nnod     = max(EL2NOD(:));
nip      = 6;
ndim     = 2;

GCOORD   = MESH.GCOORD(:,1:nnod);

% % DEBUGGING
% GCOORD      = GCOORD .* NUMSCALE.L0;
% dt          = dt * NUMSCALE.t0;
% NUMSCALE.L0 = 1;
% NUMSCALE.t0 = 1;
% % DEBUGGING

%==========================================================================
% PREPARE INTEGRATION POINTS & DERIVATIVES wrt LOCAL COORDINATES
%==========================================================================
[x_ip,w_ip] = ip_triangle(nip);
if (size(x_ip,1) > size(x_ip,2))                                           % as nip > 2, then we are using rift2ridge2D' ip_triangle()
  x_ip = x_ip';
end
[N,dNds]    = sf_dsf_tri367(x_ip,nnodel,'cell');
[N3,dN3ds]  = sf_dsf_tri367(x_ip,3,'cell');
dN3ds       = dN3ds{1}';

if isfield(SETTINGS,'MODEL')
  if SETTINGS.MODEL == "rift2ridge2D"
     for ip=1:nip
        N{ip} = N{ip}(SETTINGS.ki2mi_lnodes(1:nnodel));
        dNds{ip} = dNds{ip}(:,SETTINGS.ki2mi_lnodes(1:nnodel));
     end
  end
end
    
%==========================================================================
% DECLARE VARIABLES (ALLOCATE MEMORY)
%==========================================================================
K_all     = zeros(nnodel*(nnodel+1)/2,nel);
Rhs_all   = zeros(nnodel,nel);

%==========================================================================
% BLOCKING PARAMETERS (nelblk must be < nel)
%==========================================================================
nelblk    = min(nel, nelblk);
nblk      = ceil(nel/nelblk);
il        = 1;
iu        = nelblk;

%==========================================================================
% BLOCK LOOP - MATRIX COMPUTATION
%==========================================================================
for ib = 1:nblk
    
    %======================================================================
    % CALCULATE JACOBIAN, ITS DETERMINANT AND INVERSE
    %======================================================================
    % NOTE: For triangular elements with non-curved edges the Jacobian is
    %       the same for all integration points (i.e. calculated once
    %       before the integration loop). Further, linear 3-node shape are
    %       sufficient to calculate the Jacobian.
    els_blk    = il:iu;
    ECOORD_x   = reshape( GCOORD(1,EL2NOD(1:nvert,els_blk)), nvert, nelblk );
    ECOORD_z   = reshape( GCOORD(2,EL2NOD(1:nvert,els_blk)), nvert, nelblk );
    Jx         = (dN3ds' * ECOORD_x)';
    Jz         = (dN3ds' * ECOORD_z)';
    detJ       = Jx(:,1).*Jz(:,2) - Jx(:,2).*Jz(:,1);                      % [nelblk, 1] 
    if any(detJ<0)
        error('negative Jacobi')
    end
    invdetJ    = 1./detJ;
    invJx      = zeros(nelblk, ndim);
    invJz      = zeros(nelblk, ndim);
    invJx(:,1) = +Jz(:,2).*invdetJ;
    invJx(:,2) = -Jz(:,1).*invdetJ;
    invJz(:,1) = -Jx(:,2).*invdetJ;
    invJz(:,2) = +Jx(:,1).*invdetJ;
    
    %======================================================================
    % INTEGRATION LOOP
    %======================================================================
    T_blk     = T(EL2NOD(:,els_blk));                                      % [nnodel, nelblk]
    K_blk     = zeros(nelblk, nnodel*(nnodel+1)/2);                         % [nelblk, sum_{i=1}^nnodel i];
    Rhs_blk   = zeros(nelblk, nnodel);                                     % [nelblk, nnodel]
    
    for ip=1:nip
        %==================================================================
        % LOAD SHAPE FUNCTIONS DERIVATIVES FOR INTEGRATION POINT
        %==================================================================
        Ni        = N{ip};
        dNids     = dNds{ip};
        
        %==================================================================
        % PROPERTIES OF ELEMENTS AT ip-TH EVALUATION POINT
        %==================================================================
%         Cond_blk  = PHYSICS.K(PhaseID(els_blk));
%         RhoCp     = PHYSICS.Dens.*PHYSICS.Cp;
%         RhoCp_blk = RhoCp(PhaseID(els_blk));
%         dQdt_blk  = PHYSICS.Hp(PhaseID(els_blk));

        Cond_blk  = PHYSICS.EL_IP.K(els_blk,ip);
        RhoCp_blk = PHYSICS.EL_IP.Dens(els_blk,ip).*PHYSICS.EL_IP.Cp(els_blk,ip);
        dQdt_blk  = PHYSICS.EL_IP.Hp(els_blk,ip);

        if isfield(VAR,'Hs') && isfield(VAR,'Hs_ip')
          error('thermal2d:: heat source input defined at both vertices and integration points')
        end
        if isfield(VAR,'Hs')                                                 % kinedyn default: Hs defined at triangle vertices
            if size(VAR.Hs,2)~=nvert
                error('VAR.Hs must have size nel x 3');
            end
            Hs_blk   = VAR.Hs(els_blk,:) * N3{ip};
            dQdt_blk = dQdt_blk + Hs_blk;
        end
        if isfield(VAR,'Hs_ip')                                              % rift2ridge2D default: Hs externally defined at 6 integration points 
            if size(VAR.Hs_ip,2)~=nip
                error('VAR.Hs_ip must have size nel x 6');
            end
            Hs_blk   = VAR.Hs_ip(els_blk,ip);
            dQdt_blk = dQdt_blk + Hs_blk;
        end
        
        %==================================================================
        % DERIVATIVES wrt GLOBAL COORDINATES
        %==================================================================
        dNdx     = invJx*dNids;
        dNdz     = invJz*dNids;
        
        %==================================================================
        % NUMERICAL INTEGRATION OF ELEMENT MATRICES
        %==================================================================
        weight    = w_ip(ip)*detJ;
        
        %==================================================================
        % 3-node triangle allows lumping to improve stability
        %==================================================================
        
        NxN = Ni*Ni';
        if nnodel==3                                                       % but this does not happen here
            NxN = diag(sum(Ni*Ni',2)); % lumping
        end    
        
        % EXTENDED STIFFNESS MATRIX
        indx = 1;
        for i = 1:nnodel
            for j = i:nnodel
                K_blk(:,indx) = K_blk(:,indx) ...
                    + (NUMSCALE.t0/NUMSCALE.L0^2) ...
                      * dt*(dNdx(:,i).*dNdx(:,j) + dNdz(:,i).*dNdz(:,j)) ...             % stiffness matrix  
                      .*(weight.*Cond_blk(:)) ...
                    + NxN(i,j).*weight.*RhoCp_blk(:);                                    % mass matrix
                indx = indx + 1;
            end
        end

        % RIGHT HAND SIDE
        Rhs_blk = Rhs_blk ...
            + (NxN*(T_blk.*(weight.*RhoCp_blk(:)*ones(1,nnodel))'))'...
            + (Ni*(NUMSCALE.t0*dt*dQdt_blk(:).*weight)')';
    end
    
    %==============================================================
    % WRITE DATA INTO GLOBAL STORAGE
    %==============================================================
    K_all(:,els_blk)	 = K_blk';
    Rhs_all(:,els_blk) = Rhs_blk';
    
    %===========================================
    % READJUST START, END AND SIZE OF NEXT BLOCK
    %===========================================
    il = iu + 1;
    if ib == nblk-1
        nelblk  = nel-iu;   % Account for different number of elements in last block
    end
    iu = iu + nelblk;
end

%==========================================================================
% CREATE TRIPLET FORMAT INDICES
%==========================================================================
indx_j = repmat(1:nnodel,nnodel,1); indx_i = indx_j';
indx_i = tril(indx_i); indx_i = indx_i(:); indx_i = indx_i(indx_i>0);
indx_j = tril(indx_j); indx_j = indx_j(:); indx_j = indx_j(indx_j>0);

K_i = int32(EL2NOD(indx_i,:)); K_i = K_i(:);
K_j = int32(EL2NOD(indx_j,:)); K_j = K_j(:);

indx       = K_i < K_j;
tmp        = K_j(indx);
K_j(indx)  = K_i(indx);
K_i(indx)  = tmp;
clear indx tmp

%==========================================================================
% CONVERT TRIPLET DATA TO SPARSE MATRIX
%==========================================================================
K      = sparse3(K_i, K_j, K_all(:));
Rhs    = accumarray(EL2NOD(:), Rhs_all(:));
clear K_i K_j K_all Rhs_all;

%==========================================================================
% BOUNDARY CONDITIONS
%==========================================================================
if isfield(TBC,'ifix') && ~isempty(TBC.ifix)
    % Move prescribed temperatures to rhs
    ifree = 1:nnod; ifree(TBC.ifix) = [];
    Kbc   = K(:,TBC.ifix) + cs_transpose(K(TBC.ifix,:));
    Rhs   = Rhs - Kbc*TBC.vfix(:); clear Kbc
    T(TBC.ifix) = TBC.vfix;
else
    ifree = 1:nnod;
end

% Add heat flux to Rhs
if isfield(TBC,'NN_INTEG')
    Rhs = Rhs + (NUMSCALE.t0/NUMSCALE.L0^2)*dt*TBC.NN_INTEG*TBC.valHF(:);
end

%==========================================================================
% SOLVE SYSTEM OF EQUATIONS
%==========================================================================
[L,perm] = cholesky_factorize(K(ifree,ifree));
T(ifree) = cholesky_solve(Rhs(ifree),L,perm);

end % END OF FUNCTION thermal2d
