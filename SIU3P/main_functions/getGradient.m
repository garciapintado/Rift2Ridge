function dV = getGradient(EL2NOD, GCOORD, Vb, BC_ind, BC_val, nelblk)
  % function dV = getGradient(El2NOD, GCOORD, Vb, BC_ind, BC_val, nelblk)
  % Get gradient of a field at mesh nodes, possibly with Cauchy-type boundary conditions
  %
  % EL2NOD :: [nnodel,nel]
  % GCOORD :: [2,nnod]
  % Vb     :: [nnod,1]
  % BC_ind :: OPT, set to [] for no imposed boundary condition
  % BC_val :: OPT, set to [] for no imposed boundary condition
  % nelblk :: OPT, INTEGER, maximum block size for matrix processing
  %
  % 2020-06-19 Javier García-Pintado, MARUM

  %==========================================================================
  % MODEL INFO - HYDROTHERMAL DOMAIN
  %==========================================================================
  
  if nargin < 4
      BC_ind = [];
  end
 
  if nargin < 5
      BC_val = [];
  end
  
  if nargin < 6
      nelblk = 5000;
  end
  
  nnod     = length(unique(EL2NOD));
  [nnodel, nel] = size(EL2NOD);

  if length(GCOORD) ~= nnod
      error("getGradient: GCOORD and EL2NOD not compliant")   
  end

  if isempty(BC_val)
      BC_val = zeros(size(BC_ind));
  end
  
  Vb = Vb(:);
  if size(Vb,1) ~= nnod
      error("getGradient: GCOORD and Vb not compliant")   
  end
  
  %==========================================================================
  % CONSTANTS
  %==========================================================================
  ndim = 2;
  nip = nnodel;                                                    
  nvert = 3;                                                               % triangle vertices

  %==========================================================================
  % PREPARE INTEGRATION POINTS & DERIVATIVES wrt LOCAL COORDINATES
  %==========================================================================
  [x_ip,w_ip] = ip_triangle(nip);
  [N, dNds]   = shp_deriv_triangle(x_ip, nnodel);                            % N [nip] cell array, with [nnodel,1] size per cell. Shape function evaluated at integration points 
  Nm = cell2mat(N');                                                         % [nnodel, nip] 
  [~,dN3ds]   = sf_dsf_tri367(x_ip',3,'cell');
  dN3ds       = dN3ds{1}';                                                   % [3,2]

  %==========================================================================
  % DECLARE VARIABLES (ALLOCATE MEMORY)
  %==========================================================================
  K_all     = zeros(nnodel*(nnodel+1)/2,nel);
  Rhs_all_x = zeros(nnodel,nel);
  Rhs_all_z = zeros(nnodel,nel);
  %==========================================================================
  % BLOCKING PARAMETERS
  %==========================================================================
  nelblk = min(nel, nelblk);
  nblk   = ceil(nel/nelblk);
  il     = 1;
  iu     = nelblk;

  %==================================================================
  % i) BLOCK LOOP - MATRIX COMPUTATION
  %==================================================================
  for ib = 1:nblk
      
    els_blk    = il:iu;
    
    %==============================================================
    % i.1 GET JACOBIAN, J, J^{-1}, and their determinants
    %==============================================================
    % In these model we use linear coordinate transformation, which is 
    % subparametric with respect to the trial function (i.e. quadratic 
    % shape functions for Temperature interpolation and Galerkin weights).
    % That is, physical space triangles have straight edges, and this 
    % => constant |J| throughout each element => no need to calculate it at
    % each integration point
    
    ECOORD_x   = reshape( GCOORD(1,EL2NOD(1:nvert,els_blk)), nvert, nelblk );  % [nvert, nelblk]
    ECOORD_z   = reshape( GCOORD(2,EL2NOD(1:nvert,els_blk)), nvert, nelblk );  % [nvert, nelblk]
    Jx         = ECOORD_x' * dN3ds;                                            % [nelblk, 2] each row : [\frac{\partial x}{\partial \xi} \frac{\partial x}{\partial \eta}]
    Jz         = ECOORD_z' * dN3ds;                                            % [nelblk, 2] each row : [\frac{\partial y}{\partial \xi} \frac{\partial y}{\partial \eta}]
    detJ       = Jx(:,1).*Jz(:,2) - Jz(:,1).*Jx(:,2);                          % [nelblk, 1] 
    if any(detJ<0)
        error('negative |J|')
    end
    invdetJ    = 1./detJ;
    invJx      = zeros(nelblk, ndim);
    invJz      = zeros(nelblk, ndim);
    invJx(:,1) = +Jz(:,2).*invdetJ;
    invJx(:,2) = -Jz(:,1).*invdetJ;
    invJz(:,1) = -Jx(:,2).*invdetJ;
    invJz(:,2) = +Jx(:,1).*invdetJ;
       
    Vb_bl =  Vb(EL2NOD(:,els_blk));                                        % [nnodel,nelblk] background field
    
    %==============================================================
    % iii) INTEGRATION LOOP
    %==============================================================
   
    K_blk     = zeros(nnodel*(nnodel+1)/2, nelblk);                        % [nelblk, sum_{i=1}^nnodel i]; block element contribution to the [nnod,nnod] extended stiffness matrix 
    Rhs_blk_x = zeros(nnodel, nelblk);                                     % [nnodel, nelblk];             block element contribution to the [nnod,1] right hand side 
    Rhs_blk_z = zeros(nnodel, nelblk);

    for ip=1:nip
        %==========================================================
        % LOAD canonical SHAPE FUNCTIONS & DERIVATIVES evaluated @ this integration point
        %==========================================================
        Ni          = N{ip};                                               % [nnodel,1] 
        dNids       = dNds{ip};                                            % [nnodel,2] local coordinate derivatives
        
        NxN = Ni*Ni';
        if nnodel==3
            NxN = diag(sum(Ni*Ni',2));                                     % lumping -> improves stability
        end
        
        dNdx        = dNids * invJx';                                      % [nnodel,nelblk] global coordinate derivatives Nx  <-  Nxz = [Nx | Nz]
        dNdz        = dNids * invJz';                                      % [nnodel,nelblk] global coordinate derivatives Nz
        
        %==========================================================
        % vii) NUMERICAL INTEGRATION OF ELEMENT MATRICES
        %==========================================================
        weight_ip   = reshape(w_ip(ip) .* detJ, 1, nelblk);                % [1, nelblk]

        indx = 1;
        for i = 1:nnodel
            for j = i:nnodel                                               % [nnodel,1] block assignations
                K_blk(indx,:) = K_blk(indx,:) + NxN(i,j) .* weight_ip;     % ip block contribution to mass matrix
                indx = indx + 1;
            end
        end
        
        %RIGHT HAND SIDE
        Rhs_blk_x = Rhs_blk_x + Ni * (sum(dNdx .* Vb_bl) .* weight_ip);    % [nnodel,nelblk]  =  [nnodel,1]*[1,nelblk]
        Rhs_blk_z = Rhs_blk_z + Ni * (sum(dNdz .* Vb_bl) .* weight_ip);    % [nnodel,nelblk]
    end % for ip
    %==============================================================
    % WRITE DATA INTO GLOBAL STORAGE
    %==============================================================
    K_all(:,il:iu)	    = K_blk;
    Rhs_all_x(:,il:iu)  = Rhs_blk_x;
    Rhs_all_z(:,il:iu)  = Rhs_blk_z;
    
    il  = il + nelblk;
    if ib==nblk-1
        nelblk = nel-iu;                                                   % size and indices of next block
    end
    iu  = iu + nelblk;
  end % for ib


  %==========================================================================
  % ix) CREATE TRIPLET FORMAT INDICES
  %==========================================================================
  indx_j = repmat(1:nnodel,nnodel,1); indx_i = indx_j';
  indx_i = tril(indx_i); indx_i = indx_i(:); indx_i = indx_i(indx_i>0);
  indx_j = tril(indx_j); indx_j = indx_j(:); indx_j = indx_j(indx_j>0);

  K_i = EL2NOD(indx_i,:); K_i = K_i(:);
  K_j = EL2NOD(indx_j,:); K_j = K_j(:);

  indx       = K_i < K_j;
  tmp        = K_j(indx);
  K_j(indx)  = K_i(indx);
  K_i(indx)  = tmp;
  clear indx tmp

  %==========================================================================
  % x) CONVERT TRIPLET DATA TO SPARSE MATRIX
  %==========================================================================
  K_all  = K_all(:);
  K      = sparse2(K_i, K_j, K_all);
  Rhs_x    = accumarray(EL2NOD(:), Rhs_all_x(:));
  Rhs_z    = accumarray(EL2NOD(:), Rhs_all_z(:));
  clear K_i K_j K_all Rhs_all_x Rhs_all_z;

  dVx = zeros(nnod,1); 
  dVz = zeros(nnod,1);
  %==========================================================================
  % BOUNDARY CONDITIONS
  %==========================================================================
  PBC = [];
  if ~isempty(BC_ind)
      PBC.nod = BC_ind;                                                    % [nbdi,1]
      PBC.val = BC_val;                                                    % [nbdi,1] fluid pressure at boundary conditions
      Free    = setdiff(1:nnod,PBC.nod);                                   % [nnod-nbdi] replaces Free = 1:nnod; Free(PBC.nod)= [];
      % Move prescribed values to rhs
  
      K_L        = tril(K,-1);  % [nnod,nnod] K without diagonal
      TMP        = K(:,PBC.nod) + cs_transpose(K_L(PBC.nod,:)); clear K_L        % [nnod,nbdi]
      Rhs_x      = Rhs_x - TMP*PBC.val;
      Rhs_z      = Rhs_z - TMP*PBC.val;
      clear TMP

      K          = K(Free,Free);
      dVx(PBC.nod) = PBC.val;
      dVz(PBC.nod) = PBC.val;   
  else
      Free = 1:nnod;
  end
  
  %==========================================================================
  % Solve equation with suitesparse cholesky direct solver
  %==========================================================================
  perm = amd(K); % amd reorder
  %==========================================================================
  % FACTORIZATION - ideally L = lchol(K, perm)
  %==========================================================================
  K = cs_transpose(K);
  K = cs_symperm(K,perm);
  K = cs_transpose(K);
  L = lchol(K);
  %==========================================================================
  % SOLVE
  %==========================================================================
  dVx(Free(perm)) = cs_ltsolve(L,cs_lsolve(L,Rhs_x(Free(perm))));
  dVz(Free(perm)) = cs_ltsolve(L,cs_lsolve(L,Rhs_z(Free(perm))));
  dV = [dVx';dVz'];

  % example call and demonstration of error in non-solved approximation
  % dPfb = getGradient(HMESH.EL2HNO, HMESH.HCOORD, HVAR.Pfb, [], [], SOLVER.nelblo);
  % figure(); plot_tF(HVAR.Pfb,  HMESH.HCOORD, HMESH.EL2HNO, [], "Pfb [Pa]",0);
  % figure(); plot_tF(dPfb(1,:),  HMESH.HCOORD, HMESH.EL2HNO, [-10.E04 5.E04], "\partial Pfb / \partial x [Pa m^{-1}]",0);
  % figure(); plot_tF(dPfb(2,:),  HMESH.HCOORD, HMESH.EL2HNO, [-1.5E05 1.E05], "\partial Pfb / \partial y [Pa m^{-1}]",0);
  % [dPfbx,dPfby] = gradient_field(HMESH.HCOORD, HVAR.Pfb, HMESH.EL2HNO, 6); 
  % figure(); plot_tF(dPfbx,  HMESH.HCOORD, HMESH.EL2HNO, [-10.E04 5.E04], "Miguel: \partial Pfb / \partial x [Pa m^{-1}]",0);
  % figure(); plot_tF(dPfby,  HMESH.HCOORD, HMESH.EL2HNO, [-1.5E05 1.E05], "Miguel: \partial Pfb / \partial y [Pa m^{-1}]",0);
  % hnnod6 = numel(dPfby);
  % figure(); plot_tF(dPfbx-dPfb(1,1:hnnod6),HMESH.HCOORD, HMESH.EL2HNO, [-1.5E04 1.5E04], "Miguel-solved: \partial Pfb / \partial x [Pa m^{-1}]",0);
  % figure(); plot_tF(dPfby-dPfb(2,1:hnnod6),HMESH.HCOORD, HMESH.EL2HNO, [-.5E05 .5E05], "Miguel-solved: \partial Pfb / \partial y [Pa m^{-1}]",0);
   
end % function
