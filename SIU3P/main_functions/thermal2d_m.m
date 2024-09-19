function T = thermal2d_m(ELEM2NODE, GCOORD, Temp, ...                      % T,Temp: [1,nnod7]
                         RhoCp, Cond, dQdt, ...                           % [nel,nip] variables
                         Bc_ind, Bc_val, dt, reorder, nelblk)
% THERMAL2D Two dimensional finite element thermal problem solver of MILAMIN

%   Part of MILAMIN: MATLAB-based FEM solver for large problems, Version 1.0
%   Copyright (C) 2007, M. Dabrowski, M. Krotkiewski, D.W. Schmid
%   University of Oslo, Physics of Geological Processes
%   http://milamin.org
%   See License file for terms of use.

%  1. ELEM2NODE :: INT  [nnodel,nel]
%  2. GCOORD    :: REAL [2,nnod]
%  3. Temp      :: REAL [nnod,1] or [1,nnod]
%  4. RhoCp     :: REAL [nel, nip], density (Rho) times specific heat at constant pressure (Cp) defined at integration points
%  5. Hp        :: REAL [nel, nip], radiogenic heat production
%  6. Cond      :: REAL [nel, nip], thermal conductivity defined at integration points 
%  7. Hs        :: REAL [nel, nip], additional heat source
%  8. Bc_ind    :: INT  [nbct], vector with node i indices in GCOORD(:,i)
%  9. Bc_val    :: REAl [nbct], ºC values corresponding to Bc_ind
% 10. dt        :: REAL [1] (s)
% 11. reorder   :: CHAR, either 'amd' or 'metis' 
% 12. nelblk    :: INT  [1] number of element per computation block
if nargin ~= 11
    error('thermal2d_m:: incorrect number of input arguments')
end

if nargin < 11
  nelblk = 5000;
end

%==========================================================================
% MODEL INFO
%==========================================================================
nnod   = max(ELEM2NODE,[],'all');                                          % number of nodes
nnodel = size(ELEM2NODE,1);                                                % number of nodes per element [default=7 in rift2ridge]
nel    = size(ELEM2NODE,2);                                                % number of elements
nip = nnodel;                                                              % number of integration points

%==========================================================================
% CONSTANTS & error checks
%==========================================================================
ndim   = 2;
nvert  = 3;                                                                % triangle  

if ~isequal(size(RhoCp), [nel, nip])
    error('thermal2d_m:: size of RhoCp ~= [nel, nip]');
end
if ~isequal(size(Cond), [nel, nip])
    error('thermal2d_m:: size of K ~= [nel, nip]');
end

Bc_inmesh = ismember(Bc_ind,1:nnod);                                       % in case {Bc_ind, Bc_val} include central nodes
if ~all(Bc_inmesh) 
    % error("thermal2d_m :: ~all(Bc_inmesh)")                              
    Bc_ind = Bc_ind(Bc_inmesh); 
    Bc_val = Bc_val(Bc_inmesh);
end

%==========================================================================
% PREPARE INTEGRATION POINTS & DERIVATIVES wrt LOCAL COORDINATES
%==========================================================================
[x_ip, IP_w] = ip_triangle(nip);                                           % x_ip [nip,2], IP_w [1,nip] with nip=6; 
[N, dNdu]    = shp_deriv_triangle(x_ip, nnodel);                           % N [nip] cell array, with [nnodel,1] size per cell. Shape function evaluated at integration points 
[~, dN3ds]   = sf_dsf_tri367(x_ip',3,'cell');                                                                      % dNdu [nip] cell array, with [nnodel,2] size per cell. Derivative of shape functions at integration points
dN3ds        = dN3ds{1}';                                                  % [nnodel3,2]
                                                                           
%==========================================================================
% DECLARE VARIABLES (ALLOCATE MEMORY)
%==========================================================================
K_all   = zeros(nnodel*(nnodel+1)/2,nel);                                      % nnodel*(nnodel+1)/2==nnodel+(nnodel-1)+...+1
                                                                               % element contribution to the [nnod,nnod] extended stiffness matrix 
Rhs_all = zeros(nnodel,nel);                                                   % element contribution to the [nnod,1] right hand side 

%==================================================================
% DEFINE BLOCK PARAMETERS
%==================================================================
nelblk = min(nel, nelblk);
nblk   = ceil(nel/nelblk);
il     = 1;
iu     = nelblk;

%==================================================================
% i) BLOCK ELEMENT LOOP - MATRIX COMPUTATION
%==================================================================
%         fprintf(1, 'MATRIX ASSEMBLY:    '); tic;
for ib = 1:nblk
    els_blk    = il:iu;                                                        % elements in blok                                                  
    %==============================================================
    % i.1 GET JACOBIAN, J, J^{-1}, and their determinants
    %==============================================================
    % In these model we use linear coordinate transformation, which is 
    % subparametric with respect to the trial function (i.e. quadratic 
    % shape functions for Temperature interpolation and Galerkin weights).
    % That is, physical space triangles have straight edges, and this 
    % => constant |J| throughout each element => no need to calculate it at
    % each integration point

    ECOORD_x = reshape( GCOORD(1,ELEM2NODE(1:nvert,els_blk)), nvert, nelblk);     % REAL [nvert, nelblk]
    ECOORD_y = reshape( GCOORD(2,ELEM2NODE(1:nvert,els_blk)), nvert, nelblk);     % REAL [nvert, nelblk]
    Jx       = ECOORD_x' * dN3ds;                                                 % [nelblk, 2] each row : [\frac{\partial x}{\partial \xi} \frac{\partial x}{\partial \eta}]
    Jy       = ECOORD_y' * dN3ds;                                                 % [nelblk, 2] each row : [\frac{\partial y}{\partial \xi} \frac{\partial y}{\partial \eta}]
    detJ = Jx(:,1).*Jy(:,2) - Jy(:,1).*Jx(:,2);                                   % [nelblk, 1] \frac{\partial x}{\partial \xi} * \frac{\partial y}{\partial \eta} - \frac{\partial y}{\partial ji} * \frac{\partial x}{\partial eta}
    if any(detJ<0)
        error('negative |J|')
    end  
    invdetJ    = 1./detJ;                                                         % [nelblk, 1]
    invJx      = zeros(nelblk, ndim);                                             % invJx:: [nelblk, 2] :: [\xi/x ; \eta/x]
    invJy      = zeros(nelblk, ndim);                                             % invJy:: [nelblk, 2] :: [\xi/y ; \eta/y]
    invJx(:,1) = +Jy(:,2).*invdetJ;                                               % $\frac{\partial \xi}{x}$
    invJx(:,2) = -Jy(:,1).*invdetJ;                                               % $\frac{\partial \eta}{x}$
    invJy(:,1) = -Jx(:,2).*invdetJ;                                               % $\frac{\partial \xi}{y}$
    invJy(:,2) = +Jx(:,1).*invdetJ;                                               % $\frac{\partial \eta}{y}$
   
    T_blk   = reshape(Temp(ELEM2NODE(:,els_blk)), nnodel, nelblk);         % REAL [nnodel, nelblk] 
    K_blk   = zeros(nnodel*(nnodel+1)/2, nelblk);                          % block element contribution to the [nnod,nnod] extended stiffness matrix 
    Rhs_blk = zeros(nnodel, nelblk);                                       % block element contribution to the [nnod,1] right hand side 
    %==============================================================
    % i.2) INTEGRATION POINT LOOP
    %==============================================================
    for ip=1:nip    
        %==========================================================
        % i.2.1) LOAD canonical SHAPE FUNCTIONS & DERIVATIVES evaluated @ this integration point 
        %==========================================================
        Ni        = N{ip};                                                 % [nnodel, 1]
        dNdui     = dNdu{ip};                                              % [nnodel, 2] $\bN_{\xi\eta}$ canonical partial derivatives [$1/\partial{\xi}$, /$1/partial{\eta}$] for all shape functions evaluated at this integration point
                
        %==========================================================
        % i.2.2) NUMERICAL INTEGRATION OF ELEMENT MATRICES
        %==========================================================
        % ip terms for mass matrix - constant for all elements 
        % NxN    = diag(sum(Ni*Ni',2)); it can work faster (only possible for nnodel=3)
        NxN  = Ni*Ni';                                                     % [nnodel,nnodel] integration point contribution to mass matrix
        % ip terms for stiffness matrix - convert derivatives of the shape functions from canonical to physical coordinates
        dNdx = dNdui * invJx';                                             % [nnodel,nelblk]  <- [nelblk, 2][2,nnodel] $\frac{\partial N}{\partial x}$
        dNdy = dNdui * invJy';                                             %   "           "                "          $\frac{\partial N}{\partial y}$

        weight_ip = reshape(IP_w(ip)*detJ,1,nelblk);                       % REAL [1,nelblk]  
        
        Cond_blk  = Cond(els_blk,ip)';                                     % REAL [1,nelblk] thermal conductivity
        RhoCp_blk = RhoCp(els_blk,ip)';                                    % REAL [1,nelblk] 
        dQdt_blk = dQdt(els_blk,ip)';                                       % heat sources
        %dQdt_blk  = Hp(els_blk,ip)' + Hs(els_blk,ip)';                     % REAL [1,nelblk] radiogenic + shear heat production [J.s-1]
        
        indx = 1;                                                          % indx: column in extended stiffness matrix
        for i = 1:nnodel
            for j = i:nnodel
                K_blk(indx,:) = K_blk(indx,:) ...                                               % extended stiffness
                    + NxN(i,j) .* RhoCp_blk .* weight_ip ...                                    % mass matrix term           rho_i * cp_i * NxX  * |J_i| * w_i        
                    + dt*(dNdx(i,:).*dNdx(j,:) + dNdy(i,:).*dNdy(j,:)) .* Cond_blk .* weight_ip; % stiffness matrix term      dt * cond_i * dNxdN * |J_i| * w_i
                    %     Diffusion::computeQpResidual() { return _grad_test[i][_qp] * _grad_u[_qp] }
                    indx = indx + 1;
            end
        end
        
        %RIGHT HAND SIDE
        Rhs_blk  = Rhs_blk ...                                             % [nnodel, nelblk]
            + Ni * ((Ni'*T_blk) .* RhoCp_blk .* weight_ip) ...             % [nnodel, nelblk] previous temperature
            + Ni * (dt * dQdt_blk .* weight_ip);                           % [nnodel, nelblk] <- [nnodel,1].[1,nelblk] heat source
    
    end
    %==============================================================
    % ix) WRITE DATA INTO GLOBAL STORAGE
    %==============================================================
    K_all(:,els_blk)    = K_blk;
    Rhs_all(:,els_blk)  = Rhs_blk;
    %===========================================
    % READJUST START, END AND SIZE OF NEXT BLOCK
    %===========================================
    il = iu + 1;
    if ib == nblk-1
        nelblk = nel - iu;
    end
    iu = iu + nelblk;
end

%==========================================================================
% ix) CREATE TRIPLET FORMAT INDICES
%==========================================================================
% fprintf(1, 'TRIPLET INDICES:    '); tic
indx_j = repmat(1:nnodel,nnodel,1); indx_i = indx_j';                        % [nnodel,nnodel]
indx_i = tril(indx_i); indx_i = indx_i(:); indx_i = indx_i(indx_i>0);        % [21,1] <- 21=nnodel+(nnodel-1)+...+1
indx_j = tril(indx_j); indx_j = indx_j(:); indx_j = indx_j(indx_j>0);        % [21,1]

K_i = ELEM2NODE(indx_i,:); K_i = K_i(:);                                   % [21,nel] -> [21*nel,1]
K_j = ELEM2NODE(indx_j,:); K_j = K_j(:);                                   % [21,nel] -> [21*nel,1]

indx       = K_i < K_j;                                                    % LOGICAL [21*nel,1]
tmp        = K_j(indx);
K_j(indx)  = K_i(indx);
K_i(indx)  = tmp;
% fprintf(1, [num2str(toc),'\n']);

%==========================================================================
% x) CONVERT TRIPLET DATA TO SPARSE MATRIX
%==========================================================================
% fprintf(1, 'SPARSIFICATION:     '); tic
K_all  = K_all(:);
K      = sparse2(K_i, K_j, K_all);
Rhs    = accumarray(ELEM2NODE(:), Rhs_all(:));
clear K_i K_j K_all Rhs_all;
% fprintf(1, [num2str(toc),'\n']);

%==========================================================================
% BOUNDARY CONDITIONS
%==========================================================================
% fprintf(1, 'BDRY CONDITIONS:    '); tic;
Free        = 1:nnod;
Free(Bc_ind)= [];

TMP         = K(:,Bc_ind) + cs_transpose(K(Bc_ind,:));                     % [nnod, nbdy]
Rhs         = Rhs - TMP * Bc_val(:);                                       % [nnod,1]
K           = K(Free,Free);                                                % standard indexing
% fprintf(1, [num2str(toc),'\n']);

%==========================================================================
% REORDERING
%==========================================================================
% fprintf(1, 'REORDERING:         '); tic;
switch reorder
    case 'metis'
        perm = metis(K);
    case 'amd'
        perm = amd(K);
    otherwise
        error('Unknown reordering')
end
% fprintf(1, [num2str(toc),'\n']);

%==========================================================================
% FACTORIZATION - ideally L = lchol(K, perm)
%==========================================================================
% fprintf(1, 'FACTORIZATION:      '); tic;
K = cs_transpose(K);
K = cs_symperm(K,perm);
K = cs_transpose(K);
L = lchol(K);
% fprintf(1, [num2str(toc,'%8.6f'),'\n']);

%==========================================================================
% SOLVE
%==========================================================================
% fprintf(1, 'SOLVE:              '); tic;
T             = zeros(1,nnod);
T(Bc_ind)     = Bc_val;
T(Free(perm)) = cs_ltsolve(L,cs_lsolve(L,Rhs(Free(perm))));
% fprintf(1, [num2str(toc,'%8.6f'),'\n']);
