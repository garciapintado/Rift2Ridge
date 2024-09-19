function Z = topographical_diffusion_smooth(Dx, Z, S, sedthk, dt, SP)
    % Z = topographical_diffusion(Dx, Z, S, sedthk, dt, SP, nnodel)
    %
    % +++ purpose +++
    % Calculates a new topography by solving sediment transport PDE. 
    % Processes include: linear hill-slope diffusion, fluvial
    % sediment transport and submarine sediment transport.
    %
    % INPUT
    % Dx       :: REAL [nnod-1,1] x-increment [m] between consecutive nodes 
    % Z        :: REAL [nnod,1] topography at         : [0,cumsum(Dx)] coordinates
    % S        :: REAL [nnod,1] [m/s] external source at nodes: [0,cumsum(Dx)] coordinates [source and Neuman's boundary conditions added together]
    % sedthk   :: REAL [nnod,1] sediment thickness at : [0,cumsum(Dx)] coordinates
    % dt       :: REAL, scalar time step
    % SP       :: STRUCT
    %   .De       :: REAL scalar, fluvial transport coefficient (c x alpha^nexp)
    %   .kappa    :: REAL scalar,
    %   .kdecay   :: REAL scalar, controlling the decay of diffusivity with the water depth
    %   .kdistalb :: REAL scalar,
    %   .kdistals :: REAL scalar,
    %   .kseddump :: REAL scalar, exponent for the Type I erosion depth-dependent diffusion modulator. Normally > 1. 0 => no modulation.
    %   .ksedthk0 :: REAL scalar, nominal erodible layer thickness
    %   .ktidal   :: REAL scalar,
    %   .lriver   :: REAL [2], left/right river length 
    %   .nexp     :: REAL scalar, exponent of water flux for fluvial transport
    %   .q_bc     :: REAL [2], [m/s] left/right river sediment load from the boundaries, additional to S
    %   .sealevel :: REAL scalar, sea level
    %   .LIVDcomp :: LOGICAL, true to use dynamic qc_bc boundary as source/sink to nudge toward initial condition coordinates
    % suzo :: [m] target altitude of upper corners for surce/sink nudging 
    
    % OUTPUT 
    % Z        :: REAL [nnod] topography after the time step

    % Details:
    % This solver supersedes function_sealand() by Armitage, which used dense matrices and halted for long submerged zones.

    % Also, this function considers a limited erodible layer approximation (Type I erosion),
    % and a continuous coastal-distal diffusivity parametrization.
    % Results vary considerably with respect to function_sealand() 
    %
    % FEM is done here with linear shape functions and 2 integration points. This function cannot be
    % used at its current state with different shape functions and integration
    % points. 
    % 
    % Author: J. García-Pintado, MARUM, 2020-10-30
    %         2020-10-30 : first version 

    nnodel = 2;
    nvert  = 2;
    nelblk = 5000;

    nnod = length(Z);                    % number of nodes   

    if any([size(Dx,1)+1, size(Z,1), size(S,1), size(sedthk,1)] ~= nnod)           % matrix dimension check
        error("topographical diffusion: non-compliant matrices")
    end

    nel = ceil(nnod/(nnodel-1)) - 1;     % number of elements
    nip = nnodel;
    ELEM2NODE = repmat((1:nnodel)',1,nel) + repmat((nnodel-1)*(0:nel-1),nnodel,1); % [nnodel,nel]

    % Integration points
    [ipx, ipw] = ip_line(nip);
    [N, dNdu]  = shp_deriv_line(ipx, nnodel);                                  % N [nip] cell array, with [nnodel,1] size per cell. Shape function evaluated at integration points 
    Nm = cell2mat(N');                                                         % [nnodel, nip]
    [~, dN2du]   = shp_deriv_line(ipx, nvert);                                 % dNdu [nip] cell array, with [nnodel,2] size per cell. Derivative of shape functions at integration points
    dN2du        = dN2du{1};                                                   % [nnodel2,1]

    if any(Dx <= 0.)
        error("topographical diffusion:: non-monotonic x increments")
    end

    GCOORD = [0. cumsum(Dx)'];                                                 % [1,nnod] x-coordinates
                                                

    %==========================================================================
    % CALCULATE THE SLOPE DIRECTION AND WATER FLUX
    %==========================================================================

    Dx_sign = Dx.*(sign(diff(Z)));                                            % [nnod-1,1] positive => upward slope from left to right

    %                       . Z4 = 5
    %                      /          Z       : topography defined in the nodes
    %           .Z1 = 1   /           diff(Z) : defined at each node but the last one
    %          / \       . Z3 = 0
    %         /   \     /
    %              \   /
    %               \.Z2 = -1
    %
    % diff(Z) = [-2   1   5]
    % sign(diff(Z)) = [-1 1 1]

    % Find top of the hill and valley nodes
    indx = [find(diff(diff(Z)>0)~=0)+1; nnod]; % hill+valley node indices + final node

    X_d = zeros(nnod,1);                        % [nnod,1] distances to each node from top of corresponding hill node
    ind0 = 1;                                   % start of monotonic slope node index
    for n = 1:length(indx)
        inde            = indx(n);              % end of monotonic slope node index
        dx_sign         = Dx_sign(ind0:inde-1); % [nnodtr-1,1] where nnodtr is the number of nodes in this slope including the final node
        length_river    = abs(sum(dx_sign));    % length of the proyection of the slope in the horizontal plane
        
        % Calculate vector of horizontal distances from top of the hill to nodes in transect [monotonic slope] 
        if all(dx_sign > 0.)                            % upward left2right slope
            x_dn = length_river - cumsum([0; dx_sign]); % [nnodtr,1] accumulated positive distance from hill node
        elseif all(dx_sign < 0.)                        % downward left2right slope  
            x_dn = cumsum([0; -dx_sign]);               % [nnodtr,1] accumulated positive distance from hill node
        else
            error('Bad calculation of distance to the drainage divide')
        end
        
        % Check rivers coming out from the domain
        % ------------------
        if ind0==1 && dx_sign(1) < 0.                   % outer left river flows into domain 
            x_dn = x_dn + SP.lriver(1);
        end
        if inde==nnod && dx_sign(1) > 0.   
            x_dn = x_dn + SP.lriver(2);                 % outer right river flows into domain
        end
        
        X_d(ind0:inde) = x_dn;
        ind0 = inde;
    end

    S([1 end]) = S([1 end]) + SP.q_bc(:);                                        % add source from boundaries

    % transport diffusivity: init with fluvial (aerial) values
    Dn = SP.kappa + SP.De * (X_d.^SP.nexp);                                   % [nnod,1] transport diffusivity
    clear("X_d")

    % transport diffusivity: substitute submarine nodes
    ismarine = Z < SP.sealevel;                                               % [nnod,1]
    if (any(ismarine))
        Dnsea = SP.kistalb + SP.ktidal * exp(-SP.kdecay * (SP.sealevel - Z(ismarine)));
        Dn(ismarine) = Dnsea;                                                 %+1/365.25/24/3600;
    end

    % allocate memory                                                            [nnodel*(nnodel+1)/2==nnodel+(nnodel-1)+...+1]
    K_all   = zeros(nnodel*(nnodel+1)/2,nel);                                  % element contribution to the [nnod,nnod] extended stiffness matrix 
    Rhs_all = zeros(nnodel,nel);

    % block parameters
    nelblk = min(nel, nelblk);
    nblk   = ceil(nel/nelblk);
    il     = 1;
    iu     = nelblk;

    for ib = 1:nblk
        %disp("ib: "+ib)
        els_blk = il:iu;
        
        % get Jacobian - [subparametric transformation if nnodel=2]
        ECOORD_x = reshape( GCOORD(ELEM2NODE([1 nnodel],els_blk)), nvert, nelblk);       % REAL [nvert, nelblk]
        Jx = ECOORD_x' * dN2du;                                                      % [nelblk,1]
        detJ = Jx;                                                                   % [nelblk,1]

        %xyv0 = GCOORD(ELEM2NODE(1,els_blk));                                        % [1,nel] element start corner 
        %xyv1 = GCOORD(ELEM2NODE(nnodel,els_blk));                                   % [1,nel] element end corner
        %detJ0 = (xyv1-xyv0) / 2.;                                                    % [1,nel] |J| = element length / canonical vector length
        
        if any(detJ<0)
            error('negative |J|')
        end
        invdetJ = 1./detJ;                                                         % [nelblk,1]
        invJx   = invdetJ;

        Z_blk    = reshape(Z(ELEM2NODE(:,els_blk)), nnodel, nelblk);               % REAL [nnodel, nelblk] 
        dSdt_blk = Nm' * reshape(S(ELEM2NODE(:,els_blk)), nnodel, nelblk);                  % REAL [nip, nelblk]
        D_blk    = Nm' * reshape(Dn(ELEM2NODE(:,els_blk)), nnodel, nelblk);        % REAL [nip, nelblk]
        

        K_blk   = zeros(nnodel*(nnodel+1)/2, nelblk);                              % [sum_{i=1}^nnodel i, nelblk]; block element contribution to the [nnod,nnod] extended stiffness matrix 
        Rhs_blk = zeros(nnodel, nelblk);                                           % [nnodel, nelblk];             block element contribution to the [nnod,1] right hand side 

        % integration point loop
        for ip = 1:nip
            % canonical shape function @ ip
            Ni    = N{ip};                                                         % [nnodel,1]
            dNdui = dNdu{ip};                                                      % [nnodel,ndim], ndim=1
            
            % TODO FROM HERE
            % numerical integration of element matrices
            NxN = Ni*Ni';                                                          % [nnodel,nnodel] integration point contribution to mass matrix
            dNdx = dNdui * invJx';                                                 % [nnodel,nelblk] global coordinate derivatives
            weight_ip = reshape(ipw(ip)*detJ,1,nelblk);                            % [1,nelblkn]

            indx = 1;
            for i = 1:nnodel
                for j = i:nnodel
                    K_blk(indx,:) = K_blk(indx,:) ...                              % [1,nelblk]
                    + dt * dNdx(i,:).*dNdx(j,:) .* D_blk(ip,:) .* weight_ip ...    % stiffness matrix component
                    + NxN(i,j) .* weight_ip;                                       % mass matrix component
                    indx = indx + 1;
                end
            end
            Rhs_blk = Rhs_blk ...                                                  % [nnodel,nelblk] 
            + Ni * ((Ni' * Z_blk) .* weight_ip) ...
            + Ni * (dt * dSdt_blk(ip,:) .* weight_ip);
        end % ip
        K_all(:,els_blk)    = K_blk;
        Rhs_all(:,els_blk)  = Rhs_blk;
        
        il = iu + 1;
        if ib == nblk-1
            nelblk = nel - iu;
        end
        iu = iu + nelblk;
    end % block
 
    % create triplet 
    indx_j = repmat(1:nnodel,nnodel,1); indx_i = indx_j';                  % [nnodel,nnodel]
    indx_i = tril(indx_i); indx_i = indx_i(:); indx_i = indx_i(indx_i>0);  % [3,1] <- 3=nnodel+(nnodel-1)+1
    indx_j = tril(indx_j); indx_j = indx_j(:); indx_j = indx_j(indx_j>0);  % [3,1] 
    
    K_i = ELEM2NODE(indx_i,:); K_i = K_i(:);                               % [3,nel] -> [3*nel,1]
    K_j = ELEM2NODE(indx_j,:); K_j = K_j(:);                               % [3,nel] -> [3*nel,1]
    
    indx       = K_i < K_j;                                                % LOGICAL [3*nel,1]
    tmp        = K_j(indx);
    K_j(indx)  = K_i(indx);
    K_i(indx)  = tmp;
    
    % sparsification
    K = sparse2(K_i, K_j, K_all(:));
    Rhs    = accumarray(ELEM2NODE(:), Rhs_all(:));
    clear("K_i", "K_j", "K_all", "Rhs_all");
     
    % boundary conditions
    
    Free = 1:nnod;
    if SP.bc_type == "Dirichlet"
        K(1,:) = 0.0;
        K(1,1) = 1.0;
        Rhs(1) = SP.bc_val(1);
        K(end,:) = 0.0;
        K(end,end) = 1.0;
        Rhs(end) = SP.bc_val(2);
    end

    % reordering 
    perm = amd(K);
    % factorization
    K = cs_transpose(K);
    K = cs_symperm(K,perm);
    K = cs_transpose(K);
    L = lchol(K);

    % solve
    Z = zeros(nnod,1);
    Z(perm) = cs_ltsolve(L,cs_lsolve(L,Rhs(perm)));
    
end % function
