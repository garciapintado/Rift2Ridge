function P = fluidPressure2D(MESH,SETTINGS,VAR,PHYSICS,PBC,dt)

% THERMAL2D Two dimensional finite element thermal problem solver of MILAMIN

%   Part of MILAMIN: MATLAB-based FEM solver for large problems, Version 1.0
%   Copyright (C) 2007, M. Dabrowski, M. Krotkiewski, D.W. Schmid
%   University of Oslo, PHYSICS of Geological Processes
%   http://milamin.org
%   See License file for terms of use.

t=tic;
%==========================================================================
% MODEL INFO
%==========================================================================
GCOORD   = MESH.GCOORD;
EL2NOD   = MESH.ELEM2NODE;
nnod     = length(unique(EL2NOD));
%nnod     = size(GCOORD,2); %[nnod7]
nnodel   = size(EL2NOD,1); 
nel      = size(EL2NOD,2);
PhaseID  = MESH.PhaseID(:);

if length(VAR.Mu_f)==nnod
    if isfield(MESH,'upwind_nodes') && length(MESH.upwind_nodes)==nel
        f_prop_method = 'upwind_nodes';
        upwind_nodes  = MESH.upwind_nodes;
    else
        f_prop_method = 'average_nodes';
    end
elseif length(VAR.Mu_f)==nel
    f_prop_method = 'elcenter';
else
    error('Length of VAR.Mu_f must be either ==nnod of ==nel.');
end

%==========================================================================
% CONSTANTS
%==========================================================================
ndim         =   2;
nelblo       = 760;
nip          = SETTINGS.nip;

%==========================================================================
% BLOCKING PARAMETERS (nelblo must be < nel)
%==========================================================================
nelblo       = min(nel, nelblo);
nblo         = ceil(nel/nelblo);

%==========================================================================
% PREPARE INTEGRATION POINTS & DERIVATIVES wrt LOCAL COORDINATES
%==========================================================================
[IP_X,IP_w]  = ip_triangle(nip);
[N,dNds,~]   = shp_deriv_triangle(IP_X,nnodel);

%==========================================================================
% DECLARE VARIABLES (ALLOCATE MEMORY)
%==========================================================================
K_all        = zeros(nnodel*(nnodel+1)/2,nel);
Rhs_all      = zeros(nnodel,nel);

%==================================================================
% DECLARE VARIABLES (ALLOCATE MEMORY)
%==================================================================
K_block     = zeros(nelblo,nnodel*(nnodel+1)/2);
Rhs_block   = zeros(nelblo, nnodel);
invJx       = zeros(nelblo, ndim);
invJy       = zeros(nelblo, ndim);
il          = 1;
iu          = nelblo;

%==================================================================
% i) BLOCK LOOP - MATRIX COMPUTATION
%==================================================================
for ib = 1:nblo
    %==============================================================
    % ii) FETCH DATA OF ELEMENTS IN BLOCK
    %==============================================================
    ECOORD_x    = reshape( GCOORD(1,EL2NOD(:,il:iu)) , nnodel, nelblo);    % [nnodel, nelblo]
    ECOORD_y    = reshape( GCOORD(2,EL2NOD(:,il:iu)) , nnodel, nelblo);
    
    % FLUID PROPERTIES (stored at nodes, averaged over each element)
    switch f_prop_method
        case 'average_nodes'
            Mu_f_bl     = mean( VAR.Mu_f   (EL2NOD(:,il:iu)) )';           % [nelblo,1] fluid viscosity
            Rho_f_bl    = mean( VAR.Rho_f  (EL2NOD(:,il:iu)) )';           % [nelblo,1] fluid density
            Beta_f_bl   = mean( VAR.Beta_f (EL2NOD(:,il:iu)) )'; % fluid compressiility
            Alpha_f_bl  = mean( VAR.Alpha_f(EL2NOD(:,il:iu)) )'; % fluid thermal expansion coefficient
            Kr_bl  = mean( VAR.Kr(EL2NOD(:,il:iu)) )'; % relative permeability
           
        case 'upwind_nodes'
            Mu_f_bl     = VAR.Mu_f   (upwind_nodes(il:iu)); % fluid viscosity
            Rho_f_bl    = VAR.Rho_f  (upwind_nodes(il:iu)); % fluid density
            Beta_f_bl   = VAR.Beta_f (upwind_nodes(il:iu)); % fluid compressiility
            Alpha_f_bl  = VAR.Alpha_f(upwind_nodes(il:iu)); % fluid thermal expansion coefficient
            Kr_bl       = VAR.Kr(upwind_nodes(il:iu)); % fluid thermal expansion coefficient
            if SETTINGS.upwind_weight<1
                a          = SETTINGS.upwind_weight;
                Mu_f_bl    = a*Mu_f_bl    + (1-a)*VAR.Mu_f   (upwind_nodes(il:iu)); % fluid viscosity
                Rho_f_bl   = a*Rho_f_bl   + (1-a)*VAR.Rho_f  (upwind_nodes(il:iu)); % fluid density
                Beta_f_bl  = a*Beta_f_bl  + (1-a)*VAR.Beta_f (upwind_nodes(il:iu)); % fluid compressiility
                Alpha_f_bl = a*Alpha_f_bl + (1-a)*VAR.Alpha_f(upwind_nodes(il:iu)); % fluid thermal expansion coefficient
                Kr_bl       = a*Kr_bl + (1-a)*VAR.Kr(upwind_nodes(il:iu)); % fluid thermal expansion coefficient
            end
        case 'elcenter'
            Mu_f_bl     = VAR.Mu_f   (il:iu); % fluid viscosity
            Rho_f_bl    = VAR.Rho_f  (il:iu); % fluid density
            Beta_f_bl   = VAR.Beta_f (il:iu); % fluid compressiility
            Alpha_f_bl  = VAR.Alpha_f(il:iu); % fluid thermal expansion coefficient
            Kr_bl       = VAR.Kr(il:iu); % fluid thermal expansion coefficient
    end
    
    % MATRIX PROPERTIES (calculated using element propterty "PhaseID")
    Phi_bl      = PHYSICS.Phi (PhaseID(il:iu));
    if ~isfield(VAR,'Perm')
        Perm_bl = PHYSICS.Perm(PhaseID(il:iu));
    else
        Perm_bl = mean( VAR.Perm(EL2NOD(:,il:iu)) )';
    end
    if size(Perm_bl,2) > 1
       Perm_bl = Perm_bl(:);
    end 
    
    % PRESSURE, TEMPERATURE, SATURATION TIME DERIVATIVE
    P_bl        = VAR.P   (EL2NOD(:,il:iu)); % pressure
    DTDt_bl     = mean( VAR.DTDt(EL2NOD(:,il:iu)) )'; % temperature time derivative

    disp("JGP: pressure2d.L123:: warning, check why Lars changed this") % the day we visited Lars he changed the above to the following: Why?
    %Dsat_bl     = dt .* Phi_bl.*(Rho_f_bl - PHYSICS.rho_ah).*mean( VAR.DSDt(EL2NOD(:,il:iu)) )';
    Rphi_bl = PHYSICS.Rphi(PhaseID(il:iu));
    if size(Rphi_bl,2) > 1
        Rphi_bl = Rphi_bl(:);
    end  
    Dsat_bl     = dt .* Rho_f_bl .* Rphi_bl;
    %     
    ED          = dt.*Kr_bl.*Perm_bl.*Rho_f_bl./Mu_f_bl;
    ET          = Rho_f_bl.*Beta_f_bl.*Phi_bl;    
    Alpha_term  = dt.*DTDt_bl.*Rho_f_bl.*Alpha_f_bl.*Phi_bl;
    
    %==============================================================
    % iii) INTEGRATION LOOP
    %==============================================================
    K_block(:)      = 0;
    Rhs_block(:)    = 0;
    
    for ip=1:nip
        %==========================================================
        % iv) LOAD SHAPE FUNCTIONS DERIVATIVES FOR INTEGRATION POINT
        %==========================================================
        Ni          =  N{ip};                                              % [3,1]
        dNids       = dNds{ip};                                            % [3,2]
        
        if ip==1 % For triangular elements with non-curved edges the
                 % Jacobian is the same for all integration points
                 
            %==========================================================
            % v) CALCULATE JACOBIAN, ITS DETERMINANT AND INVERSE
            %==========================================================
            Jx          = ECOORD_x'*dNids;
            Jy          = ECOORD_y'*dNids;
            detJ        = Jx(:,1).*Jy(:,2) - Jx(:,2).*Jy(:,1);

            invdetJ     = 1./detJ;
            invJx(:,1)  = +Jy(:,2).*invdetJ;
            invJx(:,2)  = -Jy(:,1).*invdetJ;
            invJy(:,1)  = -Jx(:,2).*invdetJ;
            invJy(:,2)  = +Jx(:,1).*invdetJ;
        end
        
        %==========================================================
        % vi) DERIVATIVES wrt GLOBAL COORDINATES
        %==========================================================
        dNdx        = invJx*dNids';                                        % [nnodel,3]
        dNdy        = invJy*dNids';                                        % [nnodel,3]
        
        %==========================================================
        % vii) NUMERICAL INTEGRATION OF ELEMENT MATRICES
        %==========================================================
        weight      = IP_w(ip)*detJ;                                       % [nnodel]
        
        % WHEN USING CYLINDER COORDINATES: SCALE ELEMENT AREA BY RADIUS 
        if isfield(PHYSICS,'use_cyl_coords') && PHYSICS.use_cyl_coords==1
            radius_ip = (Ni' * ECOORD_x)';
            weight    = weight .* radius_ip;
        end
        
        N_lumped    = diag(sum(Ni*Ni',2));

        indx = 1;
        for i = 1:nnodel
            for j = i:nnodel                                               % [nnodel,1] block assignations
                K_block(:,indx)  =   K_block(:,indx) + ...
                    (dNdx(:,i).*dNdx(:,j)+ dNdy(:,i).*dNdy(:,j)).*(weight.*ED)... %...    ED \in [nnodel,nnodel] : ERROR!
                    +  N_lumped(i,j).*weight.*ET;
                indx = indx + 1;
            end
        end
        
        %RIGHT HAND SIDE                                                   [nelblo,nnodel]
        Rhs_block  = Rhs_block - dNdy .* repmat((Rho_f_bl.*weight.*ED.*PHYSICS.g),1,nnodel)...
            + (N_lumped*(P_bl.*repmat(weight'.*ET', [nnodel,1])))' ...
            + (Ni*(Alpha_term.*weight)')'...
            - (Ni*(Dsat_bl.*weight)')';
    end
    
    %==============================================================
    % ix) WRITE DATA INTO GLOBAL STORAGE
    %==============================================================
    K_all(:,il:iu)	  = K_block';
    Rhs_all(:,il:iu)  = Rhs_block';
    
    %==============================================================
    % READJUST START, END AND SIZE OF BLOCK. REALLOCATE MEMORY
    %==============================================================
    il  = il+nelblo;
    if(ib==nblo-1)
        nelblo 	  = nel-iu;
        K_block	  = zeros(nelblo, nnodel*(nnodel+1)/2);
        Rhs_block = zeros(nelblo, nnodel);
        invJx     = zeros(nelblo, ndim);
        invJy     = zeros(nelblo, ndim);
    end
    iu  = iu + nelblo;
end


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
Rhs    = accumarray(EL2NOD(:), Rhs_all(:));
clear K_i K_j K_all Rhs_all;

%==========================================================================
% BOUNDARY CONDITIONS
%==========================================================================
Free       = setdiff(1:nnod,PBC.nod); % replaces Free = 1:nnod; Free(PBC.nod)= [];
% Move prescribed temperatures to rhs
K_L        = tril(K,-1); % K without diagonal
TMP        = K(:,PBC.nod) + cs_transpose(K_L(PBC.nod,:)); clear K_L
Rhs        = Rhs - TMP*PBC.val'; clear TMP

% Add heat flux to Rhs
if isfield(PBC,'NN_INTEG')
    Rhs = Rhs + dt*PBC.NN_INTEG*PBC.valHF(:);
end



K          = K(Free,Free);
P          = VAR.P;
P(PBC.nod) = PBC.val;

%==========================================================================
% Solve equation
%==========================================================================
switch SETTINGS.solver
    case 'lchol'
        %==========================================================================
        % REORDERING
        %==========================================================================
        switch SETTINGS.reorder
            case 'metis'
                perm = metis(K);
            case 'amd'
                perm = amd(K);
            otherwise
                error('Unknown reordering')
        end

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
        P(Free(perm)) = cs_ltsolve(L,cs_lsolve(L,Rhs(Free(perm))));
        if SETTINGS.disp_profiling
            fprintf(1, ' PRESSURE (lchol)               : %7.2f sec\n',toc(t));
        end
        
    case 'agmg'
        %==========================================================================
        % Solve equation using Algebraic Multigrid algorithm by Yvan Notay
        %==========================================================================
        K       = K + tril(K,-1)';
        tol     = SETTINGS.OPT_P.rtol;
        itmax   = SETTINGS.OPT_P.itmax;
        verbose = 0;
        restart = 10; 
        [P(Free),flag,relres,nit] = agmg(K,Rhs(Free),restart,tol,itmax,verbose,P(Free));
        if flag
            error('AGMG did not copnverge.');
        end
        if SETTINGS.disp_profiling
            fprintf(1, ' PRESSURE (AGMG,tol=%.1e)    : %7.2f sec (%4i iterations)\n',...
                relres(end),toc(t),nit);
        end

    case 'chol'
        %==========================================================================
        % Solve equation using Matlab's built-in functions
        %==========================================================================
        [L,~,perm]    = chol(K,'lower','vector');
        P(Free(perm)) = L' \ (L\ Rhs(Free(perm)));
        if SETTINGS.disp_profiling
            fprintf(1, ' PRESSURE (chol)                : %7.2f sec\n',toc(t));
        end
        
    case 'pcg'
        %==========================================================================
        % Solve equation using Conjugate Gradient solver
        %==========================================================================
        K = K + tril(K,-1)';
        [P(Free),profiling] = pcg_solver(K,Rhs(Free),P(Free),SETTINGS.OPT_P);
        if SETTINGS.disp_profiling
            fprintf(1, ' PRESSURE (PCG,%s,tol=%.1e): %7.2f sec (%4i iterations)\n',...
                profiling.PC,profiling.tol,toc(t),profiling.itCG);
        end
        
    case '\'
        %==========================================================================
        % Solve equation using Matlab's backslash
        %==========================================================================
        K       = K + tril(K,-1)';
        P(Free) = K \ Rhs(Free);
        if SETTINGS.disp_profiling
            fprintf(1, ' PRESSURE (backslash)           : %7.2f sec\n',toc(t));
        end
end

end
