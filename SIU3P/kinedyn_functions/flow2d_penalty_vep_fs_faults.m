 function VAR = flow2d_penalty_vep_fs_faults(VAR,FAULTS,MESH,UBC,SETTINGS,PHYSICS,NUMSCALE,dt,FigNo0,SERP)
% Usage: VAR = flow2d_penalty_vep_fs(VAR, FAULTS, MESH, UBC, SETTINGS, PHYSICS, NUMSCALE, dt, FigNo0, SERP)
%
% Purpose: Visco-elasto-plastic flow solver; fixed or free top surface
%
% This is a version (aimed to be backward-compatible) of flow2d_penalty_vep_fs_faults() in kinedyn
% Basically the aim is to improve comments and check consistency with Miguel's mechanical2d_m.m 
% 
% Input:
% 1. VAR      : [structure] :     variable fields: (Ux           IO [MESH.nU, 1]
%                                                   Uz           IO [MESH.nU, 1]
%                                                   P            IO [MESH.nP, 1]
%              stresses:                            Tau_xx       I  [nel, nip_stress], with nip_stress=3 by default 
%                                                   Tau_xz       I  [nel, nip_stress]
%                                                   Tau_zz       I  [nel, nip_stress] 
%                                                   Tau_II_totl  I  [nel, nip_stress] max shear stress
%              strain rates:                        Er_xx        I  [nel, nip_stress]
%                                                   Er_xz        I  [nel, nip_stress]
%                                                   Er_zz        I  [nel, nip_stress]
%                                                   Er_II_totl   I  [nel, nip_stress]
%              elasticity args:                     Vort_xz      IO [nel, nip_stress] vorticity; used for elasticity; used and then recalculated here for NEXT time step
%              plasticity args:                     Tau_yld      I  [nel, nip_stress] OPTIONAL, yield stress; this group only used if 'is_plastic'
%              if plastic & strain softening args:  Eacc_visc    I  [nel, nip_stress] OPTIONAL,                                                                                   
%                                                   Eacc_plastic I  [nel, nip_stress] OPTIONAL,
%                                                   Eacc_II      I  [nel, nip_stress] OPTIONAL,
%                                                   Cohesion     I  [nel, nip_stress] OPTIONAL, cohesion, if not given, initialized from PHYSICS.Cohesion0(MESH.PhaseID)
%                                                   PhiFric      I  [nel, nip_stress] OPTIONAL, friciton angle, if not given, initialized from PHYSICS.PhiFric0(MESH.PhaseID)
%                                                   F_xx         I  [nel, nip_stress] OPTIONAL,
%                                                   F_zz         I  [nel, nip_stress] OPTIONAL,
%                                                   F_zx         I  [nel, nip_stress] OPTIONAL,  
%                                                   Visc_v          ?
% 1. FAULTS   : array of structures, where each element array contains:
%                     FAULTS(i):                    activated     CHARACTER ['yes' or 'no']                   FAULTS is only used by SUBROUTINE calc_accumulated_strains       
%                                                   deactivated   CHARACTER ['yes' or 'no']                   If no faults are considered, a dummy FAULTS input is a length-one array as:
%                                                   PointID       INTEGER [], identifier of fault nodes        FAULTS(1).activated = 'no'
%                                                   dipdir        CHARACTER ['left' or 'right']                FAULTS(1).deactivated = 'no'
%  
% 3. MESH     : [structure] : FE mesh parameters
%                                                   nel           INTEGER, number of elements 
%                                                   nU            INTEGER, number of velocity nodes
%                                                   nP            INTEGER, number of pressure nodes
%                                                   PointID       INTEGER, [:]  
%                                                   PointID_bot   INTEGER, [:] point identifier, within 1:nnod of nodes at the bottom of the domain
%                                                   PointID_top   INTEGER, [:] point identifier, within 1:nnod of nodes at the top of the domain 
%                                                   GCOORD        INTEGER, [2,nnod] 
%                                                   EL2NOD        INTEGER, [nnodel, nel] 
%                                                   EL2NODP       INTEGER, [3, nel]  
%  
% 4. UBC      : [structure] : velocity boundary conditions
%                                                   ifix          INTEGER, [:] indices (within 1:nnod?) of nodes with fixed velocity  
%                                                   vfix          REAL,    [:] corresponding velocity values (vertical velocity only?)
%                                                   bottom_bc     CHARACTER, OPTIONAL ['solve_for': solve for bottom inflow]
%                                                   iPdof_ignore_divU INTEGER [:], or LOGICAL [:], OPTIONAL, used to subset divU values to calculate norm_divU in pressure iterations      
% 5. SETTINGS : [structure] : model parameters
%                                                   visc_in_el    CHARACTER ['min', 'mean', 'linear' or 'any']
%                                                   is_elastic    CHARACTER ['yes' or 'no']
%                                                   is_plastic    CHARACTER ['yes' or 'no']
%                                                   MODEL         CHARACTER, OPTIONAL ['kinedyn' or 'rift2ridge2D']. Defaults to 'kinedyn'
% 6. PHYSICS  : [structure] : physical properties
%                                                   g             REAL, gravitational constant
%                                                   minVisc       REAL, lower cut-off for viscosity
%                                                   maxVisc       REAL, upper cut-off for viscosity
%                                                   ShearG        REAL, [nphases]
%                                                   SS            ?
% 7. NUMSCALE : [structure] : numerical scaling parameters
% 8. dt       : [double]    : time step
% 8. FigNo0
% 10.SERP
%   
% Output:
%   VAR      : [structure] : major variable fields
%   this function additional fields:                Visc
%                                                   divU   divergence of velocity field%
%                                                   Tau_xx_old : [nel, nip_stress], generated by local function rotate_stresses(). In rift2ridge2D these are calculated after the mechanical within main() 
%                                                   Tau_zz_old : [nel, nip_stress],    "              "            " 
%                                                   Tau_xz_old : [nel, nip_stress],    "              "            " 
%
%   
% Part of M2TRI_LGR - 2D FINITE ELEMENT CODE
%
% Original (pure-viscous) version developed by J.Hasenclever & J.P. Morgan,
% 2007-2010. Matrix assembly-style based on MILAMIN.
% Email contact: jhasenclever@geomar.de
% For numerical methods see online Ph.D. thesis of J.Hasenclever
% (http://www.sub.uni-hamburg.de/opus/volltexte/2010/4873/)
%
% Details on this flow solver:
% - Consistent penalty method, Powell-Hestenes (Uzawa) iterations
% - elements: 2D Crouzeix-Raviert triangles
%             (7 velocity nodes, 3 discontinuous pressure nodes)
% - Cholesky direct solver
% - element assembly using standard method OR much faster block-assembly
%   (MILAMIN style)
% - rheology either defined at nodes or per element
% - current version is visco-elastic and has either free or fixed surface
% - plasticity added but not yet benchmarked
%
% JH Jan 2011
% JH/JPM Feb2011
% JH Jan 2012
% JH Feb 2014: added free surface (algorithm by Miguel & Jason)
% JH Feb 2014: added visco-elastic formulation
%              (based o work by Lars, Tine, Jason, Miguel, Marta)
% JH Apr 2014: cleaned up; benchmarked
% JH Aug 2014: restructured to improve iterative updates of stiffness
%              matrix (e.g. for powerlaw, strain-softening, plasticity)
% JH May 2017: added powerlaw & plasticity
%              (based on work by Miguel, Albert, Elena, Lars)
% JH Jul 2017: added strain softening
%              (based on work by Miguel, Albert, Elena)
%

t0=tic;
% ====

  if isfield(SETTINGS,'MODEL')
    MODEL = SETTINGS.MODEL
  else
    MODEL = "kinedyn"
  end
  if ~ismember(MODEL,["kinedyn","rift2ridge2D"])
    error('flow2d_penalty :: MODEL should be either "kinedyn" or "rift2ridge2D"')
  end
  
% =========================================================================
% SETUP FOR VISCO-ELASTO-PlASTIC FLOW SOLVER
% =========================================================================
if nargin>7 && FigNo0>0
    FigNo_PenaltyConv = 1000+FigNo0; % if ==0 --> no plot
      % number of figure showing convergence of Powell-Hestenes iterationms
    FigNo_StrainConv  = FigNo0; % if ==0 --> no plot
      % number of figure showing convergence of strain rate iterationms
else
    FigNo_PenaltyConv = 0; % if ==0 --> no plot
      % number of figure showing convergence of Powell-Hestenes iterationms
    FigNo_StrainConv  = 0; % if ==0 --> no plot
      % number of figure showing convergence of strain rate iterationms
end
itminStrainRate = 5;
  % minimum number of strain rate iterarions
itmaxStrainRate = 50;
  % maximum number of strain rate iterarions
itmaxPowHest    = 5;
  % maximum number of Powell-Hestenes (pressure) iterations
  % no tolerance is defined for leaving the iterations - instead iterations
  % are aborted when convergence stalls at machine precision
dU_tol     = 1e-6; % 0 --> disabled, using rX_tol instead
  % max relative change in velocity (to exist strain rate iterations)
rX_tol     = 0; % 0 --> disabled, using dU_tol instead
  % max norm of non-linear residual (to exist strain rate iterations)
fpenalty   = 1e4;
  % penalty factor for Uzawa solution (typically 1e4-1e8)
  % This factor must not include max(viscosity), as it will be multiplied
  % later by max(viscosity)!

% Options for matrix assembly
% (based on MILAMIN method; see Dabrowski et al. 2008 (doi:10.1029/2007GC001719)
nelblk        = 1200;
  % number of elements assembled at once (block size)
nelblk_stress = 1200;
  % number of elements for which stresses are calculated at once (block size)

nip_flow   = 6; % 6 or 7
  % number of integration points

%==========================================================================
% Check that the rheology is defined correctly & set nip_stress
%==========================================================================
if ~isfield(SETTINGS,'is_elastic') || ...
  (~strcmp(SETTINGS.is_elastic,'yes') && ...
   ~strcmp(SETTINGS.is_elastic,'no'))
    error(' SETTINGS.is_elastic undefined. Set to "yes" or "no".');
end
if ~isfield(SETTINGS,'is_plastic') || ...
  (~strcmp(SETTINGS.is_plastic,'yes') && ...
   ~strcmp(SETTINGS.is_plastic,'no'))
    error(' SETTINGS.is_plastic undefined. Set to "yes" or "no".');
end
if ~isfield(SETTINGS,'visc_in_el')
    if strcmp(SETTINGS.is_plastic,'yes')
        error(' SETTINGS.visc_in_el undefined. Allowed are "min", "mean", "linear", or "any".');
    else
        SETTINGS.visc_in_el = 'any';
    end
end
if strcmp(SETTINGS.is_elastic,'no') && nargin < 8
    dt = 0;
end

if strcmp(SETTINGS.visc_in_el, 'any')
    nip_stress = nip_flow;
else
    nip_stress = 3; % 3, 6, 7
      % number of integration points for stress calculation
      % Use 3 integration points for stress integration, because these
      % points have been found to be the most accurate ones.
      % S Rajendran, KM Liew: "Optimal stress sampling points of plane
      % triangular elements for patch recovery of nodal stresses",
      % International Journal for Numerical Methods in Engineering, 2003
end

% =========================================================================
% Prepare velocity and pressure solution vectors
% =========================================================================
if ~isfield(VAR,'Ux')
    VAR.Ux = zeros(MESH.nU,1);
end
if ~isfield(VAR,'Uz')
    VAR.Uz = zeros(MESH.nU,1);
end
if ~isfield(VAR,'P')
    VAR.P = zeros(MESH.nP,1);
end
U       = [VAR.Ux(:)';VAR.Uz(:)'];
U       = U(:);                                   % [nU*2,1]
P       = VAR.P;                                  % [nP,1]
nU      = length(U);
nP      = length(P);
ifix    = UBC.ifix;                               % list of constrained (fixed) velocity degrees of freedom
vfix    = UBC.vfix;                               % list of values at fixed velocity degrees of freedom
if ~issorted(ifix)
    [ifix,perm] = sort(ifix);
    vfix        = vfix(perm);
end

% Initialize variables related to calculation of strain rate, stress,
% deformation, etc. This is only done in the first time step, and only if
% the fields are not yet allocated (or after a remesh).
VAR = allocate_fields(VAR,MESH,SETTINGS,PHYSICS,nip_stress);     % *LOCAL FUNCTION*

% if strcmp(SETTINGS.is_elastic,'yes')
    %======================================================================
    % ROTATE OLD ELASTIC STRESSES (JAUMANN ROTATION)
    %======================================================================
    VAR = rotate_stresses(VAR, dt);                              % *LOCAL FUNCTION*
% end

if strcmp(SETTINGS.is_plastic,'yes') && strcmp(SETTINGS.use_strainsoft,'yes')
    % =========================================================================
    % ACCUMULATED STRAIN FOR STRAIN SOFTENING
    % =========================================================================
    VAR = calc_accumulated_strains(VAR,FAULTS,MESH,U,dt,nelblk); % *LOCAL FUNCTION*
    
%     fault_nodes=find(ismember(MESH.PointID,MESH.PointID_fault));
%     [~,fault_elements]=find(ismember(MESH.EL2NOD(1:3,:),fault_nodes));
%     VAR.Eacc_II(fault_elements,1:3)=1;
end

% =========================================================================
% Begin of flow solver iterations (only when including plasticiy,
% Non-Newtonian rheology (power law), or strain-softening)
% =========================================================================
converged      = 0;
itPowHest_totl = 0;
for itStrainRate = 1:itmaxStrainRate
    %======================================================================
    % BOUNDARY CONDITIONS (DIFFERENT WAYS FOR BOTTOM INFLOW)
    %======================================================================
    if isfield(UBC,'bottom_bc') && strcmp(UBC.bottom_bc,'solve_for')
        if itStrainRate == 1
            % IN FIRST ITERATION ONLY:
            % (1) disable gravity
            g           = PHYSICS.g;
            PHYSICS.g   = 0;
            % (2) make free inflow bottom boundary and solve for the
            %     vertical inflow profile
            idof_Uz_bot = 2*find(ismember(MESH.PointID,MESH.PointID_bot));
            indx_Uz_bot = find(ismember(ifix,idof_Uz_bot));
            ifix(indx_Uz_bot) = [];
            vfix(indx_Uz_bot) = [];
%             tmp.ifix = ifix; tmp.vfix = vfix; plot_velocity_bc(661,MESH,tmp);
        elseif itStrainRate == 2
            % IN SECOND ITERATION ONLY:
            % (1) enable gravity again (was disabled above)
            PHYSICS.g = g;
            % (2) set vertical bottom velocity to the values calculated in
            %     the first iteration
            ifix      = UBC.ifix; % list of constrained (fixed) velocity degrees of freedom
            vfix      = UBC.vfix; % list of values at fixed velocity degrees of freedom
            vfix(indx_Uz_bot) = U(idof_Uz_bot);
%             tmp.ifix = ifix; tmp.vfix = vfix; plot_velocity_bc(662,MESH,tmp);
%             Ux_plot  = U(1:2:end);
%             Uz_plot  = U(2:2:end);
%             plot_2d_fedata(172, MESH.GCOORD, MESH.EL2NOD, Ux_plot, Ux_plot, Uz_plot, 'none', true, MODEL);title('Ux');
%             plot_2d_fedata(173, MESH.GCOORD, MESH.EL2NOD, Uz_plot, Ux_plot, Uz_plot, 'none', true, MODEL);title('Uz');
        end
    end
    
    U(ifix)     = vfix; % prescribed velocity values
    ifree       = 1:nU;
    ifree(ifix) = [];   % list of free velocity degrees of freedom
    
    % =====================================================================
    % ASSEMBLY OF GLOBAL MATRICES
    % =====================================================================
    if itStrainRate == 1
        [KK,ViscEl_ip,VAR,GG,invMM,Rhs0,DensEl] = matrix_assembly(VAR,MESH,...
            SETTINGS, PHYSICS, NUMSCALE, dt, fpenalty, nelblk, nip_flow, itStrainRate); % *LOCAL FUNCTION*
        % =====================================================================
        % Winkler
        % =====================================================================
        nodes_bot   = find(ismember(MESH.PointID,MESH.PointID_bot));
        surf_elem   = find(sum(ismember(MESH.EL2NOD(1:6,:),nodes_bot))==3);
        local_nodes = ismember(MESH.EL2NOD(:,surf_elem),nodes_bot);
        F_ip_ast    = SETTINGS.P_lithstat_w*ones(3,length(surf_elem));
        F_ext0 = int_f_surf(surf_elem,local_nodes,F_ip_ast,2,MESH.GCOORD,MESH.EL2NOD); % external function
        F_ext0 = F_ext0(:);
        Rhs0=Rhs0+F_ext0;
    
    else
        [KK,ViscEl_ip,VAR] = matrix_assembly(VAR,MESH,...
            SETTINGS,PHYSICS,NUMSCALE,dt,fpenalty,nelblk,nip_flow,itStrainRate); % *LOCAL FUNCTION*
    end     
    
    % =====================================================================
    % FREE SURFACE TERMS
    % =====================================================================
    switch SETTINGS.top_surface
        case 'free'
            % Free surface
            if itStrainRate==1
                alpha = SETTINGS.fs_alpha; % weighting begin/end of time step
                beta  = SETTINGS.fs_beta;  % 1 --> incl. slope effects
                [KKfs_zz,~,RHSfs_z] = free_surface_bc(MESH, PHYSICS, NUMSCALE, ...
						      DensEl, alpha, beta, dt);
                
%                 G=[0, PHYSICS.g];
%                 [RHSfs_z,KKfs_zz,~] = fssa(MESH,NUMSCALE,VAR.Dens,G,dt,3, ...
%                 alpha,beta,MESH.PointID_top,[],PHYSICS.Dens_water);
            
            end
            KK            = KK - tril(KKfs_zz);
            itmaxFreeSurf = 10;
            tol_FS        = 1e-6; % tolerance for leaving free surface iterations
            
        case 'fix'
            % Fixed surface
            itmaxFreeSurf = 1;
            beta          = 0;
            RHSfs_z       = 0;
            
        otherwise
            error(' SETTINGS.top_surface must be either "free" or "fix"');
    end  % switch
    
    % =====================================================================
    % Boundary conditions
    % =====================================================================
    KKbc    = KK(:,ifix) + KK(ifix,:)';
    bc_rhs  = KKbc*vfix(:); clear KKbc
    %beta is always zero so dont change it    
    Rhs     = Rhs0 - bc_rhs + GG*P;% + RHSfs_z; % build right-hand-side vector %HERE
    
    if itStrainRate <= 2
        try
            [LL,perm] = cholesky_factorize(KK(ifree,ifree)); % Factorize KK
        catch
            disp('ERROR: Cholesky factorization failed!! Check mesh in figure(666)!');
            plot_2d_fedata(666, MESH.GCOORD, MESH.EL2NOD, MESH.PhaseID, [], [], 'k', true, MODEL);
            title('Mesh for which Cholesky factorizatiuon failed');
            [LL,perm] = cholesky_factorize(KK(ifree,ifree)); % Factorize KK
        end
    else
        % can re-use perm; saves a bit of time
        LL = cholesky_factorize(KK(ifree,ifree),perm); % Factorize KK
    end
    
    % =========================================================================
    % CALCULATE RESIDUALS OF MOMENTUM AND CONTINUITY EQUATIONS
    % =========================================================================
    if itStrainRate > 1
        % Calculate residuals
        KK2         = KK(ifree,ifree) + tril(KK(ifree,ifree),-1)'; % KK2 at current iteration
        Rhs_Gp      = Rhs(ifree);
        r1          = Rhs_Gp - KK2*U(ifree);
        r2          = GG'*U;
        norm_rX     = norm([r1(:); r2(:)]);

        % Change in velocity and pressure during previous iteration
        dU_conv     = norm(U-U_last)/norm(U);
        dP_conv     = norm(P-P_last)/norm(P);
        
        % Check for convergence of strain rate iterations
        converged   = 0;
        if ~isempty(dU_tol) && dU_tol > 0
            if dU_conv<dU_tol
                converged(1) = 1;
            else
                converged(1) = 0;
            end
        else
            converged(1) = 1;
        end
        if ~isempty(rX_tol) && rX_tol>0 
            if norm_rX<rX_tol
                converged(2) = 1;
            else
                converged(2) = 0;
            end
        else
            converged(2) = 1;
        end
        converged = min(converged);
        
        if FigNo_StrainConv
            sfigure(FigNo_StrainConv);
            if itStrainRate == 2
                clf
            end
            plot(itStrainRate,log10(norm_rX),'r.'); hold on
            plot(itStrainRate,log10(dU_conv),'k.');
            plot(itStrainRate,log10(dP_conv),'b.');
            if itStrainRate == 2
                if ~isempty(dU_tol) && dU_tol>0
                    line([0 itmaxStrainRate],log10([dU_tol dU_tol]),'Color','k');
                end
                if ~isempty(rX_tol) && rX_tol>0
                    line([0 itmaxStrainRate],log10([rX_tol rX_tol]),'Color','r');
                end
                xlabel('strain rate iteration');
                ylabel('log_{10}');
                legend('norm(rX)','dU','dP');
                set(gca,'YLim',[-8 5],'XLim',[0 itmaxStrainRate]); hold all
            end
            drawnow
        end
    end % if strainRate > 1
    
    % =====================================================================
    % Begin of free surface iterations (only when including slope-effect)
    % =====================================================================
    U_last   = U; % needed to check convergence for plasticity
    P_last   = P;
    Rhs_nofs = Rhs;
    for itFreeSurf = 1:itmaxFreeSurf
        if beta==0
            if itFreeSurf==2
                break % no slope effect considered ---> no need to iterate
            end
        else
            error('Free surface iterations must be benchmarked');
            if itFreeSurf>1 && norm(U0-U)/norm(U)<tol_FS
                break
            end
            U0  = U;
            Rhs = Rhs_nofs + KKfs_xz * U;
            P   = P_last;
        end
        
        % =================================================================
        % Begin of Powell-Hestenes pressure iterations
        % =================================================================
        norm_divU = zeros(1,itmaxPowHest);
        for itPowHest = 1:itmaxPowHest
            U(ifree) = cholesky_solve(Rhs(ifree),LL,perm);
            
            % M2TRI:
            divU     = -GG'*U;
            dP       = (fpenalty*PHYSICS.maxVisc) * (invMM * divU);
            Rhs      = Rhs + GG*dP;
            P        = P + dP;
            
            if isfield(UBC,'iPdof_ignore_divU')
                tmp = divU;  tmp(UBC.iPdof_ignore_divU) = [];
                norm_divU(itPowHest) = norm(tmp); % check if solution converged
            else
                norm_divU(itPowHest) = norm(divU); % check if solution converged
            end
            
% %             % for comparison: equivalent MILAMIN lines
% %             if itPowHest==1
% %                 GG = -GG;
% %                 PF = (fpenalty*PHYSICS.maxVisc);
% %             end
% %             divU     = invMM * (GG'*U);
% %             Rhs      = Rhs - PF*(GG*divU);
% %             P        = P + PF*divU;
% %             div_max  = max(abs(divU(:)));
% %             disp([' PH_ITER: ', num2str(itPowHest), ' ', num2str(div_max)]);
% %             norm_divU(itPowHest) = norm(invMM * divU); % check if solution converged

            if itPowHest>2 && (norm_divU(itPowHest)/norm_divU(itPowHest-1)>0.8)
                break % leave if covergence stalls at machine precision
            end
        end % END OF POWELL-HESTENES ITERATIONS
        itPowHest_totl = itPowHest_totl + itPowHest;
        
        if FigNo_PenaltyConv % Plot convergence
            sfigure(FigNo_PenaltyConv); if itStrainRate==1 && itFreeSurf==1; clf; end
            plot(log10(norm_divU(1:itPowHest)),'r.-'); hold on
            % axis([0 itmaxPowHest 1e-14 1e2]); hold all
            set(gca,'FontSize',8); grid on
            xlabel('Pressure iteration'); ylabel('Norm of Np-integrated divergence vector');
            drawnow;
        end
    end % for itFreeSurf (FREE SURFACE ITERATIONS)
    
    % Check solution
    if any(isnan(U)) || any(isnan(P))
        error('Solution containes NaNs !!!!');
    end
    
    % Store velocity/pressure solution in structure VAR
    VAR.Ux   = U(1:2:end); % horizontal velocity solution
    VAR.Uz   = U(2:2:end); % vertical velocity solution
    VAR.P    = P;          % pressure solution
    VAR.divU = divU;       % divergence of velocity field

    % =========================================================================
    % POST PROCESSING (STRESSES, etc)
    % =========================================================================
    VAR = calc_stresses(VAR,ViscEl_ip,MESH,U,SETTINGS,PHYSICS,NUMSCALE,dt,nelblk_stress,SERP); % *LOCAL FUNCTION*

%     if strcmp(SETTINGS.is_plastic,'yes')
%         VAR.Tau_exc = max(0,VAR.Tau_II_totl-VAR.Tau_yld);
%         if itStrainRate>1
%             isPlastic_old = isPlastic_now;
%         end
%         isPlastic_now = all(VAR.Tau_exc>0,2);
%         nel_Plastic   = length(find(isPlastic_now));
%         fprintf(' %1i elements (==%.1f%%) with brittle failure.\n',...
%             nel_Plastic,100*nel_Plastic/MESH.nel);
%         if itStrainRate>1
%             nel_Plastic_new    = length(find(isPlastic_now & ~isPlastic_old));
%             nel_Plastic_nomore = length(find(isPlastic_old & ~isPlastic_now));
%             fprintf(' %1i new elements failed in this iteration.\n',nel_Plastic_new);
%             fprintf(' %1i elements no longer fail in this iteration.\n',nel_Plastic_nomore);
%         end
%     end
    
    if strcmp(SETTINGS.is_plastic,'yes') || ...
       strcmp(SETTINGS.method_eval_visc,'powlaw')
        if itStrainRate > itminStrainRate && converged
            fprintf(' Convergence after %1i iterations :-)\n',itStrainRate);
            break
        end
        
        if itStrainRate == itmaxStrainRate
            fprintf(' No convergence after %1i iterations :-(\n',itmaxStrainRate);
        end
    else
        break % DONE (NO NEED TO ITERATE)
    end
end % for itStrainRate (STRAIN RATE ITERATIONS)

% Save the effective viscosity
VAR.Visc = ViscEl_ip;

VAR.Er_visc=VAR.Tau_II_totl./VAR.Visc_v./2;
VAR.Er_plastic=VAR.Tau_II_totl./VAR.Visc_p./2;

positive=find(VAR.Tau_II_totl>VAR.Tau_yld);
% VAR.Er_plastic(VAR.Tau_II_totl<VAR.Tau_yld)=1e-30;

VAR.Eacc_visc=VAR.Eacc_visc+VAR.Er_visc.*dt;
VAR.Eacc_plastic(positive)=VAR.Eacc_plastic(positive)+VAR.Er_plastic(positive).*dt;




fprintf(1, ' FLOW SOLVER (Penalty, it=%3i/%3i)   : %7.2f sec\n',...
    itStrainRate,itPowHest_totl,toc(t0));

% if strcmp(SETTINGS.is_plastic,'yes')
%     filename = [SETTINGS.outdir '/StrainRateConvergence.png'];
%     sfigure(FigNo_StrainConv);
%     if isfield(SETTINGS,'name_model')
%         ttxt = SETTINGS.name_model; ttxt(ttxt=='_') = ' '; title(ttxt);
%     end
%     print('-dpng','-r200',filename);
% end

end % END OF FUNCTION flow2d_penalty_vep_fs

% #########################################################################
%                               LOCAL FUNCTIONS
% #########################################################################

function VAR = allocate_fields(VAR,MESH,SETTINGS,PHYSICS,nip_stress)

  nel = MESH.nel;

  % Stresses
  if ~isfield(VAR,'Tau_xx')
    VAR.Tau_xx = zeros(nel,nip_stress);
  end
  if ~isfield(VAR,'Tau_xz')
    VAR.Tau_xz = zeros(nel,nip_stress);
  end
  if ~isfield(VAR,'Tau_zz')
    VAR.Tau_zz = zeros(nel,nip_stress);
  end
  if ~isfield(VAR,'Tau_II_totl') % max_shear stress
    VAR.Tau_II_totl = zeros(nel,nip_stress);
  end

  % Strain rates
  if ~isfield(VAR,'Er_xx')
    VAR.Er_xx = zeros(nel,nip_stress);
  end
  if ~isfield(VAR,'Er_xz')
    VAR.Er_xz = zeros(nel,nip_stress);
  end
  if ~isfield(VAR,'Er_zz')
    VAR.Er_zz = zeros(nel,nip_stress);
  end
  if ~isfield(VAR,'Er_II_totl') % max strain rate
    VAR.Er_II_totl = zeros(nel,nip_stress);
  end

  % if strcmp(SETTINGS.is_elastic,'yes')
    % Vorticity
    if ~isfield(VAR,'Vort_xz')
        VAR.Vort_xz = zeros(nel,nip_stress);
    end
  % end

  if strcmp(SETTINGS.is_plastic,'yes')
    % yield stress
    if ~isfield(VAR,'Tau_yld')
        VAR.Tau_yld = zeros(nel,nip_stress);
    end
    if strcmp(SETTINGS.use_strainsoft,'yes')
        if ~isfield(VAR,'Eacc_visc')
            VAR.Eacc_visc = zeros(nel,nip_stress);
        end
        if ~isfield(VAR,'Eacc_plastic')
            VAR.Eacc_plastic = zeros(nel,nip_stress);
        end
        if ~isfield(VAR,'Eacc_II')
            VAR.Eacc_II = zeros(nel,nip_stress);
        end
        if ~isfield(VAR,'Cohesion')
            VAR.Cohesion = repmat(PHYSICS.Cohesion0(MESH.PhaseID)',1,nip_stress);
        end
        if ~isfield(VAR,'PhiFric')
            VAR.PhiFric = repmat(PHYSICS.PhiFric0(MESH.PhaseID)',1,nip_stress);
        end
        if ~isfield(VAR,'F_xx')
            VAR.F_xx = ones(nel,nip_stress);
        end
        if ~isfield(VAR,'F_zz')
            VAR.F_zz = ones(nel,nip_stress);
        end
        if ~isfield(VAR,'F_xz')
            VAR.F_xz = zeros(nel,nip_stress);
        end
        if ~isfield(VAR,'F_zx')
            VAR.F_zx = zeros(nel,nip_stress);
        end
    end
  end 

end % END OF LOCAL FUNCTION allocate_fields

% #########################################################################

function VAR = calc_accumulated_strains(VAR,FAULTS,MESH,U,dt,nelblk)

    % =========================================================================
    % MODEL PARAMETERS
    % =========================================================================
    EL2NOD  = MESH.EL2NOD;    % connectivity matrix for velocity problem
    GCOORD  = MESH.GCOORD;    % node coordinates
    nnodel  = size(EL2NOD,1); % number of nodes in each element
    nel     = size(EL2NOD,2); % number of elements
    ndim    = size(GCOORD,1); % number of spatial dimensions
    nUdofel = ndim*nnodel;  % number of velocity dofs in each element
    EL2DOF                   = zeros(nUdofel, nel);
    EL2DOF(1:ndim:nUdofel,:) = ndim*(EL2NOD-1)+1;
    EL2DOF(2:ndim:nUdofel,:) = ndim*(EL2NOD-1)+2;

    %==========================================================================
    % PREPARE INTEGRATION POINTS & DERIVATIVES wrt LOCAL COORDINATES
    %==========================================================================
    nip_stress = size(VAR.Tau_xx,2);
    % local coordinates of points where stresses are calculated
    [IP_X,~]   = ip_triangle(nip_stress);
    % derivative of velocity shape functions at points where strain rates
    % and stresses are calculated
    [~,dNUds] = sf_dsf_tri367(IP_X,nnodel,'cell');

    % =========================================================================
    % BEGIN OF POST PROCESSING (STRESSES etc)
    % =========================================================================
    nelblk  = min(nelblk,nel);
    nblk    = ceil(nel/nelblk);
    il      = 1;
    iu      = nelblk;

    %==========================================================================
    % BLOCK LOOP - MATRIX COMPUTATION
    %==========================================================================
    for iblk = 1:nblk % LOOP OVER ELEMENT BLOCKS
        ECOORD_x = reshape( GCOORD(1,EL2NOD(:,il:iu)), nnodel, nelblk );
        ECOORD_z = reshape( GCOORD(2,EL2NOD(:,il:iu)), nnodel, nelblk );
        
        %==============================================================
        % STRESS CALCULATION LOOP
        %==============================================================
        Vel_blk = reshape(U(EL2DOF(:,il:iu))',nelblk,nUdofel);
        for ip = 1:nip_stress % LOOP OVER STRESS EVALUATION POINTS
            %==================================================================
            % CALCULATE JACOBIAN, ITS DETERMINANT AND INVERSE
            %==================================================================
            Jx       = ECOORD_x'*dNUds{ip}';
            Jz       = ECOORD_z'*dNUds{ip}';
            detJ     = Jx(:,1).*Jz(:,2) - Jx(:,2).*Jz(:,1);
            if any(detJ<0)
                error('negative Jacobi')
            end
            invdetJ    = 1./detJ;
            invJx      = zeros(nelblk,ndim); % storage for x-components of Jacobi matrix
            invJz      = zeros(nelblk,ndim); % storage for z-components of Jacobi matrix
            invJx(:,1) =  Jz(:,2).*invdetJ;
            invJx(:,2) = -Jz(:,1).*invdetJ;
            invJz(:,1) = -Jx(:,2).*invdetJ;
            invJz(:,2) =  Jx(:,1).*invdetJ;
            
            %==============================================================
            % ACCUMULATED STRAINS AT ip-TH EVALUATION POINT
            %==============================================================
            F0_xx_blk  = VAR.F_xx(il:iu,ip);
            F0_zz_blk  = VAR.F_zz(il:iu,ip);
            F0_xz_blk  = VAR.F_xz(il:iu,ip);
            F0_zx_blk  = VAR.F_zx(il:iu,ip);

            %==========================================================
            % DERIVATIVES wrt GLOBAL COORDINATES
            %==========================================================
            dNUdx      = invJx*dNUds{ip};
            dNUdz      = invJz*dNUds{ip};
            
            %==========================================================
            % STRAIN RATES & ACCUMULATED STRAIN
            % (based on Malvern, 1969, p. ~117)
            %==========================================================
            dUxdx      = sum(Vel_blk(:,1:2:end-1).*dNUdx,2);
            dUxdz      = sum(Vel_blk(:,1:2:end-1).*dNUdz,2);
            dUzdx      = sum(Vel_blk(:,2:2:end  ).*dNUdx,2);
            dUzdz      = sum(Vel_blk(:,2:2:end  ).*dNUdz,2);
            
            % Calculates the rate of change of the deformation gradient
            Fr_xx_blk  = dUxdx.*F0_xx_blk + dUxdz.*F0_zx_blk;
            Fr_xz_blk  = dUxdx.*F0_xz_blk + dUxdz.*F0_zz_blk;
            Fr_zx_blk  = dUzdx.*F0_xx_blk + dUzdz.*F0_zx_blk;
            Fr_zz_blk  = dUzdx.*F0_xz_blk + dUzdz.*F0_zz_blk;
            
            % Calculates the deformation gradient
            F_xx_blk   = F0_xx_blk + dt*Fr_xx_blk;
            F_xz_blk   = F0_xz_blk + dt*Fr_xz_blk;
            F_zx_blk   = F0_zx_blk + dt*Fr_zx_blk;
            F_zz_blk   = F0_zz_blk + dt*Fr_zz_blk;

            % Calculates the accumulated strains
            %   [F]^T      *       [F]       -        I 
            % |Fxx  Fzx|        |Fxx  Fxz|          |1  0|
            % |Fxz  Fzz| times  |Fzx  Fzz|  minus   |0  1|
            Eacc_xx_blk = (1/2)*(F_xx_blk.^2 + F_zx_blk.^2 - 1);
            Eacc_xz_blk = (1/2)*(F_xx_blk.*F_xz_blk + F_zx_blk.*F_zz_blk - 0);
            Eacc_zz_blk = (1/2)*(F_xz_blk.^2 + F_zz_blk.^2 - 1);
            
            % Calculate 2nd invariant of the accumulated strains
            Eacc_II_blk = calc_invariant_II(Eacc_xx_blk, Eacc_zz_blk, Eacc_xz_blk);
            
            %==============================================================
            % WRITE DATA INTO GLOBAL STORAGE
            %==============================================================
            VAR.F_xx   (il:iu,ip) = F_xx_blk;
            VAR.F_zz   (il:iu,ip) = F_zz_blk;
            VAR.F_xz   (il:iu,ip) = F_xz_blk;
            VAR.F_zx   (il:iu,ip) = F_zx_blk;
            VAR.Eacc_II(il:iu,ip) = Eacc_II_blk; % 2nd invariant of accumulated strain
        end % for ip = 1:nip_stress (END OF STRESS EVALUATION LOOP)
        
        %==============================================================
        % READJUST START, END AND SIZE OF NEXT BLOCK
        %==============================================================
        il  = il + nelblk;
        if(iblk==nblk-1)
            nelblk = nel-iu;
        end
        iu  = iu + nelblk;
    end % END OF ELEMENT BLOCK LOOP

    uc_el=find(ismember(MESH.PhaseID,5));
    VAR.Eacc_plastic(uc_el,:)=0;

    for i = 1:length(FAULTS)
        if strcmp(FAULTS(i).activated,'yes') && strcmp(FAULTS(i).deactivated,'no')
            
            fault_nodes=find(ismember(MESH.PointID,FAULTS(i).PointID));
            [fault_nodes2,fault_elements]=find(ismember(MESH.EL2NOD(1:3,:),fault_nodes));
            
            fault_x=MESH.GCOORD(1,fault_nodes);
            [fault_x,order]=sort(fault_x);
            fault_z=MESH.GCOORD(2,fault_nodes);
            fault_z=fault_z(order);
            
            if strcmp(FAULTS(i).dipdir,'left')
                for j=1:length(fault_elements)
                    mean_xz=[mean(MESH.GCOORD(1,MESH.EL2NOD(1:3,fault_elements(j))))...
                        mean(MESH.GCOORD(2,MESH.EL2NOD(1:3,fault_elements(j))))];
                    
                if(mean_xz(1) < interp1(fault_z,fault_x,mean_xz(2),'linear','extrap'))
                    for k=1:3
                            % if(FAULTS(i).PointID~=MESH.PointID(MESH.EL2NOD(k,fault_elements(j))))
                            VAR.Eacc_plastic(fault_elements(j),k)=1;
                            % end
                    end
                end
                
                end
                
            else % FAULTS(i),dipdir=="right"
                for j=1:length(fault_elements)
                    mean_xz=[mean(MESH.GCOORD(1,MESH.EL2NOD(1:3,fault_elements(j))))...
                        mean(MESH.GCOORD(2,MESH.EL2NOD(1:3,fault_elements(j))))];
                    
                if(mean_xz(1)> interp1(fault_z,fault_x,mean_xz(2),'linear','extrap'))
                    for k=1:3
                            % if(FAULTS(i).PointID~=MESH.PointID(MESH.EL2NOD(k,fault_elements(j))))
                                VAR.Eacc_plastic(fault_elements(j),k)=1;
                            % end
                    end
                end
                
                end % for j
            end % FAULTS(i),dipdir=="right"
        
        end
    end % for FAULTS 

end % END LOCAL FUNCTION calc_accumulated_strains

% #########################################################################

function VAR = calc_stresses(VAR,ViscEl_ip,MESH,U,SETTINGS,PHYSICS,NUMSCALE,dt,nelblk,SERP)

  % =========================================================================
  % MODEL PARAMETERS
  % =========================================================================
  EL2NOD  = MESH.EL2NOD;    % connectivity matrix for velocity problem
  GCOORD  = MESH.GCOORD;    % node coordinates
  PhaseID = MESH.PhaseID;   % element Phase-IDs
  if length(PhaseID)==1
    PhaseID = PhaseID * ones(1,MESH.nel,'int32');
  end
  nnodel  = size(EL2NOD,1); % number of nodes in each element
  nel     = size(EL2NOD,2); % number of elements
  ndim    = size(GCOORD,1); % number of spatial dimensions
  EL2NODP = MESH.EL2NODP;   % connectivity matrix for pressure problem
  nvertx  = 3;
  nPdofel = size(EL2NODP,1);
  nUdofel = ndim*nnodel;  % number of velocity dofs in each element
  EL2DOF                   = zeros(nUdofel, nel);
  EL2DOF(1:ndim:nUdofel,:) = ndim*(EL2NOD-1)+1;
  EL2DOF(2:ndim:nUdofel,:) = ndim*(EL2NOD-1)+2;

  %==========================================================================
  % PREPARE INTEGRATION POINTS & DERIVATIVES wrt LOCAL COORDINATES
  %==========================================================================
  nip_stress = size(VAR.Tau_xx,2);
  % local coordinates of points where stresses are calculated
  [IP_X,~]   = ip_triangle(nip_stress);
  % derivative of velocity shape functions at points where strain rates
  % and stresses are calculated
  [NU,dNUds] = sf_dsf_tri367(IP_X,nnodel,'cell');
  % pressure shape functions (their derivatives are not needed)
  [NP,~]     = sf_dsf_tri367(IP_X,nPdofel,'cell');
  % linear (3-node) shape functions; can be used to interpolate density
  [NL,~]     = sf_dsf_tri367(IP_X,nvertx,'cell');

  % =========================================================================
  % BEGIN OF POST PROCESSING (STRESSES etc)
  % =========================================================================
  nelblk  = min(nelblk,nel);
  nblk    = ceil(nel/nelblk);
  il      = 1;
  iu      = nelblk;

  %==========================================================================
  % BLOCK LOOP - MATRIX COMPUTATION
  %==========================================================================
  for iblk = 1:nblk % LOOP OVER ELEMENT BLOCKS
    ECOORD_x = reshape( GCOORD(1,EL2NOD(:,il:iu)), nnodel, nelblk );
    ECOORD_z = reshape( GCOORD(2,EL2NOD(:,il:iu)), nnodel, nelblk );
    
    %==============================================================
    % STRESS CALCULATION LOOP
    %==============================================================
    Vel_blk = reshape(U(EL2DOF(:,il:iu))',nelblk,nUdofel);                           % [nelblk, ndim*nnodel] :: (:, v_1_x v_1_z ... v_nnodel_z)
    for ip = 1:nip_stress % LOOP OVER STRESS EVALUATION POINTS
        %==============================================================
        % VISCOSITY OF ELEMENTS AT ip-TH EVALUATION POINT
        %==============================================================
        ViscIP_blk = ViscEl_ip(il:iu,ip);
        
        %==================================================================
        % CALCULATE JACOBIAN, ITS DETERMINANT AND INVERSE
        %==================================================================
        Jx       = ECOORD_x'*dNUds{ip}';                    % [nelblk,2]
        Jz       = ECOORD_z'*dNUds{ip}';                    % [nelblk,2]
        detJ     = Jx(:,1).*Jz(:,2) - Jx(:,2).*Jz(:,1);     % [nelblk,1]
        if any(detJ<0)
            error('negative Jacobian')
        end
        invdetJ    = 1./detJ;
        invJx      = zeros(nelblk,ndim); % storage for x-components of Jacobi matrix
        invJz      = zeros(nelblk,ndim); % storage for z-components of Jacobi matrix
        invJx(:,1) =  Jz(:,2).*invdetJ;
        invJx(:,2) = -Jz(:,1).*invdetJ;
        invJz(:,1) = -Jx(:,2).*invdetJ;
        invJz(:,2) =  Jx(:,1).*invdetJ;
        
        %==========================================================
        % DERIVATIVES wrt GLOBAL COORDINATES
        %==========================================================
        dNUdx      = invJx*dNUds{ip};                                 % [nelblk,nnodel]
        dNUdz      = invJz*dNUds{ip};                                 % [nelblk,nnodel]
        dUxdx      = sum(Vel_blk(:,1:2:end-1).*dNUdx,2);              % [nelblk,1]  == Er_xx_block
        dUxdz      = sum(Vel_blk(:,1:2:end-1).*dNUdz,2);              % [nelblk,1]
        dUzdx      = sum(Vel_blk(:,2:2:end  ).*dNUdx,2);              % [nelblk,1]
        dUzdz      = sum(Vel_blk(:,2:2:end  ).*dNUdz,2);              % [nelblk,1]  == Er_zz_block
        Er_xz_blk  = 0.5*(dUxdz + dUzdx);                             % [nelblk,1]  == Er_xz_block = Er_zx_block
        % Second strain rate invariant / effective total strain rate
        Er_II_blk   = calc_invariant_II(dUxdx, dUzdz, Er_xz_blk);
        
        %==================================================================
        % VISCOUS STRESSES
        %==================================================================
        Tau_xx_blk =  4/3.*ViscIP_blk.*dUxdx - 2/3.*ViscIP_blk.*dUzdz;
        Tau_zz_blk = -2/3.*ViscIP_blk.*dUxdx + 4/3.*ViscIP_blk.*dUzdz;
        Tau_xz_blk =  2  .*ViscIP_blk.*Er_xz_blk; % *2 because Er_xz_blk = 0.5*...
        
        if strcmp(SETTINGS.is_elastic,'yes')
            %==============================================================
            % CALC VORTICITY
            % (required for Jaumann stress rotation in NEXT time step)
            %==============================================================
            Vort_xz_blk = 0.5*(dUxdz - dUzdx);         % clockwise positive. JGP: standard angular velocity (positive counterclockwise) would be omega = 0.5(dUzdx - dUxdz)

            %==============================================================
            % TOTAL STRESSES ARE VISCOUS STRESSES PLUS REMAINING OLD
            % ELASTIC STRESSES (decay facor "Xi")
            %==============================================================
            Xi         = ViscIP_blk./(PHYSICS.ShearG(PhaseID(il:iu))'.*dt);
            Tau_xx_blk = Tau_xx_blk + Xi.*VAR.Tau_xx_old(il:iu,ip);
            Tau_zz_blk = Tau_zz_blk + Xi.*VAR.Tau_zz_old(il:iu,ip);
            Tau_xz_blk = Tau_xz_blk + Xi.*VAR.Tau_xz_old(il:iu,ip);
        end
        
        % Maximum shear stresses (2nd invariant of deviatoric stress tensor)
        Tau_II_blk  = calc_invariant_II(Tau_xx_blk, Tau_zz_blk, Tau_xz_blk);
        
        if strcmp(SETTINGS.is_plastic,'yes')
            %==============================================================
            % PLASTICITY [STRAIN SOFTENING] PART
            %==============================================================
            
            % Plasticity parameters for lithologic phases
            Phi_blk = PHYSICS.PhiFric0(PhaseID(il:iu))';  % friction angle
            Coh_blk = PHYSICS.Cohesion0(PhaseID(il:iu))'; % cohesion
            
            if strcmp(SETTINGS.use_strainsoft,'yes')
                SS         = PHYSICS.SS;
                I2         = VAR.Eacc_II(il:iu,ip);

                I2p       = VAR.Eacc_plastic(il:iu,ip);
                I2v       = VAR.Eacc_visc(il:iu,ip);
                
                % =============== Friction angle shoftening ===============
                Phi_ss_blk = Phi_blk; % make a copy for testing

                % Load limits for friction angle softening
                Phi_ss_max = PHYSICS.PhiFric0(PhaseID(il:iu))';
                Phi_ss_min = min(SS.PhiFric_lims);
                I2_Phi_min = min(SS.Eac_II_lims_PhiFric);
                I2_Phi_max = max(SS.Eac_II_lims_PhiFric);

                % Calculate "Phi" at integration points where "I2" is 
                % within the limits; i.e. I2_Phi_min <= I2 <= Phi_ss_max
                ind_ss_Phi = I2p >= I2_Phi_min & I2p <= I2_Phi_max;

                % Phi = (I2-rangeI2_1)*rangePhi/rangeI2+rangePhi_1
                % A   = (B -C        )*D       /E      +F
                Phi_ss_blk(ind_ss_Phi) = ...              % A
                    (-I2p(ind_ss_Phi)+I2_Phi_max) .* ...    % B-C
                    (Phi_ss_max(ind_ss_Phi)-Phi_ss_min) ./ ...        % D
                    (I2_Phi_max-I2_Phi_min) + Phi_ss_min; % E+F

                % Calculate Phi at integration points where "I2" is larger
                % than the upper limit; set to minimum value "Phi_ss_min"
                Phi_ss_blk(I2p > I2_Phi_max) = Phi_ss_min;


                % ================== Cohesion shoftening ==================
                Coh_ss_blk = Coh_blk; % make a copy for testing

                % Load limits for friction angle softening
                Coh_ss_max = PHYSICS.Cohesion0(PhaseID(il:iu))';
                Coh_ss_min = min(SS.Cohesion_lims);
                I2_Coh_min = min(SS.Eac_II_lims_Cohesion);
                I2_Coh_max = max(SS.Eac_II_lims_Cohesion);

                % Cslculate "Coh" at integration points where "I2" is 
                % within the limits; i.e. I2_Phi_min <= I2 <= Phi_ss_max
                ind_ss_Coh = I2p >= I2_Coh_min & I2p <= I2_Coh_max;

                % Cohesion = (I2-rangeI2_1)*rangeC/rangeI2+rangeC_1
                % A        = (B -C        )*D     /E      +F
                
                Coh_ss_blk(ind_ss_Coh) = ...              % A
                    (-I2p(ind_ss_Coh)+I2_Coh_max) .* ...    % B-C
                    (Coh_ss_max(ind_ss_Coh)-Coh_ss_min) ./ ...        % D
                    (I2_Coh_max-I2_Coh_min) + Coh_ss_min; % E+F

                % Calculate Coh at integration points where "I2" is larger
                % than the upper limit; set to minimum value "Coh_ss_min"
                Coh_ss_blk(I2p > I2_Coh_max) = Coh_ss_min;
                
                
                % =========================================================
                % Only update the phases selected for strain softening 
                % =========================================================
                iel_SS          = ismember(PhaseID(il:iu),SS.PhasesID_ss)';
                Phi_blk(iel_SS) = Phi_ss_blk(iel_SS);
                Coh_blk(iel_SS) = Coh_ss_blk(iel_SS);
                %==========================================================
                % update Phi for serpentinization mantle
                %==========================================================
                if strcmp(SETTINGS.SERPENTINISATION,'yes')
                    Dserp = VAR.Dserp(MESH.EL2NOD(1:3,il:iu))' * NP{ip};
                    LA = find(PhaseID(il:iu)==1 | PhaseID(il:iu)==2);
                    Phi_blk(LA,:) = Phi_blk(LA,:).*(1-SERP.Phi_alpha.*Dserp(LA,:));
                    Phi_blk = max(Phi_blk,SERP.Phi_min);
                end
            end % if use_strainsoft

            % Lithostatic pressure at IP (must be in units of NUMSCALE.P0!
            if isfield(VAR,'P_lithstat')
                Pl_blk = VAR.P_lithstat(EL2NODP(:,il:iu))' * NP{ip};
                %                 Z_ip  = reshape(GCOORD(2,EL2NOD(1:3,il:iu)),3,[])' * NP{ip}; % depth of IP
                %                 figure(88);clf;plot(Pl_blk,Z_ip)
            else
                error('Lithostatic pressure must be provided.');
            end
            
            % DEFINITION OF YIELD FUNCTION; Spiegelman, May, Wilson (2016)
            % Tau_yld = A + B * (P_lith + alpha*P_dyn)
            switch SETTINGS.plasticity_model
                case 'VM' % "von Mises" plasticity model
                    A      = Coh_blk;
                    B      = 0;
                    alpha  = 0;
                    Pl_blk = 0;
                    Pd_blk = 0;
                    
                case 'DDM' % depth-dependent "von Mises"
                    A      = Coh_blk .* cosd(Phi_blk);
                    B      = sind(Phi_blk);
                    alpha  = 0;
                    Pd_blk = 0;

                case 'DP'
                    A     = Coh_blk .* cosd(Phi_blk);
                    B     = sind(Phi_blk);
                    alpha = 1;
                    
                    if NUMSCALE.Dens0~=0
                        error('This code part needs to be modified for NUMSCALE.Dens0~=0');
                    end
                    if ~strcmp(SETTINGS.top_surface,'free')
                        error('Can only use dynamic pressure when top surface is free.');
                    end
                    
                    % Interpolate dynamic pressure integration point
                    Pd_blk = VAR.P(EL2NODP(:,il:iu));
                    Pd_blk = Pd_blk' * NP{ip};
                    %Pd_blk = sum(Pd_blk' .* NP{ip},2);

                otherwise
                    error('SETTINGS.plasticity_model must be "VM", "DDM", or "DP".');
            end  % SETTINGS.plasticity_model
            
            % Yield stress
            Tau_yld_blk = A + B.*((1-alpha).*Pl_blk + alpha.*Pd_blk);
            
            if isfield(PHYSICS,'Tau_yld_min') && ~isempty(PHYSICS.Tau_yld_min)
                Tau_yld_blk = max(Tau_yld_blk,PHYSICS.Tau_yld_min); % min cut-off
            end
            if isfield(PHYSICS,'Tau_yld_max') && ~isempty(PHYSICS.Tau_yld_max)
                Tau_yld_blk = min(Tau_yld_blk,PHYSICS.Tau_yld_max); % max cut-off
            end
        end % if SETTINGS.is_plastic
        
        %==============================================================
        % WRITE DATA INTO GLOBAL STORAGE
        %==============================================================
        VAR.Er_xx      (il:iu,ip) = dUxdx;
        VAR.Er_zz      (il:iu,ip) = dUzdz;
        VAR.Er_xz      (il:iu,ip) = Er_xz_blk;
        VAR.Er_II_totl (il:iu,ip) = Er_II_blk; % 2nd strain rate invariant / effective total strain rate
        VAR.Tau_xx     (il:iu,ip) = Tau_xx_blk;
        VAR.Tau_zz     (il:iu,ip) = Tau_zz_blk;
        VAR.Tau_xz     (il:iu,ip) = Tau_xz_blk;
        VAR.Tau_II_totl(il:iu,ip) = Tau_II_blk; % 2nd stress invariant
        if strcmp(SETTINGS.is_elastic,'yes')
            VAR.Vort_xz(il:iu,ip) = Vort_xz_blk;
        end
        if strcmp(SETTINGS.is_plastic,'yes')
            VAR.Tau_yld(il:iu,ip) = Tau_yld_blk; % Current yield stress
            if strcmp(SETTINGS.use_strainsoft,'yes')
                VAR.Cohesion(il:iu,ip) = Coh_blk;
                VAR.PhiFric (il:iu,ip) = Phi_blk;
            end
        end
    end % for ip = 1:nip_stress END OF STRESS EVALUATION LOOP
    
    %==============================================================
    % READJUST START, END AND SIZE OF NEXT BLOCK
    %==============================================================
    il  = il + nelblk;
    if(iblk==nblk-1)
        nelblk = nel-iu;
    end
    iu  = iu + nelblk;
  end % END OF ELEMENT BLOCK LOOP

end % END LOCAL FUNCTION calc_stresses

% #########################################################################

function DensEl_ip = element_density(VAR, SETTINGS, PHYSICS, EL2NOD, PhaseID, els, Nip)
  % obtain density at one specific integration point for a subset of elements
  %
  % VAR.Dens        :: REAL, OPTIONAL [nnodel,nel], required for 'interp_nodal'
  %                                   [nel],          "          'elem_var'
  %                                   [nel,nnodel]    "          'nodal_disc', in this case, it is assumed that nnodel=3
  % SETTINGS.method_eval_dens:: CHARACTER in {'interp_nodal','elem_phases','elem_var','nodal_disc'}
  % PHYSICS.Dens:: REAL, OPTIONAL [nphases], required for 'elem_phases'
  %        .minDens :: REAL, OPTIONAL, lower density cut-off
  %        .maxDens :: REAL, OPTIONAL, upper density cut-off
  % EL2NOD          :: REAL [nnodel,nel]
  % PhaseID         :: INTEGER [nel]
  % els             :: INTEGER [nelblk] vector, indicating a subset of elements within EL2NOD
  % Nip             :: REAL [nnodel,1] local shape functions evaluated at the integration point

  % output
  % DensEl_ip       :: REAL [nelblk,1] density evaluated at the integration point for all elements in the block
  
  nnodel = length(Nip);
  switch SETTINGS.method_eval_dens
    case 'interp_nodal'
        % interpolation from nodal values to the integration point
        DensEl_ip = VAR.Dens( EL2NOD(1:nnodel,els) )'*Nip;              % [nelblk,1]
        
    case 'elem_phases'
        % use the constant densities defined in PHYSICS.Dens
        DensEl_ip = PHYSICS.Dens(PhaseID(els))';                        % [nelblk,1]
        
    case 'elem_var'
        DensEl_ip = VAR.Dens(els);                                      % [nelblk,1]
        
    case 'nodal_disc'
        % Density is defined at vertex nodes of each element
        DensEl_ip = VAR.Dens(els,:)*Nip;                                % [nelblk,1]
        
    otherwise
        error(' Unknown "case" for evaluating density at integration points.');
  end

  % Make sure a column-vector is returned
  DensEl_ip = DensEl_ip(:);

  % DENSITY CUT-OFFS (IF DEFINED)
  if isfield(PHYSICS,'minDens') && ~isempty(PHYSICS.minDens)
    DensEl_ip = max(DensEl_ip,PHYSICS.minDens); % lower cut-off for density
  end
  if isfield(PHYSICS,'maxDens') && ~isempty(PHYSICS.maxDens)
    DensEl_ip = min(DensEl_ip,PHYSICS.maxDens); % upper cut-off for density
  end

end % END LOCAL FUNCTION element_density

% #########################################################################

function [KK,ViscEl_sip,VAR,GG,invMM,Rhs,DensEl] = matrix_assembly...
    (VAR,MESH,SETTINGS,PHYSICS,NUMSCALE,dt,fpenalty,nelblk,nip_flow,itStrainRate)

    % =========================================================================
    % MODEL PARAMETERS
    % =========================================================================
    EL2NOD     = MESH.EL2NOD;    % connectivity matrix for velocity problem
    GCOORD     = MESH.GCOORD;    % node coordinates
    PhaseID    = MESH.PhaseID;   % element Phase-IDs
    nnodel     = size(EL2NOD,1); % number of nodes in each element
    nel        = size(EL2NOD,2); % number of elements
    ndim       = size(GCOORD,1); % number of spatial dimensions
    nnod       = size(GCOORD,2); % number of nodes
    EL2NODP    = MESH.EL2NODP;   % connectivity matrix for pressure problem
    nUdofel    = ndim*nnodel;  % number of velocity dofs in each element
    nPdofel    = size(EL2NODP,1); % number of pressure dofs in each element
    nP         = max(EL2NODP(:));
    nU         = ndim*nnod;
    nip_stress = size(VAR.Tau_xx,2);
    penalty    = PHYSICS.maxVisc * fpenalty;

    %==========================================================================
    % PREPARE INTEGRATION POINTS & DERIVATIVES wrt LOCAL COORDINATES
    %==========================================================================
    [IP_X,IP_w]      = ip_triangle(nip_flow);
    % local coordinates and weights of points for integration of
    % velocity/pressure matrices
    [NU,dNUds]       = sf_dsf_tri367(IP_X,nnodel,'cell');
    % velocity shape functions and their derivatives
    [NP,~]           = sf_dsf_tri367(IP_X,nPdofel,'cell');
    % pressure shape functions (their derivatives are not needed)
    [IP_X_stress,IP_w_stress] = ip_triangle(nip_stress);
    % local coordinates of points where stresses are calculated
    [~,dNUds_stress] = sf_dsf_tri367(IP_X_stress,nnodel,'cell');
    % derivative of velocity shape functions at points where strain rates
    % and stresses are caclculated
    [NL,~]           = sf_dsf_tri367(IP_X,3,'cell');
    % linear (3-node) shape functions; can be used to interpolate density
    [N3_stress,~]    = sf_dsf_tri367(IP_X_stress,3,'cell');
    % linear (3-node) shape functions at stress evaluation points
    [N6_stress,~]    = sf_dsf_tri367(IP_X_stress,6,'cell');
    % quadratic (6-node) shape functions at stress evaluation points

    %==========================================================================
    % STORAGE FOR DATA OF ALL ELEMENT MATRICES/VECTORS
    %==========================================================================
    K_all = zeros(nUdofel*(nUdofel+1)/2, nel); % storage for stiffness matrix data, incl penalty terms
    if itStrainRate == 1
        G_all    = zeros(nPdofel*nUdofel, nel); % storage for gradient matrix data
        Rhs_all  = zeros(nUdofel        , nel); % storage for global force vector
        invM_all = zeros(nPdofel*nPdofel, nel); % storage for inverse of global mass matrix
        DensEl   = zeros(nel,1);
    end
    
    %==========================================================================
    % BLOCKING PARAMETERS (nelblk must be < nel)
    %==========================================================================
    nelblk  = min(nel,nelblk);  % in case nel<nelblk
    nblk    = ceil(nel/nelblk); % number of blocks
    il      = 1;
    iu      = nelblk;

    if isfield(VAR,'DilEl') && any(VAR.DilEl~=0)
        error(' This version does not support dilation terms. Remove VAR.DilEl or set to zero.');
    end

    %==========================================================================
    % VISCOSITY AT ALL STRESS-EVALUATION POINTS OF ALL ELEMENTS
    %==========================================================================
    [ViscEl_sip,VAR.Visc_v,VAR.Visc_e,VAR.Visc_p] = calc_visc_vep_all_ip(VAR,SETTINGS,PHYSICS,NUMSCALE,...
        GCOORD,EL2NOD,EL2NODP,PhaseID,1:nel,N3_stress,N6_stress,dt);
    % plot_2d_fedata(71,GCOORD,EL2NOD,log10(NUMSCALE.Visc0*ViscEl_sip),[],[],'none',true,MODEL);title('log Visc');

    switch SETTINGS.visc_in_el
        case 'min'
            % Minimum value
            ViscEl_fip = min(ViscEl_sip,[],2);
            % viscosity at integration points used during matrix assembly
            ViscEl_sip = repmat(ViscEl_fip,1,nip_stress);
            % viscosity at stress evaluation points
            
        case 'mean'
            % Harmonic mean
            ViscEl_fip = nip_stress./sum(1./ViscEl_sip,2);
            % viscosity at integration points used during matrix assembly
            ViscEl_sip = repmat(ViscEl_fip,1,nip_stress);
            % viscosity at stress evaluation points
            
        case 'linear'
            % Linear variation over element
            % (1) Map from stress-evaluation points to vertex nodes
            NN          = [N3_stress{:}];
            ViscEl_vnod = (NN \ ViscEl_sip')';
            % (2) Interpolate from vertex nodes to integration points
            ViscEl_fip  = zeros(nel,nip_flow);
            for ip=1:nip_flow
                ViscEl_fip(:,ip) = ViscEl_vnod*NL{ip};
                % viscosity at integration points used during matrix assembly
            end
            % ViscEl_sip = ViscEl_sip;
            % viscosity at stress evaluation points (does not change)
            
            %         % FOR CHECKING THE ABOVE NN \ ...
            %         ViscEl_nd  = zeros(nelblk,nvertx);
            %         for iel=1:nelblk
            %             ViscEl_nd(iel,:) = (NN \ ViscEl_ip(il+iel-1,:)')';
            %         end
            % 
            %         % FOR PLOTTING VALUES AT INTEGRATION POINTS AND NODES
            %         figure(666);clf;
            %         iel = 46;
            %         patch(ECOORD_x(1:3,iel),ECOORD_z(1:3,iel),ViscEl_nd(iel,:)');
            %         hold on
            %         for ip_stress=1:nip_stress
            %             X_ip  = N3_stress{ip_stress}' * ECOORD_x(1:3,iel);
            %             Z_ip  = N3_stress{ip_stress}' * ECOORD_z(1:3,iel);
            %             scatter(X_ip,Z_ip,200,ViscEl_ip(iel,ip_stress),'filled','MarkerEdgeColor','k');
            %         end
            %         for ip_flow=1:nip_flow
            %             X_ip  = NL{ip_flow}' * ECOORD_x(1:3,iel);
            %             Z_ip  = NL{ip_flow}' * ECOORD_z(1:3,iel);
            %             scatter(X_ip,Z_ip,200,ViscEl_nd(iel,:)*NL{ip_flow},'filled','Marker','s','MarkerEdgeColor','k');
            %         end
            %         colorbar

            % MIN/MAX CUT-OFFS ARE NEEDED WHEN MAPPING FROM 
            % 3 IPs --> NODES --> 7 IPs
            if isfield(PHYSICS,'minVisc') && ~isempty(PHYSICS.minVisc)
                ViscEl_fip = max(ViscEl_fip,PHYSICS.minVisc); % lower cut-off for viscosity
            end
            if isfield(PHYSICS,'maxVisc') && ~isempty(PHYSICS.maxVisc)
                ViscEl_fip = min(ViscEl_fip,PHYSICS.maxVisc); % upper cut-off for viscosity
            end

        case 'any'
            % Any variation over element (simply use stresses evaluated
            % at velocity integration points)
            ViscEl_fip = ViscEl_sip;

        otherwise
            error('SETTINGS.visc_in_el must be one of "min", "mean", "linear", or "any"');
    end % switch SETTINGS.visc_in_el 

  %==========================================================================
  % BLOCK LOOP - MATRIX COMPUTATION
  %==========================================================================
  for iblk = 1:nblk
    ECOORD_x = reshape( GCOORD(1,EL2NOD(:,il:iu)), nnodel, nelblk );
    ECOORD_z = reshape( GCOORD(2,EL2NOD(:,il:iu)), nnodel, nelblk );
    
    %======================================================================
    % STORAGE FOR DATA OF ELEMENTS IN BLOCK
    %======================================================================
    K_blk     = zeros(nelblk, nUdofel*(nUdofel+1)/2);
      % symmetric stiffness matrix, dim=vel dofs, but only upper triangle
    M_blk     = zeros(nelblk, nPdofel*(nPdofel+1)/2);
      % symmetric mass matrix, dim=pressure dofs, but only upper triangle
    invM_blk  = zeros(nelblk, nPdofel*nPdofel);
      % storage for inv(M) matrix
    invMG_blk = zeros(nelblk, nPdofel*nUdofel);
      % storage for inv(M)*G' matrix product
    G_blk     = zeros(nelblk, nPdofel*nUdofel);
      % asymmetric gradient matrix, pressure dofs x vel dofs
    if itStrainRate == 1
        Rhs_blk = zeros(nelblk, nUdofel);
          % storage for right-hand-side vector
    end
    
    %==============================================================
    % INTEGRATION LOOP: DEAL WITH A BLOCK OF ELEMENTS AT ONCE
    %==============================================================
    for ip = 1:nip_flow
        %==============================================================
        % PROPERTIES OF ELEMENTS AT ip-TH INTEGRATION POINT
        %==============================================================
        if itStrainRate==1
            Dens_blk = element_density...                                  % [nelblk,1]
                (VAR, SETTINGS, PHYSICS, EL2NOD, PhaseID, il:iu, NL{ip});  % *LOCAL FUNCTION*
            % Calculate average density in each element
            % (used for free surface boundary condition)
            DensEl(il:iu) = DensEl(il:iu) + Dens_blk(:)./nip_flow;         % [nelblk,1] :: JGP: why this approximation? It should better be a weighted mean from the integration point values 
            
            % Gravitational force at ip-th integration point
            Fg_blk = PHYSICS.g * NUMSCALE.Bscale .* (Dens_blk-NUMSCALE.Dens0);
        end
        
        if size(ViscEl_fip,2)==nip_flow
            Visc_blk = ViscEl_fip(il:iu,ip);
        else
            Visc_blk = ViscEl_fip(il:iu);
        end
        
        %==================================================================
        % CALCULATE JACOBIAN, ITS DETERMINANT AND INVERSE
        %==================================================================
        Jx       = ECOORD_x'*dNUds{ip}';
        Jz       = ECOORD_z'*dNUds{ip}';
        detJ     = Jx(:,1).*Jz(:,2) - Jx(:,2).*Jz(:,1);
        
        if any(detJ<0)
            plot_mesh(GCOORD,EL2NOD,666,'k'); hold on
            ind_bad = find(detJ<0);
            for iel=ind_bad
                scatter(mean(ECOORD_x(:,iel)),mean(ECOORD_z(:,iel)),10,'r');
                scatter(mean(ECOORD_x(:,iel)),mean(ECOORD_z(:,iel)),200,'r');
                for inod=1:nnodel
                    text(ECOORD_x(inod,iel),ECOORD_z(inod,iel),num2str(inod));
                end
            end
            error('negative Jacobian')
        end
        
        invdetJ    = 1./detJ;
        invJx      = zeros(nelblk,ndim); % storage for x-components of Jacobi matrix
        invJz      = zeros(nelblk,ndim); % storage for z-components of Jacobi matrix
        invJx(:,1) =  Jz(:,2).*invdetJ;
        invJx(:,2) = -Jz(:,1).*invdetJ;
        invJz(:,1) = -Jx(:,2).*invdetJ;
        invJz(:,2) =  Jx(:,1).*invdetJ;
        
        %==========================================================
        % DERIVATIVES wrt GLOBAL COORDINATES
        %==========================================================
        dNUdx    = invJx*dNUds{ip};
        dNUdz    = invJz*dNUds{ip};
        
        %==================================================================
        % PRESSURE SHAPE FUNCTIONS
        %==================================================================
        NP_blk = repmat(NP{ip}',nelblk,1);
        
        %==========================================================
        % NUMERICAL INTEGRATION OF ELEMENT MATRICES
        %==========================================================
        w_blk   = detJ*IP_w(ip); % integration weight x area of triangle
        
        % STIFFNESS MATRICES FOR ALL ELEMENTS IN BLOCK
        C1      =  4/3; % Used instead of 'D' matrix in standard assembly
        C2      = -2/3; % see Zienkiewicz book, Vol 2, 4th edition, p 519
        w2_blk  = Visc_blk.*w_blk;
          % "weight" times viscosity at integration point in all elements
        indx = 1;
        for i=1:nnodel
            % x-velocity equation (1st, 3rd, 5th,... row of stiffness matrices)
            for j=i:nnodel
                % x-velocity (1st, 3th, 5th,... column)
                K_blk(:,indx) = K_blk(:,indx) + w2_blk .* ...
                    (C1 * dNUdx(:,i) .* dNUdx(:,j) + dNUdz(:,i) .* dNUdz(:,j));
                indx          = indx + 1;
                % z-velocity (2nd, 4th, 6th,... column)
                K_blk(:,indx) = K_blk(:,indx) + w2_blk .* ...
                    (C2 * dNUdx(:,i) .* dNUdz(:,j) + dNUdz(:,i) .* dNUdx(:,j));
                indx          = indx + 1;
            end
            % z-velocity equation (2nd, 4th, 6th,... row of stiffness matrices)
            for j=i:nnodel
                if j>i
                    %  x-velocity (3rd, 5th, 7th,... column)
                    K_blk(:,indx) = K_blk(:,indx) + w2_blk .* ...
                        (C2 * dNUdz(:,i) .* dNUdx(:,j) + dNUdx(:,i) .* dNUdz(:,j));
                    indx          = indx + 1;
                end
                % z-velocity (2nd, 4th, 6th,... column)
                K_blk(:,indx) = K_blk(:,indx) + w2_blk .* ...
                    (C1 * dNUdz(:,i) .* dNUdz(:,j) + dNUdx(:,i) .* dNUdx(:,j));
                indx          = indx + 1;
            end
        end
        
        % MASS MATRICES FOR ALL ELEMENTS IN BLOCK
        indx = 1;
        for i=1:nPdofel
            for j=i:nPdofel
                M_blk(:,indx) = M_blk(:,indx) + w_blk .* NP_blk(:,i).*NP_blk(:,j);
                indx = indx + 1;
            end
        end
        
        % GRADIENT MATRICES FOR ALL ELEMENTS IN BLOCK
        for i=1:nPdofel
            tmp1        = w_blk.*NP_blk(:,i);
            tmp2        = tmp1(:,ones(1,nnodel));
            ii          = (i-1)*nUdofel + (1:2:nUdofel);
            G_blk(:,ii) = G_blk(:,ii) + tmp2.*dNUdx; % MILAMIN has minus sign
            ii          = (i-1)*nUdofel + (2:2:nUdofel);
            G_blk(:,ii) = G_blk(:,ii) + tmp2.*dNUdz; % MILAMIN has minus sign
        end
        
        if itStrainRate==1
            % ASSEMBLY OF BUOYANCY FORCES FOR ALL ELEMENTS IN BLOCK
            Rhs_blk(:,2:2:nUdofel) = Rhs_blk(:,2:2:nUdofel) + (w_blk .* Fg_blk) * NU{ip}';
        end
    end % END OF INTEGRATION LOOP
    
    %======================================================================
    % CALCULATE R.H.S.-TERMS FROM ELASTIC STRESSES (done in separate loop 
    % because another number of integration points is used)
    %======================================================================
    if strcmp(SETTINGS.is_elastic,'yes') && dt>0 && itStrainRate==1
        Rhs_e_blk = zeros(nelblk, nUdofel);
        for ip = 1:nip_stress
            %==============================================================
            % CALCULATE JACOBIAN, ITS DETERMINANT AND INVERSE
            %==============================================================
            Jx       = ECOORD_x'*dNUds_stress{ip}';
            Jz       = ECOORD_z'*dNUds_stress{ip}';
            detJ     = Jx(:,1).*Jz(:,2) - Jx(:,2).*Jz(:,1);
            if any(detJ<0)
                figure(SETTINGS.Fig_number);hold on
                scatter(mean(ECOORD_x(:,detJ<0)),mean(ECOORD_z(:,detJ<0)),50,'w','filled');
                error('negative Jacobian')
            end
            invdetJ    = 1./detJ;
            invJx      = zeros(nelblk,ndim); % storage for x-components of Jacobi matrix
            invJz      = zeros(nelblk,ndim); % storage for z-components of Jacobi matrix
            invJx(:,1) =  Jz(:,2).*invdetJ;
            invJx(:,2) = -Jz(:,1).*invdetJ;
            invJz(:,1) = -Jx(:,2).*invdetJ;
            invJz(:,2) =  Jx(:,1).*invdetJ;
            
            %==============================================================
            % DERIVATIVES wrt GLOBAL COORDINATES
            %==============================================================
            dNUdx      = invJx*dNUds_stress{ip};
            dNUdz      = invJz*dNUds_stress{ip};
            
            % new weight for integration as number if ips changed
            ws_blk     = detJ*IP_w_stress(ip);
            
            %==============================================================
            % ADD ELASTIC STRESSES FROM PREVIOUS TIME STEP TO RHS
            %==============================================================
            % Decay-factor "Xi" for elastic stresses
            Visc_blk = ViscEl_sip(il:iu,ip);
            Xi_blk   = Visc_blk./(PHYSICS.ShearG(PhaseID(il:iu))'.*dt);

            Rhs_e_blk(:,1:2:nUdofel) = Rhs_e_blk(:,1:2:nUdofel) + ...
                - dNUdx.*repmat(Xi_blk.*VAR.Tau_xx_old(il:iu,ip).*ws_blk,[1,nnodel])...
                - dNUdz.*repmat(Xi_blk.*VAR.Tau_xz_old(il:iu,ip).*ws_blk,[1,nnodel]);
            Rhs_e_blk(:,2:2:nUdofel) = Rhs_e_blk(:,2:2:nUdofel) + ...
                - dNUdz.*repmat(Xi_blk.*VAR.Tau_zz_old(il:iu,ip).*ws_blk,[1,nnodel])...
                - dNUdx.*repmat(Xi_blk.*VAR.Tau_xz_old(il:iu,ip).*ws_blk,[1,nnodel]);
        end % END OF ELASTIC STRESS LOOP
        
        % Add elastic part to right-hand-side 
        Rhs_blk = Rhs_blk + Rhs_e_blk;
    end % if SETTINGS.is_elastic && dt > 0 && itStrainRate==1    
    
    %======================================================================
    % CALCULATE INVERSE OF ELEMENT MASS MATRICES
    %======================================================================
    % --------------------------invM-------------------------------
    if nPdofel==1
        invM_blk = 1./M_blk;
    else
        M_blk = M_blk .* invdetJ(:,ones(1,nPdofel*(nPdofel+1)/2));
        % How to calculate the determinante of several symmetric 3x3 matrices:
        % det(M) =   M11*M22*M33 + M12*M23*M31 + M13*M21*M32
        %          - M13*M22*M31 - M12*M21*M33 - M11*M23*M32
        %   written in 1-index notation:
        %        =   M1 *M5 *M9  + M4 *M8 *M3  + M7 *M2 *M6
        %          - M7 *M5 *M3  - M4 *M2 *M9  - M1 *M8 *M6
        %   considering symmetry and using lower triangular part only
        %   (M4=M2, M7=M3, M8=M6)
        %        =   M1 *M5 *M9  + M2 *M6 *M3  + M3 *M2 *M6
        %          - M3 *M5 *M3  - M2 *M2 *M9  - M1 *M6 *M6
        %   and knowing where the 6 different values are stored in M_blk
        %   (1-->M1, 2-->M2, 3-->M3, 4-->M5, 5-->M6, 6-->M9)
        %        =   M_blk1*M_blk4*M_blk6 + M_blk2*M_blk5*M_blk3 + M_blk3*M_blk2*M_blk5
        %          - M_blk3*M_blk4*M_blk3 - M_blk2*M_blk2*M_blk6 - M_blk1*M_blk5*M_blk5
        %   re-arranging
        %        =   M_blk1 * (M_blk4*M_blk6 - M_blk5*M_blk5)
        %          + M_blk2 * (M_blk5*M_blk3 - M_blk2*M_blk6)
        %          + M_blk3 * (M_blk2*M_blk5 - M_blk4*M_blk3)
        detM_blk = M_blk(:,1) .* (M_blk(:,4).*M_blk(:,6) - M_blk(:,5).*M_blk(:,5)) + ...
                M_blk(:,2) .* (M_blk(:,5).*M_blk(:,3) - M_blk(:,2).*M_blk(:,6)) + ...
                M_blk(:,3) .* (M_blk(:,2).*M_blk(:,5) - M_blk(:,4).*M_blk(:,3));
        detM_blk = detJ .* detM_blk;
        
        % The determinante is used to calculate the inverse of the symmetric
        % 3x3 element mass matrices. The same logic as above is used.
        invM_blk(:,1) = (M_blk(:,4).*M_blk(:,6) - M_blk(:,5).*M_blk(:,5))./detM_blk;
        invM_blk(:,2) = (M_blk(:,5).*M_blk(:,3) - M_blk(:,2).*M_blk(:,6))./detM_blk;
        invM_blk(:,3) = (M_blk(:,2).*M_blk(:,5) - M_blk(:,4).*M_blk(:,3))./detM_blk;
        invM_blk(:,4) = invM_blk(:,2);
        invM_blk(:,5) = (M_blk(:,1).*M_blk(:,6) - M_blk(:,3).*M_blk(:,3))./detM_blk;
        invM_blk(:,6) = (M_blk(:,2).*M_blk(:,3) - M_blk(:,1).*M_blk(:,5))./detM_blk;
        invM_blk(:,7) = invM_blk(:,3);
        invM_blk(:,8) = invM_blk(:,6);
        invM_blk(:,9) = (M_blk(:,1).*M_blk(:,4) - M_blk(:,5).*M_blk(:,5))./detM_blk;
    end
    
    % --------------------------invM*G'----------------------------
    for i=1:nPdofel
        for j=1:nUdofel
            for k=1:nPdofel
                invMG_blk(:,(i-1)*nUdofel+j) = invMG_blk(:,(i-1)*nUdofel+j) ...
                    + invM_blk(:,(i-1)*nPdofel+k).*G_blk(:,(k-1)*nUdofel+j);
            end
        end
    end
    
    % -----------------Kp = Kp + penalty*G*invM*G'-------------------
    indx  = 1;
    for i=1:nUdofel
        for j=i:nUdofel
            for k=1:nPdofel
                K_blk(:,indx) = K_blk(:,indx) ...
                    + penalty*G_blk(:,(k-1)*nUdofel+i).*invMG_blk(:,(k-1)*nUdofel+j);
            end
            indx = indx + 1;
        end
    end
    
    %==============================================================
    % STORE DATA OF ALL ELEMENTS IN BLOCK FOR ASSEMBLY
    %==============================================================
    K_all(:,il:iu) = K_blk';
    if itStrainRate==1
        G_all   (:,il:iu) = G_blk';
        invM_all(:,il:iu) = invM_blk';
        Rhs_all (:,il:iu) = Rhs_blk';
    end
    
    %==============================================================
    % READJUST START, END AND SIZE OF NEXT BLOCK
    %==============================================================
    il = il + nelblk;
    if iblk==nblk-1
        % Account for different number of elements in last block
        nelblk = nel-iu; % number of remaining elements
    end
    iu = iu  + nelblk;
  end % END OF ELEMENT BLOCK LOOP

    %==========================================================================
    % ASSEMBLY OF GLOBAL SPARSE MATRICES AND RHS-VECTOR
    %==========================================================================
    % global stiffness matrix
    OPTS_MUTILS.symmetric  = 1;
    OPTS_MUTILS.n_node_dof = 2;
    OPTS_MUTILS.nthreads   = 1; % >1 is very slow !
    KK = sparse_create(EL2NOD,K_all,OPTS_MUTILS); clear K_all

    if itStrainRate==1
        % global inverse mass matrix
        OPTS_MUTILS.symmetric  = 0;
        OPTS_MUTILS.n_node_dof = 1;
        OPTS_MUTILS.nthreads   = 1; % >1 is very slow !
        invMM = sparse_create(EL2NODP,invM_all,OPTS_MUTILS); clear invM_all

        % global gradient matrix
        EL2DOF                   = zeros(nUdofel,nel,'int32');
        EL2DOF(1:ndim:nUdofel,:) = ndim*(EL2NOD-1)+1;
        EL2DOF(2:ndim:nUdofel,:) = ndim*(EL2NOD-1)+2;

        G_i    = repmat(EL2DOF,nPdofel,1);
        indx_j = repmat(1:nPdofel,nUdofel,1);
        G_j    = EL2NODP(indx_j,:);
        GG     = sparse3(G_i(:) , G_j(:) , G_all(:) , nU , nP); clear G_all

        % global r.h.s.-vector
        Rhs = accumarray(double(EL2DOF(:)),Rhs_all(:));
    end

end % END OF LOCAL FUNCTION matrix_assembly

% #########################################################################

function VAR = rotate_stresses(VAR, dt, jaumann)
    % rotation of elastic stresses for moving Lagrangian points
    % JGP: augmented for more precise rotations, further than Jaumann.
    % Precise rotations set as default
    if nargin < 3
        jaumann = true;
    end
    phi = Var.Vort_xz * dt;                                                         % [nel, nip_stress] rotation angle 
    
    if jaumann		                                                          % Jaumann' simplified co-rotation formulas
        VAR.Tau_xx_old = VAR.Tau_xx + 2 * phi .* VAR.Tau_xz;                          % [nel, nip_stress] as Vort_xz sign is contrary to that in Gerya, actually this matches Gerya, his Eq. (12.24)
        VAR.Tau_zz_old = VAR.Tau_zz - 2 * phi .* VAR.Tau_xz;                          % [nel, nip_stress]                      "                   "                                     Eq. (12.25)
        VAR.Tau_xz_old = VAR.Tau_xz - (VAR.Tau_xx - VAR.Tau_zz) .* phi;               % [nel, nip_stress]                      "                   "                                     Eq. (12.26) 
    else % added by JGP
        VAR.Tau_xx_old = VAR.Tau_xx .* cos(phi).^2 + VAR.Tau_zz .* sin(phi).^2 + VAR.Tau_xz .* sin(2*phi); % sin(2*phi) == 2*sin(phi)*cos(phi)
        VAR.Tau_zz_old = VAR.Tau_xx .* sin(phi).^2 + VAR.Tau_zz .* cos(phi).^2 - VAR.Tau_xz .* sin(2*phi);
        VAR.Tau_xz_old = (VAR.Tau_zz - VAR.Tau_xx) .* sin(phi).*cos(phi) + VAR.Tau_xz .* cos(2*phi);       % cos(2*phi) == cos(phi)^2 - sin(phi)^2
    end
end % END OF LOCALFUNCTION rotate_stresses
