function [Vel, Pressure, PRES_IP, TAU_xx, TAU_yy, TAU_xy, ...
    STRAIN_xx, STRAIN_yy, STRAIN_xy, ...
    Gamma, Yield_T2, YC, E2all, Mu_all, Mu_dis_all, Mu_dif_all, Mu_b_all, Dx_e, ...
    F_xx, F_xy, F_yx, F_yy, I2, GIP_x_all, GIP_y_all, ...
    THETA_all, W_xy, nw_it, Residual] = mechanical2d_m(ELEM2NODE, Phases, GCOORD, Point_id, Temp, E2all, ...
    Mu_all, RHEOL, RHEOLvar, Phi0, Cohesion0, R, Shearm, Rho, G, Bc_ind, ...
    Bc_val, nip, reorder, ext_erate, top_surface, fs_alpha, fs_beta, dt, ...
    Bc_ind_fs, Bc_val_fs, bc_t, SS, PHY, F_xx, F_xy, F_yx, F_yy, I2, THETA_all, ...
    TAU_xx_old, TAU_yy_old, TAU_xy_old, is_elastic, nw_it, ...
    F_ext, rho_ext, Load_el, PLASTICITY, SCALE, PLOT, SOLVER)

% MECHANICAL2D Two dimensional finite element mechanical problem solver of MILAMIN
%
% INPUT [I]
% ELEM2NODE  :: [nnodel,nel]
% Phases     :: [1,nel]
% CGOORD     :: [2,nnod7]
% Point_id   :: [1,nnod6]
% Temp       :: [1,nnod7]
% E2_all     :: [nel,nip], IO strain rate
% Mu_all     :: [nel,nip], IO viscosity
% F_xx       :: [nel,nip], IO historic (or accumulated) change of deformation gradient (Malvern, 1969)
% F_xy       :: [nel,nip], IO   " "
% F_yx       :: [nel,nip], IO   " "
% F_yy       :: [nel,nip], IO   " "
% I2         :: IO, struct ['c','f','p'], with all keys in [nel,nip]
% Rho        :: [nel,nip], I density 
% THETA_all  :: [nel,nip], IO
% nw_it      :: IO, vector augmented at each timestep with the number of solver iterations
%
% OUTPUT [O]
% DISPL      :: [2,nnod]
% Pressure   :: [np,nel], where np=3
% PRES_IP    :: [nel,nip]
% TAU_xx     :: [nel,nip]
% TAU_yy     :: [nel,nip]
% TAU_xy     :: [nel,nip]
% STRAIN_xx  :: [nel,nip]
% STRAIN_xy  :: [nel,nip]
% STRAIN_yy  :: [nel,nip]
% Gamma      :: [nel,nip]
% Yield_T2   :: [nel,nip]
% YC         :: [nel,nip] LOGICAL, yield criterium
% E2_all     :: [nel,nip], IO
% Mu_all     :: [nel,nip], IO
% Mu_b_all   :: [nel,nip]
% Mu_dif_all :: [nel,nip]
% Mu_dis_all :: [nel,nip]
% Dx_e       :: diagnostic of surface node vertical displacement for courant condition
% F_xx       :: [nel,nip], IO
% F_xy       :: [nel,nip], IO
% F_yx       :: [nel,nip], IO
% F_yy       :: [nel,nip], IO
% I2         :: IO, struct ['c','f','p'], with all keys in [nel,nip]
% GIP_x_all  :: [nel,nip]
% GIP_y_all  :: [nel,nip]
% THETA_all  :: [nel,nip], IO
% W_xy       :: [nel,nip] REAL half vorticity [i.e. angular velocity]. Positive for counterclockwise rotation
% nw_it      :: IO, vector augmented at each timestep with the number of solver iterations
% Residual   :: [1,niter], giving residual at each solver iteration

%   Part of MILAMIN: MATLAB-based FEM solver for large problems, Version 1.0
%   Copyright (C) 2007, M. Dabrowski, M. Krotkiewski, D.W. Schmid
%   University of Oslo, Physics of Geological Processes
%   http://milamin.org
%   See License file for terms of use.

%   Edited by Miguel Andres-Martinez, Elena Ros-Bernabeu and Marta 
%   Perez-Gussinye. Earth Sciences Department. This version includes:
%       - Free surface algorithm
%       - Fixed diagonal terms of the Stiffnes matrix (the value was
%           previously duplicated by the code.
%       - Strain softening
%       - Elasticity
%       - Varying rheologies along a phase
%       - Plasticity based on Moresi, 2003
%       - Yield based on discontinues pressures
%       - Density dependent on temperature
%       - Dislocation and diffusion softening controlled by different
%           parameters
%       - Viscous strain softening dependent on Arrhenius equation
%           (temperature dependent)
%       - Random weak noise
%       - Sea level pressures

% TODO Linear shape functions for stress, strain and viscosities
% TODO Continuous pressure
% TODO Put viscosities in a structure and change plots
% TODO Integration point coordinate calculation into a function and remove
%   GIP from the output
% TODO scaling of the strain softening? (coded in 
%   mechanical2d_cohesion_co.m)
% TODO yielding criterium before strain softening in phi and cohesion? 
%   (coded in mechanical2d_cohesion.m)

% 05/01/2016 ER
    % Activation volume added to the flow law
% 01/05/2016 DD & MA
    % Fixed error in which the coefficients for transforming rheological
    % parameters from uniaxial/triaxial to second square invariants was
    % accumulated along the rheol-var loop
% 08/06/2016 MA
    % Strain softening in an independent function to speed up the code
% 29/11/2016 MA
    % Distributed water load taken out of the function for generalization
    % purposes
% 03/11/2020 JGP: Note: the solution provides approximate but not necessarily equal module of the normal deviatoric stresses
%                       this is derived from the fact that it is not guaranteed STRAIN_xx + STRAIN_yy = 0., and deviatoric stress calculations at L1181... 
%==========================================================================
% MODEL INFO
%==========================================================================
nnod        = size(GCOORD,2);
nel         = size(ELEM2NODE,2);

%==========================================================================
% CONSTANTS
%==========================================================================
ndim        = 2;
nnodel      = 7;
nedof       = nnodel*ndim;
sdof        = 2*nnod;
np          = 3;

%DEV   = [ 4/3 -2/3 0;...
%         -2/3  4/3 0;...
%           0    0  1];

mu_max      = SOLVER.mu_max;
mu_min      = SOLVER.mu_min;
PF          = SOLVER.fpen;

C1 = 4/3.;
C2 = 2/3.;

%==========================================================================
% BLOCKING PARAMETERS (nelblo must be < nel)
%==========================================================================
nelblo          = min(nel, SOLVER.nelblo);
nblo            = ceil(nel/nelblo);

%==========================================================================
% i) PREPARE INTEGRATION POINTS & DERIVATIVES wrt LOCAL COORDINATES
%==========================================================================
[IP_X, IP_w]    = ip_triangle(nip);
[   N, dNdu]    = shp_deriv_triangle(IP_X, nnodel);                       % x_ip [nip,2], IP_w [1,nip] with nip=6; 

%==========================================================================
% DECLARE VARIABLES (ALLOCATE MEMORY)
%==========================================================================
% A_all       = zeros(nedof*(nedof+1)/2,nel);
% Q_all       = zeros(nedof*np,nel);
% invM_all    = zeros(np*np,nel);
% Rhs_all     = zeros(nedof,nel);
% Tauxx       = zeros(nel,nip);
Vel         = zeros(sdof,1);                     % O
% ELEM_DOF = zeros(nedof, nel,'int32');
% ELEM_DOF(1:ndim:end,:) = ndim*(ELEM2NODE-1)+1;
% ELEM_DOF(2:ndim:end,:) = ndim*(ELEM2NODE-1)+2;

GIP_x_all   = zeros(nel,nip);                    % O      
GIP_y_all   = zeros(nel,nip);                    % O
STRAIN_xx   = zeros(nel,nip);                    % O
STRAIN_yy   = zeros(nel,nip);                    % O
STRAIN_xy   = zeros(nel,nip);                    % O
TAU_xx      = zeros(nel,nip);                    % O
TAU_yy      = zeros(nel,nip);                    % O
TAU_xy      = zeros(nel,nip);                    % O
E2it        = ones(nel,nip) * ext_erate;         % L
Mu_b_all    = zeros(nel,nip);                    % O 
Mu_dif_all  = zeros(nel,nip);                    % O
Mu_dis_all  = zeros(nel,nip);                    % O

Mu_c_all    = zeros(nel,nip);                    % L
%MAX_TAU     = zeros(nel,nip);                   % currently not used
PRES_IP     = zeros(nel,nip);                    % O
%YIELD       = zeros(nel,nip);                    % currently not used
F0_xx      = F_xx;% 
F0_xy      = F_xy;
F0_yx      = F_yx;
F0_yy      = F_yy;
W_xy       = zeros(nel,nip);                     % O, Rotation
YC         = repelem(false,nel,nip);             % O, LOGICAL
YC_cum     = YC;
Gamma      = zeros(nel,nip);                     % O
Yield_T2   = zeros(nel,nip);                     % O
I2p0       = I2.p;                               % L, [nel,nip]
I2c0       = I2.c;                               % L, [nel,nip]
Residual   = Inf;                                % O

% Calculate parameters for strain softening
[Phi_all,Cohesion_all,Pef_dif_all,Pef_dis_all] = strain_softening_SST_rand ...
    (nip,nel,Phi0,Cohesion0,I2,PRES_IP,SS, PHY, Temp, ELEM2NODE, Phases,GCOORD,Point_id);

% Calculate lithostatic pressures if needed
switch PLASTICITY.yc
    case 'vonmissesP'
        [xip,zip] = ip_coord(GCOORD,ELEM2NODE,nel,nip);
        PLITHOST = plithost(GCOORD,ELEM2NODE,Point_id,Phases,Rho, ...
            G,xip,zip);
end

% Plastic phases
P_Phases = ismember(Phases,PLASTICITY.Phases);

% Set up plotting function
plot_handle = str2func(PLOT.func);
%==========================================================================
% INDICES EXTRACTING LOWER PART
%==========================================================================
indx_l = tril(ones(nedof)); indx_l = indx_l(:); indx_l = indx_l==1;
ei = 0;
track_pli = 0;

while true %strain rate iteration
    ei = ei + 1;
    
    %======================================================================
    % FLOW SOLUTION OPTIMIZED VERSION
    %======================================================================
    %==================================================================
    % DECLARE VARIABLES (ALLOCATE MEMORY)
    %==================================================================
    nelblo      = min(nel, SOLVER.nelblo);
    A_block     = zeros(nelblo, nedof*(nedof+1)/2);
    Q_block     = zeros(nelblo, np*nedof);
    M_block     = zeros(nelblo, np*(np+1)/2);
    invM_block  = zeros(nelblo, np*np);
    invMQ_block = zeros(nelblo, np*nedof);
    Pi_block    = zeros(nelblo, np);
    Rhs_block   = zeros(nelblo, nedof);

    A_all       = zeros(nedof*(nedof+1)/2,nel);
    Q_all       = zeros(nedof*np, nel);
    invM_all    = zeros(np*np, nel);
    Rhs_all     = zeros(nedof, nel);
        
    il          = 1;
    iu          = nelblo;
    %==================================================================
    % i) BLOCK LOOP - MATRIX COMPUTATION
    %==================================================================

    fprintf(1, 'MATRIX COMPUTATION: '); tic;
    for ib = 1:nblo
        els_blk = il:iu;

        invJx   = zeros(nelblo, ndim);
        invJy   = zeros(nelblo, ndim);
        
        %==============================================================
        % ii) FETCH DATA OF ELEMENTS IN BLOCK
        %==============================================================
        ECOORD_x    = reshape( GCOORD(1,ELEM2NODE(:,els_blk)), nnodel, nelblo);
        ECOORD_y    = reshape( GCOORD(2,ELEM2NODE(:,els_blk)), nnodel, nelblo);
        Temp_bl     = reshape( Temp(ELEM2NODE(:,els_blk)), nnodel, nelblo);
        
        Adis_block  = RHEOL.Adis(Phases(il:iu),:);
        Ndis_block  = RHEOL.Ndis(Phases(il:iu),:);
        Qdis_block  = RHEOL.Qdis(Phases(il:iu),:);
        Vdis_block  = RHEOL.Vdis(Phases(il:iu),:);
        
        Adif_block  = RHEOL.Adif(Phases(il:iu),:);
        Ndif_block  = RHEOL.Ndif(Phases(il:iu),:);
        Qdif_block  = RHEOL.Qdif(Phases(il:iu),:);
        Vdif_block  = RHEOL.Vdif(Phases(il:iu),:);
        
        Shear_block = Shearm(Phases(il:iu));                                   % [nelblo,1]
        
        % Plastic phases of the block
        P_Phases_block = P_Phases(il:iu);
        
        Var_block   = zeros(nelblo,size(RHEOLvar,2));
        
        %ERho        = Rho(Phases(il:iu));
 
        %==============================================================
        % iii) INTEGRATION POINT LOOP
        %==============================================================
        A_block(:)      = 0;
        Q_block(:)      = 0;
        M_block(:)      = 0;
        invM_block(:)   = 0;
        invMQ_block(:)  = 0;
        Rhs_block(:)    = 0;

        a23   = ECOORD_x(2,:).*ECOORD_y(3,:) - ECOORD_x(3,:).*ECOORD_y(2,:);
        a31   = ECOORD_x(3,:).*ECOORD_y(1,:) - ECOORD_x(1,:).*ECOORD_y(3,:);
        a12   = ECOORD_x(1,:).*ECOORD_y(2,:) - ECOORD_x(2,:).*ECOORD_y(1,:);
        area  = a23 + a31 + a12;

        for ip=1:nip
            
            Ni = N{ip};                                                        % [nnodel, 1]
            %Temp at integration point
            Temp_ip = (Ni'*Temp_bl)';                                          % [nelblo,1]         ([1,nnodel]*[nnodel,nelblk])'
            ED = zeros(nelblo,1);                                              % [nelblo,1]
            Mu_c_block = zeros(nelblo,1);                                      % [nelblo,1]
            
            % Load cohesion, friction angle and viscous softening factors
            Cohesion = Cohesion_all(il:iu,ip);                                 % [nelblo,1]
            Phi = Phi_all(il:iu,ip);                                           % [nelblo,1]
            Pef_dis = Pef_dis_all(il:iu,ip);                                   % [nelblo,1]
            Pef_dif = Pef_dif_all(il:iu,ip);                                   % [nelblo,1]
            for nvar = 1:length(RHEOLvar) 
                Var_block(:,nvar) = RHEOLvar{nvar}(il:iu,ip);
            end
            
            ERho = Rho(il:iu,ip);
            
            %==============================================================
            % VISCOSITY UPDATE
            %==============================================================
            if ei == 1                                                         % viscosity update: first iteration
                if all(Mu_all(il:iu,ip)==0)
                    E2               = mean(E2all(il:iu,:),2);
                    Sc_dis = zeros(length(il:iu),1);
                    for n = 1:size(Ndis_block,2)
                        Sc_dis = ...
                            1./(2.^((Ndis_block(:,n)-1)./Ndis_block(:,n)).* ...
                            3.^((Ndis_block(:,n)+1)./(2*Ndis_block(:,n))));
                        
                        ED = ED + Var_block(:,n).*...                          % [nelblo,1]
                            (Sc_dis.*Adis_block(:,n).^(-1./Ndis_block(:,n)).* ...
                            E2.^(1./Ndis_block(:,n)-1).* ...
                            exp(Qdis_block(:,n)./ ...
                            (Ndis_block(:,n).*R.*(Temp_ip+273))));
                    end
                    if is_elastic
                       ED = (ED .* Shear_block.*dt) ./ (ED + Shear_block.*dt);   % total resistance of // resistors
                    end
                    ED(ED>=mu_max)   = mu_max;
                    ED(ED<=mu_min)   = mu_min;
                    Mu_all(il:iu,ip) = ED;
                    %                     Mu_ve = ED;
                    %                     Mu_ve_all(il:iu,ip) = Mu_ve;
                else
                    ED = Mu_all(il:iu,ip);
                    %                     Mu_ve = Mu_ve_all(il:iu,ip);
                end
            else  % viscosity update: > 1st iteration
                % Load the strain rate of the previous iteration
                E2 = E2all(il:iu,ip);
                
                % Initialize Mu_dis_block and Mu_dif_block
                Mu_dis_block = zeros(length(il:iu),1);
                Mu_dif_block = zeros(length(il:iu),1);
                Sc_dis = zeros(length(il:iu),1);
                Sc_dif = zeros(length(il:iu),1);
                
                % Loop for rheologic variation in the same phase
                for n = 1:size(Ndis_block,2)
                    % Factors for scaling triaxial and uniaxial experiment
                    % parameters (GERYA 2010)
                    Sc_dis = ...
                        1./(2.^((Ndis_block(:,n)-1)./Ndis_block(:,n)).* ...
                        3.^((Ndis_block(:,n)+1)./(2*Ndis_block(:,n))));
                    Sc_dif = 1/3;
                    
                    % Dislocation creep
                    Mu_dis_block = Mu_dis_block + Var_block(:,n).* ...
                        (Sc_dis.*(Pef_dis.*Adis_block(:,n)).^(-1./Ndis_block(:,n)) ...
                        .* E2.^(1./Ndis_block(:,n)-1) ...
                        .* exp((Qdis_block(:,n) + PRES_IP(il:iu,ip).*Vdis_block(:,n)) ...
                        ./(Ndis_block(:,n).*R.*(Temp_ip+273))));
                    % Diffusion creep
                    Mu_dif_block = Mu_dif_block + Var_block(:,n).* ...
                        (Sc_dif.*(Pef_dif.*Adif_block(:,n)).^(-1./Ndif_block(:,n)) ...
                        .* (RHEOL.Grain/RHEOL.Burger).^0 ...
                        .* exp((Qdif_block(:,n)+PRES_IP(il:iu,ip).*Vdif_block(:,n)) ...
                        ./(Ndif_block(:,n).*R.*(Temp_ip+273))));
                end
                                                
                % Effective viscosity
                EDI = 1./Mu_dis_block;                                         % [nelblo,1]
                EDI(Ndif_block(:,1)~=0) = EDI(Ndif_block(:,1)~=0) + 1./Mu_dif_block(Ndif_block(:,1)~=0);
                if is_elastic
                    EDI = EDI + 1./ (Shear_block.*dt);
                end
                % + 1./Mu_brittle(Ndif_block(:,1)==0)
                ED = 1./EDI;
                clear EDI;
                       
                %                 % Visco-elastic
                %                 Mu_ve = ED;

                Mu_c_block(Ndif_block(:,1)~=0) = ...
                    (1./Mu_dis_block(Ndif_block(:,1)~=0) ...
                    + 1./Mu_dif_block(Ndif_block(:,1)~=0)).^-1;
                Mu_c_block(Ndif_block(:,1)==0) = ...
                    (1./Mu_dis_block(Ndif_block(:,1)==0)).^-1;
                
                % Limit viscosities
                ED(ED <= mu_min) = mu_min;               
                ED(ED >  mu_max) = mu_max;
                Mu_c_block(Mu_c_block<=mu_min) = mu_min;
                Mu_c_block(Mu_c_block>mu_max)  = mu_max;

                %==========================================================
                % PLASTICITY
                %==========================================================
                % Different yield [limit elastic behaviour] criterions (formulation from Spiegelman, 2016)
                switch PLASTICITY.yc
                    case 'druckerprager'                  % default
                        Ayc = Cohesion.*cos(Phi);
                        Byc = sin(Phi);
                        P_IP = PRES_IP(il:iu,ip);         % pressure at this ip for all elements in the block
                    case 'vonmisses'
                        Ayc = Cohesion;
                        Byc = 0;
                        P_IP = zeros(nelblo,1);
                    case 'vonmissesP'
                        Ayc = Cohesion.*cos(Phi);
                        Byc = sin(Phi);
                        P_IP = PLITHOST(il:iu,ip);
                end % switch PLASTICITY.yc
                T2_yield = Ayc + Byc.*P_IP;
                
                Smaller_c = T2_yield < Ayc;
                %                 if sum(sum(Smaller_c))  % TODO check why this happens
                %                     warning("COHESION LARGER THAN YIELD!")
                %                 end
                T2_yield(Smaller_c) = Ayc(Smaller_c);
                                
                switch PLASTICITY.type
                    case {'moresi','moresi_cum'} % default:"moresi"
                        %ED_old = ED;
                        T2old = sqrt(0.5*(TAU_xx_old(il:iu,ip).^2 + TAU_yy_old(il:iu,ip).^2) + TAU_xy_old(il:iu,ip).^2);              % IInd invariant of deviatoric stresses
                                                                                                                                      % currently not used
                        % Forecasting new stresses with previous iteration strain
                        % rates and current viscosity
                        Txx_block_f = 4/3.*ED.*STRAIN_xx(il:iu,ip) ...
                            -2/3.*ED.*STRAIN_yy(il:iu,ip) + ...
                            2.*ED.*THETA_all(il:iu,ip).*TAU_xx_old(il:iu,ip);
                        Tyy_block_f = -2/3.*ED.*STRAIN_xx(il:iu,ip) ...
                            + 4/3.*ED.*STRAIN_yy(il:iu,ip) ...
                            + 2.*ED.*THETA_all(il:iu,ip).*TAU_yy_old(il:iu,ip);
                        Txy_block_f = 2.*ED.*STRAIN_xy(il:iu,ip) + ...
                            2.*ED.*THETA_all(il:iu,ip).*TAU_xy_old(il:iu,ip);
                        T2forecast = sqrt(0.5*(Txx_block_f.^2 + Tyy_block_f.^2) ...                                                   % IInd invariant of deviatoric stress forecast
                            + Txy_block_f.^2);
                        
                        % Yield criterium
                        % No memory for plasticity:
                        YC(il:iu,ip) = T2forecast > T2_yield;
                        % Memory of plasticity (uncomment):
                        %YC(il:iu,ip) = T2forecast>T2_yield | YC(il:iu,ip);
                        YC_block = YC(il:iu,ip);
                        
                        % Ductile viscosity
                        Yield_b           = T2_yield;
                        
                        % Calculate E2eff as in Eq. 35, Moresi, 2003
                        E_eff_xx = 2.*STRAIN_xx(il:iu,ip) ...
                            + 1./(Shear_block.*dt).*TAU_xx_old(il:iu,ip);
                        E_eff_yy = 2.*STRAIN_yy(il:iu,ip) ...
                            + 1./(Shear_block.*dt).*TAU_yy_old(il:iu,ip);
                        E_eff_xy = 2.*STRAIN_xy(il:iu,ip) ...
                            + 1./(Shear_block.*dt).*TAU_xy_old(il:iu,ip);
                        E2eff = sqrt(0.5*(E_eff_xx.^2 + E_eff_yy.^2) + E_eff_xy.^2);
                        
                        Mu_brittle = Yield_b ./ E2eff;
                        
                        Gamma(il:iu,ip) =  Yield_b.*(1./Mu_brittle - 1./ED);
                        Yield_T2(il:iu,ip)      = Yield_b;
                        
                        if PLASTICITY.type == "moresi_cum" && ei>10
                            YC_block = YC_block | YC_cum(il:iu,ip);
                            YC_cum(il:iu,ip) = YC_block;
                        end
                        
                        %ED(YC_block~=0 & T2old~=0) = Mu_brittle(YC_block~=0 & T2old~=0);
                        ED(YC_block & P_Phases_block') = ...
                            Mu_brittle(YC_block & P_Phases_block');
                        
                    case 'moresi_no_fc'
                        % Yield criterium
                        % No memory for plasticity:
                        YC(il:iu,ip) = T2all(il:iu,ip) > T2_yield; %  LOGICAL
                        YC_block = YC(il:iu,ip);                   % LOGICAL
                        
                        % Ductile viscosity
                        Yield_b           = T2_yield;
                        
                        % Calculate E2eff as in Eq. 35, Moresi, 2003
                        E_eff_xx = 2.*STRAIN_xx(il:iu,ip) ...
                            + 1./(Shear_block.*dt).*TAU_xx_old(il:iu,ip);
                        E_eff_yy = 2.*STRAIN_yy(il:iu,ip) ...
                            + 1./(Shear_block.*dt).*TAU_yy_old(il:iu,ip);
                        E_eff_xy = 2.*STRAIN_xy(il:iu,ip) ...
                            + 1./(Shear_block.*dt).*TAU_xy_old(il:iu,ip);
                        E2eff = sqrt(0.5*(E_eff_xx.^2+E_eff_yy.^2)+E_eff_xy.^2);
                        
                        Mu_brittle = Yield_b./E2eff;
                        
                        ED(YC_block~=0 & P_Phases_block') = ...
                            Mu_brittle(YC_block & P_Phases_block');
                    case 'maxwell'
                        Mu_brittle = T2_yield./(2.*E2);
                        ED(P_Phases_block') = ED(P_Phases_block').* ...
                            Mu_brittle(P_Phases_block')./ ...
                            (ED(P_Phases_block')+ ...
                            Mu_brittle(P_Phases_block'));
                    case 'bullshit'
                        Mu_brittle = T2_yield./(2.*E2);
                        ED(P_Phases_block') = ED(P_Phases_block').* ...
                            Mu_brittle(P_Phases_block')./ ...
                            (ED(P_Phases_block')+ ...
                            Mu_brittle(P_Phases_block'));
                    case 'harmonic'
                        Mu_brittle = T2_yield./(2.*E2);
                        % Number of viscosities averaged taking checking if
                        % dislocation, diffusion viscosities and elasticity
                        % term is in place. The + 1 is for the brittle
                        % effective viscosity which is always included in
                        % this case
                        n_harm = (Ndif_block(:,1)~=0) + ...
                                 (Ndis_block(:,1)~=0) + 1.;
                        if is_elastic
                            n_harm = n_harm + ones(nelblo,1);
                        end
                        ED(P_Phases_block') = n_harm(P_Phases_block').* ...
                            ED(P_Phases_block').* ...
                            Mu_brittle(P_Phases_block')./ ...
                            (ED(P_Phases_block')+ ...
                            Mu_brittle(P_Phases_block'));
                    case 'no'
                        Mu_brittle = 0;
                end % switch plasticity
                % Limit viscosities
                ED(ED<=mu_min)   = mu_min;
                ED(ED>mu_max)   = mu_max;
                
                % Save viscosities
                Mu_all(il:iu,ip) = ED;
                %                 Mu_ve_all(il:iu,ip) = Mu_ve;
                Mu_b_all(il:iu,ip) = Mu_brittle;
                Mu_dis_all(il:iu,ip) = Mu_dis_block;
                Mu_dif_all(il:iu,ip) = Mu_dif_block;
                Mu_c_all(il:iu,ip) = Mu_c_block;
            end % viscosity update
            
            % Calculates elasticity factors
            if is_elastic
                Etheta = 1. ./ (2. * Shear_block * dt);                         % [nelblo,1]
            else
                Etheta = repmat(0., nelblo, 1);
            end
            THETA_all(il:iu,ip) = Etheta;
            
            %==========================================================
            % iv) SHAPE FUNCTIONS DERIVATIVES FOR INTEGRATION POINT
            %==========================================================
            dNdui = dNdu{ip};
            GIP_x = Ni'*ECOORD_x;
            GIP_y = Ni'*ECOORD_y;

            tmp   = ECOORD_x(3,:).*GIP_y - GIP_x.*ECOORD_y(3,:);
            eta1  = a23 + tmp + ECOORD_y(2,:).*GIP_x - GIP_y.*ECOORD_x(2,:);
            eta2  = a31 - tmp + ECOORD_x(1,:).*GIP_y - GIP_x.*ECOORD_y(1,:);

            Pi_block(:,1) = eta1./area;
            Pi_block(:,2) = eta2./area;
            Pi_block(:,3) = 1 - Pi_block(:,1) - Pi_block(:,2);

            %==========================================================
            % v) CALCULATE JACOBIAN, ITS DETERMINANT AND INVERSE
            %==========================================================
            Jx          = ECOORD_x'*dNdui;
            Jy          = ECOORD_y'*dNdui;
            detJ        = Jx(:,1).*Jy(:,2) - Jx(:,2).*Jy(:,1);

            invdetJ     = 1.0./detJ;
            invJx(:,1)  = +Jy(:,2).*invdetJ;
            invJx(:,2)  = -Jy(:,1).*invdetJ;
            invJy(:,1)  = -Jx(:,2).*invdetJ;
            invJy(:,2)  = +Jx(:,1).*invdetJ;

            %==========================================================
            % vi) DERIVATIVES wrt GLOBAL COORDINATES
            %==========================================================
            dNdx        = invJx*dNdui';
            dNdy        = invJy*dNdui';

            %==========================================================
            % vii) NUMERICAL INTEGRATION OF ELEMENT MATRICES
            %==========================================================
            weight      = IP_w(ip)*detJ;
            weightD     =    weight.*ED;
            % ------------------------A matrix-------------------------
            indx  = 1;
            for i = 1:nnodel
                % x-velocity equation
                for j = i:nnodel
                    A_block(:,indx) = A_block(:,indx) + ( C1.*dNdx(:,i).*dNdx(:,j) + dNdy(:,i).*dNdy(:,j)).*weightD;
                    indx = indx+1;
                    A_block(:,indx) = A_block(:,indx) + (-C2.*dNdx(:,i).*dNdy(:,j) + dNdy(:,i).*dNdx(:,j)).*weightD;
                    indx = indx+1;
                end
                % y-velocity equation
                for j = i:nnodel
                    if(j>i)
                        A_block(:,indx) = A_block(:,indx) + (-C2.*dNdy(:,i).*dNdx(:,j) + dNdx(:,i).*dNdy(:,j)).*weightD;
                        indx = indx+1;
                    end
                    A_block(:,indx) = A_block(:,indx) + ( C1.*dNdy(:,i).*dNdy(:,j) + dNdx(:,i).*dNdx(:,j)).*weightD;
                    indx = indx+1;
                end
            end

            % ------------------------Q matrix-------------------------
            for i=1:np
                TMP1 = weight.*Pi_block(:,i);
                TMP2 = TMP1(:,ones(1,nnodel));
                Q_block(:,(i-1)*nedof + (1:2:nedof)) =  Q_block(:,(i-1)*nedof + (1:2:nedof)) - TMP2.*dNdx;
                Q_block(:,(i-1)*nedof + (2:2:nedof)) =  Q_block(:,(i-1)*nedof + (2:2:nedof)) - TMP2.*dNdy;
            end

            % ------------------------M matrix-------------------------
            indx = 1;
            for i = 1:np
                for j = i:np
                    M_block(:,indx) = M_block(:,indx) + weight.*Pi_block(:,i).*Pi_block(:,j);
                    indx = indx + 1;
                end
            end

            % -----------------------Rhs vector------------------------
            Rhs_block(:,1:2:nedof) = Rhs_block(:,1:2:nedof) + G(1)*(ERho.*weight)*Ni' ...    % [nelblo,nnodel]
                - dNdx.*repmat(2.*ED.*Etheta.*TAU_xx_old(il:iu,ip).*weight,[1,nnodel])...
                - dNdy.*repmat(2.*ED.*Etheta.*TAU_xy_old(il:iu,ip).*weight,[1,nnodel]);
            Rhs_block(:,2:2:nedof) = Rhs_block(:,2:2:nedof) + G(2)*(ERho.*weight)*Ni' ...
                - dNdy.*repmat(2.*ED.*Etheta.*TAU_yy_old(il:iu,ip).*weight,[1,nnodel])...
                - dNdx.*repmat(2.*ED.*Etheta.*TAU_xy_old(il:iu,ip).*weight,[1,nnodel]);

        end % integration point loop

        %==============================================================
        % viii) STATIC CONDENSATION
        %==============================================================

        % --------------------------invM-------------------------------
        TMP     = 1./area';
        M_block = M_block.*TMP(:,ones(1,np*(np+1)/2));

        detM_block = M_block(:,1).*(M_block(:,4).*M_block(:,6) - M_block(:,5).*M_block(:,5)) + ...
            M_block(:,2).*(M_block(:,5).*M_block(:,3) - M_block(:,2).*M_block(:,6)) + ...
            M_block(:,3).*(M_block(:,2).*M_block(:,5) - M_block(:,4).*M_block(:,3));

        detM_block = detM_block./TMP;
        invM_block(:,1) = (M_block(:,4).*M_block(:,6) - M_block(:,5).*M_block(:,5))./detM_block;
        invM_block(:,2) = (M_block(:,5).*M_block(:,3) - M_block(:,2).*M_block(:,6))./detM_block;
        invM_block(:,3) = (M_block(:,2).*M_block(:,5) - M_block(:,4).*M_block(:,3))./detM_block;
        invM_block(:,4) = invM_block(:,2);
        invM_block(:,5) = (M_block(:,1).*M_block(:,6) - M_block(:,3).*M_block(:,3))./detM_block;
        invM_block(:,6) = (M_block(:,2).*M_block(:,3) - M_block(:,1).*M_block(:,5))./detM_block;
        invM_block(:,7) = invM_block(:,3);
        invM_block(:,8) = invM_block(:,6);
        invM_block(:,9) = (M_block(:,1).*M_block(:,4) - M_block(:,5).*M_block(:,5))./detM_block;

        % --------------------------invM*Q'----------------------------
        for i=1:np
            for j=1:nedof
                for k=1:np
                    invMQ_block(:,(i-1)*nedof+j) = invMQ_block(:,(i-1)*nedof+j) + invM_block(:,(i-1)*np+k).*Q_block(:,(k-1)*nedof+j);
                end
            end
        end

        % -------------------A = A + PF*Q'*invM*Q'---------------------
        indx = 1;
        for i=1:nedof
            for j=i:nedof
                for k=1:np
                    A_block(:,indx) = A_block(:,indx) + PF*Q_block(:,(k-1)*nedof+i).*invMQ_block(:,(k-1)*nedof+j);
                end
                indx = indx + 1;
            end
        end

        %==============================================================
        % ix) WRITE DATA INTO GLOBAL STORAGE
        %==============================================================
        A_all(:,il:iu)  = A_block';
        Q_all(:,il:iu)  = Q_block';
        invM_all(:,il:iu)  = invM_block';
        Rhs_all(:,il:iu)  = Rhs_block';

        %==============================================================
        % READJUST START, END AND SIZE OF NEXT BLOCK. REALLOCATE MEMORY
        %==============================================================
        il = iu + 1;
        if ib == nblo-1
            nelblo 	= nel - iu;
            A_block = zeros(nelblo, nedof*(nedof+1)/2);
            Q_block = zeros(nelblo, np*nedof);
            M_block = zeros(nelblo, np*(np+1)/2);
            invM_block = zeros(nelblo, np*np);
            invMQ_block = zeros(nelblo, np*nedof);
            Pi_block  = zeros(nelblo, np);
            Rhs_block = zeros(nelblo, nedof);
        end
        iu = iu + nelblo;
    end  % for ib = 1:nblk
    fprintf(1, [num2str(toc),'\n']);

    %======================================================================
    % ix) CREATE TRIPLET FORMAT INDICES
    %======================================================================
    tic; fprintf(1, 'TRIPLET INDICES:    ');
    %A matrix
    ELEM_DOF = zeros(nedof, nel,'int32');
    ELEM_DOF(1:ndim:end,:) = ndim*(ELEM2NODE-1)+1;                           % x
    ELEM_DOF(2:ndim:end,:) = ndim*(ELEM2NODE-1)+2;                           % y
    indx_j = repmat(1:nedof,nedof,1); indx_i = indx_j';
    indx_i = tril(indx_i); indx_i = indx_i(:); indx_i = indx_i(indx_i>0);
    indx_j = tril(indx_j); indx_j = indx_j(:); indx_j = indx_j(indx_j>0);

    A_i = ELEM_DOF(indx_i(:),:);
    A_j = ELEM_DOF(indx_j(:),:);

    indx       = A_i < A_j;
    tmp        = A_j(indx);
    A_j(indx)  = A_i(indx);
    A_i(indx)  = tmp;

    %Q matrix
    Q_i = repmat(int32(1:nel*np),nedof,1);
    Q_j = repmat(ELEM_DOF,np,1);

    %invM matrix
    indx_j = repmat(1:np,np,1); indx_i = indx_j';
    invM_i = reshape(int32(1:nel*np),np, nel);
    invM_j = invM_i(indx_i,:);
    invM_i = invM_i(indx_j,:);

    fprintf(1, [num2str(toc),'\n']);

    %======================================================================
    % x) CONVERT TRIPLET DATA TO SPARSE MATRIX
    %======================================================================
    %clear block variables
    clear 'Q_block' 'A_block' 'invMQ_block' 'Rhs_block'
    fprintf(1, 'SPARSIFICATION:     '); tic
    A    = sparse2(A_i(:)   ,    A_j(:),    A_all(:));
    Q    = sparse2(Q_i(:)   ,    Q_j(:),    Q_all(:));
    invM = sparse2(invM_i(:), invM_j(:), invM_all(:));
    Rhs  = accumarray(ELEM_DOF(:), Rhs_all(:));
    clear A_i A_j A_all Q_i Q_j Q_all invM_i invM_j invM_all Rhs_all;
    fprintf(1, [num2str(toc),'\n']);
    
    %     ind_Asparse_example = find(A~=0); % Uncomment this for plotting the
    %         % elements of the matrix A. Set a break point in 'Break point
    %         % here'.
    %     
    %     [iAsparse,jAsparse] = ind2sub(size(A),ind_Asparse_example);
    %     figure(3)
    %     plot(jAsparse,iAsparse,'.');
    %     xlabel('J')
    %     ylabel('I')
    %     axis ij 
    %     % Break point here
    
    %======================================================================
    % ADD EXTERNAL LOAD
    %======================================================================
    if any(F_ext~=0)
        Rhs = Rhs + F_ext;
    end
    
    %======================================================================
    % CONVERGENCE CRITERIUM PREVIUS STEP (SPIEGELMAN, 2016)
    %======================================================================
    if ei > 1
        if ei == 2
            Residual = [];
        end
        
        % CALCULATE AND SCALE RESIDUAL FOR NON-SCALED MODEL (SI)
        % ------------------------------------------------------
        % Load the scaling factors
        L0 = SCALE.L0;
        V0 = SCALE.V0;
        t0 = L0/V0;
        mu0 = SCALE.mu0;
        
        % Remove penalty-related terms from A matrix
        AR = A - tril(Q'*PF*invM*Q);
        AR_UPPER = cs_transpose(AR);
        AR_UPPER(speye(sdof)==1) = 0;
        AA = AR + AR_UPPER;
        
        % Update right-hand-side Dirichlet boundary contidions
        % using new viscosities
        TMP = AR(:,Bc_ind_fs) + AR_UPPER(:,Bc_ind_fs);
        Rhs_update = Rhs - TMP*Bc_val_fs';
        
        % Calculate residue for Stokes: R_v = Rhs - (A*v + Q'*P)
        ResidualV = Rhs_update(Free) ...
            - (AA(Free,Free)*Vel(Free) + Q(:,Free)'*Pressure(:));                        % note: Pressure(:) as column vector
        % Calculate residue for incompressibility cond.: R_p = Q*v = div(v)
        ResidualP = Q(:,Free)*Vel(Free);
        switch SCALE.t
            case 'residual'
                ResidualV = ResidualV*t0/(mu0*L0);
                ResidualP = ResidualP*t0/(mu0*L0);
            case 'before_back_solve'
                error(['mechanical2d_m has not the scaling option ', ...
                    '"before_back_solve"'])
        end
        
        % 2-limit norm of the scaled residual
        Residual(ei-1) = norm([ResidualV; ResidualP]);
        
        if PLOT.make
            PLOT.sec = "residual";
            plot_handle();
        end
    end % if ei > 1

    %======================================================================
    % SWITCH FOR THE BOUNDARY CONDITIONS & POWELL-HESTENES ITERATIONS - UZAWA-like - TO SOLVE FOR PRESSURE AT VERTICES
    %======================================================================
    switch top_surface
        case 'fix'    
            %==============================================================
            % BOUNDARY CONDITIONS FOR A FIXED TOP SURFACE
            %==============================================================
            fprintf(1, 'BDRY CONDITIONS:    '); tic;
            Free        = 1:sdof;
            Free(Bc_ind)= [];
            A_BC_UPPER = transpose(A);
            A_BC_UPPER = A_BC_UPPER - diag(diag(A_BC_UPPER));
            TMP         = A(:,Bc_ind) + A_BC_UPPER(:,Bc_ind);
            %TMP         = A(:,Bc_ind) + cs_transpose(A(Bc_ind,:));
            Rhs         = Rhs - TMP*Bc_val';
            A           = A(Free,Free);
            fprintf(1, [num2str(toc),'\n']);
            Dx_e = 0;
            
            %==============================================================
            % REORDERING
            %==============================================================
            fprintf(1, 'REORDERING:         '); tic;
            switch reorder
                case 'metis'
                    perm = metis(A);
                case 'amd'
                    perm = amd(A);
                otherwise
                    error('Unknown reordering')
            end
            fprintf(1, [num2str(toc),'\n']);

            %==============================================================
            % FACTORIZATION - ideally L = lchol(K, perm)
            %==============================================================
            fprintf(1, 'FACTORIZATION:      '); tic;
            A = cs_transpose(A);
            A = cs_symperm(A,perm);
            A = cs_transpose(A);
            L = lchol(A);
            fprintf(1, [num2str(toc,'%8.6f'),'\n']);

            %==============================================================
            % POWELL-HESTENES ITERATIONS - UZAWA-like - PRESSURE
            %==============================================================
            tic
            div_max_uz  = 1e-16; div_max     = realmax;
            uz_iter     =     0; uz_iter_max =       5;

            Pressure    = zeros(nel*np, 1);
            Vel         = zeros(sdof  , 1);
            Vel(Bc_ind) = Bc_val;

            while (uz_iter < uz_iter_max)       %div_max>div_max_uz  && 
                uz_iter         = uz_iter + 1;
                Vel(Free(perm)) = cs_ltsolve(L,cs_lsolve(L,Rhs(Free(perm))));          %BACK & FORWARD SUBS
                Div             = invM*(Q*Vel);                                        %COMPUTE QUASI-DIVERGENCE
                Rhs             = Rhs - PF*(Q'*Div);                                   %UPDATE RHS
                Pressure        = Pressure + PF*Div;                                   %UPDATE TOTAL PRESSURE (negative sign convention)
                div_max         = max(abs(Div(:)));                                    %CHECK INCOMPRESSIBILITY
                disp([' PH_ITER: ', num2str(uz_iter), ' ', num2str(div_max)]);
            end
            Pressure = reshape(Pressure, np, nel);                                     % [3,nel] discontinuous
            fprintf(1, 'P-H ITERATIONS:     ');
            fprintf(1, [num2str(toc,'%8.6f'),'\n']);   
            
        case 'free' % default
    
            %==============================================================
            % BOUNDARY CONDITIONS FOR A FREE SURFACE
            %==============================================================
            fprintf(1, 'BDRY CONDITIONS:    '); tic;
            
            % Free surface stabilisation algorithm
            % ------------------------------------
            MESH.Point_id   = Point_id;
            MESH.EL2NOD     = ELEM2NODE;
            MESH.GCOORD     = GCOORD;
            MESH.nnod       = nnod;
            [AX_FS_SPARSE,AY_FS_SPARSE,Dx_e] = fssa(MESH,Rho,G,dt,nip, ...
                fs_alpha,fs_beta, max(Point_id)-1,Load_el,rho_ext);
            
            A = A - AY_FS_SPARSE; % Subtracts from A the symmetric part of the free surface correction.
            
            if bc_t == "cdr"
                [~,AYY_WINK,~] = fssa(MESH,Rho,G,dt,nip, ...
                    0,0,1,[],0);
                A = A - AYY_WINK;
            end
                
            Free = 1:sdof;
            Free(Bc_ind_fs) = []; % Boundary conditions without taking into
                % account the free surface.
            A_BC_UPPER = cs_transpose(A);
            % If removed both following lines it is left as MILAMIN was
            % with the duplicated diagonal. The duplicated diagonal turns
            % out not to be a problem since TMP is only used to subtract 
            % the boundary conditions from the Rhs. In this case duplicated
            % diagonal terms are not used for such a thing:
            %
            % Example:
            %  / k11 k12         \   / v1 \     / F1 \
            % |  k21 k22 k23      | |  v2  | = |  F2  |
            % |      k32 k33 k34  | |  v3  |   |  F3  |
            %  \         k43 k44 /   \ v4 /     \ F4 /
            %
            % If v1 is fixed:
            % 
            %  / k22 k23     \   / v2 \     / F2 - k21v1 \
            % |  k32 k33 k34  | |  v3  | = |      F3      |
            %  \     k43 k44 /   \ v4 /     \     F4     /
            %
            % Note that only off-diagonal terms are used for the boundary
            % conditions (k21), therefore we can save computing time by not
            % substracting the duplicated diagonal term, and lines below
            % can be negleted as MILAMIN does. See in the solver when 
            %         Vel(Free(perm)) =
            %             cs_ltsolve(L,cs_lsolve(L,Rhs(Free(perm))));
            % the Free inside the Rhs gets only the Rhs that are not part
            % of the boundary conditions and therefore it does not get
            % duplicated diagonal terms.
            %-----------------------------------------------------------------|
            % Old way of doing it                                             |
            %A_BC_UPPER(1:sdof,1:sdof) = A_BC_UPPER - diag(diag(A_BC_UPPER));%|
            % New way setting 0 the diagonal without substraction             |
            %A_BC_UPPER(speye(sdof)==1) = 0;                                 %|
            %-----------------------------------------------------------------|
            TMP         = A(:,Bc_ind_fs) + A_BC_UPPER(:,Bc_ind_fs);
            Rhs         = Rhs - TMP*Bc_val_fs';
            A           = A(Free,Free);
                
            fprintf(1, [num2str(toc),'\n']);
            
            %     ind_Asparse_example = find(A~=0); % Uncomment this for plotting the
            %         % elements of the matrix A. Set a break point in 'Break point
            %         % here'.
            %     
            %     [iAsparse,jAsparse] = ind2sub(size(A),ind_Asparse_example);
            %     figure(4)
            %     plot(jAsparse,iAsparse,'.');
            %     xlabel('J')
            %     ylabel('I')
            %     axis ij 
            %     % Break point here
            
            %==============================================================
            % REORDERING
            %==============================================================
            fprintf(1, 'REORDERING:         '); tic;
            switch reorder
                case 'metis'
                    perm = metis(A);
                case 'amd'
                    perm = amd(A);
                otherwise
                    error('Unknown reordering')
            end
            fprintf(1, [num2str(toc),'\n']);

            %==============================================================
            % FACTORIZATION - ideally L = lchol(K, perm)
            %==============================================================
            fprintf(1, 'FACTORIZATION:      '); tic;
            A = cs_transpose(A);
            A = cs_symperm(A,perm);
            A = cs_transpose(A);
            L = lchol(A);
            fprintf(1, [num2str(toc,'%8.6f'),'\n']);

            %==============================================================
            % POWELL-HESTENES ITERATIONS - UZAWA-like - PRESSURE
            %==============================================================
            if fs_beta > 0
                tic
                div_max_uz  = 1e-16; div_max     = realmax;
                uz_iter     =     0; uz_iter_max =       5;
                fs_iter     =     0; fs_iter_max =       10;
                
                Pressure    = zeros(nel*np, 1);
                Vel         = zeros(sdof  , 1);
                Vel(Bc_ind_fs) = Bc_val_fs;
                Vel_sav = zeros(sdof,fs_iter_max);
                vel_dif = 0;
                
                Rhs0 = Rhs;
                
                while fs_iter<fs_iter_max
                    fs_iter = fs_iter+1;
                    Rhs = Rhs0 + AX_FS_SPARSE*Vel;
                    uz_iter = 0;
                    
                    while (uz_iter<uz_iter_max) %div_max>div_max_uz  &&
                        uz_iter         = uz_iter + 1;
                        Vel(Free(perm)) = cs_ltsolve(L,cs_lsolve(L,Rhs(Free(perm))));
                        %BACK & FORWARD SUBS
                        Div             = invM*(Q*Vel);
                        %COMPUTE QUASI-DIVERGENCE
                        Rhs             = Rhs - PF*(Q'*Div);
                        %UPDATE RHS
                        Pressure        = Pressure + PF*Div;
                        %UPDATE TOTAL PRESSURE (negative sign convention)
                        div_max         = max(abs(Div(:)));
                        %CHECK INCOMPRESSIBILITY
                        %disp([' PH_ITER: ', num2str(uz_iter), ' ', num2str(div_max)]);
                        % In order to check convergence uncomment the
                        % following:
                        if uz_iter==uz_iter_max;
                            Vel_sav(:,fs_iter) = Vel;
                        end
                    end
                    if fs_iter>1
                        vel_dif = max(abs(Vel_sav(:,fs_iter)-Vel_sav(:,fs_iter-1)));
                    end
                    disp([' FS_ITER: ', num2str(fs_iter), ' ', num2str(div_max),' Vel diff: ',num2str(vel_dif)]);
                end
                Pressure = reshape(Pressure,np, nel);
                fprintf(1, 'P-H ITERATIONS:     ');
                fprintf(1, [num2str(toc,'%8.6f'),'\n']);
                
            else % if fs_beta = 0
                tic
                div_max_uz  = 1e-16; div_max     = realmax;
                uz_iter     =     0; uz_iter_max =       5;
                
                Pressure    = zeros(nel*np, 1);
                Vel         = zeros(sdof  , 1);
                Vel(Bc_ind_fs) = Bc_val_fs;
                
                while (uz_iter<uz_iter_max) %div_max>div_max_uz  &&
                    uz_iter         = uz_iter + 1;
                    Vel(Free(perm)) = cs_ltsolve(L,cs_lsolve(L,Rhs(Free(perm))));
                    %BACK & FORWARD SUBS
                    Div             = invM*(Q*Vel);
                    %COMPUTE QUASI-DIVERGENCE
                    Rhs             = Rhs - PF*(Q'*Div);
                    %UPDATE RHS
                    Pressure        = Pressure + PF*Div;
                    %UPDATE TOTAL PRESSURE (negative sign convention)
                    div_max         = max(abs(Div(:)));
                    %CHECK INCOMPRESSIBILITY
                    disp([' PH_ITER: ', num2str(uz_iter), ' ', num2str(div_max)]);
                    % In order to check convergence uncomment the
                    % following:
                    %                         if uz_iter==uz_iter_max;
                    %                             Vel_sav(:,fs_iter) = Vel;
                    %                         end
                end
                Pressure = reshape(Pressure, np, nel);
                fprintf(1, 'P-H ITERATIONS:     ');
                fprintf(1, [num2str(toc,'%8.6f'),'\n']);
            end  % if fs_beta = 0
        otherwise
            error("mechanical2d_m:: non accepted top_surface value")    
    end % switch top_surface
    
    %======================================================================
    % UPDATE STRESS AND STRAINS
    %======================================================================

    tic
    
    %======================================================================
    % SMOOTH PRESSURES
    %======================================================================
    % Pc = call_el2nod_pressure(GCOORD,ELEM2NODE,nel,Pressure);
    
    %======================================================================
    % DECLARE VARIABLES (ALLOCATE MEMORY)
    %======================================================================
    nelblo      = SOLVER.nelblo;
    nelblo      = min(nel, nelblo);

    il          = 1;
    iu          = nelblo;

    % T2all(il:iu,ip)      = zeros(nel, nip);
    % YIELD(il:iu,ip)      = zeros(nel, nip);

    %==================================================================
    % i) BLOCK LOOP -
    %==================================================================
    for ib = 1:nblo
        %==============================================================
        % ii) FETCH DATA OF ELEMENTS IN BLOCK
        %==============================================================
        invJx  = zeros(nelblo, ndim);
        invJy  = zeros(nelblo, ndim);

        ECOORD_x        = reshape( GCOORD(1,ELEM2NODE(:,il:iu)), nnodel, nelblo);
        ECOORD_y        = reshape( GCOORD(2,ELEM2NODE(:,il:iu)), nnodel, nelblo);
        Vel_block       = reshape(Vel(ELEM_DOF(:,il:iu))',  nelblo, 2*nnodel);       % [nelblk, ndim*nnodel] :: (:, v_1_x v_1_z ... v_nnodel_z)
        %PRESSURE
        Pi_block        = zeros(nelblo, np);  % pressure shape functions

        a23   = ECOORD_x(2,:).*ECOORD_y(3,:) - ECOORD_x(3,:).*ECOORD_y(2,:);
        a31   = ECOORD_x(3,:).*ECOORD_y(1,:) - ECOORD_x(1,:).*ECOORD_y(3,:);
        a12   = ECOORD_x(1,:).*ECOORD_y(2,:) - ECOORD_x(2,:).*ECOORD_y(1,:);
        area  = a23 + a31 + a12;

        for ip=1:nip
            %==========================================================
            % iv) LOAD SHAPE FUNCTIONS DERIVATIVES FOR INTEGRATION POINT
            %==========================================================
            ED     = Mu_all(il:iu,ip);
            %             Mu_ve    = Mu_ve_all(il:iu,ip);
            Etheta = THETA_all(il:iu,ip);
            
            Ni     =        N{ip};
            dNdui  =     dNdu{ip};
            GIP_x  = Ni'*ECOORD_x; %x and y coordinate of the integrations point
            GIP_y  = Ni'*ECOORD_y;

            tmp   = ECOORD_x(3,:).*GIP_y - GIP_x.*ECOORD_y(3,:);
            eta1  = a23 + tmp + ECOORD_y(2,:).*GIP_x - GIP_y.*ECOORD_x(2,:);
            eta2  = a31 - tmp + ECOORD_x(1,:).*GIP_y - GIP_x.*ECOORD_y(1,:);

            %pressure shape functions
            Pi_block(:,1) = eta1./area;
            Pi_block(:,2) = eta2./area;
            Pi_block(:,3) = 1 - Pi_block(:,1) - Pi_block(:,2);

            %pressure at the integration point
            Pres_ip_block = sum(Pi_block .* Pressure(:,il:iu)',2); 
            %==========================================================
            % v) CALCULATE JACOBIAN, ITS DETERMINANT AND INVERSE
            %==========================================================
            Jx          = ECOORD_x'*dNdui;                            % [nelblk,2] 
            Jy          = ECOORD_y'*dNdui;                            % [nelblk,2]
            detJ        = Jx(:,1).*Jy(:,2) - Jx(:,2).*Jy(:,1);        % [nelblk,1]

            invdetJ     = 1.0./detJ;
            invJx(:,1)  = +Jy(:,2).*invdetJ;
            invJx(:,2)  = -Jy(:,1).*invdetJ;
            invJy(:,1)  = -Jx(:,2).*invdetJ;
            invJy(:,2)  = +Jx(:,1).*invdetJ;

            %==========================================================
            % vi) DERIVATIVES wrt GLOBAL COORDINATES
            %==========================================================
            dNdx            = invJx*dNdui';                          % [nelblk,nnodel]
            dNdy            = invJy*dNdui';                          % [nelblk,nnodel]

	        B               = zeros(nelblo, 2*nnodel);
            B(:,1:2:end-1)  = dNdx;
            Erx_block       = sum(B.*Vel_block,2);
            B(:)            = 0;
            B(:,2:2:end)    = dNdy;
            Ery_block       = sum(B.*Vel_block,2);
            B(:)            = 0;
            B(:,1:2:end-1)  = dNdy;
            B(:,2:2:end)    = dNdx;
            Erxy_block      = sum(B.*Vel_block,2);
            
            Wxy_block       = 0.5.*(sum(dNdx.*Vel_block(:,2:2:end),2) - sum(dNdy.*Vel_block(:,1:2:end-1),2)); % angular velocity (= half vorticity)
            
            Txx_block       =  4/3.*ED.*Erx_block -2/3.*ED.*Ery_block + 2.*ED.*Etheta.*TAU_xx_old(il:iu,ip);           % Thus TAUxx_old == TAUyy_old = 0 => Txx + Tyy = 0 if and only if Erx + Ery = 0. 
            Tyy_block       = -2/3.*ED.*Erx_block +4/3.*ED.*Ery_block + 2.*ED.*Etheta.*TAU_yy_old(il:iu,ip);           % " "
            Txy_block       = ED.*Erxy_block + 2.*ED.*Etheta.*TAU_xy_old(il:iu,ip);
            
            if isfield(SOLVER,"symmetric_tau")
                if SOLVER.symmetric_tau              % impose theoretical constraint than sum of normal deviatoric stresses is 0.
                    Txx_block = sign(Txx_block) .* (abs(Txx_block) + abs(Tyy_block)) / 2.;
                    Tyy_block = - Txx_block;
                end
            end

            Erxy_block      = Erxy_block./2;
            % IInd invariant of (effective?) strain rates
            E2it_block      = sqrt(0.5*(Erx_block.^2+Ery_block.^2) + Erxy_block.^2);
            
            %2nd invariant of the deviatoric stress
            T2it_block    = sqrt(0.5*(Txx_block.^2+Tyy_block.^2)+Txy_block.^2);
            
            %TAU(il:iu,ip) = sqrt(1/2*(TAU_xx_old(il:iu,ip).^2 + TAU_yy_old(il:iu,ip).^2) + TAU_xy_old(il:iu,ip).^2);
            
            %Difference in principal stresses (2 ways of calculating)
            
            %             Diff_ps_block1  = 2*sqrt((Txx_block - Tyy_block).^2 / 4 + Txy_block.^2);
            Max_Tau_block   = 2*ED.*E2it_block; % maximum shear stress
            %equal to difference in principal stresses divided by two /Turcotte Schubert page 83
            % principal stresses and
            % strain rate invariants
            % page 36 and 76 in
            % Ranalli's book
            
            % currently not used:
            %Yield_block                         = -GIP_y'*9.81*2700*sin(Phi0) + Cohesion0; %sin(Phi).*Pres_ip_block + Cohesion;% 200*1e6;%sin(Phi).*Pres_ip_block + Cohesion;
            %Yield_block(Yield_block<Cohesion0)  = Cohesion0;
            
            %--------------------------------------------------------------
            % STRAIN SOFTENING: Calculating and accumulating I2 MODIFY
            %--------------------------------------------------------------
            % Calculates the spatial gradients of the velocity, L (Malvern p. 146) 
            
            dVxdx = sum(Vel_block(:,1:2:end-1).*dNdx,2);                       % [nelblk,1] == Er_xx_block       L matrix := | Vxx Vxy | :: spatial velocity gradient
            dVxdy = sum(Vel_block(:,1:2:end-1).*dNdy,2);                       % [nelblk,1]                                  | Vyx Vyy |
            dVydx = sum(Vel_block(:,2:2:end  ).*dNdx,2);                       % [nelblk,1]
            dVydy = sum(Vel_block(:,2:2:end  ).*dNdy,2);                       % [nelblk,1] == Er_zz_block
            % update the historic change of deformation gradient F (Malvern, 1969) for the current strain iteration
            %                @previous time  + dt * \dot(F) (where \dot(F) is the rate of change (i.e .time derivative) of the deformation gradient F; \dot{F}= L * F) [Eq 4.5.14 Malvern]
            F_xx(il:iu,ip) = F0_xx(il:iu,ip) + dt * (dVxdx.*F0_xx(il:iu,ip) + dVxdy.*F0_yx(il:iu,ip)); %  dt * | Vxx Vxy | | Fxx Fxy |, where F is the "material deformation gradient tensor"
            F_xy(il:iu,ip) = F0_xy(il:iu,ip) + dt * (dVxdx.*F0_xy(il:iu,ip) + dVxdy.*F0_yy(il:iu,ip)); %       | Vyx Vyy | | Fyx Fyy | 
            F_yx(il:iu,ip) = F0_yx(il:iu,ip) + dt * (dVydx.*F0_xx(il:iu,ip) + dVydy.*F0_yx(il:iu,ip)); %
            F_yy(il:iu,ip) = F0_yy(il:iu,ip) + dt * (dVydx.*F0_xy(il:iu,ip) + dVydy.*F0_yy(il:iu,ip)); %

            % Calculates the accumulated strain tensors :: E = 1/2 (F'*F - I) [Eq. 4.5.8 Malvern]
            E_xx = (1/2)*(F_xx(il:iu,ip).^2 + F_yx(il:iu,ip).^2 - 1);
            E_xy = (1/2)*(F_xx(il:iu,ip).*F_xy(il:iu,ip) + F_yx(il:iu,ip).*F_yy(il:iu,ip));
            E_yy = (1/2)*(F_xy(il:iu,ip).^2 + F_yy(il:iu,ip).^2 - 1);
            
            %             % Calculates the deviatoric strain (not needed because incompressibility)
            %             Ed_xx = E_xx - (1/2).*(E_xx+E_yy);
            %             Ed_yy = E_yy - (1/2).*(E_xx+E_yy);
            
            % Calculates the second invariant of the accumulated deviatoric strain I2
            I2.f(il:iu,ip) = sqrt((1/2)*(E_xx.^2 + E_yy.^2) + E_xy.^2);
            clear("E_xx","E_xy","E_yy")

            %--------------------------------------------------------------
            % STRAIN SOFTENING: Calculating and accumulating PI2 MODIFY
            %-------------------------------------------------------------- 
            if ei~=1 && PHY.SS.opt=="partI"
                yield_blk    = YC(il:iu,ip) & P_Phases(il:iu)';
                gamma_blk    = Gamma(il:iu,ip);
                Yield_blk_T2 = Yield_T2(il:iu,ip);
                Strain_xx_pl = zeros(length(il:iu),1);                                % only the II invariant of Strain_pl [plastic strainrate] is relevant here
                Strain_xy_pl = zeros(length(il:iu),1);
                Strain_yy_pl = zeros(length(il:iu),1);
                Strain_II_pl = zeros(length(il:iu),1);
                
                Strain_xx_pl(yield_blk~=0) = 0.5 * gamma_blk(yield_blk~=0)...
                    .*Txx_block(yield_blk~= 0)./Yield_blk_T2(yield_blk~= 0);
                Strain_xy_pl(yield_blk~=0) = 0.5 * gamma_blk(yield_blk~=0)...
                    .*Txy_block(yield_blk~=0)./Yield_blk_T2(yield_blk~= 0);
                Strain_yy_pl(yield_blk~=0) = 0.5 * gamma_blk(yield_blk~= 0)...
                    .*Tyy_block(yield_blk~=0)./Yield_blk_T2(yield_blk~=0);
                
                Strain_II_pl    = sqrt(0.5*(Strain_xx_pl.^2+ ...
                    Strain_yy_pl.^2)+Strain_xy_pl.^2);
                I2.p(il:iu,ip)  = I2p0(il:iu,ip) + dt .* Strain_II_pl; 
                
                Ec_block        = 0.5.*T2it_block./Mu_c_all(il:iu,ip);
                I2.c(il:iu,ip)  = I2c0(il:iu,ip) + dt .* Ec_block;
            end
            %==============================================================
            % ix) WRITE DATA INTO GLOBAL STORAGE
            %==============================================================
            STRAIN_xx(il:iu,ip)  = Erx_block;                                                      % strain rates
            STRAIN_yy(il:iu,ip)  = Ery_block;                                                      % "       "
            STRAIN_xy(il:iu,ip)  = Erxy_block;                                                     % "       "
            TAU_xx(il:iu,ip)     = Txx_block;                                                      % deviatoric stresses   
            TAU_yy(il:iu,ip)     = Tyy_block;                                                      % "       "
            TAU_xy(il:iu,ip)     = Txy_block;                                                      % "       "
            E2it(il:iu,ip)       = E2it_block;                                                     % 2nd invariant of strain rates above 
            T2all(il:iu,ip)      = T2it_block;                                                     % 2nd invariant of the deviatoric stresses above [used for yield evaluation]
            % MAX_TAU(il:iu,ip)    = Max_Tau_block;
            PRES_IP(il:iu,ip)    = Pres_ip_block;                                                  % pressure at integration points
            GIP_x_all(il:iu,ip)  = GIP_x';                                                         % x coordinate of integration points
            GIP_y_all(il:iu,ip)  = GIP_y';                                                         % y coordinate of integration points
            %YIELD(il:iu,ip)      = Yield_block;
            W_xy(il:iu,ip)       = Wxy_block;                                                      % angular velocity (= half vorticity)
            
        end % for ip


        %==============================================================
        % READJUST START, END AND SIZE OF BLOCK. REALLOCATE MEMORY
        %==============================================================
        il = iu + 1;
        if ib == nblo-1
            nelblo = nel - iu;
        end
        iu = iu + nelblo;
    end % block loop
    
    %======================================================================
    % Plots
    %======================================================================
    if PLOT.make
        PLOT.sec = 'postprocessing';
        plot_handle();
    end
    
    %======================================================================
    % Closes the loop
    %======================================================================
    strain(ei) = max(abs(E2all(:) - E2it(:)))./max(E2it(:));
    %strain(ei) = max((abs(E2all(:) - E2it(:)))./(E2it(:)))
    E2all      = E2it;
    %     if (strain(end)<1e-3 && ei>1  || ei>SOLVER.nitmax || ...
    %         (all(RHEOL.Ndis(:)==1) && PLASTICITY.type == "no"))
    %         ei
    %         break
    %     end
    disp(['RESIDUAL:           ',num2str(Residual(end))])
    if (Residual(end)<SOLVER.tol_Exit && ei>1  || ei>SOLVER.nitmax || ...
        (all(RHEOL.Ndis(:)==1) && PLASTICITY.type == "no"))
        ei
        break
    end
    
    % ISBRITTLE = (MAX_TAU - YIELD);%zeros(size(YIELD));
    % %ISBRITTLE(MAX_TAU>=YIELD) = 1;
    %     %PLOTTING FOR TESTING
    %     EL2N      = zeros(nel,3);
    %     GCOORD_N  = zeros(2,nel*3);
    %     TAU_xxn   = zeros(3,nel);
    %     Mu_n      = zeros(3,nel);
    %     Pres_n    = zeros(3,nel);
    %     Brittle   = zeros(3,nel);
    %     EL2N(1,:) = 1:3;
    %     nip1 = 6;
    %     nnodel1 = 6;
    %     [IP_X, IP_w]    = ip_triangle(nip1);
    %     [   Nbig]    = shp_triangle(IP_X, nnodel1);
    %     for i=1:nel
    %         is         = (i-1)*3+1; ie = (i-1)*3+3;
    %         GCOORD_N(:,is:ie) = GCOORD(:,ELEM2NODE([1 2 3],i));
    %         EL2N(i,:) = is:ie;
    %         Dummy      = Nbig'\YIELD(i,:)';
    %         TAU_xxn(:,i)= Dummy(1:3);
    %         Dummy      = Nbig'\Mu_all(i,:)';
    %         Mu_n(:,i)= Dummy(1:3);
    %         Dummy      = Nbig'\PRES_IP(i,:)';
    %         Pres_n(:,i)= Dummy(1:3);
    %         Dummy      = Nbig'\ISBRITTLE(i,:)';
    %         Brittle(:,i)= Dummy(1:3);
    %
    % 
    %     end
    %     figure(1), clf,
    %     patch('faces',EL2N,'vertices',GCOORD_N','facevertexcdata',Pres_n(:),'FaceColor','flat')
    %     shading interp
    %     colorbar
    % %     hold on
    % %     quiver(GCOORD(1,:), GCOORD(2,:), Vel(1:2:end-1)', Vel(2:2:end)','w');
    %     axis tight
    %     drawnow
    %  
    % %     axis([-6 6 -6 6])
    %     title('Pressure')
    %     drawnow
    %     figure(2),clf
    %     patch('faces',EL2N,'vertices',GCOORD_N','facevertexcdata',TAU_xxn(:),'FaceColor','flat')
    %     shading interp
    %     colorbar
    %     title('yield stress')
    %     drawnow
    %     figure(3),clf
    %     patch('faces',EL2N,'vertices',GCOORD_N','facevertexcdata',log(Mu_n(:)),'FaceColor','flat')
    % %     shading interp
    %     colorbar
    %     title('Mu')
    %     drawnow
    % 
    %     figure(4),clf
    %     patch('faces',EL2N,'vertices',GCOORD_N','facevertexcdata',Brittle(:),'FaceColor','flat')
    %     shading interp
    %     colorbar
    %     title('diff yield stress')
    %     drawnow
    % 
    
    %     plot(GIP_x_all(YC==1 & ~YC_old)/1000,GIP_y_all(YC==1 & ~YC_old)/1000,'.','Color',rand(1,3))
    %     axis([-250 500 -400 5])
    %     hold on
    disp(['Number of yielding ips: ',num2str(sum(YC(:)))])
    disp(' ')
    %     YC_old = YC;

    fprintf(1, 'STRESS STRAINS UPDATE:     ');
    fprintf(1, [num2str(toc,'%8.6f'),'\n']);
    
    %CLEAR UNUSED VARIABLES
    clear 'perm' 'invM' 'Q' 'L' 'A' 'EL2N' 'GCOORD_N';

end  % while true strain iteration 
  % JGP: would it be better to enforce deviatoric stress to sum 0 (as the analytical solution) as?:
  % TAU_xx = sign(TAU_xx) .* (abs(TAU_xx) + abs(TAU_yy)) / 2.;
  % TAU_yy = - TAU_xx;


  nw_it = [nw_it ei];
  Vel = reshape(Vel,2,nnod);
end % function
