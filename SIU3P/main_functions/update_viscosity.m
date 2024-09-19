function [MU,dMu_dU,Mu_up,Mu_do] = ...
    update_viscosity(RHEOL,R,plast,Shear,dt,es,P,Er_II,ToII, ...
    Phi,Cohesion,Pef_dis,Pef_dif,T,SOLVER,Er_xx,Er_zz,Er_xz,dNUdx,Ue,dU)
% [MU,MU_C,MU_P,dMU_DU,MU_UP,MU_DO] = UPDATE_VISCOSITY_V1(VISC0,PLAST,P,
% ER_II,PHI,COHESION,DNUDX,UE,DU) calculates effective, Newtonian and
% apparent plastic viscosities MU, MU_C, MU_P for a given Newtonian
% viscosity VISC0, pressure P, second invariant of the strain rate ER_II,
% friction angle PHI, and COHESION. Additionally, if global derivatives of
% the shape functions DNUDZ and velocities UE are given it also calculates
% partial derivatives of the viscosities respect the velocities DMU_DU
% needed for the calculation of the Jacobian for Newton iterations. It also
% uses a phases type variable PLAST to determine if the layer is plastic or
% only viscous. DU is an increment in velocity used for calculating the
% viscous partial derivative numerically instead of analytically, by
% calculating MU_UP and MU_DO which difference is the numerical derivative
% of the viscosity (this lines are uncommented). This function works either
% for a single element or for a block of elements.
%
% Author: Joerg Hasenclever, University of Bremen
%         joerg3@uni-bremen.de
%
% 19/04/2017 MA
    % Added plastic term into the viscosities and derivatives
    % Vectorisation of the calculation so that the function can deal with
    % element blocks

% TODO Include again the numerical derivatives (uncomment and switch)
% TODO Include again non-Newtonian term together with plasticity
% TODO Include elasticity term into the viscosities and its derivatives
% TODO Include all viscosities and derivatives into a structure

%==========================================================================
% VARIABLES
%==========================================================================
nelblo  = size(Er_II,1);
mu_max  = SOLVER.mu_max;
mu_min  = SOLVER.mu_min;

if exist('T','var')
    T   = T+273;
else
    T   = zeros(nelblo,1);
end

Adis    = RHEOL.Adis_block;
Ndis    = RHEOL.Ndis_block;
Qdis    = RHEOL.Qdis_block;
Vdis    = RHEOL.Vdis_block;

Adif    = RHEOL.Adif_block;
Ndif    = RHEOL.Ndif_block;
Qdif    = RHEOL.Qdif_block;
Vdif    = RHEOL.Vdif_block;

Var     = RHEOL.Var_block;

nelblo  = size(Er_II,1);

% Initialize viscosities
MU.dis  = zeros(nelblo,1);
MU.dif  = zeros(nelblo,1);
MU.p    = zeros(nelblo,1);
MU.eff  = zeros(nelblo,1);
Sc_dis  = zeros(nelblo,1);
Sc_dif  = zeros(nelblo,1);

%==========================================================================
% CREEP
%==========================================================================
% Loop for rheologic variation in the same phase
for n = 1:size(Ndis,2)
    % Factors for scaling triaxial and uniaxial experiment
    % parameters (GERYA 2010)
    Sc_dis = ...
        1./(2.^((Ndis(:,n)-1)./Ndis(:,n)).* ...
        3.^((Ndis(:,n)+1)./(2*Ndis(:,n))));
    Sc_dif = 1/3;
    
    % Dislocation creep
    MU.dis = MU.dis + Var(:,n).* ...
        (Sc_dis.*(Pef_dis.*Adis(:,n)).^(-1./Ndis(:,n)) ... (Sc_dis.*(Pef_dis.*Adis(:,n)).^(-1./Ndis(:,n)) ...
        .* Er_II.^(1./Ndis(:,n)-1) ...
        .* exp((Qdis(:,n) + P.*Vdis(:,n)) ...
        ./(Ndis(:,n).*R.*T)));
    
    % Diffusion creep
    MU.dif = MU.dif + Var(:,n).* ...
        (Sc_dif.*(Pef_dif.*Adif(:,n)).^(-1./Ndif(:,n)) ... (Sc_dif.*(Pef_dif.*Adif(:,n)).^(-1./Ndif(:,n)) ...
        .* exp((Qdif(:,n)+P.*Vdif(:,n)) ...
        ./(Ndif(:,n).*R.*T)));
end
Large_dis           = MU.dis > mu_max;
MU.dis(Large_dis)   = mu_max;
Large_dif           = MU.dif > mu_max;
MU.dif(Large_dif)   = mu_max;

Small_dis           = MU.dis < mu_min;
MU.dis(Small_dis)   = mu_min;
Small_dif           = MU.dif < mu_min;
MU.dif(Small_dif)   = mu_min;

% % Effective viscosity
% dif = Ndif(:,1)~=0;
% MU.c(dif) = (1./MU.dis(dif) ...
%     + 1./MU.dif(dif) ...
%     ...+ elasticity_s*(1./(Shear_block(dif).*dt))...
%     ).^-1;
% 
% MU.c(~dif) = (1./MU.dis(~dif) ...
%     ...+ elasticity_s*(1./(Shear_block(~dif).*dt))...
%     ).^-1;

%==========================================================================
% PLASTICITY
%==========================================================================
Tyield  = P.*sin(Phi) + Cohesion.*cos(Phi);
MU.p    = Tyield./(2.*Er_II + ToII./(dt*Shear));
Large_p = MU.p > mu_max;
MU.p(Large_p)   = mu_max;

%==========================================================================
% ELASTICITY
%==========================================================================
MU.e    = Shear.*dt.*ones(size(MU.p));

%==========================================================================
% EFFECTIVE VISCOSITY
%==========================================================================
difb        = Ndif(:,1)~=0;
pb          = plast==1;
eb          = ones(size(pb))==es;

MU_dis          = MU.dis;
MU_dif          = MU.dif;
MU_dif(~difb)   = 1;
MU_p            = MU.p;
MU_p(~pb)       = 1;
MU_e            = MU.e;
MU_e(~eb)       = 1;

MU.eff           = MU_dis.*MU_dif.*MU_p.*MU_e./ ...
                    (MU_dif.*MU_p.*MU_e + ...
                     MU_dis.*MU_p.*MU_e.*difb + ...
                     MU_dis.*MU_dif.*MU_e.*pb + ...
                     MU_dis.*MU_dif.*MU_p.*eb);

% if Mu < 1e18
%     Mu = 1e18;
% end
% if Mu > 1e24
%     Mu = 1e24;
% end
if nargin==15
    dMu_dU = [];
    Mu_up  = [];
    Mu_do  = [];
    return
end

% % syms Visc0 Ndis dNdx dNdz Ux Uz 
% % Er_xx   = dNdx*Ux
% % Er_zz   = dNdz*Uz
% % Er_xz   = 0.5*(dNdz*Ux + dNdx*Uz)
% % Er_II   = sqrt( 0.5*(Er_xx.^2 + Er_zz.^2) + Er_xz.^2 )
% % Mu      = Visc0 .* Er_II.^(1./Ndis-1)
% % dMu_dUx = diff(Mu,Ux)
% % dMu_dUz = diff(Mu,Uz)
% % dMu_dUx = (Visc0*(1/Ndis-1)*(dNdz*((Ux*dNdz)/2 + (Uz*dNdx)/2) + Ux*dNdx^2)*(((Ux^2*dNdx^2)/2 + (Uz^2*dNdz^2)/2 + ((Ux*dNdz)/2 + (Uz*dNdx)/2)^2)^(1/2))^(1/Ndis-2))/(2*((Ux^2*dNdx^2)/2 + (Uz^2*dNdz^2)/2 + ((Ux*dNdz)/2 + (Uz*dNdx)/2)^2)^(1/2))
% % dMu_dUz = (Visc0*(1/Ndis-1)*(dNdx*((Ux*dNdz)/2 + (Uz*dNdx)/2) + Uz*dNdz^2)*(((Ux^2*dNdx^2)/2 + (Uz^2*dNdz^2)/2 + ((Ux*dNdz)/2 + (Uz*dNdx)/2)^2)^(1/2))^(1/Ndis-2))/(2*((Ux^2*dNdx^2)/2 + (Uz^2*dNdz^2)/2 + ((Ux*dNdz)/2 + (Uz*dNdx)/2)^2)^(1/2))

%==========================================================================
% VISCOSITY DERIVATIVES
%==========================================================================
if ~iscell(dNUdx)
    % Element-by-element
    nnode   = size(Ue,2);
%     [Er_II,Er_xx,Er_zz,Er_xz] = calc_strain_rates(dNUdx,Ue);
    dNdx     = dNUdx(1,:);
    dNdz     = dNUdx(2,:);
    % Second invariant derivatives
    dEr_IIdx = 1/(2*Er_II).*(Er_xx*dNdx + Er_xz*dNdz);
    dEr_IIdz = 1/(2*Er_II).*(Er_zz*dNdz + Er_xz*dNdx);
    % Plastic derivatives
    if plast
        dMu_pdx  = -Tyield/(2*Er_II.^2).*dEr_IIdx;
        dMu_pdz  = -Tyield/(2*Er_II.^2).*dEr_IIdz;
    else
        dMu_pdx  = zeros(1,nnode);
        dMu_pdz  = zeros(1,nnode);
    end
    % Derivatives of the effective viscosities
    dMu_dUx  = (1/Mu_c+1/Mu_p)^(-2)*(1/Mu_c^2*dMu_cdx + 1/Mu_p^2*dMu_pdx);
    dMu_dUz  = (1/Mu_c+1/Mu_p)^(-2)*(1/Mu_c^2*dMu_cdz + 1/Mu_p^2*dMu_pdz);
    dMu_dU_a = [dMu_dUx;dMu_dUz];
    dMu_dU_a = dMu_dU_a(:)';
else
    % Vectorised
    nelblo   = size(Ue,1);
    nnode    = size(Ue,2)/2;
    rep_mat  = ones(1,nnode);
%     [Er_II,Er_xx,Er_zz,Er_xz] = calc_strain_rates(dNUdx,Ue);
    dNdx     = dNUdx{1};
    dNdz     = dNUdx{2};
    
    % Second invariant derivatives
    %-----------------------------
    dEr_IIdx = 1./(2*Er_II(:,rep_mat)).* ...
        (Er_xx(:,rep_mat).*dNdx + Er_xz(:,rep_mat).*dNdz);
    dEr_IIdz = 1./(2*Er_II(:,rep_mat)).* ...
        (Er_zz(:,rep_mat).*dNdz + Er_xz(:,rep_mat).*dNdx);
    
    % Dislocation derivatives
    %------------------------
    dMu_dis_ct  = zeros(nelblo,1);
    % Loop for rheologic variation in the same phase
    for n = 1:size(Ndis,2)
        % Factors for scaling triaxial and uniaxial experiment
        % parameters (GERYA 2010)
        Sc_dis = ...
            1./(2.^((Ndis(:,n)-1)./Ndis(:,n)).* ...
            3.^((Ndis(:,n)+1)./(2*Ndis(:,n))));
        
        % Dislocation creep
        dMu_dis_ct = dMu_dis_ct + Var(:,n).* ...
            (Sc_dis.*(Pef_dis.*Adis(:,n)).^(-1./Ndis(:,n)) ... (Sc_dis.*(Pef_dis.*Adis(:,n)).^(-1./Ndis(:,n)) ...
            .* Er_II.^(1./Ndis(:,n)-2) ...
            .* exp((Qdis(:,n) + P.*Vdis(:,n)) ...
            ./(Ndis(:,n).*R.*T))).*(1./Ndis(:,n)-1);
    end
    dMu_dis_ct(Large_dis)   = 0;
    dMu_dis_ct(Small_dis)   = 0;
    % Dislocation derivatives
    %------------------------
    dMu_disdx   = dMu_dis_ct(:,rep_mat).*dEr_IIdx;
    dMu_disdz   = dMu_dis_ct(:,rep_mat).*dEr_IIdz;
    
    % Plastic derivatives
    %--------------------
    dMu_pdx     = -Tyield(:,rep_mat)./((2*Er_II(:,rep_mat) + ...
        ToII(:,rep_mat)./(dt.*Shear)).^2).*dEr_IIdx;
    dMu_pdx(Large_p)    = 0;
    dMu_pdx(~plast,:)   = 0;
    dMu_pdz     = -Tyield(:,rep_mat)./((2*Er_II(:,rep_mat) + ...
        ToII(:,rep_mat)./(dt.*Shear)).^2).*dEr_IIdz;
    dMu_pdz(Large_p)    = 0;
    dMu_pdz(~plast,:)   = 0;
    
    % Derivatives of the effective viscosities
    %-----------------------------------------
    MU_EFF  = MU.eff(:,rep_mat);
    MU_D    = MU_dis(:,rep_mat);
    MU_P    = MU_p(:,rep_mat);
    dMu_dUx    = MU_EFF.^2 .* ...
                    (dMu_disdx./(MU_D.^2) + dMu_pdx./(MU_P.^2));
    dMu_dUz    = MU_EFF.^2 .* ...
                    (dMu_disdz./(MU_D.^2) + dMu_pdz./(MU_P.^2));
    
    % Reordering the derivatives of the viscosity
    dMu_dU_a = zeros(nelblo,nnode*2);
    dMu_dU_a(:,1:2:end-1)     = dMu_dUx;
    dMu_dU_a(:,2:2:end)       = dMu_dUz;
end

% nUe     = numel(Ue);
% Mu_up   = zeros(1,nUe);
% Mu_do   = zeros(1,nUe);
% for i=1:nUe
%     Ue_up    = Ue;
%     Ue_up(i) = Ue_up(i) + dU;
%     % Strain rates and strain rate invariant
%     Er_II_up = calc_strain_rates(dNUdx,Ue_up);
%     Mu_up(i) = Visc0 .* Er_II_up.^(1./Ndis-1);
%     
%     Ue_do    = Ue;
%     Ue_do(i) = Ue_do(i) - dU;
%     % Strain rates and strain rate invariant
%     Er_II_do = calc_strain_rates(dNUdx,Ue_do);
%     Mu_do(i) = Visc0 .* Er_II_do.^(1./Ndis-1);
% end
% 
% dMu_dU_n = (Mu_up-Mu_do)./(2*dU');
dMu_dU   = dMu_dU_a;

end