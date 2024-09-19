function Hs = shear_heat(TAU_xx,TAU_yy,TAU_xy,TAU_xx_old,TAU_yy_old, ...
    TAU_xy_old,STRAIN_xx,STRAIN_yy,STRAIN_xy,Phases,Shearm,dt)

% HS = SHEAR_HEAT(TAU_XX,TAU_YY,TAU_XY,TAU_XX_OLD,TAU_YY_OLD,TAU_XY_OLD,
% STRAIN_XX,STRAIN_YY,STRAIN_XY,SHEARM,DT) calculates the shear heating HS
% for given deviatoric stresses, old stresses and strains TAU, TAU_OLD and
% STRAIN respectively. Old stresses TAU_OLD are needed to calculate the
% unelastic strain needed for the shear heating formulation

%--------------------------------------------------------------------------
% Function written by Miguel Andres-Martinez, PhD student at Royal Holloway
% University of London, 25-11-2014. Email: mandresmartinez87@gmail.com
%--------------------------------------------------------------------------

% Calculate elastic strain (formula from Moresi, 2007)
CTS = repmat(1./(2*Shearm(Phases)*dt),1,size(TAU_xx,2));
E_xx_e = CTS.*(TAU_xx-TAU_xx_old);
E_yy_e = CTS.*(TAU_yy-TAU_yy_old);
E_xy_e = CTS.*(TAU_xy-TAU_xy_old);

% Calculate unelastic strain
E_xx_ue = STRAIN_xx - E_xx_e;
E_yy_ue = STRAIN_yy - E_yy_e;
E_xy_ue = STRAIN_xy - E_xy_e;

% Calculate shear heating (formula from Gerya, 2010)
Hs = TAU_xx.*E_xx_ue + TAU_yy.*E_yy_ue + 2.*TAU_xy.*E_xy_ue;