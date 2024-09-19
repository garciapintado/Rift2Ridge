function [Ec_xx,Ec_zz,Ec_xz,Ec_II] = calculate_c_strains(Exx,Ezz,Exz, ...
    TAUxxo,TAUzzo,TAUxzo,Shearm,dt)

Ec_xx = 2.*Exx + TAUxxo./(Shearm.*dt);
Ec_zz = 2.*Ezz + TAUzzo./(Shearm.*dt);
Ec_xz = 2.*Exz + TAUxzo./(Shearm.*dt);

Ec_II = sqrt(0.5.*(Ec_xx.^2 + Ec_zz.^2) + Ec_xz.^2);