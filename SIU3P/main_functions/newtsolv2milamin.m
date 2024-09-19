% Script to clean up the variables from Newton's solver

Vel         = U;
Pressure    = P;

TAU_xx      = TAU.xx;
TAU_yy      = TAU.zz;
TAU_xy      = TAU.xz;
W_xy        = TAU.W_xz;
T2all       = TAU.II;
STRAIN_xx   = ER.xx;
STRAIN_yy   = ER.zz;
STRAIN_xy   = ER.xz;
STRAIN_xx_ip= ER.xx_ip;
STRAIN_yy_ip= ER.zz_ip;
STRAIN_xy_ip= ER.xz_ip;
E2all       = ER.II;
F_xx        = SS.F_xx;
F_xy        = SS.F_xz;
F_yx        = SS.F_zx;
F_yy        = SS.F_zz;
Mu_all      = MU.Visc;
Mu_dif_all  = MU.Viscd;
Mu_dis_all  = MU.ViscD;
Mu_b_all    = MU.ViscP;

PRES_IP     = zeros(nel,nip);
[IP_X,~]    = ip_triangle_m2tri(nip);
[NP,~]      = sf_dsf_tri367_N(IP_X,3,'cell');
MESH.EL2NODP    = reshape(1:(3*MESH.nel),3,[]);
for ip = 1:nip
    PRES_IP(:,ip)   = P(MESH.EL2NODP(1:3,:))'*NP{ip};
end

clear MU YIELD g dP_dz MESH nnodel ifix vfix EL2NODP nUg nPg U P TAU ER ...
    SS.F_xx SS.F_xz SS.F_zx SS.F_zz IP_X NP