% Initializing variables
Tp = 0:1:1500;
Ep = 0:0.01:1.25;
TP = repmat(Tp',1,length(Ep));
EP = repmat(Ep,length(Tp),1);

% Calculate theoretical functions
% -------------------------------

% Diffusion
A=log(SS.Pef_dif(2))*SS.dif_upp_lim/(SS.dif_low_lim-SS.dif_upp_lim);
Pmax_dif = exp(A*(TP-SS.dif_upp_lim)/SS.dif_upp_lim);
Pmax_dif(TP<=SS.dif_low_lim) = SS.Pef_dif(2);
Pmax_dif(TP>=SS.dif_upp_lim) = SS.Pef_dif(1);
Pp_dif = (EP-SS.I2_pef_dif(1)).*(Pmax_dif-SS.Pef_dif(1)) ...
    ./(SS.I2_pef_dif(2)-SS.I2_pef_dif(1))+SS.Pef_dif(1);
Pp_dif(EP<=SS.I2_pef_dif(1)) = SS.Pef_dif(1);
Pp_dif(EP>=SS.I2_pef_dif(2)) = Pmax_dif(EP>=SS.I2_pef_dif(2));

% Dislocation
A=log(SS.Pef_dis(2))*SS.dis_upp_lim/(SS.dis_low_lim-SS.dis_upp_lim);
Pmax_dis = exp(A*(TP-SS.dis_upp_lim)/SS.dis_upp_lim);
Pmax_dis(TP<=SS.dis_low_lim) = SS.Pef_dis(2);
Pmax_dis(TP>=SS.dis_upp_lim) = SS.Pef_dis(1);
Pp_dis = (EP-SS.I2_pef_dis(1)).*(Pmax_dis-SS.Pef_dis(1)) ...
    ./(SS.I2_pef_dis(2)-SS.I2_pef_dis(1))+SS.Pef_dis(1);
Pp_dis(EP<=SS.I2_pef_dis(1)) = SS.Pef_dis(1);
Pp_dis(EP>=SS.I2_pef_dis(2)) = Pmax_dis(EP>=SS.I2_pef_dis(2));

% Check results with theoretical
% ------------------------------

% Difussion
SS.F_dif = TriScatteredInterp(EP(:),TP(:),Pp_dif(:));

% Dislocation
SS.F_dis = TriScatteredInterp(EP(:),TP(:),Pp_dis(:));

SS.TP = TP;
SS.EP = EP;
SS.Pp_dif = Pp_dif;
SS.Pp_dis = Pp_dis;