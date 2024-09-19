function [Phi,Cohesion,Pef_dif,Pef_dis] = strain_softening_SST_rand ...
    (nip,nel,Phi0,Cohesion0,I2,PRES_IP,SS,PHY,Temp,ELEM2NODE,Phases,GCOORD,Point_id)

% TODO scaling of the strain softening? (coded in mechanical2d_cohesion_co.m)
% TODO yielding criterium before strain softening in phi and cohesion? (coded in mechanical2d_cohesion.m)

switch PHY.SS.opt
    case 'fullI'
        I2p = I2.f;
        I2c = I2.f;
    case 'partI'       % default
        I2p = I2.p;
        I2c = I2.c;
    otherwise
        error('Strain softening method not recognised')
end

% Initialize variables
% ====================
% Friction angle for strain softening
Phi = Phi0*ones(nel,nip);
% Cohesion for strain softening
Cohesion = Cohesion0*ones(nel,nip);
Pef_dis = PHY.SS.Pef_dis(1)*ones(nel,nip);       % Preexponential factor for strain softening dislocation
Pef_dif = PHY.SS.Pef_dif(1)*ones(nel,nip);       % Preexponential factor for strain softening difussion

% Switch for the random damage
if PHY.SS.rand_s
    % Assign the random parameters to the phi-softening function
    Phi_1 = SS.Phi1_rand;
    Phi_2 = SS.Phi2_rand;
    I2_phi_1 = SS.I2_phi1_rand;
    I2_phi_2 = SS.I2_phi2_rand;
    
    % Assign the random parameters to the dislocation factor function
    Pef_dis1 = SS.Pef_dis1_rand;
    Pef_dis2 = SS.Pef_dis2_rand;
    I2_pef_dis1 = SS.I2_pef_dis1_rand;
    I2_pef_dis2 = SS.I2_pef_dis2_rand;
    
    % Assign the random parameters to the diffusion factor function
    Pef_dif1 = SS.Pef_dif1_rand;
    Pef_dif2 = SS.Pef_dif2_rand;
    I2_pef_dif1 = SS.I2_pef_dif1_rand;
    I2_pef_dif2 = SS.I2_pef_dif2_rand;
else
    % Assign the constant parameters to the phi-softening function
    Phi_1 = PHY.SS.Phi(1)*ones(size(Phi));
    Phi_2 = PHY.SS.Phi(2)*ones(size(Phi));
    I2_phi_1 = PHY.SS.I2_phi(1)*ones(size(Phi));
    I2_phi_2 = PHY.SS.I2_phi(2)*ones(size(Phi));
    
    % Assign the constant parameters to the dislocation factor function
    Pef_dis1    = PHY.SS.Pef_dis(1)*ones(size(Pef_dis));
    Pef_dis2    = PHY.SS.Pef_dis(2)*ones(size(Pef_dis));
    I2_pef_dis1 = PHY.SS.I2_pef_dis(1)*ones(size(Pef_dis));
    I2_pef_dis2 = PHY.SS.I2_pef_dis(2)*ones(size(Pef_dis));
    
    % Assign the constant parameters to the diffusion factor function
    Pef_dif1    = PHY.SS.Pef_dif(1)*ones(size(Pef_dif));
    Pef_dif2    = PHY.SS.Pef_dif(2)*ones(size(Pef_dif));
    I2_pef_dif1 = PHY.SS.I2_pef_dif(1)*ones(size(Pef_dif));
    I2_pef_dif2 = PHY.SS.I2_pef_dif(2)*ones(size(Pef_dif));
end

% Switch for the viscous-strain-softening dependance on temperature
if PHY.SS.ss_dep_t % default
    [Pef_max_dif,Pef_max_dis,TEMP_IP] = ...                                % {'Pef_max_dif','Pef_max_dis','TEMP_IP'} \in [nel,nip]
        ss_dep_temp(Temp,PHY,ELEM2NODE,nel,nip);
else
    Pef_max_dif = Pef_dif2.*ones(nel,nip);
    Pef_max_dis = Pef_dis2.*ones(nel,nip);
end

Pef_dis1(Pef_dis1>Pef_max_dis) = Pef_max_dis(Pef_dis1>Pef_max_dis);
Pef_dif1(Pef_dif1>Pef_max_dif) = Pef_max_dif(Pef_dif1>Pef_max_dif);

% Calculate softening
% ===================

% Friction angle shoftening
% -------------------------
% Phi where I2 is in the range of PHY.SS.I2_phi or the random interval
% SS.I2_phi#_rand
ind_ss_phi = I2p >= I2_phi_1 & I2p <= I2_phi_2;
Phi(ind_ss_phi) = ...                               % A
    (I2p(ind_ss_phi)-I2_phi_1(ind_ss_phi)).* ...    % B-C
    (Phi_2(ind_ss_phi)-Phi_1(ind_ss_phi))./ ...     % D
    (I2_phi_2(ind_ss_phi)-I2_phi_1(ind_ss_phi))+ ...% E
    Phi_1(ind_ss_phi);                              % F
%       Phi = (I2-rangeI2_1)*rangePhi/rangeI2+rangePhi_1
%       A   = (B -C        )*D       /E      +F

% For I2 bigger than the highest value of I2_phi assigns
% the minimum value of the range Phi
Phi(I2p > I2_phi_2) = Phi_2(I2p > I2_phi_2);

% Cohesion softening
% -------------------
% Cohesion where I2 is in the range of PHY.SS.I2_c
Cohesion(I2p >= PHY.SS.I2_c(1) & I2p <= PHY.SS.I2_c(2)) = ... % A
    (I2p(I2p >= PHY.SS.I2_c(1) & ...                      % B
    I2p <= PHY.SS.I2_c(2))-PHY.SS.I2_c(1))* ...              % B-C
    (PHY.SS.C(2)-PHY.SS.C(1))/ ...                          % D
    (PHY.SS.I2_c(2)-PHY.SS.I2_c(1))+PHY.SS.C(1);                % E+F
%       Cohesion = (I2-rangeI2_1)*rangeC/rangeI2+rangeC_1
%       A        = (B -C        )*D     /E      +F

% For I2 bigger than the highest value of PHY.SS.I2_c assigns
% the minimum value of the range PHY.SS.C
Cohesion(I2p > PHY.SS.I2_c(2)) = PHY.SS.C(2);

% Dislocation shoftening
% ----------------------
ind_ss_dis = I2c >= I2_pef_dis1 & I2c <= I2_pef_dis2;
Pef_dis(ind_ss_dis) = ...                                       % A
    (I2c(ind_ss_dis)-I2_pef_dis1(ind_ss_dis)).* ...                      % B-C
    (Pef_max_dis(ind_ss_dis)-Pef_dis1(ind_ss_dis))./ ...                % D
    (I2_pef_dis2(ind_ss_dis)-I2_pef_dis1(ind_ss_dis))+Pef_dis1(ind_ss_dis);          % E+F
%       Pef_dis = (I2-rangeI2_1)*rangePef_dis/rangeI2+rangePef_dis_1
%       A   = (B -C        )*D       /E      +F

% For I2 bigger than the highest value of PHY.SS.I2_pef assigns
% the maximum value of the range PHY.SS.Pef_dis
Pef_dis(I2c > I2_pef_dis2) = Pef_max_dis(I2c > I2_pef_dis2);

% Diffusion shoftening
% --------------------
ind_ss_dif = I2c >= I2_pef_dif1 & I2c <= I2_pef_dif2;
Pef_dif(ind_ss_dif) = ...                                       % A
    (I2c(ind_ss_dif)-I2_pef_dif1(ind_ss_dif)).* ...                      % B-C
    (Pef_max_dif(ind_ss_dif)-Pef_dif1(ind_ss_dif))./ ...                % D
    (I2_pef_dif2(ind_ss_dif)-I2_pef_dif1(ind_ss_dif))+Pef_dif1(ind_ss_dif);          % E+F
%       Pef_dif = (I2-rangeI2_1)*rangePef_dif/rangeI2+rangePef_dif_1
%       A   = (B -C        )*D       /E      +F

% For I2 bigger than the highest value of PHY.SS.I2_pef_dif assigns
% the maximum value of the range PHY.SS.Pef_dif
Pef_dif(I2c > I2_pef_dif2) = Pef_max_dif(I2c > I2_pef_dif2);

% % Plot (uncomment)
% if PHY.SS.ss_dep_t
%     % Check
%     check_pef_dep_T % For this to work pef_softening_function_check has to
%         % used in the main
%     %plot_pef % Uncomment if plot is not needed but the check is needed
% end

% Softening phases
% ================

% Prepare a boolean vector to remove strain softening from
% the phases that do not need strain softening (the ones
% that are not in SS.Phases_'param')
Phases_bool_phi = ismember(Phases,PHY.SS.Phases_phi);
Phases_bool_c = ismember(Phases,PHY.SS.Phases_c);
Phases_bool_pef_dis = ismember(Phases,PHY.SS.Phases_pef_dis);
Phases_bool_pef_dif = ismember(Phases,PHY.SS.Phases_pef_dif);

% Remove strain softening
Phi(~Phases_bool_phi,:) = Phi0;
Cohesion(~Phases_bool_c,:) = Cohesion0;
Pef_dis(~Phases_bool_pef_dis,:) = PHY.SS.Pef_dis(1);
Pef_dif(~Phases_bool_pef_dif,:) = PHY.SS.Pef_dif(1);
