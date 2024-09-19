function [Phi,Cohesion,Pef_dif,Pef_dis] = strain_softening_SST ...
    (nip,nel,Phi0,Cohesion0,I2,PRES_IP,PHY,Temp,ELEM2NODE,Phases)

% TODO scaling of the strain softening? (coded in mechanical2d_cohesion_co.m)
% TODO yielding criterium before strain softening in phi and cohesion? (coded in mechanical2d_cohesion.m)

% Initialize variables
% ====================
% Friction angle for strain softening
Phi = Phi0*ones(nel,nip);
% Cohesion for strain softening
Cohesion = Cohesion0*ones(nel,nip);
% Preexponential factor for strain softening dislocation
Pef_dis = ones(nel,nip);
% Preexponential factor for strain softening difussion
Pef_dif = ones(nel,nip);

if PHY.SS.ss_dep_t
    [Pef_max_dif,Pef_max_dis,TEMP_IP] = ...
        ss_dep_temp(Temp,PHY,ELEM2NODE,nel,nip);
else
    Pef_max_dif = PHY.SS.Pef_dif(2)*ones(nel,nip);
    Pef_max_dis = PHY.SS.Pef_dis(2)*ones(nel,nip);
end

% Calculate softening
% ===================

% Friction angle shoftening
% -------------------------
% Phi where I2 is in the range of PHY.SS.I2_phi
Phi(I2 >= PHY.SS.I2_phi(1) & I2 <= PHY.SS.I2_phi(2)) = ...  % A
    (I2(I2 >= PHY.SS.I2_phi(1) & ...                    % B
    I2 <= PHY.SS.I2_phi(2))-PHY.SS.I2_phi(1))* ...          % B-C
    (PHY.SS.Phi(2)-PHY.SS.Phi(1))/ ...                      % D
    (PHY.SS.I2_phi(2)-PHY.SS.I2_phi(1))+PHY.SS.Phi(1);          % E+F
%       Phi = (I2-rangeI2_1)*rangePhi/rangeI2+rangePhi_1
%       A   = (B -C        )*D       /E      +F

% For I2 bigger than the highest value of PHY.SS.I2_phi assigns
% the minimum value of the range PHY.SS.Phi
Phi(I2 > PHY.SS.I2_phi(2)) = PHY.SS.Phi(2);

% Cohesion shoftening
% -------------------
% Cohesion where I2 is in the range of PHY.SS.I2_c
Cohesion(I2 >= PHY.SS.I2_c(1) & I2 <= PHY.SS.I2_c(2)) = ... % A
    (I2(I2 >= PHY.SS.I2_c(1) & ...                      % B
    I2 <= PHY.SS.I2_c(2))-PHY.SS.I2_c(1))* ...              % B-C
    (PHY.SS.C(2)-PHY.SS.C(1))/ ...                          % D
    (PHY.SS.I2_c(2)-PHY.SS.I2_c(1))+PHY.SS.C(1);                % E+F
%       Cohesion = (I2-rangeI2_1)*rangeC/rangeI2+rangeC_1
%       A        = (B -C        )*D     /E      +F

% For I2 bigger than the highest value of PHY.SS.I2_c assigns
% the minimum value of the range PHY.SS.C
Cohesion(I2 > PHY.SS.I2_c(2)) = PHY.SS.C(2);

% Dislocation shoftening
% ----------------------
ind_ss_dis = I2 >= PHY.SS.I2_pef_dis(1) & I2 <= PHY.SS.I2_pef_dis(2);
Pef_dis(ind_ss_dis) = ...                                       % A
    (I2(ind_ss_dis)-PHY.SS.I2_pef_dis(1)).* ...                      % B-C
    (Pef_max_dis(ind_ss_dis)-PHY.SS.Pef_dis(1))./ ...                % D
    (PHY.SS.I2_pef_dis(2)-PHY.SS.I2_pef_dis(1))+PHY.SS.Pef_dis(1);          % E+F
%       Pef_dis = (I2-rangeI2_1)*rangePef_dis/rangeI2+rangePef_dis_1
%       A   = (B -C        )*D       /E      +F

% For I2 bigger than the highest value of PHY.SS.I2_pef assigns
% the maximum value of the range PHY.SS.Pef_dis
Pef_dis(I2 > PHY.SS.I2_pef_dis(2)) = Pef_max_dis(I2 > PHY.SS.I2_pef_dis(2));

% Diffusion shoftening
% --------------------
ind_ss_dif = I2 >= PHY.SS.I2_pef_dif(1) & I2 <= PHY.SS.I2_pef_dif(2);
Pef_dif(ind_ss_dif) = ...                                       % A
    (I2(ind_ss_dif)-PHY.SS.I2_pef_dif(1)).* ...                      % B-C
    (Pef_max_dif(ind_ss_dif)-PHY.SS.Pef_dif(1))./ ...                % D
    (PHY.SS.I2_pef_dif(2)-PHY.SS.I2_pef_dif(1))+PHY.SS.Pef_dif(1);          % E+F
%       Pef_dif = (I2-rangeI2_1)*rangePef_dif/rangeI2+rangePef_dif_1
%       A   = (B -C        )*D       /E      +F

% For I2 bigger than the highest value of PHY.SS.I2_pef_dif assigns
% the maximum value of the range PHY.SS.Pef_dif
Pef_dif(I2 > PHY.SS.I2_pef_dif(2)) = Pef_max_dif(I2 > PHY.SS.I2_pef_dif(2));

% Plot (uncomment)
if PHY.SS.ss_dep_t
    % Check
    check_pef_dep_T % For this to work pef_softening_function_check has to
        % used in the main
    plot_pef % Uncomment if plot is not needed but the check is needed
end

% Softening phases
% ================

% Prepare a boolean vector to remove strain softening from
% the phases that do not need strain softening (the ones
% that are not in PHY.SS.Phases_'param')
Phases_bool_phi = ismember(Phases,PHY.SS.Phases_phi);
Phases_bool_c = ismember(Phases,PHY.SS.Phases_c);
Phases_bool_pef_dis = ismember(Phases,PHY.SS.Phases_pef_dis);
Phases_bool_pef_dif = ismember(Phases,PHY.SS.Phases_pef_dif);

% Remove strain softening
Phi(~Phases_bool_phi) = Phi0;
Cohesion(~Phases_bool_c,:) = Cohesion0;
Pef_dis(~Phases_bool_pef_dis,:) = PHY.SS.Pef_dis(1);
Pef_dif(~Phases_bool_pef_dif,:) = PHY.SS.Pef_dif(1);
