function SS = random_damage_WS(PHY, RDWS_factorIP, nel, nip, splot)
% SS = RANDOM_DAMAGE_WS(PHY, RDWS_factorIP, NEL, NIP) modifies softening functions defined in SS
% to start from a random value in the range of SS.Param(1)-SS.Param_var/2
% to SS.Param(1)+SS.Param_var/2 defined outside of the function (i.e.
% instead of starting at SS.Phi(1) the function will start in a range
% between SS.Phi(1)-SS.Phi_var/2 to SS.Phi(1)+SS.Phi_var/2). The new
% softening functions will preserve their slope so that the maximum strain
% before softening stops SS.I2_param also needs to be modified
% SS.I2_param_rand.
%
%--------------------------------------------------------------------------
% Function written by Elena Ross 25-06-2015.
%--------------------------------------------------------------------------
% TODO Calculate randomization also for cohesion and viscous softening

if (nargin < 4)
    splot = false;
end

% Change random parameters so randomized distribution is not repeated
if isfield(PHY.SS,'rng')
    if ~isempty(PHY.SS.rng)
        rng(PHY.SS.rng)
    else
        rng('shuffle');
    end
else
    rng('shuffle');
end

SS = [];
SS.RDWS_factorIP = RDWS_factorIP;

% Calculate random field
SS.Random  = rand(nel,1);
SS.Random  = repmat(SS.Random,1,nip);                                      % [nel,nip]

% Calculate friction angle ranges for softening functions
SS.Phi1_rand = RDWS_factorIP{1}.* (SS.Random * (PHY.SS.Phi_var)) + ...     % [nel,nip]
    (PHY.SS.Phi(1)- PHY.SS.Phi_var*(RDWS_factorIP{1}-1/2));
SS.Phi2_rand = repmat(PHY.SS.Phi(2), size(SS.Random));                     %   "  "

% Calculate accumulated strain ranges for softening functions
SS.I2_phi1_rand = repmat(PHY.SS.I2_phi(1),size(SS.Random));                % [nel,nip]
SS.I2_phi2_rand = (SS.Phi2_rand-SS.Phi1_rand) .* (PHY.SS.I2_phi(2)-PHY.SS.I2_phi(1)) ...
                  ./(PHY.SS.Phi(2)-PHY.SS.Phi(1)) + PHY.SS.I2_phi(1);      % [nel,nip]

% Calculate dislocation ranges for softening functions
SS.Pef_dis1_rand = (RDWS_factorIP{2}-1).*SS.Random + PHY.SS.Pef_dis(1);    % [nel, nip]
SS.Pef_dis2_rand = repmat(PHY.SS.Pef_dis(2), size(SS.Pef_dis1_rand) );     % [nel, nip]

% Calculate accumulated strain ranges for dislocation softening functions
SS.I2_pef_dis1_rand = repmat(PHY.SS.I2_pef_dis(1),size(SS.Random));        % [nel, nip]
SS.I2_pef_dis2_rand = (SS.Pef_dis2_rand-SS.Pef_dis1_rand) ...              % [nel, nip]
    * ((PHY.SS.I2_pef_dis(2)-PHY.SS.I2_pef_dis(1)) ...
      /(PHY.SS.Pef_dis(2) - PHY.SS.Pef_dis(1))) + PHY.SS.I2_pef_dis(1);

% Calculate diffusion ranges for softening functions
SS.Pef_dif1_rand = (RDWS_factorIP{2}-1).*SS.Random + PHY.SS.Pef_dif(1);    % [nel, nip]
SS.Pef_dif2_rand = repmat(PHY.SS.Pef_dif(2), size(SS.Pef_dif1_rand));      % [nel, nip]

% Calculate accumulated strain ranges for diffusion softening functions
SS.I2_pef_dif1_rand = repmat(PHY.SS.I2_pef_dif(1),size(SS.Random));        % [nel, nip]
SS.I2_pef_dif2_rand = (SS.Pef_dif2_rand-SS.Pef_dif1_rand) ...              % [nel, nip]
    * ((PHY.SS.I2_pef_dif(2)-PHY.SS.I2_pef_dif(1)) ...
      /(PHY.SS.Pef_dif(2) - PHY.SS.Pef_dif(1))) + PHY.SS.I2_pef_dif(1);

% % Plot of the element strain softening functions (uncomment for plotting)
if splot
    figure()
    plot([SS.I2_phi1_rand(:,1) SS.I2_phi2_rand(:,1)]', ...
        ([SS.Phi1_rand(:,1) SS.Phi2_rand(:,1)].*180./pi)')
    title('Friction angle strain softening function')
    xlabel('Friction angle [Degrees]')
    ylabel('Accumulated strain')

% % Plot of the element strain softening functions (uncomment for plotting)
    figure()
    plot([SS.I2_pef_dis1_rand(:,1) SS.I2_pef_dis2_rand(:,1)]', ...
        ([SS.Pef_dis1_rand(:,1) SS.Pef_dis2_rand(:,1)])')
    title('Dislocation strain softening function')
    xlabel('Preexponential softening')
    ylabel('Accumulated strain')

% % Plot of the element strain softening functions (uncomment for plotting)
    figure()
    plot([SS.I2_pef_dif1_rand(:,1) SS.I2_pef_dif2_rand(:,1)]', ...
        ([SS.Pef_dif1_rand(:,1) SS.Pef_dif2_rand(:,1)])')
    title('Diffusion strain softening function')
    xlabel('Preexponential softening')
    ylabel('Accumulated strain')
end
