function SSo = random_damage(SS,nel,nip)
% SS = RANDOM_DAMAGE(SS,NEL,NIP) modifies softening functions defined in SS
% to start from a random value in the range of SS.Param(1)-SS.Param_var/2
% to SS.Param(1)+SS.Param_var/2 defined outside of the function (i.e.
% instead of starting at SS.Phi(1) the function will start in a range
% between SS.Phi(1)-SS.Phi_var/2 to SS.Phi(1)+SS.Phi_var/2).
% The new softening functions will preserve their slope so that the maximum strain
% before softening stops SS.I2_param also needs to be modified SS.I2_param_rand.
%
% input is SS structure with elements
% SS.Phi_var      : REAL
% SS.Phi          : REAL [2], range of friction angle for strain softening
%
% output:
% SSo structure with keys ['Random','Phi1_rand','Phi2_rand','I2_phi1_rand'], all REAL [nel,nip].
%--------------------------------------------------------------------------
% Function written by Elena Ross 25-06-2015.
%--------------------------------------------------------------------------
% TODO Calculate randomization also for cohesion and viscous softening

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

SSo = [];

% Calculate random field
SSo.Random  = repmat(rand(nel,1), 1, nip);                                      % [nel,nip]

% Friction angle ranges for softening functions
SSo.Phi1_rand = SSo.Random .* (SS.Phi_var) + (SS.Phi(1)-SS.Phi_var/2);          % [nel,nip]
SSo.Phi2_rand = repmat(SS.Phi(2),[size(SSo.Random)]);                           % [nel,nip]

% Accumulated strain ranges for softening functions
SSo.I2_phi1_rand = repmat(SS.I2_phi(1),size(SSo.Random));                        % [nel,nip]
SSo.I2_phi2_rand = (SSo.Phi2_rand-SSo.Phi1_rand) ...                              % [nel,nip]
    .*(SS.I2_phi(2)-SS.I2_phi(1))./(SS.Phi(2)-(SS.Phi(1))) + SS.I2_phi(1);

% % Plot of the element strain softening functions (uncomment for plotting)
% plot([SSo.I2_phi1_rand(:,1) SSo.I2_phi2_rand(:,1)]', ...
%     ([SSo.Phi1_rand(:,1) SSo.Phi2_rand(:,1)].*180./pi)')
% title('Friction angle strain softening function')
% xlabel('Friction angle [Degrees]')
% ylabel('Accumulated strain')
