 function [dT_sensible_heat, dE_latent_heat, dT_ele_affected] = heat_release(GCOORD, ELEM2NODE, Temp, ign_body, dt, Bct_ind, EL_IP, PHY, nelblk, GEO)
    % function [Temp,dT] = heat_release(GCOORD, ELEM2NODE, ign_body, area_melt, Tmelt)
    % 
    % +++ purpose +++
    % heat release at a number of dikes and sills represented by tracking points
    %
    % OUTPUT
    % Temp :: [1,nnod7] temperature field updated with (possibly diffused)
    %                   temperature increments from emplaced igneous bodies (basalt dikes + gabbro sills)
    %
    % Details:
    %  consider a number of weights linearly decreasing for increasing distance from the center of the melting body 
    %  the sum of weights should sum 1
    %  then the sum of weights for all melting bodies is a scalar, which
    %  need to be fitted so that the integral of incremental heat input
    %  matches the heat relase by the emplaced magma. 
    %  It is assumed that the specific heat times the density is equal at
    %  both sides of the equation. Thus, these are neglected in the
    %  temperature updating
    %
    % Javier Garcia-Pintado, MARUM, 2021
    
    Tmelt = PHY.temperature_release; % [ºC] temperature or allocated melting bodies
    dT_cutoff_d = PHY.dT_cutoff_distance;
    
    [nnodel, nel] = size(ELEM2NODE);
    if nnodel == 7
      ELEM2NODE = ELEM2NODE(1:6,:);
    end
    nnod = max(ELEM2NODE,[],'all');
    GCOORD = GCOORD(:,1:nnod);
    
    if length(Temp) == (length(GCOORD)+nel)
        isTri7 = true;
    else
        isTri7 = false;
    end 
    Temp = Temp(1:nnod);                                                       % Tri6 background temperature field

    % a] get total energy needed to raise the system to the melting body
    % temperature [assume density and specific heat capacity of the rock are equal at both sides of equations]
    deltaT_x_volume_sen = 0.; % sensible T x volume [i.e. sensible heat divided by Cp x rho]
    deltaEnergy_sen = 0.;    % sensible heat
    deltaEnergy_lat = 0.;     % latent crystallization heat
    coo = [];
    nodesin = false(1,nnod);
    if ~isempty(ign_body.dike)
       try 
          tempbody = sample2D(GCOORD, ELEM2NODE, Temp, ign_body.dike.coo);
       catch
          error("heat_release: dike emplaced on top of surface") % figure(); hold on; plotN(Temp, GCOORD, ELEM2NODE); hold on; plot(ign_body.dike.coo(1,:)/1000,ign_body.dike.coo(2,:)/1000,'o')
       end
       deltaT_x_volume_sen = deltaT_x_volume_sen + ign_body.dike.area_tp * sum(Tmelt - tempbody); % [m3.K] [energy / (cp*rho)] | Note: ign_body.dike.area_tp is a scalar representing the volume of each individual tracking point 
       %deltaEnergy_sen = deltaEnergy_sen + 1200 * PHY.Rho_melt * ign_body.dike.area_tp * sum(Tmelt - tempbody);                      % [J] = [J.kg-1.K-1][kg.m-3][m3][K]
       deltaEnergy_lat = deltaEnergy_lat + PHY.Rho_melt * PHY.melt_crystallization_heat * ign_body.dike.area_tp * length(tempbody);  % [J] = [m3][kg.m-3][J.kg-1]
       coo = [coo ign_body.dike.coo];
       nodesin(ign_body.dike.nodesin(1:nnod)) = true;
    end
    if ~isempty(ign_body.sills)
       for i=1:length(ign_body.sills)
           tempbody = sample2D(GCOORD, ELEM2NODE, Temp, ign_body.sills(i).coo);
           deltaT_x_volume_sen = deltaT_x_volume_sen + ign_body.sills(i).area_tp * sum(Tmelt - tempbody);                 % [m3.K] [energy / (cp*rho)]
           %deltaEnergy_sen = deltaEnergy_sen + 1200 * PHY.Rho_melt * ign_body.sills(i).area_tp * sum(Tmelt - tempbody); % [J] = [J.kg-1.K-1][kg.m-3][m3][K]
           deltaEnergy_lat = deltaEnergy_lat + ign_body.sills(i).area_tp * PHY.Rho_melt * PHY.melt_crystallization_heat * length(tempbody); % [J]
           coo = [coo ign_body.sills(i).coo];
           nodesin(ign_body.sills(i).nodesin(1:nnod)) = true;
       end
    end

    % b] get weight field
    d2melt = pdist2(coo',GCOORD','euclidean','Smallest',1);
    w = 1.0 - min(1.0, d2melt / dT_cutoff_d);                              % [1,nnod6]  figure(); plot_tF(w, GCOORD, ELEM2NODE);
    w(nodesin) = 1.0;
    while sum(w) < 1.0
        dT_cutoff_d = dT_cutoff_d + 500;                                   % [m] increase cut-off until enough sampling shares are found                      % hold on; plot_meshF(ELEM2NODE, GCOORD)
        w = 1.0 - min(1.0, d2melt / dT_cutoff_d);
    end
    
    % c] get increment of temperature field for sensible heat and energy field for latent heat
    % dTmax = deltaEnergy_sen / (integrate2D(GCOORD, ELEM2NODE, w) * 1200 * PHY.Rho_melt);         % [ºC] maximum increment of temperature by sensible heat
    dTmax = deltaT_x_volume_sen / integrate2D(GCOORD, ELEM2NODE, w);       % [ºC] maximum increment of temperature by sensible heat
    dT_sensible_heat = dTmax * w;                                          % [1,nnod] [ºC] temperature increment field by sensible heat
    dT_ele_affected = any(dT_sensible_heat(ELEM2NODE) > 0);                % elements affected by increments in sensible heat
    dEmax = deltaEnergy_lat / integrate2D(GCOORD, ELEM2NODE, w);           % [J.m-3] maximum increment of energy by latent heat
    dE_latent_heat = dEmax * w;                                            % [1,nnod] [J.m-3] increment in the energy density field by latent heat

    % figure();  plotN(Temp + dT, GCOORD, ELEM2NODE); hold on; plot_meshF(ELEM2NODE, GCOORD); hold on; plotGEO(GEO,1:10,'')
    % hold on; plot_meshF(ELEM2NODE, GCOORD, dT_ele_affected, "red");
    % d] heat diffusion [consider the increments as linear adition to the background field and allow to diffuse]
    if PHY.dT_diffuse_timeratio > 0.0
        if PHY.dT_diffuse_timeratio > 1.0
            error("PHY.dT_diffuse_timeratio > 1.0")
        end
        if isfield(EL_IP,'K_hyd')                                          % active hydrothermally-enhanced thermal conductivity [only for parametrized cooling
            Cond = EL_IP.K_hyd;
        else
            Cond = EL_IP.K;
        end
        %dTb = dT;
        % these two can alternatively be placen into the fluid circulation module
        % figure(); plot_ip_valF(EL_IP.Hp, GCOORD, ELEM2NODE,GEO(9).coo) % radiogenic heat production
        % figure(); plot_ip_valF(EL_IP.Hs, GCOORD, ELEM2NODE,GEO(9).coo) %  shear heating
        % figure(); plot_ip_valF(EL_IP.Rho_ip .* EL_IP.Cp, GCOORD, ELEM2NODE,GEO(9).coo) %  thermal capacity
        % dTq = thermal2d_m(ELEM2NODE, GCOORD, zeros(1,nnod)', EL_IP.Rho_ip .* EL_IP.Cp, EL_IP.Hp, Cond, EL_IP.Hs, Bct_ind, Bct_ind*0.0, dt*PHY.dT_diffuse_timeratio, 'amd', nelblk);
        dT_sensible_heat = thermal2d_m(ELEM2NODE, GCOORD, dT_sensible_heat, EL_IP.Rho_ip .* EL_IP.Cp, Cond, EL_IP.Hp*0, ...                        % [nel,nip] variables       nnod6 =  max(ELEM2NODE(1:6,:),[],'all');
                 Bct_ind, Bct_ind*0.0, dt*PHY.dT_diffuse_timeratio, 'amd', nelblk);
    end
    if isTri7
        dT_sensible_heat = addCentroidVal(ELEM2NODE(1:6,:), dT_sensible_heat); % [1,nnod7]
    end
end
