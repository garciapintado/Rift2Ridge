function RHO = dens_temp_Dpl(GCOORD, ELEM2NODE, Phases, GEOn, Basement, ...
       Temp, Dpl, nip, PHY, SP, ocrust)
    % +++ purpose +++
    % calculate density at integration points
    % initial densities (in PHY) are influenced here by thermal expansion and depletion due to melting
    %
    % ELEM2NODE :: [nnodel, nel]
    % Phases    :: [1, nel]
    % Temp      :: [1,nnod6]
    % PHY
    %    .Ce        :: [nphases,1] thermal expansion
    %    .Dens      :: [nphases,1] density
    %    .Dplf      :: [1] depletion factor
    % nip       :: [1] number of integration points
    % Dpl       :: [1,nnod]
    % ocrust    :: [2,:] polyline representing the base of the oceanic crust. Empty for no oceanic crust
    %
    % Details:
    %         This function assumes 3 phases [1:mantle, 2:lower crust, 3:upper crust]
    %         After evaluation of the global density, density is recalculated at oceanic crust and sediment
    %         The reference temperature for sediment and oceanic crust is 0 ºC  
    %         IPs above Basement are considered sediment, with its specific density,
    %
    % last version: 2021-10 JGP
    %               some speed up
    
    [nnodel,nel] = size(ELEM2NODE);

    RHO = zeros(nel,nip);

    Tref = 0.0;                                                               % reference temperature

    [IP_X, ~] = ip_triangle(nip);                                              % [nip,2] normalised integration point coordinates
    [N, ~]    = shp_deriv_triangle(IP_X, nnodel);                              % {nip}[np,1] normalised shape function
    Nm = cell2mat(N');                                                         % [nnodel,nip]       

    if size(N{1},1) == 1
        for ip = 1:nip
            N{ip} = N{ip}(:);                                                  % into column vectors
        end
    end
    T_nel  = reshape(Temp(ELEM2NODE),nnodel,[]);                            % [nnodel,nel] temperature  
    dT_nel = T_nel - Tref;                                                  % [nnodel,nel] delta temperature
    Dp_nel = reshape(Dpl(ELEM2NODE), nnodel,nel);                           % [nnodel,nel] depletion

    Rho_el = reshape(PHY.Dens(Phases),1,nel); % [1,nel]
    Ce_el  = reshape(PHY.Ce(Phases),1,nel);   % [1,nel]
    
    for ip = 1:nip
        dT         = N{ip}' * dT_nel;                                           % [1,nel]
        Dp_elip    = N{ip}' * Dp_nel;                                           % [1,nel]
        RHO(:,ip)  = Rho_el .* (1 - Ce_el .* dT - PHY.Dplf * Dp_elip);          % [1,nel] (Ros_al2017, supporting material)  
    end
    
    % phases defined by IPs [oceanic crust and sediment]
    if ~isempty(ocrust) || SP.make
        T_ip = (Nm' * T_nel)';                                                 % [nel,nip] 
        [gipx,gipy] = ip_coord(GCOORD, ELEM2NODE, nel, nip);                   % [nel,nip] global coordinates of integration points
    end
    
    % oceanic crust density
    if ~isempty(ocrust)
        ocrust_poly = [ocrust fliplr(GEOn(3).coo)];                                   
        isocr = inpoly2([gipx(:) gipy(:)], ocrust_poly');                    % [nel*nip,1]
        if any(isocr)
            if isfield(PHY.OC,'thermal_expansivity')
                alpT = PHY.OC.thermal_expansivity;
            else
                alpT = 0.0;
            end
            RHO(isocr) = PHY.OC.Dens * (1-alpT*T_ip(isocr));                % assumes Tref=0.0 for the oceanic crust density input
        end
    end
    
    % sediment density 
    % TODO: add Athys (1930) parametrization for compaction - but then we should explicitly compact the mesh.
    %       [doable but not a priority here]
    if SP.make                                          
        sedim_poly = [Basement fliplr(GEOn(end-1).coo)];
        issed = inpoly2([gipx(:) gipy(:)], sedim_poly');                     % [nel*nip,1] LOGICAL; whether the IP is sediment
        if any(issed)
            if isfield(SP,'thermal_expansivity')
                alpT = SP.thermal_expansivity;
            else
                alpT = 0.0;
            end
            RHO(issed) = SP.Dens * (1-alpT*T_ip(issed));                     % assumes Tref=0.0 for the sediment density input
        end
    end % if SP.make
end % function
