function [EL_IP, EL_EL, Phases_ip, Phases_el] = updateMaterial(EL_IP, GCOORD, EL2NOD, GEOn, Phases, ocrust, Basement, SETTINGS, PHY, SP)
    % +++ purpose +++
    % update material properties evaluated at quadrature points:
    % .Cp :: [J.kg-1.K-1] specific heat capacity
    % .K  :: [W.m-1.K-1] thermal conductivity
    % .Hp :: [W.m-3] radiogenic heat production 

    [nel,nip] = size(EL_IP.Hs);
    EL_EL = [];

    [gipx,gipy] = ip_coord(GCOORD, EL2NOD, size(EL2NOD,2), nip);

    Phases_el = Phases;                                                         % [1,nel] phases for the thermal solver
    Phases_ip = repmat(Phases(:),1,nip);                                        % [nel,nip] for new hydrothermal parameterization (porosity, permeability)
    if ~isempty(ocrust)
        ocrust_poly = [ocrust fliplr(GEOn(3).coo)];
        isocr_qp = reshape(inpoly2([gipx(:) gipy(:)], ocrust_poly'),nel,nip);                 % [nel,nip]
        isocr_el = sum(isocr_qp,2) >= nip/2;                                    % figure(); plot_meshF(EL2NOD, GCOORD, []); hold on; plot_meshF(EL2NOD, GCOORD, isocr_el,'red')
        if any(isocr_qp(:))                                                        % hold on; plotGEO(GEOn,1:10,'')
            Phases_ip(isocr_qp) = SETTINGS.phaseoc;
        end
        if any(isocr_el)
            Phases_el(isocr_el) = SETTINGS.phaseoc;                          % [1,nel] new phase code for sediment [considers as sediment only elements fully above tracked basement]
        end
    end
    if SP.make
        sedim_poly = [Basement fliplr(GEOn(end-1).coo)];
        issed_qp = reshape(inpoly2([gipx(:) gipy(:)], sedim_poly'),nel,nip);                  % [nel,nip]
        issed_el = sum(issed_qp,2) >= nip/2;
        if any(issed_qp(:))
            Phases_ip(issed_qp) = SETTINGS.phasesed;
            %Ks   = [PHY.K;  SP.K];
            %Cps  = [PHY.Cp; SP.Cp];
            %Hps  = [PHY.Hp; SP.Hp];
        end
        if any(issed_el)
            Phases_el(issed_el) = SETTINGS.phasesed;                          % [1,nel] new phase code for sediment [considers as sediment only elements fully above tracked basement]
        end
    end                                                                    % figure(); plot_phasesF(ELEM2NODE, GCOORD, Phases_el)
    % augment nominal values per phase class 
    Cp = PHY.Cp;                                                           % [nphases,1]
    Cp(SETTINGS.phaseoc)  = PHY.OC.Cp;
    Cp(SETTINGS.phasesed) = SP.Cp; 
    K = PHY.K;                                                             % [nphases,1]
    K(SETTINGS.phaseoc)  = PHY.OC.K;
    K(SETTINGS.phasesed) = SP.K;
    Hp = PHY.Hp;
    Hp(SETTINGS.phaseoc)  = PHY.OC.Hp;                                     % [nphases,1]
    Hp(SETTINGS.phasesed) = SP.Hp;

    EL_IP.Cp = Cp(Phases_ip);                                              % [nel,nip] 
    EL_IP.K  = K(Phases_ip);                                               % [nel,nip]
    EL_IP.Hp = Hp(Phases_ip);                                              % [nel,nip]

    %EL_EL.Cp = Cp(Phases_el);                                             % [nel,1] in case it is needed for archaea
    EL_EL.K  = K(Phases_el);                                               % [nel,1] needed by the conductive heat flux functions 
    %EL_EL.Hp = Hp(Phases_el);                                             % [nel,1] in case it is needed for archaea
end