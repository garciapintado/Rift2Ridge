
function [Crust_thickness, TRACKP_melt, dT_sensible, ign_body, tp_melt, dE_latent, dT_ele_affected] = emplacementDikeSill(...
          Crust_thickness, TRACKP_melt, Temp, area_melt, GCOORD, ELEM2NODE, GEO, Phases, ext_rate, dt, ...
          E2all, ErP, dF, Dserp, PHY, SETTINGS, istep, Bct_ind, EL_IP, nelblk, BDTat, ocrust)

    %This function calculates the place where the dike or the sill is located.
    %DIKE: It will be located below the Moho and centered where the maximum
    %strain rate takes place.
    %SILL: It will be located below the Moho over a distance equal to the width
    %of the melting area.
    % ext_rate: [m/s] full spreading velocity
    %
    % SETTINGS :: struct with components
    %         .heat_release :: LOGICAL      
    %         .igneous_body :: STRING in {"dike","sill"}
    %
    %
    % GEO :: STRUCT, OPTIONAL for plotting purposes and debugging
    % PHY :: STRUCT [global PHY.MELT]. 
    %    .T_melt   :: for heat release, currently not used
    %    .Rho_melt :: for heat release, currently not used
    %
    % OUTPUT
    %
    % ign_body :: STRUCT with fields
    %         .x :: REAL
    %         .y :: REAL
    %         .i :: INTEGER
    % tp_melt :: LOGICAL true if melting has ocurred in this time step
    %
    % Author: - Elena Ros: 1st version
    %         - Javier GP, 2021-09
    %           modified to account for 
    %           a] brittle-ductile transition [BDT] and dual dike-sill emplacement
    %           d] heat-release for multiple igneous bodies
    
    if nargin ~= 22
        error("emplacementDikeSill:: all 22 arguments are needed");
    end
    
    dT_ele_affected = false(1,size(ELEM2NODE,2));
    dT_sensible = [];
    dE_latent = [];
    
    if area_melt == 0.
        tp_melt = false;
        Crust_thickness(istep) = 0.;
        ign_body = [];
    else
        tp_melt = true;
        dz_crust = area_melt/(ext_rate*dt);                                    % [m] estimated crustal thickness, only as diagnostic (assumes all melt at this time step accumulates to create new crust)

        %===================== 1st part ====================================
        % Calculate the column of dike where maximum strain rate takes place
        %===================================================================
        Crust_thickness(istep) = dz_crust;                                 % [m]
        %nnod6 = max(ELEM2NODE(1:6,:),[],'all');
        [dike, sills, TP_xmelt, TP_ymelt, x_shallow, y_shallow, radius_tpmelt_all, intrusive_type] = ...
            placeDikeSill_belowserp(GCOORD, ELEM2NODE, GEO, Phases, area_melt, ext_rate, dt, ...
            E2all, ErP, Dserp, Temp, PHY, SETTINGS, BDTat, dF, ocrust);

        %TRACKP_melt{istep} = [TP_xmelt(:)'; TP_ymelt(:)']; 
        TRACKP_melt = [TRACKP_melt,[TP_xmelt;...
                                    TP_ymelt;...
                                    radius_tpmelt_all;...
                                    repmat(istep,[1 length(TP_xmelt)]);...
                                    intrusive_type]];
        ign_body.dike = dike;
        ign_body.sills = sills;
        
        %===================== 2nd part ====================================
        % heat release 
        %===================================================================
        if SETTINGS.heat_release
            % figure(); plotGEO(GEO,1:10,'');
            % Tisoc = getIsolines(GCOORD,Temp,[300, 325, 350],true);
            % hold on; plot(Tisoc(2).coo(1,:)/1000, Tisoc(2).coo(2,:)/1000); % 325ºC
            % hold on; plot(ign_body.dike.coo(1,:)/1000, ign_body.dike.coo(2,:)/1000,'o')
            % for i=1:length(ign_body.sills)
            %     hold on; plot(ign_body.sills(i).coo(1,:)/1000, ign_body.sills(i).coo(2,:)/1000,'o','color',[.4 .4 .0])
            % end
            [dT_sensible, dE_latent, dT_ele_affected] = heat_release(GCOORD, ELEM2NODE, Temp, ign_body, dt, Bct_ind, EL_IP, PHY, nelblk, GEO);    % [1,nnod7] fields of updated temperature and increment ocurred by heat of the melting body emplacement
            % [Temp] = heat_release_dike(GCOORD,Temp,Geo_id,dz_crust,K,...
            %                            Cp,ext_rate,dt, PHY.T_melt, PHY.Rho_melt,x_shallow,y_shallow);
            % Temp = heat_release_dike_errorfunction(GCOORD,Temp,dz_crust,K,...
            %         Cp,ext_rate,dt, PHY.T_melt, PHY.Rho_melt, x_shallow,y_shallow);
        end
    end % area_melt > 0.   
end % function 
