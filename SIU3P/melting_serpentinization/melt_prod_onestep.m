function  [Dpl, Temp, dF, dFaxis, dFaxisfac] = melt_prod_onestep(Dpl, Temp, PHY, PRESS_CONT, GCOORD, ELEM2NODE, Phases, GEOn, E2all, nodtopids)
    % +++ purpose +++
    % Estimates depletion and temperature at melting origin,
    % at equilibrium with an updated solidus. That is, updated temperature
    % at areas producing new melting lie at the updated solidus resulting from 
    % the updated depletion.
    %
    % INPUT
    % ---
    % Dpl  :: REAL, DIM(1,nnod) depletion at nodes
    % Temp :: REAL, DIM(1,nnod) temperature at nodes
    % PHY  :: local name for PHY.MELT, struct with at least the parameters:
    %    .z2p      
    %    .Ts0
    %    .dTs_dP
    %    .dTs_dF
    %    .QL  :: REAL [ºC] enthalpy of fusion approximately mapped into degrees
    % PRESS_CONT :: REAL, DIM(1,nnod)
    %               Note: atmospheric pressure is used as
    %                     reference for the solidus calculation. Only for
    %                     consistency with general atmospheric load into the model
    % GCOORD    :: REAL, OPTIONAL, DIM(2,nnod) global coordinates only for debugging
    % ELEM2NODE :: INTEGER, OPTIONAL, DIM(nnodel,nel]) mesh topology
    %
    % OUTPUT
    % --- 
    % Dpl       :: REAL, DIM(1,nnod) output depletion
    % Temp      :: REAL, DIM(1,nnod) output temperature
    % dF        :: REAL, DIM(1,nnod) depletion increment
    % dFaxis    :: REAL, DIM(1,nnod) depletion increment contributing to feed the magma emplaced below the spreading center 
    % dFaxisfac :: REAL, DIM(1,nnod) porportion of dF going into dFcenter [only for used for diagnosis]  
    %
    % Details : The solution is consistent in the sense that the updated
    %   Temperature lies at the new solidus resulting from the updated
    %   depletion, which is consistent with the decrease in temperature
    %   for the nodes that lie above the initial solidus. This direct
    %   one-step solution was solved by Javier GP (unpublished as by 2021). 
    %
    %   An allowance is included to allow than only a proportion of the 
    %   melt is considered for later melt emplacement under the spreading center 
    %  [idea of Leila Mezri after Behnand & Grove (2015)]
    
    % Javier Garcia-Pintado, 2021
    
    atm2pascal = 101325.;
    
    Pref = atm2pascal;                                                 % [Pa]
      
    nnod = size(Temp,2);
    if size(Dpl,2) ~= nnod || size(PRESS_CONT,2) ~= nnod
        error("melt_prod_new_press_cont:: inconsistent dimensions in input variables")
    end
     
    Ts0    = PHY.Ts0;           % unpack physical parameters for melting
    dTs_dP = PHY.dTs_dP;
    dTs_dF = PHY.dTs_dF;
    QL     = PHY.QL;
                      
    Pnew = PRESS_CONT;                                                     % (1,nnod)

    Ts_dry = Ts0 + dTs_dP.*(Pnew-Pref) + dTs_dF.*Dpl;                      % [ºC] (1,nnod) (Elena Phd, Eq. 2.68) commonly pressure gardient has higher influence than previous depletion
    
    a = Ts0 + dTs_dP.*(Pnew-Pref);                                         % (1,nnod)
    b = dTs_dF;                                                            % scalar
    k = 1 / (QL + b);                                                      % scalar
    
    boo = Temp > Ts_dry;
    Temp1 = Temp;                                                          % (1,nnod) updated temperature    
    Temp1(boo) = (a(boo) + b * Dpl(boo) + b*k*Temp(boo)) / (1 + b*k); 
    dF = max(0.0, - k * (Temp1-Temp));                                     % figure(); plot_tF(Ts_dry0, GCOORD, ELEM2NODE, [], "", 0);

    notmantleboo = ~ismember(1:nnod, unique(ELEM2NODE(:, Phases == 1)));
    dF(notmantleboo) = 0.0;                                                 % prevent crustal node depletion - nodes @ moho allowed to deplete

    Dpl = Dpl + dF;                                                         % figure(); plot_tF(dF, GCOORD, ELEM2NODE, [], "", 0); 
    Temp = Temp1;
    
    % updated solidus:
    % Ts_dry1 = Ts0 + dTs_dP.*(Pnew-Pref) + dTs_dF.* Dpl;                   % figure(); plot_tF(Ts_dry0, GCOORD, ELEM2NODE, [], "", 0);
    % figure(); plot_tF(Ts_dry1 - Temp, GCOORD, ELEM2NODE, [], "", 0);
    
    % partition of dF to be emplaced in the ridge axis    
    dFaxisfac = 1.0;
    dFaxis = dF;
    
    hwidth = PHY.axis_feeding_hwidth;
    decran = PHY.axis_feeding_decay_range;
    if hwidth < Inf
        alt = 2;
        switch alt
            case 1 % maximum surface E2all as xcenter criterion
                topelids = find(sum(ismember(ELEM2NODE,nodtopids),1) >= 1); % any triangle touching the  surface  
                EL2NODtop = ELEM2NODE(:,topelids);
                [~,maxintopids] = max(max(E2all(topelids,:),[],2));
                xcenter = GCOORD(1,EL2NODtop(7,maxintopids));                         % hold on; plot(xcenter/1000, GCOORD(2,EL2NODtop(7,maxintopids))/1000,'x','color','cyan','markersize',20)
            case 2 % self-centered magma chamber as xcenter criterion [proposed by Leila]
                dx = 500;                                                              % [m] approximate sampling spacing [discretization] of the horizontal estimation for serpentinization 
                incx = diff(minmax(GCOORD(1,:)));                                       % [m] horizontal size of domain
                nx =  ceil(incx/dx);                                                    % number of horizontal discretization nodes                                      
                dFIx = approxXmarginalIntegrate2D(GCOORD, ELEM2NODE, dF, GEOn, nx);     % figure(); plot(dFIx.seq,dFIx.int)
                [~,xid] = max(dFIx.int);
                xcenter = dFIx.seq(xid);
        end
        dFaxisfac = 1 - variomodels(max(0.0,abs(GCOORD(1,:)-xcenter)-hwidth), ...
                                     0, 1, decran, 'gau');
        dFaxis = dFaxisfac .* dF;                        
        % figure(); plot_tF(dFaxisfac, GCOORD, ELEM2NODE); xline((xcenter-hwidth)/1000,'linewidth',2); xline((xcenter+hwidth)/1000,'linewidth',2);
        % figure(); plot_tF(dFaxis, GCOORD, ELEM2NODE);
    end
    
end % function
