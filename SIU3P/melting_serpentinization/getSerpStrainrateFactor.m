function srfac = getSerpStrainrateFactor(ELEM2NODE, Phases, inH, ErP, PHY, eids)
    % function srfac = getSerpStrainrateFactor(GCOORD, ELEM2NODE, Phases, inH, ErP, PHY, eids)
    % 
    % +++ purpose +++
    % Get a REAL [1,nnod] vector with values between 0 and 1, which will be eventually used as a multiplier of the kinetics of the serpentinization reaction 
    % Both, the phase (only allowed for mantle) anf an input hydrothermal mask are used to set srfac=0.0 when any of these masks if false. 
    %
    % INPUT
    % GCOORD    ::
    % ELEM2NODE ::
    % Phases [1,nel]
    % inH       :: [1,nel] hydrothermal criterion. Only elements within this are allowed to be serpentinised
    % Temp      :: [nel,nip] temperature, for temperature criterion
    % ErP       :: [nel,nip] II invariant of the plastic strain rate, as additional criterion 
    % PHY       :: STRUCT, corresponding to the global PHY
    %
    % Details:
    % inH: hydrothermal domain: as minimum, it includes includes a geothermal gradient criterion
    
    % Author : Javier García-Pintado, MARUM, 2021
 
    % plastic strainrate factor in (0,1) masked by hydroelboo (i.e. 0 where hydroelboo=false) and Phase

    ErPnod = ipval2nodval(ELEM2NODE, ErP, eids, true)'; % [1,nnod7]
    nnod = size(ErPnod,2);

    if isempty(PHY.E2fac)
        ErPnom = PHY.ErPnom;                                               % fixed-input nominal plastic strainrate to allow for full serpentinization kinetics
    else
        if PHY.E2fac < 0. || PHY.E2fac > 1.
            error("PHY.E2fac should be in [0,1]")
        end
        ErPnom = PHY.E2fac * (max(ErPnod,[],'all') -  min(ErPnod,[],'all')) + min(ErPnod,[],'all');
    end
    srfac = variomodels(ErPnod, 0, 1, ErPnom, PHY.sfun);                   % [nnod,1] \in (0,1)
    
    inHn = false(1,nnod);                                                  % hydrothermal domain
    inHn(unique(ELEM2NODE(:,inH))) = true;
    inPh = false(1,nnod);                                                  % mantle domain
    inPh(unique(ELEM2NODE(:,Phases==1))) = true;

    srfac(~inHn | ~inPh) = 0.0;

                                                                           % inH = env['hydroelboo'] == 1; inP = env['Phases'] == 1
    %TempV = Temp(ELEM2NODE(1:3,:));                                        % TempV = env['Temp'][env['ELEM2NODE'][0:3,:]-1]
    %inT = any(TempV >= PHY.Tmin & TempV <= PHY.Tmax);                      % Tmin = env['PHY']['SERP'].item()['Tmin'].item()
                                                                           % Tmax = env['PHY']['SERP'].item()['Tmax'].item()
    %serpelboo = inH & inE & inT & ismember(Phases, 1);                     % inT = np.any((TempV >= Tmin) & (TempV <= Tmax),0)
                                                                            % np.sum(inH & inE & inT & inP)    
end % function
