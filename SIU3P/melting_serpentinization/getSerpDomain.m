function serpelboo = getSerpDomain(GCOORD, ELEM2NODE, Phases, inH, Temp, ErP, PHY)
    % function serpelboo = getSerpDomain(GCOORD, ELEM2NODE, Phases, inH, Temp, ErP, PHY)
    % 
    % +++ purpose +++
    % Get a LOGICAL vector indicating elements subject to hydrothermal cooling.
    % Elements are considered, as a whole, as subject/not-subject to hydrothermal domain
    %
    % INPUT
    % GCOORD    ::
    % ELEM2NODE ::
    % Phases [1,nel]
    % inH       :: [1,nel] hydrothermal criterion. Only elements within this are allowed to be serpentinised
    % Temp      :: [nel,nip] temperatue, for temperature criterion
    % ErP       :: [nel,nip] II invariant of the plastic strain rate, as additional criterion 
    % PHY       :: STRUCT, corresponding to the global PHY
    %
    % Details:
    % inH: hydrothermal domain: as minimum, it includes includes a) geothermal gradient criterion, b) thickness criterion
    
    % Author : Javier García-Pintado, MARUM, 2020
 
    % plastic strainrate criterion [not necessarily used to obtain the hydrothermal domain]
    if isempty(PHY.E2fac)
        ErPmin = PHY.ErPmin;                                               % ErPmin = env['PHY']['SERP'].item()['ErPmin'].item()
    else
        if PHY.E2fac < 0. || PHY.E2fac > 1.
            error("PHY.E2fac should be in [0,1]")
        end
        ErPmin = PHY.E2fac * (max(ErP,[],'all') -  min(ErP,[],'all')) + min(ErP,[],'all');
    end         
    inE = any(ErP' >= ErPmin);                                             % LOGICAL [1,nel]  inE = np.any(env['ErP'] >= ErPmin,1)
                                                                           % inH = env['hydroelboo'] == 1; inP = env['Phases'] == 1
    TempV = Temp(ELEM2NODE(1:3,:));                                        % TempV = env['Temp'][env['ELEM2NODE'][0:3,:]-1]
    inT = any(TempV >= PHY.Tmin & TempV <= PHY.Tmax);                      % Tmin = env['PHY']['SERP'].item()['Tmin'].item()
                                                                           % Tmax = env['PHY']['SERP'].item()['Tmax'].item()
    serpelboo = inH & inE & inT & ismember(Phases, 1);                     % inT = np.any((TempV >= Tmin) & (TempV <= Tmax),0)
                                                                           % np.sum(inH & inE & inT & inP)
end % function
