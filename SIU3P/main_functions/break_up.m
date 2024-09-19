function BREAKUP = break_up(BREAKUP, SETTINGS, GEO, ntime)
% BREAKUP = break_up(BREAKUP, GEO, ntime) 
% +++ purpose +++
% checks for break-up conditions
%
% GEO :: defines the geometry of the mesh
% BREAKUP.mid       :: corresponds to the id of the moho in GEO
% BREAKUP.min_thick :: to the crustal thickness for the break-up to be assumed
% BREAKUP.bool :: set to True (by this function) when the break-up condition is met

%--------------------------------------------------------------------------
% Function written by Miguel Andres-Martinez, PhD student at Royal Holloway
% University of London, 15-05-2015.
%
% Javier GP, 2021: modified: integration time is now not modified and
% BREAKUP records the break time [s]
%--------------------------------------------------------------------------

    
    if BREAKUP.bool                                                        % break up condition already met
        return
    end
    
    if BREAKUP.sp_s
        Surface = BREAKUP.Basement;                                        % active surface processes: top basement node coordinates
    else
        Surface = GEO(end-1).coo;
    end
    Moho = GEO(SETTINGS.mohoid).coo;
    
    M2S = interp1(Surface(1,:),Surface(2,:),Moho(1,:),'linear');           % Interpolate surface y-coordinates at the moho node
    S2M = interp1(Moho(1,:),Moho(2,:),Surface(1,:),'linear');              % Interpolate moho y-coordinates at the surface node
    
    M_Sbool = abs(M2S-Moho(2,:))    < SETTINGS.minthick;                   % moho nodes meeting the break-up criterium
    S_Mbool = abs(S2M-Surface(2,:)) < SETTINGS.minthick;                   % surface nodes meeting the break-up criterium
    
    if any(M_Sbool) || any(S_Mbool)                                        % any node meets the break-up criterium
        disp(['BREAK UP CRITERIUM MET: ', ...
            num2str(sum(M_Sbool)+sum(S_Mbool)),' nodes'])
        
        BREAKUP.bool = true;
        BREAKUP.breaktime = ntime;
        %time_int = floor((ntime+2*BREAKUP.break_time)/ ...                % Calculate the new finishing time
        %    BREAKUP.break_time)*BREAKUP.break_time;
    end
    
    %     NOT RUN
    %     figure(); plotGEO(GEOn); hold on; plot(BREAKUP.Basement(1,:)/1000,BREAKUP.Basement(2,:)/1000)
    %     plot(Moho(1,:)/1000,M2S/1000,'.r')
    %     hold on
    %     plot(Surface(1,S_Mbool)/1000,Surface(2,S_Mbool)/1000,'or')
    %     plot(Surface(1,:)/1000,S2M/1000,'.b')
    %     plot(Moho(1,M_Sbool)/1000,Moho(2,M_Sbool)/1000,'ob')
end