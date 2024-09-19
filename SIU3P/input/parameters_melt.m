function MPAR = parameters_melt(Rho)
% function MPAR = parameters_melt(Rho)
%  
% +++ purpose +++
% define melting parameters
%

MPAR = [];                                  % melting parameter struct
% igneous body emplacement parameters
MPAR.BDT_temp = 500;                        % [ºC] initiation of BDT. Taken as limit for boundary for dike-sill emplacement type. [450ºC] is for Diabase formation. 500ºX is from Maher et al (2021)
MPAR.axis_feeding_hwidth = 50000.;          % [m] [after Behnand & Grove (2015)] distance from the magma changer x-center that is fully included as feeder of axial magma emplacement
MPAR.axis_feeding_decay_range = 5000.;      % [m] decay range out of the window above
MPAR.Dserp_blocking_threshold = 0.75;       % [0/1] serpentinization degree considered as blocking for upwelling melt. Values >=1 deactivate blocking by serpentinization
MPAR.ErP_blocking_threshold = 1.0E-13;      % IInd invariant of plastic strain rate to be used as blocking value for upwelling melt 
MPAR.sills_a2b = 4.;                        % a/b axes ratio for ellipses emulating sill emplacement
MPAR.emplace_sdev = 500.;                   % [m] standard deviation to randomise the x-location of the dike-sills along the moho polyline. 0 to not randomise
MPAR.Rho_melt = 2300.;                      % melt density [kg/m^3]             
MPAR.temperature_release = 1100.;           % [ºC] temperature of the allocated melting bodies% [ºC] (Tm'= Tm + L/(Cp*rho)) Ros PhD pp.41 (Wangen, 2010) 
MPAR.dT_cutoff_distance = 1000.;            % [m] temperature increment around allocated body decays linearly with distance, falling to 0.0 at this distance
MPAR.dT_diffuse_timeratio = 1.0;            % [-] temperature increments are allowed to diffuse during dt*dT_diffuse_timeratio
MPAR.melt_crystallization_heat = 4e5;       % [J.kg-1] [Sleep & Warren 2014, after Kojitani and Akaogi 1997] we assume crystallization heat is uniformly distributed in time

% PETROLOGICAL CONSTANTS FOR MELTING (I am doing everything for one single component)
% Some constants are not used but are kept for later implementations.
        %nc =1; % number of petrological components in mantle
        %[PETRO] = mantle_composition_marta(nc,Rho);
        
        % Ts0    :: solidus temperature of rock at surface pressure
        % dTs_dP :: change in solidus temperature with increasing pressure
        % dTs_dF :: change in solidus temperature with increasing depletion
        % --> Ts = Ts0 + dTs_dP*P + dTs_dF*F
        
        % Solidus function pool (see JPM 2000)
        Ts(1,:) = [1081  132*1e-9  350]; % [Elena PhD pp.31] (2) fertile peridotite (use Dpl0 to make refractory PD!!) THIS IS HARZBURGITE IN JOERG'S PHD
        iTs  = 1;  % pick solidus function(s) from pool (==1 here)
        Vol0 = 1;  % define volume fraction of each component (MUST SUM UP TO 1)
        Dpl0 = 0;  % define initial depletion of each component
        X0   = 200;% define initial water content for each component
        
       % Rho    = Rho.cnst;
        nc     = length(iTs); % number of melting components
        Ts0    = Ts(iTs,1)';  % [ºC]    solidus function/s: solidus temperature at zero depletion and surface pressure 
        dTs_dP = Ts(iTs,2)';  % [ºC/Pa] solidus function/s: (positive) dependence on pressure
        dTs_dF = Ts(iTs,3)';  % [ºC]    solidus function/s: (positive) dependence on depletion.
        
        %==========================================================================
        % MORE MELTING RELATED PARAMETERS
        %==========================================================================
        % Definition of latent heat QL(enthalpy of fusion)
        % empty ([]) --> QL = (1/3) * T
        % a number   --> latent heat kJ kg^-1; e.g. QL = 400
        % PETRO.QL     = 400;
        % This is going to be T*Increment S/Cp = 550 C % Phipps Morgan,
        % G-cubed, 2001
        QL     = 550;  % [degrees Celsius] (L/Cp) enthalpy of fusion
        
        
        % Maximum upwelling (decompression) interval; code starts to iterate if
        % decompression is larger in a single time step; use values less than 1 km
        dzmax  = 500; % [m]; 0.5 is a good value
        
        % Define adiabatic temperature increase with depth (K/km)
        % (set to zero to use potential temperatures)
        ad_K_km = 0; %0.3;
        
        % Volume fraction below which a component will be removed (i.e. treated as
        % exhausted). Good value: 0.001
        Vol_ex = 0.001;
        
        % Bulk partition coefficient for water between solid and melt
        Dbulk  = 0.01;
        
        % Coefficients for parameterization of wet melting (Katz et al. 2003)
        % dTs(X) = aX * X^bX, with dTs being the decrease in solidus temperature
        % due to the presence of water content X
        aX     = -43;
        bX     = 0.75;
        
        p2z      = 1/(Rho(2)*9.81); % conversion pressure (Pa) --> depth (m)
        z2p      = (Rho(2)*9.81); % conversion depth (m) --> pressure (Pa)
        ad_K_GPa = ad_K_km * p2z;

        
        

        MPAR.Ts = Ts;
        MPAR.iTs = iTs;
        MPAR.Vol0 = Vol0;
        MPAR.Dpl0 = Dpl0;
        MPAR.X0 = X0;
        MPAR.nc = nc;
        MPAR.Ts0 = Ts0;
        MPAR.dTs_dP = dTs_dP;
        MPAR.dTs_dF = dTs_dF;
        MPAR.QL = QL;
        MPAR.aX = aX;
        MPAR.bX = bX;
        MPAR.dzmax = dzmax;
        MPAR.ad_K_km = ad_K_km;
        MPAR.Vol_ex = Vol_ex;
        MPAR.Dbulk = Dbulk;
        MPAR.p2z = p2z; 
        MPAR.z2p = z2p;
        MPAR.ad_K_GPa = ad_K_GPa;
end % function

