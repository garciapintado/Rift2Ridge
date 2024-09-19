function NUMSCALE = numerical_eq_scaling(scaling,U0)

%==========================================================================
% SCALING OF EQUATIONS TO REDUCE NUMERICAL ROUND-OFF PROBLEMS
%
% This is needed if not all input parameters are in SI-units (e.g. length
% in km, time in Myr, viscosity in multiples of 1e19 Pa-s,...).
% Scaling of the viscous flow equations results in a single term "Bscale"
% that multiplies the buoyancy terms in the equations. Thermal diffusivity
% is also scaled so that it matches the time and length units t0 and L0, 
% respectively. 
%
% JH Jul 2014

if nargin==0
    error('Define scaling method as input argument.');
end

switch scaling
    case 'SI' % SI units
        NUMSCALE.L0     = 1;  % m; length unit
        NUMSCALE.t0     = 1;  % s; time unit
        NUMSCALE.Visc0  = 1;  % Pa-s; reference viscosity
        NUMSCALE.Dens0  = 0;  % kg m^-3; g*Dens0 will be subtracted as lithostatic pressure
        NUMSCALE.unit_t = 's';
        NUMSCALE.unit_L = 'm';
    
    case 'mantle'
        Myr             = 1e6*365.25*24*60*60;
        NUMSCALE.L0     = 1000;  % m; length unit
        NUMSCALE.t0     = Myr;   % s; time unit
        NUMSCALE.Visc0  = 1e19;  % Pa-s; reference viscosity
        NUMSCALE.Dens0  = 3300;  % kg m^-3; g*Dens0 will be subtracted as lithostatic pressure
        NUMSCALE.unit_t = 'Myr';
        NUMSCALE.unit_L = 'km';
        
    case 'magma'
        yr              = 365.25*24*60*60;
        NUMSCALE.L0     = 1000; % m; length unit
        NUMSCALE.t0     = yr;   % s; time unit
        NUMSCALE.Visc0  = 1e3;  % Pa-s; reference viscosity
        NUMSCALE.Dens0  = 2800; % kg m^-3; g*Dens0 will be subtracted as lithostatic pressure
        NUMSCALE.unit_t = 'yr';
        NUMSCALE.unit_L = 'km';
        
    case 'syrup'
        NUMSCALE.L0     = 0.001; % m; length unit
        NUMSCALE.t0     = 1;     % s; time unit
        NUMSCALE.Visc0  = 1;     % Pa-s; reference viscosity
        NUMSCALE.Dens0  = 1400;  % kg m^-3; g*Dens0 will be subtracted as lithostatic pressure
        NUMSCALE.unit_t = 's';
        NUMSCALE.unit_L = 'mm';

    case 'Spieg2016'
        NUMSCALE.L0     = 120000;         % m; length unit
        NUMSCALE.t0     = NUMSCALE.L0/U0; % s; time unit
        NUMSCALE.Visc0  = 1e22;  % Pa-s; reference viscosity
        NUMSCALE.Dens0  = 2700;  % kg m^-3; g*Dens0 will be subtracted as lithostatic pressure
        NUMSCALE.unit_t = 'L0/U0';
        NUMSCALE.unit_L = 'x 120 km';
        
    case 'Kaus2010'
        NUMSCALE.L0     = 10000; % m; length unit
        NUMSCALE.t0     = 1e15;  % s; time unit
        NUMSCALE.Visc0  = 1e20;  % Pa-s; reference viscosity
        NUMSCALE.Dens0  = 0;     % kg m^-3; g*Dens0 will be subtracted as lithostatic pressure
        NUMSCALE.unit_t = 'x 1e15 s';
        NUMSCALE.unit_L = 'x 10 km';

    otherwise
        error(' Unknown scaling method. See function "numerical_eq_scaling" for available methods.');
end

%==========================================================================
% Do not edit the next lines...

NUMSCALE.U0     = NUMSCALE.L0 / NUMSCALE.t0;
NUMSCALE.unit_U = [NUMSCALE.unit_L '/' NUMSCALE.unit_t];
    % velocity unit
NUMSCALE.Kappa0 = NUMSCALE.L0^2 / NUMSCALE.t0;
    % scaling for diffusivity
NUMSCALE.P0     = NUMSCALE.Visc0 / NUMSCALE.t0;
    % scaling for pressure and shear modulus

% Calculate buoyancy scaling factor
% This factor multiplies the gravitaional force in the element assembly
NUMSCALE.Bscale = NUMSCALE.L0^2 / (NUMSCALE.Visc0 * NUMSCALE.U0);

end % END OF FUNCTION numerical_eq_scaling