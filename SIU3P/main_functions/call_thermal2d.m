function T =  call_thermal2d(ELEM2NODE, GCOORD, Phases, Temp, Cp, Hp, ...
                             K, RHO, Hs, ...                                        % [nel,nip] variables
                             Bct_ind, Bct_val, dt, reorder, nelblo, SETTINGS)
  % purpose: wrapper for Jorg version of thermal solver thermal2d()
  % rift2ridge2D has thermal2d_m() as solver for thermal difussion. This script allows for direct comparison with kinedyn' thermal2d()
  %          to be callable from Miguedyn
  % arguments in this function are exactly as those for rift2ridge2D' thermal2d_m()

  if nargin < 14
    nelblo = 5000;
  end
  
  % nnodel = 6;                                  % forced in thermal2d()
  nip = 6;                                       % forced in thermal2d()

  VAR.Hs_ip = Hs;                                % Heat source. Optional argument to thermal2d()

  MESH.EL2NOD  = ELEM2NODE;                      % possible bubble (7th) node neglected in thermal2d()
  MESH.GCOORD  = GCOORD;                         % possible bubble (7th) node neglected in thermal2d()
  % MESH.PhaseID                                 % not used
                                                 % The following PHYSICS are defined in kinedyn as input parameters
                                                 % at each phase for the thermal solver (except for the density),
                                                 % while kinedyn' thermal2d() considers them as input defined at each integration point
                                                 % TODO: Look upstream, if these should not be better defined at IPs as in Kinedyn.
  PHYSICS.EL_IP.Cp   = repmat(Cp(Phases),1,nip); % [nel,nip] specific heat at constant pressure
  PHYSICS.EL_IP.Dens = RHO;                      % [nel,nip] density
  PHYSICS.EL_IP.Hp   = repmat(Hp(Phases),1,nip); % [nel,nip] radiogenic heat production
  PHYSICS.EL_IP.K    = K;                        % [nel,nip] thermal conductivity 
  
  %NUMSCALE = numerical_eq_scaling('mantle');     % as current kinedyn setup: a soft link to this file is needed within rift2ridge2D
  NUMSCALE = numerical_eq_scaling('SI'); 
  
  SETTINGS.nelblk = nelblo;
  SETTINGS.nnodel = 7;                           
  
  [TBC.ifix,p] = unique(Bct_ind);                 % Cauchy-type Temperature boundary conditions; sort and guarantee uniqueness
  TBC.vfix     = Bct_val(p);
  %TBC.NN_INTEG = 
  %TBC.valHF    = 

  T = thermal2d(Temp, VAR, MESH, SETTINGS, PHYSICS, NUMSCALE, TBC, dt);   
 
end % function call_thermal2d
