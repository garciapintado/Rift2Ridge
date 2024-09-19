function [SETTINGS, PHY, SOLVER, SP, PLOT] = param_defaults(dsninput, COMPSET)
  % function [SETTINGS, PHY, SOLVER, SP, PLOT] = param_defaults(dsninput, COMPSET)
  % +++ purpose +++
  % rift2ridge2D default model input
  % 
  % INPUT
  % dsninput :: path where to fin specific external input (e.g. thermodynamic tables)
  %
  % Authors: Javier GP, Marta PG. MARUM, 2020

  km      = 1000.;                                                           % [m]
  year    = 365.25*24*60*60;                                                 % [s]
  ma      = 1.e6*year;                                                       % [s]
  bar2pascal = 100000.;                                                      % [Pa/bar]
  atm2pascal = 101325.;                                                      % [Pa/atm]
  %---------------
  %% SETTINGS
  %---------------
  SETTINGS = [];
  SETTINGS.COMPSET = COMPSET;
  SETTINGS.eids = "645";                         % [rift2ridge2D counterclockwise edge node convention]
  %SETTINGS.MODEL = "rift2ridge2D";

  %% saving options
  SETTINGS.save_by = 10;                        % number of steps between saves. TODO: change to time unit [e.g. 10yr]
  SETTINGS.save_part = true;                    % if true, continuous save of:
  SETTINGS.sav_var_part = {'GCOORD','ELEM2NODE','Point_id','dt','Vel', ...
      'istep','ntime','ma','km','ISOCHRONS','Basement', ...
      'E2all','Topography','GEOn'};

  % rift2ridge2D - kinedyn local node transformation
  SETTINGS.ki2mi_lnodes = [1,2,3,5,6,4,7]; % sf_dsf_tri367() output local nodes to milamin
  SETTINGS.mi2ki_lnodes = [1,2,3,6,4,5,7]; % shape_drv_fun() output local nodes to kinedyn

  % plotting options
  PLOT = [];
  PLOT.make = false;
  PLOT.func = "plot_b_sb";            % name of plotting function inside the mechanical
  PLOT.FigNo = 1;                     % plot figure number to start off

  %---------------
  %% SOLVER
  %---------------
  SOLVER = [];
  SOLVER.mu_min = 1.e18;                                                     % Lower viscosity limit
  SOLVER.mu_max = 1.e24;                                                     % Upper viscosity limit
  SOLVER.nip_stress  = 6;                                                    % Nº integration points for viscosity and strains
  SOLVER.nip_elast   = 6;
  SOLVER.mu_loc      = "ips";                                                % Location of the viscosity calculation
  SOLVER.fpen     = 1.e6*SOLVER.mu_max;                                      % Penalty factor
  SOLVER.iteration_type = "MILAMIN-penalty";                                 % Iterative solver: 'MILAMIN-penalty' or 'Newton'
  SOLVER.nitmax   = 50;                                                      % Maximum number of iterations
  SOLVER.mat_method  = "block";                                              % Method to build matrices: 'element' or 'block'
  SOLVER.symmetric_tau = false;                                              % force normal deviatoric stresses in 2D to be tau_x + tau_y = 0. (as should theoretically be) by simple averaging
  if ismember(computer, ["MACI64","GLNXA64"])
      SOLVER.nelblo     = 5000;                                              % Nº elements per block if 'block' option is chosen above
  else                                                                       % (also passed to other functions)
      SOLVER.nelblo     = 5000;
  end

  % iterations
  SOLVER.tol_Exit = 1.e-4;                                                   % Tolerance for leaving iterations
  SOLVER.tol_E_rel= 1.e-7;                                                   % " " " relative to the initial residual
  SOLVER.tol_Newt = 1.e-3;                                                   % Tolerance for switching from Picard to Newton
  SOLVER.tol_N_rel= 1.e-5;                                                   % " " " relative to the initial residual
  SOLVER.beta = 1.00;                                                        % Weighting factor in u(iter+1) = u(iter) + beta*du
  SOLVER.method_UP = "Coupled";                                              % Solving method:: 'Penalty': penalty-picard iteration | 'Couple': a first 'Picard' iteration, then 'Newton'
  SOLVER.method = "Picard";                                                  % Iteration method: always start with Picard

  %---------------
  %% PHYSICS
  %---------------
  PHY = [];
  % constants
  PHY.R = 8.314463;                                                         % [J / (mol K)] Gas constant
  PHY.G = [0; -9.81];                                                       % Gravity

  % 
  PHY.dt  = 10000*year;                                                      % [s] reference time step
  PHY.time_int = 50*ma;                                                      % [s] Time interval for the model

  % Break up conditions
  % -------------------
  SETTINGS.BRKU.mohoid = 3;                                                 % Moho id
  SETTINGS.BRKU.minthick  = 150;                                            % [m] Crustal thickness for which the model is considered to be on break-up

  %--------------------------------------------------------------------------
  %% THERMODYNAMIC PROPERTIES OF RHEOLOGICAL UNITS
  %--------------------------------------------------------------------------
  PHY.K    = [3.3;      2.5;            2.1];                          % [W.m-1.K-1] thermal conductivity [mantle, LC, UC]
  PHY.Cp   = [1200.;  1200.;          1200.];                          % $\C_p$ heat capacity             [  ",     ",  "]
  PHY.Hp   = [0.;    0.2e-6;         1.3e-6];                          % $H_r$ radiogenic heat production      ./[1000; 1000; 1000]./Rho;
  PHY.Ce   = [3.;       2.4;            2.4] * 1.e-5;                  % $\alpha_T$ thermal expansion coefficient  [mantle, LC, UC]
  PHY.Dplf = 0.04;                                                     % $\beta$ parameter for change of rho due to melting (value for garnet, below ~60 km at mid-ocean ridges)

  %--------------------------------------------------------------------------
  %%      - Rheologies
  %          * Olivine preexponential factors -
  %--------------------------------------------------------------------------
  PHY.RHEOL = [];
    % WET DIFFUSION AND DISLOCATION PARAMETERS FOR OLIVINE
    Ndif_w        = 1;                      % Power law exponent
    Adif_w        = 1.e6*(1.e6)^(-Ndif_w);  % Pa^-1*s^-1
    Qdif_w        = 335.e3;                 % Activation energy
    Vdif_w        = 4.e-6;
    
    Ndis_w        = 3.5;                    % Power law exponent
    Adis_w        = 90.*(1.e6)^(-Ndis_w);   % Pa^-n*s^-1
    Qdis_w        = 480.e3;                 % Activation energy
    Vdis_w        = 10.e-6;
    
    % DRY DIFFUSION AND DISLOCATION PARAMETERS
    Ndif_d        = 1;                      % Power law exponent
    Adif_d        = 1.5e9*(1.e6)^(-Ndif_d); % Pa^-1*s^-1
    Qdif_d        = 375.e3;                 % Activation energy
    Vdif_d        = 6.e-6;
    
    Ndis_d        = 3.5;                    % Power law exponent
    Adis_d        = 1.1e5*(1.e6)^(-Ndis_d); % Pa^-n*s^-1
    Qdis_d        = 530.e3;                 % Activation energy
    Vdis_d        = 13.e-6;
    
    m = 3;
    f = 500.;                               % fugacity ppm ; H/Si
    d = 6000.;                              % 6 mm in microns

  PHY.RHEOL.Bdis_w = Adis_w*f^(1.2);          % used by load_rheol()
  PHY.RHEOL.Bdif_w = Adif_w.*d.^(-m)*f^1;     % "   " 
  PHY.RHEOL.Bdis_d = Adis_d;                  % "   "
  PHY.RHEOL.Bdif_d = Adif_d.*d.^(-m);         % "   "
  PHY.RHEOL.type   = "hk_lc_gr";              % "   "             {"hk_lc_gr": mafic granulite lower crust (e.g. Gakkel ridge) | "hk_lc_an": wet anorthite lower crust  
  PHY.RHEOL.Grain  = 1.;                      %                   | "hk_lc_wq": wet quartzite lower crust (e.g. angola margin) | "all_peridotite": emmulate all mantle with 3 layers}
  PHY.RHEOL.Burger = 1.;                      %
  PHY.Shearm       = [74.e9;           40.e9;           36.e9];
                                                  
  %--------------------------------------------------------------------------
  %%      - Viscosity limits and non-linear iterations
  %--------------------------------------------------------------------------
  PHY.RHEOL.mu_min = 1.e18;                                                       % Lower viscosity limit
  PHY.RHEOL.mu_max = 1.e24;                                                       % Upper viscosity limit
  PHY.RHEOL.max_it = Inf;                                                         % maximum number of iterations
  PHY.RHEOL.var_f  = "rheol_var_Dpl_ip";                                          % rheology changes inside a phase (need in case >1 column in A, N and Q
                                                                                  % 1 for Hirth because is included inside A
  SETTINGS.is_elastic = true;                                                     % enable/disable ELASTICITY for kinematic solver
  PHY.sh_factor  = [1.; 1.; 1.];                                                  % shear heating multiplier at each layer. 0.: no heating | 1.: full heating

  %--------------------------------------------------------------------------
  %%      - Yielding criteria and strain weakening parameters
  %--------------------------------------------------------------------------
  PHY.PLASTICITY.type   = "moresi";                   % "maxwell", "moresi", "moresi_no_fc" (no forecasting of the stresses), "harmonic" or "no" plasticity
  PHY.PLASTICITY.yc     = "druckerprager";            % yield criterion: "vonmisses", "vonmissesP" or "druckerprager"
  PHY.PLASTICITY.Phases = [1 2 3];                    % phases where plasticity is active

  PHY.Phi = 30. * pi/180.;                            % Friction angle [0-30]                        note: kinedyn has this defined for each phase
  PHY.Cohesion = 10.E06;                              % [Pa] Make sure this value is equal to SS.C(1)            "               "

  % STRAIN SOFTENING
  % ================
  % Smoothing of the finite strain before strain softening parameter
  % calculation to choose between:
  %   - 'elem_avrg'   - average inside an element, constant strain along the
  %                     element
  %   - 'linear       - linear variation of the strain rate inside the
  %                     element
  %   - 'linear_avrg' - linear variation of the strain rate inside the
  %                     element, continue strain rate across the elements by
  %                     averaging at the nodes
  %   - 'no'          - no smoothing
  PHY.SS = [];
  PHY.SS.e_smooth = "no";                             % used by mechanical2dN()
  PHY.SS.rI2 = "linear_N";                            % strain remeshing ["triscatteredinterp"]

  % Strain softening option to choose among 'allI' for softening all different
  % deformation mechanisms with the full accummulated strain or, 'partI' for
  % softening each mechanism with its correspondent deformation. Note that in
  % 'partI' the partial I2 are calculated without using the change of
  % gradient of deformation.
  PHY.SS.opt = "partI";                                  % used by strain_softening_SST_rand() and mechanical2d_m()
 
  % Friction angle reduction
  % ------------------------
  % Layers where friction angle strain softening is applied, 1 is the lower
  % layer and 4 is the upper layer. Set and empty vector [] for no softening
  PHY.SS.Phases_phi = [1 2 3];
  PHY.SS.Phi = [PHY.Phi; 15.*pi/180.];                       % Range of the friction angle for strain softening
  PHY.SS.I2_phi = [0. 1.];                                   % Range of accumulated II invariant of the strain where strain softening is applied
  
  % Cohesion reduction
  % ------------------
  % Layers where cohesion strain softening is applied, 1 is the lower
  % layer and 4 is the upper layer. Set and empty vector [] for no softening
  PHY.SS.Phases_c = [];                                     % TODO: where is this used?
  PHY.SS.C = [PHY.Cohesion; 10.E06];                        % [Pa] Range of the cohesion for the strain softening
  PHY.SS.I2_c = [0. 1.];                                    % Range of accumulated II invariant of strain where strain softening is applied
                                                
  % Viscosity preexponential factor dislocation
  % --------------------------------------------
  % Layers where viscous softening is applied, 1: lower layer; 4: upper layer. Set empty vector [] for no softening
  PHY.SS.Phases_pef_dis = [1 2 3];
  PHY.SS.Pef_dis    = [1 30];             % range of preexponential factor for strain softening
  PHY.SS.I2_pef_dis = [0 1];              % range of accumulated II invariant of strain where strain softening is applied
  
  % Viscosity preexponential factor diffusion
  % ------------------------------------------
  % Layers where viscous softening is applied, 1: lower layer; 4: upper layer. Set empty vector [] for no softening
  PHY.SS.Phases_pef_dif = [1 2 3];
  PHY.SS.Pef_dif = [1 30];                % range of preexponential factor for strain softening
  PHY.SS.I2_pef_dif = [0 1];              % range of accumulated II invariant of strain where strain softening applied
 
  % Viscous softening dependency of the temperature
  % -----------------------------------------------
  % ss_dep_t to: 0 if softening dependence on temperature is not needed, 1 otherwise
  PHY.SS.ss_dep_t = true;
  PHY.SS.dif_low_lim = 800;                     % Diffusion in ss_dep_temp()
  PHY.SS.dif_upp_lim = 1200;
  PHY.SS.dis_low_lim = 800;                     % Dislocation in ss_dep_temp()
  PHY.SS.dis_upp_lim = 1200;
  
  % Run to check that the softening factor for the pef is correct:
  % pef_softening_function_check
    
  %--------------------------------------------------------------------------
  %%      - Random damage
  %--------------------------------------------------------------------------
  PHY.SS.rand_s = true;                                                          % At the moment it only supports random damage at the friction angle
  % Random seeding: same random parameters for every rand if it has a value,
  % otherwise if empty or non-existing randomisation is shuffled. IMPORTANT:
  % different seeds imply different integers, if for example 2 seeds are
  % given in two different models and their values are 0.3 and 0.9 rng
  % function considers this as the same seed (1). Use instead for example 3 and 9.
  PHY.SS.rng = 965;
  PHY.SS.Phi_var = 4. * pi/180.;        % Variation from Phi in the randomization SS.Phi_var/2 above and SS.Phi_var/2 below

  % Parameters for the random damage weak seed (SETTINGS.ini_deformation = 6)
  PHY.RDWS.function = "gaussian_band_rand_damage";
  PHY.RDWS.sigma    = 100. *km;                                                    % [m] wave length
  PHY.RDWS.factor   = [2. 5.];                                                     % Factor that multiplies the randomization (amplitude)

  %%      - Deformation seed
  %--------------------------------------------------------------------------
  % ini_deformation is a switch to define different initial topographies for the interfaces:
  %       0 - Flat topography
  %       1 - Sinusoidal [negative] topography for Moho only
  %       2 - Sinusoidal [negative] topography for Moho and resulting topography for surface
  %       3 - Flat topography
  %           Increment of temperature on the model center with vertical and horizontal gaussian decay
  %       4 - Flat topography
  %           Boundary conditions for the temperature elevated for the center of the model (necking of the lithosphere)
  %       5 - Flat topography. Damage weak seed in terms of strain rate
  %       6 - Flat topography. Random damage seed. To use this option, activate random damage (PHY.SS.rand_s = true) 
  %           and select the RDWS function parameters
  SETTINGS.ini_deformation = 3;

  PHY.Dens = [3360.; 2850.; 2700.];                                               % [kg m^-3] Initial density at each phase

  %% BOUNDARY CONDITIONS
  %%      - Temperature                                               
  %--------------------------------------------------------------------------
  PHY.temp_bc_depth = -120.*km;                                                   % 1) BC temp: depth [m]
  PHY.temp_bc = 1300.;                                                            % 2) BC temp: value [ºC]
  PHY.Ttop = 0.;                                                                  % 3) BC temp: value at surface [ºC]

  PHY.temp_bc_depth_ini = -120.*km;                                               % depth BC for the initial conditions
  PHY.temp_bc_ini = 1300.;                                                        % temperature BC for the initial conditions
  PHY.load_T = false;                                                             % true: loads initial temperatures from regular grid (calculated with SteadyTempGUI.m and TransientGUI.m thermal solvers)
  PHY.Tini_dir = '';                                                              % file including the regular mesh and the temperature if PHY.load_T = true
  PHY.Tini_load = '';                                                             % name of function to load temperature

  %% BOUNDARY CONDITIONS
  % FREE SURFACE SETUP
  SETTINGS.top_surface = "free";                   % "free"| "fix": top model surface as a x-fixed=0 surface (free-slip)
  SETTINGS.fs_alpha    = 2./3.;                    % Miguel TODO: evaluate alpha. 0 is more accurate but less stable. A value of 0.5 is less accurate but more stable. ~2/3 seems optimal
  SETTINGS.fs_beta     = 0.;                       % Parameter to control the influence of the free surface x-term in the calculation of the average of h for the next step for the free surface.

  %% INITIAL CONDITIONS
  %%      - Geometry                                                  
  %--------------------------------------------------------------------------
  PHY.Lx_model = 400.;                            % [km] initial domain full width
  PHY.xmin =  -PHY.Lx_model/2 * km;
  PHY.xmax =   PHY.Lx_model/2 * km;
  PHY.ymin1 = -150.*km;                 % layer 1 [mantle]      - bottom [== bottom of domain]
  PHY.ymax1 =  -35.*km;                 % layer 1 [mantle]      - top [== layer.2 bottom == zmoho] 
  PHY.ymax2 =  -20.*km;                 % layer 2 [lower crust] - top [== layer.3 bottom]
  PHY.ymax3 =   0.*km;                  % layer 3 [upper crust] - top
  PHY.zlab  = -125.*km;                 % lithosphere-asthenosphere-boundary - only to be tracked and used for initialization of depletion

  PHY.sigma_moho = 20.*km;         % IC anomaly length:  initial moho anomaly [for SETTINGS.ini_deformation \in {1,2}]
  PHY.topo_moho  = -4.*km;         % IC anomaly ampli.:     "      "      "      "      "      "
  
  % - Temperature anomaly as week seed initializer of the deformation with spatial Gaussian decay (ini_deformation = 3 and 4)
  PHY.t_seed_max = 100.;                                                        % (ºC) maximum temperature anomaly
  PHY.t_seed_x = mean([PHY.xmin,PHY.xmax]);                                     % (m)  maximum x (horizontal) location
  PHY.t_seed_y = -30.*km;                                                       % (m)  maximum y (depth) location 
  PHY.t_sigma_x = 10.*km;                                                       % (m) sdev for the x Gaussian decay
  PHY.t_sigma_y = 20.*km;                                                       % (m) sdev for the y Gaussian decay
  
  % Parameters for the damage weak seed (ini_deformation = 5)
  PHY.WS = [];
  PHY.WS.damage = 1;                                                             % In terms of strain rate
  PHY.WS.size   = 2.5*km;                                                        % Radius of the ares affected by the damage [m]
  PHY.WS.coord  = [0*km -42.5*km];                                               % [x,y] center of damage area                                              
  PHY.WS.time   = Inf;                                                           % time step to stop
  PHY.WS.shape  = "circle";                                                      % used by damage_seed()   
  PHY.WS.res    = 200.;                                                          % [m]

  % MESH RESOLUTION
  % true: randomize increment for crustal resolution varying from 50 to -50 m to generate different meshes
  % JGP: I've taken randomization of (r_inc) out of here. If desired, this should included in readMC() Monte Carlo function
  % note: it is possibly advantageous to use a more flexible resolution input as in kinedyn 
  resmul = 1.;
  SETTINGS.hres = resmul * [5000. 4000. 2000. 1000.];                    % [m] 'soft' resolution for horizontal interfaces from bottom to top
  %                        | mantle |  lc  |  uc  |                       mantle includes dry and wet olivine mantle 
  SETTINGS.geomf = "make_geometry3l";                                   % name of function to make geometry. Use this function as template to check returned values
  
  % Switch for the different functions used to deal with too thin layers and intersections of interfaces::
    %  'none'           no rupture of interfaces - use for option for compressional settings
    %  'dis_layers'     for real discontinuous layers (not working)
    %  'layer_collapse' for a simpler concept for real discontinuous layers (working)
    %  'lineintersect' for layers no thinner than shift
  SETTINGS.intersect_s = "layer_collapse";
  % Smaller distance between two interfaces. When the distance is smaller
    % either 'dis_layers' removes the layer in this area or 'lineintersect'
    % fixes the thickness of the layer at a distance shift
    % 
    % .tolerance: smaller values reduce numerical diffusion by resampling but leads to more nodes and computation needs

  SETTINGS.layshifts = [1000., 200.];                                            % [m] pre-/post-breakup clipping distance between GEO quasi-horizontal interfaces
  SETTINGS.tolerancef = [.2 .1 .1 .1];                                           % [%] clipping distance for each horizontal layer (later mapped to vertical sides) = tolerancef * GEOres  [m] in simplifyGEO() within a GEOMETRY interface
  SETTINGS.temp_remesh_type = "shpf_cuadr";                                      % Type of temperature remeshing \in {'shpf_cuadr','interp2d_cubic','triscatter'}

  %--------------------------------------------------------------------------
  %%      - Mechanical                                                
  %--------------------------------------------------------------------------
  PHY.bc_t = "winkler";                                                          % mechanical boundary condition type: {'pureshear','winkler','winkler+LIVD'}, where XIVD stands for Lateral Infinite Virtual Domain
  PHY.ext_rate = 10. / 1000. / year;                                             % [m/s] = [mm/yr * m/mm * yr/s] full extension rate [km/My == mm/yr]

  %--------------------------------------------------------------------------
  %% SCALING                                                          
  %--------------------------------------------------------------------------
  % Scaling type to choose among scaling the 'residual', 'full' or 
  % 'before_back_solve' ('before_back_solve' only works with Newton solver)
  SETTINGS.NUMSCALE = [];
  SETTINGS.NUMSCALE.t = "residual";  

  %% - mechanical BDT location
  PHY.BDT_Tcent = 450.; % [ºC] center : temperature criterion to constraint the location of the mechanical BDT. The likelihood of its emplacement is 
  PHY.BDT_Trang = 150.; % [ºC] range  : assumed to be maximum in the range [Tmean-Tsdev,Tmean+Tsdev], and then decays with a Gaussian variogram model with range=Trange 

  %% - MELTING  -  note: this version contains diking as diagnostic. It does not include sills nor dikes have any feedback into the mesh                                               
  SETTINGS.MELT = [];                              % general melting switches - melting physical parameters go into PHY.MELT
  SETTINGS.MELT.make = false;
  if contains(COMPSET,"melt")
      SETTINGS.MELT.make = true;                    % true: compute & allocate melting
      SETTINGS.MELT.emplace_at = "E2allmax_moho"; % {"E2allmax_moho":maximum deviatoric strainrate at moho|"BDT_nn_moho":nearest neighbour in moho to BDT}
      SETTINGS.MELT.emplace_idwmagma = true;      % weights horizontal distance to the maximum of the dF value (magma chamber at the corresponding timestep)
      SETTINGS.MELT.heat_release = false;                  % release of heat from allocated melting. WARNING: Elena code does not seem to work fine in this part
      SETTINGS.MELT.sillmat = [1,3];                       % dimensions [nrow,ncol] for the array of gabbro sills emplaced at each timestep
      SETTINGS.MELT.dila_type = "";                        % {"":no dilation,"linear":linear decrease of dilation with depth, "deep_ml": deep melt lens} 
      SETTINGS.MELT.eids = SETTINGS.eids;
  end

  %% - OCEANIC CRUST DOMAIN ESTIMATION
  SETTINGS.OC.make = true;                 % consider serpentinization and emplaced melting to get a subdomain considered as oceanic crust with OC.Dens density
  PHY.OC.serp_dx =  500.;                  % [m] horizontal sampling spacing for serpentinization integral over vertical sections
  PHY.OC.melt_dx = 5000.;                  % [m] width of blocks for piece-wise integral of emplaced melting bodies - the wider, the smoother the estimation
  % properties to augment phase-based material properties
  PHY.OC.Cp   = 1200.;                     % [J.kg-1.K-1] $\C_p$ specific heat capacity
  PHY.OC.Dens = 2900.;                     % [kg.m-3] average density assigned to oceanic crust integration points
  PHY.OC.Hp   = 0.46E-06;                  % [W.m-3] $H_r$ radiogenic heat production [value for subalkalic gabbro in Hasterok and Webb (2017)] 
  PHY.OC.K    = 3.1;                       % [W.m-1.K-1] thermal conductivity [after Kardell et al. (2021)]  
  PHY.OC.thermal_expansivity = 0.0;        % [K-1] a realistic value is 3.5E-05 [see Afonso et al., 2004]
  PHY.OC.Cohesion = [10.0E06, 10.0E06];    % initial and strain-softened cohesion [decay along the same accumulated plastic strain than the rest of the domain]

  %% - SERPENTINIZATION
  SETTINGS.SERP = [];
  SETTINGS.SERP.make = false;              % true requires SETTINGS.HYDRO.make = true to get the hydrothermal domain;
  if contains(COMPSET,"serp")
      SETTINGS.SERP.make = true; 
      PHY.SERP = [];                      % serpentinization physical parameters
      PHY.SERP.Tmin = 100.;               % [ºC] minimum temperature for serpentinization to occur [e.g. from Wenner and Taylor 1971]
      PHY.SERP.Tmax = 450.;               % [ºC] maximum temperature for serpentinization to occur [e.g. from Wenner and Taylor 1971]
      PHY.SERP.k0 = 1.0E-10;              % [s-1] maximum kinetic coefficient for the serpentinization reaction
      PHY.SERP.equation = "SKGAU"         % {'E2006','M2012','SKGAU'} see serp_underplate() for details
      %PHY.SERP.Hserp  = 300.;            % [ºC] exotermic heat of olivine -> serpentinites transform ([ºC] : i.e. normalized by heat capacity) Ros_al2017 (after Iyer_al2012)
      PHY.SERP.Hserp_heat = 2.9e5;        % [J.kg-1] exothermic heat of olivine (MacDonald and Fyfe, 1985). TODO: replace Hserp to introduce this as source term in the heat equation
                                          % H = 19.5 kcal.mol-1 serpentinite => 19.5*4.184/(2*140.69)*1000 = 290 kJ.kg-1 Mg2SiO4
      PHY.SERP.E2fac = [];                % serpentinization domain: value in the range [0,1]. If not empty, E2min is calculated at each timestep as E2min = E2fac * (max(E2) - min(E2)) + min(E2); as Ros_al2017
      PHY.SERP.ErPnom = 5.0E-14;          % nominal II plastic strain invariant for which the kinetics of serpentinization are only controlled by temperature
      PHY.SERP.sfun = "gau";
      PHY.SERP.olivine_density = 3300;      % [kg.m-3]
      PHY.SERP.olivine_H2_molar_ratio = 15; % Olivine -> Serp + 15*H2          Based on Malvoisin 2012; should be different for different isochemical EQs
      PHY.SERP.olivine_molar_mass = 146.7;  % [gr.mol-1] x=0.91; 2*(x*Mg+(1-x)*Fe)+SiO4 = 2*(x*24.305+(1-x)*55.845)+ (28.44+4*15.999) g/mol  
      %PHY.SERP.PARAM.dxmin = 3.  *km; 
      %PHY.SERP.PARAM.dzmax = 15. *km;    % maximum depth [below the sea floor] at which serpentinization can be found %cut off  km depth for generation of serpentinite (generally lowest point is at ~3km depth)
  end

  %% - SURFACE PROCESSES    
  SP = [];
  SP.make = false;
  if contains(COMPSET,"spro")
      SP.make = true;                                  % true: compute surface processes

      %% - SP physics                                                  
      % ----------------------------------------------- 
      SP.tp_isoc = true;                               % true: track isochrons (at the expense of slowing model down)
      SP.iso_sav = 10;                                 % Time steps between isochrons to save. TODO: switch into [Ma] criterion as Kinedyn
      SP.func = "diff_topo_sealand";            % Function to calculate erosion and sedimentation
      %% - SP physics                                                
      %------------------------------------------------ %
      SP.dt    = 1000.*year;                           % [s] Time step for the surface processes
      SP.kappa = 0.25/year;                            % Hill-slope diffusion [m^2/s]
      SP.c     = 1.e-3;                                % Discharge transport coefficient
      SP.nexp  = 1;                                    % Inverse Hack's law
      SP.alpha_sed = 1./year;                          % Precipitation rate
      SP.pelagic_t = "cnst";                           % Pelagic sedimentation type: 'cnst' for constant or 'var' for a function. Currently not used
      SP.pelagic_rate = 2.e-05/year;                   % [m s-1] Pelagic sedimentation rate [m/s]. Empty (i.e. = []) for no external source of sediment
      SP.pelagic_Hdep = false;                         % true: make a linearly increasing depth-dependent decay: pelagic_rate = H/Hmax * SP.pelagic_rate
      SP.pelagic_Hmax = [];                            % [m] if not empty, water depth Hmax. Otherwise, calculated as model maximum depth
      SP.pelagic_rate_below_ccd = 1.e-05/year;         % [m s-1] Pelagic sedimentation rate below carbonate compensation depth     
      SP.ccd = 4000.;                                  % [m] carbonate compensation depth - make (= Inf) to not consider it at all
      SP.kdistalb = 1.0e-15/year;                      % basement-top, set as minimum diffusivity everywhere [m^2/s]
      SP.ktidal  = 1.e3/year;                          % Hill-slope diffusion under the sea [m^2/s] (later multiplied by the decay)
      SP.kdecay = 2.5e-02;                             % Decay of the tidal diffusion under the sea
      SP.scale = false;                                % true: scale variables
      SP.K     = 2.1;                                  % [W.m-1.K-1] Thermal conductivity for the sediments. If []; i.e. isempty(.) , equals that from the upper crust
      SP.Cp    = 1200.;                                % $\C_p$ heat capacity [  ",     ",  "]
      SP.Hp    = 1.3e-6;                               % $H_r$ radiogenic heat production for sediments;
      SP.lriver = [0. 0.] * km;                        % [m] (>=0) river length for aerial diffusivity parametrization [left and right rivers coming out of the domain] 
      SP.q_bc = [0. 0.] / year;                        % [m/s] Sediment flux into the model from boundaries. Positive for source. Negative for sink.
      SP.LIVDcomp = false;                             % if true, a dynamic SP.q_bc will nudge topography toward initial conditions
      %% - Neumann BC                                                
      SP.bc_type = "Neumann";
      SP.bc_val = [0. 0.];                             % [m] if bc_type='Dirichlet'. Not implemented for Neumann.     
      SP.Dens = 2450.;                                 % mean from Shulgin and Artemieva (2019)
      SP.thermal_expansivity = 0.0;                    % [K-1] a realistic value is maybe ~2.0E-05. Used in density calculation.  
  end % if contains(COMPSET,"spro")
  %--------------------------------------------------------------------------
  %% - Sea level                                                 
  %--------------------------------------------------------------------------
  PHY.rho_w = 1028.;                               % [kg m^-3] density of sea water  (set 0 for no sea level load)
  PHY.c_sealevel = 30. *km;                        % Typical crustal thickness for sea level [m]
  PHY.sea_level = [];                              % Automatically calculated. Set to specific value only for debugging.
  PHY.add_atmPres = true;                          % whether to add atmospheric pressure on top of free surfaces and for winkler pressures

  %% TRACKING POINTS
  SETTINGS.tracking = true;                        % true: track points along time steps
  SETTINGS.res_tp = 1. *km;                        % [m] initial resolution of tracking point lines
  SETTINGS.sedtocrust = true;                      % make thick sediments upper crust phase for the mechanical
  SETTINGS.phasesed = 4;                           % TODO: include all crustal-like properties into the mechanical via the corresponding properties at IPs 
  SETTINGS.phaseoc  = 5;

  %% KINETIC FAULTS
  SETTINGS.FAUL = [];
  SETTINGS.FAUL.make = false;

  %% - HYDROTHERMAL FLOW
  SETTINGS.HYDR = []; 
  SETTINGS.HYDR.make = false;
  
  if contains(COMPSET,"hydr") || contains(COMPSET,"serp")
      SETTINGS.HYDR.make  = true;
      SETTINGS.HYDR.model = false;       % false for parameterised hydrothermal cooling, true for Darcy flow

      if SETTINGS.HYDR.model % settings for the hydrothermal model
          % FLUID PROPERTY SETTINGS - only used for SETTINGS.HYDR.model = true
          % 'PROST'    - properties calculated using tables generated by PROST [PROperties of Water and STeam "PROST 4.1"]
          % 'tplinear' - density function of constant alpha and beta
          % 'PROSTmex' - alpha_f & beta_f calculated using tables generated by PROST
          %              Rho_f, Mu_f & Cp_f calculated using mex function "mexsteam"
          %              (EXPERIMENTAL VERSION)
          % reaction 
          SETTINGS.HYDR.reaction = false;
          SETTINGS.HYDR.feedback_k = false;
          SETTINGS.HYDR.start_from_seawater = true; % concentration of Ca, So4 initialized as seawater
          % fluid properties
          SETTINGS.HYDR.f_prop = "PROST";  % fluid property switch (PFR, HYD), where the file containing the lookup table is:
          SETTINGS.HYDR.H2O_table_file = fullfile(dsninput,"prost_tables_P2kbar_T1200_falko.mat"); % contains {'ALPHA','BETA','CP','H','MU','PP','RHO','TT'} all REAL [500,300]
          
          SETTINGS.HYDR.adv_scheme = "FVM";  % advection scheme {"FVM","SLM"}
          
          switch SETTINGS.HYDR.adv_scheme
            case 'FVM'
              SETTINGS.HYDR.FVM_order = 1; % first or second order finite volumes
            case 'SLM'
              SETTINGS.HYDR.FVM_order = 1; % first or second order finite volumes
              SETTINGS.HYDR.SLM_BT = 'RK4';     % 'EU1' == Euler (very inaccurate !!!!)
              SETTINGS.HYDR.dt_limit = 1e-2;    % limit dt so we dont move more than this frac of domain size
          end % switch HYDR.SETTINGS.adv_scheme
          SETTINGS.HYDR.dt_max = 10. * year; % [s] max timestep for hydrothermal model
          SETTINGS.HYDR.dt_ini = 60;         % [s] max timestep for first step
          SETTINGS.HYDR.dt_maxincr = 1.2;    % factor that limits increase of dt between time steps
          SETTINGS.HYDR.disp_profiling = false;
          
          SETTINGS.HYDR.scale  = "vel";       % {"temp","vel"} temperature or velocity advection scaling (PFR)
          SETTINGS.HYDR.vel_smooth = false;   % LOGICAL, true: velocity smoothing (PFR)
          SETTINGS.HYDR.compress = false;     % LOGICAL, true: compresible fluid (PFR)
          SETTINGS.HYDR.ftop = "open";        % {"closed","open"} top convection style (PFR) - still not in use: automatically the "open" top boundary is selected
          SETTINGS.HYDR.use_overPfb = true;
      end % if SETTINGS.HYDR.model
      
      HYDRO = [];                     % physics for hythermal model/parameterized cooling

      HYDRO.G = PHY.G;
      % overall selection of hydrothermal domain for both the hydrothermal model and the parameterised enhanced thermal conductivity
      HYDRO.dzmax = Inf *km;          % hydrothermal domain: [m] maximum depth from the topography surface on which hydrothermal efect is active  
      HYDRO.dzcrustmax = 10. *km;     % minimum crust thickness (as approximate global measure) to initiate overall hydrothermal flow
      HYDRO.Tmax  = 600.;             % hydrothermal domain: [ºC] maximum temperature [by permeability closure]
      HYDRO.gradTmin = 300/10000.;    % (400/10000=0.04) [ºC/m] evaluated as the minimum averaged vertical thermal gradient for the surface up to a maximum HYDRO.Tmax isotherm 
      %HYDRO.E2min = 1.E-15;           % hydrothermal domain: minimum II strain rate invariant to select the hydrothermal domain   
      HYDRO.E2fac = [];               % value in the range [0,1]. If not empty, E2min is calculated at each timestep as E2min = E2fac * (max(E2) - min(E2)) + min(E2); as Ros_al2017
      HYDRO.ErPmin = 0.;              % hydrothermal domain: minimum II plastic strain rate invariant
      % larger ErPmin => delay & reduced area in the onset of hydrothermal flow (and serpentinization in its case)
      
      HYDRO.PARAM = [];
      if ~contains(COMPSET,"hydr")
        HYDRO.PARAM.cfmax = 1.;
      else
        HYDRO.PARAM.cfmax = 8.;         % maximum factor multiplying thermal conductivity for parameterized hydrothermal cooling [HYDRO.model==false]. Morgan_Chen1993b got 8 as optimal 
      end
      HYDRO.PARAM.ErPfault = 1.E-13;    % nominal II plastic strain invariant at faults matching HYDRO.param.cfmax       
      HYDRO.PARAM.vfun = "strainrate";  % {'constant','strainrate','permeability'}. Functional variable for getHThermalCondFac(). 
      HYDRO.PARAM.sfun = "exp";         % {'gau','exp','sph'} submodel for HYDRO.param.vfun. See getHThermalCondFac() for details.

      PHY.HYDR = HYDRO;
      %for fn = fieldnames(HYDRO)'
      %  PHY.HYDR.(fn{1}) = HYDRO.(fn{1});
      %end
  end % 'hydro' COMPSET
  
end % function