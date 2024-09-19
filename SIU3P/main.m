%  Rift2Ridge is a visco-elasto-plastic finite element code for geodynamics, originally branched from MILAMIN
%      
%  Two major stages of development have been carried out in Rift2Ridge
%  @ Royal Holloway University of London. Earth Sciences Department.
%      - Melting [Marta Perez-Gussinye]
%      - Friction angle, cohesion and viscous strain weakening
%      - Elasticity
%      - Varying rheologies along a phase
%      - Dike/sill emplacement and postsolidification heat release [Elena Ros]
%      - Density dependent on temperature and depletion [Elena Ros and Marta Perez-Gussinye]
%      - Free surface algorithm [Miguel Andres-Martinez]
%      - Random Pef [Miguel Andres-Martinez]
%      - Viscous strain weakening dependent on temperature [Albert de Monserrat, Leon Liu and Elena Ros]
% 
%  @ MARUM, University of Bremen, Geodynamics Research Group [2017 onwards]
%      Major code reorganization, modularization, improved parallelization, water loading, layer breakup,
%      lossless remeshing with explicit interface tracking, direct solution melting, kinematic-dynamic integration...
%
%==========================================================================
%% CLEARING AND INITIALIZATION
%==========================================================================

clc;
warning('off', 'MATLAB:triangulation:PtsNotInTriWarnId');

MODEL   = "rift2ridge2D";
VERSION = "SIU3P";
HOME    = getenv('HOME');
LOGNAME = getenv('LOGNAME');                                               % computer user
WORK    = getenv('WORK');

codebase = "siu3p";

dsnmod   = fullfile(HOME,"docs",MODEL,VERSION);                            % model code

global pypath
global year

%cwd = pwd;
addpath(fullfile(dsnmod,"input"))
addpath(fullfile(dsnmod,"hydrothermal"))
addpath(fullfile(dsnmod,"main_functions"))
addpath(fullfile(dsnmod,"melting_serpentinization"))
addpath(fullfile(dsnmod,"kinedyn_functions"))
addpath(fullfile(dsnmod,"mesh_gen_bc"))
addpath(genpath(fullfile(dsnmod,"plot_functions")))
pypath = fullfile(HOME,"docs",MODEL,"scripts","python");                   % path to python scripts [tsearchPy.py]

dsn      = fullfile(HOME,"docs",MODEL,VERSION,"data");                     % data input [scenario input parameters] - for future structued input
dsninput = fullfile(HOME,"docs",MODEL,VERSION,"input");

readCASE;                                                                  % local environment
scn = join([meshnam,nnn],".");
% dsnscn = fullfile(dsn, region, event, join([meshnam,nnn],"/"));

CASE_0 = join([codebase, COMPSET, region, event, scn], ".");               % common for all members
CASEROOT_0 = fullfile(WORK,MODEL,VERSION,"data",CASE_0);                   % common for all members
readCF;                                                                    % local environment: makes CF list
[CF.p,CF.MC,CF.MCcons] = readMC();                                         % Monte Carlo parameters

% NOTE: modify for ensemble simulation & DA
CASE = join([CASE_0,"001"],".");                                           % in an ensemble case, these 3 will be vectors of strings
CASEROOT    = fullfile(WORK,MODEL,VERSION,"data",   CASE);                 % scratch model output [mesh...]
DOUT_S_ROOT = fullfile(WORK,MODEL,VERSION,"archive",CASE);                 % e.g. could be replaced by $DATA_SHARE

t0=tic;
fprintf(1,"Start Time: %s\n",datetime("now"));                             % print start date time

%==========================================================================
% INITIALIZATION
%==========================================================================
[MUTILS_PATH, TRIANGLE_PATH, SSPARSE_PATH] = load_dir(0);
addpath(genpath(MUTILS_PATH));
addpath(genpath(SSPARSE_PATH));
addpath(fullfile(SSPARSE_PATH,"CXSparse","MATLAB","CSparse"));              % needed in some computer (e.g. cluster cluster)
triangle_mode = "text";                                                    % triangle mode (choose between "text" and "binary")
%addpath_agmg();                   % adds paths to AGMG

% SET THE DEFAULT ROOT RENDERER TO EITHER zbuffer OR OpenGL
% - zbuffer makes nicer plots than OpenGL, does not support transparency
% - zbuffer becomes very slow for larger 3D plots (e.g. a 3D scatter plot)
set(0, 'DefaultFigureRenderer', "zbuffer");                                % set default renderer

%==========================================================================
% MAKE OUTPUT DIRECTORIES
%==========================================================================
if ~isfolder(CASEROOT_0)
    mkdir(CASEROOT_0);                                                     % mesh and temporary file output
end

if ~isfolder(CASEROOT)
    mkdir(CASEROOT);
end
 
if ~isfolder(DOUT_S_ROOT)
    mkdir(DOUT_S_ROOT);
end
addpath(DOUT_S_ROOT);

%==========================================================================
% GET PARAMETERS AND SWITCHES
%==========================================================================
[SETTINGS, PHY, SOLVER, SP, PLOT] = param_defaults(dsninput, COMPSET);
[SETTINGS, PHY, SOLVER, SP] = change_model(SETTINGS, PHY, SOLVER, SP);                                   % scenario specific model parameters
SETTINGS.layshift = SETTINGS.layshifts(1);

if ~loadsave
    %==========================================================================
    %% USER INPUT:                                                      
    %==========================================================================
    
    %% SETTINGS
    %--------------------------------------------------------------------------
    %%      - General settings                                          
    %--------------------------------------------------------------------------
    % tstep_type = 'cnst' for a constant time step defined by dt and 'courant'
    % for a courant time step criteria:
    %             dt = C*(min(inc_surface_element)/max(velocity)).
    tstep_type  = "cnst";
    % Courant time step criteria
    % C           = 0.025;

    %--------------------------------------------------------------------------
    %%      - Units & unit conversions                                                    
    %--------------------------------------------------------------------------
    km      = 1000.;
    year    = 365.25*24*60*60;                                                 % [s]
    ma      = 1e6*year;                                                        % [s]
    bar2pascal = 100000.;                                                      % []
    atm2pascal = 101325.;                                                      % []
    
    %% PHYSICS
    %--------------------------------------------------------------------------
    %%      - Time                                                      
    %--------------------------------------------------------------------------
    dt          = PHY.dt;                                                      % [s] initial time step
    time_int    = PHY.time_int;                                                % [s] Time interval for the model.
    nsteps      = PHY.time_int/dt;                                             % Initial mumber of steps

    %--------------------------------------------------------------------------
    %%          * General rheological parameters                        
    %--------------------------------------------------------------------------
    PHY.RHEOL = load_rheol(PHY.RHEOL);                                            % uses PHY.RHEOL.type to load rheology parameters

    %--------------------------------------------------------------------------
    Ylim = [PHY.ymin1, PHY.ymax1, PHY.ymax2, PHY.ymax3];                                             % domain limits
    % Number of initial points at geometry boundaries ordered as follows:
    %       __7__      _PHY.ymax3
    %             |
    %       __5__|6    _PHY.ymax2
    %            |
    %       __3__|4    _PHY.ymax1
    %            |
    %       __1__|2    _PHY.ymin1
    cont_points = [floor((PHY.xmax-PHY.xmin)/SETTINGS.hres(1)) floor((PHY.ymax1-PHY.ymin1)/SETTINGS.hres(1)) ...   % 1 2   DWOM: dry and wet olivine mantle
                  floor((PHY.xmax-PHY.xmin)/SETTINGS.hres(2)) floor((PHY.ymax2-PHY.ymax1)/SETTINGS.hres(2)) ...   % 3 4
                  floor((PHY.xmax-PHY.xmin)/SETTINGS.hres(3)) floor((PHY.ymax3-PHY.ymax2)/SETTINGS.hres(3)) ...    % 5 6
                  floor((PHY.xmax-PHY.xmin)/SETTINGS.hres(4))] + 1;                                                % 7 
    Elsizes = SETTINGS.hres(2:end).^2;                                                                        % [1e8 1e7 1e7 1e7]
    
    vres = SETTINGS.hres(1:end-1);                                                           % soft bottom-top vertical interface resolution [m]
    SETTINGS.GEOres = [SETTINGS.hres(1) reshape([vres; SETTINGS.hres(2:end); vres],1,[])];   % resolution for each interface, in increasing interface order
    % error('Ya lo se, Walter')
    % GEOres = 
    
    ngids = length(SETTINGS.GEOres);
    tolerancef = repelem(SETTINGS.tolerancef(end),1,ngids);
    tolerancef([1 3:3:ngids-2 2:3:ngids 4:3:ngids]) = repmat(SETTINGS.tolerancef(1:end-1),1,3);
    SETTINGS.tolerance = SETTINGS.GEOres .* tolerancef;
    clear("vres","ngids","tolerancef")

    NUMSCALE = SETTINGS.NUMSCALE;                                                       % NUMSCALE  - see param_defaults()  - only NUMSCALE.t set
    NUMSCALE.mu0 = PHY.Shearm(1)*dt;                                                    % 'mantle' Viscosity scaling [Pa s]
    NUMSCALE.L0  = PHY.ymax3-PHY.ymin1;                                                 % Length scaling [m]
    NUMSCALE.V0  = abs(PHY.ext_rate);                                                   % Velocity scaling [m/s]
    NUMSCALE.t0  = NUMSCALE.L0/NUMSCALE.V0;                                             % Time scaling [s]
    NUMSCALE.P0  = NUMSCALE.mu0/NUMSCALE.t0;                                            % Pressure scaling [Pa]
    
    if isfield(PHY,"sea_level") && ~isempty(PHY.sea_level)                                 % fixed value for model testing
        sea_level = PHY.sea_level;
    else
        incc = (PHY.ymax3-PHY.ymax1);                      % [m] total crust thickness
        incl = (PHY.ymax2-PHY.ymax1);                      % [m] lower crust thickness
        incu = (PHY.ymax3-PHY.ymax2);                      % [m] upper crust thickness
        csea = PHY.c_sealevel;
        c_inf_sl = csea * incl / incc;           % [m] JGP: check if this OK for rho_w ~= 0.0? it's perhaps better to fix the sea level to 0  
        sea_level = ( PHY.Dens(3)*(c_inf_sl-csea) - PHY.Dens(2)*c_inf_sl ...  % sea level [m]
                    - PHY.Dens(1)*incc + PHY.Dens(1)*csea ...
                    + PHY.Dens(3)*incu + PHY.Dens(2)*incl) / PHY.Dens(1);
        clear("csea","incc","incl","incu","c_inf_sl")
        %foo = (PHY.Dens(3)*incu + PHY.Dens(3)*csea*(incl/incc - 1) + PHY.Dens(2)*incl*(1 - csea/incc) + PHY.Dens(1)*(csea - incc))/PHY.Dens(1);
    end
    SP.sealevel = sea_level;
    hydroactive = false;                                                                 % later used only for 'HYDRO' compsets
    %--------------------------------------------------------------------------
    
    % In order to initialize track points and tracers, a matrix TRACKP needs to
    % be defined. This matrix has the size (number of dimensions x number of
    % track points).
    %
    % Example:
    % To create a rectangular mesh at the surface from -5 km to 5 km at the
    % horizontal and 5 km depth:
    
    % tp_x = (-20:0.5:0)*km; % X initial distribution of the tracked points
    % tp_y = (-10:0.5:-2.3)*km; % Y initial distribution of the tracked points
    % [TP_x,TP_y] = meshgrid(tp_x,tp_y); % Calculates a regular mesh
    % TRACKP{1} = [TP_x(:)'; TP_y(:)']; % Defines track points
    % 
    % % Creates rectangular-element indexes for the tracking points
    % E2N_TP(1,:) = 1:length(tp_y):((length(tp_x)-1)*length(tp_y));
    % E2N_TP(2,:) = 2:length(tp_y):((length(tp_x)-1)*length(tp_y));
    % E2N_TP(3,:) = length(tp_y)+2:length(tp_y):((length(tp_x))*length(tp_y));
    % E2N_TP(4,:) = length(tp_y)+1:length(tp_y):((length(tp_x))*length(tp_y));
    % E2N_TP = repmat(E2N_TP,1,length(tp_y)-1);
    % ADD_INDEX = repmat((0:length(tp_y)-2),length(tp_x)-1,1);
    % ADD_INDEX = repmat(ADD_INDEX(:)',4,1);
    % E2N_TP = E2N_TP + ADD_INDEX; 
    
    % % Creates line indexes for the tracking points
    % tp_x = linspace(PHY.xmin,PHY.xmax,cont_points(7)); 
    %     % X initial distribution of the tracked points
    % tp_y = PHY.zlab*ones(size(tp_x));
    %     % Y initial distribution of the tracked points
    % TRACKP = [tp_x; tp_y]; % Defines track points
    % % LINES_TP = repmat((1:length(tp_y))',1,length(tp_x)); 
    % % LINES_TP = LINES_TP(:); 
    
    TRACKP = [PHY.xmin:SETTINGS.res_tp:PHY.xmax; ...                           % LAB TRACK
              PHY.zlab*ones(1,floor((PHY.xmax-PHY.xmin)/SETTINGS.res_tp)+1)];
    
    dzLIVD = 0.0;                                                              % [m] tracking of LIVD compensation          
    
    %% INITIAL SETUP
    BREAKUP.breaktime = Inf;
    BREAKUP.bool = false;

    if SETTINGS.FAUL.make
        FAULTS = readFAULTS();
    else
        FAULTS = [];
    end
    
    %==========================================================================
    % MESH GENERATION:
    %==========================================================================
    fprintf(1, "PREPROCESSING:      \n");
    
    %Intersect_ID = [];
    
    if SETTINGS.geomf == "make_geometry3l"                        
        [GEOMETRY, Geo_id] = make_geometry3l(PHY.xmin, PHY.xmax, PHY.ymin1, PHY.ymax1, ...
        PHY.ymax2, PHY.ymax3, cont_points, SETTINGS.ini_deformation, PHY.Dens, ...
        PHY.sigma_moho, PHY.topo_moho);    
    else
        fprintf("making non-default geometry\n");
        make_geom_handle = str2func(SETTINGS.geomf);
        [GEOMETRY, Geo_id] = make_geom_handle(PHY.xmin, PHY.xmax, PHY.ymin1, PHY.ymax1, ...
        PHY.ymax2, PHY.ymax3, cont_points, SETTINGS.ini_deformation, PHY.Dens, ...
        PHY.sigma_moho, PHY.topo_moho);  
    end
    
    GEO = geometry2GEO(GEOMETRY, Geo_id);                                   % structured domain definition. plotGEO(GEO)
    [GCOORD, ELEM2NODE, Point_id, Phases] = generate_meshGEO(GEO, ...       % GCOORD \in [2,nnod], ELEM2NODE \in [nnpe,nel], Point_id \in [nnodce], Phases \in [nel] 
    Elsizes, triangle_mode, "triangle", fullfile(CASEROOT,"mesh"));         % where nnpe is number of edge+vertex nodes per element
    
    nnod = size(GCOORD,2);                                                  % number of [vertex+edge] nodes in the continuous mesh
    ndof = size(GCOORD,1);                                                  % number of degrees of freedom for position [2D => 2]
    nel  = size(ELEM2NODE,2);                                               % plotGEO(GEO); hold on; plot_meshF(ELEM2NODE,GCOORD)
    zuco = GEO(end-1).coo(2,[1 end]);                                       % [m] altitude of upper corners 

    %augment elements and node set with 7th (central) node
    nnod6 = size(GCOORD,2);
    ELEM2NODE(7,:)  = nnod+1:nnod+nel;                                      % [nodes per el, nel] 
    GCOORD          = [GCOORD, ...                                          % [2, nnod7]
    [mean(reshape(GCOORD(1, ELEM2NODE(1:3,:)), 3, nel));...           
     mean(reshape(GCOORD(2, ELEM2NODE(1:3,:)), 3, nel))]];
    nnod7 = size(GCOORD,2);                                                 % nnod = nnodce + nel (as central node for each element is included now)
    clear("nnod")
    
    fprintf(1, ['\n Number of elems:   ', num2str(nel),'\n']);
    fprintf(1, ['\n nnod7 nodes:   ',     num2str(nnod7),'\n']);           % Note: still Point_id \in [1,nnod] 
    
    if PLOT.make
        PLOT.Point_id = Point_id;
    end
    
    corners_id = getCornerids(GCOORD, Point_id);
    ext_srate  = PHY.ext_rate / diff(GCOORD(1, corners_id(:,1)));           % [1/s] x-strain rate
    
    %==========================================================================
    % SOLVER
    %  nip                - Number of integration points (6 or higher)
    % (reorder = 'amd')   - AMD reordering
    % (reorder = 'metis') - METIS reordering
    % (method = 'std')    - Standard matrix computation
    % (method = 'opt')    - Optimized matrix computation
    %==========================================================================
    
    nip             =       6;
    reorder         =   'amd';
    method          =   'opt';
    
    %STORAGE
    extension       = zeros(1,nsteps);
    E2all           = ones(nel,SOLVER.nip_stress)*ext_srate;
    Mu_all          = zeros(nel,nip);
    Mu_b_all        = zeros(nel,nip);
    Mu_dis_all      = zeros(nel,nip);
    Mu_dif_all      = zeros(nel,nip);
    RHEOLvar        = ones(nel,1);
    
    % Initialize accumulated gradient of deformation
    F_xx = ones(nel,SOLVER.nip_stress);                                          % [nel,nip], with nip=6
    F_xy = zeros(nel,SOLVER.nip_stress);                                         % [nel,nip]
    F_yx = zeros(nel,SOLVER.nip_stress);                                         % [nel,nip]
    F_yy = ones(nel,SOLVER.nip_stress);                                          % [nel,nip]
    I2.f = zeros(nel,SOLVER.nip_stress);                                         % [nel,nip]
    I2.p = zeros(nel,SOLVER.nip_stress);                                         % [nel,nip]
    I2.c = zeros(nel,SOLVER.nip_stress);                                         % [nel,nip]
    
    % Initialize eneu and theta for elasticity
    THETA_all       = zeros(nel,nip);                                            % [nel,6] viscosity at integration points - needed for initial power law guess
    TAU_xx_old      = zeros(nel,nip);
    TAU_yy_old      = zeros(nel,nip);
    TAU_xy_old      = zeros(nel,nip);
    ttotal = 0;
    track_remesh = []; 
    nw_it = [];
    SOLVER.residue  = [];
    SOLVER.ext_srate    = ext_srate;
    % Declare isochrons rows: 1) x-coordinates, 2) y-coordinates, 3) step
    % index, 4) erosion/deposited depth
    isoclen = sum((max(Point_id)-1==Point_id));
    ISOCHRONS = [GCOORD(:,Point_id==max(Point_id)-1);                         % [x;z] - [2,ntop] nodes at the top including corners
                 zeros(1,isoclen)];                                           % [istep]  for this isochron    
    [~,indx_isoc] = sort(ISOCHRONS(1,:));                                     % guarantee monotonic increase in x coordinate - x
    ISOCHRONS(1:2,:) = ISOCHRONS(1:2,indx_isoc);                                  %    "           "                "     "      - z
    ISOCHRONS = [ISOCHRONS; ISOCHRONS(2,:)];                                  % add 4th row to preserve the isochron y-coordinates at the time it was created
    Basement = ISOCHRONS(1:2,:);                                              % extract first isochron as Basement 
    clear("isoclen")
    
    % Calculate random damage
    if PHY.SS.rand_s
        if SETTINGS.ini_deformation == 6                                      % random damage weak seed
            RDWS_handle = str2func(PHY.RDWS.function);                        % default: "gaussian_band_rand_damage"
            RDWS_factorIP = RDWS_handle(PHY.RDWS, GCOORD, ELEM2NODE, nip);    % [2] cell, where each element is an [nel,nip] array
        else
            RDWS_factorIP{1} = ones(nel,nip);                                 % for friction angles
            RDWS_factorIP{2} = ones(nel,nip);                                 % for dislocation & diffusion
        end
        SS = random_damage_WS(PHY, RDWS_factorIP, nel, nip, false);           % init SS with keys: ['RDWS_factorIP','Random','Phi1_rand','Phi2_rand','I2_phi1_rand','I2_phi2_rand'...], all REAL [nel,nip].
        % plot(GCOORD(1,ELEM2NODE(7,:))',SS.Phi1_rand(:,1)*180/pi,'.'); title('Phi values when random noise is activated')
        % xlabel('Distance [km]'); ylabel('Friction angle [^\circ]')
    else
        SS = [];
    end
    if SETTINGS.ini_deformation == 5
        WScoord = PHY.WS.coord;                  % only for tracking
    end
    
    if SETTINGS.HYDR.make
        if SETTINGS.HYDR.model && ~isfield(PHY.HYDR,'H2O_TABLES')
            PHY.HYDR.H2O_TABLES = load_H2O_table(SETTINGS.HYDR); 
        end
    end

    if SETTINGS.MELT.make
        PHY.MELT = parameters_melt(PHY.Dens);                             % melting physical parameters
        [SETTINGS, PHY, SOLVER, SP] = change_model(SETTINGS, PHY, SOLVER, SP);
    else
        PHY.MELT = [];                                                              
        PHY.MELT.Ts0    = 0.;
        PHY.MELT.dTs_dP = 0.;
        PHY.MELT.dTs_dF = 0.;
    end

    % init melt-origin crustal thickness (as diagnostic), prognostic mantle depletion (% of melted rock volume) due to melt, and prognostic serpentinization
    ocrust = []; % [2,:] polyline representing the bottom of the oceanic crust
    ocrust_melt = []; % melting component of oceanic crust, as polyline
    ocrust_serp = []; % serpentinite component layer of oceanic crust, as polyline 
    MELTING = [];                                                              % struct for melting variables                                                                  
    MELTING.Crust_thickness = zeros(1,nsteps);                                 % [nsteps]
    MELTING.Trackp = [];                                                       % melting track points
    MELTING.tp_melt = false;                                                   % init switch to track melting once in the surface
    
    Dpl   = dpl_init(GCOORD, ELEM2NODE, Phases, PHY.zlab);                     % [1,nnod7] | plot_dplF(Dpl, GCOORD, ELEM2NODE, Point_id, 0,dt,ma)
    Dserp = zeros(1,nnod7);                                                    % [1,nnod7]
    if SETTINGS.SERP.make
        SERPEN.trackt = [];                                                    % Tracking initial time step serpentinization
        SERPEN.serp = [];
        SERPEN.serp.total = zeros(1,nsteps);                                   % [m3] total volume of serpentinite
        SERPEN.serp.rate = zeros(1,nsteps);                                    % [m3.yr-1] serpentinite generation rate
        SERPEN.H2p = [];                                                       % [mol] H2 production in the serpentinization reaction
        SERPEN.H2p.total = zeros(1,nsteps);                                    % [mol] integrated H2 production [in the 1-m depth profile]
        SERPEN.H2p.rate =  zeros(1,nsteps);                                    % [mol.yr-1]
    end
    
    if size(PHY.RHEOL.Ndis,2) > 1                                              % declare rheologic variation function
        rheol_var_handle = str2func(PHY.RHEOL.var_f);                          % "rheol_var_Dpl_ip" by default. Several options seems available in ".../input" folder
        RHEOLvar = rheol_var_handle(Dpl, ELEM2NODE, Phases, ...                % 1x2 cell array, each: [nel,nip]
        SOLVER.nip_stress, GCOORD, false);
    end
    
    if SP.make                                                                 % active surface processes: declare surface process function
        surf_proc = str2func(SP.func);
    end
    BDTdepth = zeros(1,nsteps);                                                % mechanical brittle-ductile transition [shallower] depth
    BDTat = zeros(2,nsteps);
    %==========================================================================
    % CALCULATE INITIAL CONDITIONS
    %==========================================================================
    EL_IP = [];
    EL_IP.K  = repmat(PHY.K(Phases),  1, nip);       % [nel,nip] thermal conductivity at integration points
    EL_IP.Cp = repmat(PHY.Cp(Phases), 1, nip);       % [nel,nip] specific heat at constant pressure
    EL_IP.Hp = repmat(PHY.Hp(Phases), 1, nip);       % [nel,nip] radiogenic heat production
    EL_IP.Hs = zeros(nel,nip);                       % [nel,nip] shear heating
  
    if PHY.load_T                                                                % load initial temperature field from PHY.Tini_dir file 
        Tini_handle = str2func(PHY.Tini_load);
        Temp = Tini_handle(PHY.Tini_dir,GCOORD,PHY.temp_bc);
    else
        % initial temperature field as steady state solution
        [Bct_ind, Bct_val] = set_bcs_temp_ini(...
            GCOORD, Point_id, PHY.temp_bc_ini, PHY.temp_bc_depth_ini, ...
            PHY.Ttop, SETTINGS.ini_deformation, PHY.t_sigma_x, PHY.t_sigma_y);
        
        Temp = PHY.temp_bc_ini * ones(1,nnod6);                                     % [1,nnod7] (ºC) 
        EL_IP.Rho_ip = repmat(PHY.Dens(Phases),1,nip);                              % [nel, nip] (kg/m3) density at integration points [Phases is defined at element]
    
        Temp = thermal2Dss(ELEM2NODE(1:6,:), GCOORD, Temp, ...                      % Temp [1,nnod6] equilibrium solution
            EL_IP.Rho_ip .* EL_IP.Cp, EL_IP.Hp, EL_IP.K, EL_IP.Hs, ...              % [nel,nip] variables        
            Bct_ind, Bct_val, 100*ma, reorder, SOLVER.nelblo, true);                % plot_tF(Temp, GCOORD, ELEM2NODE)
        Temp = addCentroidVal(ELEM2NODE(1:6,:), Temp);                              % [1,nnod7] to make it usable by current melting routines
    end

    GEOn = getGEOn(GCOORD, ELEM2NODE, Point_id, Phases);                                      % GEO augmented with additional triangle-generated vertex nodes

    % UPDATE DENSITIES                                                                        note densities do not affect the thermal steady state
    EL_IP.Rho_ip = dens_temp_Dpl(GCOORD, ELEM2NODE, Phases, GEOn, Basement, ...
        Temp, Dpl, nip, PHY, SP, ocrust);                                           % [nel, nip=6]
        
    % CALCULATE STRESSES FOR WINKLER BOUNDARY CONDITIONS
    if contains(PHY.bc_t,"winkler")                                                 % get TAU_ast [Pa]:: lithostatic pressure - only once for all model run, unless thermodynamics included
        TAU_ast = plithost_int_node(corners_id(1,1), GCOORD, ...
            ELEM2NODE, Point_id, EL_IP.Rho_ip, PHY.G, true, SETTINGS.eids);
    end
    clear("corners_id")

    if SETTINGS.ini_deformation == 3                                                         % Flat topographies. 2D Gaussian decay temperature anomaly
        amp_t = PHY.t_seed_max * exp(-(GCOORD(2,:) - PHY.t_seed_y).^2 / PHY.t_sigma_y^2) ... % (ºC) [1,nnod7]
            .* exp(-(GCOORD(1,:) - PHY.t_seed_x).^2 / PHY.t_sigma_x^2);
        Temp = Temp + amp_t;                                                                 % (ºC) [1,nnod7] 
        Temp(Temp > PHY.temp_bc) = PHY.temp_bc;                                              % Temp0 = Temp; plot(Temp, GCOORD,ELEM2NODE)
        clear("amp_t")                                     
    end
      
    GEOn = getGEOn(GCOORD, ELEM2NODE, Point_id, Phases);                                     % GEO augmented with additional triangle-generated vertex nodes
    
    time0 = clock;
    etime = 0;
    remesh = 0;
    istep = 0;                                     % INTEGER
    ntime = 0.;                                    % REAL [s]
    isteps = 0;                                    % INTEGER [:] vector to record timestep indices [istep]
    ntimes = 0.;                                   % REAL [:]    vector to record times corresponding to "isteps": e.g. to ease exact timing of isochrons for irregular dt

else                                             % loadsave == TRUE : restart file
    if fname == "last"
        fnames = dir(DOUT_S_ROOT + "/" + CASE + ".0*.mat")
        fname = fnames(length(fnames)).name
        clear("fnames")
    end
    loadsave_file = fullfile(DOUT_S_ROOT, fname);

    pypath0 = pypath;
    load(loadsave_file)
    pypath = pypath0;

    remesh_wb = 0;
    time0 = clock;

    [SETTINGS, PHY, SOLVER, SP, CASEROOT, DOUT_S_ROOT, dsninput] = update_model(SETTINGS, PHY, SOLVER, SP, CASEROOT, DOUT_S_ROOT, dsninput); % branch from previous model
    CASE = split(CASEROOT,'/'); CASE = CASE(end);
    dsn = strrep(dsninput,'input','data');                                   % assumes similar directory tree among computers 
    dt = PHY.dt;

    if SETTINGS.HYDR.make
        if SETTINGS.HYDR.model && ~isfield(PHY.HYDR,'H2O_TABLES')            % if activated by update_model()
            PHY.HYDR.H2O_TABLES = load_H2O_table(SETTINGS.HYDR.H2O_table_file);
        end
    end
    
    if exist('MELTING') && isfield(MELTING,'make')                           % old MELTING struct - get rid of it and initialize
        MELTING = [];                                                              % struct for melting variables                                                                  
        MELTING.Crust_thickness = zeros(1,nsteps);                                 % [nsteps]
        MELTING.Trackp = [];                                                       % melting track points
        MELTING.tp_melt = false;  
    end
    if SETTINGS.FAUL.make
        if ~exist('FAULTS','var')
            FAULTS = [];
        end
        FAULTS = updateFAULTS(FAULTS);
    end
end % if loadsave

if SETTINGS.SERP.make
    if ~SETTINGS.HYDR.make
        error("set on modelled or parameterized hydrothermal flow to get the hydrothermal domain for serpentinization")
    end
    if ~exist("SERPEN")
        SERPEN.trackt = [];                                                 % tracking of serpentinization timesteps
        Dserp = zeros(1,nnod7);
    end 
end
save_step = sort([istep+1, 0:SETTINGS.save_by:(nsteps*3)]);

% MAIN LOOP
% =========

while ntime <= time_int
    % if ~exist("dtmax") && dt < PHY.dt
    %    dt = dt + 1./20. * (PHY.dt - dt);                                   % reduced dt grow back to nominal
    % end
    %if istep > 0
    %    error('ok')
    %end
    istep = istep + 1;
    ntime = ntime + dt;
    isteps = [isteps istep];
    ntimes = [ntimes ntime];
    
    disp(join(["istep: ",istep]));
    fprintf("ntime: %5.4f Myr\n", ntime/ma);
    fprintf("dt: %5.4f yr\n", dt/year);
    
    if exist("dtmax")                                                           % last saved step recovered
        save_step = sort([istep, save_step]);
        clear("dtmax")
    end
    
    % REMESHING
    switch SETTINGS.intersect_s % functions to deal with interface intersections and too thin layers
        %case 'dis_layers'
        %    error("dis_layers() has not been checked with left-right GEO")            % miguel also put a note that this is not currently working yet
        %    [GEOMETRY,Geo_id,Sgap,Intersect_ID,remesh] = ...
        %      dis_layers(GEOMETRY,Geo_id,Intersect_ID,SETTINGS.layshift,remesh);
      case 'layer_collapse'    
        if SP.make
            Base = Basement;                                                 % figure(); plot(Basement(1,:),Basement(2,:)); axis([])
        else                                                  
            Base = [];                                                       % figure(); plot_meshF(ELEM2NODE, GCOORD); hold on; plotGEO(GEOn); hold on; plotGEO(GEOlc, 1:10, 'o')
        end
        
        GEOlc = layer_collapse(GEOn, Base, SETTINGS);                        % potentially collapsed geometry
        if ~isequal([GEOlc.coo],[GEOn.coo]) && istep > 1                     % [GEOn.n; GEOlc.n]'; figure(); plot_meshF(ELEM2NODE, GCOORD)
            GEOn = GEOlc;                                                    % figure(); plot_phasesF(ELEM2NODE, GCOORD, Phases, 1:3, GEOn); hold on; plot_meshF(ELEM2NODE, GCOORD);
            remesh = true;                                                   % hold on; plot(Basement(1,:)/1000, Basement(2,:)/1000,'-','color',[.5 0 0])
        end                                                                  % hold on; plotGEO(GEOn); hold on; plotGEO(GEOlc, 1:9, 'o')
        clear("GEOlc")
      case 'none'                                                            % no action taken
      otherwise
        error("SETTINGS.intersect_s not understood") 
    end % switch SETTINGS.intersect_s

    track_remesh = [track_remesh remesh];
    if remesh                                                                        % REMESH
        
        %if (istep == 1153) % debugging
        %    GEO0 = GEOn;
        %    GCOORD0 = GCOORD; ELEM2NODE0 = ELEM2NODE; Point_id0 = Point_id; Phases0 = Phases;   % only for debugging    
        %    figure(); plot_phasesF(ELEM2NODE, GCOORD, Phases, 1:3); hold on; plotGEO(GEOn);  hold on; plot_meshF(ELEM2NODE, GCOORD)
        %end
        
        disp('main:: remeshing...')   
        [GCOORD, Point_id, Phases, ELEM2NODE, GEO, GEOn, ndof, nnod7, nel, ...
        F_xx, F_xy, F_yx, F_yy, TAU_xx_old, TAU_xy_old, TAU_yy_old, E2all, Mu_all, ... % [nel,nip] variables
        I2, SS, RHEOLvar, ...                                                          % cell/struct with [nel,nip] variables
        Temp, PRESS_CONT, Dpl, Dserp] = ...                                            % [1,nnod7] variables 
        remesh_all3l_dpl_rand_melt_I2_F_Mu_E2_linear(GCOORD, Phases, ELEM2NODE, GEOn, ...  
        Elsizes, nip, SETTINGS, PHY, "triangle", triangle_mode, fullfile(CASEROOT,"mesh"), ...
        F_xx, F_xy, F_yx, F_yy, TAU_xx_old, TAU_xy_old, TAU_yy_old, E2all, Mu_all, ... % [nel,nip] variables
        I2, SS, RHEOLvar, ...                                                          % cell/struct with [nel,nip] variables
        Temp, PRESS_CONT, Dpl, Dserp);                                                 % [1,nnod7] variables 
        nnod6 =  max(ELEM2NODE(1:6,:),[],'all');
        remesh = false;
        %  figure(); plot_phasesF(ELEM2NODE, GCOORD, Phases, 1:3, GEOn); hold on; plotGEO(GEOn,[3,6,9],'.',15); hold on;  plot_meshF(ELEM2NODE, GCOORD);
        %  hold on; plot_point_idF(GCOORD0, Point_id0, [1,2,3,4,5,6,7,8,9,10],'.',15);

        if PLOT.make
            PLOT.Point_id = Point_id;
        end
        
        if SP.make                                                                               % if active surface processes REMESH ISOCHRONS
            [ISOCHRONS,Basement] = resample_isoc(GCOORD,Point_id,ELEM2NODE, ...
                ISOCHRONS,Basement, SP.tp_isoc, SETTINGS.hres(end)/2);
        end
        
        TRACKP = resample_line(TRACKP,SETTINGS.res_tp,'spline');                                              % REMESH LAB
    end % if remesh
    
    [Topography,nodtopids] = find_topo(GCOORD,ELEM2NODE,Point_id);       % GCOORD(:,nodtopids) == Topography
    if ~isequal(nodtopids(1:2:end)',GEOn(end-1).pids)
        error('Point_id & GEOn not compliant: check last remeshing')
    end

    if ((SETTINGS.MELT.make && MELTING.tp_melt) || any(Dserp > 0.0)) && SETTINGS.OC.make
        [ocrust, ocrust_melt, ocrust_serp] = getOceanicCrust(GCOORD, ELEM2NODE, GEOn, MELTING, Dserp, PHY.OC, PHY.ext_rate);                    % hold on; plot_dikesF(MELTING, PHY, ntimes);
    end
    
    EL_IP.Rho_ip = dens_temp_Dpl(GCOORD, ELEM2NODE, Phases, GEOn, Basement, ...
        Temp, Dpl, nip, PHY, SP, ocrust);                                       % [nel, nip=6] | figure(); plot_ip_valF(EL_IP.Rho_ip, GCOORD, ELEM2NODE, GEOn(9).coo, [], "EL_IP.Rho_ip");
    
    % Damage weak seed
    if SETTINGS.ini_deformation == 5 && istep*dt < PHY.WS.time
        I2.f = damage_seed(I2.f,PHY.WS,GCOORD,ELEM2NODE,nip);
    end
    
    LOAD = [];
    LOAD.F_ext = zeros(nnod7*ndof,1);                                                % [Pa.m] terms directly added to the Rhs in the mechanical
    
    % MECHANICAL RHS FORCES BY INTEGRATION OF WINKLER & OCEAN PRESSURES
    if contains(PHY.bc_t,"winkler")                                                             % Bottom upward force [integrated pressure] (for Winkler boundary conditions)
        PwinklerG = getWinklerPressure(GCOORD, ELEM2NODE, Point_id, TAU_ast, PHY.add_atmPres);  % Galerkin mass matrix integrations              nodwin = abs(sum(PwinklerG.inn)) > 0
        LOAD.F_ext = LOAD.F_ext + PwinklerG.inn(:);                                             % [nnod7*2,1] [N] as column vector            plot(GCOORD(1,nodwin), PwinklerG.inn(1,nodwin))
    end
    if PHY.rho_w > 0. || PHY.add_atmPres % Top downward water load
        Water_thick = Topography;                                             % Water thickness
        Water_thick(2,:) = max(0., sea_level - Water_thick(2,:));
        PwaterG = getHydrostaticPressure(GCOORD, ELEM2NODE, Point_id, sea_level, PHY.rho_w, PHY.G, PHY.add_atmPres); % Galerkin mass matrix integrations: (\int N'NP) [Newtons]
        LOAD.F_ext = LOAD.F_ext + PwaterG.inn(:);                           % [nnod7*2,1] [N]    as column vector      -> RHS in mechanical
        LOAD.el = PwaterG.loadel;                                           % LOGICAL [1,nel]                          -> free surface stabilisation
        LOAD.rho = PwaterG.rhoel;                                           % [1,topnel], where topnel = sum(LOAD.el)  -> free surface stabilisation
        LOAD.rho = 0.;
    else 
        LOAD.el = zeros(1,size(ELEM2NODE,2))==1;                            % [1,nel]
        LOAD.rho = 0.;
    end % if PHY.rho_w > 0. else...
  
    % set thermomechanical boundary conditions
    [Bc_ind, Bc_val, Bc_ind_fs, Bc_val_fs, ext_srate] = ...
        set_bcs_vel(GCOORD, Point_id, PHY.ext_rate, PHY.bc_t, dt);
    [Bct_ind, Bct_val] = set_bcs_temp(GCOORD, Point_id, PHY.temp_bc, PHY.temp_bc_depth, PHY.Ttop);
  
    % update velocities & pressure
    switch SOLVER.iteration_type
      case 'MILAMIN-penalty'
        [Vel, Pressure, PRES_IP, TAU_xx, TAU_yy, TAU_xy, STRAIN_xx, STRAIN_yy, STRAIN_xy, ...                                    % O: Vel: [2,nnod7] [m/s], Pressure: [3,nel]
        Gamma,Yield_T2,YC,E2all,Mu_all,Mu_dis_all,Mu_dif_all,Mu_b_all,Dx_e,F_xx,F_xy,F_yx,F_yy,I2,GIP_x_all,GIP_y_all, ...       % IO: E2all, Mu_all, F_xx, F_xy, F_yx, F_yy, I2, THETA_all, nw_it
        THETA_all,W_xy,nw_it,residue] = mechanical2d_m(ELEM2NODE, Phases, GCOORD, Point_id, Temp, E2all, ...          % Temp: [nnod7,1] | E2all, Mu_all, Rho_ip: [nel,nip] | RHEOLvar: cell with [nel,nip] elements
        Mu_all, PHY.RHEOL, RHEOLvar, PHY.Phi, PHY.Cohesion, PHY.R, PHY.Shearm, EL_IP.Rho_ip, PHY.G, Bc_ind,Bc_val,... % 
        nip, reorder, ext_srate, SETTINGS.top_surface, SETTINGS.fs_alpha, SETTINGS.fs_beta, dt, ...                   % hold on; plot_finiteStrainEllipses(GCOORD, ELEM2NODE, F_xx, F_xy, F_yx, F_yy)
        Bc_ind_fs, Bc_val_fs, PHY.bc_t, SS, PHY, F_xx, F_xy, F_yx, F_yy, I2, THETA_all,...                            % hold on; plot_deviatoricStress(GCOORD, ELEM2NODE, TAU_xx, TAU_xy, TAU_yy,"black")
        TAU_xx_old, TAU_yy_old, TAU_xy_old, SETTINGS.is_elastic, nw_it, ...                                           % hold on; plotGEO(GEOn,1:10,'')
        LOAD.F_ext, LOAD.rho, LOAD.el, PHY.PLASTICITY, NUMSCALE, PLOT, SOLVER);
        STRAIN_xx_ip = STRAIN_xx;
        STRAIN_xy_ip = STRAIN_xy;
        STRAIN_yy_ip = STRAIN_yy;
      case 'Newton'
        milamin2newtsolv;
        [U,P,MU,residue,ER,TAU,I2,SS] = ...
        mechanical2dN(MESH, SOLVER, YIELD, PHY.RHEOL, RHEOLvar, SS, PHY, nnodel, nip, EL_IP.Rho_ip, g, PHY.R, dP_dz, ...
        ifix, vfix, ELASTICITY, NUMSCALE, U, P, Temp, E2all, TAU, I2, FSSA, dt, LOAD, PLOT);
        newtsolv2milamin;
    end % switch SOLVER.iteration_type
    
    dtmax = get_bottom_dtmax(GCOORD,ELEM2NODE, Vel, dt, Ylim(1));        % max dt for triangles not fully crossing  initial domain bottom
    
    try % dtmax
        if dt > dtmax
            error('dt > dtmax: go one saved step back')
        end
        clear("dtmax")
    catch
        fnames = dir(DOUT_S_ROOT + "/" + CASE + ".0*.mat");
        load(fullfile(DOUT_S_ROOT, fnames(length(fnames)).name));       % get last one
        clear("fnames")
        dt = dtmax;
        fprintf('Time step decreased to %5.4f Myr\n', dt/ma);
        continue;
    end % try...catch dtmax
    
    if istep ~= 1
        PRESS_CONT_old = PRESS_CONT;                                                    % [1,nnod7] store previous iteration pressure for melt_prod_onestep()
    end
    PRESS_CONT = reshape(call_el2nod_pressure(GCOORD, ELEM2NODE, nel, Pressure),1,[]); % [1,nnod7] <- [3,nel] discontinuous to continuous pressure
    if istep == 1
        PRESS_CONT_old = PRESS_CONT;
    end
    [ErP,~] = getPVstrainrateII(TAU_xx, TAU_yy, TAU_xy, YC, Yield_T2, Gamma, Mu_dif_all, Mu_dis_all); % figure(); plot_ip_valF(ErP, GCOORD, ELEM2NODE, Point_id, [0,3.5E-14]);
    [BDTdepth(istep), BDTat(:,istep)] = getBDTdepth(ErP, E2all, Temp, GCOORD, ELEM2NODE, GEOn, SETTINGS.eids, PHY);
    
    % MELTING: generation of magma
    if SETTINGS.MELT.make
        [Dpl, Temp, MELTING.dF, dFaxis, dFaxisfac] = melt_prod_onestep(Dpl, Temp, PHY.MELT, ...       % figure(); plot_tF(dF, GCOORD, ELEM2NODE, [], "\delta Depletion", ntime); hold on; plotGEO(GEOn,1:10,'')
            PRESS_CONT, GCOORD, ELEM2NODE, Phases, GEOn, E2all, nodtopids);                   % figure(); plot_tF(dFaxis, GCOORD, ELEM2NODE, [], "\delta depletion @ spreading center", ntime); hold on; plotGEO(GEOn,1:10,'')
        % hold on; plot(BDTat(1,istep)/1000,BDTat(2,istep)/1000)
        if any(dFaxis > 0.)  
            [MELTING.area_melt_axis, ~] = integrate2D(GCOORD, ELEM2NODE, dFaxis);                % [m2] total integral of produced melting feeding the spreading center magma emplacement
        else
            MELTING.area_melt_axis = 0.;
        end
    end % melting
    
    % SHEAR HEATING
    EL_IP.Hs = shear_heat(TAU_xx, TAU_yy, TAU_xy, TAU_xx_old, TAU_yy_old, ...   % [nel, nip]
    TAU_xy_old, STRAIN_xx_ip, STRAIN_yy_ip, STRAIN_xy_ip, Phases, PHY.Shearm, dt);
    SH_F = repmat(PHY.sh_factor(Phases), 1, nip);
    EL_IP.Hs =  SH_F .* EL_IP.Hs;

    % material properties at quadrature points [EL_IP{.K, .Cp. .Hp}] and Phases_ip & Phases_el indicators with oceanic crust and sediment as independent phases
    [EL_IP, ~, Phases_ip, Phases_el] = updateMaterial(EL_IP, GCOORD, ELEM2NODE, GEOn, Phases, ocrust, Basement, ...
        SETTINGS, PHY, SP);
    
    if SETTINGS.HYDR.make         % TRUE also if serpentinite is active
        if ~hydroactive
            if getDzminCrust(GEOn) <= PHY.HYDR.dzcrustmax
                save_force = true;
                hydroactive = true;                                                   % initiate and stay active hereafter
                HVAR = [];                                                            % dynamic variables for the hydrothermal model
            end
        end
        if hydroactive
            hydroelboo = getHydroDomain(GCOORD, ELEM2NODE, Point_id, Temp, ErP, PHY.HYDR, sea_level, SETTINGS.HYDR.model);      % LOGICAL [1,nel] - inherit Temp from global domain
        end
    end
    % SERPENTINIZATION
    Qserp_nodal = [];
    if SETTINGS.SERP.make && hydroactive && any(hydroelboo)
        %SERPEN.elboo = getSerpDomain(GCOORD, ELEM2NODE, Phases, hydroelboo, Temp, ErP, PHY.SERP);  % LOGICAL [1,nel] elements potentially subject to serpentinization
        SERPEN.srfac =  getSerpStrainrateFactor(ELEM2NODE, Phases, hydroelboo, ErP, PHY.SERP, SETTINGS.eids); 
        if any(SERPEN.srfac > 0.0) % figure(); plot_meshF(ELEM2NODE, GCOORD, hydroelboo, [0 0 0.8]); hold on; plot_meshF(ELEM2NODE, GCOORD, SERPEN.elboo, [0 .6 0], 2); hold on; plotGEO(GEOn)
            [dFserp, Dserp] = serp_underplate(GCOORD, ELEM2NODE, Temp, Dserp, dt, ...         % {dFserp, Dserp} in [0,1]: increment and updated serpentinization ratios
              MELTING.Trackp, SERPEN.srfac, PHY.SERP, istep);   % figure(); plot_tF(Dserp, GCOORD, ELEM2NODE, [], "serpentinization [0/1]", ntime); hold on; plotGEO(GEOn)
            % dTemp_serp + PHY.SERP.Hserp .* dFserp; % only for diagnostic | figure(); plot_tF(PHY.SERP.Hserp .* dFserp, GCOORD, ELEM2NODE, [], "dT by serpentinization [ºC]", ntime); hold on; plotGEO(GEOn)
            Qserp_nodal = PHY.SERP.Hserp_heat * PHY.SERP.olivine_density * dFserp;            % [J.m-3] = [J.kg-1].[kg.m-3].[0/1]
            if sum(Dserp > 0.05) >= 20
                SERPEN.trackt =[SERPEN.trackt ntime];
            end
            % H2 production in the Olivine -> serpentinization reaction
            olivine_to_H2 = PHY.SERP.olivine_density * PHY.SERP.olivine_H2_molar_ratio / PHY.SERP.olivine_molar_mass * 1000; % [mol.m-3] moles of H2 per m3 of serpentinized rock 
            SERPEN.serp.total(istep) = integrate2D(GCOORD, ELEM2NODE(1:6,:), Dserp(1:nnod6));             % [m3] total volume of serpentinized olivine
            SERPEN.serp.rate(istep) = integrate2D(GCOORD, ELEM2NODE(1:6,:), dFserp(1:nnod6)) / dt * year; % [m3.yr-1] volumetric rate [olivine rock] of serpentinization
            SERPEN.H2p.total(istep) = integrate2D(GCOORD, ELEM2NODE(1:6,:), Dserp(1:nnod6))  * olivine_to_H2;                % [mol] H2 accumulated in time
            SERPEN.H2p.rate(istep)  = integrate2D(GCOORD, ELEM2NODE(1:6,:), dFserp(1:nnod6)) * olivine_to_H2 / dt * year;    % [mol.yr-1]
        end
    end % serpentinisation
  
    % b] TEMPERATURE DIFFUSION + HEAT RELEASE FROM IGNEOUS BODIES
    if isfield(EL_IP,'K_hyd') 
        EL_IP = rmfield(EL_IP,'K_hyd');
    end
    
    [MELTING.Crust_thickness, MELTING.Trackp, MELTING.dT_sensible, MELTING.ign_body, MELTING.tp_melt, MELTING.dE_latent] = emplacementDikeSill(... % figure();plot(ntimes(1:end-1)/ma,MELTING.Crust_thickness(1:length(ntimes(1:end-1))),'o'); xlabel('Time [Myr]')
        MELTING.Crust_thickness, MELTING.Trackp, Temp, MELTING.area_melt_axis, GCOORD, ELEM2NODE, GEOn, Phases, PHY.ext_rate, dt, ...
        E2all, ErP, MELTING.dF, Dserp, PHY.MELT, SETTINGS.MELT, istep, Bct_ind, EL_IP, SOLVER.nelblo, BDTat(:,istep), ocrust);
        
    % HYDROTHERMAL CIRCULATION
    if SETTINGS.HYDR.make
        if hydroactive                                                                                        % figure(); plot_ip_valF(E2all, GCOORD, ELEM2NODE, GEOn(9).coo, [0,3.5E-14]);
            if any(hydroelboo)                                                                                                    % hold on; plot_meshF(ELEM2NODE, GCOORD, hydroelboo, [.8 0. 0.])
                if SETTINGS.HYDR.model                                                                        % figure(); plot_ip_valF(E2all-ErP, GCOORD, ELEM2NODE, GEOn(9).coo, []);
                    error("Hydrothermal component not yet public. Set SETTINGS.HYDR.model = false for a parameterised hydrothermal effect")
                else                             % ~SETTINGS.HYDR.model :: increase thermal conductivity for parameterised hydrothermal flow
                    disp("hydrothermal-par cooling active at " + sum(hydroelboo) + " elements")
                    Htcfac = getHThermalCondFac(ELEM2NODE, PHY.HYDR.PARAM, Temp, ErP, hydroelboo);     % [nel,nip] thermal conductivity factor due to [parameterised] hydrothermal flow
                    EL_IP.K_hyd = EL_IP.K .* Htcfac;                                                       % [nel,nip] thermal conductivity
                end % ~SETTINGS.HYDR.model                                                             % figure(); plot_ip_valF(Htcfac, GCOORD, ELEM2NODE, GEOn(9).coo); hold on; plot_meshF(ELEM2NODE,GCOORD, hydroelboo, 'red')                                                  
            end % if any(hydroelboo)                                                                   % hold on; plotGEO(GEOn); hold on; plot_meshF(ELEM2NODE, GCOORD, hydroelboo, 'red')
        end % if hydroactive
    end % if SETTINGS.HYDR.make
    
    % TEMPERATURE DIFFUSION SOLVER
    dQdt_qp = EL_IP.Hp + EL_IP.Hs;
    if ~isempty(Qserp_nodal) % exothermic heat release by serpentinization
        dQdt_qp = dQdt_qp + nodval2ipval(ELEM2NODE(1:6,:), Qserp_nodal, nip) / dt;      % [W.m-3] figure(); plotN(Qserp_nodal/dt, GCOORD, ELEM2NODE)
    end                                                                                 
    if ~isempty(MELTING.dE_latent) % latent heat of magma crystallization
        dQdt_qp = dQdt_qp + nodval2ipval(ELEM2NODE(1:6,:), MELTING.dE_latent, nip) / dt; % [W.m-3];
    end
    Temp_a = thermal2d_m(ELEM2NODE(1:6,:), GCOORD, Temp, ...                          % Temp: [1,nnod6], {'Cps','Hps'}: [nphases,1]
        EL_IP.Rho_ip .* EL_IP.Cp, EL_IP.K, dQdt_qp, ...                                 % [nel,nip] variables       nnod6 =  max(ELEM2NODE(1:6,:),[],'all');
        Bct_ind, Bct_val, dt, reorder, SOLVER.nelblo); 
    Temp_a = addCentroidVal(ELEM2NODE(1:6,:), Temp_a);                                % [1,nnod7] - for feedback into mechanical & to make it usable by current melting routines
        
    if SETTINGS.HYDR.make && hydroactive && any(hydroelboo)            % parameterised cooling by enhanced diffusion                                                         % parameterized hydrothermal influence of thermal conductivity active                 
        Temp_hyd = thermal2d_m(ELEM2NODE(1:6,:), GCOORD, Temp, ...                  % Temp: [1,nnod6], {'Cps','Hps'}: [nphases,1]
          EL_IP.Rho_ip .* EL_IP.Cp, EL_IP.K_hyd, dQdt_qp, ...                         % [nel,nip] variables       nnod6 =  max(ELEM2NODE(1:6,:),[],'all');
          Bct_ind, Bct_val, dt, reorder, SOLVER.nelblo);
        Temp_hyd = addCentroidVal(ELEM2NODE(1:6,:), Temp_hyd);                          % figure(); plot_tF(Temp_hyd, GCOORD, ELEM2NODE, [], "hydrothermal param Temp [ºC]", ntime)
        dTemp_hyd = Temp_hyd - Temp_a;                                                  % [ºC] increment on thermal state due to parameterized hydrothermal activity with respect to non-hydrothermal activity
        % dtem_hydr_int = integrate2D(GCOORD, ELEM2NODE, dTemp_hyd);                    % figure(); plot_tF(dTemp_hyd, GCOORD, ELEM2NODE, [], "hydrothermal influence on Temp [ºC]", ntime);
        Temp = Temp_hyd;    
    else % no hydrothermal effect
        Temp = Temp_a;
    end
    clear('Temp_a');
    if ~isempty(MELTING.dT_sensible)
        Temp = Temp + MELTING.dT_sensible; % delta T by sensible heat subject to their independent heat diffusion
    end
    
    % MOVE TRACKED POINTS
    if SETTINGS.tracking
        TRACKP = trackPoints(GCOORD, ELEM2NODE, GEOn, TRACKP, Vel, dt);
    end
    
    if SP.make                                                                                     % surface processes active
        if SP.tp_isoc                                                                              % track isochrons            
            ISOCHRONS(1:2,:) = trackPoints(GCOORD, ELEM2NODE, GEOn, ...
            ISOCHRONS(1:2,:), Vel, dt);
        end
        Basement = trackPoints(GCOORD, ELEM2NODE, GEOn, Basement, Vel, dt);
    end
    
    if SETTINGS.ini_deformation == 5
        WScoord = trackPoints(GCOORD, ELEM2NODE, GEOn, WScoord', Vel, dt)';
    end
    
    if SETTINGS.MELT.make && MELTING.tp_melt
        MELTING.Trackp = trackPoints(GCOORD, ELEM2NODE, GEOn, MELTING.Trackp, Vel, dt);
    end
        
    % UPDATE & KEEP OLD COORDINATES
    GCOORD_old = GCOORD;
    GCOORD = GCOORD + Vel*dt;                                                % break GCOORD <-> GEOn synchrony 
    
    % MAKE STRAIGHT EDGES
    eids = [6 4 5]; % <- rift2ridge2D | eids = [4 5 6] <- kinedyn
    GCOORD(:,ELEM2NODE(eids,:)) = 0.5*(GCOORD(:,ELEM2NODE([1 2 3],:)) + GCOORD(:,ELEM2NODE([2 3 1],:)));
    GCOORD(:,ELEM2NODE(7,:))    = 1/3*(GCOORD(:,ELEM2NODE(1,:))+GCOORD(:,ELEM2NODE(2,:))+GCOORD(:,ELEM2NODE(3,:)));
    integrate2D(GCOORD, ELEM2NODE, ones(1,nnod7));                           % area of domain [m2] : just to catch mesh errors
    
    % STRESS ROTATION :: transform stress in the elastically rotating particles to the standard xy axes
    [TAU_xx_old,  TAU_yy_old, TAU_xy_old] = rotateSymmetricTensor2D(TAU_xx, TAU_yy, TAU_xy, -W_xy*dt); % where W_xy is the angular velocity of the rotating particle
    if contains(PHY.bc_t,"winkler")
        [GCOORD, Temp, F_xx,F_xy,F_yx,F_yy,I2,TAU_xx_old,TAU_xy_old, ...
        TAU_yy_old,Mu_all,E2all,remesh] = update_bot_values_linear ...
        (GCOORD,ELEM2NODE,PHY.ymin1,Point_id, Temp, PHY.temp_bc, F_xx, ...
        F_xy,F_yx,F_yy,I2,TAU_xx_old,TAU_xy_old,TAU_yy_old,Mu_all, ...
        E2all, PHY.R, Phases, remesh);
    end
    if min(GCOORD(1,Point_id==9)) ~= min(GCOORD(1,Point_id==10))
        error("---ERR 001---")
    end    

    % SURFACE PROCESSES
    [Topography,nodtopids] = find_topo(GCOORD,ELEM2NODE,Point_id); % find topography after mechanical displacements
    if SP.make % apply surface processes
        % TODO: include time-varying pelagic rate either as input, separated from SP or calculated within surf_proc()
        New_topo = surf_proc(Basement, Topography, GCOORD, SP, dt, zuco);                     % default: diff_topo_sealand()
        % Update the topography and interpolate variables of reshaped elements by surface processes
        [GCOORD, Temp, F_xx, F_xy, F_yx, F_yy, I2, TAU_xx_old, TAU_xy_old, TAU_yy_old, ...
        Mu_all, E2all, remesh] = update_topo_values_linear( ...
            GCOORD, ELEM2NODE, New_topo, nodtopids, Temp, PHY.Ttop, F_xx, F_xy, F_yx, F_yy, I2, ...
            TAU_xx_old, TAU_xy_old, TAU_yy_old, Mu_all, E2all, PHY, PHY.R, Phases, remesh);
            
            [Topography,nodtopids] = find_topo(GCOORD,ELEM2NODE,Point_id); % find topography after mechanical displacements
            % Update isochrons
        if floor(istep/SP.iso_sav) == (istep/SP.iso_sav)
            ISOCHRONS = [ISOCHRONS [Topography; ...
            istep*ones(1,size(Topography,2)); ...
            New_topo(2,:)]];
        end
        % Update (erode) basement
        Topo_base = interp1(Topography(1,:),Topography(2,:),Basement(1,:));
        Basement = [Basement(1,:); ...
                    min([Basement(2,:); Topo_base])];
        clear("Topo_base")
    end % if SP.make
    
    for i=1:length(GEOn)
        GEOn(i).coo = GCOORD(:,GEOn(i).pids);                                    % re-establish GCOORD <-> GEOn synchrony 
    end
    if ~isequal(minmax(GEOn(end-1).coo(1,:)),minmax(GCOORD(1,:)) )
        error("GEOn - GCOORD synchronization failed")
    end

    % start previous BDT...SERP block
    if SETTINGS.MELT.make                                                           % for next mechanical loop    
        if size(PHY.RHEOL.Ndis,2) > 1                                               % update viscosity values according to depletion
            RHEOLvar = rheol_var_handle(Dpl, ELEM2NODE, Phases, ...                 % 1x2 cell array, each: [nel,nip]
            SOLVER.nip_stress, GCOORD, false);        
        end
    end
    % end previous BDT...SERP block

    %CHECK MESH INTEGRITY
    [qn,amin,amax] = checkmesh(GCOORD,ELEM2NODE);
    areaq = el_area(GCOORD,ELEM2NODE);
    if (qn<0.2 || amin<7 || amax>170 || any(areaq<=0))
        remesh = true;
    end
    
    etime = etime  + dt;
    
    % CHECK BREAK UP
    if ~BREAKUP.bool
        BREAKUP.Basement   = Basement;
        BREAKUP.sp_s       = SP.make;
        BREAKUP = break_up(BREAKUP, SETTINGS.BRKU, GEOn, ntime);
        if BREAKUP.bool
            SETTINGS.layshift = SETTINGS.layshifts(2);
        end
    end

    %======================================================================
    % SAVING THE COORDINATES OF THE MESH, THE TEMPERATURE, THE VELOCITY, 
    % TRACERS, SOME VARIABLES AND TIME SPENDING IN THE PROCESSING
    %======================================================================
    time = clock - time0;
    if SETTINGS.save_part
        s_name = join([CASE,"p",num2str(istep,'%0.5d'),"mat"],".");
        save(fullfile(DOUT_S_ROOT,s_name), SETTINGS.sav_var_part{:});
    end
    if ismember(istep, save_step)
        s_name = join([CASE,num2str(istep,'%0.5d'),"mat"],".");
        save(fullfile(DOUT_S_ROOT,s_name));
    end
    
    %======================================================================
    % POSTPROCESSING
    %======================================================================
    
    if PLOT.make            
        clf
        drawnow
        figure(1)
        log_data = 0;
        plot_eri
        hold on
        %plot_isoc
        axis([-100 100 -40 5])
        drawnow
        hold off    
    end
    
    %SWITCH FOR THE TIME STEP CRITERIA
    if SETTINGS.top_surface == "free"
      switch tstep_type
        case 'cnst'                      % the only working one
        case 'courant'                   % Courant criteria (currently not working)
          dt = C*min(Dx_e)/max(Vel,[],'all');
      end
    end
    % fprintf(1, [num2str(toc,'%8.6f'),'\n']);

end % while ntime <= time_int
