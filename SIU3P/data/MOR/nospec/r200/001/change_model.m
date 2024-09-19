function [SETTINGS, PHY, SOLVER, SP] = change_model(SETTINGS, PHY, SOLVER, SP)
  % see param_defaults() for definitions and options 

  km      = 1000.;                                                        % [m]
  year    = 365.25*24*60*60;                                              % [s]
  ma      = 1.e6*year;                                                    % [s]

  PHY.dt  = 200*year;                                                     % [s] reference time step

  PHY.K    = [3.3;      3.3;            3.3];                             % [W.m-1.K-1] thermal conductivity [mantle, LC, UC]
  PHY.Hp   = [0.;        0.;            0.];                              % $H_r$ radiogenic heat production      ./[1000; 1000; 1000]./Rho;
  PHY.Ce   = [3.;        3.;            3.] * 1.e-5;                      % $\alpha_T$ thermal expansion coefficient  [mantle, LC, UC]

  PHY.RHEOL.type = "all_peridotite";                                      % used by load_rheol()
  PHY.Shearm = [74.e9; 74.e9; 74.e9];

  PHY.SS.I2_phi = [0. 0.3];                                               % Zhonglan values for Friction angle reduction
 
  PHY.Dens = [3300.; 3300.; 3300.];                                       % Initial density at each phase

  %% BOUNDARY CONDITIONS
  %%      - Temperature                                               
  %----------------------------------------------------------------------
  PHY.temp_bc_depth = -40.*km;                                            % 1) BC temp: depth [m]
  PHY.temp_bc = 1300.;                                                    % 2) BC temp: value [ºC]
  PHY.temp_bc_depth_ini = -15.*km;                                        % depth BC for the initial conditions
  PHY.temp_bc_ini = 1300.;                                                % temperature BC for the initial conditions

  %% INITIAL CONDITIONS
  %%      - Geometry                                                  
  %----------------------------------------------------------------------
  PHY.Lx_model = 80.;                                                     % [km] initial domain full width
  PHY.xmin =  -PHY.Lx_model/2 * km;
  PHY.xmax =   PHY.Lx_model/2 * km;
  PHY.ymin1 = -40.*km;                           % layer 1 [mantle]      - bottom [== bottom of domain]
  PHY.ymax1 =  -8.*km;                           % layer 1 [mantle]      - top [== layer.2 bottom == zmoho] 
  PHY.ymax2 =  -4.*km;                           % layer 2 [lower crust] - top [== layer.3 bottom]
  PHY.ymax3 =   0.*km;                           % layer 3 [upper crust] - top
  PHY.zlab  = -15.*km;                           % lithosphere-asthenosphere-boundary
 
  % - Temperature anomaly as week seed initializer of the deformation with spatial Gaussian decay (ini_deformation = 3 and 4)
  PHY.t_seed_max = 50.;                                                        % (ºC) maximum temperature anomaly
  PHY.t_seed_x = mean([PHY.xmin,PHY.xmax]);                                    % (m)  maximum x (horizontal) location
  PHY.t_seed_y = -15.*km;                                                      % (m)  maximum y (depth) location 
  PHY.t_sigma_x = 5.*km;                                                       % (m) sdev for the x Gaussian decay
  PHY.t_sigma_y = 5.*km;                                                       % (m) sdev for the y Gaussian decay

  % MESH RESOLUTION
  resmul = 1.;
  SETTINGS.DWOM = resmul * 200.;                                        % [m] 'soft' resolution: Dry and wet olivine mantle
  SETTINGS.LCR  = resmul * 200.;                                        % [m] 'soft' resolution: Lower crust
  SETTINGS.UCR  = resmul * 200.;                                        % [m] 'soft' resolution: Upper crust

  % - Mechanical                                                
  PHY.ext_rate = 30. / 1000 / year;                                          % [m/s] <- [km/My == mm/yr] Full extension rate -  Zhonlang: 30 mm/year for his 80km wide model

  %% - surface process parameters
  SP.pelagic_rate = 0./year;                       % Pelagic sedimentation rate [m^2/s]
  SP.kdecay = 1.e-4;                               % Decay of the Hill-slope diffusion under the sea
  SP.K = 2.5;                                      % [W.m-1.K-1] Thermal conductivity for the sediments. If []; i.e. isempty(.) , equals that from the upper crust

  SETTINGS.res_tp = 200.;                          % [m] initial resolution of tracking point lines

end % function