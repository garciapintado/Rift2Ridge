function [SETTINGS, PHY, SOLVER, SP, CASEROOT, DOUT_S_ROOT, dsninput] = update_model(SETTINGS, PHY, SOLVER, SP, CASEROOT, DOUT_S_ROOT, dsninput)
  % see param_defaults() for definitions and options 

  km   = 1000.;                                                           % [m]
  year = 365.25*24*60*60;                                                 % [s]
  ma   = 1.e6*year;                                                       % [s]

  PHY.time_int = 100*ma;
  
end % function
