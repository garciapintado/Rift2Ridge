function [SETTINGS, PHY, SOLVER, SP, CASEROOT, DOUT_S_ROOT, dsninput] = change_model(SETTINGS, PHY, SOLVER, SP, CASEROOT, DOUT_S_ROOT, dsninput)
  % COT:001 - see param_defaults() for definitions and options 

  km      = 1000.;                                                        % [m]
  year    = 365.25*24*60*60;                                              % [s]
  ma      = 1.e6*year;                                                    % [s]
  
  SETTINGS.MELT.heat_release = true;
end % function
