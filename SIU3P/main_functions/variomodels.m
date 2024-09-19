function val = variomodels(h, c0, cmax, range, vartype)
  a = 3;             % as Chiles&Delfiner 1999
  switch vartype
    case 'exp'                                                      % exponential: 0.95*cmax at range
      val = c0 + (cmax - c0) * (1 - exp(-a*h/range));
    case 'sph'
      val = c0 + (cmax - c0) * (1.5*h/range - 0.5*(h/range).^3);    % spherical: cmax at range
      val(h>range) = cmax; 
    case 'gau'
      val = c0 + (cmax - c0) * (1 - exp(-a*(h/range).^2));
    otherwise
      error("variomodels: not understood")  
  end % switch 
end % function
