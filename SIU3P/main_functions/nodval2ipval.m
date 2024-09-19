function val_ip = nodval2ipval(ELEM2NODE, valnod, nip)
  % function ipval = nodval2ipval(ELEM2NODE, valnod)
  % get values at integration points [nel,nip] from values at nodes [1,nnod]
  %
  %
  % Details: the function is hardcoded for 2D quadratic shape functions with nip=6.
  % It can use either 6 or 7 nodes, provided ELEM2NODE and valnod arrays are compliant
  %
  % Javier GP, MARUM, 2020
  
  if nargin < 3
      nip = 6;
  end
  nnod =  max(ELEM2NODE,[],'all');
  if length(valnod) > nnod
      valnod = valnod(1:nnod);    
  end
  if length(valnod) ~= nnod
      error("nodval2ipval:: variable array length do not match ELEM2NODE maximum value")   
  end
  
  [nnodel, nel] = size(ELEM2NODE);
  [ipx, ~] = ip_triangle(nip);                                             % ipx:[nip,2], ipw:[1,nip]
  N = shp_triangle(ipx, nnodel);                                           % [nnodel,nip]
  
  val_elnod = reshape(valnod(ELEM2NODE), nnodel, nel);
  val_ip = (N' * val_elnod)';                                              % [nel,nip] transpose to comply MILAMIN convention 
end % function
