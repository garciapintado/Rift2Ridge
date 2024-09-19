function V = addCentroidVal(EL2NOD, V)
  % add centroid value to a 3-/6-nnodel element by interpolating from nodes
  %
  % EL2NOD :: [3,nel] or [6,nel]
  % 
  % Author: Javier Garcia-Pintado, MARUM. 2020
  
  [nnodel, nel] = size(EL2NOD); 
  if ~ismember(nnodel, [3 6])
     error('addCentroidVal:: EL2NOD input requires nnodel in {3,6}')
  end
  if size(V,2) ~= max(EL2NOD(:))
      error('addCentroidVal:: V not compliant with EL2NOD')
  end
  
  [x_ip, ~] = ip_triangle(1);                                              % centroid as integration point
  [N, ~]    = shp_deriv_triangle(x_ip, nnodel);
  V = [V N{1}' * V(EL2NOD)];
 
end

