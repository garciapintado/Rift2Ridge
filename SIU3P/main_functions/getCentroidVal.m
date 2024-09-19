function V = addCentroidVal(EL2NOD, V)
  % add centroid value to a 3-/6-nnodel element by interpolating from nodes
   
  nnodel = size(EL2NOD,1); 
  nip = 1;                                                                 % centroid
  [x_ip, ~] = ip_triangle(nip);
  [N, ~]    = shp_deriv_triangle(x_ip, nnodel);       
  Nm = cell2mat(N)';
  V = [V Nm' * V(EL2NOD)];
 
end

