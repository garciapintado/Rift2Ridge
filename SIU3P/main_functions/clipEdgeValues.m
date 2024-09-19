function val = clipEdgeValues(EL2NOD, val, signs)
  % function val = clipEdgeValues(EL2NOD, val, signs)
  % +++ purpose+++
  % truncate values in edge triangular nodes to those in vertex nodes
  %
  % val   :: [nnod] vector or [nvar,nnod] matrix, with one variable per row
  % signs :: OPT, CHARACTER in {"-","+","-+"}. Defaults to "-+"
  
  % 2020-07 Javier García-Pintado, MARUM
  if nargin < 3
     signs = "-+"; 
  end
  
  [nnodel, nel] = size(EL2NOD);
  
  if size(val,2) == 1
      dotranspose = true;
      val = val';
  else 
      dotranspose = false;
  end
     
  for i=1:size(val,1)                                                      % every row represents a variable
      val_eln = reshape(val(i,EL2NOD),nnodel,nel);
 
      if contains(signs,"-")
          minval_eln = val_eln;                                            % [nnodel,nel]
          minval_eln(4,:) = min(val_eln(2,:),val_eln(3,:));
          minval_eln(5,:) = min(val_eln(3,:),val_eln(1,:));
          minval_eln(6,:) = min(val_eln(1,:),val_eln(3,:));   
          ids = val_eln < 0.0 & val_eln < minval_eln;
          val_eln(ids) = minval_eln(ids);
      end
      if contains(signs,"+")
          maxval_eln = val_eln;                                            % [nnodel,nel]
          maxval_eln(4,:) = max(val_eln(2,:),val_eln(3,:));
          maxval_eln(5,:) = max(val_eln(3,:),val_eln(1,:));
          maxval_eln(6,:) = max(val_eln(1,:),val_eln(3,:));   
          ids = val_eln > 0.0 & val_eln > maxval_eln;
          val_eln(ids) = maxval_eln(ids);
      end
      val(i,:) = nelval2nodval(EL2NOD, val_eln);
  end
  if dotranspose
      val = val';
  end
  
end % function
