function P = getWinklerPressure(GCOORD, ELEM2NODE, Point_id, Pbottom, plusatm)
  % P = getWinklerPressure(GCOORD, ELEM2NODE, Point_id, Pbottom, plusatm)

  atm2pascal = 101325.;                                                      % [Pa/atm]
  
  nel = length(ELEM2NODE);
  % find_bottom
  intid = 1; % bottom interface
  Surf_nodes = find(Point_id==intid);
  [TopoX,isort] = sort(GCOORD(1,Surf_nodes));                                % xsorted = x(isort) 
  TopoY = GCOORD(2,Surf_nodes);
  Topography = [TopoX; TopoY(isort)];                    
  Gcootopids = Surf_nodes(isort);

  P = GCOORD * 0.0;                                                        % Vector of pressures [dof,nnod7]

  nnoded = 3;                                                              % nodes per triangular edge
  nip = 3;                                                                 % k=5 exact polynomial
  Top_elboo = sum(ismember(ELEM2NODE, Gcootopids)) >= nnoded;              % Top element indexes with >=3 nodes (vertices + edge) on the surface
  topelids = find(Top_elboo);                                              
  [~, xordi] = sort(GCOORD(1,ELEM2NODE(1,Top_elboo)));
  topelids = topelids(xordi);                                              % figure(); plot_meshF(ELEM2NODE, GCOORD, topelids)                
                                                                           %  
  [eortho, ~] = getOrthoDpolyline(Topography(:,1:2:end), false);           % upward normal to the slope at each edge
  Press_n = Pbottom;                                                       % [Pa] scalar, uniform for all bottom nodes

  if plusatm
      Press_n = Press_n + atm2pascal;                                           % add 1 atm
  end
  % get pressure at each topographical element topo nodes
  topnel = length(topelids);
  topEL2NOD = repmat((1:nnoded)',1,topnel) + repmat((nnoded-1)*(0:topnel-1),nnoded,1); % INTEGER [nnoded,topnel]   
  Presx =  Press_n .* repmat(cosd(eortho),nnoded,1);   
  Presy =  Press_n .* repmat(sind(eortho),nnoded,1);

  [Px_inn, Px_ie] =  integrateMM1D(Topography, Presx);
  [Py_inn, Py_ie] =  integrateMM1D(Topography, Presy);
  
  P = [];
  P.inn = GCOORD * 0.0;
  P.inn(:,Gcootopids) = [Px_inn; Py_inn];
  P.ie = zeros(2,nel);
  P.ie(:,topelids) = [Px_ie; Py_ie];
  P.loadel = sum(abs(P.ie)) > 0.;
end % function
