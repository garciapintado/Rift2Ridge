function P = getHydrostaticPressure(GCOORD, ELEM2NODE, Point_id, ...
    llevel, lrho, G, plusatm)
  % P = (GCOORD, ELEM2NODE, Point_id, llevel, lrho, G, plusatm)
  % Get pressure of a hydrostatic load at nodes on the top of a mesh.
  % The nodes include corner and edge ones in a triangular mesh.
  % If integrate_basis=false, the hydrostatic pressure is evaluated at the
  % top surface nodes, if true (hardcoded as default), then basis function are integrated at each
  % element, and the returned value is a cell containing the integrals for
  % each basis and the total integral per element
  % 
  % INPUT [I]
  % GCOORD   :: REAL [2,nnod] mesh coordinates
  % ELEM2NODE:: INTEGER [nnodel,nel] element definition
  % Point_id :: INTEGER [1,nnod] indicate node classes. The function
  %             assumes that max(Point_id)-1 defines the top nodes to which pressure is applied 
  % llevel   :: load top level
  % lrho     :: load density
  % G        :: [1,2] gravity acceleration constant
  % plusatm  :: LOGICAL; whether to add atmospheric pressure
  % integrate_basis :: LOGICAL; whether (false) to return the directly evaluated
  %                    pressure values or further (true) return the integrate
  %                    cuadratic basis functions at each node and element
  % 
  % OUTPUT [O]
  % P :: [2,nnod] for integrate_basis = false
  % P :: list with {'inn','ie'} integrations for integrate_basis = true
  %
  % Details:
  % The hydrostatic pressure [Pa] is evaluated at each node for each element, and 
  % than integrated over the element [Pa.m] and weighted by 1D edge shape functions 
  % compliant with the Galerkin 2D shape function weighting. Thus  
  %. The mathematical operation at each
  % topographical edge element is $\int_{\partial_e} N'NP \partial_e$, where 
  % $\partial_e$ represents the triangular edge, N [1,nnoded], is the 1D  
  % shape function end nnoded is the number of nodes per edge.
  % The number of nodes per triagular edge is nnodel=3, which is hardcoded 
  % in this function to match edge shape functions in nnodel{6,7} quadratic
  % Lagrange polinomial in 2D triangles.
  %
  %--------------------------------------------------------------------------
  % Author: Javier García_Pintado, MARUM, 2020-05-02
  %--------------------------------------------------------------------------

  atm2pascal = 101325.;                                                    % [Pa/atm]
  integrate_basis = true;
  
  nel = size(ELEM2NODE,2);

  % find interface nodes
  intid = max(Point_id)-1;                                                 % top interface
  Surf_nodes = cast(find(Point_id==intid),'int32');
  Topography  = GCOORD(:,Surf_nodes);
  [~,isort]   = sort(Topography(1,:));                                     % xsorted = x(isort) 
  Topography  = Topography(:,isort);                    
  Gcooids = Surf_nodes(isort);
  if any(diff(GCOORD(1,Gcooids)) <= 0.)
      error("getHydrostaticPressure:: ---ERR001---")
  end

  hthick = max(0,llevel - Topography(2,:));                                % [1,ntopo] [m] hydrostatic layer thickness
  Press_n =  hthick * lrho * (-G(2));                                      % [1,ntopo] [Pa] pressure module
  if plusatm
      Press_n = Press_n + atm2pascal;                                      % [1,ntopo] add 1 atm
  end
  [eortho, ~] = getOrthoDpolyline(Topography(:,1:2:end), true);            % [1,topnel] downward normal to the slope at each interface edge

  % get pressure at each topographical element topo nodes
  nnoded = 3;                                                              % nodes per triangular edge
  topnel = numel(eortho);
  topEL2NOD = repmat((1:nnoded)',1,topnel) + repmat((nnoded-1)*(0:topnel-1),nnoded,1); % INTEGER [nnoded,topnel]   
  Presx = Press_n(topEL2NOD) .* repmat(cosd(eortho),nnoded,1);            % (3,topnel)
  Presy = Press_n(topEL2NOD) .* repmat(sind(eortho),nnoded,1);            % (3,topnel)    
  
  % integrate_basis
  % nip = 3;                                                               % k=5 exact polynomial
  %[Px_inn, Px_ie] =  integrateBasis1D(Topography, Presxy(1,:), nnodel, nip);
  %[Py_inn, Py_ie] =  integrateBasis1D(Topography, Presxy(2,:), nnodel, nip);
  [Px_inn, Px_ie] =  integrateMM1D(Topography, Presx);                    % [N, N]
  [Py_inn, Py_ie] =  integrateMM1D(Topography, Presy);                    % [N, N]
  
  elidunsort = find(any(ismember(ELEM2NODE,Gcooids(2:2:end))));            % [1,topnel]
  EL2ENOD    = ELEM2NODE(4:6,:);                                           % edge nodes in all elements
  topednod   = EL2ENOD(ismember(ELEM2NODE(4:6,:),Gcooids(2:2:end)));       % [topnel,1] topographical edge node for each topographical element
  [~, isort] = sort(GCOORD(1,topednod));                                   % [1,topnel]
  elids = elidunsort(isort);                                               % [1,topnel] left-to-right topographical elements; figure(); plot_meshF(ELEM2NODE, GCOORD, elids)
  if any(diff(GCOORD(1,ELEM2NODE(7,elids))) <= 0.) || ...
      ~all(any((ELEM2NODE(4:6,elids) - repmat(Gcooids(2:2:end),3,1)) == 0)) % check elid matches topography nodes 
      error("getHydrostaticPressure:: ---ERR002---")
  end
  [~,J] = sort(isort);                                                     % for inverse sorting:: elidunsort == elids(J)

  if 1 > 2 % tmp - check integrals
      xyv = Topography(:,1:2:end);                                        % corner nodes in topography
      xye = Topography(:,2:2:end);                                        % edge nodes in topography
      dseg = sqrt(sum((xyv(:,2:end) - xyv(:,1:end-1)).^2));               % segment length
      Py_ie_unit = Py_ie ./ dseg;                                         % [Pa]
      pltfac = 1 / max(abs(Py_ie_unit));
      plot_meshF(GCOORD,ELEM2NODE);
      for i=1:topnel
          hold on; plot(repelem(xye(1,i),2)/1000, repelem(xye(2,i),2)/1000+[0 pltfac*Py_ie_unit(i)],'LineWidth',2,'color','blue');
      end     
  end   % tmp
  
  P = [];
  P.inn = GCOORD * 0.0;                                                    % (1,nnod7)
  P.inn(:,Gcooids) = [Px_inn; Py_inn];                                     % [N]
  % xy = Topography; hold on; plot((xy(1,i)+[0. P.inn(1,Gcooids(i))*symsize])/1000, (xy(2,i)+[0. P.inn(2,Gcooids(i))*symsize])/1000,'-','color',[0 0 0.8],'linewidth',2.0);
  P.ie = zeros(2,nel);
  P.ie(:,elids) = [Px_ie; Py_ie];                                          % [N] integral over the element. Only for P.loadel evaluation
  P.loadel = sum(abs(P.ie)) > 0.;                                          % LOGICAL (1,nel)
  
  rhoel = zeros(1,topnel);                                               % [kg/m3] mean density on loaded elements
  relev   = llevel - Topography(2,:);                                      % [m] relative level with respect to the water surface
  relevel = relev(topEL2NOD);                                              % [nnoded=3,topnel]
  rhoel(relevel(1,:) > 0. & relevel(3,:) > 0.) = lrho;                     % fully submerged topographical edges
  psub = sign(relevel(1,:)) ~= sign(relevel(3,:));                         % LOGICAL (1,topnel)
  if any(psub)                                                             % partially submerged
      rhoel(psub) = lrho * max(relevel(:,psub)) ./ abs(diff(relevel([1 3],psub)));
  end
  P.rhoel = rhoel(J);                                                      % map into unsorted topo elements [as used by fssa()]
end % function
