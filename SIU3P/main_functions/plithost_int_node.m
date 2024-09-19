function PLITHOST = plithost_int_node(node, GCOORD, ELEM2NODE, Point_id,...
         Rho, G, igrid, eids)
  % PLITHOS = plithost_int_node(NODE,GCOORD,ELEM2NODE,POINT_ID,Rho,G,igrid) 
  %
  % calculates the lithostatic pressure at a given NODE using the mesh 
  % defined by GCOORD, ELEM2NODE, POINT_ID and the density Rho and gravity
  % vector G, by integrating Rho*G over a vertical regular mesh. Note that
  % NODES can be a vector of nodes and not necesarily a single value
  % variable.
  %
  % node  :: [nnode] node index within GCOORD over which lithostatic pressure is going to be calculated  
  % Rho   :: [nel,nip], with nip = 6
  % igrid :: if TRUE, vertical sampling is restricted up to topography and adapted at each node location
  %
  %--------------------------------------------------------------------------
  % Author Miguel Andres-Martinez, MARUM, 20-10-2017. Email: mandresmartinez87@gmail.com
  %--------------------------------------------------------------------------
  % last version:  JGP, 2019-09-10
  %                modified to cope with irregular vertical sampling; needed for calculations after initialization:
  %                new parameter (igrid=true default)

  if nargin < 7
    igrid = true;
  end

  [Topography,~] = find_topo(GCOORD,ELEM2NODE,Point_id);                     % [x;y] east-west domain topography [JGP warning: does not consider crustal break-up]
  res = min(diff(Topography(1,:)));                                          % minimum topography dx taken as dy resolution below 

  %==========================================================================
  % START THE IPS AND THE MESHES OVER THE ELEMENTS
  %==========================================================================

  Rhon = ipval2nodval(ELEM2NODE(1:3,:), Rho, eids);                            % Rho -> Rhon :: [nel,nip] -> [3,nel], where 3 refers to corner nodes in each element

  % Topography directly above the selected node/s
  node_topo = interp1(Topography(1,:),Topography(2,:),GCOORD(1,node));       % [nx] topography sample at each x-node location

  % bounds for integration of the density vertical profiles
  ymin = min(GCOORD(2,:));
  ymax = max(GCOORD(2,:));

  nx = numel(node);
  ny = floor((ymax-ymin)/res);

  xpmat = repmat(GCOORD(1,node),ny,1);                                       % [ny,nx] goal regular interpolation grid

  if (igrid)
    ypmat = zeros(ny,nx);
    for i=1:nx
      ypmat(:,i) = linspace(node_topo(i),GCOORD(2,node(i)),ny);  
    end
  else  
    yat = linspace(ymax,ymin,ny);                                            % [1,ny] top to bottom  
    ypmat = repmat(yat',1,nx);                                               % [ny,nx] goal regular interpolation grid
  end
    
  n = nx*ny;

  xp = reshape(xpmat,1,n);
  yp = reshape(ypmat,1,n);

  %==========================================================================
  % DENSITY INTERPOLATION FROM IPS TO COLUMNS OVER IP
  %==========================================================================
  % Find which elements contain the vertical profile coordinates
  elids = tsearch2(GCOORD,uint32(ELEM2NODE(1:3,:)), [xpmat(:)'; ypmat(:)']); % uint32, [1,n]
  Ind = find(elids==0);

  % Check of all nodes were found, and if not, it solves the problem finding
  % the closest element
  if(~isempty(Ind))
    for i=1:length(Ind)
        [~, elids(Ind(i))] = min(sqrt( ...
            (GCOORD(1,ELEM2NODE(7,:)) - xp(Ind(i))).^2 + ...
            (GCOORD(2,ELEM2NODE(7,:)) - yp(Ind(i))).^2));
    end
  end

  if(any(isnan(elids)))
    error('remeshing failed in move_contours');
  end

  % Calculates local [canonical] coordinates of the points of vertical profiles

  x = reshape(GCOORD(1,ELEM2NODE(1:3,elids)), 3, n);                         % [3,n] global coordinates of element triangle vertices for elids elements
  y = reshape(GCOORD(2,ELEM2NODE(1:3,elids)), 3, n);                         % [3,n] "       " 
  
  %xieta = zeros(2,n);                                                       % map vertical profiles into the canonical triangle
  %for i=1:n
     %xieta(:,i) = cooGlo2Loc([x(:,i)'; y(:,i)'], [xp(i);yp(i)]);
  %end
                                                                           
  xi =  (  x(1,:) .* yp     - x(1,:) .* y(3,:) + x(3,:).*y(1,:) ... 
         - y(1,:) .* xp     - x(3,:) .* yp     + y(3,:).*xp) ./ ...
        (- x(1,:) .* y(3,:) + x(1,:) .* y(2,:) ...
         - x(2,:) .* y(1,:) + x(2,:) .* y(3,:) ...
         - x(3,:) .* y(2,:) + x(3,:) .* y(1,:));                             % Local x coord for ip
  yi =  (+ y(1,:) .* xp     - y(1,:) .* x(2,:) + y(2,:) .* x(1,:)   ...
         - x(1,:) .* yp     + x(2,:) .* yp     - y(2,:) .* xp) ./ ...
         (-x(1,:) .* y(3,:) + x(1,:).*y(2,:) - x(2,:) .* y(1,:) ...
         + x(2,:) .* y(3,:) + x(3,:).*y(1,:) - x(3,:) .* y(2,:));            % Local y coord for ip
  % NN = shp_triangle(xieta',3);
  NN = shp_triangle([xi' yi'],3);                                            % [3,n] shape functions evaluated at the profile (normalised) coordinates
  Rho_xyp = reshape(sum(NN .* Rhon(:,elids)),ny,nx);                         % interpolate from nodes into vertical profile coordinates

  %==========================================================================
  % INTEGRATION OF THE DENSITIES
  %==========================================================================
  rho_dy  = -diff(ypmat);                                                    % (m)     [ny-1,nx] vertical distance between vertical profile points
  rho_val = (Rho_xyp(1:end-1,:)+Rho_xyp(2:end,:))/2;                         % (kg/m3) [ny-1,nx] average density   between vertical profile points
  pre_val = -G(2) * (rho_val .* rho_dy);                                     %         [ny-1,nx] contribution to pressure

  % Index matrix for profile coordinates below topography and above ips
  if igrid
     P = sum(pre_val);
  else
    INDXn = ypmat <= repmat(node_topo,ny,1) & ...
            ypmat >= repmat(GCOORD(2,node),ny,1);
    INDXint = (INDXn(1:end-1,:)+INDXn(2:end,:)) == 2;                        %         [ny-1,nx] only considers interval in which both bounds are within topography and above ips  
    pre_val(~INDXint) = 0.0;                                                 %         [nx] integrated pressure
    P = sum(pre_val);
  end
  % % TODO Review contributions of the upper and lower boundaries since the
  % % densities are obtained from previous intervals.
  % %--------------------------------------------------------------------------
  % % Add contribution between the upper node and the topography
  % UPPER_N = cumsum(cumsum(INDXn))==1;
  % H_up = Topo_ip - NODES_rho(UPPER_N);
  % P = P+(H_up'.*(-G(2)).*Rho_NA(UPPER_N(1:end-1,:))');
  % 
  % % Add contribution between the lower node and the integration point
  % LOWER_N = cumsum(cumsum(INDXn(end:-1:1,:)))==1;
  % LOWER_N(end:-1:1,:) = LOWER_N;
  % H_lo = NODES_rho(LOWER_N) - GCOORD(2,node);
  % P = P+(H_lo'.*(-G(2)).*Rho_NA(LOWER_N(2:end,:))');
  % %--------------------------------------------------------------------------

  PLITHOST = reshape(P',size(node,1),size(node,2));

  % % Plot (uncomment)
  % plp = sum(PLITHOST,2)/6;
  % patch('faces',ELEM2NODE(1:3,:)','vertices',GCOORD'/1000,'facevertexcdata', ...
  %     plp(:)/1e6,'FaceColor','flat')
  % shading flat
  % colorbar
  % title('Pressures [MPa]')
  % xlabel('Distance [Km]')
  % ylabel('Depth [Km]')
end % function
