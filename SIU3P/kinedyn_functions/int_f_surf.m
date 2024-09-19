function F= int_f_surf(surf_elem,local_nodes,F_n,dim,GCOORD,ELEM2NODE)
  % F = INT_F_SURF(SURF_ELEM,LOCAL_NODES,F_IP,DIM,GCOORD,ELEM2NODE)
  % 1D surface integration of loads F_IP defined at 1D surface integration 
  % points in the surface of a 2D mesh. SURF_ELEM are the indexes of the
  % elements where the load is apply and LOCAL_NODES are the indexes refering
  % to the element nodes in the surface. DIM defines for which dimension the
  % body force F is calculated. GCOORD and ELEM2NODE are the node coordinates
  % and connectivity mesh matrix respectively.
  
  
  %--------------------------------------------------------------------------
  % Function written by Miguel Andres-Martinez, Postdoc at University of
  % Bremen, 13-10-2017. Email: andresma@uni-bremen.de
  %--------------------------------------------------------------------------
  % Javier García-Pintado: modified to return a [2,nnod7] F matrix 
  %
  %==========================================================================
  % SETUP
  %==========================================================================
  
  nnod_el_s = 3;                                                           % Number of superficial nodes in the element
  
  F = GCOORD * 0.0;                                                          % Vector of the load forces [dof, nnod7]
  
  % % % Plot load (uncomment)
  % plot(Topo_load(1,:),Topo_load(2,:)+Topography(2,:))
  % hold on
  % axis equal
  % axis([10000 30000 -10000 5000])
  % plot(Topography(1,:),Topography(2,:),'k')
  % trimesh(ELEM2NODE(1:3,Top_elem3)',GCOORD(1,:),GCOORD(2,:),'Color',[0,0,0])

  %==========================================================================
  % INTEGRATION POINTS, WEIGHTS AND SHAPE FUNCTIONS
  %==========================================================================
  
  nip = 3;                                                                 % Number of integration points
  Ip = [-sqrt(3/5), 0, sqrt(3/5)];                                         % Integration points local coordinates
  Ipw = [5/9; 8/9; 5/9];                                                   % Integration weights
  N = shp_line_int(Ip,nnod_el_s);                                          % Shape functions for surface integration

  %==========================================================================
  % LOOP THROUGH THE SURFACE ELEMENTS
  %==========================================================================
  for nn = 1:size(surf_elem,2)
    % Global indexes of the surface nodes
    Gl_nodes_s = ELEM2NODE(local_nodes(:,nn),surf_elem(nn));
    
    % Order the surface elements by increasing x order
    [Xorder,Xindex] = sort(GCOORD(1,Gl_nodes_s));
    % Reorder the global node indexes
    Gl_nodes_s = Gl_nodes_s(Xindex);
    % Calculate dx
    dx_e = (max(Xorder)-min(Xorder));
    
    F_index = Gl_nodes_s;
    % %     % Element surface coordinates
    % %     Topo_local = Topography(2,ismember(Topo2nodes,Gl_nodes_s));
    % %     % Find the ip coordinates
    % %     X_ip = Ip*dx_e/2+Xorder(2);
    % %     Y_ip = interp1(Xorder,Topo_local,X_ip,'linear');
    %     
    % %     % Find the load height at the nodes
    % %     Lnode_e = Topo_load(2,ismember(Topo2nodes,Gl_nodes_s));
    % %     % Find the load height at the ips
    % %     Lip_e = interp1(Xorder,Lnode_e,X_ip,'linear');
    % %     % Correct wrong interpolations at the edges of the sea
    % %     if ~isempty(max_level)
    % %         Lip_e((Lip_e+Y_ip)>max_level) = max_level- ...
    % %             Y_ip((Lip_e+Y_ip)>max_level);
    % %         Lip_e(Lip_e<0) = 0;
    % %     end
    %     
    % %     % Calculate the increment in height
    % %     dh_e = (Topo_local(3)-Topo_local(1));
    % %     % Module of the slope vector
    % %     slope_m = sqrt(dx_e^2 + dh_e^2);
    % %     % Perpendicular-to-the-surface unitarian vector pointing downwards
    % %     n = [dh_e,-dx_e]/slope_m;
    % %     % Components of the ip load
    % %     Lip_x = Lip_e'*n(1);
    % %     Lip_y = Lip_e'*n(2);

    % Plot nodes 'o', ips 'x' and load heigth over ips '^' (uncomment)
    
    %     hold on
    %     plot(Xorder,Topo_local,'o')
    %     plot(X_ip,Y_ip,'x')
    %     plot(X_ip,Y_ip+Lip_e,'^')
    %     % Plot directed heigths over the ips
    %     plot([X_ip; X_ip-Lip_x'],[Y_ip; Y_ip-Lip_y'],'r')

    F_d =  zeros(nnod_el_s,1); % Initialize the vector F_d for the nodes in each element.
    % Integration loop
    for ip_s = 1:nip
        F_d = F_d + ...
            N{ip_s}*Ipw(ip_s)*F_n(ip_s,nn).*dx_e/2;
    end
    F(2,F_index) = F(2,F_index) + reshape(F_d,1,3); % save F
  end

  % Plot directed integrated heights over the nodes (uncomment)
  % pnodes = find(F(1:2:end-1)~=0 & F(2:2:end)~=0);
  % Fh = F/(-G(2)*rho_l);
  % plot([GCOORD(1,pnodes); GCOORD(1,pnodes)-Fh(pnodes*2-1)'], ...
  %     [GCOORD(2,pnodes); GCOORD(2,pnodes)-Fh(pnodes*2)'],'*b')
end