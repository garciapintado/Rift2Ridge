function PLITHOST = plithost_vintegr(GCOORD,ELEM2NODE,Point_id,Phases,...
    RHO,G,nip)
% PLITHOS = PLITHOS(GCOORD,ELEM2NODE,POINT_ID,PHASES,RHO,G) calculates the
% lithostatic pressure at the integration points using the mesh defined by
% GCOORD, ELEM2NODE, POINT_ID, PHASES and the density RHO and gravity
% vector G.

%--------------------------------------------------------------------------
% Function written by Miguel Andres-Martinez, Postdoc at University of
% Bremen & MARUM, 27-01-2017. Email: mandresmartinez87@gmail.com
%--------------------------------------------------------------------------

%disp('CALCULATING LITHOSTATIC PRESSURES: ')

t = tic;

% Find topography
[Topography,~] = find_topo(GCOORD,ELEM2NODE,Point_id);

% Resolution of the density mesh
res = min(diff(Topography(1,:)));

%==========================================================================
% START THE IPS AND THE MESHES OVER THE ELEMENTS
%==========================================================================

% Calculate GIPS
[GIP_x,GIP_y] = ip_coord(GCOORD,ELEM2NODE,size(ELEM2NODE,2),nip);

% Vectorize GIPs
Gip_x = GIP_x(:);
Gip_y = GIP_y(:);

% Densities to 3 nodes
RHOn = ip2nodes(RHO,GCOORD,ELEM2NODE,6,3);

% Points in topography directly above the integration points
Topo_ip = interp1(Topography(1,:),Topography(2,:),Gip_x);

% Nodes for integration of the density
max_node_rho = max(GCOORD(2,:));
min_node_rho = min(GCOORD(2,:));
depth_model = max_node_rho-min_node_rho;
Nodes_rho_y = linspace(max_node_rho,min_node_rho,floor(depth_model/res));
NODES_rho = repmat(Nodes_rho_y',1,length(Gip_x));
NODESx = repmat(Gip_x',floor(depth_model/res),1);

%==========================================================================
% DENSITY INTERPOLATION FROM IPS TO COLUMNS OVER IP
%==========================================================================
% Find which elements contain the NODES_rho
Tris_F = tsearch2(GCOORD,uint32(ELEM2NODE(1:3,:)), ...
    [NODESx(:)'; NODES_rho(:)']);

Tris_F = reshape(Tris_F,size(NODESx));
Ind = find(Tris_F==0);

% Check of all nodes were found, and if not it solves the problem finding
% the closest element
if(~isempty(Ind))
    for i=1:length(Ind)
        [~, Tris_F(Ind(i))] = min(sqrt( ...
            (GCOORD(1,ELEM2NODE(7,:)) - NODESx(Ind(i))).^2 + ...
            (GCOORD(2,ELEM2NODE(7,:)) - NODES_rho(Ind(i))).^2));
    end
end

if(any(isnan(Tris_F)))
    error('remeshing failed in move_contours');
end

% Calculates the local coordinates of the new ip
xp = NODESx(:)'; % New ip
yp = NODES_rho(:)'; % New ip

x = reshape(GCOORD(1,ELEM2NODE(1:3,Tris_F)),3,size(xp,2)); % Reshape
y = reshape(GCOORD(2,ELEM2NODE(1:3,Tris_F)),3,size(yp,2)); % Reshape

xi = -(-x(1,:).*yp+x(1,:).*y(3,:)-x(3,:).*y(1,:)+xp.*y(1,:)+x(3,:) ...
    .*yp-xp.*y(3,:))./(-x(1,:).*y(3,:)+x(1,:).*y(2,:)-x(2,:).*y(1,:) ...
    +x(2,:).*y(3,:)+x(3,:).*y(1,:)-x(3,:).*y(2,:)); % Local x cood for ip
yi = (x(1,:).*y(2,:)-x(1,:).*yp-x(2,:).*y(1,:)+x(2,:).*yp+xp ...
    .*y(1,:)-xp.*y(2,:))./(-x(1,:).*y(3,:)+x(1,:).*y(2,:)-x(2,:).*y(1,:)...
    +x(2,:).*y(3,:)+x(3,:).*y(1,:)-x(3,:).*y(2,:)); % Local y coord for ip

% Load shape functions
[NN] = shp_triangle([xi' yi'],3);

% Interpolation
R_N = sum(NN.*RHOn(:,Tris_F));
RHO_N = reshape(R_N,size(Tris_F,1),size(Tris_F,2));

%==========================================================================
% INTEGRATION OF THE DENSITIES
%==========================================================================
% Difference in vertical distances
dNODES_rho = -diff(NODES_rho);

% Average distances between vertical nodes
RHO_NA = (RHO_N(1:1:end-1,:)+RHO_N(2:1:end,:))/2;

% Contribution to the pressure of every interval
PRES_int = RHO_NA.*(-G(2)).*dNODES_rho;

% Index matrix for NODES_rho below topography and above ips
INDXn = NODES_rho<=repmat(Topo_ip',floor(depth_model/res),1) & ...
    NODES_rho>=repmat(Gip_y',floor(depth_model/res),1);
INDXint = (INDXn(1:end-1,:)+INDXn(2:end,:))==2;

% Integrate pressures
P_int = zeros(size(INDXint));
P_int(INDXint) = PRES_int(INDXint);
P = sum(P_int);

% TODO Review contributions of the upper and lower boundaries since the
% densities are obtained from previous intervals.
%--------------------------------------------------------------------------
% Add contribution between the upper node and the topography
UPPER_N = cumsum(cumsum(INDXn))==1;
H_up = Topo_ip - NODES_rho(UPPER_N);
P = P+(H_up'.*(-G(2)).*RHO_NA(UPPER_N(1:end-1,:))');

% Add contribution between the lower node and the integration point
LOWER_N = cumsum(cumsum(INDXn(end:-1:1,:)))==1;
LOWER_N(end:-1:1,:) = LOWER_N;
H_lo = NODES_rho(LOWER_N) - Gip_y;
P = P+(H_lo'.*(-G(2)).*RHO_NA(LOWER_N(2:end,:))');
%--------------------------------------------------------------------------

PLITHOST = reshape(P',size(GIP_x,1),size(GIP_x,2));

% % Plot (uncomment)
% plp = sum(PLITHOST,2)/6;
% patch('faces',ELEM2NODE(1:3,:)','vertices',GCOORD'/1000,'facevertexcdata', ...
%     plp(:)/1e6,'FaceColor','flat')
% shading flat
% colorbar
% title('Pressures [MPa]')
% xlabel('Distance [Km]')
% ylabel('Depth [Km]')