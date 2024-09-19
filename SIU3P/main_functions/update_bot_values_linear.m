function [GCOORD,Temp,F_xx,F_xy,F_yx,F_yy,I2,TAU_xx_old,TAU_xy_old, ...
    TAU_yy_old,Mu_all,E2all,remesh] = update_bot_values_linear(GCO,ELEM2NODE, ...
    ymin,Point_id,Temp,temp_bc,F_xx,F_xy,F_yx,F_yy,I2, ...
    TAU_xx_old,TAU_xy_old,TAU_yy_old,Mu_all,E2all,R,Phases,remesh)
% UPDATE_TOPO_VALUES takes a mesh which topography has been changed to 
% NEW_TOPO (by i.e. erosion/sedimentation) and updates the coordinates of
% the nodes, calculates the new values of the FIELDS where nodes are still 
% inside of the old mesh (i.e. basement) and assign new values to nodes 
% that fall out of the old mesh (i.e. sediments). It also schedules a
% remeshing for the next time step in case some nodes of the new mesh are
% above the topography (i.e. inversion of an element).
%
% FIELDS:
% -------
%     Temperature              TEMP
%     Gradient of deformation  F*
%     Historic strain          I2
%     Rotated old stresses     TAU_*_OLD
%     Viscosity                Mu_all
%     Strain rate              E2ALL

%--------------------------------------------------------------------------
% Function written by Miguel Andres-Martinez, PhD student at Royal Holloway
% University of London, 05-09-2015. Email: mandresmartinez87@gmail.com
%--------------------------------------------------------------------------

% TODO add remesh of Rheol_var and ramdon element

%==========================================================================
% INITIALIZE VARIABLES
%==========================================================================
tic; 
GCOORD = GCO;
Bot_nodes = find(Point_id==1);
GCOORD(2,Bot_nodes) = ymin;

%==========================================================================
% RESHAPED ELEMENTS
%==========================================================================
% Nodes that changed y-coordinate
change_bot = GCO(2,Bot_nodes) ~= GCOORD(2,Bot_nodes);
NODES_CT = Bot_nodes(change_bot);
% Find reshaped elements
reshape_el = find(sum(ismember(ELEM2NODE,NODES_CT))>0);

%==========================================================================
% CALCULATE NEW EDGE AND CENTRAL NODES
%==========================================================================
% Make straight edges
GCOORD(:,ELEM2NODE([6 4 5],reshape_el)) = ...
    0.5*(GCOORD(:,ELEM2NODE([1 2 3],reshape_el)) + ...
    GCOORD(:,ELEM2NODE([2 3 1],reshape_el)));
GCOORD(:,ELEM2NODE(7,reshape_el)) = 1/3 * ...
    (GCOORD(:,ELEM2NODE(1,reshape_el)) + ...
    GCOORD(:,ELEM2NODE(2,reshape_el)) + GCOORD(:,ELEM2NODE(3,reshape_el)));
% % Plot (uncomment)
% axis(axisp)
% trimesh(ELEM2NODE(1:3,reshape_el)',GCOORD(1,:)/1000,GCOORD(2,:)/1000, ...
%     'Color',[1 0 0])
% hold on
% trimesh(ELEM2NODE(1:3,reshape_el)',GCO(1,:)/1000,GCO(2,:)/1000, ...
%     'Color',[0 0 0])
% plot(GCOORD(1,ELEM2NODE(:,reshape_el))/1000, ...
%     GCOORD(2,ELEM2NODE(:,reshape_el))/1000,'sr','MarkerFaceColor','r')
% plot(GCO(1,ELEM2NODE(:,reshape_el))/1000, ...
%     GCO(2,ELEM2NODE(:,reshape_el))/1000,'sk','MarkerFaceColor','k')

%==========================================================================
% CHECK FOR NODES BELLOW BOTTOM BOUNDARY
%==========================================================================
Not_bot = ~ismember(1:size(GCOORD,2),Bot_nodes);
[BOT_TOPO,indx_bot] = sort(GCOORD(1,Bot_nodes));
BOT_TOPO(2,:) = GCOORD(2,Bot_nodes(indx_bot));
NEW_TOPO2GCO_Y = interp1(BOT_TOPO(1,:),BOT_TOPO(2,:),GCOORD(1,Not_bot));
if sum(GCOORD(2,Not_bot)<NEW_TOPO2GCO_Y)>0
    remesh = 1;
end

%==========================================================================
% CALCULATE NEW IPS
%==========================================================================
nip = 6;
GIP_x = zeros(size(ELEM2NODE(1:6,reshape_el)'));
GIP_y = zeros(size(ELEM2NODE(1:6,reshape_el)'));
[IP_X,~] = ip_triangle(6);
[N,~] = shp_deriv_triangle(IP_X,7);
for ip=1:nip
    Ni = N{ip};
    ECOORD_x = reshape(GCOORD(1,ELEM2NODE(:,reshape_el)),7, ...
        size(ELEM2NODE(:,reshape_el),2));
    ECOORD_y = reshape(GCOORD(2,ELEM2NODE(:,reshape_el)),7, ...
        size(ELEM2NODE(:,reshape_el),2));
    GIP_x(:,ip) = Ni'*ECOORD_x;
    GIP_y(:,ip) = Ni'*ECOORD_y;
end
% % Plot (uncomment)
% axisp = [41 56 -7.5 -2.5]; % JGP: only needed for plotting below
% plot(GIP_x/1000,GIP_y/1000,'.r')

%==========================================================================
% INDEXES OF THE OLD ELEMENTS CONTAINING THE NEW NODES AND IPS
%==========================================================================
% Nodes
nindx = unique(ELEM2NODE(1:7,reshape_el));
Tris_n = tsearch2(GCO,uint32(ELEM2NODE(1:3,reshape_el)), ...
    [GCOORD(1,nindx);GCOORD(2,nindx)]);
% Ips
Tris_ip = tsearch2(GCO,uint32(ELEM2NODE(1:3,reshape_el)), ...
    [GIP_x(:)';GIP_y(:)']);

% % Plot (uncomment)
% plot(GCOORD(1,nindx(Tris_n~=0))/1000,GCOORD(2,nindx(Tris_n~=0))/1000,'or')
% plot(GIP_x(Tris_ip~=0)/1000,GIP_y(Tris_ip~=0)/1000,'or')

%==========================================================================
% INTERPOLATE FIELDS WHERE NEW NODES AND IPS ARE INSIDE THE OLD MESH
%==========================================================================
% Temperature
Temp(nindx(Tris_n~=0)) = remesh_val(Tris_n(Tris_n~=0),GCO, ...
    GCOORD(:,nindx(Tris_n~=0)),Temp,ELEM2NODE(:,reshape_el));
nnodel_r = 3;

% Gradient of deformation and rotated stresses
Fxx = F_xx(reshape_el,:);
Fxy = F_xy(reshape_el,:);
Fyx = F_yx(reshape_el,:);
Fyy = F_yy(reshape_el,:);
I.p = I2.p(reshape_el,:);
I.c = I2.c(reshape_el,:);
Txx = TAU_xx_old(reshape_el,:);
Txy = TAU_xy_old(reshape_el,:);
Tyy = TAU_yy_old(reshape_el,:);
Mu_old = Mu_all(reshape_el,:);
E_old = E2all(reshape_el,:);

% Interpolate
extrap_scheme = "ip2nodb2ip_N3";
[Fxx,Fxy,Fyx,Fyy,~,~,I,Txx,Txy,Tyy,Mu,E] = remesh_F_TAU_Mu_E2( ...
    GCOORD, ELEM2NODE(:,reshape_el), GCO, ELEM2NODE(:,reshape_el), ...
    Fxx, Fxy, Fyx, Fyy, Txx, Txy, Tyy, Mu_old, E_old, I,...
    extrap_scheme);

%==========================================================================
% ASSIGN VALUES FOR NEW MATERIALS (Tris==0)
%==========================================================================
% Temperature
Temp(nindx(Tris_n==0)) = temp_bc;

% Gradient of deformation (no deformation for new ips)
Fxx(Tris_ip==0) = 1;
F_xx(reshape_el,:) = Fxx;
Fxy(Tris_ip==0) = 0;
F_xy(reshape_el,:) = Fxy;
Fyx(Tris_ip==0) = 0;
F_yx(reshape_el,:) = Fyx;
Fyy(Tris_ip==0) = 1;
F_yy(reshape_el,:) = Fyy;

% Historic second invariant of the strain (no deformation for new ips)
I.f(Tris_ip==0) = 0;
I.p(Tris_ip==0) = 0;
I.c(Tris_ip==0) = 0;
I2.f(reshape_el,:) = I.f;
I2.p(reshape_el,:) = I.p;
I2.c(reshape_el,:) = I.c;

% Rotated stresses
Txx(Tris_ip==0) = 0;
TAU_xx_old(reshape_el,:) = Txx;
Txy(Tris_ip==0) = 0;
TAU_xy_old(reshape_el,:) = Txy;
Tyy(Tris_ip==0) = 0;
TAU_yy_old(reshape_el,:) = Tyy;

% % Viscosity (dislocation)
% Adis_block  = RHEOL.Adis(Phases(reshape_el),:);
% Ndis_block  = RHEOL.Ndis(Phases(reshape_el),:);
% Qdis_block  = RHEOL.Qdis(Phases(reshape_el),:);
% Var_block   = RHEOLvar(Phases(reshape_el),:);
% 
% E2 = mean(E2all(reshape_el,:),2);
% Temp_ip = temp_surf;
% Sc_dis = zeros(size(E2));
% ED = zeros(size(E2,1),6);
% for ip = 1:nip
%     for n = 1:size(Ndis_block,2)
%         Sc_dis = Sc_dis + Var_block(:,n).* ...
%             1./(2.^((Ndis_block(:,n)-1)./Ndis_block(:,n)).* ...
%             3.^((Ndis_block(:,n)+1)./(2*Ndis_block(:,n))));
%         
%         ED(:,ip) = ED(:,ip) + Var_block(:,n).*(Sc_dis.*Adis_block(:,n).^ ...
%             (-1./Ndis_block(:,n)).*E2.^(1./Ndis_block(:,n)-1).* ...
%             exp(Qdis_block(:,n)./(Ndis_block(:,n).*R.*(Temp_ip+273))));
%     end
% end
% ED(ED<=1e18)   = 1e18;               
% ED(ED>1e24)    = 1e24;
% Mu(Tris_ip==0) = ED(Tris_ip==0);
% Mu_all(reshape_el,:) = Mu;
% 
% % Strain rate invariant (average of the element)
% E22 = repmat(E2,1,6);
% E(Tris_ip==0) = E22(Tris_ip==0);
% E2all(reshape_el,:) = E;

% Find which elements the outside ips belong to
% Indexes of the ips
Ip_i = Tris_ip==0;
% Make matrix of element indexes
IP2EL = repmat((1:length(reshape_el))',1,nip);

% Viscosity (average of the rest of the element)
Mu(Tris_ip==0) = sum(Mu_old(IP2EL(Ip_i),:),2)/nip;
Mu(Mu<=1e18)   = 1e18;               
Mu(Mu>1e24)    = 1e24;
Mu_all(reshape_el,:) = Mu;

% Strain rate invariant (average of the rest of the element)
E(Tris_ip==0) = sum(E_old(IP2EL(Ip_i),:),2)/nip;
E2all(reshape_el,:) = E;

% %==========================================================================
% % PLOT FINAL VALUES (UNCOMMENT)
% %==========================================================================
% np = GCOORD(1,:)/1000>=axisp(1) & GCOORD(1,:)/1000<=axisp(2) & ...
%     GCOORD(2,:)/1000>=axisp(3) & GCOORD(2,:)/1000<=axisp(4);
% ipp = GIP_x/1000>=axisp(1) & GIP_x/1000<=axisp(2) & ...
%     GIP_y/1000>=axisp(3) & GIP_y/1000<=axisp(4);
% % Temperature
% scatter(GCOORD(1,np)/1000,GCOORD(2,np)/1000,50,Temp(np),'filled')
% colormap('jet')
% colorbar
% 
% % Fxx
% scatter(GIP_x(ipp)/1000,GIP_y(ipp)/1000,50,Fxx(ipp),'filled')
% colorbar
% % Fxy
% scatter(GIP_x(ipp)/1000,GIP_y(ipp)/1000,50,Fxy(ipp),'filled')
% colorbar
% % Fyx
% scatter(GIP_x(ipp)/1000,GIP_y(ipp)/1000,50,Fyx(ipp),'filled')
% colorbar
% % Fyy
% scatter(GIP_x(ipp)/1000,GIP_y(ipp)/1000,50,Fyy(ipp),'filled')
% colorbar
% 
% % I
% scatter(GIP_x(ipp)/1000,GIP_y(ipp)/1000,50,I(ipp),'filled')
% colorbar
% 
% % Txx
% scatter(GIP_x(ipp)/1000,GIP_y(ipp)/1000,50,Txx(ipp),'filled')
% colorbar
% % Txy
% scatter(GIP_x(ipp)/1000,GIP_y(ipp)/1000,50,Txy(ipp),'filled')
% colorbar
% % Tyy
% scatter(GIP_x(ipp)/1000,GIP_y(ipp)/1000,50,Tyy(ipp),'filled')
% colorbar
% 
% % Mu
% scatter(GIP_x(ipp)/1000,GIP_y(ipp)/1000,50,log10(Mu(ipp)),'filled')
% colorbar
% 
% % E
% scatter(GIP_x(ipp)/1000,GIP_y(ipp)/1000,50,log10(E(ipp)),'filled')
% colorbar

fprintf(1,'BOT REMESH:   ');
fprintf(1,[num2str(toc),'\n']);
