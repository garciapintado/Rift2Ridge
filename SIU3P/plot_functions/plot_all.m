% PLOT_ General plots for the models

%MAPPING FROM INTEGRATION POINTS TO NODES - IP CANNOT BE PLOTTED
EL2N      = zeros(nel,3); % new connectivity matrix
GCOORD_N  = zeros(2,nel*3);
TAU_xxn   = zeros(3,nel);
TAU_yyn   = zeros(3,nel);
Pres_n    = zeros(3,nel);
Mu_n      = zeros(3,nel);
E2_n      = zeros(3,nel);
T_n      = zeros(3,nel);
EL2N(1,:) = 1:3;
nip = 6;
nnodel = 6;
[IP_X, IP_w]    = ip_triangle(nip);
[   Nbig]    = shp_triangle(IP_X, nnodel);
for i=1:nel
    is         = (i-1)*3+1; ie = (i-1)*3+3;
    GCOORD_N(:,is:ie) = GCOORD(:,ELEM2NODE([1 2 3],i));
    EL2N(i,:) = is:ie;
    T_n(is:ie) = Temp(ELEM2NODE([1 2 3],i));
    Dummy      = Nbig'\TAU_xx(i,:)';
    TAU_xxn(:,i)= Dummy(1:3);
    Dummy      = Nbig'\TAU_yy(i,:)';
    TAU_yyn(:,i)= Dummy(1:3);
    Dummy      = Nbig'\Mu_all(i,:)';
    Mu_n(:,i)= Dummy(1:3);
    Dummy      = Nbig'\PRES_IP(i,:)';
    Pres_n(:,i)= Dummy(1:3);
    Dummy      = Nbig'\E2all(i,:)';
    E2_n(:,i)= Dummy(1:3);
    
end

% Calculates the density field
DENSITY = Rho(Phases)';
DENSITY = repmat(DENSITY,3,1);

% BUILDS THE BOX AND THE INTERFACES FOR THE PLOT
n_interf = (max(Point_id)-4)/3;
interf_npoint = 0;
for i = 1:n_interf;
    interf_npoint = interf_npoint + sum(Point_id==i*3);
end
Boxx = zeros(1,sum(Point_id~=0,2)-3-interf_npoint);
Boxy = zeros(1,sum(Point_id~=0,2)-3-interf_npoint);
Boxord = zeros(1,sum(Point_id~=0,2)-3-interf_npoint);
ind0 = 2;
ind1 = 1+sum(Point_id==1);

Boxx(1) = GCOORD(1,Corner_id(1));
Boxy(1) = GCOORD(2,Corner_id(1));

[Boxx(ind0:ind1),Boxord(1:sum(Point_id==1))] = sort(GCOORD(1,Point_id==1));
Tmp = GCOORD(2,Point_id==1);
Boxy(ind0:ind1) = Tmp(Boxord(1:sum(Point_id==1)));

ind0 = ind1+1;
ind1 = ind0;
Boxx(ind0) = GCOORD(1,Corner_id(2));
Boxy(ind0) = GCOORD(2,Corner_id(2));

point_id = 2;
max_point_id = max(Point_id);

while point_id < max_point_id
    ind0 = ind1+1;
    ind1 = ind0-1+sum(Point_id==point_id);
    [Boxy(ind0:ind1),Boxord(ind0:ind1)] = sort(GCOORD(2,Point_id==point_id));
    Tmp = GCOORD(1,Point_id==point_id);
    Boxx(ind0:ind1) = Tmp(Boxord(ind0:ind1));
    point_id = point_id+3;
end

ind0 = ind1+1;
ind1 = ind0;
Boxx(ind0) = GCOORD(1,Corner_id(3));
Boxy(ind0) = GCOORD(2,Corner_id(3));

point_id = point_id-2;
ind0 = ind1+1;
ind1 = ind0-1+sum(Point_id==point_id);
[Boxx(ind0:ind1),Boxord(ind0:ind1)] = sort(GCOORD(1,Point_id==point_id),'descend');
Tmp = GCOORD(2,Point_id==point_id);
Boxy(ind0:ind1) = Tmp(Boxord(ind0:ind1));

ind0 = ind1+1;
ind1 = ind0;
Boxx(ind0) = GCOORD(1,Corner_id(4));
Boxy(ind0) = GCOORD(2,Corner_id(4));

point_id = point_id+1;

while point_id > 1
    ind0 = ind1+1;
    ind1 = ind0-1+sum(Point_id==point_id);
    [Boxy(ind0:ind1),Boxord(ind0:ind1)] = sort(GCOORD(2,Point_id==point_id),'descend');
    Tmp = GCOORD(1,Point_id==point_id);
    Boxx(ind0:ind1) = Tmp(Boxord(ind0:ind1));
    point_id = point_id-3;
end

ind0 = ind1+1;
ind1 = ind0;
Boxx(ind0) = GCOORD(1,Corner_id(1));
Boxy(ind0) = GCOORD(2,Corner_id(1));

Boxx = Boxx/km;
Boxy = Boxy/km;

%Chapuza!!!!
while Boxx(end)==0 && Boxy(end)==0
    Boxx(end)=[];
    Boxy(end)=[];
end
%End chapuza

% Interface 3
Point_id(Cornin_id(1:2)) = 3;
Interf_nodes3 = find(Point_id==3);
Interf_elem3 = find(sum(ismember(ELEM2NODE,Interf_nodes3),1)==3);

% Interface 6
Point_id(Cornin_id(3:4)) = 6;
Interf_nodes6 = find(Point_id==6);
Interf_elem6 = find(sum(ismember(ELEM2NODE,Interf_nodes6),1)==3);

% Interface 9
Point_id(Cornin_id(5:6)) = 9;
Interf_nodes9 = find(Point_id==9);
Interf_elem9 = find(sum(ismember(ELEM2NODE,Interf_nodes9),1)==3);

clf

subplot(231),patch('faces',EL2N,'vertices',GCOORD_N'/1000,'facevertexcdata', ...
    TAU_xxn(:),'FaceColor','flat')
shading interp
axis tight
colormap(jet)
freezeColors
colorbar
cbfreeze(colorbar)
hold on
plot(Boxx,Boxy,'k')
hold on

% Plots interface 3
for i = 1:size(Interf_elem3,2)
    Interf_node_el = Point_id(ELEM2NODE(1:6,Interf_elem3(i)))==3;
    Interf_gl = ELEM2NODE(Interf_node_el,Interf_elem3(i));
    [X_int,X_int_ord] = sort(GCOORD(1,Interf_gl));
    Interf_gl = Interf_gl(X_int_ord);
    Y_int = GCOORD(2,Interf_gl);
    plot(X_int/km,Y_int/km,'k')
    hold on
end

% Plots interface 6
for i = 1:size(Interf_elem6,2)
    Interf_node_el = Point_id(ELEM2NODE(1:6,Interf_elem6(i)))==6;
    Interf_gl = ELEM2NODE(Interf_node_el,Interf_elem6(i));
    [X_int,X_int_ord] = sort(GCOORD(1,Interf_gl));
    Interf_gl = Interf_gl(X_int_ord);
    Y_int = GCOORD(2,Interf_gl);
    plot(X_int/km,Y_int/km,'k')
    hold on
end

% Plots interface 9
for i = 1:size(Interf_elem9,2)
    Interf_node_el = Point_id(ELEM2NODE(1:6,Interf_elem9(i)))==9;
    Interf_gl = ELEM2NODE(Interf_node_el,Interf_elem9(i));
    [X_int,X_int_ord] = sort(GCOORD(1,Interf_gl));
    Interf_gl = Interf_gl(X_int_ord);
    Y_int = GCOORD(2,Interf_gl);
    plot(X_int/km,Y_int/km,'k')
    hold on
end

xlabel('Distance [km]')
ylabel('Depth [km]')

title('Horizontal Deviatoric Stress [Pa]')

drawnow

subplot(232),patch('faces',EL2N,'vertices',GCOORD_N'/1000,'facevertexcdata',TAU_yyn(:),'FaceColor','flat')
shading interp
axis tight
colormap(jet)
freezeColors
colorbar
cbfreeze(colorbar)
hold on
plot(Boxx,Boxy,'k')
hold on

% Plots interface 3
for i = 1:size(Interf_elem3,2)
    Interf_node_el = Point_id(ELEM2NODE(1:6,Interf_elem3(i)))==3;
    Interf_gl = ELEM2NODE(Interf_node_el,Interf_elem3(i));
    [X_int,X_int_ord] = sort(GCOORD(1,Interf_gl));
    Interf_gl = Interf_gl(X_int_ord);
    Y_int = GCOORD(2,Interf_gl);
    plot(X_int/km,Y_int/km,'k')
    hold on
end

% Plots interface 6
for i = 1:size(Interf_elem6,2)
    Interf_node_el = Point_id(ELEM2NODE(1:6,Interf_elem6(i)))==6;
    Interf_gl = ELEM2NODE(Interf_node_el,Interf_elem6(i));
    [X_int,X_int_ord] = sort(GCOORD(1,Interf_gl));
    Interf_gl = Interf_gl(X_int_ord);
    Y_int = GCOORD(2,Interf_gl);
    plot(X_int/km,Y_int/km,'k')
    hold on
end

% Plots interface 9
for i = 1:size(Interf_elem9,2)
    Interf_node_el = Point_id(ELEM2NODE(1:6,Interf_elem9(i)))==9;
    Interf_gl = ELEM2NODE(Interf_node_el,Interf_elem9(i));
    [X_int,X_int_ord] = sort(GCOORD(1,Interf_gl));
    Interf_gl = Interf_gl(X_int_ord);
    Y_int = GCOORD(2,Interf_gl);
    plot(X_int/km,Y_int/km,'k')
    hold on
end

xlabel('Distance [km]')
ylabel('Depth [km]')

title('Vertical Deviatoric Stress [Pa]')

drawnow

subplot(233)
quiver(GCOORD(1,:)/km, GCOORD(2,:)/km, Vel(1:2:end-1)', Vel(2:2:end)','k');
title('Velocity and Density Fields')
xlabel('Distance [km]')
ylabel('Depth [km]')
hold on
patch('faces',EL2N,'vertices',GCOORD_N'/1000,'facevertexcdata',DENSITY(:),'FaceColor','flat')
axis tight
shading interp
colormap(jet)
freezeColors
colorbar
cbfreeze(colorbar)
hold on
plot(Boxx,Boxy,'k')
hold on

% Plots interface 3
for i = 1:size(Interf_elem3,2)
    Interf_node_el = Point_id(ELEM2NODE(1:6,Interf_elem3(i)))==3;
    Interf_gl = ELEM2NODE(Interf_node_el,Interf_elem3(i));
    [X_int,X_int_ord] = sort(GCOORD(1,Interf_gl));
    Interf_gl = Interf_gl(X_int_ord);
    Y_int = GCOORD(2,Interf_gl);
    plot(X_int/km,Y_int/km,'k')
    hold on
end

% Plots interface 6
for i = 1:size(Interf_elem6,2)
    Interf_node_el = Point_id(ELEM2NODE(1:6,Interf_elem6(i)))==6;
    Interf_gl = ELEM2NODE(Interf_node_el,Interf_elem6(i));
    [X_int,X_int_ord] = sort(GCOORD(1,Interf_gl));
    Interf_gl = Interf_gl(X_int_ord);
    Y_int = GCOORD(2,Interf_gl);
    plot(X_int/km,Y_int/km,'k')
    hold on
end

% Plots interface 9
for i = 1:size(Interf_elem9,2)
    Interf_node_el = Point_id(ELEM2NODE(1:6,Interf_elem9(i)))==9;
    Interf_gl = ELEM2NODE(Interf_node_el,Interf_elem9(i));
    [X_int,X_int_ord] = sort(GCOORD(1,Interf_gl));
    Interf_gl = Interf_gl(X_int_ord);
    Y_int = GCOORD(2,Interf_gl);
    plot(X_int/km,Y_int/km,'k')
    hold on
end

drawnow

Surfx = [GCOORD(1,Corner_id(3)) GCOORD(1,Corner_id(4)) GCOORD(1,Point_id==12)];
Surfy = [GCOORD(2,Corner_id(3)) GCOORD(2,Corner_id(4)) GCOORD(2,Point_id==12)];
[Surfx, Surford] = sort(Surfx);
Surfy = Surfy(Surford);
Velx_surf = Vel(1:2:end-1);
Velx_surf = Velx_surf(Point_id==max(Point_id)-1);
Vely_surf = Vel(2:2:end);
Vely_surf = Vely_surf(Point_id==max(Point_id)-1);

subplot(234)
plot(Surfx/km,Surfy,'k');
hold on
quiver(GCOORD(1,Point_id==max(Point_id)-1)/km, GCOORD(2,Point_id==max(Point_id)-1), Velx_surf', Vely_surf',0.2,'k');
axis tight
title('Surface Velocities')
ylabel('Topography [m]')
xlabel('Distance [km]')

drawnow

% Viscosity plot
subplot(235)
title('Viscosity Field')
xlabel('Distance [km]')
ylabel('Depth [km]')
h = patch('faces',EL2N,'vertices',GCOORD_N'/1000,'facevertexcdata',log10(Mu_n(:)),'FaceColor','flat');
axis tight
shading interp
colormap(flipud(jet))
freezeColors
colorbar
cbfreeze(colorbar)
hold on
plot(Boxx,Boxy,'k')
hold on

% Plots interface 3
for i = 1:size(Interf_elem3,2)
    Interf_node_el = Point_id(ELEM2NODE(1:6,Interf_elem3(i)))==3;
    Interf_gl = ELEM2NODE(Interf_node_el,Interf_elem3(i));
    [X_int,X_int_ord] = sort(GCOORD(1,Interf_gl));
    Interf_gl = Interf_gl(X_int_ord);
    Y_int = GCOORD(2,Interf_gl);
    plot(X_int/km,Y_int/km,'k')
    hold on
end

% Plots interface 6
for i = 1:size(Interf_elem6,2)
    Interf_node_el = Point_id(ELEM2NODE(1:6,Interf_elem6(i)))==6;
    Interf_gl = ELEM2NODE(Interf_node_el,Interf_elem6(i));
    [X_int,X_int_ord] = sort(GCOORD(1,Interf_gl));
    Interf_gl = Interf_gl(X_int_ord);
    Y_int = GCOORD(2,Interf_gl);
    plot(X_int/km,Y_int/km,'k')
    hold on
end

% Plots interface 9
for i = 1:size(Interf_elem9,2)
    Interf_node_el = Point_id(ELEM2NODE(1:6,Interf_elem9(i)))==9;
    Interf_gl = ELEM2NODE(Interf_node_el,Interf_elem9(i));
    [X_int,X_int_ord] = sort(GCOORD(1,Interf_gl));
    Interf_gl = Interf_gl(X_int_ord);
    Y_int = GCOORD(2,Interf_gl);
    plot(X_int/km,Y_int/km,'k')
    hold on
end

drawnow

% 2nd strain rate invariant
subplot(236)
title('2nd strain rate invariant')
xlabel('Distance [km]')
ylabel('Depth [km]')
patch('faces',EL2N,'vertices',GCOORD_N'/1000,'facevertexcdata',E2_n(:),'FaceColor','flat')
axis tight
shading interp
colormap(jet)
freezeColors
colorbar
cbfreeze(colorbar)
hold on
plot(Boxx,Boxy,'k')
hold on

% Plots interface 3
for i = 1:size(Interf_elem3,2)
    Interf_node_el = Point_id(ELEM2NODE(1:6,Interf_elem3(i)))==3;
    Interf_gl = ELEM2NODE(Interf_node_el,Interf_elem3(i));
    [X_int,X_int_ord] = sort(GCOORD(1,Interf_gl));
    Interf_gl = Interf_gl(X_int_ord);
    Y_int = GCOORD(2,Interf_gl);
    plot(X_int/km,Y_int/km,'k')
    hold on
end

% Plots interface 6
for i = 1:size(Interf_elem6,2)
    Interf_node_el = Point_id(ELEM2NODE(1:6,Interf_elem6(i)))==6;
    Interf_gl = ELEM2NODE(Interf_node_el,Interf_elem6(i));
    [X_int,X_int_ord] = sort(GCOORD(1,Interf_gl));
    Interf_gl = Interf_gl(X_int_ord);
    Y_int = GCOORD(2,Interf_gl);
    plot(X_int/km,Y_int/km,'k')
    hold on
end

% Plots interface 9
for i = 1:size(Interf_elem9,2)
    Interf_node_el = Point_id(ELEM2NODE(1:6,Interf_elem9(i)))==9;
    Interf_gl = ELEM2NODE(Interf_node_el,Interf_elem9(i));
    [X_int,X_int_ord] = sort(GCOORD(1,Interf_gl));
    Interf_gl = Interf_gl(X_int_ord);
    Y_int = GCOORD(2,Interf_gl);
    plot(X_int/km,Y_int/km,'k')
    hold on
end

drawnow

%     figure(1), clf
%     subplot(121),patch('faces',EL2N,'vertices',GCOORD_N'/1000,'facevertexcdata',T_n(:),'FaceColor','flat')
%     shading interp
%             axis tight
%     title('Temperature [�C]')
%     xlabel('Distance [km]')
%         ylabel('Depth [km]')
%     colorbar
%     
%     
%    subplot(122),patch('faces',EL2N,'vertices',GCOORD_N'/1000,'facevertexcdata',TAU_xxn(:),'FaceColor','flat')
%     shading interp
%         axis tight
%     colorbar
%         xlabel('Distance [km]')
%         ylabel('Depth [km]')
% 
%         title('Horizontal Deviatoric Stress [Pa]')
%         
%         drawnow
%         
%      subplot(223),patch('faces',EL2N,'vertices',GCOORD_N','facevertexcdata',E2_n(:),'FaceColor','flat')
%     shading interp
%         axis tight
%     colorbar
%        title('E2')
       
%      subplot(224),patch('faces',EL2N,'vertices',GCOORD_N','facevertexcdata',TAU_xxn(:),'FaceColor','flat')
%     shading interp
%     colorbar
%     hold on
%     quiver(GCOORD(1,:), GCOORD(2,:), Vel(1:2:end-1)', Vel(2:2:end)','w');
%     axis tight
%         title('Tau.')
%     drawnow
% 
%       figure(2), clf
%     subplot(211),patch('faces',EL2N,'vertices',GCOORD_N','facevertexcdata',Pres_n(:),'FaceColor','flat')
% %     shading interp
%         axis tight
%     colorbar
    
    
%    subplot(212)
%     plot(GCOORD(1,Point_id==9), GCOORD(2,Point_id==9),'bx', GCOORD(1,Point_id==6), GCOORD(2,Point_id==6),'rx') 
%         axis tight
%     drawnow
%     
%     figure(111),plot(GEOMETRY(1,:), GEOMETRY(2,:), 'rx')
%     drawnow
%     
