%MAPPING FROM INTEGRATION POINTS TO NODES - IP CANNOT BE PLOTTED
EL2N      = zeros(nel,3); % new connectivity matrix
GCOORD_N  = zeros(2,nel*3);
TAU_xxn   = zeros(3,nel);
TAU_yyn   = zeros(3,nel);
Pres_n    = zeros(3,nel);
Mu_n      = zeros(3,nel);
E2_n      = zeros(3,nel);
Dpl_n      = zeros(3,nel);
dP_melting_n = zeros(3,nel);
T_n      = zeros(3,nel);
dF_n      = zeros(3,nel);
dFup_n      = zeros(3,nel);
dFdiff_n  = zeros(3,nel);
EL2N(1,:) = 1:3;
nip = 6;
nnodel = 6;
[IP_X, IP_w]    = ip_triangle(nip);
[   Nbig]    = shp_triangle(IP_X, nnodel);
for i=1:nel
    is         = (i-1)*3+1; ie = (i-1)*3+3; %is=1,4,7,.. %ie=3,6,9,..
    GCOORD_N(:,is:ie) = GCOORD(:,ELEM2NODE([1 2 3],i)); %GCOORD(ndim,nnod)=> %GCOORD(x e y, nodos 1 2 y 3, de todos los elementos i)
    EL2N(i,:) = is:ie;
    Dpl_n(is:ie) = Dpl(ELEM2NODE([1 2 3],i));  
    %dF_n(is:ie) = dF(ELEM2NODE([1 2 3],i));  
    %dFup_n(is:ie) = dF_up(ELEM2NODE([1 2 3],i));  
    %dFdiff_n(is:ie) = dF_diffusive(ELEM2NODE([1 2 3],i));  
%     dP_melting_n(is:ie) = dP_melting(ELEM2NODE([1 2 3],i));  
    %T_n(is:ie) = Temp(ELEM2NODE([1 2 3],i));  

end

% % Calculates the density field
% DENSITY = Rho(Phases)';
% DENSITY = repmat(DENSITY,3,1);
% 
% % BUILDS THE BOX AND THE INTERFACES FOR THE PLOT
% n_interf = (max(Point_id)-4)/3;
% interf_npoint = 0;
% for i = 1:n_interf;
%     interf_npoint = interf_npoint + sum(Point_id==i*3);
% end
% Boxx = zeros(1,sum(Point_id~=0,2)-3-interf_npoint);
% Boxy = zeros(1,sum(Point_id~=0,2)-3-interf_npoint);
% Boxord = zeros(1,sum(Point_id~=0,2)-3-interf_npoint);
% ind0 = 2;
% ind1 = 1+sum(Point_id==1);
% 
% Boxx(1) = GCOORD(1,Corner_id(1));
% Boxy(1) = GCOORD(2,Corner_id(1));
% 
% [Boxx(ind0:ind1),Boxord(1:sum(Point_id==1))] = sort(GCOORD(1,Point_id==1));
% Tmp = GCOORD(2,Point_id==1);
% Boxy(ind0:ind1) = Tmp(Boxord(1:sum(Point_id==1)));
% 
% ind0 = ind1+1;
% ind1 = ind0;
% Boxx(ind0) = GCOORD(1,Corner_id(2));
% Boxy(ind0) = GCOORD(2,Corner_id(2));
% 
% point_id = 2;
% max_point_id = max(Point_id);
% 
% while point_id < max_point_id
%     ind0 = ind1+1;
%     ind1 = ind0-1+sum(Point_id==point_id);
%     [Boxy(ind0:ind1),Boxord(ind0:ind1)] = sort(GCOORD(2,Point_id==point_id));
%     Tmp = GCOORD(1,Point_id==point_id);
%     Boxx(ind0:ind1) = Tmp(Boxord(ind0:ind1));
%     point_id = point_id+3;
% end
% 
% ind0 = ind1+1;
% ind1 = ind0;
% Boxx(ind0) = GCOORD(1,Corner_id(3));
% Boxy(ind0) = GCOORD(2,Corner_id(3));
% 
% point_id = point_id-2;
% ind0 = ind1+1;
% ind1 = ind0-1+sum(Point_id==point_id);
% [Boxx(ind0:ind1),Boxord(ind0:ind1)] = sort(GCOORD(1,Point_id==point_id),'descend');
% Tmp = GCOORD(2,Point_id==point_id);
% Boxy(ind0:ind1) = Tmp(Boxord(ind0:ind1));
% 
% ind0 = ind1+1;
% ind1 = ind0;
% Boxx(ind0) = GCOORD(1,Corner_id(4));
% Boxy(ind0) = GCOORD(2,Corner_id(4));
% 
% point_id = point_id+1;
% 
% while point_id > 1
%     ind0 = ind1+1;
%     ind1 = ind0-1+sum(Point_id==point_id);
%     [Boxy(ind0:ind1),Boxord(ind0:ind1)] = sort(GCOORD(2,Point_id==point_id),'descend');
%     Tmp = GCOORD(1,Point_id==point_id);
%     Boxx(ind0:ind1) = Tmp(Boxord(ind0:ind1));
%     point_id = point_id-3;
% end
% 
% ind0 = ind1+1;
% ind1 = ind0;
% Boxx(ind0) = GCOORD(1,Corner_id(1));
% Boxy(ind0) = GCOORD(2,Corner_id(1));
% 
% Boxx = Boxx/km;
% Boxy = Boxy/km;
% 
% %Chapuza!!!!
% while Boxx(end)==0 && Boxy(end)==0
%     Boxx(end)=[];
%     Boxy(end)=[];
% end
% %End chapuza
% 
% % Interface 3
% Point_id(Cornin_id(1:2)) = 3;
% Interf_nodes3 = find(Point_id==3);
% Interf_elem3 = find(sum(ismember(ELEM2NODE,Interf_nodes3),1)==3);
% 
% % Interface 6
% Point_id(Cornin_id(3:4)) = 6;
% Interf_nodes6 = find(Point_id==6);
% Interf_elem6 = find(sum(ismember(ELEM2NODE,Interf_nodes6),1)==3);
% 
% % Interface 9
% Point_id(Cornin_id(5:6)) = 9;
% Interf_nodes9 = find(Point_id==9);
% Interf_elem9 = find(sum(ismember(ELEM2NODE,Interf_nodes9),1)==3);

%==========================1 DEPLETION==============================================% 
figure(12)
rect = [-1200, -1200, 1000, 1700];%rect = [left, bottom, width, height]
set(gcf,'OuterPosition',rect)
patch('faces',EL2N,'vertices',GCOORD_N'/1000,'facevertexcdata',Dpl_n(:),'FaceColor','flat')
shading interp
axis tight
colorbar
hold on
plot(Boxx,Boxy,'k')
hold on
caxis([0 0.5])  %define el min y max de los colores de la paleta que pinta

% % Plots interface 3
% for i = 1:size(Interf_elem3,2)
%     Interf_node_el = Point_id(ELEM2NODE(1:6,Interf_elem3(i)))==3;
%     Interf_gl = ELEM2NODE(Interf_node_el,Interf_elem3(i));
%     [X_int,X_int_ord] = sort(GCOORD(1,Interf_gl));
%     Interf_gl = Interf_gl(X_int_ord);
%     Y_int = GCOORD(2,Interf_gl);
%     plot(X_int/km,Y_int/km,'k')
%     hold on
% end
% 
% % Plots interface 6
% for i = 1:size(Interf_elem6,2)
%     Interf_node_el = Point_id(ELEM2NODE(1:6,Interf_elem6(i)))==6;
%     Interf_gl = ELEM2NODE(Interf_node_el,Interf_elem6(i));
%     [X_int,X_int_ord] = sort(GCOORD(1,Interf_gl));
%     Interf_gl = Interf_gl(X_int_ord);
%     Y_int = GCOORD(2,Interf_gl);
%     plot(X_int/km,Y_int/km,'k')
%     hold on
% end
% 
% % Plots interface 9
% for i = 1:size(Interf_elem9,2)
%     Interf_node_el = Point_id(ELEM2NODE(1:6,Interf_elem9(i)))==9;
%     Interf_gl = ELEM2NODE(Interf_node_el,Interf_elem9(i));
%     [X_int,X_int_ord] = sort(GCOORD(1,Interf_gl));
%     Interf_gl = Interf_gl(X_int_ord);
%     Y_int = GCOORD(2,Interf_gl);
%     plot(X_int/km,Y_int/km,'k')
%     hold on
% end

xlabel('Distance [km]')
ylabel('Depth [km]')
%axis ([-240 240 -100 5])
title('Depletion ')
hold on;


drawnow
hold on;


