% Elena Ros

function [Temp] = heat_release_dike_errorfunction(GCOORD,Temp,dz_crust,K,Cp,ext_rate,dt,T_melt,Rho_melt,x_shallow,y_shallow)

% T_n has to be the interpolation of Temp into the y_points_release
% (background temperature approximation 1D)
length_x= 1000;
y_points_release = linspace(y_shallow,y_shallow-dz_crust,50); %100
[xq,yq] = meshgrid(x_shallow,y_points_release);
T_n = griddata(GCOORD(1,:),GCOORD(2,:),Temp,xq,yq);
%scatter(xq/1000,yq/1000,50,T_n,'fill');hold on

% Postsolidification dike.
% Apply error function solution (Magnus-Wangen) in direction X
%dist1d = (-500:50:500);  %distance in which apply the  solution
dist1d = (0:50:length_x);  %distance in which apply the  solution
%plot(xq(1)/1000,y_points_release(1)/1000,'*r')

%kappa=1e-6
kappa = K(1)/(Rho_melt*Cp(1));
a=ext_rate*dt;
to = a^2/kappa;

xq_centered=0;
%          x_hat =(dist1d-xq(1))/a    %z_hat =(z-z1)/a
x_hat =(dist1d-xq_centered)/a;    %z_hat =(z-z1)/a
%plot(xq(1)/1000+dist1d/1000,repmat(y_points_release(1)/1000,[1,length(x_hat)]),'*')
t=1e2*365*24*60*60; %100 years
t_hat = (t-0)/to;

Temp_hat = (1/2)*(erf((x_hat+1)./(2*sqrt(t_hat)))  - erf((x_hat-1)/2*sqrt(t_hat)));
A = fliplr(Temp_hat);
A(:,end)=[];
Temp_hat2 = [A,Temp_hat];

dist1d2 = (-length_x:50:length_x);
%figure;scatter(dist1d2(:)/1000,repmat(y_points_release(1)/1000,1,length(dist1d2)),50,Temp_hat2(:),'fill');hold on

%T_cooled = To*T_hat;
%Temp_total = T_n + (T_melt - T_n)*Temp_cooled
Temp_total= zeros(length(T_n),length(Temp_hat2));
for i=1:length(T_n)%50
    %for j=1:length(Temp_hat) %21
    Temp_total(i,:) = T_n(i) + (T_melt - T_n(i))*Temp_hat2;
    %end
end

x_m_coord = x_shallow+dist1d2;
[x_mesh,y_mesh] = meshgrid(x_m_coord,y_points_release);
%scatter(x_mesh(:)/1000,y_mesh(:)/1000,50,Temp_total(:),'fill');hold on
               
       
% Interpolation of new temperatures into the nodes enclosed by the area

  size_mesh=size(x_mesh');
  x_mesh_one_colum = reshape(x_mesh',1,size_mesh(1)*size_mesh(2)); 
                
  %size_y=size(y_mesh');
  y_mesh_one_colum = reshape(y_mesh',1,size_mesh(1)*size_mesh(2));
                
  Temp_from_dike_one_row = reshape(Temp_total',1,size_mesh(1)*size_mesh(2));%20100
%plot_t_optimized_series % BEFORE
  F = scatteredInterpolant(x_mesh_one_colum', y_mesh_one_colum',Temp_from_dike_one_row'); 
                
  found_nodes_indx = GCOORD(1,:)>=(x_mesh(1,1)) & GCOORD(1,:)<=(x_mesh(1,end))&...
                                GCOORD(2,:)>=(y_mesh(end,1)) & GCOORD(2,:)<=y_mesh(1,1);      
                            
  found_nodes = find(found_nodes_indx==1);
  Temp_new = F(GCOORD(1,found_nodes), GCOORD(2,found_nodes));
  Temp(found_nodes) = Temp_new;

% plot_t_optimized_series % AFTER
%figure; scatter(GCOORD(1,found_nodes)/1000,GCOORD(2,found_nodes)/1000,50,Temp(found_nodes),'fill');colorbar          
end
   
   
   
   