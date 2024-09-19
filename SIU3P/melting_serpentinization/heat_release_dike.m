function [Temp,TRACKP_melt,cont_melt,tp_melt, xdike, ydike]  = heat_release_dike(TRACKP_melt,cont_melt,GCOORD,...
                        Temp,Geo_id,dz_crust,K,Cp,ext_rate,dt,...
                        Point_id,ELEM2NODE,Phases,E2all,km)


if(dz_crust==0)
  TRACKP_melt{end+1} = [NaN'; NaN']; 
  xdike=NaN;
  ydike=NaN;                      
else
  cont_melt = cont_melt+1;
  [column_melt,tp_melt,TP_xmelt,TP_ymelt,x_shallow,y_shallow] = column_melt_maxE2all(GCOORD,ext_rate,...
            dt,dz_crust,Point_id,ELEM2NODE,Phases,E2all);
  if cont_melt==1
     TRACKP_melt{end+1} = [TP_xmelt(:)'; TP_ymelt(:)']; 
  else
     TRACKP_melt{end+1} = [TRACKP_melt{end}(1,:), TP_xmelt(:)'; TRACKP_melt{end}(2,:), TP_ymelt(:)'];   
  end
  
 % Magnitudes
   Rho_melt=2900; %kg/m3
   %L_Cp= 320;     %ratio L/Cp in grades C 
   T_melt=1100.;  %Assumption for Temperature of melt once is upwelled
   %L_melt=334e3;
   %kappa=K(1)/(Rho_melt*Cp(1)); %thermal diffusivity=thermal conductivity/(density*heat capacity)
               
 % Area of heat release
   line_bottom = [linspace(x_shallow-ext_rate*dt,x_shallow+ext_rate*dt,...
                        floor(ext_rate*dt/(2*10))); (y_shallow-dz_crust)...
                        *ones(1,floor(ext_rate*dt/(2*10)))];
   line_right=[(x_shallow+ext_rate*dt)*ones(1,floor(dz_crust/(2*10)));...
                         linspace(y_shallow,y_shallow-dz_crust,floor(dz_crust/(2*10)))];
   line_top = [linspace(x_shallow-ext_rate*dt,x_shallow+ext_rate*dt,...
                        floor(ext_rate*dt/(2*10))); y_shallow*ones(1,floor(ext_rate*dt/(2*10)))];
   line_left=[(x_shallow-ext_rate*dt)*ones(1,floor(dz_crust/(2*10)));...
                         linspace(y_shallow,y_shallow-dz_crust,floor(dz_crust/(2*10)))];  
   
%%%%%%%%%%%%%%%%    
    xdike = [ line_left(1,:),line_top(1,:),line_right(1,:),fliplr(line_bottom(1,:))] ;
    ydike = [ line_left(2,:),line_top(2,:),line_right(2,:),fliplr(line_bottom(2,:))];
    %indike = inpolygon(GCOORD(1,:),GCOORD(2,:),xdike,ydike);

  clear 'line_bottom' 'line_right' 'line_top' 'line_left'
%%%%%%%%%%%%%%%% 


% T_o has to be the interpolation of Temp into the y_points_release
% (background temperature approximation 1D)
               
    y_points_release = linspace(y_shallow,y_shallow-dz_crust,50);%100
    [xq,yq] = meshgrid(x_shallow,y_points_release);
    T_o = griddata(GCOORD(1,:),GCOORD(2,:),Temp,xq,yq);
    
% Solidification dike laplace in direction X and direction Y
   ifig = 3;
   dist1d = -400:40:400;  %distance in which apply the laplace solution
   %TEMP FUNCTION VS X            
   for nn=1:length(T_o)
       [T_x,x_m] = solidification_dike_laplace(K,Rho_melt,Cp,T_o(nn),T_melt,Geo_id,ext_rate*dt,dist1d,ifig);
       T_xx{nn} = T_x;
   end
%    ifig = 4;
%    dist1d = -10000:1000:10000;
%    % TEMP FUNCTION VS Y
%    T_oy=T_o(length(T_o)/2);
%    [T_y,y_m] = solidification_dike_laplace(K,Rho_melt,Cp,T_oy,T_melt,Geo_id,dz_crust/2,dist1d,ifig);
   x_m_coord = x_shallow+x_m;
   %y_m_coord = y_shallow-(dz_crust/2)+y_m;
               
   [x_mesh,y_mesh] = meshgrid(x_m_coord,y_points_release);
       
   for i=1:length(T_xx) %100%
       for j=1:length(T_xx{end}) %201%
               Temp_from_dike(i,j)=T_xx{i}(j);
       end
   end
       
% Interpolation of new temperatures into the nodes enclosed by the area

  [g,h]=size(x_mesh');
  tot=g*h;
  x_mesh_column = x_mesh';
  x_mesh_one_colum = reshape(x_mesh_column,1,tot) ; 
                
  [g,h]=size(y_mesh');
  tot=g*h;
  y_mesh_column = y_mesh';
  y_mesh_one_colum = reshape(y_mesh_column,1,tot) ;
                
  Temp_from_dike_column = Temp_from_dike';
  Temp_from_dike_one_colum = reshape(Temp_from_dike_column,1,tot) ;%20100

  F = scatteredInterpolant(x_mesh_one_colum', y_mesh_one_colum',Temp_from_dike_one_colum'); 
                
  found_nodes_indx= GCOORD(1,:)>=(x_mesh(1,1)) & GCOORD(1,:)<=(x_mesh(1,end))&...
                                GCOORD(2,:)>=(y_mesh(end,1)) & GCOORD(2,:)<=y_mesh(1,1);      
                            
  found_nodes=find(found_nodes_indx==1);
  Temp_super_new = F(GCOORD(1,found_nodes), GCOORD(2,found_nodes));
  Temp(found_nodes)=Temp_super_new;
               


           
end
end
   
   
   
   