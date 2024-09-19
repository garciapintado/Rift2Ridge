% Elena Ros 

function [sill,TP_xmelt,TP_ymelt,dF_column3,xsort,ysort] = columns_melt_sill...
    (GCOORD,dF,Point_id,area_melt,dF_interp,ext_rate,dt,ELEM2NODE,Phases,E2all,istep, GEO)

    % ext_rate :: [m/s] full extension rate
    km = 1000.;
    
    % Define DIMENSIONS of irregular grid where the sill will be located
    
    % JGP: ?? these are the new melt generation coordinates - not sill allocation
    xminmelt = min(GCOORD(1,dF>0)) - 0.5*km;                               % figure(); plot_tF(dF, GCOORD, ELEM2NODE, [], "dF", 0); 
    xmaxmelt = max(GCOORD(1,dF>0)) + 0.5*km;                               % this does not make sense!!! Why assymetric?
    yminmelt = min(GCOORD(2,dF>0)) - 10*km;
    ymaxmelt = max(GCOORD(2,dF>0)) +  5*km;
        
% Select nodes in the Moho
       Moho_nodes = find(Point_id==6);  %FOR POINT_ID==6 distances among nodes are different

        
% Identify the ends of the melting area to relate them to the nodes in the
% Moho

       [~,xa_node] = min(abs(GCOORD(1,Moho_nodes)-xminmelt));
       [~,xb_node] = min(abs(GCOORD(1,Moho_nodes)-xmaxmelt));
       
% Select the nodes in the Moho between the two ends Xa and Xb
 
       moho_indx_xa_xb= GCOORD(1,Moho_nodes)>=GCOORD(1,Moho_nodes(xa_node)) & ...
                        GCOORD(1,Moho_nodes)<=GCOORD(1,Moho_nodes(xb_node));
         
% Define x and y coordinates of the mesh over the melting to compute the sill

       xmesh = [sort(GCOORD(1,Moho_nodes(moho_indx_xa_xb)))] ;  %!!! to create the mesh I sort the nodes of the Moho to get the related dF_column in the same order
       ymesh = linspace(yminmelt,ymaxmelt,length(xmesh));
       [Xmesh,Ymesh] = meshgrid(xmesh,ymesh);
       
% Interpolate dF for cases with more than a melt point, otherwise it is not
% posible to do the calculus of the area_melt 
% (CHECK THIS AGAIN!!!!)

      is_dF = find(dF>0);
      if length(is_dF) == 1
%         area_melt = area_melt; % the one of the melt column
         dF_column3 = dF_interp(is_dF); % the one of the melt column
%          [Crust_thickness,xdike,ydike,indike,tp_melt,TP_xmelt,TP_ymelt] = column_melt_dike_maxE2all(GCOORD,...
%           area_melt,ext_rate,dt,Point_id,ELEM2NODE,Phases,E2all,istep,Crust_thickness);
         [dike,TP_xmelt,TP_ymelt] = column_melt_dike_maxE2all_newmesh(GCOORD,...
          area_melt,ext_rate,dt,Point_id,ELEM2NODE,Phases,E2all,istep);
      
         sill = dike;
         
         xsort = NaN;
         ysort = NaN;
      
      elseif length(is_dF)>1         
         dF_interp = griddata(GCOORD(1,:),GCOORD(2,:),dF,Xmesh,Ymesh,'linear');
         boolean = isnan(dF_interp);
         dF_interp(boolean) = 0 ;  
         
% Areas
         x_left = xmesh(1:end-1);
         x_right = xmesh(2:end);
         delta_x = x_right-x_left; 
         delta_y = ymesh(end)-ymesh(end-1);
         delta_xx = repmat(delta_x,length(xmesh),1);
         delta_yy = repmat(delta_y,length(xmesh),length(delta_x));
         area_xy_node = delta_xx.*delta_yy ; 
         
% dF in each column in relation to the corresponding Moho node
% (between the two ends Xa and Xb)
         dF_=dF_interp;
         dF_(:,end)=[]; 
         dF_column = sum(area_xy_node.*dF_,1);
         
% Crustal thickness
         area_melt = sum(dF_column);    
         %dF_column2 = dF_column/(ext_rate*dt);% like in the dike
         dF_column3 = dF_column./delta_x;
    
      
        dz_crust = area_melt/(ext_rate*dt);
        Crust_thickness(istep) = dz_crust/2;
      

% Index of the nodes in moho_indx_xa_xb==1 and remove the node 
% in the end Xb to preserve size

       mo_nodes = find(moho_indx_xa_xb==1); 
    
% The sorted indexes of xsill are needed to order coordinate y 
% as well, and remove the corresponding node Xb

       [xsort,xsill_indx] = sort(GCOORD(1,Moho_nodes(moho_indx_xa_xb)));
       ysort = GCOORD(2,Moho_nodes(mo_nodes(xsill_indx)));
       ysort(:,end) = [];
    
       xsort(:,end) = [];  
% Calculate the thickness the columns of the sill

    %dz_thick = ysort-dF_column;  %dz_crust in each Moho node
       dz_thick = ysort-dF_column3;  %dz_crust in each Moho node
          
% Calculate the the points forming each column of dz_thick
   
     for i=1:length(mo_nodes)-1
       
       if (dF_column3(i))<10
          scale = 1;
       else
          scale = 10;
       end    
       
       column_melt = [xsort(i)*ones(1,ceil(dF_column3(i)/scale));...
                     linspace(ysort(i),dz_thick(i),ceil(dF_column3(i)/scale))];
                 
       [TP_xmelt,TP_ymelt] = meshgrid(column_melt(1,1), linspace(ysort(i),...
                             dz_thick(i),ceil(dF_column3(i)/scale)));
                         
       if i==1
         TP_xm = TP_xmelt';
         TP_ym = TP_ymelt';
       else
         TP_xm = [TP_xm,TP_xmelt'];
         TP_ym = [TP_ym,TP_ymelt'];  
       end

      xsill = [xsort(1)',xsort,xsort(end)',fliplr(xsort)] ;
      ysill= [ysort(1)',ysort,ysort(end)',fliplr(dz_thick)];
      insill = inpolygon(GCOORD(1,:),GCOORD(2,:),xsill,ysill);
   
      TP_xmelt = TP_xm;
      TP_ymelt = TP_ym;

   end
  end
end