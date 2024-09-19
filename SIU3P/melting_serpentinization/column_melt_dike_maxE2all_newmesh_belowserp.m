function [dike, TP_xmelt, TP_ymelt, x_shallow, y_shallow, radius_tpmelt_all] = ...
                column_melt_dike_maxE2all_newmesh_belowserp(GCOORD, ELEM2NODE, GEO, Phases, area_melt, dz_crust, ext_rate, dt,  ...
                E2all, Dserp)
    % Author: Elena Ros
    %
    % ext_rate: full spreading rate           
    %
    % GCOORD :: REAL, OPTIONAL, DIM(2,nnod) node physical coordinate
    % 
    % OUTPUT
    % ---
    % dike :: STRUCT with elements
    %     .coo ::     REAL, DIM(2,:) x, y coordinates of polygon defining the dike
    %     .nodesin :: LOGICAL, DIM(1,:) true for mesh nodes inside the
    %                 polygon defined by the dike
    % TP_xmelt  :: REAL, DIM(2,ntrk) where 'ntrk' is the number of returned tracking points
    % TP_xmelt  :: REAL, DIM(2,ntrk)
    % x_shallow :: REAL x of the top center of allocated dike
    % y_shallow :: REAL y of the top center of allocated dike
    % radius_tpmelt_all :: REAL (1,ntrk) radius of tracking point
    %                      associated melting areas, such that
    %                      sum(radius_tpmelt_all) == area_melt
    %
    % WARNING: this function is adapted to three layers: so moho is
    % identified as GEO(3).pids
    
    % Javier Garcia-Pintado: simplified, documented  and minor modifications. Same logic
    
    
    % Find nodes in the Moho and select the elements of the upper mantle in which the Moho is present.
    Moho_nodes = GEO(3).pids;
    mohoe = find(sum(ismember(ELEM2NODE,Moho_nodes),1) > 1 & Phases==1);

    % Find element index in mesh with the maximum strain rate, touching the Moho and highest vertical coordinate 
    [~,maxeli] = max(max(E2all(mohoe,:),[],2));
    elid = mohoe(maxeli); 
    inodes = ELEM2NODE(:, elid);                                               % nodes of this element

    [~,ilnode] = max(GCOORD(2,inodes));                                        % local node index
    inode = inodes(ilnode);                                                    % mesh node index where melt is going to be allocated
    xy = GCOORD(:,inode);                                                      % inode coordinates

    x_shallow = xy(1);                                                     % [m]
    y_shallow = xy(2);                                                     % [m] figure(); plotGEO(GEO); hold on; plot(xy(1)/1000,xy(2)/1000,'o','markersize',20)
    
    %dwidmax = ext_rate * dt; % maximum dike width
    
    %xy = [x_shallow + 10* dwidmax/2*[-1 1]; ...
    %      y_shallow   y_shallow]; 
    %dwid = diff(griddata(GCOORD(1,:),GCOORD(2,:),Vel(1,:) * dt,xy(:,1),xy(:,2),'linear'));
    
    if max(Dserp) >= 0.001
            
        Dserp_test = Dserp;
        Dserp_test(Dserp < 0.001) = 0.;                                    %  
        
        %%Dserp_test(nodes_Moho_maxE2all(i_coord)) % can happen that the point itself does not have Dserp
        %any(Dserp(ELEM2NODE(:,max_eri_i))~=0)  % if equals 1 means that this triangle has Dserp
        %Dserp_nodes= find(Dserp_test>0);
        %Dserp_elems= find(sum(ismember(ELEM2NODE,Dserp_nodes),1)& Phases==1);
        % %plot(GCOORD(1,Dserp_nodes)/1000,GCOORD(2,Dserp_nodes)/1000,'*y')
        % %Tris   = tsearch2(GCOORD,uint32(ELEM2NODE(1:3,:)),GCOORD(:,Dserp_nodes));
        % %plot(GCOORD(1,ELEM2NODE(:,Tris))/1000,GCOORD(2,ELEM2NODE(:,Tris))/1000,'*m')
        
        % Dserp_nodes_uni= unique(ELEM2NODE(:,Dserp_elems));

        % % find(abs(GCOORD(1,Dserp_nodes_uni)-x_shallow)==min(abs(GCOORD(1,Dserp_nodes_uni)-x_shallow)) & ...
        % %      max(abs(GCOORD(2,Dserp_nodes_uni)-y_shallow)))
        % 
        % % [val,i_node]= min(GCOORD(2,nodes_elem_maxE2all))
        % % plot(GCOORD(1,nodes_elem_maxE2all(i_node))/1000,GCOORD(2,nodes_elem_maxE2all(i_node))/1000,'*m','MarkerSize',14)

        
        % %  [ma,~] = find(E2all(Dserp_elems,:)==max(max(E2all(Dserp_elems,:))));
        % %  mai = Dserp_elems(ma); 
        
        [X_sha,Y_sha_depth] = meshgrid(x_shallow, y_shallow + (0:-10:-20000));
        
        %Dserp_nodes_no=find(Dserp_test==0);
        % [X_serp1,Y_serp1] = meshgrid(linspace(x_shallow-1000,x_shallow+1000), linspace(y_shallow,-20000,100));
        %%%%Dserp_test_interp = griddata(GCOORD(1,Dserp_nodes_uni),GCOORD(2,Dserp_nodes_uni),Dserp_test>0,X_serp1,Y_serp1,'linear');
        %Dserp_test_interp = griddata(GCOORD(1,:),GCOORD(2,:),Dserp_test,X_serp1,Y_serp1,'linear')
        %plot(X_serp1(Dserp_test_interp>0)/1000,Y_serp1(Dserp_test_interp>0)/1000,'y*')

        Dserp_test_interp = griddata(GCOORD(1,:),GCOORD(2,:),Dserp_test,X_sha,Y_sha_depth,'linear');
        % figure(); plot_tF(Dserp, GCOORD, ELEM2NODE, [], "Dserp", 0); hold on; plotGEO(GEO)
        % hold on;  plot(X_sha(Dserp_test_interp>0)/1000,Y_sha_depth(Dserp_test_interp>0)/1000,'go')
        
        if max(Dserp_test_interp) ~= 0
            y_shallow = min(Y_sha_depth(Dserp_test_interp>0));
        end
        % [val,ii]=find(abs(GCOORD(1,Dserp_nodes_no)-x_shallow)== ...
        %      min(abs(GCOORD(1,Dserp_nodes_no)-x_shallow)))

        % for i=1:length(Y_sha_depth)
        % [val,ii]=find(abs(GCOORD(2,Dserp_nodes_no)-Y_sha_depth(i))== ...
        %      max(abs(GCOORD(2,Dserp_nodes_no)-Y_sha_depth(i))))
        % end
        %     for i=1:length(Dserp_nodes)
        %         [val, els(Dserp_nodes(i))] = min(sqrt((x_shallow - GCOORD(1,Dserp_nodes(i))).^2 ));
        %     end
    end % max(Dserp) >= 0.001



    y_bottom = y_shallow - dz_crust;                                        % [m] y-coordinate of the bottom of dike
    
    % a] dike modelling as tracking points
    vres = min(50, dz_crust);                                              % [m] approximate vertical spacing between tracking nodes representing the dike 
    
    ext_hrate = ext_rate / 2.;                                             % [m/s] half extension rate 
    dwid = ext_rate * dt;                                                  % total dike width; (area_melt == dwid*dz_crust)
    
    ny = round(dz_crust/vres);                                             % number of vertical nodes to represent the dike
    ntrk = 2 * ny;                                                         % number of total tracking points representing the dike 
    
    radius_tpmelt = sqrt(area_melt/ntrk / pi);                             % [m] radius of each circle, such that the sum of areas represents area_melt
    radius_tpmelt_all = repmat(radius_tpmelt,1,ntrk);
    
    ulx = (x_shallow - dwid/2);  % x upper-left corner of dike
    urx = ulx + dwid;            % x upper-right corner of dike
    dike_left =  [ulx * ones(1,ny);...              % [m] (2,:) position of left nodes of dike
                  linspace(y_shallow, y_bottom, ny)];
    dike_right = [urx * ones(1,ny);...              % [m] (2,:) position of left nodes of dike
                  linspace(y_shallow, y_bottom, ny)];

    TP_xmelt = [dike_left(1,:) dike_right(1,:)];                           % (ntrk,1)
    TP_ymelt = [dike_left(2,:) dike_right(2,:)];                           % (ntrk,1)
    
    % b] dike modelling for detecting mesh nodes inner to dike polygon
    nx = floor(dwid/(2*10));                                               % number horizontal nodes in modelled dike
    ny = floor(dz_crust/(2*10));                                           % number vertical nodes in modelled dike
    
    line_bottom = [linspace(x_shallow-ext_hrate*dt,x_shallow+ext_hrate*dt,nx); ... % left -> right
                   y_bottom * ones(1,nx)];
    
    line_top = line_bottom;                                                        
    line_top(2,:) = line_bottom(2,:) + dz_crust;
    line_top = fliplr(line_top);                                           % right -> left
               
    line_left = [(x_shallow - dwid/2) * ones(1, ny);...
                 linspace(y_shallow, y_bottom, ny)];                        % up -> bottom
             
    line_right = line_left;                                                
    line_right(1,:) = line_left(1,:) + dwid;
    line_right = fliplr(line_right);                                       % bottom -> up         
    
    dike.coo     = [ line_bottom(:,1:end-1) line_right(:,1:end-1) line_top(:,1:end-1) line_left(:,1:end-1)];
    dike.nodesin = inpolygon(GCOORD(1,:),GCOORD(2,:), dike.coo(1,:), dike.coo(2,:)); % plotGEO(GEO); hold on; plot(dike.coo(1,:)/1000,dike.coo(2,:)/1000,'x')
end % hold on; plot_meshF(ELEM2NODE,GCOORD); hold on; plot(GCOORD(1,dike.nodesin)/1000, GCOORD(2,dike.nodesin)/1000, 'o','markersize',20)


