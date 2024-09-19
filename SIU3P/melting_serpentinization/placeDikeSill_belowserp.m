function [dike, sills, TP_xmelt, TP_ymelt, x_shallow, y_shallow, radius_tpmelt_all, intrusive_type] = ...
          placeDikeSill_belowserp(GCOORD, ELEM2NODE, GEO, Phases, area_melt, ext_rate, dt,  ...
          E2all, ErP, Dserp, Temp, PHY, SETTINGS, BDTat, dF, ocrust)
    % +++ purpose +++
    % place a dike, sill or a dike-sill combination. The functions assumes
    % - the top of the emplacement is at the location of highest total
    %   strainrate in the mesh
    % - the top of the emplacement is no higher than the Moho
    % - serpentinite act as a seal so the emplacement is not allow to pass
    %    vertically through it.
    %
    % INPUT   
    % ---
    % GCOORD    :: REAL, DIM(2,nnod) node physical coordinate
    % ELEM2NODE :: INTEGER DIM(nnodel,nel), mesh connectivity matrix
    % Phases    :: INTEGER, DIM(1,nel), phase indicator
    % area_melt :: REAL, [m2] total amount of melt that need to be allocated
    % ext_rate  :: full spreading rate [m/s]           
    % dt        :: REAL, timestep [s]
    % E2ll      :: REAL, DIM(nel,nip) strain rate field at quadrature points
    % Dserp     :: REAL, DIM(1,nnod) total serpentinization degree in (0,1) 
    % Dserp_blocking_threshold :: REAL serpentinization degree (0,1) that
    %              is considered for blocking upwelling magma
    % Temp      :: REAL, DIM(1,nnod) temperature field
    % BDT_Temp  :: REAL temperature at which obtain an isotherm consider as the BDT
    
    % OUTPUT
    % ---
    % dike :: empty, or STRUCT with elements
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
    
    % Javier Garcia-Pintado, 2021.
    %
    % Details
    % This function is branched from a function, by Elena Ros, which was to allocate dikes.
    % Here, I consider that a brittle-ductile transition acts as boundary for dike-sill emplacement
    % type. All that is above the BDT places as a dike. Below the BDT, the
    % remaining magma which has not been placed a dike is then placed as an
    % elliptical bodies, representing gabbro sills.
    %
    % Leila Mezri has brought the idea of using ErP as blocking barrier in
    % addition/alternatively to serpentinization blocking
    
    % Find nodes in the Moho and select the elements of the upper mantle in which the Moho is present.

    nnod = size(GCOORD,2);
    sills_a2b                = PHY.sills_a2b;                              % ellipse ratio betwwen axes for sill bodies
    BDT_temp                 = PHY.BDT_temp;                               % [ºC] initiation of BDT. taken as limit for boundary for dike-sill emplacement type
    Dserp_blocking_threshold = PHY.Dserp_blocking_threshold;               % [0/1] serpentinization degree considered as blocking for upwelling melt                                                      
    ErP_blocking_threshold   = PHY.ErP_blocking_threshold; 
    
    dz_crust = area_melt/(ext_rate*dt);  % [m]
    
    ELEM2NODE = uint32(ELEM2NODE);
    mohon = uint32(GEO(3).pids);                                           % moho nodes [only triangle vertices]
    mohoe = find(sum(ismember(ELEM2NODE,mohon),1) > 1 & Phases==1);        % mantle elements with an edge at the Moho

    % Find element index in mesh with the maximum strain rate, touching the Moho and highest vertical coordinate 
    % [~,elid0] = max(max(E2all,[],2));
    moho = GCOORD(:,mohon);
    nmoho = size(moho,2);
    E2alln = ipval2nodval(ELEM2NODE,E2all,SETTINGS.eids,true);
    E2allmoho = E2alln(mohon);
    wE2all = reshape(E2allmoho/max(E2allmoho),1,[]);               % [1,nmoho] weight due to E2all \in [0,1]
    
    if SETTINGS.emplace_idwmagma
        [~,cid] = max(dF(1:nnod));                                     % get center of magma chamber
        xym = GCOORD(:,cid);                                           % location of maximum melt generated at this step
        idist2m = 1 ./ sqrt((moho(1,:)-xym(1)).^2 + (moho(2,:)-xym(2)).^2);
        wdism = idist2m/max(idist2m);                                  % [1,nmoho]
    else
        wdism = zeros(1,nmoho);
    end 

    % hold on; plot(moho(1,:)/1000, moho(2,:)/1000,'color','white')
    % hold on; plot(ocrust(1,:)/1000, ocrust(2,:)/1000,'color','red')
    switch SETTINGS.emplace_at 
        case 'E2allmax_moho' % maximum deviatoric strainrate at moho
            wcost = (wE2all + wdism); % figure(); plot_ip_valF(E2all, GCOORD, ELEM2NODE, GEO(9).coo, [0,3.5E-14], "II invariant of strain rate");
            [~,mid] = max(wcost);
            xy = moho(:,mid);          % hold on; plot(xy(1)/1000,xy(2)/1000,'x','color','cyan','markersize',15)
                                       % figure(); plot(GCOORD(1,mohon),wcost,'color','blue')
                                       % hold on; plot(GCOORD(1,mohon), wdism,'color','red')
                                       % hold on; plot(GCOORD(1,mohon), wE2all,'color','green')
        case 'BDT_nn_moho'    % nearest node in moho to estimated BDT coordinates
            idist2bdt = 1 ./ sqrt((moho(1,:)-BDTat(1)).^2 + (moho(2,:)-BDTat(2)).^2); % hold on; plot(BDTat(1)/1000,BDTat(2)/1000,'o','color','red','markersize',15)
            wdisbdt = idist2bdt/max(idist2bdt);                                  % [1,nmoho]
            wcost = (wdisbdt + wdism);  
            [~,mid] = max(wcost);
            xy = moho(:,mid);           % hold on; plot(xy(1)/1000,xy(2)/1000,'x','color','red','markersize',15)
                                        
        %case 'BDT_ortho_moho' % orthogonal projection from estimated BDT at moho - does not cling onto a node, perturbation of emplacement does not apply
        %    r = pointPolylineD(moho',BDTat');
        %    xy = [r.xc,r.yc];
        otherwise
            error("SETTINGS.emplace_at unknown case")
    end
    if PHY.emplace_sdev > 0 && exist('mid','var') % exact at a moho-node - perturb along the moho polyline
        chainage = [0 cumsum(diff2D(moho))];
        chainxy = chainage(mid) + normrnd(0,PHY.emplace_sdev);
        xy(1) = interp1(chainage,moho(1,:),chainxy);
        xy(2) = interp1(chainage,moho(2,:),chainxy);                       %       hold on; plot(xy(1)/1000,xy(2)/1000,'x','color','white','markersize',10)
    end
                    
    x_shallow = xy(1);                                                     % [m]
    y_shallow = xy(2);                                                     % [m] figure(); plotGEO(GEO,1:10,''); hold on; plot(xy(1)/1000,xy(2)/1000,'o','markersize',20)
    
    dwid = ext_rate * dt;                                                  % [m] total dike width; (area_melt == dwid*dz_crust)
    
    x_shallow_corners = x_shallow + dwid/2*[-1 1];
    y_shallow_corntop = interp1(GEO(end-1).coo(1,:),GEO(end-1).coo(2,:),x_shallow_corners);
    y_shallow = min(y_shallow, min(y_shallow_corntop) - 10.0);              
    
    %dwidmax = ext_rate * dt; % maximum dike width
    
    %xy = [x_shallow + 10* dwidmax/2*[-1 1]; ...
    %      y_shallow   y_shallow]; 
    %dwid = diff(griddata(GCOORD(1,:),GCOORD(2,:),Vel(1,:) * dt,xy(:,1),xy(:,2),'linear'));
    
    % BLOCKING CRITERIA
    [X_sha,Y_sha_depth] = meshgrid(x_shallow, y_shallow + (0:-10:-20000));         % top to bottom 
    
    % serpentinites
    if max(Dserp) >= 0.001
        Dserp = Dserp(1:nnod);
        Dserp(Dserp < 0.001) = 0.;                                         % just locally 
        
        % % find(abs(GCOORD(1,Dserp_nodes_uni)-x_shallow)==min(abs(GCOORD(1,Dserp_nodes_uni)-x_shallow)) & ...
        % %      max(abs(GCOORD(2,Dserp_nodes_uni)-y_shallow)))
        % 
        % % [val,i_node]= min(GCOORD(2,nodes_elem_maxE2all))
        % % plot(GCOORD(1,nodes_elem_maxE2all(i_node))/1000,GCOORD(2,nodes_elem_maxE2all(i_node))/1000,'*m','MarkerSize',14)

        
        % %  [ma,~] = find(E2all(Dserp_elems,:)==max(max(E2all(Dserp_elems,:))));
        % %  mai = Dserp_elems(ma); 
        
        Dserp_interp = griddata(GCOORD(1,:),GCOORD(2,:),Dserp,X_sha,Y_sha_depth,'linear'); % top-to-bottom serpentinization column for the location
        
        % figure(); plot(Y_sha_depth, Dserp_interp)
        
        % figure(); plot_tF(Dserp, GCOORD, ELEM2NODE, [], "Dserp", 0); hold on; plotGEO(GEO,1:10,'')
        % hold on;  plot(X_sha/1000,Y_sha_depth/1000,'go')
        % hold on;  plot(X_sha(Dserp_interp>0)/1000,Y_sha_depth(Dserp_interp>0)/1000,'go')
        
        if max(Dserp_interp) >  Dserp_blocking_threshold
            y_shallow = min(Y_sha_depth(Dserp_interp > Dserp_blocking_threshold)); % underplate serpentinized column
        end
    end % max(Dserp) >= 0.001                                              % figure(); plot_tF(Dserp, GCOORD, ELEM2NODE, [], "Dserp", 0); hold on; plotGEO(GEO,1:10,'')
    
    % faults
    if max(ErP,[],'all') > ErP_blocking_threshold
        [ipx,ipy] = ip_coord(GCOORD, ELEM2NODE(1:3,:), size(ELEM2NODE,2), 6);
        ErP_interp = griddata(ipx(:),ipy(:),ErP(:),X_sha,Y_sha_depth,'linear'); % top-to-bottom ErP column for the location
        if max(ErP_interp) > ErP_blocking_threshold
            y_shallow = min(y_shallow, min(Y_sha_depth(ErP_interp > ErP_blocking_threshold))); % underplate serpentinized column
        end
    end
    
    y_bottom = y_shallow - dz_crust;                                       % [m] y-coordinate of the bottom of dike | hold on; plot(repelem(x_shallow,1,2)/1000,[y_shallow,y_bottom]/1000,'color','green','linewidth',3)
    
    isoT = getIsolines(GCOORD, Temp(1:nnod), BDT_temp, true, true);        % hold on; plot(isoT.coo(1,:)/1000, isoT.coo(2,:)/1000,'--','color','white')
    %isoT900 = getIsolines(GCOORD, Temp(1:nnod), 900, true, true);          % hold on; plot(isoT900.coo(1,:)/1000, isoT900.coo(2,:)/1000,'--','color','yellow')
    y_BDT = interp1(isoT.coo(1,:),isoT.coo(2,:),x_shallow);                % [m] y-coordinate of the parametrized BDT | hold on; plot(x_shallow/1000,y_BDT/1000,'o','color','red')
    %Leila's alternative
    %[isoTX,isort] = sort(unique(isoT.coo(1,:)));                                % xsorted = x(isort)
    %isoTcoo = [isoTX; isoT.coo(2,isort)];
    %y_BDT = interp1(isoTcoo(1,:),isoTcoo(2,:),x_shallow);
    
    % alt. 1] gabro dike/basalt sills divide at BDT level  
    % y_bottom_dike = min(y_shallow, max(y_BDT, y_bottom));                  % [m] y-coordinate for the bottom of the dike
    % al1. 2] no gabbro/basalt division: all in a column [top would
    % represent dikes and lower would be vertically distributed dikes -
    % only the part that would exceed the base of the estimated ocrust at
    % the previous timestep
    if isempty(ocrust)
        y_bottom_dike = min(y_shallow, max(y_BDT, y_bottom));                % use alt. 1
    else
      y_ocrust = interp1(ocrust(1,:),ocrust(2,:),x_shallow);                 % ocrust criterion to improve visual matching of estimated oceanic
      y_bottom_dike = min(y_shallow, max(min(y_BDT,y_ocrust), y_bottom));
    end   
   
    dz_dike = y_shallow - y_bottom_dike;                                   % [m] >=0
    
    area_melt_dike = dwid * dz_dike;
    area_melt_sills = area_melt - area_melt_dike;
    
    radius_tpmelt_all = [];
    TP_xmelt = [];
    TP_ymelt = [];
    intrusive_type = [];
    dike = [];
    sills = [];
    
    if area_melt_dike > 0.0   
        % a] dike modelling as tracking points
        vres = min(50, dz_dike);                                               % [m] approximate vertical spacing between tracking nodes representing the dike 
    
        ny = round(dz_dike/vres);                                              % number of vertical nodes to represent the dike
        ntrk = 2 * ny;                                                         % number of total tracking points representing the dike at this time step
    
        radius_tpmelt = sqrt(area_melt_dike/ntrk / pi);                        % [m] radius of each tracking point (assumed as circle), such that the sum of areas represents area_melt
        radius_tpmelt_all = repelem(radius_tpmelt,1,ntrk);                 % Note: layer 1 in oceanic crust is on average 0.4 km thick, consisting of unconsolidated or semiconsolidated sediments - this is considered elsewhere
        intrusive_type = repelem(2,1,ntrk);                                % 2 : dikes (assumed as layer 2 in oceanic crust; commonly 2A [~0.5km pillow basalt] plus 2B [~1.5km diabase dikes])
        
        ulx = (x_shallow - dwid/2);                                        % x upper-left corner of dike
        urx = ulx + dwid;                                                  % x upper-right corner of dike
        dike_left =  [ulx * ones(1,ny);...                                 % [m] (2,:) position of left nodes of dike
                  linspace(y_shallow, y_bottom_dike, ny)];
        dike_right = [urx * ones(1,ny);...                                 % [m] (2,:) position of left nodes of dike
                  linspace(y_shallow, y_bottom_dike, ny)];

        TP_xmelt = [dike_left(1,:) dike_right(1,:)];                           % (ntrk,1)
        TP_ymelt = [dike_left(2,:) dike_right(2,:)];                           % (ntrk,1) hold on; plot(TP_xmelt/1000, TP_ymelt/1000,'o','color',[.7 .1 .1])
                                                                           % cooids =   % indices within the polygon
        % b] dike modelling for detecting mesh nodes inner to dike polygon
        nx = floor(dwid/(2*10));                                               % number horizontal nodes in modelled dike
        ny = floor(dz_dike/(2*10));                                           % number vertical nodes in modelled dike
        
        line_bottom = [linspace(x_shallow-ext_rate/2*dt,x_shallow+ext_rate/2*dt,nx); ... % left -> right
            y_bottom_dike * ones(1,nx)];
        
        line_top = line_bottom;
        line_top(2,:) = line_bottom(2,:) + dz_dike;
        line_top_topo = interp1(GEO(end-1).coo(1,:),GEO(end-1).coo(2,:),line_top(1,:));
        line_top(2,:) = min(line_top(2,:), line_top_topo - 10.0);          % assure no location is above topography 
        line_top = fliplr(line_top);                                           % right -> left
        
        line_left = [(x_shallow - dwid/2) * ones(1, ny);...
            linspace(y_shallow, y_bottom_dike, ny)];                     % up -> bottom
        
        line_right = line_left;
        line_right(1,:) = line_left(1,:) + dwid;
        line_right = fliplr(line_right);                                       % bottom -> up
        
        dike.xyc     = [x_shallow, y_shallow - dz_dike/2];
        dike.coo     = [ line_bottom(:,1:end-1) line_right(:,1:end-1) line_top(:,1:end-1) line_left(:,1:end-1)];
        dike.nodesin = inpoly2(GCOORD', dike.coo'); % plotGEO(GEO); hold on; plot(dike.coo(1,:)/1000,dike.coo(2,:)/1000,'x')
        dike.area_tp = pi * radius_tpmelt^2;
    end % dike modelling
    if area_melt_sills > 1.0E-01                                                   % neglect otherwise
        % c] sill modelling as [ellipse] tracking points
        nsills = prod(SETTINGS.sillmat);
        area_melt_sill = area_melt_sills/nsills;
        a = sqrt(area_melt_sill * sills_a2b / pi);                                 % ellipse semi-mayor axis A_ellipse = a*b*pi
        b = a / sills_a2b;
        
        ntrk = max(4, ceil(ellipsePerimeter(a,b)/50));                             % number of tracking points representing each ellipsoidal sill
        theta = linspace(0,2*pi,ntrk);                                             % [1,ntrk]
        XY = [cos(theta); sin(theta)];                                             % [2,ntrk] standarised generic circle
        
        nrow = SETTINGS.sillmat(1);
        ncol = SETTINGS.sillmat(2);
        x_centers = x_shallow + linspace(-1,1,ncol)*((ncol-1)*(a+25));             % 50-m spacing between ellipses
        y_centers = y_bottom_dike - b - 50 - linspace(0,2,nrow)*((nrow-1)*(b+25)); % 50 meters below the bottom of the dike + 50-m spacing between ellipses
        x_centers = repmat(x_centers,nrow,1);
        y_centers = repmat(y_centers(:),1,ncol);

        radius_tpmelt = sqrt(area_melt_sill/ ntrk / pi);                                % [m] radius of each tracking point (assumed as circle), such that the sum of areas represents area_melt
        radius_tpmelt_all = [radius_tpmelt_all, repelem(radius_tpmelt,1,ntrk*nsills)]; % [1,ntrk*nsills]
        intrusive_type = [intrusive_type, repelem(3,1,ntrk*nsills)];                   % [1,ntrk*nsills] 3: gabbros [assumed as layer 3 in oceanic crust]

        for i=1:nsills
            Xell = x_centers(i) + a(:) * XY(1,:);                                           % [1,ntrk] parametric equation for ellipse without rotation
            Yell = y_centers(i) + b(:) * XY(2,:);                                           % [1,ntrk] hold on; plot(Xell/1000, Yell/1000,'+')
            TP_xmelt = [TP_xmelt, Xell];                                          %
            TP_ymelt = [TP_ymelt, Yell];                                          %
            sills(i).xyc  = [x_centers(i),y_centers(i)];
            sills(i).coo  = [Xell; Yell];
            sills(i).nodesin = inpoly2(GCOORD', sills(i).coo');
            sills(i).area_tp = pi * radius_tpmelt^2;                              % [m2] area melt for each individual tracking point in the sill 
        end
    end % sill modelling
end % hold on; plot_meshF(ELEM2NODE,GCOORD); hold on; plot(GCOORD(1,dike.nodesin)/1000, GCOORD(2,dike.nodesin)/1000, 'o','markersize',20)

