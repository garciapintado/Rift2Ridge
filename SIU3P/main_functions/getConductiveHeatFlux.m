function r = getConductiveHeatFlux(GCOORD, ELEM2NODE, GEOn, Temp, cond_el, xy, h, dc, downscale)
    % CONDUCTIVE HEAT FLOW ALONG A POLYLINE
    % 
    % +++ purpose +++
    % This functions get the conductive heat flux across a polyline ($\mathbf{F}_c^T = - \lambda \nabla T$)
    % 
    % The nodes along the polyline a resampled, so that the returning
    % structure indicate a new resampled polyline along with the
    % corresponding heat flux. This function obtains a smooth kernel regression
    % along polyline locations.
    % The functions uses quadratic shape functions to sample the
    % temperature field to estimate the thermal gradients
    %
    % GCOORD    :: [2,nnod] 
    % ELEM2NODE :: [nnodel,nel]
    % GEOn      :: GEOMETRY struct
    % Temp      :: [1,nnod] [ºC] Temperature field
    % cond_el   :: [1,nel] [W.m-1.K-1] thermal conductivity
    % xy        :: [2,n] : input polyline
    % h         :: bandwidth for Gaussian kernel. Empty for automated bandwidth detection
    % dc        :: resampling spacing for downscaling. Empty for no kernel smoothing
    % downscale :: downscale [approximate] resolution to increase sampling previous to upscaling.
    %
    % OUTPUT
    % integral of heat flux along the polyline [W] (assumes the polyline represents a 1m-depth area in a 3D domain)
    % Javier GP, MARUM 2021
    
    % figure(); plot_meshF(ELEM2NODE, GCOORD);
    % figure(); plotGEO(GEOn,1:10,'')
    % hold on; plot(xy(1,:)/1000,xy(2,:)/1000,'.-','markersize',10,'color','red')
    
    if nargin < 7
        h = 500;
    end
    if nargin < 8
        dc = 500;
    end
    if nargin < 9
        downscale = [];
    end
    
    if ~isempty(downscale)    
        xy = resampleLine(xy, downscale, "linear");
    end
    
    % guarantee polyline is below top of domain
    ytop = interp1(GEOn(end-1).coo(1,:),GEOn(end-1).coo(2,:),xy(1,:));     % [1,n] assume GEO is monotonic
    xy = [xy(1,:); min(xy(2,:), ytop-0.001)];                              % [2,n] - 0.001 to avoid numerical errors in sample 2D | xy_chain = [0 cumsum(diff2D(xy))];
    grad_T6_xy = sampleGrad2D(GCOORD, ELEM2NODE(1:6,:), Temp, xy);         % [2,n] sample the gradient of the temperature field   | figure; plot(xy_chain, grad_T6_xy(2,:), 'color', 'black')
                                                                           
    % ALT 1]
    %[~, northo] = getOrthoDpolyline(xy, false);                            % upward-pointing orthogonal vectors at xy locations
    %grad_T6x_xyn = grad_T6_xy(1,:) .* cosd(northo);                        % project x-component onto normal to the surface
    %grad_T6y_xyn = grad_T6_xy(2,:) .* cosd(northo - 90);
    %grad_T6n_xy = grad_T6x_xyn + grad_T6y_xyn;                             % gradient projected onto the [upward] normal to the basement [at basement nodes] 
    if 1 > 2 % NOT RUN
        symsize = 500;
        figure(); plot(xy(1,:),xy(2,:),'.-','markersize',10);
        xy2 = [xy(1,:) + cosd(northo)*symsize; ...
               xy(2,:) + sind(northo)*symsize];
        for i=1:length(northo)                                              % orthogonal angles at edges
            hold on; plot([xy(1,i) xy2(1,i)],[xy(2,i) xy2(2,i)],'-','color','green'); 
        end
        %text(xy(1,:), xy(2,:)+0.01,string(round(northo,2)),'color','red'); 
        figure(); plot(xy(1,:),grad_T6_xy(1,:)); hold on; plot(xy(1,:),grad_T6x_xyn,'color',[.6 .0 .0]); yline(0,'color',[.3 .3 .3]);
        figure(); plot(xy(1,:),grad_T6_xy(2,:)); hold on; plot(xy(1,:),grad_T6y_xyn,'color',[.6 .0 .0]); yline(0,'color',[.3 .3 .3]);
        figure(); plot(xy(1,:),grad_T6_xy(2,:)); hold on; plot(xy(1,:),grad_T6n_xy,'color',[.6 .0 .0]);  yline(0,'color',[.3 .3 .3]);
    end
    
    % ALT 2]
    [~,alphasd] = getOrthoDpolyline(xy,1,0);                                % [1,n] average segment directions at vertices nodes in the xy polyline
    [~,grad_T6n_xy] = rotateParallel(grad_T6_xy(1,:)',grad_T6_xy(2,:)',-alphasd*pi/180); 
    grad_T6n_xy = grad_T6n_xy';                                            % [1,n]
    
    cond_nod = elval2nodval(ELEM2NODE(1:6,:), cond_el);
    cond_xy = sample2D(GCOORD, ELEM2NODE(1:6,:), cond_nod, xy);            % [1,n]
    
    heatflux_xy = - cond_xy .* grad_T6n_xy;                                % [W.m-2] = [W.m-1.K-1].[K.m-1] diffusive thermal energy per unit area per unit time
    
    chainage = [0. cumsum(diff2D(xy))];
    if isempty(dc)
       r = [];
       r.xy = xy;                                                          % [m]     [2,n]
       r.chainage = chainage;                                              % [m]     [1,n]
       r.heatflux = heatflux_xy;                                           % [W-m-2] [1,n]
       r.sensible_heat = sum(diff2D(xy) .* (heatflux_xy(1:end-1) + heatflux_xy(2:end))/2); % [W] total integral over the transect conducted over the un-smoothed samples
       return
    end
    
    chat = unique([0:dc:chainage(end), chainage(end)]); 
    chat = linspace(min(chat),max(chat),length(chat));
    xyat = [interp1(chainage,xy(1,:),chat); ...
            interp1(chainage,xy(2,:),chat)];
    
    if isempty(h)                                                          % estimate bandwidth [over input locations]
        rh = ksr(chainage, heatflux_xy);
        h = rh.h;
    end
    
    r = ksr(chainage, heatflux_xy, h, length(chat));                                          % smooth by kernel regression
    r.sensible_heat = sum(diff2D(xy) .* (heatflux_xy(1:end-1) + heatflux_xy(2:end))/2);       % [W] total integral over the transect conducted over the un-smoothed samples
    r.xy = xyat;                                                                              % [2,nresampled]
    %r.sensible_heat = sum(diff(r.x).* (r.f(1:end-1)+r.f(2:end)) / 2);                        % [W] total conductive sensible heat transfer along the transect [for each 1m width cross section]
    r = renameStructField(r,'x','chainage');                                                  % [m]     [1,nresamppled] 
    r = renameStructField(r,'f','heatflux');                                                  % [W-m-2] [1,nresampled] kernel-smoothed heat flux [for plotting purposes] 
    r.h = h; 
    % hold on; plot(xyat(1,:), r.f, 'color',[.8 .5 .0],'linewidth',2)
    % hold on; plot(xy(1,:), grad_T3y_xy , 'color','cyan','linewidth',2)
    % 
end
