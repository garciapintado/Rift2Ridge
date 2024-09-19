function TopoXY = diff_topo_sealand(Basement, TopoXY, GCOORD, SP, dt, zuco)
    % [TOPOXY,SP] = DIFF_TOPO_ARMITAGE(TOPOXY,GCOORD,SP,DT,MA) takes current
    % topography TOPOXY and calculates new topographies by diffusion taking
    % into account the parameters defined by GCOORD, SP, DT and MA. This
    % algorithm uses the finite element diffusion function
    % function_nl_erosion_fixedelevation developed by John Armitage.
    % ALL units must be in International System.

    maxslope    = max(atand(abs(diff(TopoXY(2,:))./diff(TopoXY(1,:)))));       % [\deg] absolute maximum slope
    cwide       = max(10., (500-10)/(25-4)*maxslope - 160.);                   % [m]
    sres        = 10;                                                          % [m]
    highan      = 20;
    highan_air  = 25;

    %==========================================================================
    % SCALING
    %==========================================================================
    model_long  = diff(minmax(GCOORD(1,:)));                                   % [m] longitude of the model
    year = 365.25*24*60*60;
    ma   = 1e6*year;                                                        % [s]
    
    % Make variables adimensional
    if SP.scale                                                                % local SP. Scaling factors: model_long [L], kappa...
        error('scaling not operational')
        kappa = SP.kappa;
        SP.kappa    = 1.0;
        SP.kdecay   = SP.kdecay   * model_long;                                % ?
        SP.kdistalb = SP.kdistalb * year / kappa;
        SP.kdistals = SP.kdistals * year / kappa;
        SP.ktidal   = SP.ktidal   * year / kappa;
        SP.lriver   = SP.lriver   / model_long;                                % [m] -> []
        SP.q_bc     = SP.q_bc     / model_long;                                % [m/s] -> [/s]
        SP.sealevel = SP.sealevel / model_long;                                % [m] -> []
        
        SP.De = SP.c * SP.alpha_sed^SP.nexp * (model_long / kappa);
        Dx_sed = diff(TopoXY(1,:));
        Dx_sed = Dx_sed'/model_long;
        dt_sed = SP.dt * kappa / model_long^2;
        
        Topography =  TopoXY(2,:)' / model_long;                               % [m] -> []
        
        clear("kappa")
    else
        SP.De = SP.c * SP.alpha_sed^SP.nexp;
        dt_sed = SP.dt;
    end
    
    Old_topo = TopoXY;
    MM_SED = [1/3 1/6; 1/6 1/3]; %1/9*([16/5 4/5; 4/5 16/5]+[5/8 5/8; 5/8 5/8]); %[1/3 1/6; 1/6 1/3]; % [1/4 1/4; 1/4 1/4];
    
    ntsteps_sed = ceil(dt/SP.dt);                                              % Number of sediment time steps to reach a mechanical time step
    last_dts = dt - (ntsteps_sed-1)*SP.dt;

    % Loop through the erosion/sedimentation time steps
    for j = 1:ntsteps_sed
        % disp("diss_topo_sealand:: it:" + j);
        % Remesh steep slopes in the sea
        %-------------------------------
        if SP.ktidal ~= 0
            % Find submarine nodes
            subm_nodes = TopoXY(2,:) < SP.sealevel;                                % LOGICAL [1,ntopo]
            % Find the steep slopes
            steep_slopes = ...
                atand(abs(diff(TopoXY(2,:))./diff(TopoXY(1,:))))>=highan;
            steep_slopes = ([steep_slopes 0]+[0 steep_slopes]) ~= 0;               % LOGICAL [1,ntopo]
            % Find subaerial steep slopes
            steep_slopes_air = ...
                atand(abs(diff(TopoXY(2,:))./diff(TopoXY(1,:)))) >= highan_air;
            steep_slopes_air = ([steep_slopes_air 0]+[0 steep_slopes_air]) ~= 0;   % LOGICAL [1,ntopo]
            % Find shallow water and near sea subaerial topographies
            near_sea = TopoXY(2,:) >= SP.sealevel-cwide/2 & ...                    % LOGICAL [1,ntopo]
                TopoXY(2,:) <= SP.sealevel+cwide/2;
            
            % x-resolution check
            need_remesh = diff(TopoXY(1,:)) > 2*sres;
            need_remesh = ([need_remesh 0]+[0 need_remesh]) ~= 0;                  % LOGICAL [1,ntopo]
            % Submarine steep slopes
            sbm_ss = ...                                                           % LOGICAL [1,ntopo]
                ((subm_nodes & steep_slopes) | near_sea | steep_slopes_air) ...
                & need_remesh;
            % Calculate steep slope segment indexes
            counter_sbm_ss = 1:length(sbm_ss);
            sbm_ss_seg = sbm_ss*1;
            sbm_ss_seg(sbm_ss) = ...
                [1 cumsum(diff(counter_sbm_ss(sbm_ss))>1)+1];
            % Number of segments
            nseg = unique(sbm_ss_seg(sbm_ss_seg~=0));
            
            %         % Plot
            %         figure(2); clf
            %         plot([TopoXY(1,1) TopoXY(1,end)],[SP.sealevel SP.sealevel],'--')
            %         hold on
            %         plot(TopoXY(1,:),TopoXY(2,:))
            %         plot(TopoXY(1,subm_nodes),TopoXY(2,subm_nodes),'x')
            
            % Loop through steep slope submarine segments
            for n = 1:size(nseg,2)
                % Indexes of the nodes in the current segment
                seg_nodes = find(sbm_ss_seg==n);
                xleft   = TopoXY(1,seg_nodes(1));
                xright  = TopoXY(1,seg_nodes(end));
                % Remesh x
                xslope = linspace(xleft,xright,round((xright-xleft)/sres));
                % Remesh y
                yslope = interp1(TopoXY(1,seg_nodes),TopoXY(2,seg_nodes),xslope);
                % Add remeshed nodes to topography
                ileft = 1:length(sbm_ss_seg)<seg_nodes(1);
                iright = 1:length(sbm_ss_seg)>seg_nodes(end);
                TopoXY = [TopoXY(:,ileft) [xslope; yslope] TopoXY(:,iright)];
                sbm_ss_seg = [sbm_ss_seg(ileft) (n*ones(length(xslope),1))' ...
                    sbm_ss_seg(iright)];
            end
            %         plot(TopoXY(1,sbm_ss_seg~=0),TopoXY(2,sbm_ss_seg~=0),'o')
        end
        
        % Check if corner nodes are present and if not add them back
        deleted_corners = TopoXY(1,[1 end])~=Old_topo(1,[1 end]);
        if deleted_corners(1)
            TopoXY = [Old_topo(:,1) TopoXY];
        end
        if deleted_corners(2)
            TopoXY = [TopoXY Old_topo(:,end)];
        end
        
        Dx_sed = diff(TopoXY(1,:))';                                           % [ntopo-1,1]
        Topography = TopoXY(2,:)';                                             % [ntopo,1]
        ntopo = length(Topography);
        
        %Topo_base = interp1(TopoXY(1,:),TopoXY(2,:),Basement(1,:)); %
        %Basement = [Basement(1,:); ...
        %            min([Basement(2,:); Topo_base])];                         % [m] possibly eroded basement level in its original x-coordinates
        sedthkXZ = get_sed_thickness(Basement, TopoXY);                        % [m] [2,ntopo] sediment thickness at topography-x coordinates
        % plot(TopoXY(1,:)/1000,TopoXY(2,:)/1000,'Color',"blue");
        % hold on;
        % plot(Basement(1,:)/1000,Basement(2,:)/1000,'Color',"red");
        % plot(sedthkXZ(1,:)/1000,sedthkXZ(2,:)/1000,'Color',"orange");
        
        source = zeros(ntopo,1);                                               % [ntopo,1]
        
        if ~isempty(SP.pelagic_rate)
            belowSea = TopoXY(2,:)' < SP.sealevel;                             % [ntopo,1]
            Hw = reshape(SP.sealevel - TopoXY(2,belowSea),[],1);               % [nbelowsea,1] water depth (thickness)
            Hwini = 10;                                                        % [m] minimum depth for pelagic sources to reach SP.pelagic_rate. This is mostly to avoid a step with respect to aerial nodes
            if isfield(SP,'pelagic_Hmax') && ~isempty(SP.pelagic_Hmax)         % and to avoid replenishment of developing basins. 
                Hwmax = SP.pelagic_Hmax;
            else
                Hwmax = max(Hw);
            end
            Hwmax = max(Hwini,Hwmax);
            if SP.pelagic_Hdep
                Beta_sea = min(Hw ./ Hwmax,1.0);                           % [nbelowsea,1]
            else
                Beta_sea = min(Hw ./ Hwini,1.0);                           % [nbelowsea,1]
            end
            belowCCD = Hw >= SP.ccd;                                       % LOGICAL [nbelowsea,1], deeper than the carbonate compensation depth for belowSea nodes
            pelagic_rate = repelem(SP.pelagic_rate,sum(belowSea),1);       % [nbelowsea,1]
            pelagic_rate(belowCCD) = SP.pelagic_rate_below_ccd;
            source(belowSea) = ...
                source(belowSea) + pelagic_rate .* Beta_sea;               % [ntopo,1] [m/s]
        end
        clear("belowSea","Hw","Beta_sea","belowCCD","pelagic_rate")
        
        if j == ntsteps_sed
            dt_sed = last_dts;
            if SP.scale
                dt_sed = last_dts*SP.kappa/(model_long*model_long);
            end
        end
        
        % Calculating new topography
        if SP.LIVDcomp
            SP.q_bc = - (Topography([1 end]) - zuco(:))/dt - source([1 end]);
        end
        
        Topography = topographical_diffusion_smooth(Dx_sed, Topography, source, sedthkXZ(2,:)', dt_sed, SP);
        % figure(); plot(TopoXY(1,:), Topography, 'color',[.5 .5 .5])
        % hold on;  plot(TopoXY(1,:), Topography0, 'color','blue')
        % hold on;  plot(TopoXY(1,:), Topography1, 'color','yellow')
        % hold on;  plot(TopoXY(1,:), source/max(source), 'color','brown')
        % hold on;  plot(TopoXY(1,:), sedthkXZ(2,:)/max(sedthkXZ(2,:)), 'color','green')
        TopoXY = [TopoXY(1,:); Topography'];                               % hold on; plot(TopoXY(1,:)/1000,Topography/1000,'r')
        
    end % time loop

    % Back to old mesh
    if SP.ktidal ~= 0
        TopoXY = [Old_topo(1,:); ...                                       % [2,:]
                  interp1(TopoXY(1,:),TopoXY(2,:),Old_topo(1,:))]; 
    end

    % Redimensioning topography
    if SP.scale
        Topography = Topography * model_long;
    end
end % function


% % Plot old and new topographies (uncomment)
% figure(2)
% plot(TopoXY(1,:)/1000,Old_topo/1000,'k')
% hold on
% plot(TopoXY(1,:)/1000,Topography/1000,'r')
