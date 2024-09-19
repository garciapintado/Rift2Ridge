function [GCOORD, Point_id, Phases, ELEM2NODE, GEO, GEOn, ndof, nnod7, nel, ...
          F_xx, F_xy, F_yx, F_yy, TAU_xx_old, TAU_xy_old, TAU_yy_old, E2all, Mu_all, ...      % [nel,nip] variables  
          I2, SS, RHEOLvar, ...                                                               % cell/struct with [nel,nip] variables
          Temp, PRESS_CONT, Dpl, Dserp] = ...                                                        % [1,nnod7] variables
   remesh_all3l_dpl_rand_melt_I2_F_Mu_E2_linear(GCOORD, Phases, ELEM2NODE, GEO, ...
          Elsizes, nip, SETTINGS, PHY, triangle_path, triangle_mode, meshname, ...
          F_xx, F_xy, F_yx, F_yy, TAU_xx_old, TAU_xy_old, TAU_yy_old, E2all, Mu_all, ...      % [nel,nip] variables
          I2, SS, RHEOLvar, ...                                                               % cell/struct with [nel,nip] variables
          Temp, PRESS_CONT, Dpl, Dserp)                                                              % [1,nnod7] variables
    %      
    %
    % REMESH_ALL generates a new mesh and remeshes the temperature, the
    % accumulated gradient of the deformation and the stress
    %
    % INPUT [I]
    % SETTINGS.intersect_s      :: switch for function dealing with thin quasi-horizontal layers
    %         .layshift         :: clipping distance between interfaces 
    %         .tolerance        :: clipping distance within interface
    %         .temp_remesh_type :: switch for temperature remeshing
    % PHY                       :: struct with physical parameters. Melting parameter as needed, and input as
    % PHY.MELT
    % 
    % OUTPUT [O]
    % GEO  : structured GEOMETRY as input to generate_meshGEO()
    % GEOn : structured GEOMETRY as output from generate_meshGEO()
    
    % Author: Miguel Andres-Martinez, branched from MILAMIN
    %         Javier Garcia-Pintado, 2020-03 adapted to the use of GEO structs

    eids = "645";                                      % rift2ridge convention [required for kinedyn-shared trimesh_p2_to_p1()]
    check_monotony = false; % do not check or x-monotony in resampleGEO()
    
    E2N = ELEM2NODE;
    GCO = GCOORD;
    Pha = Phases;

    mi = 2 - [GEO.horizontal];                                               % monotonicity checking index
    
    if SETTINGS.intersect_s == "layer_collapse" % tmp: raise error for debugging - monotony should be dealt within layer_collapse()
        for i=1:length(GEO)
            if any(diff(GEO(i).coo(mi(i),:)) <= 0.0)
                error("remesh_all:: non monotonic interface i:" + i)
            end
        end
    end

    GEO = simplifyGEO(GEO, SETTINGS.tolerance, SETTINGS.layshift);
    GEO = resampleGEO(GEO, SETTINGS.GEOres, SETTINGS.layshift, check_monotony); 
    
    if SETTINGS.intersect_s == "layer_collapse" % tmp: raise error for debugging - monotony should be dealt within layer_collapse()
        for i=1:length(GEO)
            if any(diff(GEO(i).coo(mi(i),:)) <= 0.0)
                error("remesh_all:: simplify+resample made non monotonic interface i:" + i)
            end
        end
    end

    genmesh_t = tic;                             
    [GCOORD,ELEM2NODE, Point_id, Phases] = generate_meshGEO(GEO, ...       % figure(); plot_phasesF(ELEM2NODE, GCOORD, Phases);
            Elsizes, triangle_mode, triangle_path, meshname);              % hold on; plot_meshF(ELEM2NODE, GCOORD); hold on; plotGEO(GEO)
    
    nnod    = size(GCOORD,2);
    ndof    = size(GCOORD,1);
    nel     = size(ELEM2NODE,2);
    disp(['remesh_all3l_dpl_rand_melt_I2_F_Mu_E2_linear:: genmesh: ',num2str(toc(genmesh_t))]);

    %add 7th node
    ELEM2NODE(7,:)  = nnod+1:nnod+nel;
    GCOORD          = [GCOORD, [...
        mean(reshape(GCOORD(1, ELEM2NODE(1:3,:)), 3, nel));...
        mean(reshape(GCOORD(2, ELEM2NODE(1:3,:)), 3, nel))]];

    nnod7    = size(GCOORD,2);                                                                     % after remeshing

    ngeo = length(GEO);
  
    if SETTINGS.intersect_s == "layer_collapse"   % give id to undefined phases because of discontinuous layers
        if any(Phases==0)
            nboxes = floor(ngeo/3);
            ELnoPhase = find(Phases==0);
            centralNode = ELEM2NODE(7, ELnoPhase);
                
            for ib=1:nboxes
                intb = max(1, (ib-1)*3);                                              % bottom interface
                intr = ib*3 - 1;                                                      % right    " "
                intt = intr + 1;                                                      % top      " "
                intl = intt + 1;                                                      % left     " "
                BOX_n = [GEO(intb).coo GEO(intr).coo fliplr(GEO(intt).coo) fliplr(GEO(intl).coo)]; % BOX nodes
                InLayer = inpolygon(GCOORD(1,centralNode), GCOORD(2,centralNode),...
                                    BOX_n(1,:), BOX_n(2,:));
                Phases(ELnoPhase(InLayer)) = ib;
            end
        end                                                                  % figure(); plot_phasesF(GCOORD, ELEM2NODE, Phases);
    end % "layer_collapse"                                                   % hold on;  plot_meshF(ELEM2NODE, GCOORD);
                                                                             % hold on; plotGEO(GEOr);
    GEOn = getGEOn(GCOORD, ELEM2NODE, Point_id, Phases);                     % GEO definition including triangle-added vertices

    %[~,Geo_mark] = ismember(round(GEOMETRY'*1e3),round(GCOORD'*1e3),'rows');  % old geometry markers

    % Remesh variables in the ips (accumulated gradient of the deformation,
    % accumulated strain invariant, old stresses, viscosities and strain rates)
    if PHY.SS.rI2 == "triscatteredinterp"
        % Old ips
        [GIPold_x,GIPold_y] = ip_coord(GCO,E2N,size(E2N,2),6);
        % New ips
        [GIPnew_x,GIPnew_y] = ip_coord(GCOORD,ELEM2NODE,size(ELEM2NODE,2),6);
        Fc      = scatteredInterpolant(GIPold_x(:),GIPold_y(:),I2.c(:), ...
            'natural','linear');
        Itsi_c  = Fc(GIPnew_x,GIPnew_y);
        Fp      = scatteredInterpolant(GIPold_x(:),GIPold_y(:),I2.p(:), ...
            'natural','linear');
        Itsi_p  = Fp(GIPnew_x,GIPnew_y);
    end

    nnodel_r = 3;
    extrap_scheme = "ip2nodb2ip_N3";                                           % fine, with some overshooting
    [F_xx, F_xy, F_yx, F_yy, GIP_xF, GIP_yF, I2, TAU_xx_old, ...
        TAU_xy_old, TAU_yy_old, Mu_all, E2all] = remesh_F_TAU_Mu_E2( ...
        GCOORD, ELEM2NODE, GCO, E2N, F_xx, F_xy, F_yx, F_yy, ...
        TAU_xx_old, TAU_xy_old, TAU_yy_old, Mu_all, E2all, I2, extrap_scheme);

    Mu_all(Mu_all<1e18) = 1.e18;
    Mu_all(Mu_all>1e24) = 1.e24;

    if PHY.SS.rI2 == "triscatteredinterp"
        I2.p = Itsi_p;
        I2.c = Itsi_c;
    else       % "linear_N"
        I2.f(I2.f<0) = 0;
        I2.p(I2.p<0) = 0;
        I2.c(I2.c<0) = 0;
    end
    %Mu_all(Mu_all<0) = -Mu_all(Mu_all<0);
    
    % Remesh random
    if PHY.SS.rand_s
        SS = remesh_random_element(GCOORD, ELEM2NODE, Phases, GCO, E2N, Pha, SS, PHY);
    end

    %Convert mesh to 4-times as many linear elements
    EL2NOD4Tri = trimesh_p2_to_p1(E2N(1:6,:), Phases, eids) ;  %this function(trimesh_p2_to_p1) needs only 6 nodes
    nVnod = max(max(EL2NOD4Tri(1:3,:))); % number of vertex nodes

    els = tsearch2(GCO(:,1:nVnod), uint32(EL2NOD4Tri(1:3,:)),GCOORD);
    Ind22 = find(els==0); %check of all elements were found
    if(~isempty(Ind22))
        for i=1:length(Ind22)
            [val, els(Ind22(i))] = min(sqrt((GCO(1,E2N(7,:)) - GCOORD(1,Ind22(i))).^2 + (GCO(2,E2N(7,:)) - GCOORD(2,Ind22(i))).^2));
        end
    end
    if any(els==0)
        error('Remeshing failed in move_contours');
    end

    %Interpolate continuous pressure to new nodes
    LCOORD = local_coords_2d(GCO(:,1:nVnod),EL2NOD4Tri,els,GCOORD);
    PRESS_CONT = interp2d_cubic(GCO,EL2NOD4Tri,els,LCOORD,PRESS_CONT,'nelblk');
    % PRESS_CONT(Ind22)=0; % I think this is wrong (Miguel's note)
    PRESS_CONT = reshape(PRESS_CONT,1,[]);                                     % [1,nnod7]

    %Interpolate Depletion to new nodes
    Dpl = interp2d_cubic(GCO,EL2NOD4Tri,els,LCOORD,Dpl,'nelblk');
    Dpl = reshape(Dpl,1,[]);                                                   % [1,nnod7]
    Dpl(Ind22)=0.1;
    Dpl(unique(ELEM2NODE(:,Phases ~= 1))) = 0.0;                               % nodes out of mantle set to their original values

    if isfield(PHY,'MELT') && PHY.MELT.Ts0 ~= 0
        Ts_dry = PHY.MELT.Ts0 + PHY.MELT.dTs_dP .* PRESS_CONT + PHY.MELT.dTs_dF .* Dpl; % [1,nnod7]
    end

    % Interpolate temperatures
    switch SETTINGS.temp_remesh_type
        case 'shpf_cuadr'
            % old element indexes in which the new node coordinates are contained
            trisT = tsearch2(GCO,uint32(E2N(1:3,:)),GCOORD);
            IndT = find(trisT==0);
            
            % If not all envolving elements found, find closest element
            if(~isempty(IndT))
                for i=1:length(IndT)
                    [~,trisT(IndT(i))] = min(sqrt( ...
                        (GCO(1,E2N(7,:)) - GCOORD(1,IndT(i))).^2 + ...
                        (GCO(2,E2N(7,:)) - GCOORD(2,IndT(i))).^2));
                end
            end
            
            if(any(isnan(trisT)))
                error('remeshing failed in move_contours');
            end
            Temp = remesh_val(trisT,GCO,GCOORD,Temp,E2N);
            Temp(IndT) = PHY.Ttop;
        case 'triscatter'
            FT = scatteredInterpolant(GCO(1,:)',GCO(2,:)',Temp','natural', ...
                'linear');
            Temp = FT(GCOORD(1,:)',GCOORD(2,:)');
    end
    Temp = reshape(Temp,1,[]);

    if isfield(PHY,'MELT') && PHY.MELT.Ts0 ~= 0          % correct temperatures about old solidus
        Temp = min(Temp, Ts_dry);                              % parallel min
    end

    %Interpolate serpentinization
    if all(Dserp == 0.)
        Dserp = zeros(1,nnod7);
    else
        Dserp = interp2d_cubic(GCO,EL2NOD4Tri,els,LCOORD,Dserp,'nelblk');
        Dserp = reshape(Dserp,1,[]);
        Dserp(Dserp < 0.) = 0.;
        to0boo = ~ismember(1:nnod7, unique(ELEM2NODE(:, Phases == 1)));
        Dserp(to0boo) = 0.;                                                 % serpentine only allowed in nodes touching mantle
    end

    %==========================================================================
    % REMESH FACTORS FOR RHEOLOGIC CHANGES INSIDE PHASES
    %==========================================================================
    if size(RHEOLvar,2) > 1
        %     RHEOLvar = ...
        %         remesh_rheol_var(GCOORD,ELEM2NODE,Phases,GCO,E2N,Pha,RHEOLvar);
        RHEOLvar = rheol_var_Dpl_ip(Dpl, ELEM2NODE, Phases, nip, GCOORD);           % Dpl \in [1,nnod7], GCOORD \in [2,nnod7], ELEM2NODE \in [7,nel], Phases \in [1,nel], nip=6 here
    else
        RHEOLvar = {ones(nel,nip)};
    end

end
