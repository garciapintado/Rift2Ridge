function [AX_FS_SPARSE,AY_FS_SPARSE,Dx_e] = fssa(MESH,RHO,G,dt,nip, ...
    fs_alpha,fs_beta, int_id, Load_el, rho_ext)
    % [AX_FS_SPARSE,AY_FS_SPARSE,DX_E] = FSSA(MESH,RHO,G,DT,NIP,ALPHA,BETA,
    % INT_ID,LOAD_EL,RHO_EXT) calculates AX_FS_SPARSE and AY_FS_SPARSE
    % penalisation terms for the free surface stabilisation algorithm described
    % in Andres-Martinez, 2015. The stabilisation algorithm is applied to the
    % interface INT_ID, with BETA (for x) and ALPHA (for y) penalisation
    % factors. MESH is a structure containing information about the FEM mesh.
    % RHO is a density matrix defined in the ips. DT is the time step, NIP is
    % number of integration points per 2D element. In case an external load is
    % applied to the surface and this load corresponds to an external body with
    % a given density (i.e. water load), LOAD_EL indicates which elements of
    % this surface are loaded and RHO_EXT is the density of the external body.
    % Note that the fssa doesn't include the force of an external body, so that
    % this parameters are only used for adjusting the penalisation terms to the
    % density difference. DX_E are the elements loaded.
    %
    %
    %                 Top_nodes x
    %     x-----x-----x-----x-----x-----x-----x-----x-----x Surface
    %                  \Top_elem3/ \Top_elem3/
    %                   \       /   \       /
    %                    o     o     o     o
    %                     \   /       \   /
    %                      \ /         \ /
    %                       o           o
    %

    %--------------------------------------------------------------------------
    % Function written by Miguel Andres-Martinez, Postdoc at University of
    % Bremen, 20-10-2017. Email: andresma@uni-bremen.de
    % Algorithm described in Andres-Martinez, 2015.
    %--------------------------------------------------------------------------

    %==========================================================================
    % DECLARE VARIABLES
    %==========================================================================
    Point_id    = MESH.Point_id;
    ELEM2NODE   = MESH.EL2NOD;
    GCOORD      = MESH.GCOORD;
    sdof        = 2*MESH.nnod;

    if isempty(Load_el)
        Load_el = zeros(1,size(ELEM2NODE,2))==1;                               % [1,nel]
    end

    %==========================================================================
    % CALCULATE FREE SURFACE STABILISATION TERMS
    %==========================================================================
    nip_fs = 3;     % Number of integration points.
    nnod_el_fs = 3; % Number of superficial nodes in the element.
    Top_nodes = find(Point_id==int_id);
    % Indexes of the top nodes.
    Top_elem3 = find(sum(ismember(ELEM2NODE,Top_nodes),1)==...
        nnod_el_fs);
    topnel = numel(Top_elem3);

    % Top element indexes with 3 nodes on the surface.
    % Density of the external material vector with value 0 where
    % material is not in place, and rho_ext where the element is
    % loaded (i.e. water load)
    if length(rho_ext) == 1
        rho_ext = repelem(rho_ext,1,topnel);
    end
    if ~isequal(size(rho_ext),size(Top_elem3))
        error("fssa:: length(rho_ext) does not match length(topel)")
    end

    if any(Load_el)
        if ~isequal(find(Load_el),Top_elem3)
             error("fssa:: input Load_el does not match top element identification")
        end
    end

    Rho_ext = Load_el(Top_elem3)' .* rho_ext(:);             % [topnel,1]

    inc_rho = sum(RHO(Top_elem3,:),2)/nip - Rho_ext;
    % Increment of the density.
    Ip_fs = [-sqrt(3/5), 0, sqrt(3/5)]; % Integration points
    Ipw_fs = [5/9; 8/9; 5/9]; % Integration weight
    [N_fs] = shp_line_int(Ip_fs, nnod_el_fs); % Shape functions

    AX_I_FS = zeros(nnod_el_fs*nnod_el_fs,size(Top_elem3,2));     % Matrix of indexes i for free surface elements in the x component.
    AX_J_FS = AX_I_FS;                                            % Matrix of indexes j for free surface elements in the x component. i.e. (3*3, elements on top)
    AX_FS   = AX_I_FS;                                            % Values for correction of matrix A for the free surface in the x component.
    AY_I_FS = AX_I_FS;                                            % Matrix of indexes i for free surface elements in the y component.
    AY_J_FS = AX_I_FS;                                            % Matrix of indexes j for free surface elements in the y component.
    AY_FS   = AX_I_FS;                                            % Values for correction of matrix A for the free surface in the y component.

    Dx_e = zeros(size(Top_elem3,2),1); % dx of surface elements for the Courant time step criterion.

    F_fs = zeros(sdof, 1);             % Vector of the values for the correction of the Rhs for the free surface.

    for i = 1:size(Top_elem3,2)                                      % Loop for each superficial element
        Lo_nodes_fs = Point_id(ELEM2NODE(1:6,Top_elem3(i)))==int_id; % Local  indexes of superficial nodes.
        Gl_nodes_fs = ELEM2NODE(Lo_nodes_fs,Top_elem3(i));           % Global indexes of superficial nodes.
        
        [Xorder, Xindex] = sort(GCOORD(1,Gl_nodes_fs)); % order surface elements from left to right
        Gl_nodes_fs = Gl_nodes_fs(Xindex);              % reorder global node indexes
        Dx_e(i) = (max(Xorder)-min(Xorder));            % calculate dx.
        h_fs = GCOORD(2,Gl_nodes_fs);                   % Heights of the nodes.
        dhdx = (h_fs(end)-h_fs(1))/Dx_e(i);             % Slope at each node.
        
        Ax_ie = (2*double(Gl_nodes_fs))*ones(1,nnod_el_fs);    % indexes i for the matrix Ax_fse of the evaluated element for the x component.
        Ay_ie = (2*double(Gl_nodes_fs))*ones(1,nnod_el_fs);    % indexes i for the matrix Ay_fse of the evaluated element for the y component.
        Ax_ie = Ax_ie(:);
        Ay_ie = Ay_ie(:);
        Ax_je = ones(nnod_el_fs,1)*(2*double(Gl_nodes_fs)-1)'; % indexes j for the matrix Ax_fse of the evaluated element for the x component.
        Ay_je = ones(nnod_el_fs,1)*2*double(Gl_nodes_fs)';     % indexes j for the matrix Ay_fse of the evaluated element for the y component.
        Ax_je = Ax_je(:);
        Ay_je = Ay_je(:);
        F_index = 2*double(Gl_nodes_fs);                       % indexes for the vector F
        
        % Initialize the matrices for the correction of A and F for each element.
        Ax_fse = zeros(nnod_el_fs,nnod_el_fs);
        Ay_fse = zeros(nnod_el_fs,nnod_el_fs);
        F_fse =  zeros(nnod_el_fs,1);
        
        for ip_fs = 1:nip_fs             % Integration loop
            Ax_fse = Ax_fse + N_fs{ip_fs}*N_fs{ip_fs}'*...
                Ipw_fs(ip_fs)*Dx_e(i)/2; % Ni*Nj*w*dx
            Ay_fse = Ay_fse + N_fs{ip_fs}*N_fs{ip_fs}'*...
                Ipw_fs(ip_fs)*Dx_e(i)/2; % Ni*Nj*w*dx
            F_fse = F_fse + N_fs{ip_fs}*Ipw_fs(ip_fs)*Dx_e(i)/2;
            % Ni*w*dx
        end
        
        Ax_fse_col = Ax_fse(:)*inc_rho(i)*G(2)*fs_alpha*fs_beta * -(dhdx)*dt; % Afs for x = inc_density*g* fs_beta*dhdx*dt*integral(Ni*Nj*dx)
        Ay_fse_col = Ay_fse(:)*inc_rho(i)*G(2)*fs_alpha *dt;                  % Afs for y = inc_density*g* fs_alpha    *dt*integral(Ni*Nj*dx)
        
        % Fills the global matrices for the free surface correction
        AY_I_FS(:,i) = Ay_ie;
        AY_J_FS(:,i) = Ay_je;
        AY_FS(:,i) = Ay_fse_col;
        
        AX_I_FS(:,i) = Ax_ie;
        AX_J_FS(:,i) = Ax_je;
        AX_FS(:,i) = Ax_fse_col;
        
        F_fs(F_index) = 0; %F_fse(:).*h_fs'*inc_rho*G(2);
        % Ffs = inc_density*g*h*integral(Ni*dx)
        
    end % for

    % Delete the non-diagonal terms of the local matrices
    % matrices for the 3 surficial nodes of each surficial element
    %
    %  / k1  k4  k7 \   Remove           / k1         \
    % |  k2  k5  k8  |  non-diagonal => |  k2  k5      |
    %  \ k3  k6  k9 /   terms            \ k3  k6  k9 /

    V_rm = AY_J_FS(:) > AY_I_FS(:); % Vector for removing the elements of the matrix A_FS which j>i.
    AY_FS(V_rm)   = [];
    AY_I_FS(V_rm) = [];
    AY_J_FS(V_rm) = [];

    AX_FS_SPARSE = sparse2(AX_I_FS, AX_J_FS, AX_FS, sdof, sdof); % sparsification of AX_FS.
    AY_FS_SPARSE = sparse2(AY_I_FS, AY_J_FS, AY_FS, sdof, sdof); % sparsification of AY_FS
end % function