function Pc = el2nod_pressure(MESH, Pd, MODEL, method)

    % WRITTEN BY JOERG. ADAPTED BY ELENA
    % Unified by JGP
    
    % MESH.GCOORD = GCOORD;
    % MESH.EL2NOD = ELEM2NODE;
    % MESH.nel    = nel;

    if nargin < 3
        MODEL = "kinedyn";
    end
    switch MODEL
        case "kinedyn"
            eid = [4 5 6];
        case "rift2ridge2D" 
            eid = [6 4 5];
        otherwise
        error('el2nod_pressure:: wrong model specification')
    end

    if nargin < 4
        method = 'std';  % e.g. default for kinedyn
    end
    
    FigNo  = 0;

    Pd= Pd(:);

    if length(Pd)~=3*MESH.nel
        error(' Pressure vector must have length of 3 x nel');
    end

    switch method
        case 'std'
            [Pc,MMc,MMd,D2C] = standard_method(MESH,Pd);
            
        case 'opt'
            [Pc,MMc,MMd,D2C] = optimized_method(MESH,Pd);
    %         [Pc2,MMc2,MMd2,D2C2] = standard_method(MESH,Pd);
    end

    Pc(MESH.EL2NOD(eid(1),:)) = 0.5 .* (Pc(MESH.EL2NOD(1,:)) + Pc(MESH.EL2NOD(2,:)));
    Pc(MESH.EL2NOD(eid(2),:)) = 0.5 .* (Pc(MESH.EL2NOD(2,:)) + Pc(MESH.EL2NOD(3,:)));
    Pc(MESH.EL2NOD(eid(3),:)) = 0.5 .* (Pc(MESH.EL2NOD(1,:)) + Pc(MESH.EL2NOD(3,:)));

    Pc(MESH.EL2NOD(7,:)) = (1/3) .* ...
        (Pc(MESH.EL2NOD(1,:)) + Pc(MESH.EL2NOD(2,:)) + Pc(MESH.EL2NOD(3,:)));

    if FigNo
        EL2NODP = int32(reshape(1:(3*MESH.nel),3,[]));
        meshcol = 'none';
        visible = 1;
        plot_2d_pressure(FigNo,MESH.GCOORD,MESH.EL2NOD,EL2NODP,Pd,[],[],meshcol,visible);
    %     plot_2d_tridata(FigNo+1,MESH.GCOORD,MESH.EL2NOD,Pc,[],[],meshcol,visible);
    %     plot_2d_tridata(FigNo+2,MESH.GCOORD,MESH.EL2NOD,Pc2,[],[],meshcol,visible);
    end

end % END OF FUNCTION el2nod_pressure

% #########################################################################
%                               SUBFUNCTIONS
% #########################################################################

function [Pc,MMc,MMd,D2C] = optimized_method(MESH,Pd)

    nip     = 3;
    nPdofel = 3;
    nvert   = 3;            % number of vertex nodes
    nelblk  = 1200;
    nel     = MESH.nel;
    GCOORD  = MESH.GCOORD;
    EL2NOD  = MESH.EL2NOD(1:3,:);
    EL2P    = reshape(1:nPdofel*nel,nPdofel,nel);

    [IP_X,IP_w] = ip_triangle_m2tri(nip);
        % local coordinates of points where stresses are calculated
    [N,dN]       = sf_dsf_tri367(IP_X,nvert,'cell');
    %[N,dN]       = shp_deriv_triangle(IP_X,nvert);
        % shape function values at integration points
        
    M_all   = zeros(nPdofel*(nPdofel+1)/2,nel); % storage for global mass matrix
    nelblk  = min(nelblk,nel);
    nblk    = ceil(nel/nelblk);
    il      = 1;
    iu      = nelblk;
    for ib=1:nblk % loop over element blocks
        %==============================================================
        % CALCULATE JACOBIAN, ITS DETERMINANT AND INVERSE
        %==============================================================
        % NOTE: For triangular elements with non-curved edges the Jacobian is
        %       the same for all integration points (i.e. calculated once
        %       before the integration loop). Further, linear 3-node shape are
        %       sufficient to calculate the Jacobian.
        ECOORD_x   = reshape( GCOORD(1,EL2NOD(1:nvert,il:iu)), nvert, nelblk );
        ECOORD_z   = reshape( GCOORD(2,EL2NOD(1:nvert,il:iu)), nvert, nelblk );
        %     Jx         = ECOORD_x'*dN{1}';
        %     Jz         = ECOORD_z'*dN{1}';
        Jx         = ECOORD_x'*dN{1}';
        Jz         = ECOORD_z'*dN{1}';

        detJ       = Jx(:,1).*Jz(:,2) - Jx(:,2).*Jz(:,1);
        
        M_blk      = zeros(nelblk,nPdofel*(nPdofel+1)/2);
        for ip=1:nip
            % ASSEMBLE RHS VECTOR FOR LSQ FIT
            
            % MASS MATRICES FOR ALL ELEMENTS IN BLOCK
            weight = detJ.*IP_w(ip);
            indx   = 1;
            for i=1:nPdofel
                for j=i:nPdofel
                    M_blk(:,indx) = M_blk(:,indx) + weight .* N{ip}(i)*N{ip}(j);
                    indx = indx + 1;
                end
            end
        end % END OF INTEGRATION LOOP
        
        M_all  (:,il:iu) = M_blk';
        
        %==============================================================
        % READJUST START, END AND SIZE OF BLOCK. REALLOCATE MEMORY
        %==============================================================
        il  = il + nelblk;
        if(ib==nblk-1)
            nelblk = nel-iu;
        end
        iu  = iu + nelblk;
    end % END OF BLOCK LOOP

    % FORM GLOBAL MASS MATRIX
    indx_j     = repmat(1:nPdofel,nPdofel,1);
    indx_i     = indx_j';
    indx_i     = tril(indx_i); indx_i = indx_i(:); indx_i = indx_i(indx_i>0);
    indx_j     = tril(indx_j); indx_j = indx_j(:); indx_j = indx_j(indx_j>0);

    Mc_i       = EL2NOD(indx_i,:); Mc_i = Mc_i(:);
    Mc_j       = EL2NOD(indx_j,:); Mc_j = Mc_j(:);
    indx       = Mc_i < Mc_j;
    tmp        = Mc_j(indx);
    Mc_j(indx) = Mc_i(indx);
    Mc_i(indx) = tmp; clear tmp indx
    %MMc        = sparse3(Mc_i, Mc_j, M_all(:));
    MMc        = sparse2(Mc_i, Mc_j, M_all(:));
    MMc        = MMc + tril(MMc,-1)';

    Md_i       = EL2P(indx_i,:); Md_i = Md_i(:);
    Md_j       = EL2P(indx_j,:); Md_j = Md_j(:);
    indx       = Md_i < Md_j;
    tmp        = Md_j(indx);
    Md_j(indx) = Md_i(indx);
    Md_i(indx) = tmp; clear tmp indx
    %MMd        = sparse3(Md_i, Md_j, M_all(:));
    MMd        = sparse2(Md_i, Md_j, M_all(:));
    MMd        = MMd + tril(MMd,-1)';

    % pointer continuous --> discontinuous
    indx_j     = repmat(1:nPdofel,nPdofel,1);
    indx_i     = indx_j';
    Mc_i       = EL2NOD(indx_i,:); Mc_i = Mc_i(:);
    Md_i       = EL2P  (indx_i,:); Md_i = Md_i(:);
    %D2C        = sparse3(Mc_i(:),Md_i(:),ones(size(Md_i(:))));
    D2C        = sparse2(Mc_i(:),Md_i(:),ones(size(Md_i(:))));
    D2C(D2C>0) = 1; % MMc = D2C*MMd*D2C' (I checked that!!!!)
    %Pc         = MMc \ ((D2C'*MMc)'*Pd);
    Pc         = MMc \ (D2C*MMd*Pd);

end % function

% #########################################################################

function [Pc,MMc,MMd,D2C] = standard_method(MESH,Pd)

    nip     = 3;
    nPdofel = 3;
    nvert   = 3;            % number of vertex nodes
    nel     = MESH.nel;
    GCOORD  = MESH.GCOORD;
    EL2NOD  = MESH.EL2NOD(1:3,:);
    EL2P    = reshape(1:nPdofel*nel,nPdofel,nel);

    [IP_X,IP_w] = ip_triangle_mtri(nip);
        % local coordinates of points where stresses are calculated
    [N,dN]       = sf_dsf_tri367(IP_X,nvert,'cell');
        % shape function values at integration points

    % Global matrices
    rhs  = zeros(nPdofel*nel,1);
    MMDv = zeros(nPdofel*nPdofel,nel); % storage for mass matrix data
    MMDi = zeros(nPdofel*nPdofel,nel); % storage for mass matrix i indices
    MMDj = zeros(nPdofel*nPdofel,nel); % storage for mass matrix j indices
    MMCv = zeros(nPdofel*nPdofel,nel); % storage for mass matrix data
    MMCi = zeros(nPdofel*nPdofel,nel); % storage for mass matrix i indices
    MMCj = zeros(nPdofel*nPdofel,nel); % storage for mass matrix j indices

    for iel = 1:nel % main element loop
        elnod  = double(EL2NOD(1:nvert,iel));
        elnodP = double(EL2P(1:nvert,iel));    
        
        %==============================================================
        % CALCULATE JACOBIAN, ITS DETERMINANT AND INVERSE
        %==============================================================
        % NOTE: For triangular elements with non-curved edges the Jacobian is
        %       the same for all integration points so that it can be 
        %       calculated once before the integration loop. Linear 3-node 
        %       shape are sufficient to calculate the Jacobian.
        ECOORD = GCOORD(:,elnod);       % element node vertex coordinates
        J      = dN{1}*ECOORD';         % Jacobian: reference element ==> current element
        detJ   = J(1)*J(4) - J(2)*J(3); % Determinat of Jacobian
        % Check element deformation
        if(detJ<=0)
            err = ['negative jacobian in element ' iel];
            error(err);
        end

        % Initialize element arrays
        MMe    = zeros(nPdofel);   % element mass matrix
        Pe     = zeros(nPdofel,1);
        for ip=1:nip
            weight = detJ*IP_w(ip);
            Pe     = Pe  + weight * N{ip} * (N{ip}'*Pd(elnodP)); % 
            MMe    = MMe + weight * (N{ip}*N{ip}');  % element mass matrix
        end
        
        % Assemble discontinuous global mass matrices / vectors
        rhs(elnodP) = rhs(elnodP) + Pe;
        rows        = elnodP*ones(1,nPdofel);
        cols        = ones(nPdofel,1)*elnodP';
        MMDv(:,iel)  = MMe(:);  % values
        MMDi(:,iel)  = rows(:); % i indices
        MMDj(:,iel)  = cols(:); % j indices
        
        % Assemble continuous global mass matrices
        rows        = elnod*ones(1,nPdofel);
        cols        = ones(nPdofel,1)*elnod';
        MMCv(:,iel)  = MMe(:);  % values
        MMCi(:,iel)  = rows(:); % i indices
        MMCj(:,iel)  = cols(:); % j indices
    end
    MMc = sparse2(MMCi(:),MMCj(:),MMCv(:)); % continuous global mass
    MMd = sparse2(MMDi(:),MMDj(:),MMDv(:)); % discontinuous global mass

    % pointer continuous --> discontinuous
    % d2c            = zeros(nPnod,1);
    % d2c(nodesP(:)) = reshape(nodes(:,1:3),[],1);
    D2C        = sparse2(MMCi(:),MMDi(:),ones(size(MMDi(:))));
    D2C(D2C>0) = 1; % MMc = D2C*MMd*D2C' (I checked that!!!!)
    %Pc         = MMc \ ((D2C'*MMc)'*Pd);
    Pc         = MMc \ (D2C*MMd*Pd);

    % CC        = DC' * MMc * DC; % continous mass matrix is CC;
    % Rhs       = DC' * MMc * P;% rhs for pressure smoothing ('force' vector)
    % Pcont     = CC\Rhs; % Solve for LSQ continuous Pressure, could also use lumped mass iterative smoothing

end % function
