function [VALnod,GCO,E2N] = ip2node_smooth(VALip,GCOORD,ELEM2NODE,nip)

FigNo  = 0;
method = 'std';

% if length(VALip)~=3*nel
%     error(' Pressure vector must have length of 3 x nel');
% end

MESH.EL2NOD = ELEM2NODE;
MESH.GCOORD = GCOORD;
MESH.nel = size(ELEM2NODE,2);
MESH.nip = nip;

switch method
    case 'std'
        [VALnod,MMc,Rhs,GCO,E2N] = standard_method(MESH,VALip);
        
    case 'opt'
        [VALnod,MMc,Rhs,GCO,E2N] = optimized_method(MESH,VALip);
%         [Pc2,MMc2,MMd2,D2C2] = standard_method(MESH,VALip);
end

if FigNo
%     patch('faces',E2N','vertices',GCO'/1000,'facevertexcdata', ...
%         log10(VALnod(:)),'FaceColor','flat')
%     shading flat
    X = [100000 110000];
    Y = [-100000 1000];
    
    ind = GCO(1,:)>X(1) & GCO(1,:)<X(2) & GCO(2,:)>Y(1) & GCO(2,:)<Y(2);
    scatter(GCO(1,ind)/1000,GCO(2,ind)/1000,70,'w^','MarkerFace','flat')
    hold on
    scatter(GCO(1,ind)/1000,GCO(2,ind)/1000,50,log10(VALnod(ind)'),'^','MarkerFace','flat')
    hold on
    
    GIP_x_all = zeros(size(ELEM2NODE,2),6);
    GIP_y_all = zeros(size(ELEM2NODE,2),6);
    
    [IP_X, IP_w]    = ip_triangle(nip);
    [   N, dNdu]    = shp_deriv_triangle(IP_X, 6);
    ECOORD_x    = reshape( GCOORD(1,ELEM2NODE(1:6,:)), 6, size(ELEM2NODE,2));
    ECOORD_y    = reshape( GCOORD(2,ELEM2NODE(1:6,:)), 6, size(ELEM2NODE,2));
    for ip = 1:nip
        Ni = N{ip};
        GIP_x_all(:,ip) = Ni'*ECOORD_x;
        GIP_y_all(:,ip) = Ni'*ECOORD_y;
    end
    indx = GIP_x_all(:,6)>X(1) & GIP_x_all(:,6)<X(2) & ...
        GIP_y_all(:,6)>Y(1) & GIP_y_all(:,6)<Y(2);
    GX = GIP_x_all(indx,:)';
    GY = GIP_y_all(indx,:)';
    VALp = VALip(indx,:);
    scatter(GX(:)/1000,GY(:)/1000,70,'w','MarkerFace','flat')
    scatter(GX(:)/1000,GY(:)/1000,50,log10(VALp(:)),'MarkerFace','flat')
end

end % END OF FUNCTION ip2nod_smooth

% #########################################################################
%                               SUBFUNCTIONS
% #########################################################################

function [VALnod,MMc,Rhs,GCO,E2N] = optimized_method(MESH,VALip)

nip     = MESH.nip;
nelblk  = 1200;
ndofel = 3;
nvert   = 3;            % number of vertex nodes
nel     = size(MESH.EL2NOD,2);
GCOORD  = MESH.GCOORD;
EL2NOD  = MESH.EL2NOD;
EL2P    = reshape(1:ndofel*nel,ndofel,nel);
% Resize connectivity matrix and coordinates
GCO     = GCOORD(:,ismember(1:size(GCOORD,2),EL2NOD(1:3,:)));
EL2NOD3 = EL2NOD(1:nvert,:);
[E2Nsort,ind] = sort(EL2NOD3(:));
uni = unique(E2Nsort);
[~,newind] = ismember(E2Nsort,uni);
sortback = 1:length(E2Nsort);
sortback(ind) = sortback;
E2N = reshape(newind(sortback),size(EL2NOD3,1),size(EL2NOD3,2));

[IP_X,IP_w] = ip_triangle(nip);
    % local coordinates of points where stresses are calculated
% [N,dN]       = sf_dsf_tri367(IP_X,nvert);
[N,dN]    = shp_deriv_triangle(IP_X,nvert);
    % shape function values at integration points   
    
M_all   = zeros(ndofel*(ndofel+1)/2,nel); % storage for global mass matrix
Rhs_all = zeros(nel,ndofel);
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
    ECOORD_x   = reshape( GCO(1,E2N(1:nvert,il:iu)), nvert, nelblk );
    ECOORD_z   = reshape( GCO(2,E2N(1:nvert,il:iu)), nvert, nelblk );
    Jx         = (dN{1}'*ECOORD_x)';
    Jz         = (dN{1}'*ECOORD_z)';
    detJ       = Jx(:,1).*Jz(:,2) - Jx(:,2).*Jz(:,1);
    
    M_blk      = zeros(nelblk,ndofel*(ndofel+1)/2);
    F_blk      = zeros(nelblk,ndofel);
    
    for ip=1:nip
        % MASS MATRICES FOR ALL ELEMENTS IN BLOCK
        weight = detJ.*IP_w(ip);
        indx   = 1;
        
        % ASSEMBLE RHS VECTOR FOR LSQ FIT
        F_blk(:,1:ndofel) = F_blk(:,1:ndofel) + (VALip(il:iu,ip).*weight)*N{ip}';
        
        for i=1:ndofel
            for j=i:ndofel
                M_blk(:,indx) = M_blk(:,indx) + weight .* N{ip}(i)*N{ip}(j);
                indx = indx + 1;
            end
        end
    end % END OF INTEGRATION LOOP
    
    Rhs_all(il:iu,:) = Rhs_all(il:iu,:) + F_blk;
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
indx_j     = repmat(1:ndofel,ndofel,1);
indx_i     = indx_j';
indx_i     = tril(indx_i); indx_i = indx_i(:); indx_i = indx_i(indx_i>0);
indx_j     = tril(indx_j); indx_j = indx_j(:); indx_j = indx_j(indx_j>0);

Mc_i       = E2N(indx_i,:); Mc_i = Mc_i(:);
Mc_j       = E2N(indx_j,:); Mc_j = Mc_j(:);
indx       = Mc_i < Mc_j;
tmp        = Mc_j(indx);
Mc_j(indx) = Mc_i(indx);
Mc_i(indx) = tmp; clear tmp indx
MMc        = sparse2(Mc_i, Mc_j, M_all(:));
MMc        = MMc + tril(MMc,-1)';
E2N_DOF    = E2N';
Rhs  = accumarray(E2N_DOF(:), Rhs_all(:));

VALnod = MMc\Rhs;

end

% #########################################################################

function [VALnod,MMc,Rhs,GCO,E2N] = standard_method(MESH,VALip)

nip     = MESH.nip;
ndofel = 3;
nvert   = 3;            % number of vertex nodes
nel     = size(MESH.EL2NOD,2);
GCOORD  = MESH.GCOORD;
EL2NOD  = MESH.EL2NOD;
EL2P    = reshape(1:ndofel*nel,ndofel,nel);
% Resize connectivity matrix and coordinates
GCO     = GCOORD(:,ismember(1:size(GCOORD,2),EL2NOD(1:3,:)));
EL2NOD3 = EL2NOD(1:nvert,:);
[E2Nsort,ind] = sort(EL2NOD3(:));
uni = unique(E2Nsort);
[~,newind] = ismember(E2Nsort,uni);
sortback = 1:length(E2Nsort);
sortback(ind) = sortback;
E2N = reshape(newind(sortback),size(EL2NOD3,1),size(EL2NOD3,2));

[IP_X,IP_w] = ip_triangle(nip);
    % local coordinates of points where stresses are calculated
% [N,dN]       = sf_dsf_tri367(IP_X,nvert);
[N,dN]    = shp_deriv_triangle(IP_X,nvert);
    % shape function values at integration points

% Global matrices
Rhs  = zeros(size(GCO,2),1);
MMCv = zeros(ndofel*ndofel,nel); % storage for mass matrix data
MMCi = zeros(ndofel*ndofel,nel); % storage for mass matrix i indices
MMCj = zeros(ndofel*ndofel,nel); % storage for mass matrix j indices

for iel = 1:nel % main element loop
    elnod  = double(E2N(1:nvert,iel));   
    
    %==============================================================
    % CALCULATE JACOBIAN, ITS DETERMINANT AND INVERSE
    %==============================================================
    % NOTE: For triangular elements with non-curved edges the Jacobian is
    %       the same for all integration points so that it can be 
    %       calculated once before the integration loop. Linear 3-node 
    %       shape are sufficient to calculate the Jacobian.
    ECOORD = GCO(:,elnod); % element node vertex coordinates
    J      = ECOORD*dN{1};   % Jacobi matrix: reference element ==> current element
    detJ   = J(1)*J(4) - J(2)*J(3); % Determinate of Jacobi matrix
    % Check element deformation
    if(detJ<=0)
        err = ['negative jacobian in element ' iel];
        error(err);
    end

    % Initialize element arrays
    MMe    = zeros(ndofel);   % element mass matrix
    Fe     = zeros(ndofel,1);
    for ip=1:nip
        weight = detJ*IP_w(ip);
        Fe     = Fe  + weight * N{ip} * VALip(iel,ip); % 
        MMe    = MMe + weight * (N{ip}*N{ip}');  % element mass matrix
    end
    
    % Assemble continuous global mass matrices
    Rhs(elnod)  = Rhs(elnod) + Fe;
    rows        = elnod*ones(1,ndofel);
    cols        = ones(ndofel,1)*elnod';
    MMCv(:,iel)  = MMe(:);  % values
    MMCi(:,iel)  = rows(:); % i indices
    MMCj(:,iel)  = cols(:); % j indices
end
MMc = sparse2(MMCi(:),MMCj(:),MMCv(:)); % continuous global mass

VALnod = MMc\Rhs;

end
