function plot_val(MESH,Val,nel,nip,alphac,Phases)

GCOORD = MESH.GCOORD;
ELEM2NODE = MESH.EL2NOD;

if isfield(MESH,'sign')
    if MESH.sign==-1
        GCOORD(1,:) = -GCOORD(1,:);
        sc = -1;
    else
        sc = 1;
    end
else
    sc = 1;
end

EL2N      = zeros(nel,3); % new connectivity matrix
GCOORD_N  = zeros(2,nel*3);
VV_n      = zeros(3,nel);
EL2N(1,:) = 1:3;
[IP_X, ~]    = ip_triangle_m2tri(nip);
[   Nbig,~]    = sf_dsf_tri367_N(IP_X, 3,'matrix');
% [IP_X, ~]    = ip_triangle(nip);
% [   Nbig,~]    = shp_deriv_triangle(IP_X, 3);
% 
% U = zeros(3,length(Nbig));
% for n = 1:length(Nbig)
%     U(:,n) = Nbig{n};
% end
% Nbig = U;

for i=1:nel
    is         = (i-1)*3+1; ie = (i-1)*3+3;
    GCOORD_N(:,is:ie) = GCOORD(:,ELEM2NODE([1 2 3],i));
    EL2N(i,:) = is:ie;
    Dummy      = Nbig'\Val(i,:)';
    VV_n(:,i)= Dummy(1:3);
end

GCOORD_N(1,:) = GCOORD_N(1,:)*sc;

if nargin==6
    % Phases
    Pha = unique(Phases);
    nphases = length(Pha);
    % Loop through phases
    for n = 1:nphases
        % Averaging at nodes
        el2n = ELEM2NODE(1:3,Phases==Pha(n));
        [~,~,rep_e2n] = unique(el2n);
        EL2Nt = 1:length(rep_e2n);
        EL2Nt = EL2Nt(rep_e2n);
        VV_np = VV_n(:,Phases==Pha(n));
        Vector_smooth = accumarray(EL2Nt(:),VV_np(:),[],@mean);
        % Find vector to repeat nodes
        [~,~,rep_order] = unique(EL2Nt);
        % Organize in matrix form to plot
        VV_n(:,Phases==Pha(n)) = ...
            reshape(Vector_smooth(rep_order),size(VV_np,1),size(VV_np,2));
    end
end

if nargin<5
    patch('faces',EL2N,'vertices',GCOORD_N','facevertexcdata',VV_n(:), ...
        'FaceColor','flat')
    shading interp
    return
elseif isempty(alphac)
    patch('faces',EL2N,'vertices',GCOORD_N','facevertexcdata',VV_n(:), ...
        'FaceColor','flat')
    shading interp
else
    patch('faces',EL2N,'vertices',GCOORD_N','FaceVertexAlphaData', ...
        VV_n(:),'AlphaDataMapping','scaled','FaceAlpha','interp',...
        'FaceColor',alphac,'EdgeColor','none')
end
end