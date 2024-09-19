function [VAL] = remesh_val_ip(GCO,E2N,VAL_old,nip,nnodel,GCOORD,ELEM2NODE)

% nip          = 6;
% nnodel       = 6;
% Brings back F to the nodes
nel_old = size(E2N,2);
nel = size(ELEM2NODE,2);

%EL2N      = zeros(nel_old,nnodel);

[IP_X, IP_w]    = ip_triangle(nip);
[   Nbig]    = shp_triangle(IP_X, nnodel);

GIP_xF = zeros(nel,nip);
GIP_yF = zeros(nel,nip);

VAL_oldn = zeros(nnodel,nel_old);

for i=1:nel_old
    is         = (i-1)*nnodel+1; ie = (i-1)*nnodel+nnodel;
    %GCO_N(:,is:ie) = GCO(:,E2N([1:nnodel_r],i));
    %EL2N(i,:) = is:ie;
    Dummy      = Nbig'\VAL_old(i,:)';
    VAL_oldn(:,i)= Dummy(1:nnodel)';    
end

% Uncomment the following to plot the F_xx at the nodes
% scatter(GCO(1,E2N(1:nnodel_r,:)),GCO(2,E2N(1:nnodel_r,:)),5,F_xx_n(:))

% Average the values at the nodes???

% Find new ips

[   N, dNdu]    = shp_deriv_triangle(IP_X, size(E2N,1));

ECOORD_x = reshape(GCOORD(1,ELEM2NODE), size(ELEM2NODE,1), nel);
ECOORD_y = reshape(GCOORD(2,ELEM2NODE), size(ELEM2NODE,1), nel);

for ip=1:nip
    
    Ni      =        N{ip};

    GIP_x   = Ni'*ECOORD_x;
    GIP_y   = Ni'*ECOORD_y;
    
    GIP_xF(:,ip)   = GIP_x;
    GIP_yF(:,ip)   = GIP_y;

end

%size_gip = size(GIP_xF);
Tris = tsearch2(GCO,uint32(E2N(1:3,:)),[GIP_xF(:)';GIP_yF(:)']);

Tris = reshape(Tris,size(GIP_xF));
Ind = find(Tris==0);

%check of all elements were found, and if not it solves the problem finding
% the closest element
if(~isempty(Ind))
    for i=1:length(Ind)
        [~,Tris(Ind(i))] = min(sqrt((GCO(1,E2N(7,:)) - GIP_xF(Ind(i))).^2 + (GCO(2,E2N(7,:)) - GIP_yF(Ind(i))).^2));
    end
end

if(any(isnan(Tris)))
    error('remeshing failed in move_contours');
end

% For plotting the new mesh with the new integration points uncomment the
% next lines
%
% figure(265)
% trimesh(ELEM2NODE(1:3,:)',GCOORD(1,:),GCOORD(2,:))
% hold on
% plot(GIP_x_all(:),GIP_y_all(:),'o')

% Calculates the local coordinates of the new ip

xp = GIP_xF(:)'; % New ip
yp = GIP_yF(:)'; % New ip

x = reshape(GCO(1,E2N(1:3,Tris)),3,size(xp,2)); % Reshape old coord
y = reshape(GCO(2,E2N(1:3,Tris)),3,size(yp,2)); % Reshape old coord

xi = -(-x(1,:).*yp+x(1,:).*y(3,:)-x(3,:).*y(1,:)+xp.*y(1,:)+x(3,:) ...
    .*yp-xp.*y(3,:))./(-x(1,:).*y(3,:)+x(1,:).*y(2,:)-x(2,:).*y(1,:) ...
    +x(2,:).*y(3,:)+x(3,:).*y(1,:)-x(3,:).*y(2,:)); % Local x cood for ip
yi = (x(1,:).*y(2,:)-x(1,:).*yp-x(2,:).*y(1,:)+x(2,:).*yp+xp ...
    .*y(1,:)-xp.*y(2,:))./(-x(1,:).*y(3,:)+x(1,:).*y(2,:)-x(2,:).*y(1,:)...
    +x(2,:).*y(3,:)+x(3,:).*y(1,:)-x(3,:).*y(2,:)); % Local y coord for ip

% Load shape functions

[NN] = shp_triangle([xi' yi'], nip);

% Calculate the value of F for the new ips

VAL = sum(NN.*VAL_oldn(1:nnodel,Tris));
VAL = reshape(VAL,size(Tris,1),nip);