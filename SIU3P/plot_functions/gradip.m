function [gradm,Vn] = gradip(Vip,MESH)

nip     = size(Vip,2);
nel     = size(Vip,1);
GCOORD  = MESH.GCOORD;
EL2NOD  = MESH.EL2NOD;
if nip==6
    nnodel = 6;
elseif nip==3
    nnodel = 3;
end

EC_x    = reshape(GCOORD(1,EL2NOD(1:3,:)),3,nel);
EC_z    = reshape(GCOORD(2,EL2NOD(1:3,:)),3,nel);

[IP_X,IP_w] = ip_triangle(nip);
[N, dNds]   = shp_deriv_triangle(IP_X,nnodel);
[N3,dN3ds]  = shp_deriv_triangle(IP_X,3);
nn  = [];
for ip = 1:nip
    nn  = [nn; N{ip}'];
end
N = nn;

Vn      = zeros(nel,nnodel);
gradx   = zeros(nel,nip);
grady   = zeros(nel,nip);
E2N     = zeros(nel,3);
GC      = zeros(2,nel*3);
% Element loop
for el = 1:nel
    E2N(el,:)   = 3*(el-1)+(1:3);
    GC(:,3*(el-1)+(1:3))    = GCOORD(:,EL2NOD(1:3,el));
    Vn(el,:)  = N\Vip(el,:)';
    for ip = 1:nip
        % Jacobi matrix
        J      = dN3ds{ip}'*[EC_x(:,el) EC_z(:,el)];
        % Determinant of J
        detJ   = J(1)*J(4) - J(2)*J(3);
        
        % global derivative of shape function
        dNUdx  = J\dNds{ip}';
        
        gradx(el,ip)    = dNUdx(1,:)*Vn(el,:)';
        grady(el,ip)    = dNUdx(2,:)*Vn(el,:)';
    end
end

gradm = sqrt(gradx.^2+grady.^2);

maxgm = repmat(max(gradm,[],2),1,3)';

gradm3 = gradm(:,1:3)';
Vn3 = Vn(:,1:3)';

h1 = subplot(2,2,1);
patch('faces',E2N,'vertices',GC','facevertexcdata',log10(maxgm(:)),'FaceColor','flat','LineStyle','none')
%patch('faces',E2N,'vertices',GC','facevertexcdata',log10(gradm3(:)),'FaceColor','flat','LineStyle','none')
colorbar
title('Gradients of the viscosity [max(log10(grad(\eta)))]')
xlabel('Distance [m]')
ylabel('Depth [m]')

h2 = subplot(2,2,2);
patch('faces',E2N,'vertices',GC','facevertexcdata',log10(abs(Vn3(:))),'FaceColor','flat','LineStyle','none')
colorbar
shading interp
title('Viscosity [log10(\eta)]')
xlabel('Distance [m]')
ylabel('Depth [m]')

linkaxes([h1 h2])