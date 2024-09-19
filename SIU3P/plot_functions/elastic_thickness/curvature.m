function [MESH] = curvature(MESH)

inod_surface = find(ismember(MESH.PointID,MESH.PointID_top));

x_surface=MESH.GCOORD(1,inod_surface);
z_surface=MESH.GCOORD(2,inod_surface);

[x_surface,order]=sort(x_surface);
z_surface=z_surface(order);

[x_surface,p]=unique(x_surface);
z_surface=z_surface(p);

x_surface=x_surface.*1000;
z_surface=z_surface.*1000;
der1=diff(z_surface)./diff(x_surface);
der1=[der1(1) der1];
der2=diff(der1)./diff(x_surface);
der2=[der2 der2(end)];

if isfield(MESH,'surf_curv') 
    clearvars MESH.surf_curv;
end

MESH.surf_curv=der2;


%debug
figure(666); clf
plot(x_surface,MESH.surf_curv)


end

