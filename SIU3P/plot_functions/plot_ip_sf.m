% PLOT_IP_SF plots linear shape functions and integration points

%% Input
% Number of integration points
nip = 6;
% Number of nodes and shape functions
nnodel = 3;


%% Ips and shape functions
[IP_X, IP_w]    = ip_triangle(nip);
[   Nbig]    = shp_triangle(IP_X, nnodel);

%% Plot
clf
plot(IP_X(:,1),IP_X(:,2),'.')
hold on
axis equal
axis([0 1 0 1])
plot([0 1 0 0],[0 0 1 0])
text(IP_X(:,1)+0.01,IP_X(:,2)+0.01,{{1},{2},{3},{4},{5},{6}})
g1 = patch('Faces',[1 2 3],'Vertices',[0 0 1; 0 1 0; 1 0 0],'FaceColor','green');
g2 = patch('Faces',[1 2 3],'Vertices',[0 0 0; 0 1 1; 1 0 0],'FaceColor','blue');
g3 = patch('Faces',[1 2 3],'Vertices',[0 0 0; 0 1 0; 1 0 1],'FaceColor','red');
g1.FaceAlpha = 0.2;
g2.FaceAlpha = 0.2;
g3.FaceAlpha = 0.2;
plot3(IP_X(:,1),IP_X(:,2),Nbig(1,:),'og','MarkerFaceColor','green')
plot3(IP_X(:,1),IP_X(:,2),Nbig(2,:),'ob','MarkerFaceColor','blue')
plot3(IP_X(:,1),IP_X(:,2),Nbig(3,:),'or','MarkerFaceColor','red')