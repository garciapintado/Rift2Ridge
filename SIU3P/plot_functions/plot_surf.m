% Plot surface with velocities

Surfx = [GCOORD(1,Corner_id(3)) GCOORD(1,Corner_id(4)) GCOORD(1,Point_id==12)];
Surfy = [GCOORD(2,Corner_id(3)) GCOORD(2,Corner_id(4)) GCOORD(2,Point_id==12)];
[Surfx, Surford] = sort(Surfx);
Surfy = Surfy(Surford);
Velx_surf = Vel(1:2:end-1);
Velx_surf = Velx_surf(Point_id==max(Point_id)-1);
Vely_surf = Vel(2:2:end);
Vely_surf = Vely_surf(Point_id==max(Point_id)-1);

figure(15)
plot(Surfx/km,Surfy/km,'k');
hold on
% quiver(GCOORD(1,Point_id==12)/km, GCOORD(2,Point_id==12), Velx_surf', Vely_surf',0.2,'k');
% axis tight
title('Surface Velocities')
ylabel('Topography [m]')
xlabel('Distance [km]')