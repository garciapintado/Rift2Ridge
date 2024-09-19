% Plots the surface of the model together with the layer below

GCO_p = GCOORD(:,Point_id==max(Point_id-1));
[GCO_p(1,:),indxs] = sort(GCO_p(1,:));
GCO_p(2,:) = GCO_p(2,indxs);

plot(GCO_p(1,:)/km,GCO_p(2,:),'k','LineWidth',2);
title('Model surface')
xlabel('Distance [Km]')
ylabel('Topography [m]')
