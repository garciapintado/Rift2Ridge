% GEOMETRY PLOT
% Supports both MILAMIN_rift and M2TRI input

%--------------------------------------------------------------------------
% Function written by Miguel Andres-Martinez, PhD student at Royal Holloway
% University of London, 17-02-2015. Email: mandresmartinez87@gmail.com
%--------------------------------------------------------------------------

% Load variables from the correct format
if size(GEOMETRY,1) == 1
    geometry_p = GEOMETRY.bnd;
    Geo_id = GEOMETRY.id;
else
    geometry_p = GEOMETRY;
end

% Generate random colors for the different interfaces
color_r = rand(max(Geo_id),3);

% Loop to plot the different interfaces
for layer_plot = 1:max(Geo_id)
    plot(geometry_p(1,Geo_id==layer_plot)/km, ...
        geometry_p(2,Geo_id==layer_plot)/km,'.-', ...
        'Color',color_r(layer_plot,:),'LineWidth',1,'MarkerSize',20)
    if true
      text(mean(GEOMETRY(1,Geo_id==layer_plot)/km), mean(GEOMETRY(2,Geo_id==layer_plot)/km), ...
           string(layer_plot),'Fontsize',20)
    end
    
    hold on
end

% Headings
title('Geometry')
xlabel('Distance [Km]')
ylabel('Depth [Km]')

hold off