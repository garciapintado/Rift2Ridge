% GEOMETRY PLOT
% Supports both MILAMIN_rift and M2TRI input

%--------------------------------------------------------------------------
% Function written by Miguel Andres-Martinez, PhD student at Royal Holloway
% University of London, 17-02-2015. Email: mandresmartinez87@gmail.com
%--------------------------------------------------------------------------
function plot_geomF(GEOMETRY, Geo_id, geoids, pch, cex)
% Load variables from the correct format
  if nargin < 3
      geoids = unique(Geo_id);
  end  
  if nargin < 4
      pch = ".";
  end
  if nargin < 5
      cex = 20;
  end
  
  if size(GEOMETRY,1) == 1                       % kinedyn: GEOMETRY structure
    geometry_p = GEOMETRY.bnd;
    Geo_id = GEOMETRY.id;
  else
    geometry_p = GEOMETRY;                       % rift2ridge2D
  end

  km = 1000.;

  % Generate random colors for the different interfaces
  color_r = rand(max(Geo_id),3);
  linestyle = pch + "-"; 

  % Loop to plot the different interfaces
  for layer_plot = 1:max(Geo_id)
    if ismember(layer_plot, geoids)
        plot(geometry_p(1,Geo_id==layer_plot)/km, ...
            geometry_p(2,Geo_id==layer_plot)/km, linestyle, ...
            'Color',color_r(layer_plot,:),'LineWidth',1,'MarkerSize',cex)
        if true
            text(mean(geometry_p(1,Geo_id==layer_plot)/km), mean(geometry_p(2,Geo_id==layer_plot)/km), ...
                 string(layer_plot),'Fontsize',20);
        end
    end
    hold on;
  end

  % Headings
  title('Geometry')
  xlabel('Distance [Km]')
  ylabel('Depth [Km]')

  hold off
end % function plot_geomF