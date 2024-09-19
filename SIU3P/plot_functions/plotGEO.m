function plotGEO(GEO, gids, pch, cex, add_labels)
  % function plotGEO(GEO, gids, pch, cex)
  % +++ purpose +++
  %  GEOMETRY PLOT
  % 
  % Supports both MILAMIN_rift and M2TRI triangle node conventions
  %
  % INPUT
  % ---
  % GEO  :: STRUCT
  % gids :: INTEGER, OPTIONAL, DIM(:) vector indicating a subset og GEO.gid interfaces to be
  %         ploted
  % pch  :: CHARACTER, OPTIONAL symbol plot at nodes. Defaults to "."
  % cex  :: REAL, OPTIONAL symbol scaling. Defaults to 20
  
  %--------------------------------------------------------------------------
  % Author: Javier Garcia-Pintado, MARUM, 2020-03
  %--------------------------------------------------------------------------

  if nargin < 2
      gids = [GEO.gid];
  end
  
  if nargin < 3
      pch = ".";
  end
  
  if nargin < 4
      cex = 20;
  end
  
  if nargin < 5
      add_labels = true;
  end
  
 
  km = 1000.;

  gidsall = [GEO.gid];
  iboo = ismember(gidsall,gids);
  
  % Generate random colors for the different interfaces
  color_r = rand(length(GEO),3);
  linestyle = pch + "-"; 

  % Loop to plot the different interfaces
  for i = 1:length(GEO)
    if iboo(i)
        plot(GEO(i).coo(1,:)/km,GEO(i).coo(2,:)/km, linestyle, ...
            'Color',color_r(i,:),'LineWidth',1,'MarkerSize',cex)
        if add_labels
            text(mean(GEO(i).coo(1,:)/km), mean(GEO(i).coo(2,:)/km), ...
                 string(GEO(i).gid),'Fontsize',20);
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