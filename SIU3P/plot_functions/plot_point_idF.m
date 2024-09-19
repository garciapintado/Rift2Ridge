% GEOMETRY PLOT
% Supports both MILAMIN_rift and M2TRI input
%--------------------------------------------------------------------------
% Author: Javier GP: MARUM 2020-02-15
%--------------------------------------------------------------------------
function plot_point_idF(GCOORD, Point_id, pids, pch, cex)
    
    punique = unique(Point_id);
    if nargin < 3
        pids = punique;
    end  
    if nargin < 4
        pch = ".";
    end
    if nargin < 5
        cex = 20;
    end
     
    km = 1000.;
    
    % Generate random colors for the different interfaces
     
    color_r = rand(length(punique),3);
    % linestyle = pch + "-"; 
    
    % Loop to plot the different interfaces
    for i = 1:length(punique)
        if ismember(punique(i),pids)
            plot(GCOORD(1,Point_id==punique(i))/km, ...
                 GCOORD(2,Point_id==punique(i))/km, pch, ...
                 'Color',color_r(i,:),'MarkerSize',cex)
            %if true
            %    text(mean(geometry_p(1,Geo_id==layer_plot)/km), mean(geometry_p(2,Geo_id==layer_plot)/km), ...
            %         string(layer_plot),'Fontsize',20);
            % end
            hold on;
        end
    end

    % Headings
    xlabel('Distance [Km]')
    ylabel('Depth [Km]')

    hold off;
end % function plot_point_idF