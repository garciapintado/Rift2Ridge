% VELOCITY QUIVER PLOT IN A REGULAR MESH
%   Optional input variables:
%
%               USE                     OPTIONS                 DEFAULTS
%
% c_ar          Color of the arrows     #interfaces x 3 vector  Black
%                                       with values [0 1]
%
% scale_v       Scale of the arrows     scalar                  1e3x3600x
%                                                               24x365.25
%
% auto_scale    Use autoscale in quiver 0 no                    0
%                                       1 yes
%
% xx, yy        Vectors for regular     vectors of equally      An equally 
%               mesh                    spaced points along     spaced mesh
%                                       x- and y-axis           of 100 
%                                                               points in x
%                                                               and y axis
%
% PhaseP        Phases in which         vector                  all the
%               velocities are needed                           phases
%
% Rinterf       X-velocities relative   0: absolut velocity     0
%               to the defined          IDint: relative to
%               interface               interface with id IDint

%--------------------------------------------------------------------------
% Function written by Miguel Andres-Martinez. 22-01-2018. 
% Email: mandresmartinez87@gmail.com
%--------------------------------------------------------------------------

%==========================================================================
% PLOT
%==========================================================================
% Define the arrow color for the plot in case color_int doesn't exist
function plot_velF(GCOORD, ELEM2NODE, Phases, Vel, scale_v, auto_scale)
    km = 1000.;
    year = 3600*24*365.25;
    
    % If a vector scale_v doesn't exist creates it
    if nargin < 5
        scale_v = 1e3*year; % 1 kyr
    end
  
      % If a auto_scale doesn't exist creates it
    if nargin < 6
        auto_scale = false;
    end

    if ~exist('c_ar','var')
        c_ar = 'k';
    end

    

    % Check type of velocity
    if ~exist('Rinterf','var')
        Rinterf = 0;
    end

    % Check if relative vy choice vyc exist or not
    if ~exist('vyc','var')
        vyc = 1;
    end

    % If a vector XX doesn't exist creates a default mesh
    if ~exist('xx','var')
        xx = linspace(min(GCOORD(1,:)),max(GCOORD(1,:)),500);
        yy = linspace(min(GCOORD(2,:)),max(GCOORD(2,:)),500);
        XX = repmat(xx,length(yy),1);
        YY = repmat(yy',1,length(xx));
    else
        XX = repmat(xx,length(yy),1);
        YY = repmat(yy',1,length(xx));
    end
    XX = XX(:)';
    YY = YY(:)';

    % If a vector PhaseP doesn't exist creates a PhaseP
    if ~exist('PhaseP','var')
        PhaseP = unique(Phases);
    end

    % Find points of the regular mesh in the triangular mesh
    Tris                = tsearch2(GCOORD,uint32(ELEM2NODE(1:3,:)),[XX; YY]);
    % Find points of the regular mesh in the appropriate phase
    Phase_tris          = zeros(size(Tris));
    Phase_tris(Tris~=0) = Phases(Tris(Tris~=0));
    Points2plot = Tris~=0 & ismember(Phase_tris,PhaseP);

    % Interpolate velocities at the trackpoints
    Vx_tp = remesh_val(Tris(Points2plot),GCOORD, ...
        [XX(Points2plot); YY(Points2plot)],Vel(1,:),ELEM2NODE);
    Vy_tp = remesh_val(Tris(Points2plot),GCOORD, ...
        [XX(Points2plot); YY(Points2plot)],Vel(2,:),ELEM2NODE);

    switch Rinterf
        case 0
        otherwise
            int = find(Point_id==Rinterf);
            [intx,xsort] = sort(GCOORD(1,int));
            int2nodes = int(xsort);
            INT = [intx; GCOORD(2,int2nodes)];
            Vxint = Vel(1,int2nodes);
            Vyint = Vel(2,int2nodes);
            VXINT = interp1(INT(1,:),Vxint,XX(Points2plot));
            VYINT = interp1(INT(1,:),Vyint,XX(Points2plot));
            Vx_tp = Vx_tp-VXINT;
            Vy_tp = Vy_tp-vyc*VYINT;
    end

    % Plot velocity
    % -------------
    %title(['Velocity fields (',num2str(istep*dt/ma),' Myr)'])
    xlabel('Distance [km]')
    ylabel('Depth [km]')
    %patch('faces',EL2N,'vertices',GCOORD_N'/1000,'facevertexcdata',DENSITY(:),'FaceColor','flat')
    hold on
    if auto_scale
        quiver(XX(Points2plot)/km, YY(Points2plot)/km, ...
            Vx_tp, Vy_tp,scale_v,'Color',c_ar,'Linewidth',0.8);
    else
        quiver(XX(Points2plot)/km, YY(Points2plot)/km, ...
            Vx_tp*scale_v,Vy_tp*scale_v,'AutoScale','off','Color',c_ar,'Linewidth',0.8);
    end
    axis tight

    % plot_vx
    % hold on

    % Plot box and interfaces
    % -----------------------
    %plot_box

    hold off
    drawnow
end % function