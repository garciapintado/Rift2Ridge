% VELOCITY QUIVER PLOT WITH DENSITIES
%   Optional input variables:
%
%               USE                     OPTIONS                 DEFAULTS
%
% color_int     Color of the box        #interfaces x 3 vector  Black
%               and interfaces          with values [0 1]
%
% line_width    Width of the box        #interfaces x 1 vector  1
%               and interfaces          with width values
%
% c_ar          Color of the arrows     #interfaces x 3 vector  Black
%                                       with values [0 1]
%
% nod_p         Plot velocities every   #nodes scalar           1
%               nod_p nodes
%
% vel_indx      Indexes of nodes which  #nodes to plot          1:nod_p:
%               velocities need to be                           length
%               plotted. Defined using                          (GCOORD)
%               select_vel.m
%
% scale_v       Scale of the arrows     scalar                  1

%--------------------------------------------------------------------------
% Function written by Lars Rüpke. Edited by Miguel Andres-Martinez, 
% 17-02-2015. Email: mandresmartinez87@gmail.com
%--------------------------------------------------------------------------

%==========================================================================
% PLOT
%==========================================================================
% Define the interval of nodes for plotting in case int_p doesn't exist
if ~exist('nod_p','var')
    nod_p = 1;
end

% Define the arrow color for the plot in case color_int doesn't exist
if ~exist('c_ar','var')
    c_ar = 'k';
end

% If a vector vel_indx doesn't exist creates it
if ~exist('vel_indx','var')
    vel_indx = 1:nod_p:length(GCOORD);
end

% If a vector scale_v doesn't exist creates it
if ~exist('scale_v','var')
    scale_v = 1;
end

% Plot velocity
% -------------
title(['Velocity fields (',num2str(istep*dt/ma),' Myr)'])
xlabel('Distance [Km]')
ylabel('Depth [Km]')
%patch('faces',EL2N,'vertices',GCOORD_N'/1000,'facevertexcdata',DENSITY(:),'FaceColor','flat')
hold on
quiver(GCOORD(1,vel_indx)/km, GCOORD(2,vel_indx)/km, ...
    Vel(vel_indx*2-1)', Vel(vel_indx*2)',scale_v,'Color',c_ar);
axis tight

% plot_vx
% hold on

% Plot box and interfaces
% -----------------------
%plot_box

hold off
drawnow