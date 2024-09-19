% Y-VELOCITY FIELD PLOT
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

%--------------------------------------------------------------------------
% Function written by Lars Rüpke. Edited by Miguel Andres-Martinez, 
% 17-02-2015. Email: mandresmartinez87@gmail.com
%--------------------------------------------------------------------------

%==========================================================================
% MAPPING FROM INTEGRATION POINTS TO NODES - IP CANNOT BE PLOTTED
%==========================================================================
EL2N      = zeros(nel,3); % new connectivity matrix
GCOORD_N  = zeros(2,nel*3);
for i=1:nel
    is         = (i-1)*3+1; ie = (i-1)*3+3;
    GCOORD_N(:,is:ie) = GCOORD(:,ELEM2NODE([1 2 3],i));
    EL2N(i,:) = is:ie;
end

%==========================================================================
% PLOT
%==========================================================================
% Plot y-velocity field
% ---------------------
Vel_y = DISPL(2,ELEM2NODE(1:3,:));
title (['Vy [Km/Myr] (',num2str(istep*dt/ma),' Myr)'])
xlabel('Distance [Km]')
ylabel('Depth [Km]')
patch('faces',EL2N,'vertices',GCOORD_N'/1000,'facevertexcdata', ...
    Vel_y(:)*ma/km,'FaceColor','flat')
shading interp
colormap(jet)
hold on
hc = colorbar;
title(hc,'[km/Myr]','FontSize',10)

% Plot box and interfaces
% -----------------------
plot_box

hold off
drawnow