% PLASTIC FINITE STRAIN AND VISCOUS FINITE STRAIN
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
% plotbox_s     To activate or          0 no box plotting       1
%               deactivate plot_box     1 box plotting

%--------------------------------------------------------------------------
% Function written by Miguel Andres-Martinez, 12-02-2018. 
% Email: andresma@uni-bremen.de
%--------------------------------------------------------------------------

%==========================================================================
% PLOT
%==========================================================================
% Plot strain rate
% ----------------
title(['Plastic and viscous finite strains (',num2str(istep*dt/ma),' Myr)'])
xlabel('Distance [km]')
ylabel('Depth [km]')

MESH.EL2NOD = ELEM2NODE;
MESH.GCOORD = GCOORD/1000;

hold on
plot_val(MESH,I2.c,size(ELEM2NODE,2),6,[0    0.3490    1.0000],Phases)
axisp = gca;
axisp.ALim = [0 axisp.ALim(2)];
plot_val(MESH,I2.p,size(ELEM2NODE,2),6,[0.9020    0.2000         0],Phases)

axis tight

% Plot box and interfaces
% -----------------------
% Define plotbox_s if doesn't exist
if ~exist('plotbox_s','var')
    plotbox_s = 1;
end

if plotbox_s
    plot_box
end

drawnow
hold off