% FINITE PLASTIC STRAIN VS TOTAL VISCO-ELASTIC PLASTIC STRAIN
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
title(['Plastic FS vs viscous FS (',num2str(istep*dt/ma),' Myr)'])
xlabel('Distance [km]')
ylabel('Depth [km]')

MESH.EL2NOD = ELEM2NODE;
MESH.GCOORD = GCOORD/1000;
        
E_pc = I2.c+I2.p;
E_pp = I2.p./E_pc;
plot_val(MESH,E_pp,size(ELEM2NODE,2),6)
colormap('redblue')
caxis([0 1])

axis tight
hc = colorbar;

title(hc,'\frac{\varepsilon_p}{\varepsilon_p+\varepsilon_v}$', ...
    'interpreter','latex','FontSize',10)

hold on

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