function plot_p_fn(GCOORD,ELEM2NODE,PLOT,Pressure,ei,P_type)
% DYNAMIC PRESSURE PLOT
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
% Function written by Lars R??pke. Edited by Miguel Andres-Martinez, 
% 17-02-2015. Email: mandresmartinez87@gmail.com
%--------------------------------------------------------------------------

%==========================================================================
% COORDINATES
%==========================================================================
nel = size(ELEM2NODE,2);

EL2N = reshape(1:nel*3,3,nel)';
GCOORD_N = GCOORD(:,ELEM2NODE([1 2 3],:));

%==========================================================================
% PLOT
%==========================================================================
% Plot dynamic pressure
% ---------------------
Point_id = PLOT.Point_id;
Corner_id = PLOT.Corner_id;
Cornin_id = PLOT.Cornin_id;
title([P_type,' ',num2str(ei),' iteration'])
xlabel('Distance [km]')
ylabel('Depth [km]')

% Plot
hhp = patch('faces',EL2N,'vertices',GCOORD_N'/1000,'facevertexcdata', ...
    Pressure(:)/1e6,'FaceColor','flat');

colormap(jet)
axis tight
shading interp
hc = colorbar;
title(hc,'P [MPa]','interpreter','latex','FontSize',10)
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