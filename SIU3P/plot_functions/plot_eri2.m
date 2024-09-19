% 2ND STRAIN RATE INVARIANT PLOT
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
% log_data      Plot logarithmic data   0 no logarithmic data   1
%                                       1 logarithmic data
%
% plotbox_s     To activate or          0 no box plotting       1
%               deactivate plot_box     1 box plotting

%--------------------------------------------------------------------------
% Function written by Lars R??pke. Edited by Miguel Andres-Martinez, 
% 17-02-2015. Email: mandresmartinez87@gmail.com
%--------------------------------------------------------------------------

%==========================================================================
% MAPPING FROM INTEGRATION POINTS TO NODES - IP CANNOT BE PLOTTED
%==========================================================================
% Calculates values at nodes from the values at the integration points
Sxx = ip2nodes(STRAIN_xx,GCOORD,ELEM2NODE,6,3);
Sxy = ip2nodes(STRAIN_xy,GCOORD,ELEM2NODE,6,3);
Syy = ip2nodes(STRAIN_yy,GCOORD,ELEM2NODE,6,3);

E2_n = sqrt(0.5*(Sxx.^2 + Syy.^2) + Sxy.^2);

nel = size(ELEM2NODE,2);

EL2N = reshape(1:nel*3,3,nel)';
GCOORD_N = GCOORD(:,ELEM2NODE([1 2 3],:));

%==========================================================================
% PLOT
%==========================================================================
% Plot strain rate
% ----------------
title(['2nd strain rate invariant (',num2str(istep*dt/ma),' Myr)'])
xlabel('Distance [km]')
ylabel('Depth [km]')
% Check if the data needs to be plotted on logarithmic color scale
if exist('log_data','var')
    if log_data
        hhp = patch('faces',EL2N,'vertices',GCOORD_N'/1000,'facevertexcdata', ...
            log10(E2_n(:)),'FaceColor','flat');
    else
        hhp = patch('faces',EL2N,'vertices',GCOORD_N'/1000,'facevertexcdata', ...
            E2_n(:),'FaceColor','flat');
    end
else
    hhp = patch('faces',EL2N,'vertices',GCOORD_N'/1000,'facevertexcdata', ...
        log10(E2_n(:)),'FaceColor','flat');
end
colormap(parula)
axis tight
shading interp
hc = colorbar;
if exist('log_data','var')
    if log_data
        title(hc,'$log(\dot{\varepsilon}_{II})$','interpreter','latex', ...
            'FontSize',10)
    end
else
    title(hc,'$(\dot{\varepsilon}_{II})$','interpreter','latex', ...
        'FontSize',10)
end
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