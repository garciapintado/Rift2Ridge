function plot_mu_fn(GCOORD,ELEM2NODE,PLOT,Mu_all,ei,mu_min,log_data)
% VISCOSITY PLOT
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

%--------------------------------------------------------------------------
% Function written by Miguel Andres-Martinez, 
% 26-09-2016. Email: mandresmartinez87@gmail.com
%--------------------------------------------------------------------------

%==========================================================================
% MAPPING FROM INTEGRATION POINTS TO NODES - IP CANNOT BE PLOTTED
%==========================================================================
% Calculates values at nodes from the values at the integration points
Mu_n = ip2nodes(Mu_all,GCOORD,ELEM2NODE,6,3);

% Correct negative E2_n which are the result of interpolation
Mu_n(Mu_n<mu_min) = mu_min;

nel = size(ELEM2NODE,2);

EL2N = reshape(1:nel*3,3,nel)';
GCOORD_N = GCOORD(:,ELEM2NODE([1 2 3],:));

%==========================================================================
% PLOT
%==========================================================================
% Plot viscosity
% --------------
Point_id = PLOT.Point_id;
Corner_id = PLOT.Corner_id;
Cornin_id = PLOT.Cornin_id;
title(['Viscosity ',num2str(ei),' iteration'])
xlabel('Distance [km]')
ylabel('Depth [km]')
% Check if the data needs to be plotted on logarithmic color scale
if exist('log_data','var')
    if log_data
        hhp = patch('faces',EL2N,'vertices',GCOORD_N'/1000,'facevertexcdata', ...
            log10(Mu_n(:)),'FaceColor','flat');
    else
        hhp = patch('faces',EL2N,'vertices',GCOORD_N'/1000,'facevertexcdata', ...
            Mu_n(:),'FaceColor','flat');
    end
else
    hhp = patch('faces',EL2N,'vertices',GCOORD_N'/1000,'facevertexcdata', ...
        log10(Mu_n(:)),'FaceColor','flat');
end
axis tight
colormap(flipud(jet))
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
