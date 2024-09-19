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
%
% alphac        To activate plot        Any color (i.e. 'k',    No transpa-
%               transparency            'r', ...)               rency

%--------------------------------------------------------------------------
% Function written by Lars R??pke. Edited by Miguel Andres-Martinez, 
% 17-02-2015. Email: mandresmartinez87@gmail.com
%--------------------------------------------------------------------------

%==========================================================================
% MAPPING FROM INTEGRATION POINTS TO NODES - IP CANNOT BE PLOTTED
%==========================================================================
EL2N      = zeros(nel,3); % new connectivity matrix
GCOORD_N  = zeros(2,nel*3);
E2_n      = zeros(3,nel);
EL2N(1,:) = 1:3;
if exist('SOLVER','var')
    [IP_X, IP_w]    = ip_triangle(SOLVER.nip_stress);
else
    [IP_X, IP_w]    = ip_triangle(size(E2all,2));
end
[   Nbig]    = shp_triangle(IP_X, 3);
for i=1:nel
    is         = (i-1)*3+1; ie = (i-1)*3+3;
    GCOORD_N(:,is:ie) = GCOORD(:,ELEM2NODE([1 2 3],i));
    EL2N(i,:) = is:ie;
    Dummy      = Nbig'\E2all(i,:)';
    E2_n(:,i)= Dummy(1:3);
end

% Correct negative E2_n which are the result of interpolation
E2_n = abs(E2_n);

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
        E2_nl = log10(E2_n);
    else
        E2_nl = E2_n;
    end
else
    E2_nl = log10(E2_n);
end

if exist('alphac','var')
    patch('faces',EL2N,'vertices',GCOORD_N'/1000,'FaceVertexAlphaData', ...
        E2_n(:),'AlphaDataMapping','scaled','FaceAlpha','interp',...
        'FaceColor',alphac,'EdgeColor','none')
else
    patch('faces',EL2N,'vertices',GCOORD_N'/1000,'facevertexcdata', ...
        E2_nl(:),'FaceColor','flat')
    colormap(parula)
    shading interp
end
axis tight
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