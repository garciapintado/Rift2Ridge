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
% plot_gt       Plot geotherm in an     0 don't plot geotherm   1
%               additional subplot      1 plot geotherm

%--------------------------------------------------------------------------
% Function written by Lars R??pke. Edited by Miguel Andres-Martinez, 
% 17-02-2015. Email: mandresmartinez87@gmail.com
%--------------------------------------------------------------------------

%==========================================================================
% MAPPING FROM INTEGRATION POINTS TO NODES - IP CANNOT BE PLOTTED
%==========================================================================

EL2N      = zeros(nel,3); % new connectivity matrix
GCOORD_N  = zeros(2,nel*3);
T_n      = zeros(3,nel);
EL2N(1,:) = 1:3;
nip = 6;
nnodel = 6;
[IP_X, IP_w]    = ip_triangle(nip);
[   Nbig]    = shp_triangle(IP_X, nnodel);
for i=1:nel
    is         = (i-1)*3+1; ie = (i-1)*3+3;
    GCOORD_N(:,is:ie) = GCOORD(:,ELEM2NODE([1 2 3],i));
    EL2N(i,:) = is:ie;
    T_n(is:ie) = Temp(ELEM2NODE([1 2 3],i));    
end

%==========================================================================
% PLOT
%==========================================================================
% Plot temperature field
% ----------------------
if exist('plot_gt','var')
    if plot_gt
        subplot(2,1,1)
    end
end
title(['Temperature [C] (',num2str(istep*dt/ma),' Myr)'])
xlabel('Distance [Km]')
ylabel('Depth [Km]')
patch('faces',EL2N,'vertices',GCOORD_N'/1000,'facevertexcdata',T_n(:), ...
    'FaceColor','flat')
shading interp
axis tight
colormap(jet)
hc = colorbar;
hold on

% Plot box and interfaces
% -----------------------
plot_box

drawnow
hold off

% Plot geotherm
% -------------
if exist('plot_gt','var')
    if plot_gt
        subplot(2,1,2)
        plot(-GCOORD(2,:)/1000,Temp,'.')
        xlabel('Depth [Km]')
        ylabel('Temperature [C]')
    end
end