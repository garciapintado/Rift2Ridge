% DENSITY PLOT
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
% Function written by Lars R??pke. Edited by Miguel Andres-Martinez, 
% 25-04-2015. Email: mandresmartinez87@gmail.com
%--------------------------------------------------------------------------

%==========================================================================
% MAPPING FROM INTEGRATION POINTS TO NODES - IP CANNOT BE PLOTTED
%==========================================================================
EL2N      = zeros(nel,3); % new connectivity matrix
GCOORD_N  = zeros(2,nel*3);
Rho_n     = zeros(3,nel);
EL2N(1,:) = 1:3;
nip = 6;
nnodel = 3;
[IP_X, IP_w]    = ip_triangle(nip);
[   Nbig]    = shp_triangle(IP_X, nnodel);
for i=1:nel
    is         = (i-1)*3+1; ie = (i-1)*3+3;
    GCOORD_N(:,is:ie) = GCOORD(:,ELEM2NODE([1 2 3],i));
    EL2N(i,:) = is:ie;
    Dummy      = Nbig'\RHO(i,:)';
    Rho_n(:,i)= Dummy(1:3);
end

% Correct negative E2_n which are the result of interpolation
Rho_n = abs(Rho_n);

%==========================================================================
% PLOT
%==========================================================================
% Plot density
% ------------
title(['Density field (',num2str(istep*dt/ma),' Myr)'])
xlabel('Distance [Km]')
ylabel('Depth [Km]')
patch('faces',EL2N,'vertices',GCOORD_N'/1000,'facevertexcdata', ...
            Rho_n(:),'FaceColor','flat')
axis tight
colormap(flipud(jet))
shading interp
hc = colorbar;
title(hc,'Kg/m^3','interpreter','tex','FontSize',10)
hold on

% Plot box and interfaces
% -----------------------
plot_box

drawnow
hold off
