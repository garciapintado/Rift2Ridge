istep = 0.1;
ma = 0.1;
% ELASTIC FACTOR (CHI) PLOT
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
% Function written by Lars Rüpke. Edited by Miguel Andres-Martinez, 
% 17-02-2015. Email: mandresmartinez87@gmail.com
%--------------------------------------------------------------------------

%==========================================================================
% MAPPING FROM INTEGRATION POINTS TO NODES - IP CANNOT BE PLOTTED
%==========================================================================
EL2N      = zeros(nel,3); % new connectivity matrix
GCOORD_N  = zeros(2,nel*3);
Chi_n      = zeros(3,nel);
EL2N(1,:) = 1:3;
nip = 6;
nnodel = 6;
[IP_X, IP_w]    = ip_triangle(nip);
[   Nbig]    = shp_triangle(IP_X, nnodel);
for i=1:nel
    is         = (i-1)*3+1; ie = (i-1)*3+3;
    GCOORD_N(:,is:ie) = GCOORD(:,ELEM2NODE([1 2 3],i));
    EL2N(i,:) = is:ie;
    Dummy      = Nbig'\CHI(i,:)';
    Chi_n(:,i)= Dummy(1:3);
end

%==========================================================================
% PLOT
%==========================================================================
% Plot viscosity
% --------------
title(['Viscosity factor (',num2str(istep*dt/ma),' Myr)'])
xlabel('Distance [Km]')
ylabel('Depth [Km]')
patch('faces',EL2N,'vertices',GCOORD_N'/1000,'facevertexcdata', ...
    Chi_n(:),'FaceColor','flat')
axis tight
colormap(jet)
shading interp
colorbar
hold on

% Plot box and interfaces
% -----------------------
%plot_box

drawnow
hold off

nnodel = 7;
