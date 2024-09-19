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

%--------------------------------------------------------------------------
% Function written by Elena Ross
%--------------------------------------------------------------------------

%==========================================================================
% MAPPING FROM INTEGRATION POINTS TO NODES - IP CANNOT BE PLOTTED
%==========================================================================
EL2N      = zeros(nel,3); % new connectivity matrix
GCOORD_N  = zeros(2,nel*3);
Dpl_n      = zeros(3,nel);
EL2N(1,:) = 1:3;
nip = 6;
nnodel = 6;
[IP_X, IP_w]    = ip_triangle(nip);
[   Nbig]    = shp_triangle(IP_X, nnodel);
for i=1:nel
    is         = (i-1)*3+1; ie = (i-1)*3+3;
    GCOORD_N(:,is:ie) = GCOORD(:,ELEM2NODE([1 2 3],i));
    EL2N(i,:) = is:ie;
    Dpl_n(is:ie) = Dpl(ELEM2NODE([1 2 3],i));
end

% Correct negative E2_n which are the result of interpolation
E2_n = abs(E2_n);

%==========================================================================
% PLOT
%==========================================================================
% Plot strain rate
% ----------------
title(['Depletion (',num2str(istep*dt/ma),' Myr)'])
xlabel('Distance [Km]')
ylabel('Depth [Km]')
patch('faces',EL2N,'vertices',GCOORD_N'/1000,'facevertexcdata',Dpl_n(:),'FaceColor','flat')
colormap(parula)
axis tight
shading interp
hc = colorbar;
hold on

% Plot box and interfaces
% -----------------------
plot_box

drawnow
hold off