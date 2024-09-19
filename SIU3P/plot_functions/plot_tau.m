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
% Function written by Lars R??pke. Edited by Miguel Andres-Martinez, 
% 17-02-2015. Email: mandresmartinez87@gmail.com
%--------------------------------------------------------------------------

%==========================================================================
% MAPPING FROM INTEGRATION POINTS TO NODES - IP CANNOT BE PLOTTED
%==========================================================================
TAU2all = sqrt(1/2 .* (TAU_xx.^2 + TAU_yy.^2) + TAU_xy.^2);

EL2N      = zeros(nel,3); % new connectivity matrix
GCOORD_N  = zeros(2,nel*3);
TAU2_n      = zeros(3,nel);
EL2N(1,:) = 1:3;
nip = 6;
nnodel = 6;
[IP_X, IP_w]    = ip_triangle(nip);
[   Nbig]    = shp_triangle(IP_X, nnodel);
for i=1:nel
    is         = (i-1)*3+1; ie = (i-1)*3+3;
    GCOORD_N(:,is:ie) = GCOORD(:,ELEM2NODE([1 2 3],i));
    EL2N(i,:) = is:ie;
    Dummy      = Nbig'\TAU2all(i,:)';
    TAU2_n(:,i)= Dummy(1:3);
end

%==========================================================================
% PLOT
%==========================================================================
% Plot strain rate
% ----------------
title(['2nd stress invariant (',num2str(istep*dt/ma),' Myr)'])
xlabel('Distance [km]')
ylabel('Depth [km]')
% Check if the data needs to be plotted on logarithmic color scale
if exist('log_data','var')
    if log_data
        patch('faces',EL2N,'vertices',GCOORD_N'/1000,'facevertexcdata', ...
            log10(TAU2_n(:)),'FaceColor','flat')
    else
        patch('faces',EL2N,'vertices',GCOORD_N'/1000,'facevertexcdata', ...
            TAU2_n(:),'FaceColor','flat')
    end
else
    patch('faces',EL2N,'vertices',GCOORD_N'/1000,'facevertexcdata', ...
        TAU2_n(:),'FaceColor','flat')
end
colormap(jet)
axis tight
shading interp
hc = colorbar;
if exist('log_data','var')
    if log_data
        title(hc,'$log(\tau_{II})$','interpreter','latex', ...
            'FontSize',10)
    end
else
    title(hc,'$\tau_{II}$','interpreter','latex', ...
        'FontSize',10)
end
hold on

% Plot box and interfaces
% -----------------------
plot_box

drawnow
hold off