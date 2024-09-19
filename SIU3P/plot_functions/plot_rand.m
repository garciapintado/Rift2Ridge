% RANDOM FIELD PLOT
%   Optional input variables:
%
%               USE                     OPTIONS                 DEFAULTS
%
% color_int     Color of the box        #interfaces x 3 vector  Black
%               and interfaces          with values [0 1]
%
% line_width    Width of the box        #interfaces x 1 vector  1
%               and interfaces          with width values

%--------------------------------------------------------------------------
% Function written by Miguel Andres-Martinez, PhD student at Royal Holloway
% University of London, 17-02-2015. Email: mandresmartinez87@gmail.com
%--------------------------------------------------------------------------

%==========================================================================
% PLOT
%==========================================================================
% Plot historic strain
% --------------------
title(['Random field (',num2str(istep*dt/ma),' Myr)'])
xlabel('Distance [Km]')
ylabel('Depth [Km]')
patch('faces',ELEM2NODE(1:3,:)','vertices',GCOORD'/1000,'facevertexcdata',SS.Random(:,1), ...
    'FaceColor','flat')
axis tight
shading flat
colorbar
hold on

% Plot box and interfaces
% -----------------------
plot_box

drawnow
hold off