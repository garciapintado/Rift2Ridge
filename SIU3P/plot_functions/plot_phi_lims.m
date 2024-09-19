% PHI LIMITS PLOT
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
% University of London, 26-09-2016. Email: mandresmartinez87@gmail.com
%--------------------------------------------------------------------------

caxis_p = [min(min([SS.Phi1_rand; SS.Phi2_rand])) ...
    max(max([SS.Phi1_rand; SS.Phi2_rand]))]*180/pi;
%==========================================================================
% PLOT
%==========================================================================
% Plot SS.Phi1_rand
% --------------------
subplot(211)
title(['Upper \Phi limit (',num2str(istep*dt/ma),' Myr)'])
xlabel('Distance [Km]')
ylabel('Depth [Km]')
patch('faces',ELEM2NODE(1:3,:)','vertices',GCOORD'/1000,'facevertexcdata',SS.Phi1_rand(:,1)*180/pi, ...
    'FaceColor','flat')
axis tight
shading flat
caxis(caxis_p)
colorbar
hold on

% Plot box and interfaces
% -----------------------
plot_box

drawnow
hold off

% Plot SS.Phi2_rand
% --------------------
subplot(212)
title(['Lower \Phi limit (',num2str(istep*dt/ma),' Myr)'])
xlabel('Distance [Km]')
ylabel('Depth [Km]')
patch('faces',ELEM2NODE(1:3,:)','vertices',GCOORD'/1000,'facevertexcdata',SS.Phi2_rand(:,1)*180/pi, ...
    'FaceColor','flat')
axis tight
shading flat
caxis(caxis_p)
colorbar
hold on

% Plot box and interfaces
% -----------------------
plot_box

drawnow
hold off