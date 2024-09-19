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

% % Correct negative E2_n which are the result of interpolation
% Mu_n = abs(Mu_n);

nel = size(ELEM2NODE,2);

EL2N = reshape(1:nel*3,3,nel)';
GCOORD_N = GCOORD(:,ELEM2NODE([1 2 3],:));

%==========================================================================
% PLOT
%==========================================================================
% Plot viscosity
% --------------
title(['Viscosity field (',num2str(istep*dt/ma),' Myr)'])
xlabel('Distance [km]')
ylabel('Depth [km]')
% Check if the data needs to be plotted on logarithmic color scale
if exist('log_data','var')
    if log_data
        patch('faces',EL2N,'vertices',GCOORD_N'/1000,'facevertexcdata', ...
            log10(Mu_n(:)),'FaceColor','flat')
    else
        patch('faces',EL2N,'vertices',GCOORD_N'/1000,'facevertexcdata', ...
            Mu_n(:),'FaceColor','flat')
    end
else
    patch('faces',EL2N,'vertices',GCOORD_N'/1000,'facevertexcdata', ...
        log10(Mu_n(:)),'FaceColor','flat')
end
axis tight
colormap(flipud(jet))
shading interp
hc = colorbar;
hold on

% Plot box and interfaces
% -----------------------
plot_box

drawnow
hold off

% Plot contours
hold on
x_=GCOORD_N(1,:)';
y_=GCOORD_N(2,:)';
z_=log10(Mu_n(:));
xi=linspace(min(x_),max(x_),100);
yi=linspace(min(y_),max(y_),100);
[XI YI]=meshgrid(xi,yi);
ZI = griddata(x_,y_,z_,XI,YI);

%figure(12);
%[cs,h]=contour(XI/1000,YI/1000,ZI, 'LineWidth',1,'Color','w')
%[cs,h]=contour(XI/1000,YI/1000,ZI,[21.5,21.5], 'LineWidth',1,'Color','w')
[cs,h]=contour(XI/1000,YI/1000,ZI,[21,21], 'LineWidth',1,'Color','w');
%contour(XI/1000,YI/1000,ZI,[21.0,21.0], 'LineWidth',1,'Color','m')
hold on
%clabel(cs,'FontSize', 15)
cc=clabel(cs,h,'Color','w');%,'FontSize', 5);
%set(cc,'BackgroundColor','w');
set(cc,'Margin',1);
drawnow
