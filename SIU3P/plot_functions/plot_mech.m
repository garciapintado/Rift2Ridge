% DEFORMATION MECHANISMS PLOT
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
% 23-03-2015. Email: mandresmartinez87@gmail.com
%--------------------------------------------------------------------------

%==========================================================================
% PLOT
%==========================================================================
hold on
% Remove Nan
Mu_dif_allp = Mu_dif_all;
Mu_dif_allp(isnan(Mu_dif_all)) = Inf;

% Plot apparent dislocation viscosity IPs
% ---------------------------------------
dis_ip = Mu_dis_all==Mu_dis_all & Mu_dis_all<=Mu_dif_allp;
plot(GIP_x_all(dis_ip)/1000,GIP_y_all(dis_ip)/1000,'*w','linewidth',2)
plot(GIP_x_all(dis_ip)/1000,GIP_y_all(dis_ip)/1000,'*r')

% Plot apparent diffusion viscosity IPs
% -------------------------------------
dif_ip = Mu_dif_allp==Mu_dif_allp & Mu_dif_allp<=Mu_dis_all;
plot(GIP_x_all(dif_ip)/1000,GIP_y_all(dif_ip)/1000,'ow','linewidth',2)
plot(GIP_x_all(dif_ip)/1000,GIP_y_all(dif_ip)/1000,'og')

% Plot apparent brittle viscosity IPs
% -----------------------------------
brit_ip = YC==1;
plot(GIP_x_all(brit_ip)/1000,GIP_y_all(brit_ip)/1000,'.w','linewidth',2)
plot(GIP_x_all(brit_ip)/1000,GIP_y_all(brit_ip)/1000,'*b')

% Plot box and interfaces
% -----------------------
plot_box

drawnow
hold off