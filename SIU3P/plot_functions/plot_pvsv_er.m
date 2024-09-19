% PLASTIC STRAIN RATE AND VISCOUS STRAIN RATE
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
% plotbox_s     To activate or          0 no box plotting       1
%               deactivate plot_box     1 box plotting

%--------------------------------------------------------------------------
% Function written by Miguel Andres-Martinez, 12-02-2018. 
% Email: andresma@uni-bremen.de
%--------------------------------------------------------------------------

%==========================================================================
% PLOT
%==========================================================================
% Plot strain rate
% ----------------
title(['Plastic and viscous strain rates (',num2str(istep*dt/ma),' Myr)'])
xlabel('Distance [km]')
ylabel('Depth [km]')

MESH.EL2NOD = ELEM2NODE;
MESH.GCOORD = GCOORD/1000;
MESH.sign = S_pvsv;
        
% Calculate plastic strain rates
% ------------------------------
Er_xx   = zeros(size(ELEM2NODE,2),nip);
Er_xy   = zeros(size(ELEM2NODE,2),nip);
Er_yy   = zeros(size(ELEM2NODE,2),nip);
Er_xx(YC~=0)  = 0.5 * Gamma(YC~=0)...
.*TAU_xx(YC~= 0)./Yield_T2(YC~=0);
Er_xy(YC~=0)  = 0.5 * Gamma(YC~=0)...
.*TAU_xy(YC~=0)./Yield_T2(YC~= 0);
Er_yy(YC~=0)  = 0.5 * Gamma(YC~= 0)...
.*TAU_yy(YC~=0)./Yield_T2(YC~=0);
ErP     = sqrt(0.5*(Er_xx.^2 + Er_yy.^2) + Er_xy.^2);

% Calculate viscous strain rates
% ------------------------------
diff_s = ~isnan(Mu_dif_all);
Mu_c_all = Mu_dis_all;

Mu_c_all(diff_s) = (1./Mu_dis_all(diff_s)+1./Mu_dif_all(diff_s)).^(-1);

ErV = sqrt(0.5*(TAU_xx.^2+TAU_yy.^2)+TAU_xy.^2)./(2*Mu_c_all);

hold on
plot_val(MESH,ErP,size(ELEM2NODE,2),6,[0.9020    0.2000         0])
axisp = gca;
axisp.ALim = [0 axisp.ALim(2)];
plot_val(MESH,ErV,size(ELEM2NODE,2),6,[0    0.3490    1.0000])

axis tight

% Plot box and interfaces
% -----------------------dif
% Define plotbox_s if doesn't exist
if ~exist('plotbox_s','var')
    plotbox_s = 1;
end

if plotbox_s
    plot_box
end

drawnow
hold off
