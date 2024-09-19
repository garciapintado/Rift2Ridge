function [preexp_dif,preexp_dis,TEMP_IP] = ...
    ss_dep_temp(Temp,PHY,ELEM2NODE,nel,nip)
% PREEXP = SS_DEP_TEMP(TEMP,PHY) calculates maximum factors PREEXP for the
% flow-law preexponential factor used in viscous strain softening by
% setting PREEXP as an Arrhenius-like function of the temperature TEMP.

%--------------------------------------------------------------------------
% Function written by Leon Liu, PhD student at Royal Holloway
% University of London, 2015.
%--------------------------------------------------------------------------

% 2015 ER & AM
% Adaptation to RIFT code
% 2016 DD
% Improve some aspects of it
% 17/06/2016 MA
% Improve inputs, comment and add independent control for diffusion and
% dislocation

%==========================================================================
% INITIALIZATION
%==========================================================================
% Initialize variables
nnodel = size(ELEM2NODE,1);

% Calculate shape functions
%[IP_X,~]    = ip_triangle_m2tri(nip);
%[N,~]       = sf_dsf_tri367_N(IP_X,6,'cell');
[IP_X,~] = ip_triangle(nip);
[N,~]    = shp_deriv_triangle(IP_X,6);

% Reshape node temperatures to ELEM2NODE shape
TEMP_EL = Temp(ELEM2NODE(1:6,:));
% Initialize temperatures at integration points
TEMP_IP = zeros(nel,nip);
% % Uncomment to plot
% GIP_x = zeros(nel,nip);
% GIP_y = zeros(nel,nip);
% ECOORD_x = reshape(GCOORD(1,ELEM2NODE),nnodel,nel);
% ECOORD_y = reshape(GCOORD(2,ELEM2NODE),nnodel,nel);

% Integration loop to calculate temperatures at the integration points
for ip=1:nip
    Ni = N{ip};
    % Temp at integration point
    TEMP_IP(:,ip) = (Ni'*TEMP_EL)';
%     % Uncomment to plot
%     GIP_x(:,ip) = Ni'*ECOORD_x;
%     GIP_y(:,ip) = Ni'*ECOORD_y;
end

% % Uncomment to plot
% Gid = GIP_x>=-5000 & GIP_x<=5000 & GIP_y>= -60000 & GIP_y <= 0;
% scatter(GIP_x(Gid)/1000,GIP_y(Gid)/1000,5,TEMP_IP(Gid))
% colormap(jet)
% colorbar

%==========================================================================
% PREEXP CALCULATIONS
%==========================================================================
% Diffusion
% ---------
% Load limits
lower_lim = PHY.SS.dif_low_lim;
upper_lim = PHY.SS.dif_upp_lim;

% Arrhenius function
A=log(PHY.SS.Pef_dif(2))*upper_lim/(lower_lim-upper_lim);
preexp_dif = exp(A*(TEMP_IP-upper_lim)/upper_lim);

% Limit the function to the minimum and maximum Pef
preexp_dif(TEMP_IP<=lower_lim) = PHY.SS.Pef_dif(2);
preexp_dif(TEMP_IP>=upper_lim) = PHY.SS.Pef_dif(1);

% Dislocation
% -----------
% Load limits
lower_lim = PHY.SS.dis_low_lim;
upper_lim = PHY.SS.dis_upp_lim;

% Arrhenius function
A=log(PHY.SS.Pef_dis(2))*upper_lim/(lower_lim-upper_lim);
preexp_dis = exp(A*(TEMP_IP-upper_lim)/upper_lim);

% Limit the function to the minimum and maximum Pef
preexp_dis(TEMP_IP<=lower_lim) = PHY.SS.Pef_dis(2);
preexp_dis(TEMP_IP>=upper_lim) = PHY.SS.Pef_dis(1);

% % Plot (uncomment to plot)
% % ========================
% subplot(1,2,2)
% title('Diffusion')
% plot(preexp_dif,TEMP_IP,'.k')
% set(gca,'Ydir','reverse')
% subplot(1,2,1)
% title('Dislocation')
% plot(preexp_dis,TEMP_IP,'.k')
% set(gca,'Ydir','reverse')
