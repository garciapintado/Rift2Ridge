function [GEOMETRY,Geo_id,cont_points,Elsizes,x_min,x_max] = ...
    make_geometrySCS(Ginput)

%==========================================================================
% INPUT
%==========================================================================
km = Ginput.km;
H = Ginput.H*km;
D = Ginput.D*km;
d = Ginput.d*km;
M = Ginput.M*km;
Mh = Ginput.Mh*km;
t = Ginput.t*km;
At = Ginput.At*km;
Ohu = Ginput.Ohu*km;
Ohl = Ginput.Ohl*km;
To = Ginput.To*km;
Do = Ginput.Do*km;
UCR = Ginput.UCR;
LCR = Ginput.LCR;
DWOM = Ginput.DWOM;
y_min1 = Ginput.y_min1;
y_max1 = Ginput.y_max1;
y_max2 = Ginput.y_max2;
y_max3 = Ginput.y_max3;
Rho = Ginput.Rho;
Elsizes = Ginput.Elsizes;

%==========================================================================
% GEOMETRY
%==========================================================================
% Resolution
% ----------
x_min = -(At/2+D);
x_max = At/2+d+M+To+Do;
% Number of points at the boundaries
cont_points = [floor((x_max-x_min)/DWOM) floor((y_max1-y_min1)/DWOM) ...
    floor((x_max-x_min)/LCR) ...
    floor((y_max2-y_max1)/LCR) floor((x_max-x_min)/UCR) ...
    floor((y_max3-y_max2)/UCR) floor((x_max-x_min)/UCR)]+1;

% Thicknesses of regular continental crust
H1 = y_max3-y_max1;  % Total thickness of the crust
hu1 = y_max3-y_max2; % Upper
hl1 = y_max2-y_max1; % Lower

% Arc
% ---
% Surface
xx_surface = linspace(x_max,x_min,cont_points(7));
surface = y_max3 + t*exp(-xx_surface.^2./((At/5)^2));
surface(1) = y_max3;
surface(end) = y_max3;
surface(isnan(surface)) = 0;

% Calculations
delta_rho_u2 = Rho(3)-Rho(1)+hl1/hu1*(Rho(2)-Rho(1));
hu2 = (hu1*(Rho(3)-Rho(1))+hl1*(Rho(2)-Rho(1))-surface*Rho(1))./ ...
    delta_rho_u2;
hl2 = hu2*(hl1/hu1);
Lt = hu2-surface-hu1;
Mt = hl2+Lt-hl1;

% Upper-lower crust interface
xx_ulcrust = linspace(x_max,x_min,cont_points(5));
ulcrust = y_max2-interp1(xx_surface,Lt,xx_ulcrust,'spline');

% Moho
xx_moho = linspace(x_max,x_min,cont_points(3));
moho = y_max1-interp1(xx_surface,Mt,xx_moho,'spline');

% % Moho
% xx_moho = linspace(x_max,x_min,cont_points(3));
% moho_m = (surface-y_max3)*(Rho(3)/(Rho(2)-Rho(1)))+y_max1;
% moho = interp1(xx_surface,moho_m,xx_moho);

% Continent
% ---------
H2 = H;
hu2 = H2*hu1/H1;
hl2 = H2*hl1/H1;
mc = ((hu2-hu1)*Rho(3)+(hl2-hl1)*Rho(2))/Rho(1);
tc = H2-H1-mc;
lc = hu2-hu1-tc;
D_indx_s = xx_surface>=x_min & xx_surface<=x_min+D;
D_indx_ul = xx_ulcrust>=x_min & xx_ulcrust<=x_min+D;
D_indx_m = xx_moho>=x_min & xx_moho<=x_min+D;
Tc = interp1([x_min x_min+D],[tc 0],xx_surface(D_indx_s));
Lc = interp1([x_min x_min+D],[lc 0],xx_ulcrust(D_indx_ul));
Mc = interp1([x_min x_min+D],[mc 0],xx_moho(D_indx_m));
surface(D_indx_s) = surface(D_indx_s)+Tc;
ulcrust(D_indx_ul) = ulcrust(D_indx_ul)-Lc;
moho(D_indx_m) = moho(D_indx_m)-Mc;

% Proto-SCS margin
% ----------------
H2 = Mh;
hu2 = H2*hu1/H1;
hl2 = H2*hl1/H1;
mm = ((hu2-hu1)*Rho(3)+(hl2-hl1)*Rho(2))/Rho(1);
tm = H2-H1-mm;
lm = hu2-hu1-tm;
D_indx_s = xx_surface>=x_min+D+At+d & xx_surface<=x_max-To-Do;
D_indx_ul = xx_ulcrust>=x_min+D+At+d & xx_ulcrust<=x_max-To-Do;
D_indx_m = xx_moho>=x_min+D+At+d & xx_moho<=x_max-To-Do;
Tm = interp1([x_min+D+At+d x_max-To-Do],[0 tm],xx_surface(D_indx_s));
Lm = interp1([x_min+D+At+d x_max-To-Do],[0 lm],xx_ulcrust(D_indx_ul));
Mm = interp1([x_min+D+At+d x_max-To-Do],[0 mm],xx_moho(D_indx_m));
surface(D_indx_s) = surface(D_indx_s)+Tm;
ulcrust(D_indx_ul) = ulcrust(D_indx_ul)-Lm;
moho(D_indx_m) = moho(D_indx_m)-Mm;

% Proto-SCS transition
% --------------------
hu1 = hu2;
hl1 = hl2;
hu2 = Ohu;
hl2 = Ohl;
mt = (Rho(3)*(hu1-hu2)+Rho(2)*(hl1-hl2))/Rho(1);
moho_o = y_max1-mt-mm;
ul_o = y_max1-mt-mm+hl2;
surface_o = y_max1-mt-mm+hl2+hu2;
D_indx_s = xx_surface>=x_max-To-Do & xx_surface<=x_max-Do;
D_indx_ul = xx_ulcrust>=x_max-To-Do & xx_ulcrust<=x_max-Do;
D_indx_m = xx_moho>=x_max-To-Do & xx_moho<=x_max-Do;
surface(D_indx_s) = interp1([x_max-To-Do x_max-Do],[y_max1-mm+hl1+hu1 surface_o],xx_surface(D_indx_s));
ulcrust(D_indx_ul) = interp1([x_max-To-Do x_max-Do],[y_max1-mm+hl1 ul_o],xx_ulcrust(D_indx_ul));
moho(D_indx_m) = interp1([x_max-To-Do x_max-Do],[y_max1-mm moho_o],xx_moho(D_indx_m));

% Proto-SCS ocean
% ---------------
D_indx_s = xx_surface>=x_max-Do & xx_surface<=x_max;
D_indx_ul = xx_ulcrust>=x_max-Do & xx_ulcrust<=x_max;
D_indx_m = xx_moho>=x_max-Do & xx_moho<=x_max;
surface(D_indx_s) = surface_o;
ulcrust(D_indx_ul) = ul_o;
moho(D_indx_m) = moho_o;

% Plot (uncomment)
% ----------------
axis([x_min x_max y_max1-10*km y_max3+5*km])
plot(xx_surface/1000,surface/1000,'.-')
hold on
plot(xx_ulcrust/1000,ulcrust/1000,'.-')
plot(xx_moho/1000,moho/1000,'.-')
plot([-(At/2) (At/2); -(At/2) (At/2)]/1000,[y_max1 y_max1; t t]/1000, ...
    '--k')
plot([x_min x_min x_min; x_max x_max x_max]/1000, ...
    [y_max1 y_max2 y_max3; y_max1 y_max2 y_max3]/1000,'--k')
plot([x_max-M-To-Do x_max-To-Do x_max-Do; ...
    x_max-M-To-Do x_max-To-Do x_max-Do]/1000, ...
    [y_max1 y_max1 y_max1; t t t]/1000,'--k')
hold off

% % Check propotions
% % ----------------
% ulcrust_s = interp1(xx_ulcrust,ulcrust,xx_surface);
% moho_s = interp1(xx_moho,moho,xx_surface);
% ratio_uc_lc = (surface-ulcrust_s)./(ulcrust_s-moho_s);
% ratio_diff = ratio_uc_lc-hu1/hl1;
% % Check for constant uc/lc ratio along axis
% if max(abs(ratio_diff))>0.005;
%     error(['UC/LC ratios different along profile! :: ', ...
%         'Difference of: ',num2str(max(abs(ratio_diff))), ...
%         '. If this difference is still ok, change the error criterium: ', ...
%         '"if max(abs(ratio_diff))>error_criterium;" in make_geometrySCS.m'])
% end

%==========================================================================
% BUILD LAYERS
%==========================================================================
% Number of points boundaries
%   Right
r1 = round((moho(1)-y_min1)*cont_points(2)/(y_max1-y_min1));
r2 = round((ulcrust(1)-moho(1))*cont_points(4)/(y_max2-y_max1));
r3 = round((surface(1)-ulcrust(1))*cont_points(6)/(y_max3-y_max2));
%   Left
l1 = round((moho(end)-y_min1)*cont_points(2)/(y_max1-y_min1));
l2 = round((ulcrust(end)-moho(end))*cont_points(4)/(y_max2-y_max1));
l3 = round((surface(end)-ulcrust(end))*cont_points(6)/(y_max3-y_max2));

%layer1
LAYER1_1    = [ linspace(x_min,x_max,cont_points(1)); y_min1*ones(1,cont_points(1))];
Geo_id      = ones(1, size(LAYER1_1,2));
LAYER1_2    = [ x_max*ones(1,r1)    ; linspace(y_min1,moho(1),r1)];
LAYER1_2(:,[1 end]) = [];
Geo_id      = [Geo_id 2*ones(1, size(LAYER1_2,2))];
LAYER1_3    = [xx_moho; moho];
Geo_id      = [Geo_id 3*ones(1, size(LAYER1_3,2))];
LAYER1_4    = [x_min*ones(1,l1)          ; linspace(moho(end),y_min1,l1)];
LAYER1_4(:,[1 end]) = [];
Geo_id      = [Geo_id 4*ones(1, size(LAYER1_4,2))];

GEOMETRY    = [LAYER1_1 LAYER1_2 LAYER1_3 LAYER1_4];
clear 'LAYER1_1' 'LAYER1_2' 'LAYER1_3' 'LAYER1_4'

%layer2
LAYER2_2    = [ x_max*ones(1,r2)    ; linspace(moho(1),ulcrust(1),r2)];
LAYER2_2(:,[1 end]) = [];
Geo_id      = [Geo_id 5*ones(1, size(LAYER2_2,2))];
LAYER2_3    = [xx_ulcrust; ulcrust];
Geo_id      = [Geo_id 6*ones(1, size(LAYER2_3,2))];
LAYER2_4    = [x_min*ones(1,l2)          ; linspace(ulcrust(end),moho(end),l2)];
LAYER2_4(:,[1 end]) = [];
Geo_id      = [Geo_id 7*ones(1, size(LAYER2_4,2))];

GEOMETRY    = [GEOMETRY LAYER2_2 LAYER2_3 LAYER2_4];
clear 'LAYER2_2' 'LAYER2_3' 'LAYER2_4'

%layer3
LAYER3_2    = [ x_max*ones(1,r3)    ; linspace(ulcrust(1),surface(1),r3)];
LAYER3_2(:,[1 end]) = [];
Geo_id      = [Geo_id 8*ones(1, size(LAYER3_2,2))];
LAYER3_3    = [xx_surface; surface];
Geo_id      = [Geo_id 9*ones(1, size(LAYER3_3,2))];
LAYER3_4    = [x_min*ones(1,l3)          ; linspace(surface(end),ulcrust(end),l3)];
LAYER3_4(:,[1 end]) = [];
Geo_id      = [Geo_id 10*ones(1, size(LAYER3_4,2))];

GEOMETRY    = [GEOMETRY LAYER3_2 LAYER3_3 LAYER3_4];
clear 'LAYER3_2' 'LAYER3_3' 'LAYER3_4'
