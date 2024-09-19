function [GEOMETRY, Geo_id] = make_geometry3l_slope(x_min, x_max,y_min1, y_max1, y_max2, y_max3, cont_points, ini_interf, Rho, sigma_moho, topo_moho)

%interfaces
xx_moho = linspace(x_max,x_min,cont_points(3));

switch ini_interf
    case {0,5,6} % No initial topography for any interfaces
        moho = y_max1*ones(1,cont_points(3));
        
        surface = y_max3*ones(1,cont_points(7));
    case 1 % Sinusoidal negative topography only for the Moho
        moho = y_max1 + topo_moho*exp(-xx_moho.^2./(sigma_moho^2));
        moho(1) = y_max1;
        moho(end) = y_max1;
        
        surface = y_max3*ones(1,cont_points(7));
    case 2 % Sinusoidal negative topography for the Moho and resulting 
        % topography for the surface
        moho = y_max1 + topo_moho*exp(-xx_moho.^2./(sigma_moho^2));
        moho(1) = y_max1;
        moho(end) = y_max1;
        
        surface_m = (moho-y_max1)*(Rho(2)-Rho(1))/Rho(3);
        surface = interp1(linspace(x_max,x_min,cont_points(3)), ...
            surface_m, linspace(x_max,x_min,cont_points(7)));
    case {3,4} % No initial topography for any interfaces, but increment of 
        % temperature on the center of the model.
        moho = y_max1*ones(1,cont_points(3));
        
        surface = y_max3*ones(1,cont_points(7));
end

% Increment for the slope
inc = 100;

xx_surf = linspace(x_max,x_min,cont_points(7));
surface = -(xx_surf-x_min)*inc/(x_max-x_min)+inc;
incM = -Rho(3)*inc/(Rho(2)-Rho(1));
MDLB = y_max1-incM;
moho = (xx_moho-x_min)*incM/(x_max-x_min)+MDLB;

%layer1
LAYER1_1    = [ linspace(x_min,x_max,cont_points(1)); y_min1*ones(1,cont_points(1))];
Geo_id      = ones(1, size(LAYER1_1,2));
LAYER1_2    = [ x_max*ones(1,cont_points(2))    ; linspace(y_min1,y_max1,cont_points(2))];
LAYER1_2(:,[1 end]) = [];
Geo_id      = [Geo_id 2*ones(1, size(LAYER1_2,2))];
LAYER1_3    = [xx_moho; moho];
Geo_id      = [Geo_id 3*ones(1, size(LAYER1_3,2))];
LAYER1_4    = [x_min*ones(1,cont_points(2))          ; linspace(y_max1,y_min1,cont_points(2))];
LAYER1_4(:,[1 end]) = [];
Geo_id      = [Geo_id 4*ones(1, size(LAYER1_4,2))];

GEOMETRY    = [LAYER1_1 LAYER1_2 LAYER1_3 LAYER1_4];
clear 'LAYER1_1' 'LAYER1_2' 'LAYER1_3' 'LAYER1_4'

%layer2
LAYER2_2    = [ x_max*ones(1,cont_points(4))    ; linspace(y_max1,y_max2,cont_points(4))];
LAYER2_2(:,[1 end]) = [];
Geo_id      = [Geo_id 5*ones(1, size(LAYER2_2,2))];
LAYER2_3    = [linspace(x_max,x_min,cont_points(5));  y_max2*ones(1,cont_points(5))];
Geo_id      = [Geo_id 6*ones(1, size(LAYER2_3,2))];
LAYER2_4    = [x_min*ones(1,cont_points(4))          ; linspace(y_max2,y_max1,cont_points(4))];
LAYER2_4(:,[1 end]) = [];
Geo_id      = [Geo_id 7*ones(1, size(LAYER2_4,2))];

GEOMETRY    = [GEOMETRY LAYER2_2 LAYER2_3 LAYER2_4];
clear 'LAYER2_2' 'LAYER2_3' 'LAYER2_4'

%layer3
LAYER3_2    = [ x_max*ones(1,cont_points(6))    ; linspace(y_max2,y_max3,cont_points(6))];
LAYER3_2(:,[1 end]) = [];
Geo_id      = [Geo_id 8*ones(1, size(LAYER3_2,2))];
LAYER3_3    = [ xx_surf; surface];
Geo_id      = [Geo_id 9*ones(1, size(LAYER3_3,2))];
LAYER3_4    = [x_min*ones(1,cont_points(6))          ; linspace(y_max3,y_max2,cont_points(6))];
LAYER3_4(:,[1 end]) = [];
Geo_id      = [Geo_id 10*ones(1, size(LAYER3_4,2))];

GEOMETRY    = [GEOMETRY LAYER3_2 LAYER3_3 LAYER3_4];
clear 'LAYER3_2' 'LAYER3_3' 'LAYER3_4'
