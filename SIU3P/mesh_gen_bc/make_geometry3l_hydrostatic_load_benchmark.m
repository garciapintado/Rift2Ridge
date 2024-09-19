function [GEOMETRY, Geo_id] = make_geometry3l_hydrostatic_load_benchmark(x_min, x_max, y_min1, y_max1, ...
    y_max2, y_max3, cont_points, ini_deformation, Rho, sigma_moho, topo_moho)

%interfaces
xx_moho = linspace(x_min,x_max,cont_points(3));                            % x: rightward

switch ini_deformation
    case {0,3,4,5,6}                                                       % No initial topography for any interfaces
        moho    = y_max1*ones(1,cont_points(3));                           % y
        surface = y_max3*ones(1,cont_points(7));                           % y
    case {1,2}                                                             % Sinusoidal negative topography for the Moho
        moho = y_max1 + topo_moho*exp(-xx_moho.^2./(sigma_moho^2));        % WARNING: only valid for 0-centered domains
        moho([1,end]) = y_max1;
        if ini_deformation == 1                                                % 1: only affects Moho
            surface = y_max3*ones(1,cont_points(7));
        else                                                                   % 2: also affects surface
            surface_m = (moho-y_max1)*(Rho(2)-Rho(1))/Rho(3);
            surface = interp1(xx_moho, surface_m, ...
                              linspace(x_min,x_max,cont_points(7)));
            clear surface_m;
        end
end

%layer1
LAYER1_1    = [ linspace(x_min,x_max,cont_points(1));                      % rightward
                y_min1*ones(1,cont_points(1))];
Geo_id      = ones(1, size(LAYER1_1,2));
LAYER1_2    = [ x_max*ones(1,cont_points(2));
                linspace(y_min1,y_max1,cont_points(2))];                   % upward
LAYER1_2(:,[1 end]) = [];
Geo_id      = [Geo_id 2*ones(1, size(LAYER1_2,2))];
LAYER1_3    = [xx_moho; moho];                                             % rightward
Geo_id      = [Geo_id 3*ones(1, size(LAYER1_3,2))];
LAYER1_4    = [ x_min*ones(1,cont_points(2));            
                linspace(y_min1,y_max1,cont_points(2))];                   % upward
LAYER1_4(:,[1 end]) = [];
Geo_id      = [Geo_id 4*ones(1, size(LAYER1_4,2))];

GEOMETRY    = [LAYER1_1 LAYER1_2 LAYER1_3 LAYER1_4];
clear 'LAYER1_1' 'LAYER1_2' 'LAYER1_3' 'LAYER1_4'

%layer2
LAYER2_2    = [ x_max*ones(1,cont_points(4));
                linspace(y_max1,y_max2,cont_points(4))];                   % upward
LAYER2_2(:,[1 end]) = [];
Geo_id      = [Geo_id 5*ones(1, size(LAYER2_2,2))];
LAYER2_3    = [ linspace(x_min,x_max,cont_points(5));                      % rightward
                y_max2*ones(1,cont_points(5))];
Geo_id      = [Geo_id 6*ones(1, size(LAYER2_3,2))];
LAYER2_4    = [ x_min*ones(1,cont_points(4));
                linspace(y_max1,y_max2,cont_points(4))];                   % upward 
LAYER2_4(:,[1 end]) = [];
Geo_id      = [Geo_id 7*ones(1, size(LAYER2_4,2))];

GEOMETRY    = [GEOMETRY LAYER2_2 LAYER2_3 LAYER2_4];
clear 'LAYER2_2' 'LAYER2_3' 'LAYER2_4'

%layer3
LAYER3_2    = [ x_max*ones(1,cont_points(6));
                linspace(y_max2,y_max3,cont_points(6))];                   % upward
LAYER3_2(:,[1 end]) = [];
Geo_id      = [Geo_id 8*ones(1, size(LAYER3_2,2))];
LAYER3_3    = [ linspace(x_min,x_max, cont_points(7)); surface];           % rigtward
nnodl3 = length(LAYER3_3);
% central lake depression
lakeid = round([1/4 1/2 3/4]*nnodl3);                                      % lake definition indices [left center, right]
lakehwid = diff(LAYER3_3(1,lakeid(1:2)));                                  % lake half width
lakeids = lakeid(1):lakeid(3);
LAYER3_3(2,lakeids) =  y_max3  - (y_max3 - y_max2)*4/5 * exp(-0.5*(LAYER3_3(1,lakeids)-LAYER3_3(1,lakeid(2))).^2 / (0.3*lakehwid)^2 );

Geo_id      = [Geo_id 9*ones(1, size(LAYER3_3,2))];
LAYER3_4    = [ x_min*ones(1,cont_points(6));
                linspace(y_max2,y_max3,cont_points(6))];                   % upward
LAYER3_4(:,[1 end]) = [];
Geo_id      = [Geo_id 10*ones(1, size(LAYER3_4,2))];

GEOMETRY    = [GEOMETRY LAYER3_2 LAYER3_3 LAYER3_4];

end
