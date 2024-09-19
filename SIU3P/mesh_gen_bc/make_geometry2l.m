function [GEOMETRY, Geo_id] = make_geometry2l(x_min, x_max,y_min1, y_max1, y_max2, cont_points, ini_interf)

%interfaces
xx_moho = linspace(x_max,x_min,cont_points(3));

moho = y_max1*ones(1,cont_points(3));

surface = y_max2*ones(1,cont_points(5));

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
