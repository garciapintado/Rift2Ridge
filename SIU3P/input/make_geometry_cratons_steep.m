function [GEOMETRY,Geo_id,cont_points,Elsizes,x_min,x_max] = ...
    make_geometry_cratons()
% [GEOMETRY,GEO_ID,CONT_POINTS,ELSIZES,X_MIN,X_MAX] = MAKE_GEOMETRY_CRATONS
% is an input function to produce the geometries necessary to generate a
% triangular mesh. It reads geometries of the Temperature Solver developed
% by Albert de Monserrat and turns them into the MILAMIN_RIFT geometry
% language.

%--------------------------------------------------------------------------
% Function written by Miguel Andres-Martinez, PhD student at Royal Holloway
% University of London, 30-10-2014. Email: mandresmartinez87@gmail.com
%--------------------------------------------------------------------------

%==========================================================================
% INPUT
%==========================================================================

% File to load geometries
dir_geom = 'Input/GEOM_eq_thick40steep.mat';
% File to load data
dir_data = 'Input/DATA_eq_thick40steep.mat';

% Constants
km = 1000;

% General geometries of the model
% -------------------------------
x_max =  250*km;
x_min =  -250*km;
%layer 1
y_min1 = -400*km;
y_max1 = -80*km;
%layer 2
y_max2 =  -40*km;
%layer 3
y_max3 =  -20*km;
%layer 4
y_max4 =  0;

% Resolution
% ----------
% Number of points at the boundaries
cont_points = [50 50 100 15 100 7 100 10 1000]; 
% Number of points on geometry contours for ordered as follows:
%       __9__
%            |
%       __7__|8
%            |
%       __5__|6
%            |
%       __3__|4
%            |
%       __1__|2

% Resolution at the different phases
Elsizes     = [1e8 1e8 1e7 1e7]; % [1e8 1e7 1e7 1e7]

%==========================================================================
% LOAD DATA
%==========================================================================

load(dir_geom)
load(dir_data)

%==========================================================================
% CALCULATE GEOMETRIES
%==========================================================================

% Layer1
% ------
% Calculate layer
LAYER1_1 = [linspace(x_min,x_max,cont_points(1)); ...
    y_min1*ones(1,cont_points(1))];
Geo_id = ones(1,size(LAYER1_1,2));
LAYER1_2 = [x_max*ones(1,cont_points(2)); ...
    linspace(y_min1,y_max1,cont_points(2))];
LAYER1_2(:,[1 end]) = [];
Geo_id = [Geo_id 2*ones(1, size(LAYER1_2,2))];
LAYER1_3 = [linspace(x_max,x_min,cont_points(3)); ...
    y_max1*ones(1,cont_points(3))];
Geo_id = [Geo_id 3*ones(1, size(LAYER1_3,2))];
LAYER1_4 = [x_min*ones(1,cont_points(2)); ...
    linspace(y_max1,y_min1,cont_points(2))];
LAYER1_4(:,[1 end]) = [];
Geo_id = [Geo_id 4*ones(1, size(LAYER1_4,2))];

GEOMETRY = [LAYER1_1 LAYER1_2 LAYER1_3 LAYER1_4];
clear 'LAYER1_1' 'LAYER1_2' 'LAYER1_3' 'LAYER1_4'

% Layer2
% ------
% Calculate Moho
X = x_max:-DATA.NUM.dx:x_min;
Y = -GEOM.ycr_cor(end:-1:1);
CumLength = [0 cumsum(sqrt((X(2:end) - X(1:end-1)).^2 + (Y(2:end) ...
    - Y(1:end-1)).^2))];
n_cont_points = round(cont_points(5)*CumLength(end)/(x_max-x_min));
Sn = linspace(CumLength(1),CumLength(end),n_cont_points);
X_moho = interp1(CumLength,X,Sn,'linear');
Y_moho = interp1(CumLength,Y,Sn,'linear');
cont_right = round((Y_moho(1)-y_max1)*cont_points(4)/(y_max2-y_max1));
cont_left = round((Y_moho(end)-y_max1)*cont_points(4)/(y_max2-y_max1));

% Calculate layer
LAYER2_2 = [x_max*ones(1,cont_right); ...
    linspace(y_max1,Y_moho(1),cont_right)];
LAYER2_2(:,[1 end]) = [];
Geo_id = [Geo_id 5*ones(1,size(LAYER2_2,2))];
LAYER2_3 = [X_moho; Y_moho];
Geo_id = [Geo_id 6*ones(1,size(LAYER2_3,2))];
LAYER2_4 = [x_min*ones(1,cont_left); ...
    linspace(Y_moho(end),y_max1,cont_left)];
LAYER2_4(:,[1 end]) = [];
Geo_id = [Geo_id 7*ones(1, size(LAYER2_4,2))];

GEOMETRY = [GEOMETRY LAYER2_2 LAYER2_3 LAYER2_4];
clear 'LAYER2_2' 'LAYER2_3' 'LAYER2_4'

% Layer3
% ------
% Calculate Upper-Lower crust boundary
X = x_max:-DATA.NUM.dx:x_min;
Y = -GEOM.yucq_cor(end:-1:1);
CumLength = [0 cumsum(sqrt((X(2:end) - X(1:end-1)).^2 + (Y(2:end) ...
    - Y(1:end-1)).^2))];
n_cont_points = round(cont_points(7)*CumLength(end)/(x_max-x_min));
Sn = linspace(CumLength(1),CumLength(end),n_cont_points);
X_ulc = interp1(CumLength,X,Sn,'linear');
Y_ulc = interp1(CumLength,Y,Sn,'linear');
cont_right = round((Y_ulc(1)-Y_moho(1))*cont_points(6)/(y_max3-y_max2));
cont_left = round((Y_ulc(end)-Y_moho(end))*cont_points(6)/(y_max3-y_max2));

% Calculate layer
LAYER3_2 = [x_max*ones(1,cont_right); ...
    linspace(Y_moho(1),Y_ulc(1),cont_right)];
LAYER3_2(:,[1 end]) = [];
Geo_id = [Geo_id 8*ones(1,size(LAYER3_2,2))];
LAYER3_3 = [X_ulc; Y_ulc];
Geo_id = [Geo_id 9*ones(1,size(LAYER3_3,2))];
LAYER3_4 = [x_min*ones(1,cont_left); ...
    linspace(Y_ulc(end),Y_moho(end),cont_left)];
LAYER3_4(:,[1 end]) = [];
Geo_id = [Geo_id 10*ones(1,size(LAYER3_4,2))];

GEOMETRY = [GEOMETRY LAYER3_2 LAYER3_3 LAYER3_4];
clear 'LAYER3_2' 'LAYER3_3' 'LAYER3_4'

% Layer4
% ------
% Calculate topography
X = x_max:-DATA.NUM.dx:x_min;
Y = GEOM.ytopo_cor(end:-1:1);
CumLength = [0 cumsum(sqrt((X(2:end) - X(1:end-1)).^2 + (Y(2:end) ...
    - Y(1:end-1)).^2))];
n_cont_points = round(cont_points(9)*CumLength(end)/(x_max-x_min));
Sn = linspace(CumLength(1),CumLength(end),n_cont_points);
X_topo = interp1(CumLength,X,Sn,'linear');
Y_topo = interp1(CumLength,Y,Sn,'linear');
cont_right = round((Y_topo(1)-Y_ulc(1))*cont_points(8)/(y_max4-y_max3));
cont_left = round((Y_topo(end)-Y_ulc(end))*cont_points(8)/(y_max4-y_max3));

% Calculate layer
LAYER4_2 = [x_max*ones(1,cont_right); ...
    linspace(Y_ulc(1),Y_topo(1),cont_right)];
LAYER4_2(:,[1 end]) = [];
Geo_id = [Geo_id 11*ones(1,size(LAYER4_2,2))];
LAYER4_3 = [X_topo; Y_topo];
Geo_id = [Geo_id 12*ones(1,size(LAYER4_3,2))];
LAYER4_4 = [x_min*ones(1,cont_left); ...
    linspace(Y_topo(end),Y_ulc(end),cont_left)];
LAYER4_4(:,[1 end]) = [];
Geo_id = [Geo_id 13*ones(1,size(LAYER4_4,2))];

GEOMETRY = [GEOMETRY LAYER4_2 LAYER4_3 LAYER4_4];
clear 'LAYER4_2' 'LAYER4_3' 'LAYER4_4'