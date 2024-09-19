function [XX,YY,FIELD] = int2meshgrid(dir_p,reg_mesh_resol,field)
% [XX,YY,FIELD,TOPO] = INT2MESHGRID(DIR_P,REG_MESH_RESOL,FIELD) 
% interpolates a given FIELD variable defined at the integration points 
% (string) from the integration points to a regular grid which spacing is 
% given by reg_mesh_resol. The values of the FIELD are taken from the file
% given by DIR_P.

% Load the necessary variables from the file
load(dir_p,'GCOORD','ELEM2NODE','Point_id','nel','nip')
ori_field = field;
switch field
    case 'I2.p'
        FIELD_u = load(dir_p,'I2'); a = fieldnames(FIELD_u); 
        FIELD_u = FIELD_u.(a{1});
        FIELD_u = FIELD_u.p;
    otherwise
        FIELD_u = load(dir_p,field); a = fieldnames(FIELD_u); 
        FIELD_u = FIELD_u.(a{1});
end

% Find limits of the regular mesh
minx = min(GCOORD(1,:));
maxx = max(GCOORD(1,:));
miny = min(GCOORD(2,:));
maxy = max(GCOORD(2,:));

% Build regular mesh
X = minx:reg_mesh_resol:maxx;
Y = miny:reg_mesh_resol:maxy;
[XX,YY] = meshgrid(X,Y);

% Integration point calculation
[GIP_xx,GIP_yy] = int_points(GCOORD,ELEM2NODE,nip);

% FIELD function
F_FIELD = scatteredInterpolant(double(GIP_xx(:)),double(GIP_yy(:)), ...
    double(FIELD_u(:)));

% Interpolation of FIELD function to the regular mesh
FIELD = F_FIELD(XX,YY);