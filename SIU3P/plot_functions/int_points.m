function [GIPx,GIPy] = int_points(GCOORD,ELEM2NODE,nip)
% [GIPX,GIPY] = INT_POINTS(GCOORD,ELEM2NODE,NIP) calculates the coordinates of
% the NIP number of integration points associated to the triangular
% elements given by the mesh defined by the coordinates GCOORD and the
% conectivity matrix ELEM2NODE.

% Local positions of the integration points
[IP_X,~] = ip_triangle(nip);
% Shape functions
[N,~] = shp_deriv_triangle(IP_X,size(ELEM2NODE,1));

nel = size(ELEM2NODE,2);

% Initialize integration points
GIPx = zeros(size(ELEM2NODE(1:6,:)))';
GIPy = zeros(size(ELEM2NODE(1:6,:)))';

% Coordinates organized by elements
ECOORD_x = reshape(GCOORD(1,ELEM2NODE), size(ELEM2NODE,1), nel);
ECOORD_y = reshape(GCOORD(2,ELEM2NODE), size(ELEM2NODE,1), nel);

% Integration point iteration
for ip=1:nip
    % Load shape function
    Ni      =        N{ip};
    
    % Global coordinatess of the integration points
    GIP_x   = Ni'*ECOORD_x;
    GIP_y   = Ni'*ECOORD_y;
    
    GIPx(:,ip)   = GIP_x;
    GIPy(:,ip)   = GIP_y;
end