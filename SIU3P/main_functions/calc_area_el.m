function [area_el,len_el,circumcirc_el,incirc_el] = calc_area_el(GCOORD,EL2NOD)
%
% Purpose: Calculate element sizes
% Input:
%   GCOORD   :: nodal coordinates
%   EL2NOD   :: connectivity matrix
% Output:
%   area_el   :: area of each element

% JH Jan 2011
%

area_el = abs(0.5*( GCOORD(1,EL2NOD(2,:)).*GCOORD(2,EL2NOD(3,:))...
                   -GCOORD(1,EL2NOD(3,:)).*GCOORD(2,EL2NOD(2,:))...
                   -GCOORD(1,EL2NOD(1,:)).*GCOORD(2,EL2NOD(3,:))...
                   +GCOORD(1,EL2NOD(1,:)).*GCOORD(2,EL2NOD(2,:))...
                   +GCOORD(1,EL2NOD(3,:)).*GCOORD(2,EL2NOD(1,:))...
                   -GCOORD(1,EL2NOD(2,:)).*GCOORD(2,EL2NOD(1,:)))   )';

nel           = size(EL2NOD,2);
x_el_vertices = reshape(GCOORD(1,EL2NOD(1:3,:))',3,nel);
y_el_vertices = reshape(GCOORD(2,EL2NOD(1:3,:))',3,nel);
a             = sqrt( (x_el_vertices(1,:)-x_el_vertices(2,:)).^2 + ...
                      (y_el_vertices(1,:)-y_el_vertices(2,:)).^2 )';
b             = sqrt( (x_el_vertices(2,:)-x_el_vertices(3,:)).^2 + ...
                      (y_el_vertices(2,:)-y_el_vertices(3,:)).^2 )';
c             = sqrt( (x_el_vertices(3,:)-x_el_vertices(1,:)).^2 + ...
                      (y_el_vertices(3,:)-y_el_vertices(1,:)).^2 )';

% len_el = min([a b c],[],2);
len_el = mean([a b c],2);

% Diameter of the circumcircle
circumcirc_el = a.*b.*c ./ (2*area_el);

% Diameter of the incircle
incirc_el = 4.*area_el./(a+b+c);

end % END OF SUBFUNCTION calc_area_el
