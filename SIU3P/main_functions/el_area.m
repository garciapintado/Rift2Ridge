function area = el_area(GCOORD,ELEM2NODE)
% AREA = EL_AREA(GCOORD,ELEM2NODE) calculates the area of elements for a
% given mesh defined by GCOORD and ELEM2NODE. Area of elements is useful
% for the calculation of Jacobians and determinant of the Jacobians and
% also because negative areas indicate inverted elements.

%--------------------------------------------------------------------------
% Function written by Miguel Andres-Martinez, Postdoc at University of
% Bremen. Email: mandresmartinez87@gmail.com
% Function based on MILAMIN, Dabrowiski et al. 2010.
%--------------------------------------------------------------------------

ECOORD_x = reshape(GCOORD(1,ELEM2NODE),size(ELEM2NODE,1),size(ELEM2NODE,2));
ECOORD_y = reshape(GCOORD(2,ELEM2NODE),size(ELEM2NODE,1),size(ELEM2NODE,2));

a23      = ECOORD_x(2,:).*ECOORD_y(3,:) - ECOORD_x(3,:).*ECOORD_y(2,:);
a31      = ECOORD_x(3,:).*ECOORD_y(1,:) - ECOORD_x(1,:).*ECOORD_y(3,:);
a12      = ECOORD_x(1,:).*ECOORD_y(2,:) - ECOORD_x(2,:).*ECOORD_y(1,:);
area     = a23 + a31 + a12;