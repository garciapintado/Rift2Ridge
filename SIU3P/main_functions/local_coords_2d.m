function LCOORD = local_coords_2d(GCOORD, EL2NOD, els, xy)
% Usage: LCOORD = local_coords_2d(GCOORD, EL2NOD, els, xy)
%
% Purpose: Returns local coordinates of points "xy" in elements "els".
%
% Input:
%   GCOORD : [matrix]    : coordinates of all nodes in mesh
%   EL2NOD : [matrix]    : finite element connectivity matrix (nnodel x nel)
%   els    : [rowvector] : element in which each point is located
%   xy     : [2,n]    : coordinates of points to be located
% Output:
%   LCOORD : [matrix]    : local coordinates of points in each element


% JH Jan 2011
% JH Jan 2013 : handles NaN in els
% Javier GP, MARUM, 2020: error cheking. Reasonable rounding to prevent matlab
%                                        numerical errors while allocating

if size(xy,2) ~= length(els)
    error("local_coords_2d:: size(xy,2) ~= length(els)")
end

ind = ~isnan(els) & els>0;
els = els(ind);
ns  = length(els);
x   = reshape(GCOORD(1,EL2NOD(1:3,els)),3,ns);
y   = reshape(GCOORD(2,EL2NOD(1:3,els)),3,ns);

xp  = xy(1,ind);
yp  = xy(2,ind);

LCOORD = nan(2,length(ind));
LCOORD(1,ind) = -(-x(1,:).*yp+x(1,:).*y(3,:)-x(3,:).*y(1,:)+xp.*y(1,:)+x(3,:).*yp-xp.*y(3,:))./ ...
                 (-x(1,:).*y(3,:)+x(1,:).*y(2,:)-x(2,:).*y(1,:)+x(2,:).*y(3,:)+x(3,:).*y(1,:)-x(3,:).*y(2,:));
LCOORD(2,ind) = (x(1,:).*y(2,:)-x(1,:).*yp-x(2,:).*y(1,:)+x(2,:).*yp+xp.*y(1,:)-xp.*y(2,:))./ ...
               (-x(1,:).*y(3,:)+x(1,:).*y(2,:)-x(2,:).*y(1,:)+x(2,:).*y(3,:)+x(3,:).*y(1,:)-x(3,:).*y(2,:));

LCOORD = round(LCOORD,11);

end % END OF FUNCTION local_coords_2d