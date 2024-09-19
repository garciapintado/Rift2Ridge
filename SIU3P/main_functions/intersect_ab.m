function [x, y] = intersect_ab(a1,b1,a2,b2)
% INTERSECT_AB(A1,B1,A2,B2) find the common point (X,Y) between the
% intersecting lines defined by A1, B1 and A2, B2.

x = (b2-b1)./(a1-a2);
y = a1*x + b1;

end

