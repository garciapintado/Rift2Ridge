function [a,b] = line_ab(x1,y1,x2,y2)
% LINE_AB(X1,Y1,X2,Y2) returns the parameters a b for y = ax + b line
% deduced from (X1, Y1) and (X2, Y2) points.

a = (y2-y1)./(x2-x1);
b = y1-a.*x1;

end

