function diffcoo = diff2D(X)
%  
% obtain segment length for every line along a polyline
%
% X :: REAL [2,ncoo]

% Autor: Javier GP, MARUM, 2020
%--------------------------------------------------------------------------
    diffcoo = sqrt(diff(X(1,:)).^2 + diff(X(2,:)).^2);
end