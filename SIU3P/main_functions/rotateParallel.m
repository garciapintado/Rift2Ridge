function [Xr, Yr] = rotateParallel(X, Y, phi)
  % +++ purpose +++
  % rotate points counter clockwise in groups of coordinates arranged in a
  % matrix
  % 
  %
  % x   :: [n,ncoo] matrix of x-coordinates, where n is the number of
  %        independenly rotated groups and ncoo in the number of coordinates in
  %        each group
  % y   :: [n,ncoo] matrix of y-coordinates
  % phi :: [n] rotation angle, in radians. Positive is counterclockwise.
  
  % Author: Javier Garcia-Pintado. MARUM, 2020
  
  ncoo = size(X,2);
  
  s = repmat(sin(phi(:)), 1, ncoo);
  c = repmat(cos(phi(:)), 1, ncoo);
  
  Xr = X .* c - Y .* s;
  Yr = X .* s + Y .* c;

end % end function rotateParallel
