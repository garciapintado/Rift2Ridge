function xyr = rotation(xy, phi)
  % +++ purpose +++
  % rotate points counter clockwise
  %
  %
  % xy  :: [p,2] matrix of coordinates
  % phi :: rotation angle, in radians. Positive angle is counterclockwise
  %
  % xyr:    [p,2] matrix of rotated coordinates
  
  % Author: Javier Garcia-Pintado. MARUM, 2020
  
  rotmat = [cos(phi), -sin(phi); ...                                       % x transform
            sin(phi), cos(phi)];                                           % y transform
  xyr = xy * rotmat';                                                      % [p,2]
end % end function rotation
