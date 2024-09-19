function usegments = findUnderlevelSegments(Topography, sea_level)
   % function usegments = findUnderlevelSegments(Topography, sea_level)
   % 
   % +++ purpose +++
   % obtain a matrix of segment bounds for topography below a specific level (one segment per row)
   % The segments are obtained by linear interpolation so the output coordinates do not need
   % to match any specific node in the Topography array. The reference level can be a scalar (e.g. sea level) or a matrix of coordinates for a varying level
   %
   % INPUT
   % Topography :: [2,ntopo] x-monotically increasing, array of coordinates
   % sea_level  :: scalar, representing a constant level, or [2,nsea]
   %               matrix or coordinates.
   %
   % OUTPUT
   % uwsegments :: REAL [nseg,2], where each row contains the semgnet
   % bound coordinates (right, left).
   
   % Javier García_Pintado, MARUM, 2020
   
   if length(sea_level) == 1                                               % constant
     z = Topography(2,:) - sea_level;
   else
     sea_level = interp1(sea_level(1,:), sea_level(2,:), Topography(1,:), 'linear', 'extrap');  
     z = Topography(2,:) - sea_level;
   end
   if  all(z >= 0.) || all(z <= 0.)
     xo = [];
   else    
     xo = findZeroes(Topography(1,:), z); 
   end
   
   if z(1) < 0.
     xo = [Topography(1,1); xo];    
   end    
   if z(end) < 0.
     xo = [xo; Topography(1,end)];    
   end 
   usegments = reshape(xo,2,[])';                                          % underwater segment bounds
   
end % function
