function [sxx, syy, txy] = rotateSymmetricTensor2D(sigma_x, sigma_y, tau, phi, jaumann)
  % +++ purpose +++
  % transformation of 2D tensor by rotation of axes. This results from the
  % operation Q*S*Q^t, where Q = | cos(phi) sin(phi)|, Q^ is its tranpose, and S a 2D stress tensor
  %                              |-sin(phi) cos(phi)| 
  % 
  % sigma_xx ::  tensor normal component on x
  % sigma_yy ::  tensor normal component on y
  % tau_xy   ::  tensor diagonal value (e.g. shear for stress)
  % phi         :: [radians] angle to conduct counterclockwise rotation (positive angle)
  % jaumann     :: LOGICAL true to conduct the approximate Jaumann rotations
  %
  % Reference: R.M. Jones, Mechanics of Composite Materials, McGraw-Hill, 1975. 
  %
  % Details: phi indicates the angle (positive counterclokwise) of the new axes with respect to the
  % current axes. Thus, if the situation is that a block of material has rotated (let's say an angle phi>0; i.e. counterclockwise) and we use this
  % function to express its tensor state with respect to original axes, the
  % input angle phi has to be expressed as a negative value, as the
  % original axes are at -phi with respect to the block axes after the block has moved.
  % in general it is not adviced to use Jaumann rotation, unless it is
  % guaranteed that rotation angles are vey small.
  
  % Author: Javier Garcia-Pintado. MARUM, 2020
  
  if nargin < 5
      jaumann = false;
  end
  
  s = sin(phi); 
  c = cos(phi);
  
  if jaumann		                                                       % Jaumann' simplified co-rotation formulas
    sxx = sigma_x + 2 * phi .* tau;                                    
    syy = sigma_y - 2 * phi .* tau;
    txy = tau - phi .* (sigma_x - sigma_y);                                % [nel, nip_stress]                      "                   "                                     Eq. (12.26) 
  else
    sxx = sigma_x .* c.^2 + sigma_y .* s.^2 + 2.* tau .*s.*c;              % sin(2*phi) == 2*sin(phi)*cos(phi)
    syy = sigma_x .* s.^2 + sigma_y .* c.^2 - 2.* tau .*s.*c;
    txy = (sigma_y - sigma_x) .* s.*c + tau .* (c.^2 - s.^2);              % cos(2*phi) == cos(phi)^2 - sin(phi)^2
  end
  
end % end function rotateSymmetricTensor2D()
