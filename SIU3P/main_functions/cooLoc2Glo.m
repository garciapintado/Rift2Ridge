function [xy] = cooLoc2Glo(xyv, xieta)
  % transform local coordinates to global ones by simple linear transformation
  % canonical triangle assumed i,j,k -> (0,0),(0,1),(0,1)
  % input
  % xyv :  [2,3] vertices of triangle in physical (global) space
  % xieta: [2,np] local coordinates, points to be transformed

  % return:
  % xy:  [2,np] global coordinates
  
  np = size(xieta,2);
                                                      % Jacobian
  Jf = [xyv(1,2) - xyv(1,1), xyv(1,3) - xyv(1,1);     % x_xi x_eta 
        xyv(2,2) - xyv(2,1), xyv(2,3) - xyv(2,1)];    % y_xi y_eta
  xy = repmat(xyv(:,1),1,np) + Jf * xieta;

end


