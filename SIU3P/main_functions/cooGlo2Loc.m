function [xieta] = cooGlo2Loc(xyv, xy)
  % transform global coordinates to local ones by simple linear transformation
  % canonical triangle assumed i,j,k -> (0,0),(0,1),(0,1)
  % input
  % xyv : [2,3] vertices of triangle in physical (global) space
  % xy:   [2,np] global coordinates

  % return:
  % xieta: [2,np] local coordinates, points to be transformed  

  np = size(xy,2);
  idetJ = 1 / (  (xyv(1,2) - xyv(1,1)) * (xyv(2,3) - xyv(2,1)) ...
               -((xyv(2,2) - xyv(2,1)) * (xyv(1,3) - xyv(1,1))));
                                                                           % inverse of Jacobian
  Ji = idetJ * [xyv(2,3) - xyv(2,1), xyv(1,1) - xyv(1,3); ...              % x_xi x_eta 
                xyv(2,1) - xyv(2,2), xyv(1,2) - xyv(1,1)];                 % y_xi y_eta
  xieta = Ji * (xy - repmat(xyv(:,1),1,np));

  % example
  if 1 > 2
    xyv = [5, 12, 7; 3, 1, 6];
    xy  = [6, 8; 4, 4];
    xieta = cooGlo2Loc(xyv, [xy, xyv]); 
    xya   = cooLoc2Glo(xyv, xieta);                               % == [xy, xyv]  
 
   %figure(); % plot local
    %plot(nvx,nvy,'.','markersize',20);
    %hold on;
    %plot(lx, ly, '.','markersize',20,'color','red');
    %plot(lx(1), ly(1), '.','markersize',10,'color','cyan');
    %plot(lx(2), ly(2), '.','markersize',10,'color','yellow');

    %figure() % plot global
    %plot(pvx,pvy,'.','markersize',20);
    %hold on;
    %plot(gxy(1,:), gxy(2,:), '.','markersize',20,'color','red');
    %plot(gxy(1,1), gxy(2,1), '.','markersize',10,'color','cyan');
    %plot(gxy(1,2), gxy(2,2), '.','markersize',10,'color','yellow');
  
  
  
  end % example  
end


