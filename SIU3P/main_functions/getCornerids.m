function corners_id = getCornerids(GCOORD, Point_id)
% return a [2,n] matrix of indices in GCOORD, with n being the number of pseudo-horizontal layers in the
% domain, representing the corners (from bottom to top).
%
% Author: Javier García-Pintado, MARUM, 2020-03


  gids = unique(Point_id);
  hids = [1 3:3:max(gids)-1];
  corners_id = zeros(2,length(hids));
  for i=1:length(hids)
      hid = hids(i);
      locids = find(Point_id == hid);
      [~,minlid] = min(GCOORD(1,locids));
      [~,maxlid] = max(GCOORD(1,locids));
      corners_id(:,i) = locids([minlid,maxlid]);
  end 

end % end function corners_id

