function [Bc_ind, Bc_val, Point_id] = set_bcs_temp(GCOORD, Point_id, temp_bc, temp_bc_depth, temp_surf)
                                                             
  ids = find(Point_id == (max(Point_id) - 1));                             % top
  Bc_ind = ids;
  Bc_val = repmat(temp_surf,1,length(ids));

  ids = unique([find(GCOORD(2,:) <= temp_bc_depth) ...                     % bottom or nodes below boundary condition
                find(Point_id == 1)]);      
  Bc_ind = [Bc_ind ids];
  Bc_val = [Bc_val repmat(temp_bc,1,length(ids))];
end

