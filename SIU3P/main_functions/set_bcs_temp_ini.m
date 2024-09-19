function [Bc_ind, Bc_val, Point_id] = set_bcs_temp_ini(...
    GCOORD, Point_id, temp_bc, temp_bc_depth, ...
    temp_surf, ini_deformation, sigma_tx, sigma_tz)

% GCOORD          :: REAL [2,nnod]
% Point_id        :: INT  [nnod6], where nnod6 does not include the central node at each element
% temp_bc         :: REAL [1] (ºC)
% temp_bc_depth   :: REAL  (m)
% temp_surf       :: REAL [1] (ºC)
% ini_deformation :: INT  [1], flag for deformation type 
% sigma_tx        :: REAL [1] (m) sd of the Gaussian window for temperature decay - x    
% sigma_tz        :: REAL [1] (m) sd of the Gaussian window for temperature decay - z    

  ids = find(Point_id == (max(Point_id) - 1));                               % top
  Bc_ind = ids;
  Bc_val = repmat(temp_surf,1,length(ids));

  switch ini_deformation
    case {0,1,2,3,5,6}
        ids = find(GCOORD(2,:) <= temp_bc_depth);                          % nodes below boundary condition location
    case 4
        xmin = min(GCOORD(1,:));
        xmax = max(GCOORD(1,:));
        xseq = linspace(xmin,xmax,1000);
        xc   = 0.5 * (xmin + xmax);
        yisotherm = temp_bc_depth + sigma_tz*exp(-(xseq-xc).^2./((sigma_tx/2)^2)); % Gaussian shape
        ids = find(GCOORD(2,:) <= interp1(xseq,yisotherm,GCOORD(1,:)));
  end
  ids = unique([ids find(Point_id == 1)]);
  Bc_ind = [Bc_ind ids];
  Bc_val = [Bc_val repmat(temp_bc,1,length(ids))];

