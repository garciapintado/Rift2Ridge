function [Bc_ind, Bc_val, Bc_ind_fs, Bc_val_fs, ext_srate] = ...
    set_bcs_vel(GCOORD, Point_id, ext_rate, bc_t, dt)
  % set velocities 
  %
  % GCOORD    :: REAL [2,nnod] x,y coordinates of nodes
  % bc_t      :: CHARACTER, type of boundary condition
  % Point_id  :: INTEGER [1,nnodce], class with -1 for corner_id &
  %              cornin_id nodes, 1:10 for subdomain limits, and 0 for nodes
  %              interior to any domain. Note all remaining central nodes
  %              (not included in Point_id), would have a 0 if included in Point_id
  % ext_rate  :: REAL (m/s) horizontal full spread velocity
  % bv_t      :: CHARACTER, mechanical boundary condition type 

  % return::
  % Bc_ind    : INTEGER,[1,nbdy], indexes of the boundary condition location within GCOORD(:); ie.e taking into account x,y serialised indices within GCOORD  
  % Bc_val    : REAL,   [1,nbdy], velocity [m/s] values corresponding to Bc_ind
  % Point_id  : INTEGER [1,nnodce], as Point_id, but guaranteed to have a -1 value in the corner nodes
  % Bc_ind_fs : INTEGER,[1,nfsbdy], BC indices for free surface (subset of BC_ind)
  % Bc_val_fs : INTEGER,[1,nfsbdy], BC values for free surface (subset of BC_val)
  % ext_srate : REAL (1/s) horizontal strain rate as needed by the mechanical
  %
  % last version: Javier García-Pintado, MARUM, 2020-03 JGP
  %               local corner identification compliant with generate_meshGEO()

  corners_id = getCornerids(GCOORD, Point_id);                               % [2,nlayer] left,right corner for each layer from bottom to top 

  ext_srate =  ext_rate / diff(GCOORD(1, corners_id(:,1)));                  % [1/s]
  ext_srate_y = - ext_rate / (1 + ext_srate*dt);

  nlayers = length(corners_id);

  % CORNER VELOCITIES
  % bottom left and right
  Dyleft  = diff(GCOORD(2,corners_id(1,[1 nlayers]))); % domain thickness - left bound
  Dyright = diff(GCOORD(2,corners_id(2,[1 nlayers]))); % domain thickness - right bound - both should be the same

  %if Dyleft ~= Dyright
  %  error('Dyleft ~= Dyright: prepare set_bcs_flow for this case')
  %end

  il = 1;                                                                    % bottom interface
  switch bc_t                                                        
      case 'pureshear'
          Bc_ind = [2*(corners_id(1,il)-1)+1, 2*(corners_id(1,il)-1)+2, ...                           % bottom left  [i_x, i_y]  
                    2*(corners_id(2,il)-1)+1, 2*(corners_id(2,il)-1)+2];                              % bottom right [i_x, i_y]
          Bc_val = [ext_srate*GCOORD(1,corners_id(1,il)), -ext_srate_y*Dyleft, ...  % bottom left  [v_x, v_y]  [m/s]
                    ext_srate*GCOORD(1,corners_id(2,il)), -ext_srate_y*Dyright];    % bottom right [v_x, v_y]  [m/s]
      case {'winkler','winkler+LIVD'}                                                       
          Bc_ind = reshape(2*(corners_id(:,il)-1)+1,1,[]);                 % bottom left,right   [i_x]
          Bc_val = ext_srate*GCOORD(1,corners_id(:,il));                   % bottom left, right  [v_x] [m/s]
      otherwise
          error('Type of mechanical boundary conditions not recognised')
  end

  il = nlayers;                                                           % top interface
  Bc_ind = [Bc_ind reshape(2*(corners_id(:,il)-1)+1,1,[])];               % [i_x] left, right
  Bc_val = [Bc_val ext_srate*GCOORD(1,corners_id(:,il))];                 % [v_x] left, right  [m/s]

  % the following works but has been replacing by overall domain "drunken-sailor" compensation after the mechanical
  %if bc_t == "winkler+LIVD"
  %    Bc_ind = [Bc_ind reshape(2*(corners_id(:,il)-1)+2,1,[])];           % [i_y] left, right
  %    Bc_val = [Bc_val [0 0]];                                            % [v_y] left, right  [m/s]  
  %end

  %INNER CORNERS
  ils = 2:nlayers-1;
  Bc_ind = [Bc_ind reshape(2*(corners_id(:,ils)-1)+1,1,[])];              % [i_x] layer 3-left, layer 3-right, layer 6-left layer 6-right 
  Bc_val = [Bc_val ext_srate*GCOORD(1,corners_id(:,ils))];                % [v_x]    "               "             "              "
 
  % assure corner nodes are not considered for the following segments
  Point_id(corners_id) = -1;

  % SIDES
  if bc_t == "pureshear"
      ids     = find(Point_id==1);                                        % bottom line - not set for winkler
      Bc_ind  = [Bc_ind reshape(2*(ids-1)+2,1,[])];                       %              [i_y] vector                         
      Bc_val  = [Bc_val repelem(-ext_srate_y*Dyleft,1,length(ids))];      %              [v_y] vector
  end
  geomaxi = max(unique(Point_id));
  sideids = [4:3:geomaxi 2:3:geomaxi]; % left and right vertical bound indices

  ids = find(ismember(Point_id, sideids));                                % left [4,7,10...] & right [2,5,8...] vertical bound i indices in GCOORD[:,i] 
  Bc_ind  = [Bc_ind reshape(2*(ids-1)+1,1,[])];                           %              [i_x] vector
  Bc_val  = [Bc_val ext_srate*GCOORD(1,ids)];                             %              [v_x] vector

  Bc_ind_fs = Bc_ind;                                                     % free surface: everything but the upper bound set below
  Bc_val_fs = Bc_val;                                                     %   "    "         "        "   "    "     "

  ids     = find(Point_id==(geomaxi-1));                                  % top line
  Bc_ind  = [Bc_ind reshape(2*(ids-1)+2,1,[])];                           %              [i_y] vector
  Bc_val  = [Bc_val zeros(size(ids))];                                    %              [v_y] vector : fixed surface

  % Bc_tmp  = find(Point_id==6); % upper edge
  % Bc_ind  = [Bc_ind 2*(Bc_tmp-1)+1 2*(Bc_tmp-1)+2];
  % Bc_val  = [Bc_val ext_srate*abs(GCOORD(2,Bc_tmp)) zeros(1,length(Bc_tmp))]; 

end
