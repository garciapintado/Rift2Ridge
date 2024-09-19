function [GEOMETRY, Geo_id, Intersect_ID,Sgap]=  resample_interf(GCOORD,GEOMETRY,Geo_id, ...
    Corner_id, Cornin_id, cont_points, x_max, x_min, Ylim, Intersect_ID)

% GCOORD:    REAL [2,nnodes]  x,z coordinates for all nodes in the grid
% GEOMETRY:  IO, REAL [2,ngeo]    x,z coordinates defining domain and subdomain boundaries
% Geo_id  :  IO, INTEGER [1,ngeo] class of nodes in GEOMETRY
% Corner_id: [ll,lr,ur,ul]
% Cornin_id:
% cont_points: [7] Number of points in each geometry interface, see main.m for definition. These go from bottom to top an horizontal, vertical, alternatively

% JGP: on-going adaptation to avoid previous overshooting of the surface
% node, which makes the odd crash!


% Save old geometries
GEOMETRY_OLD = GEOMETRY;
Geo_id_old = Geo_id;
Intersect_ID_old = Intersect_ID;

% figure(); plot_meshF(ELEM2NODE, GCOORD);
% plot_geomF(GEOMETRY, Geo_id)
% hold on; plot_geomF(GEOMETRY1, Geo_id1); plot_geomF(GEOMETRY, Geo_id, 6, '+'); plot_geomF(GEOMETRY, Geo_id, 3, 'o')
% plot_geomF(GEOMETRY0, Geo_id0, 3, '+')


% Build Id vector for the loop of horizontal interfaces
Id = [1 3:3:max(Geo_id)];

%==========================================================================
% RESAMPLE INTERFACES
%==========================================================================
% Loop for resampling horizontal interfaces
for i = 1:length(Id)
    gid = Id(i);
    % Loads the coordinates of the i interface 
    INTERFACE     = GEOMETRY(:,Geo_id==gid);
    INTERFACE_OLD = GEOMETRY_OLD(:,Geo_id_old==gid);
    
    % Intersect_ID_old of the points below interface i
    IntIDb = Intersect_ID_old & Geo_id_old < gid;
    
    % Find parts of INTERFACE
    parts = index_parts(Intersect_ID(Geo_id==gid));                        % only one part for layyer collapse
    
    X_new = [];
    Y_new = [];
    
    % Loop through the parts of INTERFACE
    for j = 1:max(parts)
        % Save the interface part j
        INTERFACE_p = INTERFACE(:,parts==j);

        % Find intersections of lower interfaces into this part
        [Intersect_below, Intersect_below_id] = ...
            ismember(GEOMETRY_OLD(:,IntIDb)',INTERFACE_p', ...
            'rows');
        % Reorder Intersect_below_id
        Intersect_below_id = sort(Intersect_below_id);
        
        % Index the points that divide the part j into smaller parts 
        % between Intersect_below_id
        parts_parts = [1 Intersect_below_id' size(INTERFACE_p,2)];
        parts_parts(parts_parts==0) = [];
        
        % Loop for the parts between the points where lower interfaces
        % intersect with the part j
        for k = 1:size(parts_parts,2)-1
            % Save the interface part jk
            INTERFACE_pp = INTERFACE_p(:,parts_parts(k):parts_parts(k+1));
            
            % Cumulative length of the interface
            INTERFACE_pp(3,:) = [0 cumsum(sqrt((INTERFACE_pp(1,2:end) ...
                - INTERFACE_pp(1,1:end-1)).^2 + (INTERFACE_pp(2,2:end) ...
                - INTERFACE_pp(2,1:end-1)).^2))];
            
            % 0
            s_start = INTERFACE_pp(3,1);
            % Total length of the part jk
            s_end   = INTERFACE_pp(3,end);
            % Number of cont points for the part jk, for the same initial
            % resolution for the new length
            n_cont_points = round(cont_points(i*2-1)*s_end/(x_max-x_min));
            % Cumulative distance of the new points
            Sn      = linspace(s_start, s_end, n_cont_points);
            % In case the part is smaller that the resolution defined by
            % cont points set Sn as the initial point and the last point
            if length(Sn)<2
                Sn = [s_start s_end];
            end
            
            if size(INTERFACE_pp,2)>1
                % Calculation of the new coordinates
                X_new_p = interp1(INTERFACE_pp(3,:), INTERFACE_pp(1,:), ...
                    Sn, 'spline');
                Y_new_p = interp1(INTERFACE_pp(3,:), INTERFACE_pp(2,:), ...
                    Sn, 'spline');
            else
                X_new_p = INTERFACE_pp(1,:);
                Y_new_p = INTERFACE_pp(2,:);
            end
            
            % Saves new coordinates, except the last point which would be
            % added in the next iteration
            X_new = [X_new X_new_p(1:end-1)];
            Y_new = [Y_new Y_new_p(1:end-1)];
            
            % Adds the last point in case it is the last part jk of part j
            if k==sum(Intersect_below)+1
                X_new = [X_new X_new_p(end)];
                Y_new = [Y_new Y_new_p(end)];
            end

         end
        
%         if ~isempty(Inter_below)
%             % Loop for the Inter_below points in the evaluated part j
%             for k = 1:size(Inter_below,1)
%                 xo = INTERFACE_OLD_p(1,Inter_below(k));
%                 yo = INTERFACE_OLD_p(2,Inter_below(k));
%                 [a1,b1] = line_ab(X_new(1:end-1),Y_new(1:end-1), ...
%                     X_new(2:end),Y_new(2:end));
%                 b2 = yo+1./a1*xo;
%                 xi = (b2-b1)./(a1+1./a1);
%                 yi = a1.*xi+b1;
%                 inter_x = (X_new(1:end-1)>=xi & X_new(2:end)<=xi) | ...
%                     (X_new(1:end-1)<=xi & X_new(2:end)>=xi);
%                 inter_y = (Y_new(1:end-1)>=yi & Y_new(2:end)<=yi) | ...
%                     (Y_new(1:end-1)<=yi & Y_new(2:end)>=yi);
%                 ind_inter = [find(inter_x & inter_y) find(inter_x & inter_y)+1];
%                 X_new = [X_new(1:ind_inter(1)) INTERFACE_OLD_p(1,Inter_below(k)) X_new(ind_inter(2):end)];
%                 Y_new = [Y_new(1:ind_inter(1)) INTERFACE_OLD_p(2,Inter_below(k)) Y_new(ind_inter(2):end)];
%             end
%         end
        
    end % for j
    
    % Saving corners
    if Id(i) == 1
        X_new(1)    = GCOORD(1,Corner_id(1));                              % ll x
        X_new(end)  = GCOORD(1,Corner_id(2));                              % lr x
        Y_new(1)    = GCOORD(2,Corner_id(1));                              % ll y
        Y_new(end)  = GCOORD(2,Corner_id(2));                              % lr y
    elseif Id(i) == max(Id)
        X_new(1)    = GCOORD(1,Corner_id(3));                              % ur x  
        X_new(end)  = GCOORD(1,Corner_id(4));                              % ul x
        Y_new(1)    = GCOORD(2,Corner_id(3));                              % ur y 
        Y_new(end)  = GCOORD(2,Corner_id(4));                              % ur y
    else
        X_new(1)    = GCOORD(1,Cornin_id((i-1)*2));                        % r x
        X_new(end)  = GCOORD(1,Cornin_id((i-1)*2-1));                      % l x
        Y_new(1)    = GCOORD(2,Cornin_id((i-1)*2));                        % r y
        Y_new(end)  = GCOORD(2,Cornin_id((i-1)*2-1));                      % l y 
    end
    
    % Coordinates of the previous interface
    GEO_LEFT                = GEOMETRY(:,Geo_id<Id(i));
    % Coordinates of the next interface
    GEO_RIGHT               = GEOMETRY(:,Geo_id>Id(i));
    % Id of the previous interface
    Gid_left                = Geo_id(Geo_id<Id(i));
    % Id of the next interface
    Gid_right               = Geo_id(Geo_id>Id(i));
    
    % Find the new Intersect_ID for the resampled interface
    Intersect_ID_new = ismember(floor([X_new; Y_new]*100)', floor( ...
        GEOMETRY_OLD(:,Intersect_ID_old & Geo_id_old==Id(i))*100)','rows');
    % Store the new Intersect_ID
    Intersect_ID = [Intersect_ID(Geo_id<Id(i)) Intersect_ID_new' ...
        Intersect_ID(Geo_id>Id(i))];
    
    % Save new interface
    GEOMETRY                = [GEO_LEFT [X_new; Y_new] GEO_RIGHT];
    Geo_id                  = [Gid_left Id(i)*ones(1,size(X_new,2)) Gid_right];
    
end % for resampling horizontal interfaces

GEOMETRY(:,Intersect_ID) = GEOMETRY_OLD(:,Intersect_ID_old);

%==========================================================================
% BUILD NEW SGAP
%==========================================================================

% Interfaces for evaluation in terms of Geo_id
Eval_interf = [1 3:3:max(Geo_id)-1];
% Generate an id for the interfaces for evaluation
id = Geo_id(ismember(Geo_id,Eval_interf));
% Save the Intersect_ID only for the horizontal boundaries
intersect_id = Intersect_ID(ismember(Geo_id,Eval_interf));
% Changes the id from Geo_id id to Evaluation id
for j = 1:size(Eval_interf,2)
    id(id == Eval_interf(j)) = j;
end
% Calculate Sgap
[Sgap] = intid2gap(intersect_id,Geo_id,id,1);

%==========================================================================
% RESAMPLE LATERAL BOUNDARIES
%==========================================================================
% Find lateral ids
Id_lat = 1:max(Geo_id);
Id_lat(ismember(Id_lat,Id)) = [];

% Geometry indexes for the last and first points of each lateral boundary
Lim_id = repmat(1:(max(Geo_id)-1)/3,2,1);
Lim_id = Lim_id(:);

% Inner interfaces
In_id = Id(2:end-1);

% Loop through lateral boundaries
for n = 1:length(Id_lat)
    % Find corners indexes to include them as beginning and ends of the 
    % lateral boundaries on the resampling
    Cl1 = [find(Geo_id==1,1,'last') find(Geo_id==3,1,'last')];             % upward 
    Cl2 = [find(Geo_id==3,1,'first') find(Geo_id==1,1,'first')];           %
    % JGP commented
    for m = 1:length(In_id)
        Cl1 = [Cl1 find(Geo_id==In_id(m),1,'first') ...                    % 
                   find(Geo_id==In_id(m)+3,1,'last')];
        Cl2 = [Cl2 find(Geo_id==In_id(m)+3,1,'first') ...
            find(Geo_id==In_id(m),1,'last')];
    end
    % Indexes of the geometry vectors to remove the old boundary and add
    % the new resampled boundary
    Geo_left = Geo_id<Id_lat(n);
    Geo_right = Geo_id>Id_lat(n);
    
    % Load coordinates of the boundary
    X = GEOMETRY(1,Geo_id==Id_lat(n));
    X = [GEOMETRY(1,Cl1(n)) X GEOMETRY(1,Cl2(n))];
    Y = GEOMETRY(2,Geo_id==Id_lat(n));
    Y = [GEOMETRY(2,Cl1(n)) Y GEOMETRY(2,Cl2(n))];
    
    % Cumulative length of the boundary
    Cum_lat_length = [0 cumsum(sqrt((X(2:end) ...
        - X(1:end-1)).^2 + (Y(2:end) ...
        - Y(1:end-1)).^2))];
    
    s_start = Cum_lat_length(1);
    % Total length of the part jk
    s_end   = Cum_lat_length(end);
    % Number of cont points for the part jk, for the same initial
    % resolution for the new length
    n_cont_points = round(cont_points(Lim_id(n)*2)*s_end/(Ylim(Lim_id(n)+1)-Ylim(Lim_id(n))));
    % Cumulative distance of the new points
    Sn      = linspace(s_start, s_end, n_cont_points);
    % In case the part is smaller that the resolution defined by
    % cont points set Sn as the initial point and the last point
    if length(Sn)<2
        Sn = [s_start s_end];
    end
    
    if size(X,2)>1
        % Calculation of the new coordinates
        X_new = interp1(Cum_lat_length,X,Sn,'spline');
        Y_new = interp1(Cum_lat_length,Y,Sn,'spline');
    else
        X_new = X(2:end-1);
        Y_new = Y(2:end-1);
    end
        X_new = X_new(2:end-1);
        Y_new = Y_new(2:end-1);

    GEOMETRY = [GEOMETRY(1,Geo_left) X_new GEOMETRY(1,Geo_right); ...
        GEOMETRY(2,Geo_left) Y_new GEOMETRY(2,Geo_right)];
    
    Geo_id = [Geo_id(Geo_left) Id_lat(n)*ones(1,n_cont_points-2) Geo_id(Geo_right)];
end
