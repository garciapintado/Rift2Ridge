function [GEOMETRY,Geo_id,Sgap,Intersect_ID,remesh] = ...
    dis_layers(GEOMETRY,Geo_id,Intersect_ID,shift,remesh)
% [GEOMETRY,GEO_ID,SGAP,INTERSECT_ID,REMESH] = DIS_LAYERS(GEOMETRY,GEO_ID,
% INTERSECT_ID,SHIFT,REMESH) evaluates where the interfaces separating 
% layers of a model geometry (defined by GEOMETRY and GEO_ID) are 
% intersecting and, therefore, discontinuous layers are required. It adds 
% the intersection points to GEOMETRY and GEO_ID and removes the points of 
% the interfaces that are above an upper interface. SGAP is a zero vector 
% which is 1 where a segment on the segmentlist variable in generate_mesh
% function is equivalent to a gap. Intersect_ID is a zero vector which is 1
% for points in GEOMETRY that are intersections between two horizontal 
% boundaries. SHIFT stablish the minimum distance between two interfaces, 
% if the distance is smaller than SHIFT the function removes this piece of
% the layer. When a new intersection is detected REMESH is set to 1 so
% remeshing will occur in the current time step.
%
% Example: Break-up of the crust
% --------
%
%          Exhumation of the mantle!
%  SURFACE   |
% _______    |   _______ Interface 2
%  CRUST \______/ CRUST
% _______/******\_______ Discontinuous interface 1 (moho)
% ******* MANTLE *******

%--------------------------------------------------------------------------
% Function written by Miguel Andres-Martinez, PhD student at Royal Holloway
% University of London, 16-07-2014. Email: mandresmartinez87@gmail.com
%--------------------------------------------------------------------------

% 12/08/2014 MA
    % Fixed problems related with points of an interface been part of
    % another interface above without being intersection points    
% 21/08/2014 MA 
    % Fixed problems related with single intersection point produce by gaps
    % removing the old intersections (edges of the parts)
% 26/08/2014 MA 
    % Fixed dis_layers4genmesh problems when more than one interface
    % intersect with an upper interface at the same point

%==========================================================================
% LAYERS AND INTERFACES FOR EVALUATION
%==========================================================================

% Example with two layers:
%  
%  ___________._Interface 3
% |  Layer 2  | 
% |___________|_Interface 2
% |  Layer 1  |
% |___________|_Interface 1

% Number of interfaces
n_interfaces = (max(Geo_id)-1)/3+1;
% Number of layers
n_layers = (max(Geo_id)-1)/3;

% Index the different combination of interfaces
Evaluation = [];

% Builds the matrix Evaluation of the size (combinations x 2) where all the
% possible combinations of pairs of interfering interfaces are included
%
% Example for two layers:
% Evaluation = 
%   1   2 => interface 2 interseting 1
%   2   3 => interface 3 intersecting 2
%   1   3 => interface 3 intersecting 1
for i = 1:n_interfaces-1
Evaluation = [((n_interfaces-i):-1:1)' ...
    (n_interfaces-i+1)*ones(n_interfaces-i,1); Evaluation];
end

% Interfaces for evaluation in terms of Geo_id
Eval_interf = [1 3:3:max(Geo_id)-1];
% Save coordinates of the interfaces for evaluation
X = GEOMETRY(1,ismember(Geo_id,Eval_interf));
Y = GEOMETRY(2,ismember(Geo_id,Eval_interf));
% Generates an id for the interfaces for evaluation
id = Geo_id(ismember(Geo_id,Eval_interf));
% Changes the id from Geo_id id to Evaluation id
for j = 1:n_interfaces
    id(id == Eval_interf(j)) = j;
end

% In case intersect_id is empty (previous step has not discontinuous
% layers), initializes intersect_id
if isempty(Intersect_ID)
    intersect_id = zeros(size(X));
else
    intersect_id = Intersect_ID(ismember(Geo_id,Eval_interf));
end

%==========================================================================
% LOOP THROUGH DIFFERENT COMBINATIONS OF INTERFACES
%==========================================================================
for le = 1:size(Evaluation,1);
    
    % Find the continuous pieces of the lower interface and assign a
    % differnt id to each piece, storing them at parts1 vector
    parts1 = index_parts(intersect_id(id==Evaluation(le,1)));
    
    % Load the lower interface in interface1
    interface1 = [X(id==Evaluation(le,1)); Y(id==Evaluation(le,1))];
    % Load the upper interface in interface2
    interface2 = [X(id==Evaluation(le,2)); Y(id==Evaluation(le,2))];
    
    % Saves the old intersect_id of the interface1 for later use
    intersect_id_old = intersect_id(id==Evaluation(le,1));
    
    % Initialize the final interfaces, ids and other book-keeping vectors
    interface1_f = [];
    interface2_f = [];
    intersect_id_f1 = [];
    intersect_id_f2 = [];
    
    interface2_add = [];
    int_indx2_add = [];
    n_prev_int = [];
    
    intersect_bool = 0;
    
    %======================================================================
    % LOOP THROUGH THE DIFFERENT PARTS OF INTERFACE1
    %======================================================================
    for lp1 = 1:max(parts1)
        
        % Load the coordinates of the lp1 part of the interface 1
        interface1p = interface1(:,parts1==lp1);
        
        %==================================================================
        % INTERSECTION POINTS
        %==================================================================
        
        % Calculate the SHIFT of interface1p to avoid too thin layers
        if shift~=0
            SHIFT = shiftline(interface1p(1,:),interface1p(2,:),shift);
        else
            SHIFT = zeros(size(interface1p));
        end
        
        % Find the points where interface2 intersects the lp1 part of
        % interface1 (int_points) and give them an index (indx_segment)
        [int_points, indx_segment] = ...
            intersect_lines(interface1p+SHIFT,interface2);
        
        % Remove intersections within gaps of interface2
        %-----------------------------------------------
        if ~isempty(int_points)
            % Intersect id of interface2
            intersect_id2 = intersect_id(id==Evaluation(le,2));
            % Segments of interface2 that are gaps
            seg_gap = ...
                find((intersect_id2(1:end-1) + intersect_id2(2:end)) == 2);
            % Fix mistakes due to contiguous gaps
            for jj = 1:size(seg_gap,2)
                if seg_gap(jj)==1 && seg_gap(jj-1)==1
                    seg_gap(jj) = 0;
                end
            end
            % Remove intersections of within the gaps in seg_gap
            rm_intersect = ismember(indx_segment(2,:),seg_gap);
            int_points(:,rm_intersect) = [];
            indx_segment(:,rm_intersect) = [];
        end
    
        % In case there is an intersection
        if ~isempty(int_points)
            % Set the boolean intersect_bool to 1 to indicate that
            % intersection happened
            intersect_bool = 1;
            
%             % Plots the intersection points (uncomment for visualization)
%             subplot(1,3,1)
%             hold on
%             plot(int_points(1,:),int_points(2,:),'.r')
%             for i = 1:size(indx_segment,2)
%                 plot([interface1p(1,indx_segment(1,i)) ...
%                     interface1p(1,indx_segment(1,i)+1)], ...
%                     [interface1p(2,indx_segment(1,i)) ...
%                     interface1p(2,indx_segment(1,i)+1)],'--r')
%                 plot([interface2(1,indx_segment(2,i)) ...
%                     interface2(1,indx_segment(2,i)+1)], ...
%                     [interface2(2,indx_segment(2,i)) ...
%                     interface2(2,indx_segment(2,i)+1)],'--r')
%             end
%             hold off
            
            %==============================================================
            % FIND POINTS OF INTERFACE1P ABOVE INTERFACE2 (GAP POINTS)
            %==============================================================
            
            % Index for the points on interface1p above interface2
            [gap] = points_above_l2(interface1p+SHIFT,interface2);
            
            % Prevent previous intersections from being gaps          
            gap(gap==1) = [];
            gap(gap==length(interface1p)) = [];
                        
%             % Plots the gap points (uncomment for visualization)
%             subplot(1,3,1)
%             hold on
%             plot(interface1p(1,gap),interface1p(2,gap),'o', ...
%                 'MarkerEdgeColor','k','MarkerFaceColor','w')
%             hold off
            
            %==============================================================
            % ADD INTERSECTION POINTS AND REMOVE GAP POINTS
            %==============================================================
            
            % Add intersection points to both interface1p and interface2
            % and remove gap points on interface1p
            [new_interface1,new_interface2,intersect_id_new1,int_indx2] ...
                = modify_lines(interface1p,interface2,int_points, ...
                indx_segment,gap);
            
            % Store the new interface1p into the final interface1
            interface1_f = [interface1_f new_interface1];
            % Load the old intersect_id for the current part lp1
            intersect_idp = intersect_id_old(parts1==lp1);
            % Remove the ids of the gaps in intersect_idp
            intersect_idp(gap) = [];
            % Store the indexes on intersct_idp into the zero values of
            % intersect_id_new1
            intersect_id_new1(intersect_id_new1==0) = intersect_idp;
            % Store the new ids (which contain the intersection points and
            % excludes the gap points at the current lp1 part of 
            % interface1) into the final intersect_id_f1
            intersect_id_f1 = [intersect_id_f1 intersect_id_new1];
            
            % Add the intersection points to the interface2_add
            interface2_add = [interface2_add int_points];
            % Indexes of the intersection points at interface2
            int_indx2_add = [int_indx2_add int_indx2];
            % Vector that accumulates number of previous intersections
            n_prev_int = [n_prev_int size(n_prev_int,2) ...
                *ones(size(int_indx2))];
            
%             % Plot the final part and interface2 (uncomment for
%             % visualization)
%             parts = [];
%             parts = index_parts(intersect_id_new1);
%             subplot(1,3,2)
%             hold on
%             for j = 1:max(parts)
%                 plot(new_interface1(1,parts==j), ...
%                     new_interface1(2,parts==j),'.r')
%                 plot(new_interface1(1,parts==j), ...
%                     new_interface1(2,parts==j),'r')
%             end
%             plot(new_interface2(1,:),new_interface2(2,:),'.k')
%             plot(new_interface2(1,:),new_interface2(2,:),'k')
            
        % In case there is no intersection between the lp1 part of
        % interface1 and interface2
        else
            % Store the not modified interface1p to the final interface1_f
            interface1_f = [interface1_f interface1p];
            % Store the not modified intersect_id of the lp1 part into the
            % final intersect_id_f1
            int_id = intersect_id(id==Evaluation(le,1));
            intersect_id_f1 = [intersect_id_f1 int_id(parts1==lp1)];
            
        end
        
    % End of the loop through the pieces (parts) of interface1
    end
    
    %======================================================================
    % UPDATING INTERFACES
    %======================================================================
    
    % In case there was any intersection between the parts of interface1
    % and interface2 updates coordinates and indexes
    if intersect_bool
        % Set remesh 1 so remeshing will occur in the current time step
        remesh = 1;
        
        % Create a NaN matrix of the size (2 x number of points in 
        % interface2) including the intersection points
        interface2_f = nan* ...
            ones(2,sum(id==Evaluation(le,2))+size(n_prev_int,2));
        % Add the intersections to the final interface2_f
        interface2_f(:,int_indx2_add+n_prev_int) = interface2_add;
        % Add the old coordinates into the interface2_f, in the right
        % position (where interface2_f is a NaN)
        interface2_f(:,isnan(interface2_f(1,:))) = ...
            [X(id==Evaluation(le,2)); Y(id==Evaluation(le,2))];
        
        % Create a NaN vector of the size (1 x number of points in 
        % interface2) including the intersection points
        intersect_id_f2 = nan* ...
            ones(1,sum(id==Evaluation(le,2))+size(n_prev_int,2));
        % Set intersect_id_f2 to 0 for the new intersections
        intersect_id_f2(int_indx2_add+n_prev_int) = 0;
        % Save the old intersection points in the right positions (where
        % intersect_id_f2 is a NaN)
        intersect_id_f2(isnan(intersect_id_f2)) = ...
            intersect_id(id==Evaluation(le,2));
        
        % Save the new x-coordinates into X ,in the correct position
        % respect the non evaluated interfaces at this le iteration
        X = [X(id<Evaluation(le,1)) interface1_f(1,:) ...
            X(id>Evaluation(le,1) & id<Evaluation(le,2)) ...
            interface2_f(1,:) X(id>Evaluation(le,2))];
        
        % Save the new y-coordinates into Y, in the correct position
        % respect the non evaluated interfaces at this le iteration
        Y = [Y(id<Evaluation(le,1)) interface1_f(2,:) ...
            Y(id>Evaluation(le,1) & id<Evaluation(le,2)) ...
            interface2_f(2,:) Y(id>Evaluation(le,2))];
        
        % Save the new id of the intersections into intersect_id, in the 
        % correct position respect the non evaluated interfaces at this le
        % iteration
        intersect_id = [intersect_id(id<Evaluation(le,1)) ...
            intersect_id_f1 intersect_id(id>Evaluation(le,1) & ...
            id<Evaluation(le,2)) intersect_id_f2 ...
            intersect_id(id>Evaluation(le,2))];
        
        % Save the new id for the modified interfaces 1 and 2
        id = [id(id<Evaluation(le,1)) ... 
            Evaluation(le,1)*ones(size(intersect_id_f1)) ...
            id(id>Evaluation(le,1) & id<Evaluation(le,2)) ...
            Evaluation(le,2)*ones(size(intersect_id_f2)) ...
            id(id>Evaluation(le,2))];
        
    end
end

%==========================================================================
% REMOVE FLOATING PARTS ABOVE UPPER INTERFACES
%==========================================================================

% Remove segments of interfaces that are above an upper interface, and 
% disconnected to any interface
[X,Y,id,intersect_id] = rm_floating_parts(X,Y,id,intersect_id);

%==========================================================================
% FIX DISCONECTED INTERSECTIONS
%==========================================================================
% Intersection points in the upper interface might be removed by the 
% shiftline function, leaving the lower-interface intersection without
% connection to any interface

[X,Y,intersect_id,id] = connect_dots(X,Y,intersect_id,id,Geo_id);

%==========================================================================
% BUILD GLOBAL INTERSECT_ID
%==========================================================================

% Build the global intersect_id, which is ordered as GEOMETRY and Geo_id

% Initialize the Intersect_ID saving the intersect_id of the bottom
% boundary
Intersect_ID = intersect_id(id==1);
% Loop through the number of layers of the model
for j = 1:n_layers
    % Save the intersect_id of the horizontal upper boundary of the layer
    % j, and zeros for the lateral boundaries, into the Intersect_ID vector
    Intersect_ID = [Intersect_ID zeros(1,sum(Geo_id==j*3-1)) ...
        intersect_id(id==j+1) ...
        zeros(1,sum(Geo_id==j*3+1))];
end

% Make Intersect_ID a boolean vector
Intersect_ID = Intersect_ID==1;

%==========================================================================
% GEOMETRY UPDATE
%==========================================================================

% Loop through the horizontal interfaces
for j = 1:n_interfaces
    % Update coordinates with the interstection points and without the gap
    % points
    GEOMETRY = [GEOMETRY(:,Geo_id<Eval_interf(j)) [X(id==j); Y(id==j)] ...
        GEOMETRY(:,Geo_id>Eval_interf(j))];
    % Update Geo_id taking into account new intersection points and gap
    % points
    Geo_id = [Geo_id(:,Geo_id<Eval_interf(j)) ...
        Eval_interf(j)*ones(1,size(X(id==j),2)) ...
        Geo_id(:,Geo_id>Eval_interf(j))];
end

%==========================================================================
% PLOT FINAL HORIZONTAL INTERFACES (uncomment for visualization)
%==========================================================================

% color_p = rand(n_interfaces,3);
% subplot(1,3,3)
% %figure(1)
% hold on
% for n = 1:n_interfaces
%     parts = index_parts(intersect_id(id==n));
%     xp = X(id==n);
%     yp = Y(id==n);
%     for m = 1:max(parts)
%         plot(xp(parts==m),yp(parts==m),'Color',color_p(n,:))
%         plot(xp(parts==m),yp(parts==m),'.','Color',color_p(n,:))
%     end
% end
% hold off

%==========================================================================
% REMOVE REPEATED POINTS
%==========================================================================
% Initialize rm_points_gl
rm_points_gl = [];

% Loop through the interfaces
for j = 1:n_interfaces
     % Evaluated interface for j
     INTERFACE = GEOMETRY(:,Geo_id==Eval_interf(j));
     % First point index
     first_indxp = find(Geo_id==Eval_interf(j),1,'first');
     % Find the repeated points
     [~,repeated_points] = ismember(GEOMETRY',INTERFACE','rows');
     % Remove the points of the interface Eval_interf(j) and the interfaces
     % below
     repeated_points(Geo_id<=Eval_interf(j)) = 0;
     % Sort the repeated points and delete the 0s
     repeated_points = sort(repeated_points(repeated_points~=0))';
     % Remove repeated points that have repeated points at the left and
     % right (points that are not intersection points)
     Mshift = [[Inf repeated_points(1:end-1)+1]; ...
         [repeated_points(2:end)-1 Inf]];
     rml = (Mshift(1,:)==repeated_points) + (Mshift(2,:)==repeated_points);
     rm_points = repeated_points(rml==2);
     % Global points to remove
     rm_points_gl = [rm_points_gl rm_points+(first_indxp-1)];
end

% Remove rm_points_gl
GEOMETRY(:,rm_points_gl) = [];
Geo_id(rm_points_gl) = [];
Intersect_ID(rm_points_gl) = [];

% Generates a new id for the interfaces for evaluation
id = Geo_id(ismember(Geo_id,Eval_interf));
% Changes the id from Geo_id id to Evaluation id
for j = 1:n_interfaces
    id(id == Eval_interf(j)) = j;
end

% Generates a new intersect_id
intersect_id = Intersect_ID(ismember(Geo_id,Eval_interf));

%==========================================================================
% FIX TRIANGLES AT THE INTERSECTIONS
%==========================================================================
% if remesh
%     [GEOMETRY,Geo_id,Intersect_ID] = ...
%         fix_tri_int(GEOMETRY,Geo_id,Intersect_ID);
%     % Generate an id for the interfaces for evaluation
%     id = Geo_id(ismember(Geo_id,Eval_interf));
%     % Save the Intersect_ID only for the horizontal boundaries
%     intersect_id = Intersect_ID(ismember(Geo_id,Eval_interf));
%     % Changes the id from Geo_id id to Evaluation id
%     for j = 1:size(Eval_interf,2)
%         id(id == Eval_interf(j)) = j;
%     end
% end

%==========================================================================
% BUILD SGAP
%==========================================================================
% Sgap indicates which segments of segmentlist in the function generate
% mesh correspond to gaps
[Sgap] = intid2gap(intersect_id,Geo_id,id,1);

 