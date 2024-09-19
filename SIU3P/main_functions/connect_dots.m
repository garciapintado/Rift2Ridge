function [X,Y,intersect_id,id] = connect_dots(X,Y,intersect_id,id,Geo_id)
% [X,Y,INTERSECT_ID,ID] = CONNECT_DOTS(X,Y,INTERSECT_ID,ID,GEO_ID) locates
% the disconnected intersections points INTERSECT_ID, in the interfaces 
% defined by X,Y,ID and GEO_ID, and connect them to the neighbour point of
% an upper interface.

%--------------------------------------------------------------------------
% Function written by Miguel Andres-Martinez, PhD student at Royal Holloway
% University of London, 16-07-2014. Email: mandresmartinez87@gmail.com
%--------------------------------------------------------------------------

% Intersection points
IntersectP = [X(intersect_id==1);Y(intersect_id==1)];
% Intersection-points id
Int_id = id(intersect_id==1);

% Find disconected IntersectP
[~,Int_ind] = ismember([X;Y]',IntersectP','rows');
% Find how many times the intersections points are repeated
Int_rep = hist(Int_ind(Int_ind~=0),unique(Int_ind(Int_ind~=0)));
% Disconnected intersection points
Int_disc = find(Int_rep==1);

if ~isempty(Int_disc)
    
    % Find Sgap for the evaluated interfaces
    Sgap_int = intid2gap(intersect_id,Geo_id,id,0);
    
    % Initializing values for the loop
    Neighbour = [];
    Neighbour_pos = [];
    Ind_dist_global = [];
    % Loop throught the disconnected intersection points
    for k = 1:size(Int_disc,2)
        % Distance between all the points of the interfaces and the 
        % Int_disc(k)
        Dist_int = sqrt((IntersectP(1,Int_disc(k))-X).^2 ...
            + (IntersectP(2,Int_disc(k))-Y).^2);
        % Set infinite the distances corresponding to the same interface
        Dist_int(id<=Int_id(k)) = Inf;
        % Find the global index of the disconnected intersection point
        Ind_dist_global = [Ind_dist_global ...
            find(Int_ind(Int_ind==Int_disc(k))==Int_ind)];
        
        % Neighbour-point index
        Nindx = find(min(Dist_int)==Dist_int);
        % Neighbour point
        Neighbour = [Neighbour [X(Nindx);Y(Nindx)]];
        
        % Find where the neighbour point has to be located in the lower
        % disconnected interface, 1 stands left and 0 for right
        Neighbour_pos = [Neighbour_pos, ...
            Sgap_int(Ind_dist_global(k)-Int_id(Int_disc(k)))];
    end
    
    % Initializing parameters for the loop
    X_dint = [];
    Y_dint = [];
    intersect_id_dint = [];
    id_dint = [];
    c = 1;
    % Loop to add the neighbour points in the lower interface
    for p = 1:size(Int_disc,2)
        if Neighbour_pos(p)
            % Neighbour connected from the right
            
            % New coordinates
            X_dint = [X_dint X(c:Ind_dist_global(p)-1) Neighbour(1,p) ...
                X(Ind_dist_global(p))];
            Y_dint = [Y_dint Y(c:Ind_dist_global(p)-1) Neighbour(2,p) ...
                Y(Ind_dist_global(p))];
            
            % New intersect id. 1 for the neighbour that is the new
            % intersection point and 0 for the old intersection point
            intersect_id_dint = [intersect_id_dint ...
                intersect_id(c:Ind_dist_global(p)-1) 1 0];
            
            % New id
            id_dint = [id_dint id(c:Ind_dist_global(p)-1) ...
                Int_id(Int_disc(p)) id(Ind_dist_global(p))];
            
        else
            % Neighbour connected from the left
            
            % New coordinates
            X_dint = [X_dint X(c:Ind_dist_global(p)) Neighbour(1,p)];
            Y_dint = [Y_dint Y(c:Ind_dist_global(p)) Neighbour(2,p)];
            
            % New intersect id. 1 for the neighbour that is the new
            % intersection point and 0 for the old intersection point
            intersect_id_dint = [intersect_id_dint ...
                intersect_id(c:Ind_dist_global(p)-1) 0 1];
            
            % New id
            id_dint = [id_dint id(c:Ind_dist_global(p)) ...
                Int_id(Int_disc(p))];
            
        end
        
        % Updates counter
        c = Ind_dist_global(p)+1;
        
    end
    
    % Add the points after the last disconnected point
    X_dint = [X_dint X(c:end)];
    Y_dint = [Y_dint Y(c:end)];
    intersect_id_dint = [intersect_id_dint intersect_id(c:end)];
    id_dint = [id_dint id(c:end)];
    
    % Updates X,Y,intersect_id and id
    X = X_dint;
    Y = Y_dint;
    intersect_id = intersect_id_dint;
    id = id_dint;
    
end