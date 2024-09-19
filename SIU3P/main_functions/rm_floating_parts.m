function [X,Y,id,intersect_id] = rm_floating_parts(X,Y,id,intersect_id)
% [X,Y,ID,INTERSECT_ID] = RM_FLOATING_PARTS(X,Y,ID,INTERSECT_ID) removes
% the floating parts of the intersections indexed by ID, which
% coordinates are given by a X and a Y vectors. 
%
% Example:
% --------     _
%             / \ floating part of interface 2 that needs to be removed
%          ________
% ________/ /     \\_____________ interface 2
% _________/       \_____________ interface 1

%--------------------------------------------------------------------------
% Function written by Miguel Andres-Martinez, PhD student at Royal Holloway
% University of London, 28-02-2014. Email: mandresmartinez87@gmail.com
%--------------------------------------------------------------------------

% Initializes the boolean vector for removing the floating points of the
% interfaces
rm_parts = zeros(size(id));

% Loop through the different interfaces
for li = 1:max(id)
    % Set the evaluated interface for this iteration
    interface = [X(id==li);Y(id==li)];
    
    % Generates a parts vector which indexes the points of the evaluated
    % interface 
    parts = index_parts(intersect_id(id==li));
    
    % If there is more than one part
    if max(parts) > 1
        % Loop through the middle parts to check if they have 2
        % intersection points or less. If they have less it means they are
        % floating and should be removed
        for lp = 2:max(parts)-1
            % Calculates how many points are present in other interfaces
            % (how many intersections have this part)
            present_in_other_interface = ...
                ismember(interface(:,parts==lp),[X(id~=li);Y(id~=li)]);
            present_in_other_interface = ...
                present_in_other_interface(1,:) & ...
                present_in_other_interface(2,:);
            
            % If the part has less than 2 intersection points then is a
            % floating part and the part is stored in the rm_parts boolean
            % vector for later removing
            if sum(present_in_other_interface)<2
                rm_parts(id==li) = rm_parts(id==li) | parts==lp;
            end
        end
    end
end

% Changes rm_parts into the boolean format
rm_parts = rm_parts==1;

% Removes the floating parts from X, Y, ID and INTERSECT_ID
X(rm_parts) = [];
Y(rm_parts) = [];
id(rm_parts) = [];
intersect_id(rm_parts) = [];