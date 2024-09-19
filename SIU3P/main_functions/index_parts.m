function [parts] = index_parts(intersect_id)
% [PARTS] = INDEX_PARTS(INTERSECT_ID) generates a 0s PARTS vector where the
% points of the different parts are tagged with a different id for each
% part, using the INTERSECT_ID boolean vector that defines the intersection
% points:
%
% Example:
% --------
% INTERSECT_ID   1   0    0  1   x    1 0    1     1     0 0     0     1
% PARTS          1   1    1  1   x    2 2    2     3     3 3     3     3

% c   1          2    
% u   0           
%                            .   .    .
%                           /          \.____.     ._____.
%                .___.____./                              \._____._____.
%                                     |      |     |
%                |________|           |______|     |___________________|
%                    1                   2                   3



%--------------------------------------------------------------------------
% Function written by Miguel Andres-Martinez, PhD student at Royal Holloway
% University of London, 21-02-2014. Email: mandresmartinez87@gmail.com
%--------------------------------------------------------------------------

parts = zeros(size(intersect_id));
c = 1;
u = 0;
for j = 1:length(intersect_id)
    if intersect_id(j) == 0 
        parts(j) = c;
    elseif u == 0                                                          % include point in part and close part      
        parts(j) = c;
        c = c+1;                                                           % new part
        u = 1;
    else
        parts(j) = c;
        u = 0;
    end
end

% parts = zeros(size(intersect_id));
% % Points that are not intersections
% parts_no_edges = ~intersect_id;
% c = 1;
% 
% % Find the correspondent ids for the different parts
% for j = 1:length(parts_no_edges)
%     % If is not an intersection assigns PARTS(j) a value of c
%     if parts_no_edges(j) == 1
%         parts(j) = c;
%         % If j does not correspond to the las point of the line and the
%         % following point is an intersection, it means that a part is
%         % is finishing and adds 1 to c, so a new id is generted for the
%         % following part
%         if j<length(parts_no_edges)
%             if parts_no_edges(j+1) == 0
%                 c = c+1;
%             end
%         end
%     end
% end
% 
% % Stores the ids of the intersection points corresponding to the part they
% % belong to
% for k = 1:max(parts)
%     first_point = find(k==parts,1,'first');
%     last_point = find(k==parts,1,'last');
%     % The ifs are to avoid calling points out of the edges of the part
%     if first_point ~= 1
%         parts(first_point-1) = k;
%     end
%     if last_point ~= size(parts,2)
%         parts(last_point+1) = k;
%     end
% end