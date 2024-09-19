function [int_points, indx_segment] = intersect_lines(line1,line2)
% [INT_POINTS,INDX_SEGMENT2] = INTERSECT_LINES(LINE1,LINE2) finds the
% intersection points INT_POINTS between LINE1 and LINE2, where
% INDX_SEGMENT indicates the segments of LINE1 (first row) and LINE2 
% (second row) where the intersection is. LINE1, LINE2 and INT_POINTS are 
% 2xn matrices where the row 1 contains the x-coordinates, the row 2 
% contains the y-coordinates and n accounts for the number of points 
% represented by each matrix.

%--------------------------------------------------------------------------
% Function written by Miguel Andres-Martinez, PhD student at Royal Holloway
% University of London, 21-02-2014. Email: mandresmartinez87@gmail.com
%--------------------------------------------------------------------------

int_points = [[];[]];
indx_segment = [[];[]];

% Number of decimals for comparing the points
accu = 1e10;

% Define the slope and the intersection with the coordinates of all the
% segments of the lines
[a1,b1] = line_ab(line1(1,1:end-1),line1(2,1:end-1),line1(1,2:end), ...
    line1(2,2:end));
[a2,b2] = line_ab(line2(1,1:end-1),line2(2,1:end-1),line2(1,2:end), ...
    line2(2,2:end));

% Find where the different segments intersect
for j = 1:length(a1)
    % Calculates the intersection points between each segment of line1 and 
    % all the segments in line2
    if isinf(a1(j))
        % In case the segment of line1 is vertical calculates the
        % intersection with line2 segments
        x = line1(1,j)*ones(size(a2));
        y = a2*line1(1,j)+b2;
    else
        [x,y] = intersect_ab(a1(j),b1(j),a2,b2);
    end
    
    % In case some of the segments of line2 are vertical, calculates the 
    % coordinates of the intersections with line1 segment
    x(isinf(a2)) = line2(1,isinf(a2));
    y(isinf(a2)) = a1(j)*line2(1,isinf(a2))+b1(j);
    
    % Intersection coordinates with 'accu' decimals
    x_s = floor(x*accu)/accu;
    y_s = floor(y*accu)/accu;
    
    % Line1 segment coordinates with 'accu' decimals
    l1xj = floor(line1(1,j)*accu)/accu;
    l1xj1 = floor(line1(1,j+1)*accu)/accu;
    l1yj = floor(line1(2,j)*accu)/accu;
    l1yj1 = floor(line1(2,j+1)*accu)/accu;
    
    % Finds if the intersection points are inside the segment of line1
    in_segment_x1 = x_s > l1xj & x < l1xj1 | x_s < l1xj & x_s > l1xj1;
    in_segment_y1 = y_s > l1yj & y < l1yj1 | y_s < l1yj & y_s > l1yj1;
    segment2 = find(in_segment_x1 & in_segment_y1);
    
    % Take the intersection points that are in the segment of line1 and
    % evaluate if they are also in a segment of line2, which in that case
    % would mean that the two segments are intersecting and the
    % intersection point can be stored
    for k = 1:length(segment2)
        % Intersection coordinates with 'accu' decimals
        x_s = floor(x(segment2(k))*accu)/accu;
        y_s = floor(y(segment2(k))*accu)/accu;
        % Line2 segment coordinates with 'accu' decimals
        l2xk = floor(line2(1,segment2(k))*accu)/accu;
        l2xk1 = floor(line2(1,segment2(k)+1)*accu)/accu;
        l2yk = floor(line2(2,segment2(k))*accu)/accu;
        l2yk1 = floor(line2(2,segment2(k)+1)*accu)/accu;
        % Evaluate if the intersection point is inside the segment of line2
        in_segment_x2 = x_s >= l2xk & x_s <= l2xk1 | x_s <= l2xk & ...
            x_s >= l2xk1;
        in_segment_y2 = y_s >= l2yk & y_s <= l2yk1 | y_s <= l2yk & ...
            y_s >= l2yk1;
        
        % In case the intersection point is inside a segment of line1 and
        % also inside a segment of line2 then the point is stored into
        % int_points, together with the index of the segments of line1 and
        % line2.
        if in_segment_x2 && in_segment_y2
            int_points  = [int_points [x(segment2(k)); y(segment2(k))]];
            indx_segment = [indx_segment [[j];[segment2(k)]]];
        end
    end
end

% Removes repeated int_points
rm_non_unique = [];
[unique_points,order] = unique(floor(int_points*accu)'/accu,'rows');
unique_points = unique_points';
for k = 1:size(unique_points,2)
    unique_x = unique_points(1,k) == floor(int_points(1,:)*accu)/accu;
    unique_y = unique_points(2,k) == floor(int_points(2,:)*accu)/accu;
    unique_ind = find(unique_x & unique_y);
    rm_non_unique = [rm_non_unique unique_ind(2:end)];
end

int_points = int_points(:,sort(order));
indx_segment(:,rm_non_unique) = [];