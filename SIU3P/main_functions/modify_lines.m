function [new_line1,new_line2,intersect_id,int_indx2] = ...
    modify_lines(line1,line2,int_points,indx_segment,gap)
% [NEW_LINE1,NEW_LINE2,INTERSECT_ID] = MODIFY_LINES(LINE1,LINE2,INT_POINTS,
% INDX_SEGMENT,GAP) modifies LINE1 and LINE2 to add the intersection points 
% INT_POINTS in the segments indicated by INDX_SEGMENT and to remove the 
% points indicated by GAP. LINE1 and LINE2 are matrices with 2 rows in
% which row 1 corresponds to the x-coordinates and row 2 to the 
% y-coordinate. NEW_LINE1 and NEW_LINE2 are the modified LINE1 and LINE2, 
% and INTERSECT_ID is a boolean vector with a length of the LINE1, which 
% 1-value terms correspond to the points of LINE1 correspond to an 
% intersection with LINE2.

%--------------------------------------------------------------------------
% Function written by Miguel Andres-Martinez, PhD student at Royal Holloway
% University of London, 21-02-2014. Email: mandresmartinez87@gmail.com
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
%                                LINE1
%--------------------------------------------------------------------------
% Creates NaN vector with the size of LINE1 plus the size of the
% intersection points (INT_POINTS)
new_line1 = nan*(ones(2,size(line1,2)+size(int_points,2)));

% Creates a vector of 0s with the size of NEW_LINE1
intersect_id = zeros(1,size(new_line1,2));

% Takes the indexes of the segments where there are intersections
% (INDX_SEGMENT) and adds to it a vector [1 2 3 ... n] where n is the
% length of the segment index vector, so the index of the segments where
% there are intersections is obtained for the NEW_LINE1
int_indx1 = indx_segment(1,:) + (1:size(indx_segment,2));

% Adds the coordinates of the int_points in the NEW_LINE1 matrix
new_line1(:,int_indx1) = int_points;
% Adds the coordinates of the old LINE1 to the NEW_LINE1 matrix
new_line1(:,isnan(new_line1(1,:))) = line1;
% Assigns 1 to the INTERSECT_ID for the intersection points
intersect_id(int_indx1) = 1;

% Creates a vector of NaN with the size of LINE1 plus the size of the
% intersection points (INT_POINTS)
move_gap = nan*(ones(1,size(line1,2)+size(int_points,2)));
% Assigns 0 to the INTERSECT_ID for the intersection points
move_gap(int_indx1) = 0;
% Adds the old gap index vector to move_gap
move_gap(isnan(move_gap)) = 1:size(line1,2);
% Finds the new indexes of the old gap points in the move_gap vector
new_gap = ismember(move_gap,gap);

% Removes the gap points in the NEW_LINE1
new_line1(:,new_gap) = [];
intersect_id(new_gap) = [];

%--------------------------------------------------------------------------
%                                LINE2
%--------------------------------------------------------------------------

% Creates NaN vector with the size of LINE2 plus the size of the
% intersection points (INT_POINTS)
new_line2 = nan*(ones(2,size(line2,2)+size(int_points,2)));

% Takes the indexes of the segments where there are intersections
% (INDX_SEGMENT) and adds to it a vector [1 2 3 ... n] where n is the
% length of the segment index vector, so the index of the segments where
% there are intersections is obtained for the NEW_LINE2
int_indx2 = indx_segment(2,:) + (1:size(indx_segment,2));

% Adds the coordinates of the int_points in the NEW_LINE2 matrix
new_line2(:,int_indx2) = int_points;
% Adds the coordinates of the old LINE2 to the NEW_LINE2 matrix
new_line2(:,isnan(new_line2(1,:))) = line2;