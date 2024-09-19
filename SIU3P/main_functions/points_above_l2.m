function [gap] = points_above_l2(line1,line2)
% [GAP] = POINTS_ABOVE_L2(LINE1,LINE2) finds the points of LINE1 which are
% not contained inside of the area between LINE2 and the bottom of the
% domain, and store them in GAP. Both LINE1 and LINE2 have to be composed
% by consecutive points ordered in the same direction.

%--------------------------------------------------------------------------
% Function written by Miguel Andres-Martinez, PhD student at Royal Holloway
% University of London, 21-02-2014. Email: mandresmartinez87@gmail.com
%--------------------------------------------------------------------------

% Creates a sign s1 defined by the way the LINE vectors are ordered
s1 = sign(line1(1,1)-line1(1,end));
gap = [];

% Loop through the LINE1 to find which points are above LINE2
for j = 1:size(line1,2)
    % Finds which segments of LINE2 contains LINE1(J), in the x-coordinates
    in_segment2 = find( ...
        line1(1,j) <= line2(1,1:end-1) & line1(1,j) >= line2(1,2:end) | ...
        line1(1,j) >= line2(1,1:end-1) & line1(1,j) <= line2(1,2:end));
    
    % Calculates the slopes and the y-coordinate intersection of the
    % segments in_segment2
    [a2,b2] = line_ab(line2(1,in_segment2),line2(2,in_segment2), ...
        line2(1,in_segment2+1),line2(2,in_segment2+1));
    
    % Interpolates the y-coordinates at the segments of LINE2 for the LINE1
    % x-coordinate
    y2 = a2*line1(1,j)+b2;
    
    % Finds the closest y2 for the y-coordintate of LINE1(j)
    closest = find(min(abs(y2-line1(2,j)))==abs(y2-line1(2,j)));
    
    % Calculates the sign for the comparison of the y-coordinates of the 
    % LINE1 and the interpolated LINE2, so that a LINE1 point that is above
    % a LINE2 segment but inside the area between LINE2 and the bottom of
    % the domain, is not included in GAP
    %              ______
    %             /      \
    %            /        |
    %            \.LINE1   \
    %             \ point   \
    %  LINE2 _____/          \___________
    s = s1*sign(line2(1,in_segment2(closest)) ... 
        - line2(1,in_segment2(closest)+1));
    
    if ~isempty(y2)
        % Checks if the LINE1 point is outside of the domain defined by
        % LINE, and in case it is stores its index in GAP
        if (s(1)*line1(2,j))>=(s(1)*y2(closest(1)))
            gap = [gap j];
        end
    end
end