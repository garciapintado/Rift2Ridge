function [Dist] = dist_points(P1,P2)
% [DIST] = DIST_POINTS(P1,P2) calculates the distance DIST between 
% the points P1 and P2 given their coordinates P1(:,1), P1(:,2) and 
% P2(:,1), P2(:,2). It also calculates the distances between coordinate
% vectors, so distances between different pairs of points can be calculated
% at the same time.

%--------------------------------------------------------------------------
% Function written by Miguel Andres-Martinez, PhD student at Royal Holloway
% University of London, 22-06-2013. Email: mandresmartinez87@gmail.com
%--------------------------------------------------------------------------

% Difference in coordinates
X = P2(:,1)-P1(:,1);
Y = P2(:,2)-P1(:,2);

% Calculate the distance
Dist = sqrt(X.^2 + Y.^2);

% nearest neighbour from a vector to a point:
% [nnd,nni] = (dist_points(pointlist(2:3,:)',crashp))