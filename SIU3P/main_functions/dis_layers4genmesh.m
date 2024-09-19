function [regionlist,segmentlist] = dis_layers4genmesh ... 
    (pointlist,regionlist,segmentlist,Sgap,Geo_id,Elsizes)
% [REGIONLIST,SEGMENTLIST] = DIS_LAYERS4GENMESH(POINTLIST,REGIONLIST,
% SEGMENTLIST,SGAP,GEO_ID,ELSIZES) adds to the REGIONLIST variable of the
% function generate_mesh, the regionlist points corresponding to the pieces
% of the discontinuous layers of the model, using the geometry parameters
% (POINTLIST, REGIONLIST, SEGMENTLIST, GEO_ID, and ELSIZES) and the indexes
% for the gap segments SGAPS, given by the function dis_layers. It also 
% removes the gap segments from SEGMENTLIST

%--------------------------------------------------------------------------
% Function written by Miguel Andres-Martinez, PhD student at Royal Holloway
% University of London, 02-07-2014. Email: mandresmartinez87@gmail.com
%--------------------------------------------------------------------------

% 09/09/2014 MA
    % Fixed problems related with unsorted P3 that was generating a problem
    % at the calculation of the regionlist points when more than one gap
    % was happening on the model

%==========================================================================
% CREATE REGIONLIST POINTS FOR THE DISCONTINUOUS LAYERS TO BE FILLED WITH
% THE CORRESPONDING PHASE OF TRIANGLES
%==========================================================================

% Calculates the x and y coordinates of a point P inside of a triangle
% which vertex are the intersection point P1 and the adjacent 2 points of
% the interfaces intersecting P2 and P3. P is located at a line that
% connects P1 with the middle point of P2-P3, and at a distance d (1e2)
% from P1. This ensures that the region point P is between the two
% interfaces and close enough to the intersection, so that very small
% pieces of layers will not be a problem for the mesh generator
%
%              P3.__________ Interface 2
%               /
% Layer 2_____./x P
%             ^\      Phase n+1
%             | \.__________ Interface 1
%             |  ^
%  Phase n   P1  P2   Phase n
%

% First point in the gaps (right point of the gap segment)
first_p_gap = segmentlist(1,Sgap);

% Store coordinates of P1 and P2
P1 = pointlist(:,first_p_gap);
P2 = pointlist(:,first_p_gap-1);

% Indexes of P1 in pointlist. LOC will be a zero vectors where would be n
% where the n P1 is found (n = 1,2,...,number of P1s). Note that P1
% (intersection point) is repeated at interfaces 1 and 2
[~,LOC] = ismember(pointlist',P1','rows');

% Find first repetead points which are the ones on the interface evaluated
[~,Ind_rep] = unique(LOC,'first');
Ind_rep(1) = [];
% Find the common intersection points, both in the evaluated interface and
% the ones at interfaces above
P3_ind = find(LOC~=0);
% Remove repeated intersection points
P3_ind(ismember(P3_ind,Ind_rep)) = [];
% Order for P3_ind
[~,Order] = ismember(P1',pointlist(:,P3_ind)','rows');
% Reorder P3_ind
P3_ind = P3_ind(Order);
% Previous point
P3_ind = P3_ind-1;
% Store the coordinates of P3
P3 = pointlist(:,P3_ind);

% Calculate middle point P2-P3 coordinates (Mp)
Mp = (P2+P3)/2;
% Calculate the distance between P1 and Mp
D = sqrt(sum((Mp-P1).^2,1));
% Minimum size of the discontinuous layer
d = 1e-2;
% Calculate the region point P
P = repmat((d./D),2,1).*(Mp-P1) + P1;

% Find the phase of the region point P
regionlist_dis = (Geo_id(first_p_gap))/3+1;
regionlist_dis(regionlist_dis<2) = 1;

% Store the regionlist points of the layer pieces
regionlist = [regionlist [P; regionlist_dis; Elsizes(regionlist_dis)]];

%==========================================================================
% REMOVE GAPS FROM SEGMENTS
%==========================================================================
segmentlist(:,Sgap) = [];