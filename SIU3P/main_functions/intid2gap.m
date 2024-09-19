function [Sgap] = intid2gap(intersect_id,Geo_id,id,include_lat)

% [SGAP] = INTID2GAP indicates which segments of segmentlist in the 
% function generate mesh correspond to gaps. INCLUDE_LAT should be 1 if
% the lateral boundaries need to be included and 0 if they are not needed.
%
% Example:
%                                    segmentlist gap in interface 1
%                                                 |
%______        ______ Interface 2       ______    |   ______ Interface 2
%      \______/                               \....../
%      /      \                               /      \
%_____/        \_____ Interface 1       _____/        \_____ Interface 1

%--------------------------------------------------------------------------
% Function written by Miguel Andres-Martinez, PhD student at Royal Holloway
% University of London, 16-07-2014. Email: mandresmartinez87@gmail.com
%--------------------------------------------------------------------------

% Sum every intersect_id with the folowing intersect_it. This is done to
% change from coordinate ids to segments (which refers to the line between
% connected coordinates). Where the sum is 2 will mean there are 2
% intersection points and, therefore, there is a gap
Segment_gap = (intersect_id(1:end-1) + intersect_id(2:end)) == 2;

n_layers = (max(Geo_id)-1)/3;

% Fix mistakes due to contiguous gaps
%
% Example:
%
%           ..... gap
%           _____ segment              this shouldn't be a gap
%                                                 |
% Segment_gap  =    0     1     0     0     1     1     1     0
% intersect_id = 0     1     1     0     1     1     1     1     0
%                o_____o.....o_____o_____o.....o_____o.....o_____o
%
for j = 1:size(Segment_gap,2)
    if Segment_gap(j)==1 && Segment_gap(j-1)==1
        Segment_gap(j) = 0;
    end
end

% Initialize the last_segment with the last segment of the bottom boundary
last_segment = find(id==1,1,'last')-1;

% Save the Segment_gap of the bottom boundary into Sgap
Sgap = Segment_gap(1:last_segment);

% Loop through the number of layers of the model
for j = 1:n_layers
    % Find the first and the last segments of the horizontal upper boundary
    % of the j layer
    first_segment = last_segment+2;
    last_segment = find(id==j+1,1,'last')-1;
    % Add to Sgap Segment_gap of the horizontal upper boundary of the j 
    % layer, and zeros for the segments of the lateral boundaries
    Sgap = [Sgap zeros(1,(sum(Geo_id==j*3-1)+1)*include_lat) ...
        Segment_gap(first_segment:last_segment) ...
        zeros(1,(sum(Geo_id==j*3+1)+1)*include_lat)];
end

% Make Sgap a boolean vector
Sgap = Sgap == 1;