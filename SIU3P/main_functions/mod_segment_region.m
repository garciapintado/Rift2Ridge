function [segmentlist,regionlist] = mod_segment_region( ...
    segmentlist,regionlist,Mod_segment_ind,GEOMETRY,Elsizes)

% THIS FUNCTION REQUIERES COMMENTS!
%
% Function developed by Miguel Andres-Martinez 10/02/2014

% Modify regionlist 1
%--------------------
Mod_segment_ind_region = Mod_segment_ind(:,1:2:end);
% Gets indexes and ids of the new segments for the first intersection
% in each interface gap to calculate the new regionlist points
Mod_segment_ind_region(:,Mod_segment_ind_region(3,:)==0) = [];
% Removes the columns of Mod_segment_ind_region for wich the
% correspondent segment is above an overlying interface (id = 0)

% Calculates the X and Y coordinates of the center of a triangle which
% vertex correspond to the indexes Mod_segment_ind_region(1,:),
% Mod_segment_ind_region(2,:) and Mod_segment_ind_region(2,:)-1:
%
%            M2-1.__________ Layer 2
%               / **********
% Layer 2_____./x Center of the triangle
% """"""""""""^\**** Phase n+1
% """"""""""""| \.__________ Layer 1
% """"""""""""|""^""""""""""
% ""Phase n" M2" M1"" Phase n
% """"""""""""""""""""""""""

% In case there is a common intersection point
Common_in = zeros(1,size(Mod_segment_ind_region,2));
Common_c = 1;
for check_i = 1:size(Mod_segment_ind_region,2)
    if Common_in(check_i)==0
        Common_in(Mod_segment_ind_region(2,check_i) == ...
            Mod_segment_ind_region(2,:)) = Common_c;
        Common_c = Common_c+1;
    end
end

Upper_vertex = Mod_segment_ind_region(2,:)-1;

if sum(Common_in)~=0
    for check_j = 1:max(Common_in)
        check_bool = Common_in == check_j;
        ind_check = Mod_segment_ind_region(3,:);
        ind_check(~check_bool) = 0;
        for check_k = 1:sum(check_bool);
            if check_k == 1
                fmax = max(ind_check)==ind_check;
                ind_check(fmax) = 0;
            else
                Upper_vertex(max(ind_check)==ind_check) = ...
                    Mod_segment_ind_region(1,fmax);
                fmax = max(ind_check)==ind_check;
                ind_check(fmax) = 0;
            end
        end
    end
end

X_region_new = sum(reshape(GEOMETRY(1,[Mod_segment_ind_region(1,:) ...
    Mod_segment_ind_region(2,:) Upper_vertex]), ...
    size(Mod_segment_ind_region,2),3)')/3;
Y_region_new = sum(reshape(GEOMETRY(2,[Mod_segment_ind_region(1,:) ...
    Mod_segment_ind_region(2,:) Upper_vertex]), ...
    size(Mod_segment_ind_region,2),3)')/3;

regionlist_new = [X_region_new; Y_region_new; ...
    Mod_segment_ind_region(3,:)/3+1; ...
    Elsizes(Mod_segment_ind_region(3,:)/3+1)];
% Generates the new regionlist matrix. Mod_segment_ind_region(3,:)/3+1
% calculates which Phases should it be

% Modify segmentlist
%-------------------
segmentlist(:,sum(ismember(segmentlist(1:2,:), ...
    reshape(Mod_segment_ind(1,:),2,size(Mod_segment_ind,2)/2)))==2) = ...
    [];
% Remove the segments generated for interfaces that should be
% discontinuous for this segments
Mod_segment_ind(:,Mod_segment_ind(3,:)==0) = [];
% Removes the columns of Mod_segment_ind for wich the correspondent
% segment is above an overlying interface (id = 0)
segmentlist = [segmentlist Mod_segment_ind];
% Updates the segmentlist with Mod_segment_ind

% Modify regionlist 2
%--------------------
regionlist = [regionlist regionlist_new];
% Updates the regionlist with regionlist_new

