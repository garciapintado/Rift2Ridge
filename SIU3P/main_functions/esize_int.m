function esize_int(GCOORD,ELEM2NODE,Point_id,interface)
% ESIZE_INT(GCOORD,ELEM2NODE,POINT_ID,INTERFACE) returns the average, the
% maximum and the minimum of the edge longitud of the elements around an
% interface. INTERFACE is the id of the interface you want to look at,
% referenced in the same way as Point_id.
%
% Example:
%     esize_int(GCOORD,ELEM2NODE,Point_id,12)

km = 1000;
int_nodes = find(Point_id==interface);
int_el = sum(ismember(ELEM2NODE,int_nodes),1)>=1;
el_index = ELEM2NODE(1:3,int_el);
el_indp1 = [el_index(2:3,:); el_index(1,:)];
int_dist = ((GCOORD(1,el_index)-GCOORD(1,el_indp1)).^2 + ...
    (GCOORD(2,el_index)-GCOORD(2,el_indp1)).^2).^(1/2);
av_dist = sum(int_dist)/length(int_dist)/km;
max_dist = max(int_dist/km);
min_dist = min(int_dist/km);
disp(['The average distance is ', num2str(av_dist),' km'])
disp(['The maximum distance is ', num2str(max_dist),' km'])
disp(['The minimum distance is ', num2str(min_dist),' km'])

