function esize(GCOORD,ELEM2NODE,Point_id,Corner_id,Cornin_id)
% ESIZE(GCOORD,ELEM2NODE) returns the average, the
% maximum and the minimum of the edge longitud of the elements around an
% interface. 
%
% Example:
%     esize(GCOORD,ELEM2NODE)

km = 1000;
el_index = ELEM2NODE(1:3,:);
el_indp1 = [el_index(2:3,:); el_index(1,:)];
int_dist = ((GCOORD(1,el_index)-GCOORD(1,el_indp1)).^2 + ...
    (GCOORD(2,el_index)-GCOORD(2,el_indp1)).^2).^(1/2);
int_dist = reshape(int_dist,3,size(el_index,2));
av_dist = sum(int_dist)./length(int_dist)/km;
[max_dist,ii] = max(int_dist./km,[],2);
[min_dist,ij] = min(int_dist./km,[],2);
disp(['The average distance is ', num2str(min(av_dist)),' km'])
disp(['The maximum distance is ', num2str(min(max_dist)),' km'])
disp(['The minimum distance is ', num2str(min(min_dist)),' km'])

plot_mesh
hold on
for j =1:3
plot(GCOORD(1,ELEM2NODE(1:3,ij(j)))/1000,GCOORD(2,ELEM2NODE(1:3,ij(j)))/1000,'^b');
plot(GCOORD(1,ELEM2NODE(1:3,ii(j)))/1000,GCOORD(2,ELEM2NODE(1:3,ii(j)))/1000,'^r');
end

