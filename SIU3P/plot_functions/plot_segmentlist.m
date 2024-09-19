% Plots the segments sequentially. Remember to set up the axis before using
% it

for j_plot = 1:size(segmentlist,2)
    plot(pointlist(1,segmentlist(1:2,j_plot))/1000,pointlist(2,segmentlist(1:2,j_plot))/1000,'-k')
    hold on
    plot(pointlist(1,segmentlist(1:2,j_plot))/1000,pointlist(2,segmentlist(1:2,j_plot))/1000,'.k')
    pause(0.1)
end

hold off