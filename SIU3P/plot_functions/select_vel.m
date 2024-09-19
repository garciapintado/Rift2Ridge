% SELECT POINTS TO PLOT VELOCITIES 

%--------------------------------------------------------------------------
% Function written by Miguel Andres-Martinez, PhD student at Royal Holloway
% University of London, 27-09-2015. Email: mandresmartinez87@gmail.com
%--------------------------------------------------------------------------

%==========================================================================
% SELECT VELOCITIES
%==========================================================================
plot_eri
hold on
axis_p = axis;

keep_drawing = 1;
counter = 1;
X = [];
while keep_drawing
    [x,y] = ginput(1);
    keep_drawing = x>=axis_p(1) & x<=axis_p(2) & ...
        y>=axis_p(3) & y<= axis_p(4);
    if keep_drawing
        plot(x,y,'ok')
        X(1,counter) = x;
        X(2,counter) = y;
        counter = counter+1;
    end
end

% Find neighbours indexes
vel_indx = knnsearch(GCOORD'/1000,X');