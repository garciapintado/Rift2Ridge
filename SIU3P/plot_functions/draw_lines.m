function [X,Y,id] = draw_lines(nlines,ax)
% [X,Y,ID] = draw_lines(NLINES,AX) opens a GUI interface to draw a number
% of lines, defined by NLINES, in a rectangular space defined by AX, which
% is a 4-term vector which AX(1) and AX(2) are the x-coordinate minimum and
% maximum respectively, and AX(3) and AX(4) are the y-coordinate minimum
% and maximum respectively. X and Y return the coordinates of the drawn
% points and ID assigns an id to relate each point with its line. In order
% to finish a line and start a new one, you will need to mark a point
% outside of the space.

%--------------------------------------------------------------------------
% Function written by Miguel Andres-Martinez, PhD student at Royal Holloway
% University of London, 21-02-2014. Email: mandresmartinez87@gmail.com
%--------------------------------------------------------------------------

X =[];
Y = [];
id = [];
color_l = rand(nlines,3);
count=1;
for l = 1:nlines
    keep_drawing = 1;
    while keep_drawing
        plot(X(id==l),Y(id==l),'-','Color',color_l(l,:))
        plot(X(id==l),Y(id==l),'.','Color',color_l(l,:))
        axis (ax)
        hold on
        [x,y] = ginput(1);
        keep_drawing = x>=ax(1) & x<=ax(2) & y>=ax(3) & y<= ax(4);
        if keep_drawing
            X(count) = x;
            Y(count) = y;
            id(count) = l;
            count = count+1;
        end
    end
end

hold off