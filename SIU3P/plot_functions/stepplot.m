function [] = stepplot(X,Y,color)
% STEPPLOT(X,Y) plots a line point by point. Hit ENTER to continue for the
% next point

%--------------------------------------------------------------------------
% Function written by Miguel Andres-Martinez, PhD student at Royal Holloway
% University of London, 05-01-2015. Email: mandresmartinez87@gmail.com
%--------------------------------------------------------------------------

for j = 1:size(X,2)
    hold on
    plot(X(j),Y(j),color)
    pause
end