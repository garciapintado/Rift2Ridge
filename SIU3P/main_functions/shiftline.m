function [SHIFT] = shiftline(X,Y,m)
% [SHIFT] = SHIFTLINE(X,Y,M) calculates the upward shift of the profile 
% given by X and Y, where the magnitude of the shift is m.
%
% Note: vertical segments might have problems and need further debugging.

%--------------------------------------------------------------------------
% Function written by Miguel Andres-Martinez, PhD student at Royal Holloway
% University of London, 14-07-2014. Email: mandresmartinez87@gmail.com
%--------------------------------------------------------------------------

%==========================================================================
% SLOPES AND DISTANCES OF THE SEGMENTS
%==========================================================================

% Calculate slopes of the segments and constant
[a,b] = line_ab(X(1:end-1),Y(1:end-1),X(2:end),Y(2:end));
% Calculate the distance between the points of the segments
modul = sqrt(diff(X).^2+diff(Y).^2);
% Store the distances in pairs of successive segments
Modul = [modul(1:end-1);modul(2:end)];

%==========================================================================
% POINTS BEFORE AND AFTER
%==========================================================================
% Note: points before and after are at the same distance from the reference
% point

% Points before
% -------------

% Increment in x due to a slope a and the minimum distance of Modul
SlopeP1 = min(Modul)./sqrt(a(1:end-1).^2+1);
% Apply a positive increment of x
SlopeP1pos = SlopeP1+X(2:end-1);
% Apply a negative increment of x
SlopeP1neg = -SlopeP1+X(2:end-1);
% Find if SlopeP1pos is in the segment
indpos = floor(SlopeP1pos*100) >= floor(X(1:end-2)*100) & ...
    floor(SlopeP1pos*100) <= floor(X(2:end-1)*100) | ...
    floor(SlopeP1pos*100) <= floor(X(1:end-2)*100) & ...
    floor(SlopeP1pos*100) >= floor(X(2:end-1)*100);
% Find if SlopeP1neg is in the segment
indneg = floor(SlopeP1neg*100) >= floor(X(1:end-2)*100) & ...
    floor(SlopeP1neg*100) <= floor(X(2:end-1)*100) | ...
    floor(SlopeP1neg*100) <= floor(X(1:end-2)*100) & ...
    floor(SlopeP1neg*100) >= floor(X(2:end-1)*100);
% Values of x-coordinates inside the segment
SlopeP1(indpos) = SlopeP1pos(indpos);
SlopeP1(indneg) = SlopeP1neg(indneg);
% Calculate the y-coordinates of the points Slope1
SlopeP1 = [SlopeP1; a(1:end-1).*SlopeP1+b(1:end-1)];

% Points after
% -------------

% Increment in x due to a slope a and the minimum distance of Modul
SlopeP2 = min(Modul)./sqrt(a(2:end).^2+1);
% Apply a positive increment of x
SlopeP2pos = SlopeP2+X(2:end-1);
% Apply a negative increment of x
SlopeP2neg = -SlopeP2+X(2:end-1);
% Find if SlopeP2pos is in the segment
indpos = floor(SlopeP2pos*100) >= floor(X(2:end-1)*100) & ...
    floor(SlopeP2pos*100) <= floor(X(3:end)*100) | ...
    floor(SlopeP2pos*100) <= floor(X(2:end-1)*100) & ...
    floor(SlopeP2pos*100) >= floor(X(3:end)*100);
% Find if SlopeP2neg is in the segment
indneg = floor(SlopeP2neg*100) >= floor(X(2:end-1)*100) & ...
    floor(SlopeP2neg*100) <= floor(X(3:end)*100) | ...
    floor(SlopeP2neg*100) <= floor(X(2:end-1)*100) & ...
    floor(SlopeP2neg*100) >= floor(X(3:end)*100);
% Values of x-coordinates inside the segment
SlopeP2(indpos) = SlopeP2pos(indpos);
SlopeP2(indneg) = SlopeP2neg(indneg);
% Calculate the y-coordinates of the points Slope2
SlopeP2 = [SlopeP2; a(2:end).*SlopeP2+b(2:end)];

%==========================================================================
% SLOPES PERPENDICULAR TO THE PROFILE
%==========================================================================
% Slopes at the points
[Slope,~] = line_ab(SlopeP1(1,:),SlopeP1(2,:),SlopeP2(1,:),SlopeP2(2,:));

% Slopes perpendicular to the profile
Slope = -1./Slope;

% % Uncomment to plot slopes
% bp = Y(2:end-1)-X(2:end-1).*Slope;
% for k = 1:size(Slope,2)
%     hold on
%     plot([-200 200],bp(k)+Slope(k)*[-200 200],'--r')
% end

%==========================================================================
% SHIFT
%==========================================================================
% X-shift
x = m./sqrt((1+Slope.^2));
% Y-shift
y = Slope.*x;
% Add first and last point to shift
x = [0 x 0];
y = [0 y 0];
% In case the profile is horizontal, gives y_shift values of m
y(isnan(y)) = m;
% In case the profile is vertical, gives x_sift values of m
vertical = (isnan(x));
x(vertical) = m;
% Find the points below the profile
[below] = points_above_l2([X;Y],[X+[0 x(2:end-1) 0];Y+[m y(2:end-1) m]]);
% Remove the first and the last point in case is contained in the below
% vector
below(below==1|below==length(x)) = [];
% In case the point is below the profile x sholud be -x
x(below) = -x(below);
% Recalculate the y for the points below the profile
y(below) = Slope(below-1).*x(below);
% In case the profile is horizontal, gives y_shift values of m
y(isnan(y)) = -m;
y(vertical) = 0;

% Write output
SHIFT = [x;y];