function [isocx,isocy] = ero_ISOC_perp(isocx,isocy,ISOC,isoc_above)
% [ISOCX,ISOCY] = ERO_ISOC_PERP(ISOCX,ISOCY,ISOC,ISOC_ABOVE) erodes the
% sediment time-line defined by ISOCX, ISOCY coordinates with isochrons
% above defined by coordinates ISOC and indexes ISOC_ABOVE. Here, erosion
% is calculated by calculating downward perpendiculars of the time-line and
% checking for intersection between this perpendiculars and time-lines
% above. Where intersections occur erosion has taken place and therefore
% the time-line is moved to the intersection at that point. This 
% function is to solve problems when resampling of the ISOCHRONS is not
% possible due to folds with overturned limbs in the sediments, that make
% the time-line not function of x but a parametric function.

%--------------------------------------------------------------------------
% Function written by Miguel Andres-Martinez, Postdoc at University of
% Bremen, 06-06-2018. Email: andresma@uni-bremen.de
%--------------------------------------------------------------------------

% Initialisation
% --------------
% Range for intersections (distance between the most separated points)
lline = sqrt(diff(isocx([1 end])).^2 + diff(isocy([1 end])).^2);
% Resolution for detecting perpendiculars [1 m]
res = 1;
% Calculating slopes of the time-line
[a,~] = line_ab(isocx(1:end-1),isocy(1:end-1),isocx(2:end), ...
    isocy(2:end));
% Coordinate vectors of the time-line excluding the edges
x1 = isocx(2:end-1);
y1 = isocy(2:end-1);
% Calculating perpendiculars
% --------------------------
% What is called here perpendiculars to the time-line are found in the 
% time-line points. This means that the "perpendiculars" are not exactly
% perpendiculars, but lines that divide the angles in between segments of
% the time line in two equal parts.
%--------------perp2
%     //|\\
%    // | \\
%   //  |  \\
%  //   |   \\
% //    |    \\==== time-line
%     perp1
% 
% There are 2 lines that separate two segments in equal angles. In order to
% choose one of them we need to calculate a SIGN
SIGN = -sign(isocx(2:end-1)-isocx(1:end-2)) ...
    .*sign(isocx(3:end)-isocx(2:end-1));
% Takes first and second slopes of every time-line segment
a1 = a(1:end-1);
a2 = a(2:end);
% Equation to calculate the slope of a line that separates the angle
% between two lines in two equal angles. Note that SIGN is used here to get
% the right line (~perpendicular to the time-line).
A = (SIGN.*sqrt((1+a1.^2).*(1+a2.^2))+a1.*a2-1)./(a1+a2);
% Calculate the line constant
B = y1-A.*x1;

% These lines have a length of lline from the center to the two edges. In
% order to solve coordinates at the two edges we solve the following 
% cuadratic equation
c1 = 1+A.^2;
c2 = -2.*x1+2.*A.*B-2.*A.*y1;
c3 = x1.^2+B.^2-2.*B.*y1+y1.^2-(lline*1e3).^2;
% Calcualte signd depending on the slope
SIGN = sign(A);
% Add to the previous calculation a change of sign if there is an
% overturned limb in a fold
SIGN = SIGN.*sign(isocx(3:end)-isocx(1:end-2));
% Solve cuadratic equation to get coordinates of the edges of the
% perpendicular lines cutting the time-line
X = [(-c2 - SIGN.*sqrt(c2.^2 - 4.*c1.*c3))./(2.*c1); ...
    (-c2 + SIGN.*sqrt(c2.^2 - 4.*c1.*c3))./(2.*c1)];
Y = A.*X+B;
% Shorten this lines to a total length of 2*res, instead of 2*lline
incX = res*(X-[x1; x1])/lline/1e3;
incY = res*(Y-[y1; y1])/lline/1e3;
X = [x1; x1] + incX;
Y = [y1; y1] + incY;

% Consider completely horizontal lines
horz = find(isinf(A));
X(1,horz) = X(2,horz);
Y(:,horz) = [y1(horz)+sign(A(horz)).*res; y1(horz)];

% Consider completely vertical lines
vert = find(isnan(X(1,:)));
X(1,vert) = X(2,vert)+sign(y1(vert+1)-y1(vert-1)).*res;
Y(:,vert) = [y1(vert); y1(vert)];

% Redefine perpendiculars exist only relatively below the time-line
% -----------------------------------------------------------------
% Create a box from the time-line down to a thickness of lline and extra
% width of 2*lline
BOX = [isocx(1)-lline isocx isocx(end)+lline isocx(end)+lline isocx(1)-lline; ...
    isocy(1) isocy isocy(end) min(isocy)-lline-1 min(isocy)-lline-1];
% Find edges of the perpendiculars inside the box, to calculate which
% perpendiculars are below the time-line even if they are located in an
% overturned fold limb
indx = inpolygon(X,Y,BOX(1,:),BOX(2,:));

% Solve again the cuadratic equation to get coordinates of the edges of the
% perpendicular lines cutting the time-line
X = [(-c2 - SIGN.*sqrt(c2.^2 - 4.*c1.*c3))./(2.*c1); ...
    (-c2 + SIGN.*sqrt(c2.^2 - 4.*c1.*c3))./(2.*c1)];
% Limit the perpendicular line to be below exclusively relatively below the
% time-line
X = [X(indx)'; x1];
Y = A.*X+B;

% Consider completely horizontal lines
horz = find(isinf(A));
X(1,horz) = X(2,horz);
Y(:,horz) = [y1(horz)+sign(A(horz)).*lline; y1(horz)];

% Consider completely vertical lines
vert = find(isnan(X(1,:)));
X(1,vert) = X(2,vert)+sign(y1(vert+1)-y1(vert-1)).*lline;
Y(:,vert) = [y1(vert); y1(vert)];

% Control perpendiculars intersecting the time-line in other points
% -----------------------------------------------------------------
% Reoder in a vector X and Y coordinates of the perpendicular lines
x = X(:)';
y = Y(:)';

% Find intersections between perpendiculars and their time-line
[per_int_isoc,ispii] = intersect_lines(...
    [isocx(2:end-1); isocy(2:end-1)],[x; y]);
% Check if intersection occurs in any point
if ~isempty(per_int_isoc)
    % Remove intersections of segments connecting different perpendiculars 
    % (where ispii is an even number
    per_int_isoc(:,ispii(2,:)/2==floor(ispii(2,:)/2)) = [];
    ispii(:,ispii(2,:)/2==floor(ispii(2,:)/2)) = [];
    
    % Unique of intersection indexes
    ii = unique(ispii(2,:));
    % Loop through the different intersecting perpendiculars
    for q = 1:length(ii)
        % Find coordinates of the intersections points within the same
        % perpendicular
        xx_yy = per_int_isoc(:,ispii(2,:)==ii(q));
        % Get the coordinates of time-line point corresponding to the
        % perpendicular
        isocxy = [x(ispii(2,ispii(2,:)==ii(q))+1); ...
            y(ispii(2,ispii(2,:)==ii(q))+1)];
        % Calculate the distance between the time-line point and the
        % intersections
        dist_si = sqrt(sum((xx_yy-isocxy).^2));
        
%         % Plot (uncomment)
%         plot(xx_yy(1,:)/1000,xx_yy(2,:)/1000,'+b','Linewidth',3)
%         plot(isocxy(1,:)/1000,isocxy(2,:)/1000,'+k','Linewidth',3)
        
        % If the distance is larger than a tolerance
        if dist_si>1e-3
            % Create an index for the smaller distance (intersection point
            % closest to the point in the time-line
            [~,indx] = min(dist_si);
            % Half of the distance separating the time-line point from the
            % closest intersection
            x_y = (xx_yy(:,indx)+isocxy(:,indx))/2;
            % Redefine perpendicular's edge to be half way between it's
            % origin and the first intersection with the time-line
            X(ispii(2,ispii(2,:)==ii(q))) = x_y(1);
            Y(ispii(2,ispii(2,:)==ii(q))) = x_y(2);
        end
    end
end

% % Plot perpendiculars (uncomment)
% plot(X(:)'/1000,Y(:)'/1000,'--r')
% hold on
% plot(X/1000,Y/1000,'-r','linewidth',2)
% plot(isocx/1000,isocy/1000,'k','linewidth',2)

% Erosion
% -------
% Loop through interfaces for calculating erosion
for o = isoc_above
    % Find intersections between perpendiculars and the time-lines above
    [int_points,indx_segment] = ...
        intersect_lines(ISOC{o}(1:2,:),[X(:)'; Y(:)']);
    % Remove intersections of segments connecting different perpendiculars 
    % (where ispii is an even number
    int_points(:,indx_segment(2,:)/2==floor(indx_segment(2,:)/2)) = [];
    indx_segment(:,indx_segment(2,:)/2==floor(indx_segment(2,:)/2)) = [];
    
%     % Plot (uncomment)
%     plot(int_points(1,:)/1000,int_points(2,:)/1000,'or','Markerfacecolor','r')
%     plot(ISOC{o}(1,:)/1000,ISOC{o}(2,:)/1000,'r')
%     plot(x(indx_segment(2,:)+1)/1000,y(indx_segment(2,:)+1)/1000,'ob','Linewidth',3)
%     plot(isocx((indx_segment(2,:)+1)/2+1)/1000,isocy((indx_segment(2,:)+1)/2+1)/1000,'ob','Linewidth',3)
    
    % Redefine coordinates of the evaluated time-line where erosion takes
    % place (intersection with the perpendiculars)
    isocx((indx_segment(2,:)+1)/2+1) = int_points(1,:);
    isocy((indx_segment(2,:)+1)/2+1) = int_points(2,:);
    
    % Redefine coordinates of the perpendicular's origin
    X(indx_segment(2,:)+1) = int_points(1,:);
    Y(indx_segment(2,:)+1) = int_points(2,:);
    
    % Take minimum in the edges
    isocy([1 end]) = min([isocy([1 end]); ISOC{o}(2,[1 end])]);
    
%     % Plot (uncomment)
%     plot(isocx/1000,isocy/1000)
%     hold on

    % Display progress
    disp([num2str(o),'/',num2str(max(isoc_above))])
end