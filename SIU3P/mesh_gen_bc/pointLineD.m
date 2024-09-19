function xyp2xy = pointLineD(xy, xyp, orthogonal)
    % +++ purpose +++
    % obtain in 2D Euclidean space the closest point [to] and distance [from] a line segment to a set of point/s
    %
    % INPUT
    % xy  (I)    :: [2,2] matrix of coordinates defining the line segment [one 2D coordinate per row]
    % xyp (I)    :: [p,2] matrix of the external points to estimate distances to the line
    % orthogonal :: LOGICAL. If false (default) point coordinates are
    %               projected on segment and distances obtained accordingly. If true, the
    %               segment is consider to define an infinite line, on which projections
    %               are orthogonal as well as calculated distances
    %
    %               betwwen points and segment. If true, distances
    %
    % OUTPUT
    % struct with keys:
    % xyc    :: [p,2] coordinates of xyp projected into the segment/line
    % d     :: point-segment/line distance
    % dc    :: distance, along the segment, between the first point in segment and the projection of xyp onto the line
    % cross :: LOGICAL, indicating whether the projected point [xyc] falls within the segment
    % 
    % Author: Javier Garcia-Pintado. MARUM, 2021-01
    
    % example - not run
    if 1 > 2
        xy = [[1. 1.]; ...                                                 % start segment
              [2. 3.]];                                                    % end      
        xyp = (rand(100,2) + 0.1) * 3.;
        xyp2xy = pointLineD(xy, xyp);
        figure(); plot(xyp(:,1),xyp(:,2),'x','color',[.4 .4 .4]);
        hold on; plot(xy(:,1),xy(:,2),'o-','color',[.1 .1 .5]);
        hold on; plot(xyp2xy.xyc(:,1),xyp2xy.xyc(:,2),'x', 'color',[.5 .1 .1]);
        hold on; quiver(xyp(:,1), xyp(:,2),xyp2xy.xyc(:,1)- xyp(:,1), xyp2xy.xyc(:,2)-xyp(:,2), 0, 'color',[.5 .1 .1],'ShowArrowHead','off');
    end
    
    if nargin < 3
        orthogonal = false;
    end
    
    alpha = atan(abs((xy(2,2)-xy(1,2))/(xy(2,1)-xy(1,1))));                % segment rotation
    if (xy(1,1) > xy(2,1) && xy(1,2) <= xy(2,2))                           % 2nd cuadrant
      alpha = pi - alpha;
    end
    if (xy(1,1) >  xy(2,1) && xy(1,2) > xy(2,2))                           % 3rd cuadrant
      alpha = alpha + pi;
    end
    if (xy(1,1) <= xy(2,1) && xy(1,2) >= xy(2,2))                          % 4th cuadrant
      alpha = 2*pi - alpha;
    end
    if (size(xyp,1) == 2 && size(xyp,2) ~= 2)                              % undetectable for [2,2] matrices
        xyp = xyp';                                                        % [p,2]
    end
    
    p = size(xyp,1);

    % move all into a rotated space so the segment is vertical
    xyr = rotation([xy;xyp], pi/2-alpha);                                  % [p+2,2] counter clockwise rotation of segment & source points
   
    x0 = xyr(1,1);
    y0 = xyr(1,2);
    y1 = xyr(2,2);
    
    xyc = [repelem(x0,p,1), xyr(3:end,2)];                                 % [p,2]   rotated coordinates of the line points nearest to xyp points (normal projections)
    % figure(); plot(xyr(3:end,1),xyr(3:end,2),'x','color',[.4 .4 .4]); hold on; plot(xyr(1:2,1),xyr(1:2,2),'o-','color',[.1 .1 .5]);  
    % hold on; plot(xyc(:,1),xyc(:,2),'x', 'color',[.5 .1 .1]);
    d  = abs(x0 - xyr(3:end,1));                                           % [p,1] points-line distance
    dc = xyc(:,2) - y0;                                                    % [p,1] diferential chainage over [x0,y0] (> 0 if the projection goes in the segment direction)
    cross = false(p,1);      
    cross(xyc(:,2) >= y0 & xyc(:,2) <= y1) = true;                         % 1 = closest to line than to nodes
    
    if ~orthogonal
        xyc(:,2) = min(max(xyc(:,2),y0),y1);    
        d = dist2D(xyc, xyr(3:end,:));
        dc = xyc(:,2) - y0;
    end
        
    % rotate back into physical space
    %xy  = rotation(xyr, alpha-pi/2);
    xyc = rotation(xyc, alpha-pi/2);
    xyp2xy = [];                                                 % list
    xyp2xy.xy    = xy;                                           %round(xy(1:2,:),10);
    xyp2xy.xyp   = xyp;                                          %round(xy(3:(2+p),:),10);
    xyp2xy.xyc   = xyc;
    xyp2xy.d     = d;
    xyp2xy.dc    = dc;
    xyp2xy.cross = cross;
end % end function pointLineD
