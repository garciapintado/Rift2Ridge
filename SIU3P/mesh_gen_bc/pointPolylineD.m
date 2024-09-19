function anss = pointPolylineD(xy, xyp)
    % +++ purpose +++
    % obtain in 2D Euclidean space the closest point [to] and distance [from] a polyline to a set of point/s
    %
    % INPUT
    % xy  (I)    :: [ncoo,2] matrix of coordinates defining the polyline [one 2D coordinate per row]
    % xyp (I)    :: [p,2] matrix of the external points to estimate distances to the line
    % orthogonal :: LOGICAL. Argument passed to pointLineD(). See description in there.
    %
    % OUTPUT
    % struct with keys:
    % xyc   :: [p,2] coordinates of xyp projected into the segment/line
    % d     :: point-segment/line distance
    % dc    :: distance, along the segment, between the first point in segment and the projection of xyp onto the line
    % cross :: LOGICAL, indicating whether the projected point [xyc] falls within the segment
    % 
    % 
    % 
    % Author: Javier Garcia-Pintado. MARUM, 2021-01
    
    % example - not run
    if 1 > 2
        xy = [[0.5 0.5];
              [1. 1.];                                                     % start segment
              [1. 2.];
              [2. 3.];
              [4. 1.]];                                                    % end      
        xyp = (rand(100,2) + 0.1) * 3.;
        anss = pointPolylineD(xy, xyp);
        
        % mapping of xp into xy-nearest neighbour nodes
        figure(); plot(xyp(:,1),xyp(:,2),'x','color',[.4 .4 .4]); axis equal         % point cloud
        hold on; plot(xy(:,1),xy(:,2),'o-','color',[.1 .1 .5]);            % polyline
        hold on; quiver(xyp(:,1), xyp(:,2), xy(anss.inode_noortho,1) - xyp(:,1), xy(anss.inode_noortho,2) - xyp(:,2), 0, 'color',[.5 .1 .1],'ShowArrowHead','off'); 
      
        % final mapping & display a) distance, b) chainage for start of segment, and c) total chainage
        figure(); plot(xyp(:,1),xyp(:,2),'x','color',[.4 .4 .4]); axis equal         % point cloud
        hold on; plot(xy(:,1),xy(:,2),'o-','color',[.1 .1 .5]);                      % polyline
        hold on; quiver(xyp(:,1), xyp(:,2), anss.xc' - xyp(:,1), anss.yc' - xyp(:,2), 0, 'color',[.5 .1 .1],'ShowArrowHead','off'); 
        labels =  "(d:"+round(anss.dis,2)+",c0:"+round(anss.chain0,2)+",c:" +round(anss.chain0+anss.dc,2)+")"; % [1,p] 
        text(xyp(:,1), xyp(:,2), labels');
    end % example
    
    p = size(xyp,1);
    ncoo = size(xy,1);
    dchain = sqrt(diff(xy(:,1)).^2 + diff(xy(:,2)).^2);
    chainage = [0;cumsum(dchain)];                                         % [ncoo,1]
    
    d     = repelem(NaN,ncoo-1,p); % distance to each polyline segment for each point in x
    cross = d;
    for i=1:(ncoo-1)
        crlst(i) = pointLineD(xy([i i+1],:), xyp, false);                  % [nseg=ncoo-1] struct
        d(i,:) = crlst(i).d;  
        cross(i,:) = crlst(i).cross; % LOGICAL, where the ortogonal to the segment crosses it
    end
    [dsor,dord] = sort(d,1);                                               % ip = 25; dsor(:,ip) == d(dord(:,ip),ip)
    crsort = cross;                                                        % [ncoo-1, p]
    inode = zeros(1,p);
    inode_noortho = max(dord(1:2,:));                                      % [1,np] nearest node: this does not work if nearest is first or last nodes
    ip2dis_st = sqrt((xy(1,1) - xyp(:,1)).^2 + (xy(1,2) - xyp(:,2)).^2);                     % deal with exceptions:
    ip2dis_en = sqrt((xy(end,1) - xyp(:,1)).^2 + (xy(end,2) - xyp(:,2)).^2);                 % "
    ip2dis = sqrt((xy(inode_noortho,1) - xyp(:,1)).^2 + (xy(inode_noortho,2) - xyp(:,2)).^2);% "
    inode_noortho(ip2dis_st < ip2dis) = 1;                                                   % "
    ip2dis = sqrt((xy(inode_noortho,1) - xyp(:,1)).^2 + (xy(inode_noortho,2) - xyp(:,2)).^2);% "
    inode_noortho(ip2dis_en < ip2dis) = ncoo;                                                % "
 
    for ip=1:p
        crsort(:,ip) = cross(dord(:,ip),ip); % LOGICAL, whether the segment is orthogonally crossed from nearest to farthest segment
        if any(crsort(:,ip))
            inode(ip) = dord(find(crsort(:,ip) == 1,1),ip);                % if crossing normal segment is available
        else
            inode(ip) = inode_noortho(ip);                                 % nearest location is a node
        end
    end
   
    % ansn = ["inode","x0","y0","chain0","xc","yc","dc","dis"];
    anss.inode_noortho = inode_noortho;
    anss.inode = inode;                                                    % [1,p]
    anss.x0 = xy(inode,1)';                                                % [1,p]
    anss.y0 = xy(inode,2)';                                                % [1,p]
    anss.chain0 = chainage(inode)';                                        % [1,p]
    anss.xc  = zeros(1,p);
    anss.yc  = zeros(1,p);
    anss.dc  = zeros(1,p);
    anss.dis = zeros(1,p);
    for ip=1:p
        iseg = min(inode(ip),ncoo-1);
        anss.xc(ip)  = crlst(iseg).xyc(ip,1); % [1,p] mapping onto polyline - x
        anss.yc(ip)  = crlst(iseg).xyc(ip,2); % [1,p] mapping onto polyline - y 
        anss.dc(ip)  = crlst(iseg).dc(ip);    % [1,p] local segment chainage for the mapped point
        anss.dis(ip) = crlst(iseg).d(ip);     % [1,p] distance from point to mapped point
    end
    ip2dis = sqrt((xy(inode_noortho,1) - xyp(:,1)).^2 + (xy(inode_noortho,2) - xyp(:,2)).^2)'; % [1,p]
    cln = anss.dis >= ip2dis;                 % [1,p] LOGICAL; True for closer node [i.e. some node is closer than any segment or point on node]
    cln(isnan(cln)) = true;
    anss.inode(cln)  = inode_noortho(cln);
    anss.dis(cln)    = ip2dis(cln);
    anss.xc(cln)     = xy(anss.inode(cln),1);
    anss.yc(cln)     = xy(anss.inode(cln),2);
    anss.x0(cln)     = anss.xc(cln);
    anss.y0(cln)     = anss.yc(cln);
    anss.chain0(cln) = chainage(anss.inode(cln));
    anss.dc(cln)     = 0.;
   
  
end % end function pointPolylineD
