function [io, xy] = findPointsOutsideMesh(MESH, xy, getNeighbours)
  % locate points outside a mesh
  % io :: indexes, within xy of points outside the mesh
  
  %verbose = 1; t=tic;

  ncoo  = size(xy,2);
  io = true(1,ncoo);

  bdpoly =  MESH.GCOORD(:,MESH.HCObdi3([1:end, 1]));                                 % vertex-node boundary polygon 
                                                                         % figure(); plot_meshF(MESH.EL2NOD, MESH.GCOORD); hold on; 
                                                                         % plot(bdpoly(1,:)/1000,bdpoly(2,:)/1000,'LineWidth',1.5);
  [iin, ion]  = inpolygon(xy(1,:), xy(2,:), bdpoly(1,:), bdpoly(2,:));   % plot(xy(1,iin)/1000,xy(2,iin)/1000,'+','color','blue') 
                                                                         % plot(xy(1,~iin)/1000,xy(2,~iin)/1000,'+','color','red')

  io(iin)  = false; % points within domain
  io(ion)  = false; % points at the boundary of the domain

  if ~getNeighbours
      return;
  end

%% Need to separate the outside points into two groups based on which domain they originated from if required
ind_out     = find(io==1);
pt_out=xy(:,io);
% find the nearest node on boundary
for i=1:size(pt_out,2)
    [node_nearest0,index_node_nearest0]=NearestNode_Boundary(pt_out(:,i),Polypts);
%     node_nearest(:,i)=node_nearest0;
%     index_node_nearest(i)=index_node_nearest0;

% (1) the first simply option, just move the outside point to the nearest
% node on the boundary
xy(:,ind_out(i))=node_nearest0;

% (2) the second option, 
% % but for bathymetry (topography) the follow is not suitable, because top
% % boundary is not a straight line
%     num_BoundaryPolygon_ptID=length(MESH.BoundaryPolygon_ptID);
%     if(index_node_nearest0==num_BoundaryPolygon_ptID)
%         PointID_neighbour1=MESH.PointID(MESH.BoundaryPolygon_ptID(1));
%     else
%         PointID_neighbour1=MESH.PointID(MESH.BoundaryPolygon_ptID(index_node_nearest0+1));
%     end
%     PointID=MESH.PointID(MESH.BoundaryPolygon_ptID(index_node_nearest0));
%     if(index_node_nearest0==1)
%         PointID_neighbour2=MESH.PointID(MESH.BoundaryPolygon_ptID(num_BoundaryPolygon_ptID));
%     else
%         PointID_neighbour2=MESH.PointID(MESH.BoundaryPolygon_ptID(index_node_nearest0-1));
%     end
%     
%     if(PointID==PointID_neighbour1 & PointID==PointID_neighbour2) % the nearest node on the boundary edge, not a corner
% %         for outside point near a edge node, we just change one coordinate
% %         x or y according to the PointID of the nearest node
%         switch PointID
%             case 201   %horizontal wall
%                 gX_PT(2,ind_out(i))=node_nearest0(2);
%             case 203   %horizontal wall
%                 gX_PT(2,ind_out(i))=node_nearest0(2);
%             case 301   %horizontal wall
%                 gX_PT(2,ind_out(i))=node_nearest0(2);
%             case 303   %horizontal wall
%                 gX_PT(2,ind_out(i))=node_nearest0(2);
%             case 202   %vertical wall
%                 gX_PT(1,ind_out(i))=node_nearest0(1);
%             case 204   %vertical wall
%                 gX_PT(1,ind_out(i))=node_nearest0(1);
%             otherwise  %unknown PointID
%                 gX_PT(:,ind_out(i))=node_nearest0;
%         end
%     else %the nearest node is a corner
%         gX_PT(:,ind_out(i))=node_nearest0;
%     end
%     
    
end

% if ~isempty(ind_out)
%     % Find depth of top of smoker pipe
%     z_pipe      = MESH.GCOORD(2, MESH.PointID==105);
%     
%     % Do this only if pipe structure is present
%     if ~isempty(z_pipe)
%         % Top Domain
%         top_out     = find(gX_orig(2, ind_out)>= z_pipe);
%         % Now need to find pts below the bottom and interpolate it back onto the
%         % bottom
%         bt_out      = find(gX_PT(2, ind_out(top_out))<= z_pipe+1e-5);
%         gX_PT(2, ind_out(top_out(bt_out))) = z_pipe;
%         iout(ind_out(top_out(bt_out))) = 0;
%         % Now need to find pts above the top and interpolate it back onto the
%         % top
%         tp_out      = find(gX_PT(2, ind_out(top_out))>= MESH.zmax-1e-5);
%         gX_PT(2, ind_out(top_out(tp_out))) = MESH.zmax;
%         iout(ind_out(top_out(tp_out))) = 0;
%         % Now need to find pts to the right and interpolate it back onto the
%         % right boundary
%         rt_out      = find(gX_PT(1, ind_out(top_out))>= MESH.xmax-1e-5);
%         gX_PT(1, ind_out(top_out(rt_out))) = MESH.xmax;
%         iout(ind_out(top_out(rt_out))) = 0;
%         % Now need to find pts to the left and interpolate it back onto the
%         % left boundary
%         lt_out      = find(gX_PT(1, ind_out(top_out))<= MESH.xmin+1e-5);
%         gX_PT(1, ind_out(top_out(lt_out))) = MESH.xmin;
%         iout(ind_out(top_out(lt_out))) = 0;
%         
%         x_pipe_left = MESH.GCOORD(1, MESH.PointID==101);
%         x_pipe_right = MESH.GCOORD(1, MESH.PointID==102);
%         % Bottom Domain
%         bot_out     = find(gX_orig(2, ind_out)< z_pipe);
%         % Now need to find pts below the bottom and interpolate it back onto the
%         % bottom
%         bt_out      = find(gX_PT(2, ind_out(bot_out))<= MESH.zmin+1e-5);
%         gX_PT(2, ind_out(bot_out(bt_out))) = MESH.zmin;
%         iout(ind_out(bot_out(bt_out))) = 0;
%         % Now need to find pts to the right and interpolate it back onto the
%         % right boundary
%         rt_out      = find(gX_PT(1, ind_out(bot_out))>= x_pipe_right-1e-5);
%         gX_PT(1, ind_out(bot_out(rt_out))) = x_pipe_right;
%         iout(ind_out(bot_out(rt_out))) = 0;
%         % Now need to find pts to the left and interpolate it back onto the
%         % left boundary
%         lt_out      = find(gX_PT(1, ind_out(bot_out))<= x_pipe_left+1e-5);
%         gX_PT(1, ind_out(bot_out(lt_out))) = x_pipe_left;
%         iout(ind_out(bot_out(lt_out))) = 0;
%     else
%         % if no pipe structure
%         % Now need to find pts below the bottom and interpolate it back onto the
%         % bottom
%         bt_out      = find(gX_PT(2, ind_out)<= MESH.zmin+1e-5);
%         gX_PT(2, ind_out(bt_out)) = MESH.zmin;
%         iout(ind_out(bt_out)) = 0;
%         % Now need to find pts above the top and interpolate it back onto the
%         % top
%         tp_out      = find(gX_PT(2, ind_out)>= MESH.zmax-1e-5);
%         gX_PT(2, ind_out(tp_out)) = MESH.zmax;
%         iout(ind_out(tp_out)) = 0;
%         % Now need to find pts to the right and interpolate it back onto the
%         % right boundary
%         rt_out      = find(gX_PT(1, ind_out)>= MESH.xmax-1e-5);
%         gX_PT(1, ind_out(rt_out)) = MESH.xmax;
%         iout(ind_out(rt_out)) = 0;
%         % Now need to find pts to the left and interpolate it back onto the
%         % left boundary
%         lt_out      = find(gX_PT(1, ind_out)<= MESH.xmin+1e-5);
%         gX_PT(1, ind_out(lt_out)) = MESH.xmin;
%         iout(ind_out(lt_out)) = 0;
%     end
% end
% 


%
%
%
% % Fit a rectangular box into the mesh and check these boundaries
% xmin_right = min(MESH.GCOORD(1,ismember(MESH.PointID,[102 103 202])));
% xmax_left  = max(MESH.GCOORD(1,ismember(MESH.PointID,[101 104 204])));
% zmin_top   = min(MESH.GCOORD(2,ismember(MESH.PointID,[103 104 203])));
% zmax_bot   = max(MESH.GCOORD(2,ismember(MESH.PointID,[101 102 201])));
% iin        = gX_PT(1,:)<=xmin_right & gX_PT(1,:)>=xmax_left & ...
%              gX_PT(2,:)<=zmin_top   & gX_PT(2,:)>=zmax_bot;
% iout(iin)  = 0;
%
% % Deal with remaining points (if any left)
% ind        = find(iout);
% if ~isempty(ind)
%     z_PT     = gX_PT(2,ind);
%     zbnd_top = interp1(MESH.xz_top(1,:),MESH.xz_top(2,:),gX_PT(1,ind));
%     zbnd_bot = interp1(MESH.xz_bot(1,:),MESH.xz_bot(2,:),gX_PT(1,ind));
%     iout(ind(z_PT>=zbnd_bot & z_PT<=zbnd_top)) = 0;
%     if nargout==2
%         ztol               = 1e-12 * (zmin_top-zmax_bot);
%         ind2               = z_PT<zbnd_bot;
%         gX_PT(2,ind(ind2)) = zbnd_bot(ind2) + ztol;
%         ind2               = z_PT>zbnd_top;
%         gX_PT(2,ind(ind2)) = zbnd_top(ind2) - ztol;
%     end
% end

if verbose
    nout = length(find(io));
    fprintf(' Points outside domain: %1i of %1i (=%.2f%%); %.3f sec\n',...
        nout,ncoo,100*nout/ncoo,toc(t));
end

end

function [node_nearest,index_node_nearest]=NearestNode_Boundary(pt_out,Polypts_Boundary)
    % find the nearest node on boundary for the point outside 
    % to speed up this process, first set a small rectangular with center of "pt_out", and check which 
    % boundary nodes in the small rectangular, and then calculate distance
    % between "pt_out", select the minimum distance one. If there is not any
    % node in the small rectangular, make it twice bigger and check again.

    half_length_rect=1;
    x0=pt_out(1);      y0=pt_out(2);
    pts_near=[];
%     find some nodes in the small box
    count_smallbox=0;
    while(isempty(pts_near))
        xx=[x0-half_length_rect x0+half_length_rect x0+half_length_rect x0-half_length_rect];
        yy=[y0-half_length_rect y0-half_length_rect y0+half_length_rect y0+half_length_rect];
        [iin,ion]=inpolygon(Polypts_Boundary(1,:),Polypts_Boundary(2,:),xx,yy);
        ind=(iin | ion);
        pts_near=Polypts_Boundary(:,ind);
        half_length_rect=half_length_rect*2;
        count_smallbox=count_smallbox+1;
%         fprintf(1,'the %dth set small box \n',count_smallbox);
    end
%     calculate distance of each node with the outside point "pt_out"
    if(size(pts_near,2)==1)  % if only one node in the small box, it is what we want
        node_nearest=pts_near;
        index_node_nearest=find(ind==1);
    else
        ind_near=find(ind==1);
        dist=sqrt((pts_near(1,:)-x0).^2+(pts_near(2,:)-y0).^2);
        [~,id]=min(dist);
        node_nearest=pts_near(:,id(1));
        index_node_nearest=ind_near(id(1));
    end
end