% PLOT BOX AND INTERFACES
%   Optional input variables:
%
%               USE                     OPTIONS                 DEFAULTS
%
% color_int     Color of the box        #interfaces x 3 vector  Black
%               and interfaces          with values [0 1]
%
% line_width    Width of the box        #interfaces x 1 vector  1
%               and interfaces          with width values
%
% box_meth      Method to plot the box  'general' valid for     'general'
%                                       models with
%                                       discontinuous layers
%                                       'faster' for faster
%                                       plotting but with erros
%                                       for discontinuous
%                                       layers

%--------------------------------------------------------------------------
% Function written by Miguel Andres-Martinez, PhD student at Royal Holloway
% University of London, 17-02-2015. Email: mandresmartinez87@gmail.com
%--------------------------------------------------------------------------

% 01-07-2016 MA
    % Accelerated plotting
% 16-11-2016 MA
    % Ability to plot irregular boundaries fast
% 04-01-2018 MA
    % Fixed a bug by which some segments were not visible in the plot.
    % Function fully rewritten
% 11-04-2019 MA
    % Added another method for the interface processing so that it can
    % handle all type of geometries and interfaces (for example,
    % discontinuous interfaces for discontinuous layers)

%==========================================================================
% BUILDS THE BOX AND THE INTERFACES FOR THE PLOT
%==========================================================================
% Interface ids
Int_id = 3:3:max(Point_id)-1;

% Define color for the plot in case color_int doesn't exist
if ~exist('color_int','var')
    color_int = zeros(length(Int_id),3);
end

% Define the linewidth for the plot in case line_width doesn't exist
if ~exist('line_width','var')
    line_width = ones(length(Int_id),1);
end

% Define the type of plotting method in case box_meth doesn't exist
if ~exist('box_meth','var')
    box_meth = 'general';
end

switch box_meth
    % Fast method
    case 'fast'
        % Number of layers
        nlay = length(Int_id);
        % Points in the vertex of the elements
        pv = zeros(1,size(Point_id,2));
        pv(1:length(GEOMETRY)) = 1;
        pv = pv==1;
        
        switch nlay
            case 3
                % Box 1
                Point_id(Corner_id) = [1 1 9 9];
                Point_id(Cornin_id) = [3 3 6 6];
                
                Box = [GCOORD(:,ismember(Point_id,[1 2 3 4]) & pv) ...
                    GCOORD(:,Corner_id(1))];
                Box(3,:) = ones(1,size(Box,2));
                
                % Box 2
                Box2 = [GCOORD(:,Cornin_id(1)) ...
                    GCOORD(:,ismember(Point_id,[5 6 7]) & pv) ...
                    GCOORD(:,Cornin_id(2))];
                Box2(3,:) = 2*ones(1,size(Box2,2));
                Box = [Box Box2];
                
                % Box 3
                Box3 = [GCOORD(:,Cornin_id(3)) ...
                    GCOORD(:,ismember(Point_id,[8 9 10]) & pv) ...
                    GCOORD(:,Cornin_id(4))];
                Box3(3,:) = 3*ones(1,size(Box3,2));
                Box = [Box Box3];
            case 2
                % Box 1
                Point_id(Corner_id) = [1 1 6 6];
                Point_id(Cornin_id) = [3 3];
                
                Box = [GCOORD(:,ismember(Point_id,[1 2 3 4]) & pv) ...
                    GCOORD(:,Corner_id(1))];
                Box(3,:) = ones(1,size(Box,2));
                
                % Box 2
                Box2 = [GCOORD(:,Cornin_id(1)) ...
                    GCOORD(:,ismember(Point_id,[5 6 7]) & pv) ...
                    GCOORD(:,Cornin_id(2))];
                Box2(3,:) = 2*ones(1,size(Box2,2));
                Box = [Box Box2];
        end
        
        % Plot
        hold on
        for n = 1:nlay
            plot(Box(1,Box(3,:)==n)/1000,Box(2,Box(3,:)==n)/1000, ...
                'color',color_int(n,:),'LineWidth',line_width(n))
        end
        
        % Reset Point_id in the corners
        Point_id(Corner_id) = -1;
        % Reset Point_id in the interior corners
        switch nlay
            case 3
                Point_id(Cornin_id) = -1;
        end
    
    % Generalized method
    case 'general'
        % Find points in interfaces
        Interface = Point_id~=0;
        
        % TRIANGLE EDGE 1
        % ---------------
        % Load indexes from the triangle vertexes
        EL2NODed = ELEM2NODE([1 2 3],:);
        % Find which points of edge 1 are part of an interface
        isinint = ismember(ELEM2NODE([1 3 5],:),find(Interface));
        % Find in which elements edge 1 is part of the interface
        edgeEL = sum(isinint)==3;
        % Find the indexes of the vertex nodes that are part of edge 1,
        % where edge 1 belongs to an interface
        EDGES = EL2NODed([isinint(1,:); false(1,size(isinint,2)); ...
            isinint(2,:)] & repmat(edgeEL,3,1));
        % Store coordinates of edge 1
        X1 = reshape(GCOORD(1,EDGES),2,size(EDGES,1)/2);
        Y1 = reshape(GCOORD(2,EDGES),2,size(EDGES,1)/2);
        
        % TRIANGLE EDGE 2
        % ---------------
        % Load indexes from the triangle vertexes
        isinint = ismember(ELEM2NODE([2 3 4],:),find(Interface));
        % Find which points of edge 2 are part of an interface
        edgeEL = sum(isinint)==3;
        % Find the indexes of the vertex nodes that are part of edge 2,
        % where edge 2 belongs to an interface
        EDGES = EL2NODed([false(1,size(isinint,2)); isinint(1,:); ...
            isinint(2,:)] & repmat(edgeEL,3,1));
        % Store coordinates of edge 2
        X2 = reshape(GCOORD(1,EDGES),2,size(EDGES,1)/2);
        Y2 = reshape(GCOORD(2,EDGES),2,size(EDGES,1)/2);
        
        % TRIANGLE EDGE 2
        % ---------------
        % Load indexes from the triangle vertexes
        isinint = ismember(ELEM2NODE([1 2 6],:),find(Interface));
        % Find which points of edge 3 are part of an interface
        edgeEL = sum(isinint)==3;
        % Find the indexes of the vertex nodes that are part of edge 3,
        % where edge 3 belongs to an interface
        EDGES = EL2NODed([isinint(1,:); isinint(2,:); ...
            false(1,size(isinint,2))] & repmat(edgeEL,3,1));
        % Store coordinates of edge 3
        X3 = reshape(GCOORD(1,EDGES),2,size(EDGES,1)/2);
        Y3 = reshape(GCOORD(2,EDGES),2,size(EDGES,1)/2);
        
        % Store all the edges together
        X = [X1 X2 X3];
        Y = [Y1 Y2 Y3];
        
        % Plot
        hold on
        plot(X/1000,Y/1000, ...
            'color',color_int(1,:),'LineWidth',line_width(1))
end

% % Nodes at the interfaces
% Nodes_int = find(Point_id~=0);
% % Nodes at the interfaces in the shape of the connectivity matrix
% NINT2EL = ismember(ELEM2NODE(1:6,:),Nodes_int);
% % Elements with 3 nodes at the interfaces
% El_int3n = repmat(sum(NINT2EL)==3,6,1);
% % Coordinates of nodes at the 2-interface-node elements
% EL2N6 = ELEM2NODE(1:6,:);
% xp = GCOORD(1,EL2N6(NINT2EL & El_int3n));
% yp = GCOORD(2,EL2N6(NINT2EL & El_int3n));
% idp = Point_id(EL2N6(NINT2EL & El_int3n));
% 
% nelint3n = sum(El_int3n(1,:));
% 
% % Nan_v = NaN*ones(1,nelint3n);
% 
% % Interfaces without corners
% XP = reshape(xp,3,nelint3n);
% YP = reshape(yp,3,nelint3n);
% IDP = reshape(idp,3,nelint3n);
% 
% % Elements with 5 nodes at the interfaces (include corners)
% El_int5n = sum(NINT2EL)==5;
% % Nodes at the corner in the shape of connectivity matrix
% NCORNER2EL = NINT2EL(4:6,El_int5n);
% % Duplicated boolean
% DUPL_B = NCORNER2EL([1 1 2 2 3 3],:);
% % Vertex to edge
% EDGE2VER = repmat([2 3 1 3 1 2]',1,sum(El_int5n));
% %
% DUPL_el = [find(El_int5n); find(El_int5n)];
% DUPL_el = DUPL_el(:)';
% EE = ELEM2NODE(1:3,DUPL_el);
% 
% OO = reshape(EDGE2VER(DUPL_B),2,2*sum(El_int5n));
% UU = zeros(2,2*sum(El_int5n));
% for n = 1:length(EE)
%     UU(:,n) = EE(OO(:,n),n);
% end
% 
% XCP = reshape(GCOORD(1,UU),2,2*sum(El_int5n));
% YCP = reshape(GCOORD(2,UU),2,2*sum(El_int5n));
% 
% hold on
% 
% % Loop to plot interfaces
% for n = 1:length(Int_id)-1
%     % Edges
%     PID = sum(IDP==Int_id(n))==3;
%     % Corners
%     elc = find(sum(ismember(ELEM2NODE,Cornin_id([n*2-1 n*2]))==1));
%     E2NC = ELEM2NODE(1:6,elc);
%     NC = Point_id(E2NC)==Int_id(n) | Point_id(E2NC)==-1;
%     EC = E2NC(:,sum(NC)==3);
%     EC = EC(NC(:,sum(NC)==3));
%     CCX = reshape(GCOORD(1,EC),3,sum(sum(NC)==3));
%     CCY = reshape(GCOORD(2,EC),3,sum(sum(NC)==3));
%     % Plot
%     plot([XP(:,PID) CCX]/1000,[YP(:,PID) CCY]/1000, ...
%         'color',color_int(n+1,:),'LineWidth',line_width(n+1))
% end
% 
% % Plot the box
% % ------------
% % Edges
% PID = sum(IDP==Int_id(end) | IDP==1)==3;
% % Corners
% elc = find(sum(ismember(ELEM2NODE,Corner_id)==1));
% E2NC = ELEM2NODE(1:6,elc);
% NC = Point_id(E2NC)==1 | Point_id(E2NC)==Int_id(end) | Point_id(E2NC)==-1;
% EC = E2NC(:,sum(NC)==3);
% EC = EC(NC(:,sum(NC)==3));
% CCX = reshape(GCOORD(1,EC),3,sum(sum(NC)==3));
% CCY = reshape(GCOORD(2,EC),3,sum(sum(NC)==3));
% % Plot top and bot
% plot([XP(:,PID) CCX]/1000,[YP(:,PID) CCY]/1000, ...
%     'color',color_int(1,:),'LineWidth',line_width(1))
% 
% % Plot lateral boundaries
% Id_r = 2:3:max(Point_id)-2;
% Id_l = 4:3:max(Point_id);
% if exist('Corner_id','var')
%     % Corners
%     elc = find(sum((ismember(ELEM2NODE,Cornin_id) | ...
%         ismember(ELEM2NODE,Corner_id))==1));
%     E2NC = ELEM2NODE(1:6,elc);
%     NC = ismember(Point_id(E2NC),Id_r) | ismember(Point_id(E2NC),Id_l) | ...
%         Point_id(E2NC)==-1;
%     EC = E2NC(:,sum(NC)==3);
%     EC = EC(NC(:,sum(NC)==3));
%     CCX = reshape(GCOORD(1,EC),3,sum(sum(NC)==3));
%     CCY = reshape(GCOORD(2,EC),3,sum(sum(NC)==3));
% else
%     CCX = [];
%     CCY = [];
% end
% % Edges
% PID = sum(ismember(IDP,Id_r) | ismember(IDP,Id_l))==3;
% plot([XP(:,PID) CCX]/1000,[YP(:,PID) CCY]/1000, ...
%     'color',color_int(1,:),'LineWidth',line_width(1))





% %==========================================================================
% % BUILDS THE BOX AND THE INTERFACES FOR THE PLOT
% %==========================================================================
% % Interface ids
% Int_id = 3:3:max(Point_id)-1;
% 
% % Define color for the plot in case color_int doesn't exist
% if ~exist('color_int','var')
%     color_int = zeros(length(Int_id),3);
% end
% 
% % Define the linewidth for the plot in case line_width doesn't exist
% if ~exist('line_width','var')
%     line_width = ones(length(Int_id),1);
% end
% 
% hold on
% % Include corners in Point_id;
% Point_id_c = Point_id;
% Point_id_c(Corner_id) = [1 1 max(Point_id)-1 max(Point_id)-1];
% if exist('Cornin_id','var')
%     Ci_id = repmat(3:3:size(Cornin_id,2)/2*3,2,1);
%     Point_id_c(Cornin_id) = Ci_id(:)';
% end
% % Index of elements with edges
% EL_EDGE_indx = sum(ismember(Point_id_c(ELEM2NODE(1:3,:)),Int_id))==2;
% % Elements with edges
% EL_EDGE = ELEM2NODE(1:3,EL_EDGE_indx);
% % Nodes in EL_EDGE at the edges
% NODES_EDGE = ismember(Point_id_c(ELEM2NODE(1:3,EL_EDGE_indx)),Int_id);
% % 2-nodes-at-edge elements Point_id_cs
% N2_EDGE_Point_id = reshape(Point_id_c(EL_EDGE(NODES_EDGE)),2, ...
%     sum(EL_EDGE_indx));
% % Indexes for the same edge
% SAME_EDGE = repmat(N2_EDGE_Point_id(1,:)==N2_EDGE_Point_id(2,:),3,1);
% se_t = sum(SAME_EDGE(1,:));
% % Index for coordinates
% ICOORD = reshape(EL_EDGE(NODES_EDGE & SAME_EDGE),2,se_t);
% % Coordinates to plot
% XP = reshape(GCOORD(1,ICOORD),2,se_t); [XP,xs] = sort(XP);
% YP = reshape(GCOORD(2,ICOORD),2,se_t); 
% YYP(1,:) = YP(xs==1); YYP(2,:) = YP(xs==2); YP = YYP; clear YYP
% 
% % Point_id_cs for XP and YP
% PIP = Point_id_c(ICOORD(1,:));
% 
% % Loop to plot interfaces
% for n = 1:length(Int_id)-1
%     plot(XP(:,PIP==Int_id(n))/1000,YP(:,PIP==Int_id(n))/1000, ...
%         'color',color_int(n+1,:),'LineWidth',line_width(n+1))
% end
% 
% % Plot the box
% plot(XP(:,PIP==Int_id(end))/1000,YP(:,PIP==Int_id(end))/1000, ...
%         'color',color_int(1,:),'LineWidth',line_width(1))
% plot(GCOORD(1,Corner_id(1:3))/1000,GCOORD(2,Corner_id(1:3))/1000, ...
%     'color',color_int(1,:),'LineWidth',line_width(1))
% plot(GCOORD(1,Corner_id([4 1]))/1000,GCOORD(2,Corner_id([4 1]))/1000, ...
%     'color',color_int(1,:),'LineWidth',line_width(1))
% 
% hold off

% examples: add point_id to the box
	    % ncoo = size(Point_id,2); text(GCOORD(1,1:ncoo),GCOORD(2,1:ncoo),string(Point_id));
if 1 > 2
    plot_box;
    hold on;
    ids = Point_id==0;
    %  text(GCOORD(1,ids),GCOORD(2,ids),string(Point_id(1,ids)); 
    scatter(GCOORD(1,ids),GCOORD(2,ids));
    ids = Point_id == -1;
    plot(GCOORD(1,ids),GCOORD(2,ids),'.','markersize',20,'color','red');
end

				   
