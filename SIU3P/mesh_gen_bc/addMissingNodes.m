function L = addMissingNodes(L, symmetric, rmwidow)
    % get overlapping transects: this is just considered as 1D x-monotonic in both lines
    % for two partially overlapping polylines, find missing nodes in both and return the extended polylines.
    % If symmetric==false, only L{1} -> L{2} is done and linearly interpolated nodes are introduced in L{2} to remove the "isolated" edge condition 
    %
    % INPUT
    % L    :: CELL, with
    %  L{1} :: [2,n1] polyline, mutable
    %  L{2} :: [2,n2] polyline, mutable for symmetric=TRUE
    % symmetric :: true: any missing node in bot polylines will be dealt with. false: only L{1} is allowed to merge towards L{2}.
    %              That is, the mutable polyline should be L{1}. 
    % rmwidow   :: remove widow nodes in L{1} in non-overlapping L{1}-L{2} transects. Only feasible for symmetric==false
    %
    % Details: This function does not require x-monotony for quasi-horizontal interfaces   
    %         
    
    % Author: Javier García-Pintado, MARUM, 2020
    %
    if nargin < 2
        symmetric = false;
    end
    if nargin < 3
        rmwidow = true;
    end
    
    if rmwidow && symmetric
        error("addMissingNodes:: --ERR001-- widow node removal only feasible for L{1}-only mutable");
    end
    
    if length(L) ~= 2
        error("addMissingNodes:: --ERR002--");
    end

    l2r = [false, false];        % TRUE for left-to-right x input
    for i=1:2
        l2r(i) = L{i}(1,end) - L{i}(1,1) > 0.;              
        if (~l2r(i))
            L{i} = fliplr(L{i});                                           % sort both polylines from left to right
        end
        
        %while true                                                         % sanitising 1: remove isolated x-backward nodes 
        %    x = L{i}(1,:);                                                 %  this may help for function usage within clip2GEOM()
        %    dxsign = [0 sign(diff(x))];
        %    x2del = (dxsign == -1 & [sign(diff(x)) 0] ~= -1);              % isolated backward nodes
        %    if ~any(x2del)                                   
        %        break
        %    end
        %    L{i} = L{i}(:,~x2del);
        %end
        
        %if any(diff(L{i}(1,:)) < 0.)
        %    error("addMissingNodes:: L points x non-monotonic")            % plot(L{i}(1,:), L{i}(2,:), '.-')
        %end
    end
    
    %GEO(1).coo = L{1};  % local GEO for these two polylines
    %GEO(1).gid = 1;     % only local
    %GEO(2).coo = L{2};
    %GEO(2).gid = 2;
    %GEO = findJunctionsGEO(GEO);
    % figure(); plotGEO(GEO)
    % hold on; plot(L{1}(1,GEO(1).junctions)/1000,L{1}(2,GEO(1).junctions)/1000,'o','markersize',15,'color','red')
    % hold on; plot(L{2}(1,GEO(2).junctions)/1000,L{2}(2,GEO(2).junctions)/1000,'o','markersize',20,'color','green')
    while true
        TD = [];
        for i=1:2                                               % hold on; plot(L{2}(1,TD(2).tran==1)/1000,L{2}(2,TD(2).tran==1)/1000,'o')
            if i==1 
                j=2;
            else
                j=1;
            end
            TD(i).inother = ismember(L{i}(2,:), L{j}(2,:));     % LOGICAL vector: nodes existing in the other polyline
            TD(i).tran    = [0,diff(TD(i).inother)];            % INTEGER vector with values in {-1,0,1}: transect definition: 1: start connected transect, -1: start unconnected transect (possibly missing node in L2) 
            TD(i).ntr     = sum(TD(i).tran == 1);               % INTEGER vector: number of common transects
            TD(i).widows  = TD(i).tran == -1 & [TD(i).tran(2:end) 0] == 1; % LOGICAL vector, indicating widow nodes in L{i} 
            TD(i).ini     = find(TD(i).tran == 1);              % INTEGER vector: index of first node in common transects
            TD(i).ins     = find(TD(i).tran == -1);             % INTEGER vector: index of first detached node [-1] in each transect
            TD(i).ncom    = TD(i).ins - TD(i).ini;              % number of common nodes in first transect
        end
        if rmwidow && any(TD(1).widows)                                    % sanitising 2: remove widow nodes in L{1}
            L{1} = L{1}(:,~TD(1).widows);
            continue;
        end
        if TD(1).ntr == TD(2).ntr
            break 
        end
        
        trid = 1;
        while true
            if TD(1).ncom(trid) ~= TD(2).ncom(trid)
                break
            end
            trid = trid + 1;
        end
    
        if TD(2).ncom(trid) - TD(1).ncom(trid) > 0                          % modify [j] line
            i = 1;                                                          % respect
            j = 2;                                                          % mutable
        else
            i = 2;                                                          % respect
            j = 1;                                                          % mutable 
        end
        ncom = min(TD(i).ncom(trid),TD(j).ncom(trid));                       % number of common nodes in this transect beforee splitting apart
        iids = (TD(i).ini(trid)+ncom):(TD(i).ini(trid+1)-1);                 % i nodes not existing in j and to be included in the latter
        jcut = TD(j).ini(trid)+ncom-1;                                       % node before insertion at j
        if symmetric || j == 1
            L{j} = [L{j}(:,1:jcut),  L{i}(:,iids), L{j}(:,(jcut+1):end)];
        else % non-symmetric and j== 2: add a single node by linear interpolation to remove the "isolated" edge in L{1}
            newnode = mean(L{j}(:,jcut:(jcut+1)),2);
            L{j} = [L{j}(:,1:jcut), newnode, L{j}(:,(jcut+1):end)];
        end
    end % while
    
    for i=1:2
        if (~l2r(i))
            L{i} = fliplr(L{i});                                           % flip back flipped polylines
        end 
    end
end % function addMissingNodes()