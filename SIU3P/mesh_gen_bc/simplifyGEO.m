function GEO = simplifyGEO(GEO, tolerance, layshift, keepxmids)
    % +++ purpose +++
    % This function removes nodes that are "along-chainage" closer than a given tolerance in a GEO.coo
    % Both topology and surviving coordinates are preserved. 
    % Extreme values in x (peaks, and valleys) are preserved for requested interfaces.
    % Two kinds of Euclidean distances are considered for node removal:
    %   - First, 'layshift' is used as a tolerance to remove points surrounding junction nodes.
    %   - Second, 'tolerance' is a general tolerance that is used for the remaining points.
    %
    %
    % GEO:  IO. Struct, in which each (i) element represents a polyline with (at least) components 
    %    (i).coo        :: REAL [2,ni]    x,z coordinates defining domain and subdomain boundaries
    %    (i).gid        :: INTEGER, polyline identifier
    %    (i).pids       :: INTEGER [1,ni], OPT mapping of GEO(i).coo into mesh coordinates [i indices in GCOO(:,i)]
    % tolerance         :: REAL (m), scalar or [ni] vector of general tolerance
    % layshift          :: REAL (m), scalar tolerance around junction nodes
    % keepxmids         :: GEO indices for which extreme values should be preserved. Default is only the last (topographical) horizontal layer
    %
    %  Details:
    %  - If the optional input GEO.pids exists, it is recalculated and preserved in the output
    %  - The operation splits interfaces into transects, such that contiguous transects share end nodes.
    %  - This function does not require x-monotonically increasing for quasi-horizontal interfaces
    
    % Author: Javier García-Pintado [JGP], MARUM, 2020-03
    %    
    % Last modification:
    %         2021-09-03 JGP
     
    if nargin < 4
        keepxmids = 9;
    end
    
    if isfield(GEO,"pids")
        getpids = true;
    else
        getpids = false;
    end

    ngeo = length(GEO);
    if numel(tolerance) == 1
        tolerance = repelem(tolerance,1,ngeo);
    end
    if numel(tolerance) ~= ngeo
        error("simplifyGEO:: numel(tolerance) ~ ngeo")
    end
    
    GEO = fromToSidesGEO(GEO);       % add .boundids and .boundcoo slots to GEO. Interfaces with a filled .boundcoo slot indicate side ones, which do not include subdomain corner nodes
                                     % and they should never have a joint with any other GEO interface
    % GEO = findJunctionsGEO(GEO);
    for i=1:length(GEO)
        % disp("i :: " + i)
        L{i} = GEO(i).coo;
        TD = [];
        TD(i).nothers = repelem(int32(0), 1, size(GEO(i).coo,2));          % INT32 [ni] initialize: number of other polylines sharing each node in i 
        for j=1:length(GEO)
            if j==i
                continue;
            end
            nj = size(GEO(j).coo,2);
            [jini_x, iids_x] = ismember(GEO(j).coo(1,:), L{i}(1,:));
            [jini_y, iids_y] = ismember(GEO(j).coo(2,:), L{i}(2,:));
            jini = jini_x & jini_y & (iids_x == iids_y);                   % LOGICAL [nj] 
            TD(j).jini = jini;                                             % LOGICAL [nj]
            if sum(jini) == 0
                TD(j).iids = repelem(uint32(0), 1, nj);    
                continue; 
            end
            iids = uint32(iids_x);                                         % INTEGER [nj]
            iids(~jini) = 0;           
            TD(j).iids = iids;                                             % INTEGER [nj], 0 for j nodes not existing in i
            TD(i).nothers(iids(jini)) = TD(i).nothers(iids(jini)) + 1;     % figure(); plotGEO(GEO)
        end                                                                % text(L{i}(1,:)/1000,(L{i}(2,:)+500)/1000, string(TD(i).nothers),'color','red')
    
        % break into transects to preserve topology. Junctions belong to both free and shared transects
        TD(i).tran = [0,diff(TD(i).nothers)];                              % INTEGER [ni]: transect definition: 1: start connected transect, -1: start unconnected transect (possibly missing node in L2
        TD(i).ini  = find(TD(i).tran > 0);                                 % INTEGER [TD(i).ntr]: index of first node in each transects
        if isempty(TD(i).ini) || TD(i).ini(1) > 1                          % text(L{i}(1,:)/1000,(L{i}(2,:)+250)/1000, string(TD(i).tran),'color','blue')
            TD(i).ini = [1 TD(i).ini];
        end
        TD(i).ini = unique([TD(i).ini find(TD(i).tran < 0)-1]);            % INTEGER [ntr,1] junctions   hold on; plot(L{i}(1,TD(i).ini)/1000, L{i}(2,TD(i).ini)/1000, 'o','markersize',15)
        TD(i).end = [TD(i).ini(2:end) size(L{i},2)];                       % INTEGER [ntr,1]             hold on; plot(L{i}(1,TD(i).end)/1000, L{i}(2,TD(i).end)/1000, 's','markersize',15)
        TD(i).ntr  = length(TD(i).ini);                                    % INTEGER: number of common transects
                                                                           % nnt = TD(i).end - TD(i).ini + 1; text(L{i}(1,:)/1000,(L{i}(2,:)+750)/1000, string([repelem(1:TD(i).ntr,nnt-1),TD(i).ntr]),'color','green')
        iidsbL = 1:size(GEO(i).coo,2);                                     % INTEGER [ni] L{i} within line background indices
        iidsaL = iidsbL;                                                   % INTEGER [ni] L{i} within line updated indices

        if contains(GEO(i).class,"side") && TD(i).ntr > 1
            error('junction found for side GEO interface. Check GEO') 
        end
        
        for it=1:TD(i).ntr                                                 % alternating free and common transects
            
            iidsb = TD(i).ini(it):TD(i).end(it);                           % INTEGER [ntrit] this transect indices (background) within GEO(i) coordinates
            if length(iidsb) < 3
                continue
            end
            iidsa = iidsb;                                                 % INTEGER [ntrit]   "        "       "  (to be moved to)
            iidsl = repelem(false,1,length(iidsb));                        % LOGICAL [ntrit] local re-assignation conducted
            
           
            % figure(); plotGEO(GEO)
            % hold on; plot(L{i}(1,iidsb)/1000, L{i}(2,iidsb)/1000,'o-','markersize',20,'color',[0.4 0.7 0.2],'linewidth',1.5)
            
            %trlen = sum(diff2D(L{i}(:,iidsb))); hold on; plot(L{i}(1,iidsb)/1000, L{i}(2,iidsb)/1000,'o-','markersize',25,'color',[0.4 0.4 0.4])
            %
            %if trlen <= tolerance
            %   iidsa(2:end) = iidsb(1);                                     % map into first transect node
            %   iidsl(:) = true;
            %   continue; 
            %end
            slopepre = [0 diff(L{i}(2,iidsb))];                            % REAL [1,ntrit]
            slopepos = [diff(L{i}(2,iidsb)) 0];                            % REAL [1,ntrit]
            isextrem = (slopepre >= 0. & slopepos < 0) | ...               % LOGICAL [1,ntrit] peaks 
                       (slopepre < 0. & slopepos >= 0);                    %                   valleys
            isextrem([1 end]) = false;                                     % junctions bounding the transect not considered for extremes       
            
            if contains(GEO(i).class,"side")                               % only one transect allowed
                chain0 = diff2D([GEO(i).boundcoo(:,1) L{i}(:,1)]);
                chainage = chain0 + [0, cumsum(diff2D(L{i}(:,iidsb)))];         % chainage from bottom to top
                toclear = iidsb(chainage < tolerance(i));
                if ~isempty(toclear)
                    iidsa(toclear) = max(toclear) + 1;
                    iidsl(toclear) = true;
                end
                chain0 = diff2D([GEO(i).boundcoo(:,2) L{i}(:,end)]);
                chainage = chain0 + [0, cumsum(diff2D(fliplr(L{i}(:,iidsb))))]; % chainage from top to bottom
                iidsbflp = flip(iidsb);
                toclear = iidsbflp(chainage < tolerance(i));
                if ~isempty(toclear)
                    iidsa(toclear) = min(toclear) - 1;
                    iidsl(toclear) = true;
                end
            end
            
            if GEO(i).class == "layer" % junction neighbourhood 
                % forward pass
                beyondxm = repelem(false, 1,length(iidsb));                % LOGICAL [1,ntrit] beyond a extreme
                if any(isextrem)              % hold on; plot(L{i}(1,iidsb(isextrem))/1000, L{i}(2,iidsb(isextrem))/1000,'x','markersize',20,'color','green')
                    beyondxm(find(isextrem,1):end) = true; % includes first extreme
                end
                chainage = [0, cumsum(diff2D(L{i}(:,iidsb)))];                           % REAL    [1,ntrit] chainage within this transect
                toini = iidsb(chainage <= layshift & ~beyondxm & iidsb ~= max(iidsb));   % INTEGER [1,ntrit] indices in L{i} to be clipped
                toini = toini(2:end);                                                    % do not consider the same junction
                if ~isempty(toini)
                    iidsa(ismember(iidsb,toini)) = TD(i).ini(it);
                    iidsl(ismember(iidsb,toini)) = true;
                end
                % backward passs
                beyondxm = repelem(false, 1,length(iidsb));
                if any(isextrem)
                    beyondxm(find(flip(isextrem),1):end) = true;
                end
                chainagerev = [0, cumsum(diff2D(fliplr(L{i}(:,iidsb))))];
                iidsflp = flip(iidsb);
                toend = iidsflp(chainagerev <= layshift & ~beyondxm & iidsflp ~= min(iidsb));             % indices in L{i} to be clipped
                toend = toend(2:end);
                if ~isempty(toend)
                    iidsa(ismember(iidsb,toend)) = TD(i).end(it);
                    iidsl(ismember(iidsb,toend)) = true;
                end
            end % if GEO(i).class == "layer"
            
            if length(iidsb) >= 5  % now only consider transects with more than this node number
                for ii=3:length(iidsb)-2    % forward pass
                   if iidsl(ii) || (isextrem(ii) && ismember(i, keepxmids))       % disregard previously reassigned and extrema 
                       continue;
                   end
                   iiobj = find(~iidsl & (1:length(iidsb) < ii),1,'last');                  % along-chainage closest surviving node
                   if ~isempty(iiobj) && diff2D(L{i}(:,iidsb([iiobj ii]))) < tolerance(i)   % assume Euclidean distance is locally sufficient as proxi to along-chainage one
                       iidsa(ii) = iidsb(iiobj);
                       iidsl(ii) = true;
                   end
                end
                for ii=length(iidsb)-2:-1:3 % backward pass
                   if iidsl(ii) || (isextrem(ii) && ismember(i, keepxmids))       % disregard previously reassigned and extrema 
                       continue;
                   end
                   iiobj = find(~iidsl & (1:length(iidsb) > ii),1,'first');                   % along-chainage closest surviving node
                   if ~isempty(iiobj) && diff2D(L{i}(:,iidsb([iiobj ii]))) < tolerance(i)     % assume Euclidean distance is locally sufficient as proxi to along-chainage one
                       iidsa(ii) = iidsb(iiobj);
                       iidsl(ii) = true;
                   end
                end
                iidsa = sort(iidsa);
            end % if length(iidsb) >= 5
            iidsaL(iidsb) = iidsa;
        end % for it                                                             [iidsbL' iidsaL']  
        if isequal(iidsbL, iidsaL)
            continue;
        end
        GEO(i).coo = GEO(i).coo(:,iidsaL);                                     % update coordinates
        if getpids                                                             % hold on; plotGEO(GEO,i,'o')
            GEO(i).pids = GEO(i).pids(iidsaL);
        end
        
        for j=1:length(GEO)                                               
            if j==i
                continue;
            end
            if sum(TD(j).jini) > 0
                iidsj = TD(j).iids(TD(j).jini);                                % i background indices indices pointed to GEO(j)
                GEO(j).coo(:,TD(j).jini) = GEO(i).coo(:,iidsj);                % isequal(GEO(j).coo(:,TD(j).jini), GEO(i).coo(:,iidsj))
                if getpids
                    GEO(j).pids(TD(j).jini) = GEO(i).pids(iidsj);              % isequal(GEO(j).pids(TD(j).jini), GEO(i).pids(iidsj))
                end
            end
        end
        for j=1:length(GEO)                                                    % reduce
            if getpids
                [~,ia,~] = unique(GEO(j).pids,'stable');
            else    
                [~,ia,~] = unique(GEO(j).coo','stable','rows');
            end
            GEO(j).coo = GEO(j).coo(:,ia);
            if getpids
                GEO(j).pids = GEO(j).pids(ia);
            end
        end
    end % for i
    GEO = linkGEO(GEO);
end % function
