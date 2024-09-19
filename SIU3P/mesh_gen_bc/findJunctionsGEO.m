function GEO = findJunctionsGEO(GEO)
    % +++ purpose +++
    % Get nodes in a GEO interface that are junctions with any other GEO interface, and adds a "junction" element to GEO interface
    %
    %
    % GEO:  IO. Struct, in which each (i) element represents a polyline with (at least) components 
    %    (i).coo        :: REAL [2,ni]    x,z coordinates defining domain and subdomain boundaries
    %    (i).gid        :: INTEGER, polyline identifier
    %    (i).pids       :: INTEGER [1,ni], OPT mapping of GEO(i).coo into mesh coordinates [i indices in GCOO(:,i)]
   
    % Author: Javier García-Pintado, MARUM, 2020-09
   
    for i=1:length(GEO)
        L{i} = GEO(i).coo;
        TD = [];
        TD(i).nothers = repelem(int32(0), 1, size(GEO(i).coo,2));          % INT32 [ni] number of other polylines sharing this node with i 
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
        TD(i).ini = unique([TD(i).ini find(TD(i).tran < 0)-1]); 
        TD(i).end = [TD(i).ini(2:end) size(L{i},2)];
        TD(i).ntr  = length(TD(i).ini);                                    % INTEGER: number of common transects
        GEO(i).junctions = TD(i).ini; 
    end % for i
  
end % function
