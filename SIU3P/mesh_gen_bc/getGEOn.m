function GEO = getGEOn(GCOORD, ELEM2NODE, Point_id, Phases, include_edges)
    % GEO = getGEOn(GCOORD, ELEM2NODE, Point_id, Phases, include_edges)   
    %
    % +++ purpose +++
    % get mesh vertex coordinates for a specific interface matching GEOMETRY.
    % Interface definition here differs from Point_id. Mesh
    % generation (by triangle) removes duplicate coordinates so that 
    % GEOMETRY (and its Geo_id) shared interfaces between two subdomains
    % are promoted by generate_meshGEO() to the higher interfce value, as that 
    % to be surviving in the triangulation node marking. This process is
    % inverted here so that the interface returned by this function is the
    % mesh-analogous of GEOMETRY.
    % This function can, e.g., be used as workhorse to GEOMETRY resampling
    % while avoiding numerical diffussion (by forcing the resampling to be
    % snapped to the specific mesh interface as opossed to linearly
    % interpolated values between GEOMETRY nodes.
    %
    % Author: Javier Garcia-Pintado, MARUM, 2020-03
    %
    % Details: The returned interface coordinates are given in
    %          increasing-x for quasi-horizontal interfaces.
    %
    % This function is specific to rift2ridge2D quasi-horizontal interfaces and should be adapted for 
    % more complex GEOMETRYes (e.g. kinedyn) 
    %
    % 
    % interfid:: interface identifier in Point_id
    %
    % return:
    % intf_coo :: coordinates of the interface sorted from left to right in the domain
    % intf_ids :: indices whithin GCOORD corresponding to intf_coo

    if nargin < 5
        include_edges = false;
    end
    
    gids = unique(Point_id);
    gids(gids==0) = [];                                                    % 0 is for inner nodes
    hids = [1 3:3:max(gids)-1];                                            % e.g. [1 3 6 9] for 3-layer rift2ridge2D
    
    phaseids = unique(Phases);                                          
    
    % edges nodes generated by triangle have higher node order than vertices:
    if ~include_edges
        iinel = 1:3;
    else
        iinel = 1:6;
    end
    EL2NOD = ELEM2NODE(iinel,:);
    nnod = max(max(EL2NOD));
    Point_id = Point_id(1:nnod);
    GCOORD   = GCOORD(:,1:nnod);
    
    GEO = [];
    for i=1:length(gids)
        gid = gids(i);
        if ismember(gid,hids) && gid > 1
            hii = find(gid==hids);                                         % horizontal interface index within hids
            nids = hids(hii:end);                                          % searched node (Point_id) indices
            pids = phaseids(1:hii-1);                                      % searched phase indices                                       
        else
            nids = gid;                                                    % vertical and bottom horizontal interfaces
            pids = max(1,ceil((gid-1)/3));
        end
    
        elsboo = (ismember(Point_id(EL2NOD(1,:)), nids) | ismember(Point_id(EL2NOD(2,:)), nids) | ismember(Point_id(EL2NOD(3,:)), nids)) ...
                 & ismember(Phases, pids);
            
        intf_ids = intersect(find(ismember(Point_id, nids)), EL2NOD(:,elsboo));
        if ismember(gid,hids)
            horizontal = true;
            [~, sortids] = sort(GCOORD(1,intf_ids));
        else
            horizontal = false;
            [~, sortids] = sort(GCOORD(2,intf_ids));
        end
        GEO(i).gid = gid;
        GEO(i).pids = intf_ids(sortids);                                   % identifiers pointers to Point_id
        GEO(i).coo = GCOORD(:, GEO(i).pids);
        GEO(i).n = size(GEO(i).coo,2);
        GEO(i).horizontal = horizontal;
        GEO(i).to  = uint32(zeros(2,GEO(i).n));
    end
    
    ishor =  [GEO.horizontal];
    
    for i=1:length(GEO)                                                    % identifier pointers within GEO struct
        if ~ishor(i)                                                           % only quasi-horizontal interface nodes may point to parent nodes
            continue;
        end
        for k=1:length(GEO)
            if ~ishor(k) || k <= i                                         % only higher order quasi-horizontal interfaces as possible parents
                continue;
            end    
           [iink, kids] = ismember(GEO(i).pids, GEO(k).pids);
           GEO(i).to(1,iink) = kids(iink);
           GEO(i).to(2,iink) = k;
        end
    end
end % function get_interf()
