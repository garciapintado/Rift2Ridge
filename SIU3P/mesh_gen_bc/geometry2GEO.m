function GEO = geometry2GEO(GEOMETRY, Geo_id)
  % +++ purpose +++
  % transform a vector-based GEOMETRY definition into a sructured one GEO.
  % GEO has one element for each unique class in Geo_id. 
  % Each element in GEO contains a [2,nintf] polyline and additional pointer information. 
  % Quasi-horizontal polylines are written from left to right x-coordinates, and vertical interfaces in upward direction
  %
  % No attempt is done here to check monotonicity in quasi-horizontal interfaces. 
  %  
  %
  % It is advised that a GEOMETRY is first filtered by simplifyGEOM()
  % previous to application of this function
  % 
  % Each element (i)
  % has the slots:
  %      .coo         REAL [2,ncoo_i]
  %      .gid         INTEGER. Class code in Geo_id
  %      .horizontal  LOGICAL. Whether the interface is classified as quasi-horizontal
  %      .flipped     LOGICAL. Whether coordinates in GEO(i).coo are x-flipped with respect to those in GEOMETRY
  %      .n           number of coordinates
  %
  %      .to          INTEGER [2,ncoo_i]. 1st row is pointer to index in 2nd row GEO(j). [0 0] indicates no pointer. 

    gids = unique(Geo_id);
    hids = [1 3:3:max(gids)-1];

    GEO = [];
    for i=1:length(gids)
        gid = gids(i);
        GEO(i).gid = gid;
        GEO(i).coo = GEOMETRY(:,Geo_id==gid);
        GEO(i).n = size(GEO(i).coo,2);
        GEO(i).horizontal = false;
        GEO(i).flipped = false;
        GEO(i).to = uint32(zeros(2,GEO(i).n));
        if ismember(gid,hids)
            GEO(i).horizontal = true;
            l2r = GEO(i).coo(1,end) - GEO(i).coo(1,1) > 0.;
            if ~l2r
                GEO(i).coo = fliplr(GEO(i).coo);
                GEO(i).flipped = true;
            end
        else
            b2t = GEO(i).coo(2,end) - GEO(i).coo(2,1) > 0.;
            if ~b2t
                GEO(i).coo = fliplr(GEO(i).coo);
                GEO(i).flipped = true;
            end
        end
    end
GEO = linkGEO(GEO);

end % function geometry2GEO()