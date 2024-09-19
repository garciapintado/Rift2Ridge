function [GEOMETRY, Geo_id] = GEO2geometry(GEO)
  % get a vectorised version of GEO slots {'coo','gid'}
  GEOMETRY = [GEO.coo];
  Geo_id   = repelem([GEO.gid],[GEO.n]);
end

