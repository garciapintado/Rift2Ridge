function ISO = getIsolines(xy, z, at, onlyone, make_monotonic, nrows, ncols)
  % function ISO = getIsolines(xy, z, at, onlyone)
  % +++ purpose +++
  % approximate calculation of isolines by gridding and contourc()
  %
  % xy      :: REAL [ncoo,2]
  % z       :: REAL [ncoo,1] 
  % at      :: REAL, scalar, isoline levels
  % onlyone :: LOGICAL, false as default. If true, only the isoline with higher number of
  %            nodes is returned at ech level
  % make_monotonic :: LOGICAL, default to false. If true, imposes x-monotony in contour coordinates
  %            by [1~m hardcoded] lateral displacements in x coordinates 
  % nrows   :: INTEGER, OPT. Number of rows in the interpolated grid 
  % ncols   :: INTEGER, OPT. Number of colums in the interpolated grid
  %
  % OUTPUT
  %
  % ISO :: struct with isolines. Each element in the struct contains:
  %    .at  : isoline level
  %    .coo : [2,:] isoline nodes
  %
  % Details: If not given, the number of columns and rows is automatically
  % select so that the shorted domain direction in xy coordinates has 200
  % elements, and the longest is scaled to approximate a regularly squared
  % gridding
  %
  
  % Javier Garcia-Pintado, MARUM, 2020
  
  if nargin < 4
      onlyone = false;
  end
  if nargin < 5
      make_monotonic = false;
  end    
  if nargin == 6
      error("if OPT argument 'nrows' is given, 'ncols' is also required")
  end
  
  if diff(size(xy)) > 1
      xy = xy';                                                            % as [x,y]
  end
  
  if nargin < 6
      dxdy = diff(minmax(xy')');
      nrc = round(200 * dxdy/min(dxdy));
      nrows = nrc(1);
      ncols = nrc(2);
  end
  xi = linspace(min(xy(:,1)),max(xy(:,1)), nrows);
  yi = linspace(min(xy(:,2)),max(xy(:,2)), ncols);
  [XI,YI] = meshgrid(xi,yi);
  ZI = scatteredInterpolant(xy(:,1),xy(:,2), z(:));
  
  %[cout,hout] = tricontour(GCOORD',ELEM2NODE(1:3,:)',Temp',[40,45]);
  
  ats = at;
  if length(ats) == 1
      ats = [at at];                                                       % required by contourc()
  end
  isolines = contourc(xi,yi,ZI(XI,YI), ats); % [2,:] matrix, where columns with [1,i] in at indicate that 'i' is a header of a polyline, and [2,i] indicates the number of nodes in this polyline 
  
  iheads  = find(ismember(isolines(1,:),at));                               % header columns
  levels  = isolines(1,iheads);                                            % [npoly]
  ncoos   = isolines(2,iheads);
  istarts = iheads + 1;                                                    % [npoly]
  iends   = iheads + isolines(2,iheads);                                   % [npoly]
  
  ISO = [];
  j = 0;                                                                   % index in ISO
  for i=1:length(istarts)
      j = length(ISO) + 1;                                                 % append
      if i > 1 
         k = find(ismember([ISO.at], levels(i)));  
         if ~isempty(k) && onlyone 
             if ncoos(i) <= ISO(k).n
                 continue
             else
                 j = k;                                                    % substitute rather than append
             end
         end
      end
      ISO(j).at = levels(i);
      ISO(j).coo = isolines(:,istarts(i):iends(i));
      ISO(j).n = ncoos(i);
      if ISO(j).coo(1,end) < ISO(j).coo(1,1) 
         ISO(j).coo = fliplr(ISO(j).coo);                                  % sort from left to right
      end
      ISO(j).horizontal = true;
  end
  if make_monotonic
      ISO = makeMonotonicGEO(ISO,1.0);
  end
end % function
