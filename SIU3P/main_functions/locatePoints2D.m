function [xy_els, xyloc] = locatePoints2D(GCOORD, EL2NOD, xy, xy_els, WS)

%
% Purpose: Locate points with coordinates "xy" in a 2D FE mesh
% 
% Input:
%   GCOORD :: nodal coordinates
%   EL2NOD :: connectivity matrix
%   xy     :: [2,nnod] coordinates of points to be located
%   xy_els :: [nnod,1] input guess element index for each "xy" point. This
%                      is tested within the function. Indices "i", such that 
%                      xy_els(i)==0 are for xy points for which no first
%                      guess exists
%   WS     :: quadtree struct to speed up tsearch2() node-element allocation
%
% Output:
%   xy_els :: element in which each point is located
%   eX     :: local coordinates of each point in its element
%

  % JH Feb 2011
  % Javier GP, MARUM, 2020, documented, minor simplifications and error cheking

  verbose = 0;
  npt  = size(xy,2);
  nnod = size(GCOORD,2);

  % get a guess for the element in which each each point is located
  if nargin < 4
      xy_els = zeros(npt,1); 
    end 
  if nargin<5
    WS = [];
  end

  if ~isequal(size(xy_els), [npt,1])
      error("locatePoints2D:: length(xy_els) input /= number of input xy coordinates")
  end

  izero = xy_els==0;                                                      % point indices for which no first guess exist
  if any(izero)
      if npt==nnod
          [~,ind]       = unique(EL2NOD);
          xy_els(izero) = uint32(ceil(ind(izero)./size(EL2NOD,1)));
          clear ind
      else
          xy_els(izero) = uint32(ones(1,sum(izero)));
      end
  end
  % Check all local coordinates assuming points are in xy_els
  xyloc = local_coords_2d(GCOORD, EL2NOD, xy_els, xy);                     % REAL    [2,npt] Local coordinates of points in their elements
  lostids = any(xyloc>1,1) | any(xyloc<0,1) | sum(xyloc,1)>1;              % LOGICAL [1,npt] True for xy points for which the corresponding element has not been found
  lostids = lostids(:);                                                    % LOGICAL [npt,1]
  
  if verbose
      nlost = sum(lostids);
      fprintf("Need to relocate %1i (%6.2f%% ) of the points\n",...
              nlost,100*nlost/npt);
  end

  % Search for points that are in wrong elements
  if any(lostids)
    %try
        xy_els(lostids) = tsearch2(GCOORD,EL2NOD(1:3,:),xy(:,lostids), WS);         % figure(); plot_meshF(EL2NOD, GCOORD); hold on; plot(xy(1,find(lostids,1))/1000, xy(2,find(lostids,1))/1000,'o','color',[0.6 0 0])
    %catch %#ok<CTCH>
    %    scale_Lx  = 1/(max(GCOORD(1,:))-min(GCOORD(1,:)));
    %    xy_els(lostids) = tsearchn(scale_Lx*GCOORD',double(EL2NOD'),scale_Lx*xy(:,lostids)');
    %end
    
    xy_els(isnan(xy_els)) = 0;                                             % adds NaN to zeroes directly returned by tsearch2()
    newfoundids = xy_els > 0 & lostids;
    xyloc(:,lostids) = NaN;  
    xyloc(:,newfoundids) = local_coords_2d(GCOORD,EL2NOD,xy_els(newfoundids),xy(:,newfoundids));
    notfound = any(xyloc(:,newfoundids)>1,1) | any(xyloc(:,newfoundids)<0,1) | sum(xyloc(:,newfoundids),1)>1;
    if any(notfound)
        disp("number of non points outside their triangle limits: " + sum(notfound));
        error("locatePoints2D:: newly located points out of respective elements - hint: increase rounding in local_coords_2d()")
    end
  end % if any(lostids)
  
end % END OF FUNCTION locatePoints2D
