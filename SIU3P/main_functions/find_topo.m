function [Topography, topoids] = find_topo(GCOORD, ELEM2NODE, Point_id)
% [Topography, topoids] = find_topo(GCOORD, ELEM2NODE,Point_id) 
% +++ purpose +++
% finds the nodes that are part of the topography of a model mesh defined by GCOORD,
% ELEM2NODE and POINT_ID. It returns the x- and y-coordinates of the
% topographic nodes TOPOGRAPHY and the indexes of this nodes TOPO2NODES at
% the global matrix.
%
% INPUT [I]:
% GCOORD :: 
% ELEM2NODE ::
% Point_id ::
% intid :: interface identifier. Defaults to top quasi-horizontal layer 
%
% OUTPUT [O]:
% Topography :: REAL,    [ntop, 2] east to west of domain sorted x,y node coordinates
% topoids :: INTEGER, [ntop, 1] node i indices in GCOORD(:,i) corresponding to Topography
%               such that GCOORD(:,topoids) == Topography


%--------------------------------------------------------------------------
% Function written by Miguel Andres-Martinez,   RHUL  2015-09-05 
% Javier GP, simplified & adapted to meshGEO(), MARUM 2020-05-01
%--------------------------------------------------------------------------
 
  % Surface nodes
  intid = max(Point_id)-1; % top interface
  Surf_nodes = find(Point_id==intid);

  % order interface from left to right
  [TopoX,isort] = sort(GCOORD(1,Surf_nodes));                                % xsorted = x(isort) 
  TopoY = GCOORD(2,Surf_nodes);
  Topography = [TopoX; TopoY(isort)];                    
  topoids = Surf_nodes(isort);
end