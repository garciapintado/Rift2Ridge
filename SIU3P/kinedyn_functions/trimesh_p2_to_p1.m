function [EL2NOD,PhaseID] = trimesh_p2_to_p1(EL2NOD, PhaseID, eids)
% Usage: [EL2NOD,PhaseID] = trimesh_p2_to_p1(EL2NOD,PhaseID, eids)
% 
% Purpose: Returns new connectivity matrix after splitting each 6-node
%          element into four 3-node elements.
%
% Input:
%   EL2NOD  : [matrix] : finite element connectivity matrix (6 x nel)
%   PhaseID : [vector] : Phase-ID for each element (1 x nel)
%   eids    : STRING in {"456","645"}representing counterclockwise edge node indices
%
% Output:
%   EL2NOD  : [matrix] : finite element connectivity matrix (3 x 4*nel)
%   PhaseID : [vector] : Phase-ID for each element (1 x 4*nel)
%
% Part of M2TRI - 2D FINITE ELEMENT MANTLE CONVECTION CODE
%
% Developed by J.Hasenclever & J.Phipps Morgan, 2007-2010
% Email contact: joerg.hasenclever@zmaw.de, jhasenclever@geomar.de
% For numerical methods see online Ph.D. thesis of J.Hasenclever
% (http://www.sub.uni-hamburg.de/opus/volltexte/2010/4873/)
%
% JH Dec 2012
%
% JGP: 2019-10-25 adapted to include both "kinedyn" and "rift2ridge2D" node sorting input
%                 Defaults to eids == "456" [kinedyn]
  
% MODEL: kinedyn  
% Each quadratic 6-node element is split into 4 linear elements
% I.e.: element #1 has nodes 1 4 6
%       element #2 has nodes 4 2 5
%       element #3 has nodes 6 5 3
%       element #4 has nodes 4 5 6
%  3
%  | \
%  6   5
%  |     \
%  1 - 4 - 2

% This line performs the split of ALL quadratic elements in the
% mesh at once:

% MODEL: rift2ridge
% Each quadratic 6-node element is split into 4 linear elements
% I.e.: element #1 has nodes 1 6 5
%       element #2 has nodes 6 2 4
%       element #3 has nodes 5 4 3
%       element #4 has nodes 6 4 5
% 
%  3
%  | \
%  5   4
%  |     \
%  1 - 6 - 2  
%  
  
  if nargin < 3
    eids = "456";      % kinedyn convention
  end

  switch eids
    case "456"         % "kinedyn"
      EL2NOD = reshape(EL2NOD([1 4 6 ...
			       4 2 5 ...
			       6 5 3 ...
			       4 5 6],:),3,[]);
    case "645"         % rift2ridge 
      EL2NOD = reshape(EL2NOD([1 6 5 ...
			       6 2 4 ...
			       5 4 3 ...
             4 5 6],:),3,[]);
    otherwise
      error("edge node sorting not considered")
  end
  
  if nargin >= 2
    PhaseID = ones(4,1) * double(PhaseID(:)');
    PhaseID = int32(PhaseID(:)');
  end
  
end % END OF FUNCTION trimesh_p2_to_p1
