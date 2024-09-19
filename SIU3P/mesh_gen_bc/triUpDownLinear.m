function [] = triRefiner(EL2NOD, GCOORD, Point_id, rr, astriangle)
  % EL2NOD :: REAL [nnodel,nel], with nnodel \in {6,7}. The central node is not used in the refining.
  % GCOORD :: REAL [2,nnod] 

  % number of refinement recursions, each recursion divide a triangle into 4 child triangles. Default to 1 [i.e. 4 child triangles]
  % this function return a refinded triangulation and information to downscale mesh variables

  if nargin < 4
      rr = 1;
  end
  if nargin < 5
      astriangle = true;
  end

  if astriangle           % 3                   triangle convention [e.g. Miguedyn]. 
    trf = [1 6 5;         % | \                 start by local node 1. Counterclockwise triangles and nodes 
           6 2 4;         % 5  4
           5 4 3;         % |    \
           6 4 5]         % 1--6--2
  else                    % 3                   Burges [e.g. Kinedyn]
    trf = [1 4 6;         % | \      
           4 2 5;         % 6  5
           6 5 3;         % |   \
           4 5 6]         % 1--4--2
  end
  trf = trf';

  tri{1}.EL2NOD = EL2NOD;
  tri{1}.GCOO   = GCOORD;
  for ir=1:rr
      nel = length(trir{1}.EL2NOD);                           % == number of columns
      tri{ir}.EL2SEL = reshape(1:4*nel,4,nel);                % element to subelement mapping
      tri{ir+1}.EL2NOD = reshape(tri{ir}.EL2NOD(trf,:),3,[]); % [3,nel*4]
      ncoo = max(max(tri{2}.EL2NOD));                         % number of coordinates in ir [in vertices + in edges] used in ir+1 
      tri{ir+1}.GCOO =  tri{1}.GCOO(:,ncoo);
      tri{ir+1}.GCOO = [tri{ir+1}.GCOO, TODO: augment with new edges   
  end
  % figure(); plot(GCOORD(1,EL2NOD(:,1)), GCOORD(2,EL2NOD(:,1)),'x')
  % hold on; plot(GCOORD(1,EL2NOD(1,1)), GCOORD(2,EL2NOD(1,1)),'o')
  % plot(GCOORD(1,tri{2}.EL2NOD(:,1)),GCOORD(2,tri{2}.EL2NOD(:,1),'+','markersize',20);

end