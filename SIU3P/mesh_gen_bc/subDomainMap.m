function [EL2SEL, ] = subDomainMap(EL2NOD, GCOORD, elboo)
  % EL2NOD :: REAL [nnodel,nel], with nnodel \in {6,7}. The central node is not used in the refining.
  % GCOORD :: REAL [2,nnod] 
  % elboo  :: LOGICAL [1,nel] vector indicating which elements are to be included in the downscaled subdomain
  

  % downscaling
  nels = sum(elboo);
  EL2SEL = [find(elboo);                            % 1st row: element indices in parent domain
            1:nels];                                % 2nd row: element indices in children domain
  GCOO3 = 

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