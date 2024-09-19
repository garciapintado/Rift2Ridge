function rheol_var = rheol_var_Dpl_ip(Dpl, ELEM2NODE, Phases, nip, GCOORD, splot)
% rheol_var = rheol_var_Dpl_ip(GCOORD,ELEM2NODE,PHASES) 
% +++ purpose +++
% calculate mantle rheologic variation factors at integration points 
% as a linear function of depletion
%
% Dpl       :: [1,nnod]
% ELEM2NODE :: [nnodel,nel], where nnodel is the number of nodes per element
% Phases    :: [1,nel]
% nip       :: number of integration points 
% splot     :: LOGICAL, OPTIONAL, TRUE for plotting 
%
%==========================================================================
% LOAD AND INITIALIZATION
%==========================================================================
% load('Input/CeqB40_GEOM.mat')
% Number of different rheologies per phase (should be the same number as
% number of columns in A, N and Q)
% this function is adapted for two mantle rheologies (wet and dry) 

%--------------------------------------------------------------------------
% Function written by Miguel Andres-Martinez, PhD student at Royal Holloway
% University of London, 31-10-2014. Email: mandresmartinez87@gmail.com

% Javier García-Pintado - MARUM, 2020.
%                         minor error corrected and modified to use exact interpolation 
%-----------------------------------------------------------------------------
  
  if nargin < 6
    splot = false;
  end

  [nnodel,nel] = size(ELEM2NODE);

  rheol_var{1} = ones(nel, nip);       % wet mantle: default for ~mantle
  rheol_var{2} = zeros(nel, nip);      % dry mantle: default for ~mantle

  [IP_X,~] = ip_triangle(nip);
  [   N,~] = shp_deriv_triangle(IP_X, nnodel);

  Nmat = cell2mat(N');                               % [nnodel,nip]
  DplIp = (Nmat' * Dpl(ELEM2NODE))';                 % [nel,nip]

  %==========================================================================
  % CALCULATE RHEOLOGIC VARIATION FACTORS
  %==========================================================================
    
  P1  = repmat((Phases==1)',1,nip);                                   % LOGICAL [nel,nip] TRUE for mantle
  intran = P1 & DplIp > 0. & DplIp < 0.04;                            % LOGICAL [nel,nip]

  %  Wet mantle
  rheol_var{1}(P1 & DplIp <= 0.)   = 1.;                              % REAL [nel,nip]
  rheol_var{1}(intran) = 1. - (1/0.04) .* DplIp(intran);              % 1. -> 0. as depletion increases up to 0.04     
  rheol_var{1}(P1 & DplIp >= 0.04) = 0.;
  
  %  Dry mantle 
  rheol_var{2}(P1 & DplIp <= 0.)   = 0.;
  rheol_var{2}(intran) =  (1/0.04) .* DplIp(intran);                  % 0. -> 1. as depletion increases up to 0.04
  rheol_var{2}(P1 & DplIp >= 0.04) = 1.;
  
  %rheol_var{1}(~P1) = 1.; % already in initialization
  %rheol_var{2}(~P1) = 0.; % already in initialization
  
  %==========================================================================
  % PLOT
  %==========================================================================
  if splot
      figure() % plot depletion
      patch('Faces',ELEM2NODE(1:3,:)','Vertices',GCOORD', ...
            'FaceVertexCData',Dpl','FaceColor','interp')                % mean(DplEL,2) would do as well
      title('Depletion for each element')
      colorbar

      figure()
      subplot(2,1,1)
      patch('faces',ELEM2NODE(1:3,:)','vertices',GCOORD', ...
            'facevertexcdata',mean(rheol_var{1},2),'FaceColor','flat')
      title('Wet mantle-rheologic factor')
      colorbar
      shading flat
      
      subplot(2,1,2)
      patch('faces',ELEM2NODE(1:3,:)','vertices',GCOORD', ...
          'facevertexcdata',mean(rheol_var{2},2),'FaceColor','flat')
      title('Dry-mantle--rheologic factor')
      colorbar
      shading flat
  end
end % function

