function rheol_var = rheol_var_Dpl_new(Dpl,GCOORD,ELEM2NODE,Phases)
% RHEOL_VAR = RHEOL_VAR_CRATON(GCOORD,ELEM2NODE,PHASES) is a input function
% that calculates rheologic variation factors RHEOL_VAR for each element in
% case phases have varying rheologies. 
%
% This function defines lateral varying rheology from a craton on the left
% rheol_var(1,:) to a normal continental crust rheol_var(2,:)

%--------------------------------------------------------------------------
% Function written by Miguel Andres-Martinez, PhD student at Royal Holloway
% University of London, 31-10-2014. Email: mandresmartinez87@gmail.com
%--------------------------------------------------------------------------

%==========================================================================
% LOAD AND INITIALIZATION
%==========================================================================
% load('Input/CeqB40_GEOM.mat')
% Number of different rheologies per phase (should be the same number as
% number of columns in A, N and Q)
num_rheol = 2;

rheol_var = zeros(num_rheol,size(Phases,2));
%rheol_var(1,Phases~=1 &  Phases~=2) = 1;
rheol_var(1,Phases~=1) = 1;

%==========================================================================
% GEOMETRIES
%==========================================================================
% Coordinates of the vertexes of the elements
Dplnew = reshape(Dpl(1,ELEM2NODE(1:3,:)),3,size(ELEM2NODE,2));
% Calculate element centers
DplEl = [sum(Dplnew)/3];

% figure(5)
% patch('faces',ELEM2NODE(1:3,:)','vertices',GCOORD', ...
%     'facevertexcdata',DplEl(:),'FaceColor','flat')
% title('Depletion for each element')
% colorbar
% shading flat


%==========================================================================
% CALCULATE RHEOLOGIC VARIATION FACTORS
%==========================================================================
%  Wet mantle 
rheol_var(1, Phases==1  & DplEl <=0 ) = 1;
rheol_var(1, Phases==1  & DplEl >= 0.04) = 0;

%  Dry mantle 
rheol_var(2, Phases==1  & DplEl >=0.04 ) = 1;
rheol_var(2, Phases==1  & DplEl ==0 ) = 0;

% Linear transition
rheol_var(1,Phases==1  & DplEl > 0 & DplEl < 0.04) = -(1/0.04).*DplEl(Phases==1  & DplEl > 0 & DplEl < 0.04) +1;
rheol_var(2,Phases==1  & DplEl > 0 & DplEl < 0.04) = (1/0.04).*DplEl(Phases==1  & DplEl > 0 & DplEl < 0.04) ;

rheol_var(1,Phases ==3 ) = 1;
rheol_var(2,Phases ==3 ) = 0;


rheol_var = rheol_var';

%==========================================================================
% PLOT
%==========================================================================
% figure(7)
% subplot(2,1,1)
% patch('faces',ELEM2NODE(1:3,:)','vertices',GCOORD', ...
%     'facevertexcdata',rheol_var(:,1),'FaceColor','flat')
% title('Wet mantle-rheologic factor')
% colorbar
% shading flat
% 
% subplot(2,1,2)
% patch('faces',ELEM2NODE(1:3,:)','vertices',GCOORD', ...
%     'facevertexcdata',rheol_var(:,2),'FaceColor','flat')
% title('Dry-mantle--rheologic factor')
% colorbar
% shading flat