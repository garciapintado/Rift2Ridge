function rheol_var = rheol_var_craton(GCOORD,ELEM2NODE,Phases)
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
load('Input/CeqB40_GEOM.mat')
% Number of different rheologies per phase (should be the same number as
% number of columns in A, N and Q)
num_rheol = 2;

rheol_var = zeros(num_rheol,size(Phases,2));
rheol_var(1,Phases~=3) = 1;

%==========================================================================
% GEOMETRIES
%==========================================================================
% Coordinates of the vertexes of the elements
X = reshape(GCOORD(1,ELEM2NODE(1:3,:)),3,size(ELEM2NODE,2));
Y = reshape(GCOORD(2,ELEM2NODE(1:3,:)),3,size(ELEM2NODE,2));

% Calculate element centers
CENTERS = [sum(X)/3; sum(Y)/3];

% Transition points
trans1 = GEOM.trans1 + min(GCOORD(1,:));
trans2 = GEOM.trans2 + min(GCOORD(1,:));

%==========================================================================
% CALCULATE RHEOLOGIC VARIATION FACTORS
%==========================================================================
% Cratonic lower crust
rheol_var(1,Phases==3 & CENTERS(1,:)<=trans1) = 1;
% Normal lower crust
rheol_var(2,Phases==3 & CENTERS(1,:)>=trans2) = 1;
% Linear transition
rheol_var(1,Phases==3 & CENTERS(1,:)>trans1 & ...
    CENTERS(1,:)<trans2) = (CENTERS(1,Phases==3 & ...
    CENTERS(1,:)>trans1 & CENTERS(1,:)<trans2)-trans2)* ...
    1/-(trans2-trans1);
rheol_var(2,Phases==3 & CENTERS(1,:)>trans1 & ...
    CENTERS(1,:)<trans2) = (CENTERS(1,Phases==3 & ...
    CENTERS(1,:)>trans1 & CENTERS(1,:)<trans2)-trans1)* ...
    1/(trans2-trans1);

rheol_var = rheol_var';

%==========================================================================
% PLOT
%==========================================================================
subplot(2,1,1)
patch('faces',ELEM2NODE(1:3,:)','vertices',GCOORD', ...
    'facevertexcdata',rheol_var(:,1),'FaceColor','flat')
title('Craton-rheologic factor')
colorbar
shading flat

subplot(2,1,2)
patch('faces',ELEM2NODE(1:3,:)','vertices',GCOORD', ...
    'facevertexcdata',rheol_var(:,2),'FaceColor','flat')
title('Normal-crust-rheologic factor')
colorbar
shading flat