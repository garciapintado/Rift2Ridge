function [rheol_var, rheol_var_NS] = ...
    rheol_var_craton_seed(GCOORD,ELEM2NODE,Phases,Rheol_seed)
% RHEOL_VAR = RHEOL_VAR_CRATON_SEED(GCOORD,ELEM2NODE,PHASES,Rheol_seed) is 
% a input function that calculates rheologic variation factors RHEOL_VAR 
% for each element in case phases have varying rheologies. 
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
num_rheol = 3;

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

% Rheol_var without seed
rheol_var_NS = rheol_var;

%==========================================================================
% INCLUDE DEFORMATION SEED
%==========================================================================
if ~isempty(Rheol_seed)
    % Distance of the model nodest to the center of the seed
    Nodes2seed = sqrt((GCOORD(1,ELEM2NODE(7,:))-Rheol_seed(1)).^2 + ...
        (GCOORD(2,ELEM2NODE(7,:))-Rheol_seed(2)).^2);
    % Seed elements
    seedEL = Nodes2seed<Rheol_seed(3)/2;
    % Seed
    rheol_var(seedEL,3) = 1;
    rheol_var(seedEL,1:2) = 0;
end

%==========================================================================
% PLOT
%==========================================================================
RHEOL.var = rheol_var;
%plot_rheolv