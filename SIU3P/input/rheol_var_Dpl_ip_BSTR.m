function rheol_var = rheol_var_Dpl_ip_BSTR(Dpl,GCOORD,ELEM2NODE,Phases,nip)
% RHEOL_VAR = RHEOL_VAR_DPL_IP_BSTR(Dpl,GCOORD,ELEM2NODE,PHASES,NIP) is a 
% input function that calculates rheologic variation factors RHEOL_VAR for
% each element in case phases have varying rheologies. In this case BSTR
% stands for Border STRong. This means that the left crust of the model has
% a stronger or different rheology than the rest.

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

DplEl       = zeros(size(ELEM2NODE,2),nip);
rheol_var{1} = zeros(size(ELEM2NODE,2),nip);
rheol_var{2} = zeros(size(ELEM2NODE,2),nip);
%rheol_var(1,Phases~=1 &  Phases~=2) = 1;
rheol_var{1}(Phases~=1 & Phases~=2) = 1;

%==========================================================================
% GEOMETRIES
%==========================================================================
% Coordinates of the vertexes of the elements
DplEl = reshape(Dpl(1,ELEM2NODE(1:nip,:)),nip,size(ELEM2NODE,2))';

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
P1  = repmat((Phases==1)',1,nip);
rheol_var{1}(P1  & DplEl <=0) = 1;
rheol_var{1}(P1  & DplEl >= 0.04) = 0;

%  Dry mantle 
rheol_var{2}(P1  & DplEl >=0.04) = 1;
rheol_var{2}(P1  & DplEl ==0) = 0;

% Linear transition
rheol_var{1}(P1  & DplEl > 0 & DplEl < 0.04) = -(1/0.04).*DplEl(P1  & DplEl > 0 & DplEl < 0.04) +1;
rheol_var{2}(P1  & DplEl > 0 & DplEl < 0.04) = (1/0.04).*DplEl(P1  & DplEl > 0 & DplEl < 0.04) ;

P2  = repmat((Phases==2)',1,nip);
rheol_var{1}(P2) = 1;
rheol_var{2}(P2) = 0;

P3  = repmat((Phases==3)',1,nip);
rheol_var{1}(P3) = 1;
rheol_var{2}(P3) = 0;

% Calculate depletion at integration points
[IP_X,~]    = ip_triangle_m2tri(nip);
[N,~]       = sf_dsf_tri367_N(IP_X,nip,'cell');

rheol_var_n = rheol_var;
for ip = 1:nip
    rheol_var{1}(:,ip) = (N{ip}'*rheol_var_n{1}')';
    rheol_var{2}(:,ip) = (N{ip}'*rheol_var_n{2}')';
end

% STRONGER CRUST IN THE LEFT BOUNDARY
% Input
% -----
width_str = 2e4;
trans     = 1e4;

% Geometry
% --------
% Find coordinates of the integration points
[Gx,Gy] = ip_coord(GCOORD,ELEM2NODE,size(ELEM2NODE,2),nip);
% Boolean for ip in the crust
inC = repmat((Phases==max(Phases))',1,nip);
% 1) Strong
%    ------
% Find coordinate of the left boundary
leftc = min(GCOORD(1,:));
% Find beginning of transition
bt = leftc+width_str;
% Find strong indexes
inSTR = Gx>=leftc & Gx<=bt & inC;
% Change rheological factors
rheol_var{1}(inSTR) = 0;
rheol_var{2}(inSTR) = 1;

% 2) Transition
%    ----------
% Find end of transition
et = leftc+width_str+trans;
% Find transition indexes
inTRA = Gx>bt & Gx<=et & inC;
% Define rheologies
rheol_var{1}(inTRA) = (Gx(inTRA)-bt)/trans;
rheol_var{2}(inTRA) = -(Gx(inTRA)-bt)/trans+1;

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