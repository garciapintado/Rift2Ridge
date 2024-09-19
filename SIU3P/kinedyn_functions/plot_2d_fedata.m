function plot_2d_fedata(FigNo, GCOORD, EL2NOD, var, Ux, Uz, meshcol, visible, eids)
% Usage: plot_2d_fedata(FigNo, GCOORD, EL2NOD, var, Ux, Uz, meshcol, visible, eids)
% 
% Purpose: Plot any nodal scalar finite element data
%
% Input:
% 1.  FigNo   : [scalar]    : number of Matlab figure
% 2.  GCOORD  : [matrix]    : coordinates of all nodes in mesh
% 3.  EL2NOD  : [matrix]    : finite element connectivity matrix (nnodel x nel)
% 4.  var     : [colvector] : scalar nodal variable field
% 5.  Ux      : [colvector] : horizontal velocity field (optional)
% 6.  Uz      : [colvector] : vertical velocity field (optional)
% 7.  meshcol : [char]      : color of mesh being plotted ontop
%                           ('none' for don't showing mesh)
% 8.  visible : [scalar]    : 0--> Figure is invisible (not shown on screen)
% 9.  eids    : STRING, default to "456" [kinedyn node sorting convention]
%  
% Output:
%   none (creates a Matlab figure)
%
% Part of M2TRI - 2D FINITE ELEMENT MANTLE CONVECTION CODE
%
% Developed by J.Hasenclever & J.Phipps Morgan, 2007-2010
% Email contact: joerg.hasenclever@zmaw.de, jhasenclever@geomar.de
% For numerical methods see online Ph.D. thesis of J.Hasenclever
% (http://www.sub.uni-hamburg.de/opus/volltexte/2010/4873/)
%
% JH March 2011
% JH Dec 2012: Ux and Uz instead of "vectorfield"
%

if nargin < 8
  visible = true;
end   
if nargin < 9
  eids = "456";  % [kinedyn]
end
  
nnodel = size(EL2NOD,1);
nel    = size(EL2NOD,2);
if ismember(nnodel,[3 6 7])
    nvertx = 3;
else
    error(' Cannot identify element type.');
end
nVnod  = max(max(EL2NOD(1:nvertx,:)));

if nargin<7 || isempty(meshcol)
    meshcol ='none';
end
if ~strcmp(meshcol,'none')
    EL2NOD_mesh = EL2NOD(1:nvertx,:);
end

if size(var,1)==nVnod
    Vdata   = var(EL2NOD(1:nvertx,:));
    
elseif size(var,1)==nel
    if size(var,2)>1
        % Assuming values are at integration points
        
%         % Version 1
%         [IP_X, ~] = ip_triangle(size(var,2));
%         [N,~]     = sf_dsf_tri367(IP_X,6,'matrix');
%         Vdata     = zeros(nvertx,nel);
%         for iel=1:nel
%             Dummy        = N'\var(iel,:)';
%             Vdata(:,iel) = Dummy(1:nvertx);
%         end
%         EL2NOD  = EL2NOD(1:nvertx,:);

%         % Version 2
%         % --> create a discontinuous connectivity matrix and mapp from
%         %     integration points to nodes
%         tmp         = GCOORD;
%         GCOORD      = zeros(2,nvertx*nel);
%         GCOORD(1,:) = tmp(1,EL2NOD(1:nvertx,:));
%         GCOORD(2,:) = tmp(2,EL2NOD(1:nvertx,:));
%         EL2NOD = uint32(reshape(1:nvertx*nel,nvertx,nel));
%         var    = ipval_to_nodval_lsq(GCOORD,EL2NOD,var);
%         Vdata  = reshape(var,nvertx,nel);
        
        % Version 3
        Vdata  = ipval_to_nodval_disc(EL2NOD(1:nvertx,:),var)';
        EL2NOD = EL2NOD(1:nvertx,:);
    else
        EL2NOD = EL2NOD(1:nvertx,:);
        Vdata  = repmat(var(:)',nvertx,1);
    end
    
elseif size(var,2)==nel
    EL2NOD = EL2NOD(1:nvertx,:);
    switch size(var,1)
        case 1
            Vdata = repmat(var(:)',nvertx,1);
        case 3
            Vdata = var';
        otherwise
            error('Cannot handle size of variable.');
    end
    
elseif size(var,1)==3*nel
    EL2NOD = EL2NOD(1:nvertx,:);
    Vdata  = reshape(var,3,nel);
    
else
    EL2NOD = trimesh_p2_to_p1(EL2NOD, 0, eids); % PhaseID=0 as foo input
    nel    = size(EL2NOD,2);
    Vdata  = var(EL2NOD(1:nvertx,:));
end

el2V        = reshape((1:nvertx*nel)',nvertx,nel)';
Vcoord(1,:) = GCOORD(1,EL2NOD);
Vcoord(2,:) = GCOORD(2,EL2NOD);

sfigure(FigNo);clf
if  ~visible   
  set(gcf,'Visible','off');
end

patch('faces',el2V,'vertices',Vcoord','facevertexcdata',Vdata(:),...
      'FaceColor','interp','EdgeColor','none');
hold on

if ~strcmp(meshcol,'none')
    nel         = size(EL2NOD_mesh,2);
    el2V        = reshape((1:nvertx*nel)',nvertx,nel)';
    clear Vcoord
    Vcoord(1,:) = GCOORD(1,EL2NOD_mesh);
    Vcoord(2,:) = GCOORD(2,EL2NOD_mesh);
    patch('faces',el2V,'vertices',Vcoord',...
          'Facecolor','none','EdgeColor',meshcol);
end

colormap(jet(100));
colorbar
axis equal
axis([ min(GCOORD(1,:))-0.02*range(GCOORD(1,:)) ...
       max(GCOORD(1,:))+0.02*range(GCOORD(1,:)) ...
       min(GCOORD(2,:))-0.02*range(GCOORD(2,:)) ...
       max(GCOORD(2,:))+0.02*range(GCOORD(2,:)) ]);
hold all

if nargin>5 && ~isempty(Ux) && ~isempty(Uz) 
    nnod  = length(Ux);
    skip  = 27;   % skip vectors
    scale = 0.5; % length of vectors
    quiver(GCOORD(1,1:skip:nnod)',GCOORD(2,1:skip:nnod)',...
           Ux(1:skip:nnod),Uz(1:skip:nnod),scale,'k');
end

drawnow

end % END OF FUNCTION plot_2d_fedata

% #########################################################################
%                              SUB-FUNCTION
% #########################################################################

function r = range(a)
    r = max(a)-min(a);
end
