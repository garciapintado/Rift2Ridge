function [els,lc,WS] = locate_points_2d(GCOORD, EL2NOD, gX_PT, els, WS)
% Usage: [els,lc,WS] = locate_points_2d(GCOORD,EL2NOD,gX_PT,els,WS)
%
% Purpose: Locate points with coordinates "gX_PT" in a 2D FE mesh
% 
% Input:
%   GCOORD :: nodal coordinates
%   EL2NOD :: connectivity matrix
%   gX_PT  :: coordinates of points to be located
%   els    :: guess for element in which each point is located
%   WS     :: workspace for tsearch2
% Output:
%   els    :: element in which each point is located
%   lc     :: local coordinates of each point in its element
%
% Part of M3TRI 2D convection code family,
% developed by J.Hasenclever & J.Phipps Morgan, 2007-2010
% Email contact: joerg.hasenclever@zmaw.de
% For numerical methods see online Ph.D. thesis
% (http://www.sub.uni-hamburg.de/opus/volltexte/2010/4873/)

% JH Feb 2011
% JH Apr 2015 : added tsearch2
%

verbose = 0;
npt  = size(gX_PT,2);

% If no (useful) vector "els" is provided: Calculate a guess for the
% element in which each each point is located
% ================================================================
if nargin<4 || length(els)~=npt
    iloc = true(1,npt);
    els  = [];
else
    if exist('tsearch2','file')
        iloc = true(1,npt);
    else
        els(els==0)  = 1;
        nel          = size(EL2NOD,2);
        els(els>nel) = nel;
        els          = els(:)';
        
        % Check local coordinates assuming points are in els
        lc   = local_coords_2d(GCOORD,EL2NOD,els,gX_PT);
        iloc = any(lc>1,1) | any(lc<0,1) | sum(lc,1)>1;
        if all(~iloc)
            return
        end
    end
end

if nargin<5
    WS = [];
end
if verbose
    n = length(find(iloc));
    fprintf(' Need to relocate %1i (%6.2f%%) of the points\n',...
            n,100*n/npt);
end

% Search for points that are in wrong elements
nnod = max(max(EL2NOD(1:3,:)));
if exist('tsearch2','file')
    if ~isempty(els)
        [els(iloc),WS] = tsearch2(GCOORD(:,1:nnod),EL2NOD(1:3,:),gX_PT(:,iloc),WS,els(iloc));
    else
        [els,WS] = tsearch2(GCOORD(:,1:nnod),EL2NOD(1:3,:),gX_PT,WS);
    end
else
    els(iloc) = tsearchn(GCOORD(:,1:nnod)',double(EL2NOD(1:3,:)'),gX_PT(:,iloc)');
    WS        = [];
end
iloc  = ~isnan(els) & els>0;
ilost = find(~iloc);

if ~isempty(ilost)
    figure(66666);clf
    trimesh(EL2NOD(1:3,:)',GCOORD(1,:)',GCOORD(2,:)','Color','k');
    hold on
    scatter(gX_PT(1,ilost),gX_PT(2,ilost),5,'g','filled');
    title('Location of points that could not located by tsearch2.');
    iloc  = ~isnan(els) & els>0;
    ilost = find(~iloc);
end

if nargout>1
    % Local coordinates of points in their elements
    lc(:,iloc)  = local_coords_2d(GCOORD,EL2NOD,els(iloc),gX_PT(:,iloc));
    lc(:,ilost) = 0;
end

end % END OF FUNCTION locate_points_2d