function V_i = interp2d_cubic(GCOORD,EL2NOD,els,LCOORD,V,OPTS)
% Usage: V_i = interp2d_cubic(GCOORD,EL2NOD,els,LCOORD,V,OPTS)
%
% Purpose: Interpolates variables in a triangular finite element mesh using
%          quasi-cubic interpolation functions. Based on code developed by
%          Chao Shi & Jason Phipps Morgan, 2010.
%
% Input
%   GCOORD : [matrix]    : coordinates of all nodes in mesh
%   EL2NOD : [matrix]    : finite element connectivity matrix (nnodel x nel)
%   els    : [rowvector] : element in which each point is located
%   LCOORD : [matrix]    : local coordinates in each element
%   V      : [matrix]    : variables to be interpolated (nnod x nvar, i.e. 
%                          one variable field in each column)
%   OPTS   : [structure] : options for this functions
%
% Output
%   V_i    : [matrix]    : interpolated values (npt x nvar, i.e. one 
%                          interpolated field in each column)
%
% Part of M2TRI - 2D FINITE ELEMENT MANTLE CONVECTION CODE
%
% Developed by J.Hasenclever & J.Phipps Morgan, 2007-2010
% Email contact: joerg.hasenclever@zmaw.de, jhasenclever@geomar.de
% For numerical methods see online Ph.D. thesis of J.Hasenclever
% (http://www.sub.uni-hamburg.de/opus/volltexte/2010/4873/)
%
% JH Jan 2010: vectorized to improve performance
% JH Feb 2011: handles several variable fields at once (speed-up)
% JH Dec 2012: size of EL2NOD changed to dimension (nnodel x nel)
%

% DEFAULT SETTINGS (will be overwritten if options (OPTS) are provided
verbose     = 0; % display profiling data
monotonic   = 1; % cut-off for interpolated values
method_wght = 'dist'; % 'ones' :: uniform weighting
                      % 'dist' :: weighted by distance of node to element center
nelblk      = 10000; % Number of elements in block

if nargin==6
    if isfield(OPTS,'verbose'); verbose=OPTS.verbose; end
    if isfield(OPTS,'monotonic'); monotonic=OPTS.monotonic; end
    if isfield(OPTS,'method_wght'); method_wght=OPTS.method_wght; end
    if isfield(OPTS,'nelblk'); nelblk=OPTS.nelblk; end
end
if isempty(els)
    V_i = [];
    return
end
nnod = size(GCOORD,2);       % number of nodes in mesh
if size(V,2)==nnod; V=V'; end
nvar = size(V,2);
                                                                            t0=tic;
% Calculates the following spacial derivatives:
%  dV1  dV1
% ---- ----
%  dx   dy

lambda1 = 1-sum(LCOORD);
lambda2 = LCOORD(1,:);
lambda3 = LCOORD(2,:);

npt      = length(els);
N       = zeros(9,npt);
N(1,:)  = lambda1;
N(2,:)  = lambda2;
N(3,:)  = lambda3;
N(4,:)  = lambda2.*lambda1;
N(5,:)  = sqrt(2).*lambda2.*lambda3;
N(6,:)  = lambda3.*lambda1;
N(7,:)  = lambda2.*lambda1.*(1-lambda3-2*lambda2);
N(8,:)  = sqrt(2).*lambda2.*lambda3.*(lambda3-lambda2);
N(9,:)  = lambda3.*lambda1.*(1-lambda2-2*lambda3);

xel     = reshape(GCOORD(1,EL2NOD(1:3,els)),3,npt);  % three x values for one triangle in each column
yel     = reshape(GCOORD(2,EL2NOD(1:3,els)),3,npt);

return_wght = 1;
V_i = zeros(npt,nvar);
for ivar=1:nvar % loop over input variables
    B = V(:,ivar);
    
    % Calculate spatial derivatives
    [dBdx,dBdy,dBdx_el,dBdy_el,tmp] = calc_derivatives...
        (GCOORD,EL2NOD,B,nelblk,method_wght,return_wght); % *SUBFUNCTION*
    if return_wght; wght = tmp; return_wght = 0; clear tmp; end
    
    dBdx = dBdx./wght;
    dBdy = dBdy./wght;

    dBdx = reshape(dBdx(EL2NOD(1:3,els)),3,npt) - repmat(dBdx_el(els)',3,1);
    dBdy = reshape(dBdy(EL2NOD(1:3,els)),3,npt) - repmat(dBdy_el(els)',3,1);

    % also need to map the slope information to the xi,yi (xi, eta) axes
    slope_xi      = zeros(3,npt);
    slope_yi      = zeros(3,npt);
    slope_xi(1,:) = (-xel(1,:)+xel(2,:)).*dBdx(1,:)+(-yel(1,:)+yel(2,:)).*dBdy(1,:);
    slope_xi(2,:) = (-xel(1,:)+xel(2,:)).*dBdx(2,:)+(-yel(1,:)+yel(2,:)).*dBdy(2,:);
    slope_xi(3,:) = (-xel(1,:)+xel(2,:)).*dBdx(3,:)+(-yel(1,:)+yel(2,:)).*dBdy(3,:);
    slope_yi(1,:) = (-xel(1,:)+xel(3,:)).*dBdx(1,:)+(-yel(1,:)+yel(3,:)).*dBdy(1,:);
    slope_yi(2,:) = (-xel(1,:)+xel(3,:)).*dBdx(2,:)+(-yel(1,:)+yel(3,:)).*dBdy(2,:);
    slope_yi(3,:) = (-xel(1,:)+xel(3,:)).*dBdx(3,:)+(-yel(1,:)+yel(3,:)).*dBdy(3,:);

    dBdSC = sqrt(1/2)*slope_xi(3,:) - sqrt(1/2)*slope_yi(3,:); %dT/dS at point C (slope along side BC at C(0,1)):  dT/dS=dT/dxi*dxi/dS + dT/deta*deta/dS
    dBdSB = sqrt(1/2)*slope_xi(2,:) - sqrt(1/2)*slope_yi(2,:); %dT/dS at point B (slope along side BC at B(1,0))   the 'S+' is from C(0,1) to B(1,0)

    slope_info      = zeros(6,npt);
    slope_info(1,:) = (slope_xi(1,:)-slope_xi(2,:))/2; % this is to feed to N4
    slope_info(2,:) = (dBdSC-dBdSB)/2;                 % N5
    slope_info(3,:) = (slope_yi(1,:)-slope_yi(3,:))/2; % N6
    slope_info(4,:) = (slope_xi(1,:)+slope_xi(2,:))/2; % N7
    slope_info(5,:) = (dBdSC+dBdSB)/2;                 % N8
    slope_info(6,:) = (slope_yi(1,:)+slope_yi(3,:))/2; % N9

    B_i = sum( N.* [B(EL2NOD(1:3,els)) ; slope_info] )';

    if monotonic
        B_elnods = B(EL2NOD);
        B_max_el = max(B_elnods);
        B_max    = B_max_el(els); clear B_max_el
        B_i      = min(B_i,B_max(:)); clear B_max
        B_min_el = min(B_elnods);
        B_min    = B_min_el(els); clear B_min_el
        B_i      = max(B_i,B_min(:)); clear B_min
    end
    
    V_i(:,ivar) = B_i;
end

if verbose
    fprintf(' Cubic interpolation of %1i variable(s) at %1i points took %5.2sec\n',...
            nvar,npt,toc(t0));
end

end % END OF FUNCTION interp2d_cubic

% #########################################################################
%                              SUB-FUNCTIONS
% #########################################################################

function [dBdx_nod,dBdy_nod,dBdx_el,dBdy_el,wght_nod] = calc_derivatives...
             (GCOORD,EL2NOD,B,nelblk,method_wght,return_wght)

% The following lines are based on the performance-improved assembly described 
% in the paper "MILAMIN: Matlab-based FEM solver for large problems" by
% Dabrowski, Krotkiewski and Schmid; G3, VOL. 9, doi:10.1029/2007GC001719, 2008
nel    = size(EL2NOD,2);
nnod   = length(GCOORD);   % number of nodes in mesh
nelblk = min(nel,nelblk);  % in case blocksize>number of elements
nblo   = ceil(nel/nelblk); % number of blocks

invM1  = zeros(nelblk,3);  % invM1 stores the 1st row of inv(M) of all elements in block
invMx  = zeros(nelblk,3);  % invMx stores the 2nd row...
invMy  = zeros(nelblk,3);  % invMy stores the 3rd row...

dBdx_el  = zeros(nel,1);
dBdy_el  = zeros(nel,1);
dBdx_nod = zeros(nnod,1);
dBdy_nod = zeros(nnod,1);
if return_wght
    wght_nod = zeros(nnod,1);
else
    wght_nod = [];
end

il     = 1;       % 1st element in 1st block
iu     = nelblk;  % last element in 1st block
for ib=1:nblo     % loop over the element blocks
    
    elnods = EL2NOD(:,il:iu)';
    
    M1 = ones(nelblk,3);       % First column of "M" for all elements in block
    Mx = GCOORD(1,elnods);     % 2nd column of "M"... (= x-coordinates)
    Mx = reshape(Mx,nelblk,3); % 1st row = 1st element in block; 2md row...
    My = GCOORD(2,elnods);     % 3rd column of "M"... (= y-coordinates)
    My = reshape(My,nelblk,3); % 1st row = 1st element in block; 2md row...
    
    switch method_wght
        case 'ones'
            wght = ones(nelblk,3);
        case 'dist'
            Cx = repmat(sum(Mx,2)./3,1,3); % x coord of element centers
            Cy = repmat(sum(My,2)./3,1,3); % y coord of element centers
            wght = sqrt( (Mx-Cx).^2 + (My-Cy).^2 );
            wght = 1./wght;
    end
    
    % Determinate of "M" for all elements in block
    detM = M1(:,1).*Mx(:,2).*My(:,3) ...
         + M1(:,2).*Mx(:,3).*My(:,1) ...
         + M1(:,3).*Mx(:,1).*My(:,2) ...
         - M1(:,3).*Mx(:,2).*My(:,1) ...
         - M1(:,1).*Mx(:,3).*My(:,2) ...
         - M1(:,2).*Mx(:,1).*My(:,3);
     
	% 1/det(M) for all elements in block
    invdetM = 1./detM;
     
    % Inversion of M matrices for all elements in block
    % e.g.: invM1 is 1st row of inv(M) of all elements in block
    invM1(:,1) = invdetM .* ( Mx(:,2).*My(:,3) - Mx(:,3).*My(:,2) );
    invM1(:,2) = invdetM .* ( Mx(:,3).*My(:,1) - Mx(:,1).*My(:,3) );
    invM1(:,3) = invdetM .* ( Mx(:,1).*My(:,2) - Mx(:,2).*My(:,1) );
    invMx(:,1) = invdetM .* ( M1(:,3).*My(:,2) - M1(:,2).*My(:,3) );
    invMx(:,2) = invdetM .* ( M1(:,1).*My(:,3) - M1(:,3).*My(:,1) );
    invMx(:,3) = invdetM .* ( M1(:,2).*My(:,1) - M1(:,1).*My(:,2) );
    invMy(:,1) = invdetM .* ( M1(:,2).*Mx(:,3) - M1(:,3).*Mx(:,2) );
    invMy(:,2) = invdetM .* ( M1(:,3).*Mx(:,1) - M1(:,1).*Mx(:,3) );
    invMy(:,3) = invdetM .* ( M1(:,1).*Mx(:,2) - M1(:,2).*Mx(:,1) );
    
    % B_block is the temperature at the 3 nodes of all elements in block
    % ==> has same dimensions as invM1, invMx, invMy: nelblk x 3
    % Thus, we can perform a point-wise multiplication .* and sum the rows
    B_block    = B(elnods);                 % temperature at 3 nodes of all elements in block
    dBdx_elblo = sum( invMx.*B_block , 2 ); % 2nd coefficient == invM(2,:) * T_el
    dBdy_elblo = sum( invMy.*B_block , 2 ); % 3rd coefficient == invM(3,:) * T_el

    tmp1     = wght.*repmat(dBdx_elblo,1,3);
    tmp2     = accumarray( elnods(:), tmp1(:), [nnod 1] );
    dBdx_nod = dBdx_nod + tmp2;
    
    tmp1     = wght.*repmat(dBdy_elblo,1,3);
    tmp2     = accumarray( elnods(:), tmp1(:), [nnod 1] );
    dBdy_nod = dBdy_nod + tmp2;
    
    dBdx_el(il:iu) = dBdx_elblo;
    dBdy_el(il:iu) = dBdy_elblo;
    
    if return_wght
        tmp2     = accumarray( elnods(:), wght(:), [nnod 1]);
        wght_nod = wght_nod + tmp2;
    end
    
    % Go to next block of elements (increase element indices)
    il = il+nelblk;
    if (ib==nblo-1) % last block has to be scaled to the number of elements remaining
        nelblk  = nel-iu;
        invM1    = zeros(nelblk,3);
        invMx    = zeros(nelblk,3);
        invMy    = zeros(nelblk,3);
    end
    iu = iu + nelblk;
end

end % END OF SUBFUNCTION calc_derivatives