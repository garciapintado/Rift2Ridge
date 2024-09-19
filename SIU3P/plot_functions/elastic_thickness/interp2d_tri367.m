function V_i = interp2d_tri367(EL2NOD,els,lc,V)
% Usage: V_i = interp2d_tri367(EL2NOD,els,lc,V)
%
% Purpose: Interpolates variables in a triangular finite element mesh using
%          linear or quadratic interpolation functions.
%
% Input
%   EL2NOD : [matrix]    : finite element connectivity matrix (nnodel x nel)
%   els    : [rowvector] : elements in which will be interpolated
%   lc     : [matrix]    : local coordinates in each element
%   V      : [colvector] : variable field to be interpolated
%
% Output
%   V_i    : [colvector] : interpolated values
%
% Part of M2TRI - 2D FINITE ELEMENT MANTLE CONVECTION CODE
%
% Developed by J.Hasenclever & J.Phipps Morgan, 2007-2010
% Email contact: joerg.hasenclever@zmaw.de, jhasenclever@geomar.de
% For numerical methods see online Ph.D. thesis of J.Hasenclever
% (http://www.sub.uni-hamburg.de/opus/volltexte/2010/4873/)
%
% JH Jan 2011
% JH Dec 2012 : function will interpolate linearly (3-node triangles) and
%               quadratically (6- and 7-node triangles) depending on size
%               of EL2NOD
% JH Jan 2013 : handles NaN in lc
% JH Nov 2014 : calls now the function that is also used in the assembly
%

ind    = ~isnan(lc(1,:)) & ~isnan(lc(2,:));
els    = els(ind);
nnodel = size(EL2NOD,1);

% DShape function values at local coordinates
N      = sf_dsf_tri367(lc(:,ind),nnodel,'matrix');

% Interpolate temperature to points using shape functions
V_i      = nan(length(ind),1);
V_i(ind) = sum(N.*V(EL2NOD(:,els)));

end % END OF FUNCTION interp2d_tri367