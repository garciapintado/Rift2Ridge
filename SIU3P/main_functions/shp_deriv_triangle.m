function  [N, dNdu] = shp_deriv_triangle(ipx, nnodel)
%SHP_DERIV_TRIANGLE Shape functions and their derivatives with respect to local coordinates
%  Supports 3,6, and 7 node triangular elements 

%   Part of MILAMIN: MATLAB-based FEM solver for large problems, Version 1.0
%   Copyright (C) 2007, M. Dabrowski, M. Krotkiewski, D.W. Schmid
%   University of Oslo, Physics of Geological Processes
%   http://milamin.org
%   See License file for terms of use.

nip  = size(ipx,1);
N    = cell(nip,1);
dNdu = cell(nip,1);

for i=1:nip
    eta2 = ipx(i,1);
    eta3 = ipx(i,2);
    eta1 = 1-eta2-eta3;

    switch nnodel
        case 3
            SHP   = [eta1; ...
                     eta2; ...
                     eta3];
            DERIV = [-1 1 0; ...   %w.r.t eta2
                     -1 0 1];      %w.r.t eta3

        case 6                              % in kinedyn node definition this correspond to:
            SHP = [eta1*(2*eta1-1);         % N1 at coordinate (r,s)
                   eta2*(2*eta2-1);         % N2 " 
                   eta3*(2*eta3-1);         % N3 " 
                       4*eta2*eta3;         % N5
                       4*eta1*eta3;         % N6
                       4*eta1*eta2];        % N4
            %        dN1      dN2       dN3       dN5    dN6            dN4  
            DERIV = [1-4*eta1 -1+4*eta2         0 4*eta3 -4*eta3        4*eta1-4*eta2; ...   %w.r.t eta2
                     1-4*eta1         0 -1+4*eta3 4*eta2  4*eta1-4*eta3       -4*eta2];      %w.r.t eta3
             
        case 7                                        % in kinedyn node definition sf_dsf_tri7() this correspond to K.
            SHP = [eta1*(2*eta1-1)+ 3*eta1*eta2*eta3; %    K.N1
                   eta2*(2*eta2-1)+ 3*eta1*eta2*eta3; %    K.N2
                   eta3*(2*eta3-1)+ 3*eta1*eta2*eta3; %    K.N3
                     4*eta2*eta3 - 12*eta1*eta2*eta3; %    K.N5
                     4*eta1*eta3 - 12*eta1*eta2*eta3; %    K.N6
                     4*eta1*eta2 - 12*eta1*eta2*eta3; %    K.N4
                                   27*eta1*eta2*eta3];%    K.N7

            DERIV = [1-4*eta1+3*eta1*eta3-3*eta2*eta3 ... % dN1/dr
                    -1+4*eta2+3*eta1*eta3-3*eta2*eta3 ... % dN2/dr
                              3*eta1*eta3-3*eta2*eta3 ... % ...
                     4*eta3+12*eta2*eta3-12*eta1*eta3 ... %
                    -4*eta3+12*eta2*eta3-12*eta1*eta3 ... %
              4*eta1-4*eta2+12*eta2*eta3-12*eta1*eta3 ... % dN4/dr
                           -27*eta2*eta3+27*eta1*eta3;... % dN7/dr
                     1-4*eta1+3*eta1*eta2-3*eta2*eta3 ... % dN1/ds
                             +3*eta1*eta2-3*eta2*eta3 ... % ...
                    -1+4*eta3+3*eta1*eta2-3*eta2*eta3 ... %
                     4*eta2-12*eta1*eta2+12*eta2*eta3 ... %
              4*eta1-4*eta3-12*eta1*eta2+12*eta2*eta3 ... %
                    -4*eta2-12*eta1*eta2+12*eta2*eta3 ... % dN4/ds
                            27*eta1*eta2-27*eta2*eta3];   % dN7/ds
        otherwise
            error('Unknown element')

    end
	
       N{i} = SHP;
    dNdu{i} = DERIV';
 
end

