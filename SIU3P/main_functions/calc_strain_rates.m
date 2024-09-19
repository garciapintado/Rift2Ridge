function [Er_II,Er_xx,Er_zz,Er_xz] = calc_strain_rates(dNUdx,Ue)
% [ER_II,ER_XX,ER_ZZ,ER_XZ] = CALC_STRAIN_RATES(DNUDX,UE) calculates strain
% rates ER_XX, ER_ZZ, and ER_XZ and the square root from the second
% deviatoric of the strain rate ER_II from given global DNUDX derivatives
% of the shape functions and the velocity vector UE. This function works
% either for a single element or for a block of elements.
%
% Author: Joerg Hasenclever, University of Bremen
%         joerg3@uni-bremen.de

if ~iscell(dNUdx)
    % Element-by-element
    if size(Ue,2)==2
        Ue = Ue';
    end
    
    % Strain rates
    Er_xx  = sum(dNUdx(1,:).*Ue(1,:));
    Er_zz  = sum(dNUdx(2,:).*Ue(2,:));
    Er_xz  = 0.5 * sum(sum(dNUdx([2 1],:).*Ue,1));
else
    % Vectorised
    Ux = Ue(:,1:2:end-1);
    Uz = Ue(:,2:2:end);
    
    % Strain rates
    Er_xx  = sum(dNUdx{1}.*Ux,2);
    Er_zz  = sum(dNUdx{2}.*Uz,2);
    Er_xz  = 0.5 * sum(dNUdx{2}.*Ux + dNUdx{1}.*Uz,2);
end
   
    % Strain rate invariant
    Er_II  = sqrt( 0.5*(Er_xx.^2 + Er_zz.^2) + Er_xz.^2 );
   
end