function V_II = calc_invariant_II(V_xx,V_zz,V_xz)

% Second invariant
% (actually square root of second invariant)

% V_II = sqrt( 0.25*(V_xx - V_zz).^2 + V_xz.^2 );
V_II = sqrt( 0.5*(V_xx.^2 + V_zz.^2) + V_xz.^2 );

end