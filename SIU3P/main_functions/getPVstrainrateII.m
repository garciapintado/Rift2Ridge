function [ErP, ErV] = getPVstrainrateII(TAU_xx, TAU_yy, TAU_xy, YC, Yield_T2, Gamma, Mu_dif_all, Mu_dis_all)
    % function [ErP, ErV] = getStrainrates(TAU_xx, TAU_yy, TAU_xy, YC, Yield_T2, Gamma, Mu_dif_all, Mu_dis_all)
    %
    % YC :: LOGICAL [nel,nip] indicating where the yield criterion has been met
    % 
    % Javier Garcia-Pintado, MARUM, 2020 - but nothing new

    [nel,nip] = size(YC);
    
    Er_xx = zeros(nel,nip);
    Er_yy = Er_xx;
    Er_xy = Er_xx;
    
    % plastic strain rates
    Er_xx(YC) = 0.5 * Gamma(YC) .* TAU_xx(YC) ./ Yield_T2(YC);
    Er_xy(YC) = 0.5 * Gamma(YC) .* TAU_xy(YC) ./ Yield_T2(YC);  
    Er_yy(YC) = 0.5 * Gamma(YC) .* TAU_yy(YC) ./ Yield_T2(YC);
    ErP = calc_invariant_II(Er_xx, Er_yy, Er_xy);

    % viscous strain rates
    Mu_c_all = Mu_dis_all;
    isdif = ~isnan(Mu_dif_all); % LOGICAL [nel,nip] true for active diffussion
    Mu_c_all(isdif) = (Mu_dis_all(isdif) .* Mu_dif_all(isdif)) ./ (Mu_dis_all(isdif) + Mu_dif_all(isdif));
    ErV = calc_invariant_II(TAU_xx, TAU_yy, TAU_xy) ./ (2*Mu_c_all);
end
