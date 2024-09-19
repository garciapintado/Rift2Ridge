function I2f = F2I2f(F_xx, F_yy, F_xy, F_yx)
     E_xx = (1/2)*(F_xx.^2 + F_yx.^2 - 1);
     E_xy = (1/2)*(F_xx.*F_xy + F_yx.*F_yy);
     E_yy = (1/2)*(F_xy.^2 + F_yy.^2 - 1);
     I2 = [];
     I2f = sqrt((1/2)*(E_xx.^2 + E_yy.^2) + E_xy.^2);
end
