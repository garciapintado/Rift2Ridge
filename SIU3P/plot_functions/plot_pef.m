% Diffusion
figure(1)
title('Diffusion')
surf(SS.EP,SS.TP,SS.Pp_dif)
xlabel('Accumulate strain')
ylabel('Temperature')
zlabel('Preexponential factor softening')
hold on
plot3(I2.c,TEMP_IP,Pef_dif,'.r')
axis([min(min(SS.EP)) max(max(SS.EP)) min(min(SS.TP)) max(max(SS.TP)) ...
    min(min(SS.Pp_dif)) max(max(SS.Pp_dif))])

% Dislocation
figure(2)
title('Dislocation')
surf(SS.EP,SS.TP,SS.Pp_dis)
xlabel('Accumulate strain')
ylabel('Temperature')
zlabel('Preexponential factor softening')
hold on
plot3(I2.c,TEMP_IP,Pef_dis,'.r')
axis([min(min(SS.EP)) max(max(SS.EP)) min(min(SS.TP)) max(max(SS.TP)) ...
    min(min(SS.Pp_dis)) max(max(SS.Pp_dis))])