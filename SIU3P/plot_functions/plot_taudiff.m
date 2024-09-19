% Plots the diferential stress (TAU_xx - TAU_yy)
TAU_diff1 = abs(TAU_xx(Phases==1,:)-TAU_yy(Phases==1,:));
TAU_diff2 = abs(TAU_xx(Phases==2,:)-TAU_yy(Phases==2,:));
TAU_diff3 = abs(TAU_xx(Phases==3,:)-TAU_yy(Phases==3,:));
TAU_diff4 = abs(TAU_xx(Phases==4,:)-TAU_yy(Phases==4,:));

plot(-GIP_y_all(Phases==1,:)/km,TAU_diff1,'.r', ...
    -GIP_y_all(Phases==2,:)/km,TAU_diff2,'.b', ...
    -GIP_y_all(Phases==3,:)/km,TAU_diff3,'.g', ...
    -GIP_y_all(Phases==4,:)/km,TAU_diff4,'.k')

title('Differential stress')
xlabel('Depth [km]')
ylabel('Differential stress [Pa]')


