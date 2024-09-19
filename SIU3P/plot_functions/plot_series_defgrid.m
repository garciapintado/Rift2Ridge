
for j = 1:length(TRACKP)
    plot_defgrid(TRACKP{j},E2N_TP,tp_x,tp_y,km,'elements','k')
    pause ()
    clf
end