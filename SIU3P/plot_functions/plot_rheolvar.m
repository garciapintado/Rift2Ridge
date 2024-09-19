% PLOT RHEOLOGY VARIATION FACTORS

%--------------------------------------------------------------------------
% Function written by Miguel Andres-Martinez, PhD student at Royal Holloway
% University of London, 08-10-2015. Email: mandresmartinez87@gmail.com
%--------------------------------------------------------------------------

MESH.GCOORD = GCOORD/1000;
MESH.EL2NOD = ELEM2NODE;

subplot(2,1,1)
plot_val(MESH,RHEOL.var{1},nel,SOLVER.nip_stress)
title('Wet mantle-rheologic factor')
colorbar
colormap('jet')
xlabel('Distance [km]')
ylabel('Depth [km]')

subplot(2,1,2)
plot_val(MESH,RHEOL.var{2},nel,SOLVER.nip_stress)
title('Dry-mantle--rheologic factor')
colorbar
colormap('jet')
xlabel('Distance [km]')
ylabel('Depth [km]')