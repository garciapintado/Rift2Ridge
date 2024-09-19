function compare_mu(dir1,dir2)
% Compare viscosities of model with the same meshes
Mu1 = load(dir1,'Mu_all');
Mu2 = load(dir2,'Mu_all');
Mu_all = abs(Mu1.Mu_all-Mu2.Mu_all);
disp(['Maximum viscosity difference in the model: ',num2str(max(max(Mu_all)))])
load(dir1,'nel','ELEM2NODE','GCOORD','istep','dt','ma','Point_id', ...
    'Corner_id','Cornin_id','km')
plot_mu
title('Difference in viscosities')