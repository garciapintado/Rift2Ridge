function [SS] = remesh_random_element(GCOORD,ELEM2NODE,Phases,GCO,E2N,Pha,SS, PHY)

Var_phases = zeros(1,max(Pha));

% TODO this function has to be improved for not a such a great decay on the
% random factor after remeshing (maybe by remeshing using all ips?)

% Find phases with changes in rheology
for n = 1:max(Pha)
    Var_phases(n) = sum(SS.Random(Pha==n,1)==1)~=sum(Pha==n);
end
Var_p = find(Var_phases);

% Coordinates of the old vertexes of the elements
X_rand_old = reshape(GCO(1,E2N(1:3,:)),3,size(E2N,2));
Y_rand_old = reshape(GCO(2,E2N(1:3,:)),3,size(E2N,2));

% Calculate old element centers
CENTERS_old = [sum(X_rand_old)/3; sum(Y_rand_old)/3];

% Coordinates of the new vertexes of the elements
X_rand_new = reshape(GCOORD(1,ELEM2NODE(1:3,:)),3,size(ELEM2NODE,2));
Y_rand_new = reshape(GCOORD(2,ELEM2NODE(1:3,:)),3,size(ELEM2NODE,2));

% Calculate new element centers
CENTERS_new = [sum(X_rand_new)/3; sum(Y_rand_new)/3];

Random_new = ones(size(ELEM2NODE,2),size(SS.Random,2));
Random_new(:,2:end) = 0;
for m = 1:length(SS.RDWS_factorIP)
    RDWSf_new{m} = ones(size(ELEM2NODE,2),size(SS.RDWS_factorIP,2));
    RDWSf_new{m}(:,2:end) = 0;
end

% Find centers of the new triangles in the old elements
Tris_rand = tsearch2(GCO,uint32(E2N(1:3,:)),CENTERS_new);
Ind        = find(Tris_rand==0);

if(~isempty(Ind))
    for i=1:length(Ind)
        [val, Tris_rand(Ind(i))] = min(sqrt((GCO(1,E2N(7,:)) - GCOORD(1,Ind(i))).^2 + (GCO(2,E2N(7,:)) - GCOORD(2,Ind(i))).^2));
    end
end
if(any(Tris_rand==0))
    error('remeshing failed in move_contours');
end

%Random_new = repmat(SS.Random(Tris_rand,1),1,size(SS.Random,2));

% Loop through factors
for n = 1:size(SS.Random,2)
    Remesh_ind = ismember(Pha,Var_p);
    R_ind_new = ismember(Phases,Var_p);
    % Interpolation function
    F_rand = TriScatteredInterp(CENTERS_old(1,Remesh_ind)',CENTERS_old(2,Remesh_ind)', ...
        SS.Random(Remesh_ind,n),'linear');
    for m = 1:length(SS.RDWS_factorIP)                                                                  % each phase
        F_RDWS = TriScatteredInterp(CENTERS_old(1,Remesh_ind)',CENTERS_old(2,Remesh_ind)', ...
            SS.RDWS_factorIP{m}(Remesh_ind,n),'linear');
        RDWSf_new{m}(R_ind_new,n)  = F_RDWS(CENTERS_new(1,R_ind_new),CENTERS_new(2,R_ind_new));
    end
    % Interpolate factors
    Random_new(R_ind_new,n) = F_rand(CENTERS_new(1,R_ind_new),CENTERS_new(2,R_ind_new));
    % Make old factors 1 and 0 to stick to their old value
    Random_new(SS.Random(Tris_rand,n)==1 & R_ind_new',n) = 1;
    Random_new(SS.Random(Tris_rand,n)==0 & R_ind_new',n) = 0;
    Nan_ind = find(isnan(Random_new(:,n)));
    Random_new(Nan_ind,n) = SS.Random(Tris_rand(Nan_ind),n);
    for m = 1:length(SS.RDWS_factorIP)
        Nan_ind = find(isnan(RDWSf_new{m}(:,n)));
        RDWSf_new{m}(Nan_ind,n) = SS.RDWS_factorIP{m}(Tris_rand(Nan_ind),n);
    end
end

% figure(56)
% hold on
% plot_indx = X_rand_old>=-5000 & X_rand_old <=5000 & Y_rand_old>=-10000 & Y_rand_old <=0;
% scatter( X_rand_old(plot_indx)/1000,Y_rand_old(plot_indx)/1000,50,SS.Random(plot_indx),'^','fill')

SS.Random = Random_new;
SS.RDWS_factorIP = RDWSf_new;

% plot_indx = X_rand_new>=-5000 & X_rand_new <=5000 & Y_rand_new>=-10000 & Y_rand_new <=0;
% scatter( X_rand_new(plot_indx)/1000,Y_rand_new(plot_indx)/1000,50,SS.Random(plot_indx),'fill')
% colorbar

% Calculate friction angle ranges for softening functions
SS.Phi1_rand = SS.RDWS_factorIP{1}.*SS.Random.*(PHY.SS.Phi_var)+ ...
    (PHY.SS.Phi(1) - PHY.SS.Phi_var.*(SS.RDWS_factorIP{1}-1/2));
SS.Phi2_rand = repmat(PHY.SS.Phi(2), size(SS.Random));

% Calculate accumulated strain ranges for softening functions
SS.I2_phi1_rand = repmat(PHY.SS.I2_phi(1),size(SS.Random));
SS.I2_phi2_rand = (SS.Phi2_rand-SS.Phi1_rand) ...
    * ((PHY.SS.I2_phi(2)-PHY.SS.I2_phi(1))/(PHY.SS.Phi(2)-PHY.SS.Phi(1))) + PHY.SS.I2_phi(1);

% Calculate dislocation ranges for softening functions
SS.Pef_dis1_rand = (SS.RDWS_factorIP{2}-1).*SS.Random+PHY.SS.Pef_dis(1);
SS.Pef_dis2_rand = repmat(PHY.SS.Pef_dis(2), size(SS.Pef_dis1_rand));

% Calculate accumulated strain ranges for dislocation softening functions
SS.I2_pef_dis1_rand = repmat(PHY.SS.I2_pef_dis(1),size(SS.Random));
SS.I2_pef_dis2_rand = (SS.Pef_dis2_rand-SS.Pef_dis1_rand) ...
    * ((PHY.SS.I2_pef_dis(2)-PHY.SS.I2_pef_dis(1)) ...
      /(PHY.SS.Pef_dis(2)-PHY.SS.Pef_dis(1))) + PHY.SS.I2_pef_dis(1);

% Calculate diffusion ranges for softening functions
SS.Pef_dif1_rand = (SS.RDWS_factorIP{2}-1).*SS.Random+PHY.SS.Pef_dif(1);
SS.Pef_dif2_rand = repmat(PHY.SS.Pef_dif(2),size(SS.Pef_dif1_rand));

% Calculate accumulated strain ranges for diffusion softening functions
SS.I2_pef_dif1_rand = repmat(PHY.SS.I2_pef_dif(1),size(SS.Random));
SS.I2_pef_dif2_rand = (SS.Pef_dif2_rand-SS.Pef_dif1_rand) ...
    * ((PHY.SS.I2_pef_dif(2)-PHY.SS.I2_pef_dif(1)) ...
      /(PHY.SS.Pef_dif(2)-PHY.SS.Pef_dif(1))) + PHY.SS.I2_pef_dif(1);

% % Coordinates of the new vertexes of the elements
% X_rheol_new = reshape(GCOORD(1,ELEM2NODE(1:3,:)),3,size(ELEM2NODE,2));
% Y_rheol_new = reshape(GCOORD(2,ELEM2NODE(1:3,:)),3,size(ELEM2NODE,2));
% 
% % Calculate new element centers
% CENTERS_new = [sum(X_rheol_new)/3; sum(Y_rheol_new)/3];
% 
% % Find centers of the new triangles in the old elements
% Tris_rheol = tsearch2(GCO,uint32(E2N(1:3,:)),CENTERS_new);
% 
% Rheol_var_new = zeros(size(ELEM2NODE,2),size(Rheol_var,2));
% % Find the factor values for the new elements
% for n = 1:size(Rheol_var,2)
%     Rheol_var_new(:,n) = Rheol_var(Tris_rheol,n);
% end
% 
% Rheol_var = Rheol_var_new;
