function Rheol_var = remesh_rheol_var(GCOORD,ELEM2NODE,Phases,GCO,E2N, ...
    Pha,Rheol_var)

Var_phases = zeros(1,max(Pha));

% Find phases with changes in rheology
for n = 1:max(Pha)
    Var_phases(n) = sum(Rheol_var(Pha==n,1)==1)~=sum(Pha==n);
end
Var_p = find(Var_phases);

% Coordinates of the old vertexes of the elements
X_rheol_old = reshape(GCO(1,E2N(1:3,:)),3,size(E2N,2));
Y_rheol_old = reshape(GCO(2,E2N(1:3,:)),3,size(E2N,2));

% Calculate old element centers
CENTERS_old = [sum(X_rheol_old)/3; sum(Y_rheol_old)/3];

% Coordinates of the new vertexes of the elements
X_rheol_new = reshape(GCOORD(1,ELEM2NODE(1:3,:)),3,size(ELEM2NODE,2));
Y_rheol_new = reshape(GCOORD(2,ELEM2NODE(1:3,:)),3,size(ELEM2NODE,2));

% Calculate new element centers
CENTERS_new = [sum(X_rheol_new)/3; sum(Y_rheol_new)/3];

Rheol_var_new = ones(size(ELEM2NODE,2),size(Rheol_var,2));
Rheol_var_new(:,2:end) = 0;

% Find centers of the new triangles in the old elements
Tris_rheol = tsearch2(GCO,uint32(E2N(1:3,:)),CENTERS_new);

Ind22 = find(Tris_rheol==0); %check of all elements were found
if(~isempty(Ind22))
    for i=1:length(Ind22)
        [val, Tris_rheol(Ind22(i))] = min(sqrt((GCO(1,E2N(7,:)) - CENTERS_new(1,Ind22(i))).^2 + (GCO(2,E2N(7,:)) - CENTERS_new(2,Ind22(i))).^2));
    end
end
if(any(Tris_rheol==0))
    error('Remeshing failed in move_contours');
end

% Loop through factors
for n = 1:size(Rheol_var,2)
    Remesh_ind = ismember(Pha,Var_p);
    R_ind_new = ismember(Phases,Var_p);
    % Interpolation function
    F_rheol = TriScatteredInterp(CENTERS_old(1,Remesh_ind)',CENTERS_old(2,Remesh_ind)', ...
        Rheol_var(Remesh_ind,n));
    % Interpolate factors
    Rheol_var_new(R_ind_new,n) = F_rheol(CENTERS_new(1,R_ind_new),CENTERS_new(2,R_ind_new));
    % Make old factors 1 and 0 to stick to their old value
    Rheol_var_new(Rheol_var(Tris_rheol,n)==1 & R_ind_new',n) = 1;
    Rheol_var_new(Rheol_var(Tris_rheol,n)==0 & R_ind_new',n) = 0;
    Nan_ind = find(isnan(Rheol_var_new(:,n)));
    Rheol_var_new(Nan_ind,n) = Rheol_var(Tris_rheol(Nan_ind),n);
end

Rheol_var = Rheol_var_new;

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