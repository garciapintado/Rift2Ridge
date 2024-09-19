function [remesh,GCOORD,ELEM2NODE,Phases,GEOMETRY,Geo_id,Corner_id,Cornin_id, ...
    Point_id,cont_points,Elsizes,nip,x_min,x_max,ext_rate,ext_erate, ...
    boundary_condition,temp_bc,temp_bc_depth,mode,prec, ...
    Temp,F_xx,F_xy,F_yx,F_yy,TAU_xx_old,TAU_xy_old, ...
    TAU_yy_old,I2,E2all,Mu_all,ndof,nnod,nel,layer_corr,Bc_ind,Bc_val, ...
    Bc_ind_fs,Bc_val_fs,Bct_ind,Bct_val,Geo_mark,Intersect_ID,Sgap, ...
    Rheol_var,SS,Dpl,PRESS_CONT] = ...
    remesh_all2l(remesh,GCOORD,ELEM2NODE,GEOMETRY,Phases,Geo_id,dt, ...
    Corner_id,Cornin_id,cont_points,Elsizes,nip,x_min,x_max,Ylim,ext_rate, ...
    ext_erate,boundary_condition,bc_t,temp_bc,temp_bc_depth,temp_surf,mode,prec, ...
    shift,Temp,F_xx,F_xy,F_yx,F_yy,I2, ...
    TAU_xx_old,TAU_xy_old, TAU_yy_old, Intersect_ID, Sgap, ...
    intersect_s,triangle_path,meshname,Rheol_var,E2all,Mu_all,SS,Ts0, ...
    dTs_dP,dTs_dF,Dpl,PRESS_CONT,REMESH)
% REMESH_ALL generates a new mesh and remeshes the temperature, the
% accumulated gradient of the deformation and the stress

MODEL = "rift2ridge2D";
  
E2N = ELEM2NODE;
GCO = GCOORD;
Pha = Phases;

layer_corr = [0 0];

%remeshing interfaces (TODO)
[GEOMETRY, Geo_id, Intersect_ID,Sgap] =  resample_interf(GCOORD, ...
    GEOMETRY,Geo_id,Corner_id, Cornin_id,cont_points,x_max,x_min,Ylim, ...
    Intersect_ID);

% Resample can make interfaces to intersect so lineintersect or dis_layers 
% needs to be called again
inter_t = tic;
switch intersect_s
    case 'lineintersect'
        [GEOMETRY(1,Geo_id==3), GEOMETRY(2,Geo_id==3), GEOMETRY(1,Geo_id==6), GEOMETRY(2,Geo_id==6), remesh, layer_corr] = lineintersect (GEOMETRY(1,Geo_id==3), GEOMETRY(2,Geo_id==3), GEOMETRY(1,Geo_id==6), GEOMETRY(2,Geo_id==6), remesh, shift);
        Intersect_ID = zeros(1,size(GEOMETRY,2))==1;
    case 'dis_layers'
        [GEOMETRY,Geo_id,Sgap,Intersect_ID,remesh] = ...
            dis_layers(GEOMETRY,Geo_id,Intersect_ID,shift,remesh);
end
toc(inter_t)

%make mesh
[GCOORD,ELEM2NODE, Point_id, Phases] = generate_mesh2l ...
    (GEOMETRY,Geo_id,Elsizes,Sgap,mode,triangle_path,meshname);
nnod    = size(GCOORD,2);
ndof    = size(GCOORD,1);
nel     = size(ELEM2NODE,2);

% Increase the shift and remesh in case the intersections at the interfaces
% are generating very bad triangles
% switch intersect_s
%     case 'dis_layers'
%         shiftc = shift + shift;
%         mesh_c = 1;
%         cc = 1;
%         while mesh_c
%         [qn,amin,amax] = checkmesh(GCOORD,ELEM2NODE);
%             if(qn<0.2 || amin<7 || amax>170) && cc<=2
%                 [GEOMETRY,Geo_id,Sgap,Intersect_ID,remesh] = ...
%                     dis_layers(GEOMETRY,Geo_id,Intersect_ID,shiftc,remesh);
%                 [GCOORD,ELEM2NODE, Point_id, Phases] = generate_mesh ...
%                     (GEOMETRY,Geo_id,Elsizes,Sgap,mode,triangle_path);
%                 shiftc = shiftc + shift;
%                 cc = cc + 1;
%                 nnod = size(GCOORD,2);
%                 ndof = size(GCOORD,1);
%                 nel = size(ELEM2NODE,2);
%             else
%                 mesh_c = 0;
%             end
%         end
% end

%add 7th node
ELEM2NODE(7,:)  = nnod+1:nnod+nel;
GCOORD          = [GCOORD, [...
    mean(reshape(GCOORD(1, ELEM2NODE(1:3,:)), 3, nel));...
    mean(reshape(GCOORD(2, ELEM2NODE(1:3,:)), 3, nel))]];

nnod    = size(GCOORD,2);

%FIND GEOMETRY MARKERS
[~,Geo_mark] = ismember(round(GEOMETRY'*1e3),round(GCOORD'*1e3),'rows');

% CORNER NODES
Layer1_id   = find(Geo_id==1);
Layer6_id   = find(Geo_id==6);
CORNERS     = [GEOMETRY(:,Layer1_id(1)) GEOMETRY(:,Layer1_id(end)) GEOMETRY(:,Layer6_id(1))  GEOMETRY(:,Layer6_id(end))];
Corner_id   = find(ismember((round(prec*GCOORD))', (round(prec*CORNERS))','rows')==1)';
%INNER CORNERS
Layer3_id   = find(Geo_id==3);
CORNERS_IN  = [GEOMETRY(:,Layer3_id(1)) GEOMETRY(:,Layer3_id(end))];
Cornin_id   = find(ismember((round(prec*GCOORD))', (round(prec*CORNERS_IN))','rows')==1)';

if(~isequal(size(Corner_id,2),4) & ~isequal(size(Cornin_id,2),6))
    error('Corners');
end
% Layer6_id   = [find(Point_id==6) Cornin_id(4)];

%reset boundary conditions
[Bc_ind, Bc_val, Point_id, Bc_ind_fs, Bc_val_fs, ext_erate]    = ...
    set_bcs_flow2l(GCOORD, Corner_id, Cornin_id, Point_id, ext_erate, ...
    ext_rate, boundary_condition,bc_t,dt);

% Remesh variables in the ips (accumulated gradient of the deformation,
% accumulated strain invariant, old stresses, viscosities and strain rates)
if strcmp(SS.rI2,'triscatteredinterp')
    % Old ips
    [GIPold_x,GIPold_y] = ip_coord(GCO,E2N,size(E2N,2),6);
    % New ips
    [GIPnew_x,GIPnew_y] = ip_coord(GCOORD,ELEM2NODE,size(ELEM2NODE,2),6);
    Fc      = scatteredInterpolant(GIPold_x(:),GIPold_y(:),I2.c(:), ...
        'natural','linear');
    Itsi_c  = Fc(GIPnew_x,GIPnew_y);
    Fp      = scatteredInterpolant(GIPold_x(:),GIPold_y(:),I2.p(:), ...
        'natural','linear');
    Itsi_p  = Fp(GIPnew_x,GIPnew_y);
end

nnodel_r = 3;

[F_xx, F_xy, F_yx, F_yy,GIP_xF,GIP_yF, I2, TAU_xx_old, ...
    TAU_xy_old, TAU_yy_old,Mu_all,E2all] = remesh_F_TAU_Mu_E2_linear ...
    (GCOORD, ELEM2NODE, GCO, E2N, F_xx, F_xy, F_yx, F_yy, nnodel_r, ...
    TAU_xx_old, TAU_xy_old, TAU_yy_old, Mu_all, E2all, I2, nip);

Mu_all(Mu_all<1e18) = 1e18;
Mu_all(Mu_all>1e24) = 1e24;

if strcmp(SS.rI2,'triscatteredinterp')
    I2.p = Itsi_p;
    I2.c = Itsi_c;
else
    I2.f(I2.f<0) = 0;
    I2.p(I2.p<0) = 0;
    I2.c(I2.c<0) = 0;
end
%Mu_all(Mu_all<0) = -Mu_all(Mu_all<0);

% Remesh random
if SS.rand_s
    [SS] = remesh_random_element(GCOORD,ELEM2NODE,Phases,GCO,E2N,Pha,SS);
end

%Convert mesh to 4-times as many linear elements
EL2NOD4Tri = trimesh_p2_to_p1(E2N(1:6,:), 0, MODEL) ;  %this function(trimesh_p2_to_p1) needs only 6 nodes; PhaseID=0 as foo input
nVnod = max(max(EL2NOD4Tri(1:3,:))); % number of vertex nodes

els = tsearch2(GCO(:,1:nVnod), uint32(EL2NOD4Tri(1:3,:)),GCOORD);
Ind22 = find(els==0); %check of all elements were found
if(~isempty(Ind22))
    for i=1:length(Ind22)
        [val, els(Ind22(i))] = min(sqrt((GCO(1,E2N(7,:)) - GCOORD(1,Ind22(i))).^2 + (GCO(2,E2N(7,:)) - GCOORD(2,Ind22(i))).^2));
    end
end
if(any(els==0))
    error('Remeshing failed in move_contours');
end

%Interpolate continuous pressure to new nodes
LCOORD = local_coords_2d(GCO(:,1:nVnod),EL2NOD4Tri,els,GCOORD);
PRESS_CONT = PRESS_CONT';
PRESS_CONT = interp2d_cubic(GCO,EL2NOD4Tri,els,LCOORD,PRESS_CONT,'nelblk');
% PRESS_CONT(Ind22)=0; % I think this is wrong (Miguel)
PRESS_CONT = PRESS_CONT';

%Interpolate Depletion to new nodes
Dpl = interp2d_cubic(GCO,EL2NOD4Tri,els,LCOORD,Dpl,'nelblk');
Dpl(Ind22)=0.1;
Dpl = Dpl';

Ts_dry = Ts0 + dTs_dP.*PRESS_CONT + dTs_dF.*Dpl;

% Interpolate temperatures
switch REMESH.type
    case 'interp2d_cubic'
        Temp = interp2d_cubic(GCO,EL2NOD4Tri,els,LCOORD,Temp,'nelblk');
        Temp(Ind22) = temp_surf;
    case 'shpf_cuadr'
        % Old element indexes in which the new ips are contained
        trisT = tsearch2(GCO,uint32(E2N(1:3,:)),GCOORD);
        IndT = find(trisT==0);
        
        % Check of all elements were found, and if not it solves the 
        % problem finding the closest element
        if(~isempty(IndT))
            for i=1:length(IndT)
                [~,trisT(IndT(i))] = min(sqrt( ...
                    (GCO(1,E2N(7,:)) - GCOORD(1,IndT(i))).^2 + ...
                    (GCO(2,E2N(7,:)) - GCOORD(2,IndT(i))).^2));
            end
        end
        
        if(any(isnan(trisT)))
            error('remeshing failed in move_contours');
        end
        Temp = remesh_val(trisT,GCO,GCOORD,Temp,E2N);
        Temp(IndT) = temp_surf;
    case 'triscatter'
        FT = scatteredInterpolant(GCO(1,:)',GCO(2,:)',Temp','natural', ...
            'linear');
        Temp = FT(GCOORD(1,:)',GCOORD(2,:)');
end
Temp = Temp';

if Ts0~=0
    %Correct temperatures about old solidus
    correct = find(Temp >Ts_dry);
    Temp(correct) = Ts_dry(correct);
end

[Bct_ind, Bct_val] = set_bcs_temp(GCOORD, Corner_id, Cornin_id, ...
    Point_id, temp_bc, temp_bc_depth,temp_surf);

%==========================================================================
% REMESH FACTORS FOR RHEOLOGIC CHANGES INSIDE PHASES
%==========================================================================
if size(Rheol_var,2) > 1
%     Rheol_var = ...
%         remesh_rheol_var(GCOORD,ELEM2NODE,Phases,GCO,E2N,Pha,Rheol_var);
  Rheol_var = rheol_var_Dpl_ip(Dpl, ELEM2NODE, Phases, nip, GCOORD);
else
    Rheol_var = {ones(nel,nip)};
end

remesh = 0;
