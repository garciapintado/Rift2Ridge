%function [Temp] = ini_temp_conv(Tini_dir,GCOORD,d)
% TEMP = INI_TEMP_CONV(TINI_DIR,GCOORD,TEMP_BC) calculates initial
% temperatures for a model by running the mechanical and temperature
% solvers for steady mechanical boundaty conditions. Convection then starts
% and average geotherms are calculated. The function stops running when the
% average geotherm stabilizes and applies the temperature profile over the
% whole domain. TINI_DIR indicates where the model parameters are saved.
% GCOORD and TEMP_BC are variables included just to be consistent with the
% inputs of other temperature loading functions.

%%TEST FLUID2D_SOLVER
%   Part of MILAMIN: MATLAB-based FEM solver for large problems, Version 1.0
%   Copyright (C) 2007, M. Dabrowski, M. Krotkiewski, D.W. Schmid
%   University of Oslo, Physics of Geological Processes
%   http://milamin.org
%   See License file for terms of use.

%   Edited by Miguel Andres-Martinez, PhD student at Royal Holloway
%   University of London. Earth Sciences Department. This version includes:
%       - Free surface algorithm
%       - Elasticity
%       - Discontinuous layers
%       - Varying rheologies along a phase

%==========================================================================
% INITIALIZATION
%==========================================================================
Tini_dir = 'Input/temp';
load(Tini_dir)

dt = 0.01*ma;

plot_gt = 1;

cont_t = 1;
remesh = 0;
etime = 0;
istep = 1;
ntime = dt;

% Remove strain softening
SS.Phases_phi = [];
SS.Phases_c = [];
SS.Phases_pef_dis = [];
SS.Phases_pef_dif = [];

ext_rate = 0;
bc_v = 'ext_rate'; 
ext_erate = 0;

% Reset boundary conditions
[Bc_ind, Bc_val, Point_id, Bc_ind_fs, Bc_val_fs, ext_erate] = ...
    set_bcs_flow3l(GCOORD, Corner_id, Cornin_id, Point_id, ext_erate, ...
    ext_rate, 'ext_rate');

temp_bc_depth = -500000;
temp_bc = 2150;
temp_bc_ini = 2150;
temp_bc_depth_ini = -500000;
%==========================================================================
% INITIAL TEMPERATURES
%==========================================================================
% Calculates the initial temperature field by solving the thermal
% relaxation of the model
[Bct_ind, Bct_val] = set_bcs_temp_ini(...
    GCOORD, Corner_id, Cornin_id, Point_id, temp_bc_ini, temp_bc_depth_ini, ...
    ini_deformation,sigma_tx,sigma_tz);

Temp = temp_bc*ones(size(GCOORD,2),1);

% Thermal solver
nip =6;
for i=1:5
    Temp = thermal2d_rhot(ELEM2NODE, GCOORD, Phases, K, Rho, Cp, Hp, Hs, Temp, Bct_ind, Bct_val, 100*ma, nip, reorder,Ce);
end

% Temperature boundary conditions
[Bct_ind, Bct_val] = set_bcs_temp(GCOORD, Corner_id, Cornin_id, Point_id, temp_bc, temp_bc_depth);

% plot_gt = 0;
% plot_t
% hold on
% scatter(GCOORD(1,Bct_ind)/1000,GCOORD(2,Bct_ind)/1000,[],Bct_val)
% colorbar
% colormap('jet')

while cont_t
    disp(istep)
    t_it = tic;
%==========================================================================
% REMESHING
%==========================================================================

% Switch for the different methods to deal with interface intersections and
% too thin layers
inter_t = tic;
switch intersect_s
    case 'lineintersect'
        [GEOMETRY(1,Geo_id==3), GEOMETRY(2,Geo_id==3), GEOMETRY(1,Geo_id==6), GEOMETRY(2,Geo_id==6), remesh, layer_corr] = lineintersect (GEOMETRY(1,Geo_id==3), GEOMETRY(2,Geo_id==3), GEOMETRY(1,Geo_id==6), GEOMETRY(2,Geo_id==6), remesh, shift);
        [GEOMETRY(1,Geo_id==6), GEOMETRY(2,Geo_id==6), GEOMETRY(1,Geo_id==9), GEOMETRY(2,Geo_id==9), remesh, layer_corr] = lineintersect (GEOMETRY(1,Geo_id==6), GEOMETRY(2,Geo_id==6), GEOMETRY(1,Geo_id==9), GEOMETRY(2,Geo_id==9), remesh, shift);
        Intersect_ID = zeros(1,size(GEOMETRY,2))==1;
    case 'dis_layers'
        [GEOMETRY,Geo_id,Sgap,Intersect_ID,remesh] = ...
            dis_layers(GEOMETRY,Geo_id,Intersect_ID,shift,remesh);
end
toc(inter_t)
    
    if(remesh)
        disp('Remeshing..')
        
        % REMESH
        [remesh,GCOORD,ELEM2NODE,Phases,GEOMETRY,Geo_id,Corner_id,Cornin_id, ...
            Point_id,cont_points,Elsizes,nip,x_min,x_max,ext_rate,ext_erate, ...
            bc_v,temp_bc,temp_bc_depth,mode,prec, ...
            Temp,F_xx,F_xy,F_yx,F_yy,TAU_xx_old,TAU_xy_old, ...
            TAU_yy_old,I2,E2all,Mu_all,ndof,nnod,nel,layer_corr,Bc_ind,Bc_val ...
            Bc_ind_fs,Bc_val_fs,Bct_ind,Bct_val,Geo_mark,Intersect_ID,Sgap, ...
            RHEOL.var,SS,Dpl] = ...
    remesh_all3l_dpl(remesh,GCOORD,ELEM2NODE,GEOMETRY,Phases,Geo_id, ...
    Corner_id,Cornin_id,cont_points,Elsizes,nip,x_min,x_max,Ylim,ext_rate, ...
    ext_erate,bc_v,temp_bc,temp_bc_depth,mode,prec, ...
    shift,Temp,F_xx,F_xy,F_yx,F_yy, ...
    TAU_xx_old,TAU_xy_old, TAU_yy_old, plot_s,Intersect_ID,Sgap, ...
    intersect_s,triangle_path,meshname,RHEOL.var,E2all,Mu_all,SS,Dpl);
    
    end
    
    % UPDATE VELOCITIES
    [Vel,PRES_IP,TAU_xx,TAU_yy,TAU_xy,TAU_xx_old,TAU_yy_old,TAU_xy_old,STRAIN_xx,STRAIN_yy, ...
    STRAIN_xy,E2all,Mu_all,Mu_dis_all,Mu_dif_all,Mu_b_all,Dx_e,F_xx,F_xy,F_yx,F_yy,I2,GIP_x_all,GIP_y_all, ...
    THETA_all,W_xy,RHO,nw_it] = mechanical2d_PFL(ELEM2NODE,Phases,GCOORD,Temp,E2all, ...
    Mu_all,RHEOL,Phi,Cohesion,R,Grain,Burger,Shearm,Rho,G,Bc_ind, ...
    Bc_val,nip,reorder,ext_erate,top_surface,Corner_id,Point_id,dt, ...
    alpha,beta,Bc_ind_fs,Bc_val_fs,SS,F_xx,F_xy,F_yx,F_yy,I2,THETA_all, ...
    TAU_xx_old,TAU_yy_old,TAU_xy_old,elasticity_s,Ce,Dpl,Dplf,nw_it, ...
    sea_level,rho_w);

    % SHEAR HEATING
    Hs = shear_heat(TAU_xx,TAU_yy,TAU_xy,TAU_xx_old,TAU_yy_old, ...
        TAU_xy_old,STRAIN_xx,STRAIN_yy,STRAIN_xy,Phases,Shearm,dt);
    Hs =  sh_s*Hs;

    % TEMPERATURE SOLVER / INCLUDES MELT ITERATIONS / SERPENTINIZATION
    Temp = thermal2d_rhot(ELEM2NODE, GCOORD, Phases, K, Rho, Cp, Hp, Hs, Temp, Bct_ind, Bct_val, dt, nip, reorder,Ce);
    
    %INTEGRATE MELT PRODUCTIVITY
  
    % RESHAPE VELOCITIES
    DISPL       = reshape(Vel,[ndof,nnod]);
        
    %UPDATE OF COORDINATES
    GCOORD      = GCOORD + DISPL*dt;
        
    %MAKE STRAIGHT EDGES
    GCOORD(:,ELEM2NODE([6 4 5],:)) = 0.5*(GCOORD(:,ELEM2NODE([1 2 3],:)) + GCOORD(:,ELEM2NODE([2 3 1],:)));
    GCOORD(:,ELEM2NODE(7,:))       = 1/3*(GCOORD(:,ELEM2NODE(1,:))+GCOORD(:,ELEM2NODE(2,:))+GCOORD(:,ELEM2NODE(3,:)));

    %UPDATE BOUNDARY CONDITIONS
    [Bc_ind, Bc_val, Point_id, Bc_ind_fs, Bc_val_fs, ext_erate] = ...
        set_bcs_flow3l(GCOORD, Corner_id, Cornin_id, Point_id, ext_erate, ext_rate, bc_v);
    [Bct_ind, Bct_val] = set_bcs_temp(GCOORD, Corner_id, Cornin_id, Point_id, temp_bc, temp_bc_depth);
    
    %CHECK MESH INTEGRETY
    [qn,amin,amax] = checkmesh(GCOORD,ELEM2NODE);
    if(qn<0.2 || amin<7 || amax>170)
        remesh = 1;
    end
    
    %UPDATE GEOMETRY MARKERS
    GEOMETRY = GCOORD(:,Geo_mark);
    
    etime = etime  + dt;
    tint = toc(t_it);
    
    %======================================================================
    % POSTPROCESSING
    %======================================================================
    clf
    figure(1)
    plot_t
    hold on
    plot([125 125],[0 temp_bc],'--k')
    plot([0 -temp_bc_depth/1000],[1350 1350],'--k')
    hold off
    figure(2)
    plot_rho
    hold on
    plot_flow
    pause(0.1)
    %plot_eri
    istep = istep+1;
    ntime = ntime+dt;

    fprintf(1, [num2str(toc,'%8.6f'),'\n']);
    fprintf(1, ['\n']);
    
    ttotal = ttotal + tint;
    
    clc;
end