% Test interpolate sediments

loadsave_file = '/data/disk1/archiveMIGUEL/MARGIN_WIDTH_CRATON/Sh_steepCT200_GR_C20_40_DWO200_SSD1_30_RMU_WST50_moresi_HRLC/_';

load(loadsave_file)

temp_surf = 0;

%% SURFACE PROCESSES    
%%      - Settings                                                  
%--------------------------------------------------------------------------
% If surface processes are needed to be computed surf_proc_s = 1, on the
% contrary surf_proc_s = 0
surf_proc_s = 1;

% Time steps between isochrons to save
iso_sav = 1;

% Function to calculate erosion and sedimentation
surf_proc_func = 'diff_topo_armitage';
%--------------------------------------------------------------------------
%%      - Parameters                                                
%--------------------------------------------------------------------------
% Time step for the surface processes
SP.dt = 0.001*ma;
% Hill-slope diffusion [m^2/yr]
SP.kappa = 100e1;
% Discharge transport coefficient
SP.c = 1e-2;
% Inverse Hack's law
SP.nexp = 1;
% Precipitation rate
SP.alpha_sed = 1;
%--------------------------------------------------------------------------

%%
% plot_gt = 0;
% figure(2)
% plot_t
hold on
% Declare isochrons rows: 1) x-coordinates, 2) y-coordinates, 3) step
% index, 4) erosion/deposited depth
ISOCHRONS = [GCOORD(:,max(Point_id)-1==Point_id); ...
    zeros(1,sum(max(Point_id)-1==Point_id))];
[ISOCHRONS(1,:),indx_isoc] = sort(ISOCHRONS(1,:));
ISOCHRONS(2,:) = ISOCHRONS(2,indx_isoc);
ISOCHRONS = [ISOCHRONS; ISOCHRONS(2,:)];

% Declare surface process function
if surf_proc_s
    surf_proc = str2func(surf_proc_func);
end

%SURFACE PROCESSES
if surf_proc_s
    % Find topography
    [Topography,Topo2nodes] = find_topo(GCOORD,ELEM2NODE,Point_id);
    % Apply surface processes
    [New_topo,SP] = surf_proc(Topography,GCOORD,SP,dt,ma);
    % Update isochrons
    if floor(istep/iso_sav)==(istep/iso_sav)
        ISOCHRONS = [ISOCHRONS [New_topo; ...
            istep*ones(1,size(Topography,2))-1; New_topo(2,:)]];
    end

    % Update the topography and interpolate variables of reshaped elements
    % by surface processes
    [GCOORD,Temp,F_xx,F_xy,F_yx,F_yy,I2,TAU_xx_old,TAU_xy_old, ...
        TAU_yy_old,Mu_all,E2all,remesh] = update_topo_values(GCOORD, ...
        ELEM2NODE,New_topo,Topo2nodes,Temp,temp_surf,F_xx,F_xy,F_yx, ...
        F_yy,I2,TAU_xx_old,TAU_xy_old,TAU_yy_old,Mu_all,E2all, ...
        remesh);
    
    % Plot (uncomment)
    plot(Topography(1,:)/1000,Topography(2,:)/1000,'k')
    plot(New_topo(1,:)/1000,New_topo(2,:)/1000,'r')
end