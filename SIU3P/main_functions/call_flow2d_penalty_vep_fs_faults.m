function [Vel,Pressure,PRES_IP,TAU_xx,TAU_yy,TAU_xy,STRAIN_xx,STRAIN_yy, ...
    STRAIN_xy,Gamma,Yield_T2,YC,E2all,Mu_all,Mu_dis_all,Mu_dif_all,Mu_b_all,Dx_e,F_xx,F_xy,F_yx,F_yy,I2,GIP_x_all,GIP_y_all, ...
    THETA_all,W_xy,nw_it,Residual] = ...
	 mechanical2d_m(ELEM2NODE, Phases, GCOORD, Temp, E2all, ...
			Mu_all, RHEOL, Phi0, Cohesion0, R, Grain, Burger, Shearm, RHO, G, ...
			Bc_ind, Bc_val, nip, reorder, ext_erate, top_surface, Corner_id, Point_id, ...
			dt, alpha, beta, Bc_ind_fs, Bc_val_fs, bc_t, SS,F_xx, F_xy, F_yx, F_yy,I2, THETA_all, ...
                        TAU_xx_old, TAU_yy_old, TAU_xy_old, elasticity_s, nw_it, ...
                        F_ext, rho_ext, Load_el, PLASTICITY, SCALE, PLOT, SOLVER)

  % +++ purpose +++
  % wrapper for kinedyn flow2d_penalty_vep_fs_faults_j to be callable from rift2ridge2D  

  %  1. ELEM2NODE
  %  2. Phases
  %  3. GCOORD
  %  4. Temp
  %  5. E2all
  %  6. Mu_all
  %  7. RHEOL      .Adif, .Ndif, .Qdif, .Vdif, .var
  %  8. Phi0
  %  9. Cohesion0
  % 10. R
  % 11. Grain
  % 12. Burger
  % 13. Shearm
  % 14. RHO
  % 15. G
  % 16. Bc_ind
  % 17. Bc_val
  % 18. nip
  % 19. reorder
  % 20. ext_erate
  % 21. top_surface
  % 22. Corner_id
  % 23. Point_id
  % 24. dt
  % 25. alpha
  % 26. beta
  % 27. Bc_ind_fs
  % 28. Bc_val_fs
  % 29. bc_t
  % 30. SS
  % 31. F_xx
  % 32. F_xy
  % 33. F_yx
  % 34. F_yy
  % 35. I2
  % 36. THETA_all
  % 37. TAU_xx_old
  % 38. TAU_yy_old
  % 39. TAU_xy_old
  % 40. elasticity_s
  % 41. nw_it
  % 42. F_ext
  % 43. rho_ext
  % 44. Load_el
  % 45. PLASTICITY
  % 46. SCALE
  % 47. PLOT
  % 48. SOLVER
  
  % Vel [2*nnod,1]; (v_1_x v1_y ... v_nnod_x v_nnod_y) -> VAR.Ux, VAR.Uz
  % Pressure [
  
  
  VAR = [];
  VAR.Ux = Vel(1:2:end-1);
  VAR.Uz = Vel(2:2:end);

  VAR.P = 
  
  
end  

