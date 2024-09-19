function plot_tF(z, GCOORD, ELEM2NODE, ...
                 clim, title_lab, timenow, do_shading, plot_gt)
 % +++purpose+++
 % Do a patch plot of the input variable, 
 %  1. z          :: REAL [1,nnod] variable [temperature] at nodes
 %  2. GCOORD     :: REAL [2,nnod]
 %  3. ELEM2NODE  :: INTEGER [nnodel,nel]
 %  6. clim       :: REAL [2], OPTIONAL, colorscale limits. Not given, or []
 %                   for default automated limits
 %  7. title_lab  :: STRING for main lab  
 %  8. timenow    :: INTEGER, OPTIONAL, model step number, for main labeling
 %  9. do_shading :: LOGICAL, OPTIONAL, TRUE for 'shading interp' the patch (visible triangles)
 % 10. plot_gt    :: LOGICAL, OPTIONAL, whether to plot geotherm subplot
 % 
 %
 
 % +++ purpose +++
 % plot mesh using only vertex values in each triangular FEM
 
 ma = 3600*24*365.25*1.0E06; 
 
 %plot_box = true;
 %if nargin < 4
 %  Point_id = 0; 
 %  GEOMETRY = 0.;
 %  plot_box = false;
 %end   
 
 if nargin < 4
    clim = [];
 end

 if nargin < 5
   title_lab = "Temperature  [C]";
 end  
 
 if nargin < 6
   timenow = [];
 end
 
 if nargin < 7
   do_shading = true;    
 end
 
 if nargin < 8
    plot_gt = false;
 end

 nel = size(ELEM2NODE,2); 
 GCOORD_N    = zeros(2, 3*nel); 
 GCOORD_N(:) = GCOORD(:,ELEM2NODE(1:3,:));                                 % [2, 3*nel] 
 EL2N = reshape(1:nel*3,3,nel)';                                           % [nel,3] new connectivity matrix
 z_n  = z(ELEM2NODE(1:3,:));                                               % [3, nel] variable values at the corner nodes

 %==========================================================================
 % PLOT
 %==========================================================================
 % Plot temperature field
 % ----------------------
 if plot_gt
     subplot(2,1,1)
 end
 ax = gca;
 ax.FontSize = 8;
 %ax.Position(1) = ax.Position(1) - 0.05;
 %pos = get(gca, 'Position');
 %pos(1) = 0.055;
 %pos(3) = 0.9;
 %set(gca, 'Position',pos);
 
 if ~isempty(timenow)
   title_lab = strjoin([title_lab," (",num2str(timenow/ma)," Myr)"],"");
 end
 title(title_lab)
 xlabel('Distance [Km]')
 ylabel('Depth [Km]')
 patch('faces',EL2N,'vertices',GCOORD_N'/1000,'facevertexcdata',z_n(:), ...
     'FaceColor','flat')
 if do_shading
   shading interp
 end
 axis tight
 colormap(jet)
 colorbar;
 hold on

  if ~isempty(clim)
    caxis(clim);  
  end
 
  % Plot box and interfaces
  % -----------------------
  %if plot_box
  %  plot_boxF(GCOORD, ELEM2NODE, Point_id, GEOMETRY)
  %  drawnow
  %
  %end
  hold off

  if plot_gt         % plot geotherm
    subplot(2,1,2)
    plot(-GCOORD(2,:)/1000,z,'.')
    xlabel('Depth [Km]')
    ylabel('Temperature [C]')
  end
end
