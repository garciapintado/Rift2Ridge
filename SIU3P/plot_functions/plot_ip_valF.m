function plot_ip_valF(z, GCOORD, ELEM2NODE, topo, clim,  title_lab, timenow, show_triangulation, mask_ocean)
  % Direct 2D delaunay plotting of IP matrices at their locations
  %
  % z :: REAL [nel,nip]
  %
  % Javier Garcia-Pintado, MARUM, 2020
  
  if nargin < 5
    clim = [];
  end
  
  if nargin < 6
   title_lab = "";
  end  
 
  if nargin < 7
    timenow = [];
  end
 
  if nargin < 8
    show_triangulation = false;    
  end
  
  if nargin < 9
    mask_ocean = true;
  end
  
  nip = 6;
  nel = size(ELEM2NODE,2);
  
  [GIPx,GIPy] = ip_coord(GCOORD, ELEM2NODE, nel, nip);                         % both in [nel,nip]
  tri = delaunay(GIPx(:),GIPy(:));
  %plot(GIPx(:),GIPy(:),'.');
  ax = gca;
  ax.FontSize = 8;
  xlabel('Distance [Km]')
  ylabel('Depth [Km]')
  patch('Faces',tri,'Vertices',[GIPx(:),GIPy(:)]/1000,'FaceVertexCData',z(:),'FaceColor','interp');
  if ~show_triangulation
    shading interp
  end
  axis tight
  colormap(jet)
  colorbar;
 
  if mask_ocean
      %[topo, ~] = find_topo(GCOORD, ELEM2NODE, Point_id);
      topo = [topo(:,1)+[0; 2000] topo topo(:,end)+[0; 2000]];
      ps = polyshape(topo(1,:)/1000,topo(2,:)/1000, 'Simplify', false);
      hold on; plot(ps,'FaceColor','white','FaceAlpha',1,'EdgeAlpha',0);
  end
  
  if ~isempty(clim)
    caxis(clim);  
  end
  if title_lab ~= ""
      title(title_lab)
  end    
end % function