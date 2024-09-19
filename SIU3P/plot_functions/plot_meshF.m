function plot_meshF(ELEM2NODE, GCOORD, elboo, color_el, lwd, hold_on)
  % plot_meshF(ELEM2NODE, GCOORD, elboo, color_el, hold_on)
  %
  % ELEM2NODE :: INTEGER [nnodel,nel]
  % GCOORD    :: REAL [2,nnod]  
  % elboo     :: LOGICAL or INTEGER element subsetting vector. Empty '[]' to plot all elements 
  % color_el  :: element edge color, default to [0.5 0.5 0.5]
  % lwd       :: line width, default to 1
  % hold on   :: LOGICAL, default to true

  if nargin < 3
     elboo = [];
  end
  
  if nargin < 4
    color_el = [0.5 0.5 0.5];
  end
  if nargin < 5
    lwd = 1.0;
  end
  if nargin < 6
      hold_on = true; 
  end
  nel = length(ELEM2NODE);
  
  if isempty(elboo)
    elboo = repelem(true,nel);
  end

  EL2N    = ELEM2NODE(1:3,elboo)';
  GCOORDV  = GCOORD(:,1:max(max(EL2N)))';
  patch('faces',EL2N, 'vertices',GCOORDV/1000., 'FaceColor','none', 'EdgeColor',color_el,'LineWidth',lwd);
  
  %title(['Mesh (',num2str(istep*dt/ma),' Myr)'])
  xlabel('Distance [Km]')
  ylabel('Depth [Km]')

  if hold_on
      hold on;
  end
  % Plot box and interfaces
  % -----------------------
  %line_width = 2*ones(max(Phases),1);
  %plot_box
  drawnow

end
