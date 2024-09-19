% MESH PLOT
%   Optional input variables:
%
%               USE                     OPTIONS                 DEFAULTS
%
% color_int     Color of the box        #interfaces x 3 vector  Black
%               and interfaces          with values [0 1]
%
% line_width    Width of the box        #interfaces x 1 vector  1
%               and interfaces          with width values
%
% color_el      Color of the edges      1 x 3 vector with       Grey
%               of the elements         values [0 1]
%
% el_width      Width of the element    1 x 1 vector with       1
%               edges                   width values

%--------------------------------------------------------------------------
% Function written by Lars R??pke. Edited by Miguel Andres-Martinez, 
% 17-02-2015. Email: mandresmartinez87@gmail.com
%--------------------------------------------------------------------------

% Define color for the plot in case color_el doesn't exist
if ~exist('color_el','var')
    color_el = [0.5 0.5 0.5];
end

% Define the linewidth for the plot in case el_width doesn't exist
if ~exist('el_width','var')
    el_width = 1;
end

EL2N      = zeros(nel,3); % new connectivity matrix
GCOORD_N  = zeros(2,nel*3);

for i=1:nel
    is = (i-1)*3+1;
    ie = is + 2;
    GCOORD_N(:,is:ie) = GCOORD(:,ELEM2NODE([1 2 3],i));
    EL2N(i,:) = is:ie;
end

% Plot mesh
% ---------
title(['Mesh (',num2str(istep*dt/ma),' Myr)'])
xlabel('Distance [Km]')
ylabel('Depth [Km]')
patch('faces',EL2N,'vertices',(GCOORD_N')/1000,'FaceColor','none','EdgeColor',color_el)
hold on

% Plot box and interfaces
% -----------------------
line_width = 2*ones(max(Phases),1);
plot_box

drawnow
hold off