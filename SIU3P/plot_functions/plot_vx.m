% X-VELOCITY FIELD PLOT
% VELOCITY QUIVER PLOT WITH DENSITIES
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
%
% TypeV         Type of velocity                                'absolute'
%               Plots abosulte velocity 'absolute'
%               Plots velocities        'relative'
%                   relative to vel at
%                   the surface

%--------------------------------------------------------------------------
% Function written by Lars R??pke. Edited by Miguel Andres-Martinez, 
% 17-02-2015. Email: mandresmartinez87@gmail.com
%--------------------------------------------------------------------------

% Check type of velocity
if ~exist('TypeV','var')
    TypeV = 'absolute';
end

%==========================================================================
% MAPPING FROM INTEGRATION POINTS TO NODES - IP CANNOT BE PLOTTED
%==========================================================================
EL2N      = zeros(nel,3); % new connectivity matrix
GCOORD_N  = zeros(2,nel*3);
for i=1:nel
    is         = (i-1)*3+1; ie = (i-1)*3+3;
    GCOORD_N(:,is:ie) = GCOORD(:,ELEM2NODE([1 2 3],i));
    EL2N(i,:) = is:ie;
end

%==========================================================================
% PLOT
%==========================================================================
% Plot x-velocity field
% ---------------------
Vel_x = DISPL(1,ELEM2NODE(1:3,:));

% Switch for relative velocities
switch TypeV
    case 'absolute'
    case 'relative'
        [Topography,Topo2nodes] = find_topo(GCOORD,ELEM2NODE,Point_id);
        Vxtop = DISPL(1,Topo2nodes);
        VXTOP = interp1(Topography(1,:),Vxtop,GCOORD(1,ELEM2NODE(1:3,:)));
        Vel_x = Vel_x-VXTOP;
end

title (['Vx [Km/Myr] (',num2str(istep*dt/ma),' Myr)'])
xlabel('Distance Kkm]')
ylabel('Depth [Km]')
patch('faces',EL2N,'vertices',GCOORD_N'/1000,'facevertexcdata', ...
    Vel_x(:)*ma/km,'FaceColor','flat')
shading interp
colormap(jet)
hold on
hc = colorbar('location','southoutside');
xlabel(hc,'Velocity [Km/Myr]','FontSize',10)

% Plot box and interfaces
% -----------------------
plot_box

hold off
drawnow