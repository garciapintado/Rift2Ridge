% PHASES PLOT
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
% pha_type      Type of phases plot     'apatch' to plot one    apatch
%                                       element per layer
%                                       'triangles' to plot
%                                       every triangle

%--------------------------------------------------------------------------
% Function written by Miguel Andres-Martinez, PhD student at Royal Holloway
% University of London, 17-02-2015. Email: mandresmartinez87@gmail.com
%--------------------------------------------------------------------------

% Define pha_type if not existent
if ~exist('pha_type','var')
    pha_type = 'triangles';
end

% Plot phases
% -----------
switch pha_type
    case 'apatch'
        % Plot box and interfaces
        % -----------------------
        plot_box
        
        hold on 
        for np = 1:nlay
            lind = Box(3,:)==np;
            if np>1
                bot_int = [GCOORD(:,Cornin_id(np*2-3)) ...
                    GCOORD(:,ismember(Point_id,np*3-3) & pv) ...
                    GCOORD(:,Cornin_id(np*2-2))];
                BOX = [Box(1:2,lind) bot_int(1:2,end:-1:1)];
            else
                BOX = Box(1:2,lind);
            end
            sb = size(BOX,2);
            patch('faces',1:sb,'vertices',BOX'/1000, ...
                'facevertexcdata',np,'FaceColor','flat','EdgeColor','none')
        end
        
        
    case 'triangles'
        title(['Phases (',num2str(istep*dt/ma),' Myr)'])
        xlabel('Distance [km]')
        ylabel('Depth [km]')
        patch('faces',ELEM2NODE(1:3,:)','vertices',GCOORD'/1000, ...
            'facevertexcdata',Phases(:),'FaceColor','flat','EdgeColor','none')
end

load('cmap_phases')
colormap(cmap_phases)
caxis([0.7 3])
axis tight
hold on

% Plot box and interfaces
% -----------------------
plot_box

drawnow
hold off