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
% Authors: Miguel Andres-Martinez, PhD student at Royal Holloway
% University of London, 17-02-2015. Email: mandresmartinez87@gmail.com
%
% versions:
% Javier GP. 2020-03. Converted into a function 
%--------------------------------------------------------------------------
function plot_phasesF(ELEM2NODE, GCOORD, Phases, phaseids, GEOn)

    if nargin < 4
        phaseids = unique(Phases);
    end
    if nargin < 5
        pha_type = "triangles";
    else
        pha_type = "apatch";
    end
    
    % Plot phases
    % -----------
    switch pha_type
        case 'apatch'
            % Plot box and interfaces
            % -----------------------
            phase_colors = [...
                57.00, 62.00, 40.00; ...                    % RGB sampled with Gimp from plot
                89.74, 87.26, 81.16; ...
                78.00, 69.00, 61.00] / 100.0;
            reglabels = ["lithosperic mantle","lower crust","upper crust"];
            
            GEOint = [1 2 3 4;
                      3 5 6 7;
                      6 8 9 10];
                  
            for np = 1:length(phaseids)
                hold on
                pid = phaseids(np);
                GCOO = [GEOn(GEOint(pid,1)).coo, GEOn(GEOint(pid,2)).coo, ...
                       flip(GEOn(GEOint(pid,3)).coo,2) ,flip(GEOn(GEOint(pid,4)).coo,2)] /1000;
                pgon = polyshape(GCOO(1,:), GCOO(2,:),'Simplify',false);
                plot(pgon, 'FaceColor',phase_colors(pid,:), 'FaceAlpha',1);
            end
        case 'triangles'
            %title(['Phases (',num2str(istep*dt/ma),' Myr)'])
            %xlabel('Distance [km]')
            %ylabel('Depth [km]')
            phaseboo = ismember(Phases,phaseids);
            patch('faces',ELEM2NODE(1:3,:)','vertices',GCOORD'/1000, ...
                  'facevertexcdata',Phases(:),'FaceColor','flat','EdgeColor','none')
            load('cmap_phases')
            colormap(cmap_phases)
            caxis([0.7 3]) 
    end
end % function plot_phasesF
