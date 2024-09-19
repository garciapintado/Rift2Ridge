% DIKES
%   Optional input variables:
%
%               USE                     OPTIONS                 DEFAULTS
%
% color_mag     Color of the box        #interfaces x 3 vector  Red
%               and interfaces          with values [0 1]

%--------------------------------------------------------------------------
% Function written by Elena Ros
%--------------------------------------------------------------------------

% Define color for the plot in case color_mag doesn't exist
if ~exist('color_mag','var')
    color_mag = [1 0 0];
end

hold on;
%if exist('TP_xmelt')~=0
if tp_melt==1
    if Crust_thickness(istep)~=0
        %FOR DIKE_SILL_______________________________
        k = TRACKP_melt(4,:)<istep; %the same as max(find(ismember(TRACKP_melt(3,:),istep-1))); and more efficient
        %all the columns where there is answer 1 are from the previous time step, so I have to plot till there, till the maximum.
        plot(TRACKP_melt(1,k)/1000,TRACKP_melt(2,k)/1000,'.','Color',[0.5 0 0],'MarkerSize',10)
        %____________________________________________
        drawnow
    end
end
hold on;
if(Crust_thickness(istep)~=0)
    hold on
    %FOR DIKE____________________________________
    %                fill(xdike/1000,ydike/1000,'m')
    %                hold on;
    %                plot(xdike/1000,ydike/1000,'.-m')
    %FOR DIKE_SILL_______________________________
    fill(x_ign_body/1000,y_ign_body/1000,color_mag)
    plot(x_ign_body/1000,y_ign_body/1000,'-','Color',color_mag)
    %____________________________________________
end