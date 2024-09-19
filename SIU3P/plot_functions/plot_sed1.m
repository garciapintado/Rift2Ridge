function plot_sed(directory_p,name_p,t,t_horizon)

% PLOT_SED(DIRECTORY_P,NAME_P,T,T_SED,AXIS_P) plots the topography of the 
% basement and the setiments of a time T, with layers of sediments that 
% occur at T_HORIZON times (should be an horizontal vector, in brackets). 
% DIRECTORY_P refers to the directory where the data is saved and NAME to
% the name of the files. Write T and T_HORIZON as an empty vector [] if you
% wish the function to display the time interval between each save, to then 
% input manually T and T_HORIZON.

%--------------------------------------------------------------------------
% Function written by Miguel Andres-Martinez, PhD student at Royal Holloway
% University of London, 13-02-2014. Email: mandresmartinez87@gmail.com
%--------------------------------------------------------------------------

% Find all saved files in the directory
files = dir(directory_p);
% Loads the last saved file
load([directory_p,files(3).name],'time_int','n_saves','dt','ma','istep')

% Calculates save interval
save_int = time_int/n_saves/ma;
% Calculates the time of the last save
last_s = istep*dt/ma;

% Display times
disp(['The time step is ',num2str(dt/ma),' Ma'])
disp(['The time interval between each save is ',num2str(save_int),' Ma'])
disp(['The model was run until ',num2str(last_s),' Ma'])

% Checks if either T or T_HORIZON are empty and ask for input
if isempty(t)
    t = input(['At what time [Ma] do you want to plot the sediments? ']);
end

if isempty(t_horizon)
    t_horizon = input(['At what times [Ma] do you want to plot the layers?',...
        '\n','(note this should be an horizontal vector)','\n']);
end

% Calculates current step
step_p = t/(dt/ma);

% Loads the step for which we want the plot
load([directory_p,name_p,num2str(step_p),'.mat'],'GCOORD','Point_id', ...
    'Corner_id','ISOCHRONS')

%figure(36)
hold on

% Gets topography
X_basement_p = [GCOORD(1,Point_id==max(Point_id)-1) GCOORD(1,Corner_id([3 4]))];
Y_basement_p = [GCOORD(2,Point_id==max(Point_id)-1) GCOORD(2,Corner_id([3 4]))];

% Order the topography
[X_basement_pr,ind_plot] = sort(X_basement_p);
Y_basement_pr = Y_basement_p(ind_plot);

% Sorts ISOCHRONS
ISOCHRONS_sort = ISOCHRONS(:,ind_plot,:);

% Calculate current topographies of the isochrons
ISO_TOPO = (ISOCHRONS_sort(1,:,:)+repmat(Y_basement_pr,[1,1,size(ISOCHRONS_sort,3)]));

% Calculates max sediment topography
max_sed_topo = max(max(max(ISO_TOPO)));

% Loop through the different time-horizons
for jp = 1:length(t_horizon)
    
    % Calculates the step of the horizon
    step_horizon = t_horizon(jp)/(dt/ma);
    
    % Load t_horizon
    load([directory_p,name_p,num2str(step_horizon),'.mat'],'GCOORD', ...
        'Point_id','Corner_id','Sediments','istep','dt')
    
    % Store x location of sediments
    Base_sed = ...
        [GCOORD(1,Point_id==max(Point_id)-1) GCOORD(1,Corner_id([3 4]))];
    % Reorder x location of sediments
    [Base_sed,ind_sed_p] = sort(Base_sed);
    % Reorder thickness of sediments
    Sed_thickness = Sediments(ind_sed_p);
    
    % Interpolate the sediment thickness on the topography of the current
    % time step
    Sed_thickness = interp1(Base_sed,Sed_thickness,X_basement_pr);
    
    % Stores the current topography of the horizon
    Sediment_p(jp,:) = Sed_thickness + Y_basement_pr;
    
    % Removes the eroded sediments
    Xp = [X_basement_pr(1) X_basement_pr X_basement_pr(end)];
    Yp = [max(Sediment_p(jp,:)) Sediment_p(jp,:) max(Sediment_p(jp,:))];
    patch(Xp/1000,Yp/1000,jp*ones(size(Xp)),'r','EdgeColor','none')
    
    % Plot horizon
    plot3(X_basement_pr/1000,Sediment_p(jp,:)/1000,jp*ones(size(X_basement_pr)),'k')
end

% Plots topography
plot3(X_basement_pr/1000,Y_basement_pr/1000,(jp+1)*ones(size(X_basement_pr)),'k','LineWidth',2)
title(['Sediments ',num2str(t),' Ma'])
xlabel('Distance [Km]')
ylabel('Depth [Km]')