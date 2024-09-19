% TODO improve this algorithm

%% Input

directory_p = '/mnt/datadrive/archiveMIGUEL/SEDIMENTS/deleteme/';
directory_s = '/home/ross/Desktop/Sediments_John/Isochrones';

% Total time of the run
tlength = 3.41;

% Steps for the plotting
t_plot = 0.01:0.1:3.71;

% Time interval for the isochrones
t_isochrones = 0.1;

%% Plot

for n = 1:length(t_plot)
    clf 
    h=figure(1);
    plot_sed(directory_p,'_',t_plot(n),t_plot(1):t_isochrones:t_plot(n));

saveas(h,[directory_s,'/isoc',num2str(t_plot(n),'%2.2f'),'ky.fig'],'fig')
end