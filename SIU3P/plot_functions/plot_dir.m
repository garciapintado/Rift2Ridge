% Plots all time steps in a data directory

%==========================================================================
% INPUT
%==========================================================================
% Directory
dir_p = '/data/miguel';
% Figure directory
dir_f = '/data/miguel/figures';

% Step
s = 100;

% Axis
axis_p = [-200 200 -120 5];

% Exceptions
Exc = {'Test_figures','Mesh','figures'};

% Color axis
caxis_er = [];
caxis_mu = [];

% Variables to load
var_p = {'GCOORD','ELEM2NODE','istep','dt','ma','km','E2all','Mu_all', ...
    'Temp','nel','Point_id','Corner_id','GEOMETRY','Geo_id'};

addpath('/home/miguel/Dropbox/Codes/plot_functions')
% Creates dir_f
mkdir(dir_f)
% List of plotted files
list = {};

while 1
    % Find subdirectories in dir_l
    dir_l = dir(dir_p);
    % Loop through directories
    for n_p = 1:length(dir_l)-2
        if ~ismember(dir_l(n_p+2).name,Exc)
            % Create subdirectories in dir_f
            mkdir([dir_f,'/',dir_l(n_p+2).name])
            
            % Find models inside subdirectories
            dir_m = dir([dir_p,'/',dir_l(n_p+2).name]);
            
            % Loop through models
            for m_p = 1:length(dir_m)-2
                % Create model directories
                mkdir([dir_f,'/',dir_l(n_p+2).name,'/',dir_m(m_p+2).name])
                % Find steps
                model = [dir_p,'/',dir_l(n_p+2).name,'/',dir_m(m_p+2).name];
                dir_s = dir(model);
                dir_ss = dir_s;
                dir_ss(:) = [];
                dir_fstep = dir([dir_f,'/',dir_l(n_p+2).name,'/',dir_m(m_p+2).name]);
                c = 1;
                for p_p = 1:length(dir_s)-2
                    if ismember(str2double(dir_s(p_p+2).name(2:end-4)),1:s:1e6)
                        dir_ss(c) = dir_s(p_p+2);
                        c = c+1;
                    end
                end
                % Step loop
                for o_p = 1:length(dir_ss)
                    % Build the name of the timestep file so that the model
                    % can check if the file already exist so it saves time
                    % plotting
                    dot_p = strfind(dir_ss(o_p).name,'.');
                    name_p = dir_ss(o_p).name;
                    geom_name = ['g',name_p(1:dot_p-1),'.fig'];
                    % Check if figures are create already
                    if ~ismember(geom_name,{dir_fstep.name})
                        try
                            clf
                            close all
                            % Load variables to plot
                            load([model,'/',dir_ss(o_p).name],var_p{:})
                            
                            % Plot strain rate
                            h = figure('visible','off');
                            plot_eri
                            axis(axis_p)
                            saveas(h,[dir_f,'/',dir_l(n_p+2).name,'/', ...
                                dir_m(m_p+2).name,'/er',dir_ss(o_p).name(1:end-4), ...
                                '.jpg'],'jpg')
                            clf
                            close all
                            
                            % Plot viscosity
                            h = figure('visible','off');
                            plot_mu
                            axis(axis_p)
                            saveas(h,[dir_f,'/',dir_l(n_p+2).name,'/', ...
                                dir_m(m_p+2).name,'/mu',dir_ss(o_p).name(1:end-4), ...
                                '.jpg'],'jpg')
                            clf
                            close all
                            
                            % Plot temperature
                            h = figure('visible','off');
                            plot_t
                            axis(axis_p)
                            saveas(h,[dir_f,'/',dir_l(n_p+2).name,'/', ...
                                dir_m(m_p+2).name,'/t',dir_ss(o_p).name(1:end-4), ...
                                '.jpg'],'jpg')
                            clf
                            close all
                            
                            % Plot geometry
                            h = figure('visible','off');
                            plot_geom
                            axis(axis_p)
                            saveas(h,[dir_f,'/',dir_l(n_p+2).name,'/', ...
                                dir_m(m_p+2).name,'/g',dir_ss(o_p).name(1:end-4), ...
                                '.fig'],'fig')
                            clf
                            close all
                            
                            list{length(list)+1} = [model,'/',dir_ss(o_p).name];
                        end
                        clear(var_p{:})
                    end
                end
            end
        end
    end
end