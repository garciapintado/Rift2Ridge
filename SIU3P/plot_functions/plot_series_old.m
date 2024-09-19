function plot_series ...
    (directory_p,name,type_p,dir_o,format_p,save_stepC,SETp)

% PLOT_SERIES(DIRECTORY_P,NAME,TYPE_P,DIR_O,FORMAT_P,SAVE_STEPC,SETP) plots 
% and saves the series of graphs along time-steps, where DIRECTORY_P is the 
% directory where the different time-steps are, NAME is the common name of
% the saved files without number and format, TYPE_P is the type of plot 
% ('vx','vy','visco','eri','eh'), DIR_O is the output directory, FORMAT_P 
% is the format to save ('eps','pdf','movie','png','fig'), SAVE_STEPC is a 
% vector which indicates the time steps to plot that should be 0 in case 
% every time step in the directory needs to be plotted, and SETP is a
% structure containing the rest of the plot settings. You can edit and run
% PLOT_OPT script before using PLOT_SERIES to create the SETP setting
% structure. For more details on the aditional details about the settings
% in SETP read PLOT_OPT script.
%
% Example:
% plot_series('Saved/test_SS01_elast_smallt','Saved130515_SS_elast', ...
%     'eri','/home/Miguel/Files/Test_figures/','png',0,SETp)

%--------------------------------------------------------------------------
% Function written by Miguel Andres-Martinez, PhD student at Royal Holloway
% University of London, 16-09-2012. Email: mandresmartinez87@gmail.com
%--------------------------------------------------------------------------

% 31/10/2014 MA
    % Added constant color bars
% 28/01/2015 MA
    % Added invisible figures option
% 18/02/2015 MA
    % Highly improved by removing unnecessary parts, calling scripts,
    % adding options and moving all the settings to variable SETp

%==========================================================================
% LOAD DATA
%==========================================================================
files = dir(directory_p);
load([directory_p,'/',files(4).name])
mkdir(dir_o);
TRACKP = {};

% Set figure visible or invisible
switch SETp.visible
    case 'visible'
        h = figure('visible','on');
    case 'invisible'
        h = figure('visible','off');
end

% Load plotting settings
log_data = SETp.log_data;
color_int = SETp.color_int;
line_width = SETp.line_width;
Tisotherm = SETp.isot;
MUisomu = SETp.MUisomu;

% Loop for finding maximum and minimum values for then plotting constant
% colormaps

if length(SETp.ct_colors)==1 && SETp.ct_colors==1
    disp('Loading maximum and minimum variables for color bars...')
    
    % Check type of time step
    if sum(save_stepC)==0
        save_s = save_step;
    else
        save_s = save_stepC;
    end
    
    % Declare maximum and minimums of the values to be able to plot at
    % constant colormaps
    Rangec = zeros(2,length(save_s));
    neg_msg = 0;
    switch type_p
        case 'eri'
            load_var = 'E2all';
        case 'mu'
            load_var = 'Mu_all';
        case 'T'
            load_var = 'Temp';
        case {'vx','vy'}
            load_var = 'DISPL';
        case {'eh'}
            load_var = 'I2';
    end
    
    for ii = 1:length(save_s)
        L = ...
            load([directory_p,'/',name,num2str(save_s(ii)),'.mat'], ...
            load_var);
        if L.(load_var)<0
            neg_msg = 1;
        end
        switch type_p
            case {'eri','mu'}
                Rangec(1,ii) = min(min(abs(L.(load_var))));
                Rangec(2,ii) = max(max(abs(L.(load_var))));
            case {'vx'}
                load([directory_p,'/',name,num2str(save_s(ii)),'.mat'], ...
                    'ELEM2NODE');
                Rangec(1,ii) = min(DISPL(1,ELEM2NODE(1:3,:))*ma/km);
                Rangec(2,ii) = max(DISPL(1,ELEM2NODE(1:3,:))*ma/km);
            case {'vy'}
                load([directory_p,'/',name,num2str(save_s(ii)),'.mat'], ...
                    'ELEM2NODE');
                Rangec(1,ii) = min(DISPL(2,ELEM2NODE(1:3,:))*ma/km);
                Rangec(2,ii) = max(DISPL(2,ELEM2NODE(1:3,:))*ma/km);
            otherwise
                Rangec(1,ii) = min(min(L.(load_var)));
                Rangec(2,ii) = max(max(L.(load_var)));
        end
    end
    RangeC = [min(Rangec(1,:)) max(Rangec(2,:))];
    if log_data && (strcmp(type_p,'mu') || strcmp(type_p,'eri'))
        RangeC = log10(RangeC);
    end
    if neg_msg && (strcmp(type_p,'mu') || strcmp(type_p,'eri'))
        disp('Becarefull!!!! The variable you are plotting have negatives')
    end
end

clc
disp(['Plotting ',directory_p,' ....... ','0 %'])

%==========================================================================
% PLOT LOOP
%==========================================================================
for ii = 1:length(files)-3
    % Load step
    if sum(save_stepC)==0
        load([directory_p,'/',name,num2str(save_step(ii)),'.mat'])
    elseif save_stepC(ii)<=max(save_step)
        load([directory_p,'/',name,num2str(save_stepC(ii)),'.mat'])
    elseif save_stepC(ii)>max(save_step)
        break
    end
    
    clf
    
    % Switch for the different type of plots
    switch type_p
        %==================================================================
        % VERTICAL VELOCITIES
        %==================================================================
        case 'vx'
            plot_vx
            VARp = Vel_x(:)*ma/km;
        %==================================================================
        % VERTICAL VELOCITIES
        %==================================================================
        case 'vy'
            plot_vy
            VARp = Vel_y(:)*ma/km;
        %==================================================================
        % TEMPERATURE (TODO: sediments temperatures in another algorithm)
        %==================================================================
        case 'T'
            plot_gt = 0;
            plot_t
            VARp = T_n;
            
            % Plot sediment temperatures
            if exist('Sediments','var')
                if sum(Sediments) ~= 0
                    plot_t_sed
                    [X_sort,sort_p] = sort(GCOORD(1,Point_id==max(Point_id)-1));
                    Yp = GCOORD(2,Point_id==max(Point_id)-1);
                    Y_sort = Yp(sort_p);
                    plot(X_sort/km,Y_sort/km,'Color',[0.5 0.5 0.5],'LineWidth',1)
                end
            end
        %==================================================================
        % VISCOSITY
        %==================================================================
        case 'mu'
            plot_mu
            VARp = Mu_n;
        %==================================================================
        % 2ND STRAIN RATE INVARIANT
        %==================================================================
        case 'eri'
            plot_eri
            VARp = E2_n;
        %==================================================================
        % HISTORIC 2ND STRAIN INVARIANT
        %==================================================================
        case 'eh'
            plot_eh
            VARp = I2_n;
    end
    
    %======================================================================
    % PLOT WEAK SEED (TODO: remove sum(track_remesh))
    %======================================================================
    if SETp.ws && istep*dt<WS.time
        % Calculate points of the weak seed
        WSp(1,:) = WS.size*cosd([0:360])+WS.coord(1);
        WSp(2,:) = WS.size*sind([0:360])+WS.coord(2);
        hold on
        plot(WSp(1,:)/1000,WSp(2,:)/1000,'k')
    end
    
    %======================================================================
    % SEDIMENTS (TODO)
    %======================================================================
    if SETp.basement
        color_base = SETp.color_base;
        line_widthb = SETp.line_widthb;
         plot_basement
    end
    
    %======================================================================
    % ISOTHERM
    %======================================================================
    if ~isempty('Tisotherm')
        cisot = SETp.cisot;
        lwisot = SETp.lwisot;
        isotherm
    end
    
    %======================================================================
    % ISOVISCOSITY LINES
    %======================================================================
    if ~isempty('MUisomu')
        cisom = SETp.cisom;
        lwisom = SETp.lwisom;
        isomu
    end
    
    %======================================================================
    % TRACK OBJECTS
    %======================================================================
    if SETp.tp && ~isempty(TRACKP)
        hold on
        if ~exist('E2N_TP','var')
            E2N_TP = 0;
        end
        plot_defgrid(TRACKP{end},E2N_TP,tp_x,tp_y,km,SETp.dg_option,'k')
        hold on
    end
    
    %======================================================================
    % PLOT VELOCITY VECTORS
    %======================================================================
    if any(strcmp('vel_indx',fieldnames(SETp)))
        vel_indx = SETp.vel_indx;
        if any(strcmp('scale_v',fieldnames(SETp)))
            scale_v = SETp.scale_v;
        end
        plot_vel
    end
        
    %======================================================================
    % MARGINS AND GENERAL DISPLAY OPTIONS
    %======================================================================
    % Color map
    % ---------
    % If minimum and maximum values along the time steps should be used
    if length(SETp.ct_colors)==1 && SETp.ct_colors==1
        caxis(RangeC)
    % If minimum and maximum values in the current time step should be used
    elseif length(SETp.ct_colors)==1 && SETp.ct_colors==0
        switch type_p
            % For variables that logarithmic option is available
            case {'eri','mu'}
                if isempty(log_data)
                    caxis([min(min(log10(VARp))) max(max(log10(VARp)))])
                elseif log_data
                    caxis([min(min(log10(VARp))) max(max(log10(VARp)))])
                else
                    caxis([min(min(VARp)) max(max(VARp))])
                end
            % For variables that logarithmic option is not available
            otherwise
                caxis([min(min(VARp)) max(max(VARp))])
        end
    % If a fixed range of values should be used
    else
        caxis(SETp.ct_colors)
    end
 
    % Axis definition
    axis(SETp.axis)
    % Aspect ratio
    daspect(SETp.aspect)
    
    drawnow
    hold off
    
    %======================================================================
    % SAVE
    %======================================================================
    switch format_p
        case 'mov'
            if ii==1
                writerObj = VideoWriter([dir_o,'/',type_p,'.avi'], ...
                    'Motion JPEG AVI');
                set(writerObj,'FrameRate',10,'Quality',100);
                open(writerObj);
            end
            movie_(ii) = getframe;
            writeVideo(writerObj,movie_(ii));
        case 'pdf'
            orient tall
            print('-dpdf', '-loose', '-painters', [dir_o,'/',type_p,num2str(ii), '.pdf']);
        case 'eps'
            orient tall
            print('-depsc2', '-loose', '-painters', sprintf('-r%d', 600), strcat([dir_o,'/',type_p,num2str(ii)], '.eps'));
        case 'png'
            saveas(h,[dir_o,'/',type_p,num2str(istep*dt/ma,'%2.2f'),'ky.png'],'png')
        case 'pngeps'
            % Remove axis
            h = gca;
            axis_lim = axis(h);
            set(h,'visible','off')
            set(hc,'visible','off')
            
            % Remove lines
            ax = get(h,'children');
            for n = 1:length(ax)
                if strcmp(ax(n).Type,'line') || strcmp(ax(n).Type,'contour')
                    delete(ax(n))
                end
            end
            
            % Print patch in png
            print('-dpng', '-loose', '-painters', sprintf('-r%d', 600), strcat([dir_o,'/',type_p,num2str(ii)], '.png'));
            
            % Plot again axis and box
            hold on
            plot_box
            if ~isempty('Tisotherm')
                cisot = SETp.cisot;
                lwisot = SETp.lwisot;
                isotherm
                isomu
            end
            hold on
            plot([axis_lim(2) axis_lim(2)],[axis_lim(3) axis_lim(4)], ...
                'color',color_int(1,:),'LineWidth',line_width(1))
            plot([axis_lim(1) axis_lim(1)],[axis_lim(3) axis_lim(4)], ...
                'color',color_int(1,:),'LineWidth',line_width(1))
            set(gca,'visible','on')
            set(hc,'visible','on')
            
            % Remove patch
            for n = 1:length(ax)
                try
                    if strcmp(ax(n).Type,'patch')
                        delete(ax(n))
                    end
                end
            end
            
            % Print eps
            print('-depsc2', '-loose', '-painters', sprintf('-r%d', 600), strcat([dir_o,'/',type_p,num2str(ii)], '.eps'));

        case 'fig'
            saveas(h,[dir_o,'/',type_p,num2str(istep*dt/ma,'%2.2f'),'ky.fig'],'fig')
    end
    
    % Percentage display
    clc
    disp(['Plotting ',directory_p,' ....... ', ...
        num2str(floor(ii/(length(files)-3)*1000)/10),' %'])   
end

switch format_p
    case 'mov'
        close(writerObj);
end