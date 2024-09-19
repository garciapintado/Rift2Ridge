% PLOT UNCONFORMITIES AND DIAGRAM
%   Optional input variables:
%
%               USE                     OPTIONS                 DEFAULTS
%
% un_color      Color of the erosive    2x3 vector              Red and
%               and hiatus                                      blue
%               unconformities
%
%
% diag_colormap Color map of the in     colormap name           cmap_diag
%               diagram
%
% min_ero_rate  Minimum erosion rate to scalar in m/yr          0
%               be considered erosion
%               and not a depositional
%               break
%
% min_sed_rate  Minimum sedimentation   scalar in m/yr          1e-8
%               rate to be considered
%               deposition and not a
%               depositional break
%
% un_plot_type  Define the plot types:  'separate'              'separate'
%               in separate figures,    'together'
%               together in a           'unconformities'
%               subfigure or only       'diagram'
%               unconformities or
%               diagram
%
% min_segl      Minimum segment length  scalar in number of     0
%               to plot an              points
%               unconformity segment

%--------------------------------------------------------------------------
% Function written by Miguel Andres-Martinez, postdoc at University of
% Bremen, 05-11-2018. Email: mandresmartinez87@gmail.com
%--------------------------------------------------------------------------

%==========================================================================
% INPUT
%==========================================================================
% Define default color of the unconformities, both erosive and hiatus
if ~exist('un_color','var')
    un_color = [1 0 0; 0 1 0];
end

% Define default color map for the diagram
if ~exist('diag_colormap','var')
    diag_colormap = 'cmap_diag';
end

% Define the default deposition and erosion rates that are considered 
% negligible (no deposition/erosion)
if ~exist('min_ero_rate','var')
    % In m/yr
    min_ero_rate = 0;
end
if ~exist('min_sed_rate','var')
    % In m/yr
    min_sed_rate = 1e-8;
end

% Define default figure type
if ~exist('un_plot_type','var')
    % In m/yr
    un_plot_type = 'separate';
end

% Define default minimum segment length to plot as part of an unconformity
if ~exist('min_segl','var')
    min_segl = 0;
end

% Define unconformity plot
hun = gca;

% TODO fix this
u_sed_res = 1;
remove_initial_ero = 1;

%==========================================================================
% ERODE ISOCHRONS
%==========================================================================
% Remesh isochrons to have the same x-points as the topography
[ISOCHRONSp,~] = remesh_isoc(GCOORD,Point_id,ELEM2NODE,ISOCHRONS, ...
    Basement,tp_isoc);

% Calculate topography
[Topography,~] = find_topo(GCOORD,ELEM2NODE,Point_id);

% Reshape timeline matrices
isoc_i = unique(ISOCHRONSp(3,:));
ISOCHRONSpp = ISOCHRONSp(1:2,ismember(ISOCHRONSp(3,:),isoc_i));
ISOCx = reshape(ISOCHRONSpp(1,:),sum(isoc_i(1)==ISOCHRONSp(3,:)), ...
    length(isoc_i))';
ISOCy = reshape(ISOCHRONSpp(2,:),sum(isoc_i(1)==ISOCHRONSp(3,:)), ...
    length(isoc_i))';
ISOCyero = ISOCy;

% Find isochrons that have several points on the x-axis
isoc_i = unique(ISOCHRONS(3,:));
not_fx = [];
for n = 1:length(isoc_i)
    if any(diff(ISOCHRONS(1,ISOCHRONS(3,:)==isoc_i(n)))<=0)
        not_fx = [not_fx n];
    end
end

if isempty(not_fx)
    %======================================================================
    % Build erosion/sedimentation matrix
    %======================================================================
    % Initialise matrix of topographic change
    TOPOc = zeros(size(ISOCy));
    ISOCyero(1,:) = min([ISOCy; Topography(2,:)]);
    for n = 2:length(isoc_i)
        TOPOc(n,:) = ISOCy(n,:)-ISOCy(n-1,:);
        ISOCyero(n,:) = min([ISOCy(n:end,:); Topography(2,:)]);
    end
    
    % Calculate min change in heigth to be considered erosion and
    % deposition
    min_dh_sed = min_sed_rate*dt/year*iso_sav;
    min_dh_ero = min_ero_rate*dt/year*iso_sav;
    
    % Plot
    %----------------------------------------------------------------------
    % Define matrix for erosion/deposition/hiatus (-1/1/0)
    colorTC = zeros(size(TOPOc));
    colorTC(TOPOc>=min_dh_sed) = 1;
    colorTC(TOPOc<=-min_dh_ero) = -1;
    % Exclude initial erosion phase if required by the user, where initial
    % erosion terms are set to -2
    cut_ero_info = remove_initial_ero/dt*ma;
    cut_indx = find(isoc_i==cut_ero_info)-1;
    colorTCc = colorTC;
    colorTCc(1:cut_indx,:) = -2;
    % Find groups of connected types
    labelTC = label(colorTCc);
    
    %     % Plot through erosion/deposition/hiatus groups (uncomment)
    %     for nnn = 1:max(labelTC(:))
    %         imagesc(ISOCx(1,:)/km,isoc_i*dt/ma,labelTC~=nnn)
    %         if unique(colorTCc(labelTC==nnn)<0)
    %             colormap(gray)
    %         else
    %             colormap(flipud(copper))
    %         end
    %         colorbar
    %         xlabel('Distance [km]')
    %         ylabel('Time [Myr]')
    %         set(gca,'Ydir','normal')
    %         input(num2str(nnn));
    %     end
    
    %======================================================================
    % EROSIVE EVENTS
    %======================================================================
    % Define erosive events
    ero_event = colorTCc==-1;
    % Add a layer of erosive material below for correct indexing
    sed_below_mat = [-ones(1,size(colorTCc,2)); ...
        colorTCc(1:end-1,:)];
    % Find terms of the matrix over sediments
    sed_below = sed_below_mat==1 & ero_event;
    % Add a layer of erosive material above for correct indexing
    sed_above_mat = [colorTCc(2:end,:); -ones(1,size(colorTCc,2))];
    % Find terms of the matrix below sediments
    sed_above = sed_above_mat==1 & ero_event;
    % Find vertical boundaries of the erosive event
    v = ero_event & [ero_event(2:end,:); false(1,size(ero_event,2))];
    vert = (v & [colorTCc(:,2:end) -ones(size(colorTCc,1),1)]==1) | ...
        (v & [-ones(size(colorTCc,1),1) colorTCc(:,1:end-1)]==1);
    erosive = sed_above | vert;
    %erosive = (sed_above & ~sed_below) | vert;
    % Group erosive tops into unconformity groups
    UC_ero = label(erosive);
    % Remove vertical sectors
    UC_ero(vert) = 0;
    UC_ero(~sed_above) = 0;
    %UC_ero(~(sed_above & ~sed_below)) = 0;
    [C,~,ic] = unique(UC_ero);
    group_indx = 0:length(C);
    UC_ero = reshape(group_indx(ic),size(UC_ero));
    
    % Remove small unconformities
    [un_UC_ero,~,iu] = unique(UC_ero);
    counts = histc(UC_ero(:),un_UC_ero);
    size_UC_eros = reshape(counts(iu),size(UC_ero));
    UC_ero(size_UC_eros<min_segl) = 0;
    [C,~,ic] = unique(UC_ero);
    group_indx = 0:length(C);
    UC_ero = reshape(group_indx(ic),size(UC_ero));
    
%     % Plot (uncomment)
%     figure(2); clf
%     time_p = repmat(isoc_i'*dt/ma,1,size(ISOCx,2));
%     plot(ISOCx(ero_event)/km,time_p(ero_event),'.r')
%     xlabel('Distance [km]')
%     ylabel('Time [Myr]')
%     set(gca,'Ydir','normal')
%     hold on
%     plot(ISOCx(sed_below)/km,time_p(sed_below),'og')
%     plot(ISOCx(sed_above)/km,time_p(sed_above),'ob')
%     plot(ISOCx(vert)/km,time_p(vert),'.b')
%     plid = [];
%     uu.a = [];
%     uu.b = [];
%     uu.c = [];
%     for m = 1:max(UC_ero(:))
%         delete(plid)
%         delete(uu.a)
%         delete(uu.b)
%         delete(uu.c)
%         figure(2)
%         plid = plot(ISOCx(UC_ero==m)/km, ...
%             time_p(UC_ero==m),'.k','Markersize',10);
%         figure(1)
%         ero_line = [ISOCx(UC_ero==m)'; ISOCyero(UC_ero==m)'];
%         uu.a = plot(ero_line(1,:)/1000,ero_line(2,:)/1000,'k','LineWidth',5);
%         uu.b = plot(ero_line(1,:)/1000,ero_line(2,:)/1000,'w','LineWidth',2);
%         uu.c = plot(ero_line(1,:)/1000,ero_line(2,:)/1000,'--','Color', ...
%             un_color(1,:),'LineWidth',2);
%         input(num2str(m))
%     end
    
    %======================================================================
    % NON-DEPOSITIONAL EVENTS
    %======================================================================
    % Define non-depositional events
    hiat_event = colorTCc==0;
    % Add a layer of erosive material below for correct indexing
    sed_below_mat = [-ones(1,size(colorTCc,2)); ...
        colorTCc(1:end-1,:)];
    % Find terms of the matrix over sediments
    sed_below = sed_below_mat==1 & hiat_event;
    % Add a layer of erosive material above for correct indexing
    sed_above_mat = [colorTCc(2:end,:); -ones(1,size(colorTCc,2))];
    % Find terms of the matrix below sediments
    sed_above = sed_above_mat==1 & hiat_event;
    % Find vertical boundaries of the erosive event
    v = hiat_event & [hiat_event(2:end,:); false(1,size(hiat_event,2))];
    vert = (v & [colorTCc(:,2:end) -ones(size(colorTCc,1),1)]==1) | ...
        (v & [-ones(size(colorTCc,1),1) colorTCc(:,1:end-1)]==1);
    hiat = sed_above | vert;
    %hiat = (sed_above & ~sed_below) | vert;
    % Group erosive tops into unconformity groups
    UC_hiat = label(hiat);
    % Remove vertical sectors
    UC_hiat(vert) = 0;
    UC_hiat(~sed_above) = 0;
    %UC_hiat(~(sed_above & ~sed_below)) = 0;
    [C,~,ic] = unique(UC_hiat);
    group_indx = 0:length(C);
    UC_hiat = reshape(group_indx(ic),size(UC_hiat));
    
    % Remove small unconformities
    [un_UC_hiat,~,iu] = unique(UC_hiat);
    counts = histc(UC_hiat(:),un_UC_hiat);
    size_UC_hiats = reshape(counts(iu),size(UC_hiat));
    UC_hiat(size_UC_hiats<min_segl) = 0;
    [C,~,ic] = unique(UC_hiat);
    group_indx = 0:length(C);
    UC_hiat = reshape(group_indx(ic),size(UC_hiat));
    
    %======================================================================
    % PLOTTING
    %======================================================================
    % Remove unconformity sectors that are either at the basement or at the
    % surface
    rm_surf = abs(ISOCyero-repmat(ISOCyero(end,:),size(ISOCyero,1),1))< ...
        u_sed_res;
    rm_base = abs(ISOCyero-repmat(ISOCyero(1,:),size(ISOCyero,1),1))< ...
        u_sed_res;
    % Find plotting time
    time_p = repmat(isoc_i'*dt/ma,1,size(ISOCx,2));
    
    % Switches for plot type
    switch un_plot_type
        case 'separate'
            h1 = hun;
            h2 = figure(2); clf;
        case 'together'
            h1 = subplot(211);
            h2 = subplot(212);
        case 'unconformities'
            h1 = hun;
        case 'diagram'
            h2 = hun;
    end
    
    % Plot deposition/erosion vs age diagram
    if any(strcmp(un_plot_type,{'separate','together','diagram'}))
        axes(h2)
        plot(ISOCx(ero_event)/km,time_p(ero_event),'sq', ...
            'Color',[0.5 0.5 0.5],'MarkerFaceColor',[0.5 0.5 0.5],...
            'MarkerSize',10) 
        colormap(diag_colormap)
        colorbar
        xlabel('Distance [km]')
        ylabel('Time [Myr]')
        set(gca,'Ydir','normal')
    end
    
    % Plot unconformities in the selected figure
    if any(strcmp(un_plot_type,{'separate','together','unconformities'}))
        axes(h1)
        hold on
        
        % Plot disconformities (non-depositional)
        hiat_lines = [];
        for n = 1:max(UC_hiat(:))
            hiat_line = [ISOCx(UC_hiat==n)'; ISOCyero(UC_hiat==n)'];
            hiat_line(:,(rm_surf(UC_hiat==n) | rm_base(UC_hiat==n))) = NaN;
            hiat_lines = [hiat_lines NaN(2,1) hiat_line];
        end
        if ~isempty(hiat_lines)
            plot(hiat_lines(1,:)/1000,hiat_lines(2,:)/1000,'.-k','LineWidth',5)
            plot(hiat_lines(1,:)/1000,hiat_lines(2,:)/1000,'w','LineWidth',2)
            plot(hiat_lines(1,:)/1000,hiat_lines(2,:)/1000,'--','Color', ...
                un_color(2,:),'LineWidth',2)
        end
        
        % Plot erosional unconformities
        ero_lines = [];
        for n = 1:max(UC_ero(:))
            ero_line = [ISOCx(UC_ero==n)'; ISOCyero(UC_ero==n)'];
            ero_line(:,(rm_surf(UC_ero==n) | rm_base(UC_ero==n))) = NaN;
            ero_lines = [ero_lines NaN(2,1) ero_line];
        end
        if ~isempty(ero_lines)
            plot(ero_lines(1,:)/1000,ero_lines(2,:)/1000,'k','LineWidth',5)
            plot(ero_lines(1,:)/1000,ero_lines(2,:)/1000,'w','LineWidth',2)
            plot(ero_lines(1,:)/1000,ero_lines(2,:)/1000,'--','Color', ...
                un_color(1,:),'LineWidth',2)
        end
    end    

else
%     %======================================================================
%     % Slower cell-based for time-lines that are parametric functions
%     %======================================================================
%     % Isochrons to plot
%     isoc2p = isoc_i(2:isoc_int:end);
%     % Initialisation of variables
%     ISOCx = {};
%     ISOCy = {};
%     ISOC = {};
%     for n = 1:length(isoc_i)
%         ISOC{n} = ISOCHRONS(:,ISOCHRONS(3,:)==isoc_i(n));
%     end
%     
%     % Loop through eroding isochrons
%     % ------------------------------
%     if ppar==1
%         % Parallel loop
%         parfor n = 1:length(isoc2p)
%             % Load evaluated isochron
%             isocx = ISOCHRONS(1,ISOCHRONS(3,:)==isoc2p(n));
%             isocy = ISOCHRONS(2,ISOCHRONS(3,:)==isoc2p(n));
%             % Find older isochrons
%             isoc_above = find(isoc_i>isoc2p(n));
%             
%             % Run erosion algorithm in serie
%             [ISOCx{n},ISOCy{n}] = ero_ISOC_perp(isocx,isocy,ISOC,isoc_above);
%         end
%     else
%         % Serial loop
%         for n = 1:length(isoc2p)
%             % Load evaluated isochron
%             isocx = ISOCHRONS(1,ISOCHRONS(3,:)==isoc2p(n));
%             isocy = ISOCHRONS(2,ISOCHRONS(3,:)==isoc2p(n));
%             % Find older isochrons
%             isoc_above = find(isoc_i>isoc2p(n));
%             
%             % Run erosion algorithm in serie
%             [ISOCx{n},ISOCy{n}] = ero_ISOC_perp(isocx,isocy,ISOC,isoc_above);
%             %plot(ISOCx{n}/1000,ISOCy{n}/1000); drawnow; hold on
%         end
%     end
%     
%     %======================================================================
%     % Plot
%     %======================================================================
%     isoc_p = isoc_i(2:isoc_int:end)*dt/ma;
%     m = 1:length(ISOCx);
%     hold on
%     
%     switch isoc_type
%         case 'fill'
%             % Plot sediment ages
%             for n = 1:length(isoc_p)-1
%                 patch([ISOCx{m(n)} ISOCx{m(n+1)}(end:-1:1)]/1000, ...
%                     [ISOCy{m(n)} ISOCy{m(n+1)}(end:-1:1)]/1000,isoc_p(n))
%             end
%             patch([ISOCx{m(n+1)} Topography(1,end:-1:1)]/1000, ...
%                 [ISOCy{m(n+1)} Topography(2,end:-1:1)]/1000,istep*dt/ma)
%             colormap('jet')
%             hc = colorbar;
%             
%         case 'isolines'
%             % Plot isochrons
%             for n = 1:length(isoc_p)-1
%                 plot(ISOCx{m(n)}/1000,ISOCy{m(n)}/1000,'Color',isoc_color)
%             end
%     end
end
