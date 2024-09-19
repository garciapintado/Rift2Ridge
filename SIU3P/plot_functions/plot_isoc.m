% PLOT ISOCHRONS OF THE SEDIMENTS
%   Optional input variables:
%
%               USE                     OPTIONS                 DEFAULTS
%
% isoc_type     Type of plot            'isolines' to plot      'isolines'
%                                         sediment isochrons
%                                       'fill' to plot packages
%                                         of sediments with
%                                         colors relating their
%                                         ages
%
% isoc_int      Plotting interval       Interval between saved  1
%                                         steps to plot
%
% isoc_color    Color of the isochrons  1x3 vector              Black
%
% ppar          Parallelisation if any  0 no                    0
%               time-line is a          1 yes
%               parametric function

%--------------------------------------------------------------------------
% Function written by Miguel Andres-Martinez, PhD student at Royal Holloway
% University of London, 05-10-2015. Email: mandresmartinez87@gmail.com
%--------------------------------------------------------------------------

% Edited:
% MA 06-06-2018
%       Added the possibility to calculate erosion in time-lines that are
%       parametric functions and plot them

%==========================================================================
% INPUT
%==========================================================================
% Define the type of plot
if ~exist('isoc_type','var')
    isoc_type = 'isolines';
end

% Define the number of saves between each plotted isochron. Isochrons are
% one per saved file, meaning that they cannot be selected per step but per
% saved step
if ~exist('isoc_int','var')
    isoc_int = 1;
end

% Define color of the isochrons
if ~exist('isoc_color','var')
    isoc_color = [0 0 0];
end

% Define parallel processing options
if ~exist('ppar','var')
    ppar = 0;
end

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
    % Optimised matricial version for time-lines that are function of x
    %======================================================================
    for n = 1:length(isoc_i)
        ISOCy(n,:) = min([ISOCy(n:end,:); Topography(2,:)]);
    end
    
    % Plot
    %----------------------------------------------------------------------
    isoc_p = (isoc_i(2:isoc_int:end)-1) * dt/ma;
    m = 1:isoc_int:length(isoc_i);
    hold on
    
    switch isoc_type
        case 'fill'
            % Plot sediment ages
            for n = 1:length(isoc_p)-1
                patch([ISOCx(m(n),:) ISOCx(m(n+1),end:-1:1)]/1000, ...
                      [ISOCy(m(n),:) ISOCy(m(n+1),end:-1:1)]/1000, isoc_p(n))
            end
            if n~=0
                patch([ISOCx(m(n+1),:) Topography(1,end:-1:1)]/1000, ...
                    [   ISOCy(m(n+1),:) Topography(2,end:-1:1)]/1000,...
                    (istep*dt-2*dt)/ma)
            else
                patch([ISOCx(1,:) Topography(1,end:-1:1)]/1000, ...
                    [   ISOCy(1,:) Topography(2,end:-1:1)]/1000,...
                    (istep*dt-2*dt)/ma)
            end
            colormap('jet')
            hc = colorbar;
            
        case 'isolines'
            % Plot isochrons
            plot(ISOCx(m,:)'/1000,ISOCy(m,:)'/1000,'Color',isoc_color)
    end
else
    %======================================================================
    % Slower cell-based for time-lines that are parametric functions
    %======================================================================
    % Isochrons to plot
    isoc2p = isoc_i(2:isoc_int:end);
    % Initialisation of variables
    ISOCx = {};
    ISOCy = {};
    ISOC = {};
    for n = 1:length(isoc_i)
        ISOC{n} = ISOCHRONS(:,ISOCHRONS(3,:)==isoc_i(n));
    end
    
    % Loop through eroding isochrons
    % ------------------------------
    if ppar==1
        % Parallel loop
        parfor n = 1:length(isoc2p)
            % Load evaluated isochron
            isocx = ISOCHRONS(1,ISOCHRONS(3,:)==isoc2p(n));
            isocy = ISOCHRONS(2,ISOCHRONS(3,:)==isoc2p(n));
            % Find older isochrons
            isoc_above = find(isoc_i>isoc2p(n));
            
            % Run erosion algorithm in serie
            [ISOCx{n},ISOCy{n}] = ero_ISOC_perp(isocx,isocy,ISOC,isoc_above);
        end
    else
        % Serial loop
        for n = 1:length(isoc2p)
            % Load evaluated isochron
            isocx = ISOCHRONS(1,ISOCHRONS(3,:)==isoc2p(n));
            isocy = ISOCHRONS(2,ISOCHRONS(3,:)==isoc2p(n));
            % Find older isochrons
            isoc_above = find(isoc_i>isoc2p(n));
            
            % Run erosion algorithm in serie
            [ISOCx{n},ISOCy{n}] = ero_ISOC_perp(isocx,isocy,ISOC,isoc_above);
            %plot(ISOCx{n}/1000,ISOCy{n}/1000); drawnow; hold on
        end
    end
    
    %======================================================================
    % Plot
    %======================================================================
    isoc_p = isoc_i(2:isoc_int:end)*dt/ma;
    m = 1:length(ISOCx);
    hold on
    
    switch isoc_type
        case 'fill'
            % Plot sediment ages
            for n = 1:length(isoc_p)-1
                patch([ISOCx{m(n)} ISOCx{m(n+1)}(end:-1:1)]/1000, ...
                    [ISOCy{m(n)} ISOCy{m(n+1)}(end:-1:1)]/1000,isoc_p(n))
            end
            if n~=0
                patch([ISOCx{m(n+1)} Topography(1,end:-1:1)]/1000, ...
                    [ISOCy{m(n+1)} Topography(2,end:-1:1)]/1000, ...
                    (istep*dt-dt)/ma)
            else
                patch([ISOCx{1} Topography(1,end:-1:1)]/1000, ...
                    [ISOCy{1} Topography(2,end:-1:1)]/1000, ...
                    (istep*dt-dt)/ma)
            end
            colormap('jet')
            hc = colorbar;
            
        case 'isolines'
            % Plot isochrons
            for n = 1:length(isoc_p)-1
                plot(ISOCx{m(n)}/1000,ISOCy{m(n)}/1000,'Color',isoc_color)
            end
    end
end
