function r = getConductiveHeatFluxBasement(GCOORD, ELEM2NODE, GEOn, Temp, cond_el, Basement, h, dc, downscale)
    % CONDUCTIVE HEAT FLOW ALONG A POLYLINE
    % 
    % +++ purpose +++
    % This functions get the conductive heat flux across a polyline. 
    % The nodes along the polyline are resampled, so that the returning
    % structure indicate a new resampled polyline along with the
    % corresponding heat flux. This function obtains a smooth kernel regression
    % along polyline locations.
    % The functions uses quadratic shape functions to sample the
    % temperature field to estimate the thermal gradients
    %
    % GCOORD    :: [2,nnod] 
    % ELEM2NODE :: [nnodel,nel]
    % GEOn      :: GEOMETRY struct
    % Temp      :: [1,nnod] [ºC]Temperature field
    % cond_el   :: [1,nel] [W.m-1.K-1] thermal conductivity
    % Basement  :: [2,n] : Basement polyline
    % h         :: bandwidth for Gaussian kernel. Empty for automated bandwidth detection
    % dc        :: resampling spacing. Empty for no kernel smoothing
    % downscale :: downscale [approximate] resolution to increase sampling previous to upscaling.
    % 
    % Javier GP, MARUM 2021
    
    % figure(); plot_meshF(ELEM2NODE, GCOORD);
    % figure(); plotGEO(GEOn,1:10,'')
    % hold on; plot(Basement(1,:)/1000,Basement(2,:)/1000,'.-','markersize',1,'color',[.9 .2 .9])      % pink 
    
    if nargin < 7
        h = 500;
    end
    if nargin < 8
        dc = 500;
    end
    if nargin < 9
      downscale = [];
    end
    
    xy = [Basement [GEOn(9).coo(1,:); interp1(Basement(1,:),Basement(2,:),GEOn(9).coo(1,:))]]; % augment with x-topography locations
    [~,ia,~] = unique(xy(1,:));
    xy = xy(:,ia);
    xy(2,:) = min(xy(2,:), interp1(GEOn(9).coo(1,:),GEOn(9).coo(2,:)-5.0,xy(1,:)));            % set sampling polyline is at least [0.5m] deeper than the top of the model
    xy = xy(:,[true diff2D(xy) > 1.0]);                                                           % < 1 meter distance
    r = getConductiveHeatFlux(GCOORD, ELEM2NODE, GEOn, Temp, cond_el, xy, h, dc, downscale);
end
