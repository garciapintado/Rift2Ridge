% PLOT_OPT generates an structure SETp that can be used as input for
% plot_series function.

%--------------------------------------------------------------------------
% Function written by Miguel Andres-Martinez, PhD student at Royal Holloway
% University of London, 18-02-2015. Email: mandresmartinez87@gmail.com
%--------------------------------------------------------------------------

% For DEFAULTS leave it EMPTY []

% VISIBILITY
% 'visible' for visualizing the plots while generated
% 'invisible' for not visualizing the plots while generated
SETp.visible = 'visible';

% POSITION
% [] not defined position
% 1x4 position vector
SETp.position = [];

% AXIS
% 'equal'/'tight' or a [minX maxX minY maxY] vector
SETp.axis = ([-150 150 -130 5]);
    
% ASPEC RATIO
% 1x3 vector with values between 0 and one
% (i.e. [1 0.5 1] for a vertical exageration of 2)
SETp.aspect = [1 1 1];
    
% COLOR OF THE INTERFACES
% (#interfaces+1)x3 vector with values between 0 and 1
% each row defines the color of an interface
% (first row defines the color of the box of the model)
SETp.color_int = zeros(4,3); % SETp.color_int([2 4],:) = 1;
    
% WIDTH OF THE INTERFACES LINES
% (#interfaces+1)x1 vector with values for LineWidth for each interface
% (first row defines the LineWidth of the box of the model)
SETp.line_width = [2; 2; 2; 1];

% LOG DATA
% plot the logaraithm of the data
% (only for strain rate and viscosity)
SETp.log_data = 1;

% CONSTANT COLORS
% 1 in case the color scale is needed to be constant along the time steps
% taking into account the maximum and the minimum of the value along the 
% time steps
% 0 in case constant color scale is not needed, so the color scale would be
% calculated in base of the maximum and minimum at each time step
% [min_value max_value] in case the input range needs to be introduced
% manually
SETp.ct_colors = [-22 -13.5];

% PLOT ISOTHERM
% #isotherm x 1 vector with values for isotherm temperatures
% [] do not plot isotherm
SETp.isot = [1340];

% ISOTHERM COLOR
% (#interfaces+1)x3 vector with values between 0 and 1
% each row defines the color of an isotherm
SETp.cisot = zeros(1,3);

% WIDTH OF THE ISOTHERM
% (#interfaces+1)x1 vector with values for LineWidth for each isotherm
SETp.lwisot = [1];

% PLOT ISOVISCOUS LINES
% #isoviscous lines x 1 vector with values for viscosities
% [] do not plot isotherm
SETp.MUisomu = [21];

% ISOVISCOSITY LINE COLOR
% #isoviscous lines x3 vector with values between 0 and 1
% each row defines the color of an isoviscous line
SETp.cisom = zeros(1,3);

% WIDTH OF THE ISOVISCOUS LINE
% #isoviscous lines x1 vector with values for LineWidth for each isoviscous
% line
SETp.lwisom = [1];

% PLOT TRACK OBJECTS
% 1 plot track objects
% 0 do not plot track objects
SETp.tp = 0;

% WAY OF PLOTTING TRACK OBJECTS
% 'elements' to plot the track points in elements
% 'layers' to plot the trackbpoints by layers
SETp.dg_option = [];

% PLOT WEAK SEED
% 1 plot weak seed
% 2 do not plot weak seed
SETp.ws = 0;

% PLOT BASEMENT
% 1 plot basement in case model is runned with sediments
% 0 do not plot the sediments
SETp.basement = 0;

% COLOR BASEMENT/SEDIMENT INTERFACE
% 1x3 vector with values between 0 and 1
SETp.color_base = [1 0 0];

% LINE WIDTH OF THE BASEMENT/SEDIMENT INTERFACE
% Scalar With value
SETp.line_widthb = 1;

% PLOT TIME LINES
% 'isolines' to plot sediment time lines
% 'fill' to plot packages of sediments with colores relating their ages
SETp.isoc_type = 'fill';

% TIME LINES INTERVAL BETWEEN SAVED STEPS
% Scalar
SETp.isoc_int = 10;

% TIME LINES COLOR
% 1x3 vector
SETp.isoc_color = [0 0 0];

% TITLE
% 'none' no title
% 'default' default as in the given functions
% 'timestep' time step
SETp.title = 'timestep';

% COLORBAR
% 'on'
% 'off'
SETp.colorbar = 'on';

% TIMELINES IN SEDIMENTS
SETp.timelines = 1;

% SETp.visible = 'visible';
% SETp.axis = 'tight';
% SETp.aspect = [1 0.5 1];
% SETp.color_int = zeros(4,3);
% SETp.line_width = [1; 1; 1; 1];
% SETp.log_data = 1;
% SETp.ct_colors = 0;
% SETp.tp = 0;
% SETp.dg_option = [];