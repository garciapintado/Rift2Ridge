% PLOT SEA LEVEL
%   Optional input variables:
%
%               USE                     OPTIONS                 DEFAULTS
%
% sl_color      Color of the isochrons  1x3 vector              Blue
%
% sl_type       Type of line            'disc'                  'disc'

%--------------------------------------------------------------------------
% Function written by Miguel Andres-Martinez, Research Assistant at
% University of Bremen, 22-08-2017. Email: mandresmartinez87@gmail.com
%--------------------------------------------------------------------------

% Define the type of plot
if ~exist('sl_type','var')
    sl_type = 'disc';
end

% Define color of the isochrons
if ~exist('sl_color','var')
    sl_color = [0 0.4470 0.7410];
end

hold on
plot([min(GCOORD(1,:)) max(GCOORD(1,:))]/1000, ...
    [SP.sealevel SP.sealevel]/1000,'--','Color',sl_color)