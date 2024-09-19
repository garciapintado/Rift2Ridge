% PLOT ASTHENOSPHERE
%   Optional input variables:
%
%               USE                     OPTIONS                 DEFAULTS
%
% color_asth    Color of the box        1 x 3 vector with       Light brown
%               and interfaces          values [0 1]
%
% line_asth     Style of the line       style line string       -
%
% line_asth_c   Color of the line       1 x 3 vector with       Black
%                                       values [0 1]

%--------------------------------------------------------------------------
% Function written by Miguel Andres-Martinez, Postdoc at University of
% Bremen, 13-02-2018. Email: andresma@uni-bremen.de
%--------------------------------------------------------------------------

% Define color for the plot in case color_asth doesn't exist
if ~exist('color_asth','var')
    color_asth = [1.0000    0.5    0.3];
end

% Define type of line style for the plot in case line_asth doesn't exist
if ~exist('line_asth','var')
    line_asth = '-';
end

% Define color for the line in case line_asth_c doesn't exist
if ~exist('line_asth_c','var')
    line_asth_c = [0 0 0];
end

ASTH = [TRACKP(1,:) GCOORD(1,Corner_id([2 1])); ...
    TRACKP(2,:) GCOORD(2,Corner_id([2 1]))];
patch(ASTH(1,:)/1000,ASTH(2,:)/1000,ones(1,size(ASTH,2)),'FaceColor', ...
    color_asth,'LineStyle',line_asth,'EdgeColor',line_asth_c)