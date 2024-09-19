% PLOT EROSIONED MATERIAL
%   Optional input variables:
%
%               USE                     OPTIONS                 DEFAULTS
%
% color_ero     Color of the box        1 x 3 vector with       Light brown
%               and interfaces          values [0 1]
%
% line_ero      Style of the line       style line string       --
%
% line_ero_c    Color of the line       1 x 3 vector with       Red
%                                       values [0 1]

%--------------------------------------------------------------------------
% Function written by Miguel Andres-Martinez, Postdoc at University of
% Bremen, 13-02-2018. Email: andresma@uni-bremen.de
%--------------------------------------------------------------------------

% Define color for the plot in case color_ero doesn't exist
if ~exist('line_ero_c','var')
    color_ero = [0.8900    0.8450    0.8050];
end

% Define type of line style for the plot in case line_ero doesn't exist
if ~exist('line_ero','var')
    line_ero = '--';
end

% Define color for the line in case line_ero_c doesn't exist
if ~exist('line_ero_c','var')
    line_ero_c = [0.9020 0.2000 0];
end

Ero_p = [ISOCHRONS(1,ISOCHRONS(3,:)==0) GCOORD(1,Corner_id([2 1])); ...
    ISOCHRONS(2,ISOCHRONS(3,:)==0) GCOORD(2,Corner_id([2 1]))];
patch(Ero_p(1,:)/1000,Ero_p(2,:)/1000,ones(1,size(Ero_p,2)),'FaceColor', ...
    color_ero,'LineStyle',line_ero,'EdgeColor',line_ero_c)