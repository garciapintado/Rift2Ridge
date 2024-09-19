% PLOT BASEMENT
%   Optional input variables:
%
%               USE                     OPTIONS                 DEFAULTS
%
% color_base    Color of the box        1 x 3 vector with       Black
%               and interfaces          values [0 1]
%
% line_widthb   Width of the box        width value             1
%               and interfaces

%--------------------------------------------------------------------------
% Function written by Miguel Andres-Martinez, PhD student at Royal Holloway
% University of London, 08-09-2015. Email: mandresmartinez87@gmail.com
%--------------------------------------------------------------------------

% Define color for the plot in case color_int doesn't exist
if ~exist('color_base','var')
    color_base = zeros(1,3);
end

% Define the linewidth for the plot in case line_width doesn't exist
if ~exist('line_widthb','var')
    line_widthb = ones(1,1);
end

% Remesh isochrons to have the same x-points as the topography
[ISOCHRONSp,~] = remesh_isoc(GCOORD,Point_id,ELEM2NODE,ISOCHRONS, ...
    Basement,tp_isoc);

% Calculate topography
[Topography,~] = find_topo(GCOORD,ELEM2NODE,Point_id);

% Erode basement
isoc_i = unique(ISOCHRONSp(3,:));
ISOCHRONSpp = ISOCHRONSp(1:2,ismember(ISOCHRONSp(3,:),isoc_i));
ISOCx = reshape(ISOCHRONSpp(1,:),sum(isoc_i(1)==ISOCHRONSp(3,:)), ...
    length(isoc_i))';
ISOCy = reshape(ISOCHRONSpp(2,:),sum(isoc_i(1)==ISOCHRONSp(3,:)), ...
    length(isoc_i))';

Basement_ero = [ISOCx(1,:); min([ISOCy(1:end,:); Topography(2,:)])];

% Plot
hold on
plot(Basement_ero(1,:)/km,Basement_ero(2,:)/km,'color',color_base, ...
    'LineWidth',line_widthb)
