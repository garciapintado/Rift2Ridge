% PLOT BASEMENT
%   Optional input variables:
%
%               USE                     OPTIONS                 DEFAULTS
%
% color_base    Color of the box        #interfaces x 3 vector  Black
%               and interfaces          with values [0 1]
%
% line_widthb   Width of the box        #interfaces x 1 vector  1
%               and interfaces          with width values

%--------------------------------------------------------------------------
% Function written by Miguel Andres-Martinez, PhD student at Royal Holloway
% University of London, 08-09-2015. Email: mandresmartinez87@gmail.com
%--------------------------------------------------------------------------

% Define color for the plot in case color_int doesn't exist
if ~exist('color_base','var')
    color_base = zeros(n_int+1,3);
end

% Define the linewidth for the plot in case line_width doesn't exist
if ~exist('line_widthb','var')
    line_widthb = ones(n_int+1,1);
end

% Calculate topography
Corners = find(Point_id==-1,2,'last');
Topo_base = [GCOORD(:,Point_id==max(Point_id-1)) GCOORD(:,Corners)];
[Topo_base(1,:),itpb] = sort(Topo_base(1,:));
Topo_base(2,:) = Topo_base(2,itpb);

% Cut eroded basement
Basement_ero = cut(Basement,Topo_base);

% Plot
hold on
plot(Basement_ero(1,:)/km,Basement_ero(2,:)/km,'color',color_base, ...
    'LineWidth',line_widthb)
