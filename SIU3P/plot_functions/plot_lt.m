% PLOT 30 KM CRUSTAL LIMITS
% Plot  vertical lines where there is a change from a thicker than 30 Km
% crust to thinner than 30 km crust. Layer evaluated and thickness can be
% change by declaring the variables below before calling this function.
%
%   Optional input variables:
%
%               USE                     OPTIONS                 DEFAULTS
%
% layer1i       Lower interface index   value of the Point_id   6
%                                       of the interface
%
% layer2i       Upper interface index   value of the Point_id   12
%                                       of the interface
%
% thickness_change
%               Cal        value of the Point_id   6
%                                       of the interface

%--------------------------------------------------------------------------
% Function written by Miguel Andres-Martinez, PhD student at Royal Holloway
% University of London, 08-09-2015. Email: mandresmartinez87@gmail.com
%--------------------------------------------------------------------------

% Define the lower interface in case layer1i doesn't exist
if ~exist('layer1i','var')
    layer1i = 3;
end

% Define the lower interface in case layer1i doesn't exist
if ~exist('layer2i','var')
    layer2i = 9;
end

% Define the given thickness where lines should be drawn in case thickness
% change doesn't exist
if ~exist('thickness_change','var')
    thickness_change = 30*1000;
end

% Sort nodes of the interfaces along the x-axis
LAYER1 = GCOORD(:,Point_id==layer1i);
[LAYER1(1,:),ind] = sort(LAYER1(1,:));
LAYER1(2,:) = LAYER1(2,ind);
LAYER2 = GCOORD(:,Point_id==layer2i);
[LAYER2(1,:),ind] = sort(LAYER2(1,:));
LAYER2(2,:) = LAYER2(2,ind);

% Interpolate the upper interface at the x coordinates of the lower one
L21 = interp1(LAYER2(1,:),LAYER2(2,:),LAYER1(1,:));

% Thickness of the layer between the interfaces minus thickness_change.
% Where D2thick is negative the layer is thinner than thickness_change
D2thick = L21 - (LAYER1(2,:)+thickness_change);
% Indexes of the x-coordinates where the layer is thinner than
% thickness_change
Bool_thick = D2thick<0;

% Find in which nodes there is a change from thicker than thickness_change
% to thinner
Lim_thick = find(sum([Bool_thick 0 0; 0 Bool_thick 0; 0 0 Bool_thick]) ...
    ==2)-1;

% Plot
hold on
for n = 1:length(Lim_thick)
    plot(repmat(LAYER1(1,Lim_thick(n))/1000,1,2), ...
        [LAYER1(2,Lim_thick(n)) L21(Lim_thick(n))]/1000,'r')
end