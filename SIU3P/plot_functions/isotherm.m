% Plot a number of isotherms
%   Optional input variables:
%
%               USE                     OPTIONS                 DEFAULTS
%
% Tisotherm     Isotherms to plot       #isotherms x 1 vector   [1340]
%
% lwisot        Line width of the       #isotherms x 1 vector   1
%               isotherm                with width values
%
% cisot         Color of the isotherm   #isotherms  x 3 vector  Black
%                                       with values [0 1]

%--------------------------------------------------------------------------
% Function written by Elena Ros-Bernabeu. Edited by Miguel Andres-Martinez,
% 18-02-2015. Email: mandresmartinez87@gmail.com
%--------------------------------------------------------------------------

%==========================================================================
% INPUT
%==========================================================================
% Define isotherm in case Tisotherm doesn't exist
if ~exist('Tisotherm','var')
    Tisotherm = 1340;
end
% Number of isotherms
n_isot = length(Tisotherm);
% Define color for the isotherm in case cisot doesn't exist
if ~exist('cisot','var')
    cisot = zeros(n_isot,3);
end

% Define the linewidth for the isotherm in case lwisot doesn't exist
if ~exist('lwisot','var')
    lwisot = ones(n_isot,1);
end

%==========================================================================
% MAPPING FROM INTEGRATION POINTS TO NODES - IP CANNOT BE PLOTTED
%==========================================================================
EL2N      = zeros(nel,3); % new connectivity matrix
GCOORD_N  = zeros(2,nel*3);
T_n      = zeros(3,nel);
EL2N(1,:) = 1:3;
nip = 6;
nnodel = 6;
[IP_X, IP_w]    = ip_triangle(nip);
[   Nbig]    = shp_triangle(IP_X, nnodel);
for i=1:nel
    is         = (i-1)*3+1; ie = (i-1)*3+3;
    GCOORD_N(:,is:ie) = GCOORD(:,ELEM2NODE([1 2 3],i));
    EL2N(i,:) = is:ie;
    T_n(is:ie) = Temp(ELEM2NODE([1 2 3],i));    
end

%==========================================================================
% PLOT
%==========================================================================
x_=GCOORD_N(1,:)';
y_=GCOORD_N(2,:)';
z_=T_n;
xi=linspace(min(x_),max(x_),100);
yi=linspace(min(y_),max(y_),100);
[XI,YI]=meshgrid(xi,yi);
ZI = TriScatteredInterp(GCOORD(1,:)',GCOORD(2,:)',Temp(:));

hold on
for n = 1:n_isot
    [cs,h] = contour(XI/1000,YI/1000,ZI(XI,YI),[Tisotherm(n) Tisotherm(n)], ...
        '--','LineWidth',lwisot(n),'Color',cisot(n,:));
    
    %clabel(cs,'FontSize', 15)
    cc=clabel(cs,h,'Color',cisot(n,:));
end

drawnow
hold off
