% Plot isoviscous lines
%   Optional input variables:
%
%               USE                     OPTIONS                 DEFAULTS
%
% isomu        Isoviscosity lines to    #isoviscosity lines x 1 [21]
%              plot                     vector
%
% lwisom        Line width of the       #isoviscous lines x 1   1
%               isoviscous line         vector with width values
%
% cisom         Color of the isoviscous #isoviscous lines x 3   Black
%               line                    vector with values [0 1]

%--------------------------------------------------------------------------
% Function written by Elena Ros-Bernabeu. Edited by Miguel Andres-Martinez,
% 23-09-2015. Email: mandresmartinez87@gmail.com
%--------------------------------------------------------------------------

%==========================================================================
% INPUT
%==========================================================================
% Define isotherm in case Tisotherm doesn't exist
if ~exist('MUisomu','var')
    MUisomu = 21;
end
% Number of isotherms
n_isom = length(MUisomu);
% Define color for the isotherm in case cisot doesn't exist
if ~exist('cisom','var')
    cisom = zeros(n_isom,3);
end

% Define the linewidth for the isotherm in case lwisot doesn't exist
if ~exist('lwisom','var')
    lwisom = ones(n_isom,1);
end

%==========================================================================
% MAPPING FROM INTEGRATION POINTS TO NODES - IP CANNOT BE PLOTTED
%==========================================================================
EL2N      = zeros(nel,3); % new connectivity matrix
GCOORD_N  = zeros(2,nel*3);
Mu_n      = zeros(3,nel);
EL2N(1,:) = 1:3;
nip = 6;
nnodel = 6;
[IP_X, IP_w]    = ip_triangle(nip);
[   Nbig]    = shp_triangle(IP_X, nnodel);
for i=1:nel
    is         = (i-1)*3+1; ie = (i-1)*3+3;
    GCOORD_N(:,is:ie) = GCOORD(:,ELEM2NODE([1 2 3],i));
    EL2N(i,:) = is:ie;
    Dummy      = Nbig'\Mu_c_all(i,:)';
    Mu_n(:,i)= Dummy(1:3);
end

% Correct negative E2_n which are the result of interpolation
Mu_n = abs(Mu_n);

%==========================================================================
% PLOT
%==========================================================================
x_=GCOORD_N(1,:)';
y_=GCOORD_N(2,:)';
xi=linspace(min(x_),max(x_),100);
yi=linspace(min(y_),max(y_),100);
[XI,YI]=meshgrid(xi,yi);
ZI = TriScatteredInterp(GCOORD_N(1,:)',GCOORD_N(2,:)',log10(Mu_n(:)));

hold on
for n = 1:n_isom
    [cs,h] = contour(XI/1000,YI/1000,ZI(XI,YI),[MUisomu(n) MUisomu(n)], ...
        '-','LineWidth',lwisom(n),'Color',cisom(n,:));
    
    %clabel(cs,'FontSize', 15)
    %cc=clabel(cs,h,'Color',cisom(n,:));
end

drawnow
hold off
