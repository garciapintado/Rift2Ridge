% PLASTIC STRAIN RATE AND VISCOUS STRAIN RATE
%   Optional input variables:
%
%               USE                     OPTIONS                 DEFAULTS
%
% etf           etf is a value from 1   [1 0) scalar            1
%               to 0 to that multiplies 
%               the effective creep and 
%               plastic viscosities 
%               before evaluating if  
%               they are bigger or  
%               smaller than the  
%               effective elastic  
%               viscosity (i.e. material 
%               behaves mostly  
%               elastically). Note that  
%               1 is where elastic  
%               effective viscosity is  
%               smaller than other  
%               viscosities but that  
%               doesn't mean that the  
%               material is behaving  
%               exclusively elastic as  
%               the other effectives  
%               viscosities can have a  
%               similar value and then  
%               have a similar  
%               contribution to the  
%               deformation. To show  
%               material where elasticity  
%               plays most of the role etf  
%               can be defined to values  
%               lower than 1.
%
% typeTe        Defines the type of     'fill' a translucid     'fill'
%               plot                    layer covering the
%                                       elastic region
%                                       'line' a contour line
%                                       that separates material
%                                       behaving elastically
%                                       form unelastically
%                                       'both' plots both
%                                       options above
%
% colorTe       Color of the contour    1 x 3 vector with       Red
%               line that separates     values [0 1]
%               material behaving 
%               elastically from
%               unelastically in case
%               typeTe = 'line'
%
% widthTe       Width of the contour    scalar with width       1
%               line in case            values
%               typeTe = 'line'
%
% cfillTe       Color of the elastic    1 x 3 vector with       Black
%               material in case        values [0 1]
%               typeTe = 'fill'
%
% alphaTe       transparency of the     scalar with values      0.3
%               elastic patch           [0 1]
%               plot-element in case
%               typeTe = 'fill'

%--------------------------------------------------------------------------
% Function written by Miguel Andres-Martinez, 12-02-2019. 
% Email: andresma@uni-bremen.de
%--------------------------------------------------------------------------

% Options
% -------
% Define etf if it doesn't exist
if ~exist('etf','var')
    etf = 1;
end

% Define colorTe if doesn't exist
if ~exist('colorTe','var')
    colorTe = [1 0 0];
end

% Define widthTe if doesn't exist
if ~exist('widthTe','var')
    widthTe = 1;
end

% Define typeTe if doesn't exist
if ~exist('typeTe','var')
    typeTe = 'both';
end

% Define alphaTe if doesn't exist
if ~exist('alphaTe','var')
    alphaTe = 0.3;
end

% Define cfillTe if doesn't exist
if ~exist('cfillTe','var')
    cfillTe = [0 0 0];
end

resp_el = 1000;
%==========================================================================
% CALCULATE WHERE ELASTIC DEFORMATION IS DOMINANT
%==========================================================================
% Calcualte ip coordinates
[ix,iy] = ip_coord(GCOORD,ELEM2NODE,size(ELEM2NODE,2),6);
% Make a shear modules elxnip matrix
SHEAR = repmat(Shearm(Phases),1,6);

% Calculate viscosities
Mu_c_all = Mu_dis_all;
diff_s = ~isnan(Mu_dif_all);
Mu_c_all(diff_s) = (1./Mu_dis_all(diff_s)+1./Mu_dif_all(diff_s)).^(-1);
% Find where the ips deform elastically
E = double(etf*Mu_c_all>dt*SHEAR & etf*Mu_b_all>dt*SHEAR);
E(E==0) = 0;
% plot(ix(E==1)/1000,iy((E==1))/1000,'.')
hold on

switch typeTe
    case {'line','both'}
        % Generate a regular mesh
        xx = linspace(min(GCOORD(1,:)),max(GCOORD(1,:)),resp_el);
        yy = linspace(min(GCOORD(2,:)),max(GCOORD(2,:)),resp_el);
        XX = repmat(xx,length(yy),1);
        YY = repmat(yy',1,length(xx));
        
        % Interpolate
        elfunc = scatteredInterpolant(ix(:),iy(:),E(:));
        El_r = elfunc(XX,YY);
        
        % Remove interpolations over the topography
        [Topography,~] = find_topo(GCOORD,ELEM2NODE,Point_id);
        Ytopo = interp1(Topography(1,:),Topography(2,:),XX);
        El_r(YY>Ytopo) = NaN;
end

%==========================================================================
% PLOT
%==========================================================================
MESH.GCOORD = GCOORD/1000;
MESH.EL2NOD = ELEM2NODE;

switch typeTe
    case 'line'
        [cs,cel]=contour(XX/1000,YY/1000,El_r,[0.99,0.99], 'LineWidth', ...
            widthTe,'Color',colorTe);
    case 'fill'
        plot_val(MESH,E,size(ELEM2NODE,2),6,cfillTe)
        axisp = gca;
        axisp.ALim = [0 1/alphaTe];
    case 'both'
        [cs,cel]=contour(XX/1000,YY/1000,El_r,[0.99,0.99], 'LineWidth', ...
            widthTe,'Color',colorTe);
        plot_val(MESH,E,size(ELEM2NODE,2),6,cfillTe)
        axisp = gca;
        axisp.ALim = [0 1/alphaTe];
end