function [Te,TeNod] = Te_RIFT(GCOORD,ELEM2NODE,Point_id,RHEOL,Rho,Ce,G, ...
    ext_rate,km,ma,Phi,Shearm,Temp,filter_type,fwindow,acceptf)
% [TE,TENOD] = TE_RIFT(GCOORD,ELEM2NOE,POINT_ID,RHEOL,RHO,CE,G,EXT_RATE,KM,
% MA,PHI,SHEARM,TEMP,FILTER_TYPE,WINDOWF,ACCEPTF) calculates the elastic 
% thickness using Iskander Muldashev functions develop for Kinedyn (Te_all 
% and Te_Lowry), based on Lowry's codes. The first row of the TE is the 
% x-coordinate of the surface while the second row is the elastic thickness
% for the nodes at the surface. TENOD is a vector elastic thickness defined 
% for all the nodes. The input values are variables from the model except
% FILTER_TYPE,FWINDOW and ACCEPTF. This 3 variables are used to control the
% type of filter applied to the resulting Te where FILTER_TYPE is the
% filter type used in the filtering (see doc smoothdata for filter types)
% and FWINDOW is the window used for filtering. If ACCEPTF is 0 then Te is
% plotted previous to filter and after filtering and the user is asked to
% provide a new window value or accept the filter. If ACCEPTF is 1 the
% results of filtering the Te are returned without consultation.
%
% Note that this code only works for 3 layers and further generalisation
% would be needed to use it with different ammount of layers
%
% Note that this code calculates the elastic thickness using a cutoff for
% decoupled layers of 20 MPa (Burov & Diament, 1995). Currently, it is not 
% possible to change this cutoff value.
%
% Note that TE and TENOD is given in meters
%
% Note sampling is made every 50 m. If models are run with a shift in the
% mess smaller than 100 m then this needs to be changed manually in Te_all
%
% Iskander recommends using a 'rlowess' FILTER_TYPE

%==========================================================================
% DATA TRANSLATION
%==========================================================================
year = ma/1e6;
MESH.GCOORD = GCOORD/km;
MESH.EL2NOD = uint32(ELEM2NODE);
MESH.PointID = Point_id;
MESH.PointID(Point_id==-1) = [1 1 3 3 6 6 9 9];
MESH.xmax = max(GCOORD(1,:))/km;
MESH.xmin = min(GCOORD(1,:))/km;
MESH.PointID_top = 9;
MESH.PointID_bdt = 6;
MESH.PointID_moho = 3;
MESH.PointID_lab = 1;
MESH.PointID_bot = 1;
MESH.PointID_rght = [2 5 8];
MESH.PointID_left = [4 7 10];

PHYSICS.RHEOL = RHEOL;
PHYSICS.Dens = Rho;
PHYSICS.alpha = Ce;
PHYSICS.g = G(2);
PHYSICS.ext_rate = ext_rate*year*1000*2; % mm per year
PHYSICS.km2m = km;
PHYSICS.Myr2s = ma;
PHYSICS.PhiFric0 = [Phi Phi Phi Phi]*180/3.14;
PHYSICS.ShearG = Shearm; % Probably in Pa

VAR.T = Temp;

NUMSCALE = [];

%==========================================================================
% TE KINEDYN FUNCTIONS
%==========================================================================
[MESH]=Te_all(MESH,VAR,PHYSICS,NUMSCALE);
[Topography,~] = find_topo(GCOORD,ELEM2NODE,Point_id);

%==========================================================================
% FILTERING
%==========================================================================
matver = version('-release');
matver = str2double(matver(1:4));
if matver>=2018
    Te = smoothdata(MESH.Te,filter_type,fwindow);
else
    Te = smooth(MESH.Te,fwindow,filter_type);
end
if acceptf~=1
    h = figure;
    while fwindow~=0
        try
            if matver>=2018
                Te = smoothdata(MESH.Te,filter_type,fwindow);
            else
                Te = smooth(MESH.Te,fwindow,filter_type);
            end
        end
        
        clf(h)
        plot(Topography(1,:)/1000,MESH.Te)
        hold on
        plot(Topography(1,:)/1000,Te,'LineWidth',2)
        xlabel('Distance [km]')
        ylabel('Elastic thickness [km]')
        legend('Unfiltered Te','Filtered Te')
        
        fwindow = input(['Would you like to apply a new window for the', ...
            ' filter? [0 -> no, any integer for a new window]: ']);
    end
    close(h)
end

% Output setup
Te = [Topography(1,:); Te(:)'*km];
TeNod = interp1(Te(1,:),Te(2,:),GCOORD(1,:));