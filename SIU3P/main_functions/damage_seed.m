function [I2,WS] = damage_seed(I2,WS,GCOORD,ELEM2NODE,nip)
% [I2,WS] = DAMAGE_SEED(I2,WS,GCOORD,ELEM2NODE,NIP) finds the integration 
% points that belong to the weak seed defined by the parameters on WS 
% structure for a mesh defined by GCOORD and ELEM2NODE, and for nip number 
% of integration points, and outputs the historic second invariant of the
% strain that includes the weak seed damage.

%--------------------------------------------------------------------------
% Function written by Miguel Andres-Martinez, PhD student at Royal Holloway
% University of London, 14-02-2015. Email: mandresmartinez87@gmail.com
%--------------------------------------------------------------------------

%==========================================================================
% INITIALIZE WEAK SEED
%==========================================================================
if size(WS.coord,1) == 1
    switch WS.shape
        case 'cuadrilater'
            X = WS.coord(1)-WS.size(1)/2:WS.res:WS.coord(1)+WS.size(1)/2;
            Y = WS.coord(2)-WS.size(2)/2:WS.res:WS.coord(2)+WS.size(2)/2;
            POINTSX = repmat(X,length(Y),1);
            POINTSY = repmat(Y,1,length(X));
            WS.coord = [POINTSX(:) POINTSY(:)];
        case 'circle'
            % Number of points at the edge of the weak seed
            np = 120;
            % Angular distance between points
            a = 360/np;
            % Angular coordinates
            alpha = 0:a:360;
            % Cartesian coordinates
            outline = repmat(WS.coord,length(alpha),1)+[WS.size*cosd(alpha') ...
                WS.size*sind(alpha')];
            X = WS.coord(1)-WS.size:WS.res:WS.coord(1)+WS.size;
            Y = WS.coord(2)-WS.size:WS.res:WS.coord(2)+WS.size;
            POINTSX = repmat(X,length(Y),1);
            POINTSY = repmat(Y,1,length(X));
            POINTS = [POINTSX(:) POINTSY(:)];
            WS.coord = POINTS(inpolygon(POINTS(:,1),POINTS(:,2), ...
                outline(:,1),outline(:,2)),:);
    end
end

%==========================================================================
% CALCULATE INTEGRATION POINTS
%==========================================================================
% Local IPs
[IP_X,~] = ip_triangle(nip);
% Calculate shape functions
[N,~] = shp_deriv_triangle(IP_X,size(ELEM2NODE,1));

% Reoder the coordinates into element nodes
ECOORD_x = reshape(GCOORD(1,ELEM2NODE),size(ELEM2NODE,1),size(ELEM2NODE,2));
ECOORD_y = reshape(GCOORD(2,ELEM2NODE),size(ELEM2NODE,1),size(ELEM2NODE,2));

% Declare GIPs
IPx = zeros(size(ELEM2NODE,2),nip);
IPy = zeros(size(ELEM2NODE,2),nip);
% Calculate integration point coordinates (GIPs)
for ip=1:nip
    Ni = N{ip};
    GIP_x   = Ni'*ECOORD_x;
    GIP_y   = Ni'*ECOORD_y;
    IPx(:,ip)   = GIP_x;
    IPy(:,ip)   = GIP_y;
end

% Plot (uncomment)
% plot_mesh; hold on
% plot(IPx/1000,IPy/1000,'.')

%==========================================================================
% CREATE WEAK SEED
%==========================================================================
% Calculate distances from the IPs inside the weak seed region
% d2ws = sqrt((IPx-WS.coord(1)).^2 + (IPy-WS.coord(2)).^2);

% Find indexes for the weak seed
WSindx = zeros(size(IPx));
for n = 1:size(WS.coord,1)
    WSindx(sqrt((IPx-WS.coord(n,1)).^2+(IPy-WS.coord(n,2)).^2)<=WS.res) = 1;
end
WSindx = WSindx==1;
%WSindx = inpolygon(IPx,IPy,WS.coord(:,1),WS.coord(:,2));
% WSindx = d2ws<=WS.size;

% Make weak seed
I2(WSindx==1 & I2<WS.damage) = WS.damage;

% % Plot weak seed (uncomment)
% plot(WS.coord(:,1)/1000,WS.coord(:,2)/1000,'.r')
% hold on
% % plot(IPx(WSindx)/1000,IPy(WSindx)/1000,'.b')
% 
% % plot_mesh; hold on
% plot(IPx(WSindx)/1000,IPy(WSindx)/1000,'.b')
