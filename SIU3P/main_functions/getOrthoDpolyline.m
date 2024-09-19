function [eortho, northo] = getOrthoDpolyline(xy, clockwise, rotangle)
    % +++ purpose +++
    % get directions orthogonal to polyline segments. Direction are angles in [0,360] with 0 pointing to the east
    % by default direction are deviated towards the right (clockwise) with
    % respect to segment directions
    %
    % xy        :: REAL [2,ncoo] polyline coordinates 
    % clockwise :: LOGICAL, orthogonal directions to the right of segment directions
    % rotangle  :: REAL, OPTIONAL. To get directions other than normal to the polyline segments. 
    %              For example rotangle=0 just gives the direction of segments along the polyline 
    %
    % OUTPUT:
    % eortho :: orthogonal direction to polyline segments
    % northo :: orthogonal directions at polylines nodes as an average between consecutive segments
    

    % Author: Javier GP, MARUM 2020
   
    if nargin < 3
        if clockwise
            rotangle = -90;
        else
            rotangle = 90;
        end
    end
    
    diffx = diff(xy(1,:));
    diffy = diff(xy(2,:));
    alphas = atand(diffy./diffx);
    plus180 = (diffx < 0 & diffy < 0) | (diffx < 0 & diffy > 0);
    alphas(plus180) = alphas(plus180) + 180;
    alphas(alphas < 0) = alphas(alphas < 0) + 360;                          % plot(xy(1,:),xy(2,:),'.-','markersize',15)
                                                                            % 
    
    eortho = alphas + rotangle;
    eortho(eortho > 360) = eortho(eortho > 360) - 360;
    eortho(eortho < 0)   = eortho(eortho < 0)   + 360;
    
    dalphas = diff(alphas);                                                 % hold on; text(xc, yc, string(round(alphas,1)))        
    dalphas(dalphas < 0.) = dalphas(dalphas < 0.) + 360;                    % hold on; text(xy(1,1:end-1), xy(2,1:end-1),string(alphas))
    mirror = dalphas > 180.;
    nalphas = alphas(1:end-1) + dalphas/2. ;
    nalphas(nalphas > 360) = nalphas(nalphas > 360) - 360;
    nalphas(mirror) = nalphas(mirror) + 180;
    nalphas = [alphas(1) nalphas alphas(end)];                      % alphas at nodes as average between consecutive segments
    nalphas(nalphas > 360) = nalphas(nalphas > 360) - 360;
    
    northo = nalphas + rotangle;
    northo(northo > 360) = northo(northo > 360) - 360;
    northo(northo < 0)   = northo(northo < 0)   + 360;                      % hold on; text(xy(1,:), xy(2,:),string(round(northo,2)))
    
    % test
    if 1 > 2 % not run
        xy = rand(2,10);
        [eortho, northo] = getOrthoDpolyline(xy);
        diffx = diff(xy(1,:));
        diffy = diff(xy(2,:));
        
        symsize = 0.5E04;
        xc = xy(1,1:end-1) + diffx/2; yc = xy(2,1:end-1) + diffy/2;
        ecs = [cosd(eortho) ; sind(eortho)];                                % polyline edges
        ncs = [cosd(northo) ; sind(northo)];                                % polyline nodes
        figure(); plot(xy(1,:),xy(2,:),'.-','markersize',10)
        hold on; plot(xy(1,1),xy(2,1),'o','color','red')                    % start
        hold on; plot(xc,yc,'o','color','blue')   
        for i=1:length(eortho)                                              % orthogonal angles at edges
            hold on; plot(xc(i)+[0. ecs(1,i)*symsize], yc(i)+[0. ecs(2,i)*symsize],'-','color','red'); 
        end
        for i=1:length(northo)                                              % orthogonal angles at edges
            hold on; plot(xy(1,i)+[0. ncs(1,i)*symsize], xy(2,i)+[0. ncs(2,i)*symsize],'-','color','green'); 
        end
        text(xy(1,:), xy(2,:)+0.01,string(round(northo,2)),'color','red') 
    end %test
end % function
