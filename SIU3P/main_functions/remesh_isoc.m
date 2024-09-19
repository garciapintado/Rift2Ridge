function [ISOCHRONS,Basement,Topography] = remesh_isoc(GCOORD, ...
    Point_id,ELEM2NODE,ISOCHRONS,Basement,tp_isoc,high_res)

% Calculate topography
[Topography,~] = find_topo(GCOORD,ELEM2NODE,Point_id);

if nargin==7
    X = ISOCHRONS(1,:);
    X = unique(sort([X Topography(1,:)]));
    Topography = [X; interp1(Topography(1,:),Topography(2,:),X,'linear')];
end

ltopo = size(Topography,2);
if tp_isoc
    Iso_id = unique(ISOCHRONS(3,:));
    ISOCHRONS_NEW = zeros(4,length(Iso_id)*size(Topography,2));
    for n = 1:length(Iso_id)
        ISO_n = ISOCHRONS(:,ISOCHRONS(3,:)==Iso_id(n));
        Indx = (ltopo*(n-1)+1):ltopo*n;
        ISOCHRONS_NEW(1,Indx) = Topography(1,:);
        ISOCHRONS_NEW(2,Indx) = interp1(ISO_n(1,:),ISO_n(2,:), ...
            Topography(1,:),'linear');
        ISOCHRONS_NEW(3,Indx) = Iso_id(n);
        ISOCHRONS_NEW(4,Indx) = interp1(ISO_n(1,:),ISO_n(4,:), ...
            Topography(1,:),'linear');
    end
    ISOCHRONS = ISOCHRONS_NEW;
    %Basement = ISOCHRONS(1:2,ISOCHRONS(3,:)==0);
    Basement = [Topography(1,:); ....
        interp1(Basement(1,:),Basement(2,:),Topography(1,:),'linear')];
else
    Basement = [Topography(1,:); ....
        interp1(Basement(1,:),Basement(2,:),Topography(1,:),'linear')];
end
