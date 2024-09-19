function [ISOCx,ISOCy_ini,ISOCyy,ISOCpwd,ero_bool]=find_ero(ISOCHRONS,Topography)
% 

% Erode sediments
isoc_i = unique(ISOCHRONS(3,:));
ISOCHRONSpp = ISOCHRONS(1:2,ismember(ISOCHRONS(3,:),isoc_i));
ISOCHRONSpwd = ISOCHRONS(4,ismember(ISOCHRONS(3,:),isoc_i));
ISOCx = reshape(ISOCHRONSpp(1,:),sum(isoc_i(1)==ISOCHRONS(3,:)), ...
    length(isoc_i))';
ISOCy_ini = reshape(ISOCHRONSpp(2,:),sum(isoc_i(1)==ISOCHRONS(3,:)), ...
    length(isoc_i))';
ISOCpwd   = reshape(ISOCHRONSpwd,sum(isoc_i(1)==ISOCHRONS(3,:)), ...
    length(isoc_i))';

for n = 1:length(isoc_i)
    ISOCyy(n,:) = min([ISOCy_ini(n:end,:);Topography(2,:)]);
end
ero_bool = (ISOCy_ini~=ISOCyy);


