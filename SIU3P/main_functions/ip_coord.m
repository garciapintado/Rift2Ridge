function [ipx,ipy] = ip_coord(GCOORD, EL2NOD, nel, nip)
   % [GIP_X,GIP_Y] = IP_COORD(GCOORD, EL2NOD, nel, nip) 
   % Javier GP, MARUM, 2020 
 
   nvert = 3;                                                                % nvert = size(ELEM2NODE,1) if we want to use all nodes
   EL2NOD = EL2NOD(1:nvert,:);

   [ip_x, ~] = ip_triangle(nip);                                          % [nip,2], local (xi,eta) coordinates of integration points
   Nm = shp_triangle(ip_x,nvert);                             % [nvert,nip] matrix shape function evaluated at integration points                                                 
   gcoord_el_x = reshape(GCOORD(1,EL2NOD), nvert, nel);       % [nvert,nel]
   gcoord_el_y = reshape(GCOORD(2,EL2NOD), nvert, nel);       % [nvert,nel]

   % get integration point coordinates
   ipx = gcoord_el_x' * Nm;                                   % [nel,nip]= [nel,nvert]*[nvert,nip]
   ipy = gcoord_el_y' * Nm;
end % function
