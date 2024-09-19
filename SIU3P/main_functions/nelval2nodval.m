function vn = nelval2nodval(EL2NOD, vne)
    % elval2nodval(EL2NOD, vne)
    % average element-based mesh nodes into mesh nodes
    %
    % INPUT:
    % EL2NOD :: [nnodel,nel] 
    % vne    :: [nnodel,nel]
    % 
    % OUTPUT
    % vn    :: [nnod,1]
    %
    % Author: Javier GP, MARUM, 2020
    
    [nnodel, nel] = size(EL2NOD);
    nelpnod = accumarray(EL2NOD(:), ones(nnodel*nel,1));                   % [nnod,1] 
    vn      = accumarray(EL2NOD(:), vne(:)) ./ nelpnod;                       % [nnod,1)
end % function


    
