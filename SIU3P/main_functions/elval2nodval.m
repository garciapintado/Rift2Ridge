function vn = elval2nodval(EL2NOD, ve)
    % elval2nodval(EL2NOD, ve)
    % average element-based mesh nodes into mesh nodes
    %
    % INPUT:
    % EL2NOD:: [nnodel,nel] 
    % ve    :: [1,nel]
    % 
    % OUTPUT
    % vn    :: [nnod,1]
    %
    % Author: Javier GP, MARUM, 2020
    
    [nnodel, nel] = size(EL2NOD);
    nelpnod = accumarray(EL2NOD(:), ones(nnodel*nel,1));                   % [nnod,1]
    tmp     = repmat(ve(:),1,nnodel)';                                     % [nnodel,nel] 
    vn      = accumarray(EL2NOD(:), tmp(:)) ./ nelpnod;                    % [nnod,1)
end