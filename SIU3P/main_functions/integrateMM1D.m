function [Inn, Ie] = integrateMM1D(GCOO, val) 
    % numerical integration of a scalar-valued vector
    % values over a 2D polyline (e.g. a boundary line)
    % 
    % GCOO :: [2,nnod] matrix of coordinates 
    % val  :: [1,nnod] for values at polyline nodes, or [nnodel,nel], 
    %                  where 'nnodel' is the number of nodes per 1D element
    %                  and nel is the number of linear elements
    %
    % OUTPUT [O]
    % In: numerical integration of Galerkin weigthed input over the nodes [nnodel,nnod]
    % Ie: element-wise total integration 
    %
    % Details:
    % nnodel is forced to be 3, as this function is specific to be
    % a 1D integration compliant with {6,7} quadratic 2D triangular nodes  
    % GCOO(2,:) == 0 is a valid input, indicating domain is a straight line
    
    nnod = length(GCOO);                                                   % figure(); plot(GCOO(1,:), GCOO(2,:),'-o')
    if mod(nnod,2) ~= 1
        error("integrateMM1D expects nnodel=3, check polyline input")
    end
    
    nnodel = 3;
    nip = 3;
    ipx = [-sqrt(3/5), 0, sqrt(3/5)];                                      % k = 5 maximum polynomial order for exact integration          
    ipw = [5/9; 8/9; 5/9];
    
    nel = ceil(nnod/(nnodel-1)) - 1;
    ELEM2NODE = repmat((1:nnodel)',1,nel) + repmat((nnodel-1)*(0:nel-1),nnodel,1); % [nnodel,nel] in the 1D contour domain
    
    if ~isequal(size(val),size(ELEM2NODE))                                 % [nnodel,nel] values defined individually for each element nnodel
        if length(val) ~= nnod
            error("integrateMM1D:: val must be defined in all nodes or discontinuously at [nnodel,nel]")
        else
            val = val(ELEM2NODE); 
        end
    end                                     
    
    xyv0 = GCOO(:,ELEM2NODE(1,:));                                         % [2,nel] element start corner 
    xyv1 = GCOO(:,ELEM2NODE(nnodel,:));                                    % [2,nel] element end corner
    detJ = sqrt(sum((xyv1-xyv0).^2)) / 2;                                  % [1,nel] |J| = element length / canonical vector length
    
    N = shp_line_int(ipx,nnodel);                                          % cell of shape functions evaluated at integration points
    MMc = zeros(3,3);                                                      % canonical element mass matrix integration
    for ip=1:nip
      Ni = N{ip};                                                          % [nnodel,1] 
      MMc = MMc + ipw(ip) * (Ni * Ni');                                    % [nnodel,nnodel] <- [nnodel,1]*[1,nnodel]
    end                                                                    
    intBasis = MMc * (repmat(detJ,nnodel,1) .* val);                       % [nnodel,nel] <- [nnodel,nnodel]*[nnodel,nel]
    Ie = sum(intBasis);                                                    % [nel,1]
    %Inn = zeros(1,nnod);
    %for i=1:nnodel
    %    elid = ELEM2NODE(i,:);
    %    Inn(elid) =  Inn(elid) + intBasis(i,:);
    %end                                                                    % == accumarray(ELEM2NODE(:), intBasis(:))';
    Inn = accumarray(ELEM2NODE(:), intBasis(:))'; 
end

