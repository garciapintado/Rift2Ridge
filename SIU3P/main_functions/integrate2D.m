function [Id, Ie] = integrate2D(GCOORD, ELEM2NODE, val, mask, plot_error, getmean)
    % +++ purpose +++ 
    % numerical integration of field values over a 2D domain
    % Input can be defined at nodes or at integration points in a 2D FEM mesh.
    % 
    % GCOORD    :: [2,nnod] matrix of coordinates 
    % ELEM2NODE :: [nnodel,nel] mesh topology, 'nnodel' is number of nodes per element, and 'nel' is number of elements
    % val       :: [1,nnod] for values defined at nodes, or 
    %              [nnodel,nel] for [possibly] discontinuous values defined at each element, or
    %              [nel,6] is assumed to be defined at integration point
    % mask      :: OPT, LOGICAL [1,nel] If provided, only the region with true elements in mask
    %              is considered for the output. Empty array for no mask.
    %  plot_error :: LOGICAL
    %
    % getmean   :: LOGICAL, default to False. If True, the mean rather than the integral is returned
    %
    % OUTPUT [O]
    % Id: domain wide numerical integration
    % Ie: element-wise integration 
    %
    
    % Author: Javier García-Pintado, MARUM, 2020
    
    nnod = max(ELEM2NODE,[],'all');                                            % number of nodes
    nnodel = size(ELEM2NODE,1);                                                % number of nodes per element [default=7 in rift2ridge2D]
    nel    = size(ELEM2NODE,2);                                                % number of elements
    
    if nargin < 4 || isempty(mask)
        mask = true(1,nel);
    end
    
    if nargin < 5
        plot_error = true;
    end
    if nargin < 6
        getmean = false;
    end

    if length(mask) ~= nel
        error("length(mask) ~= nel")
    end

    if nnodel > 3
       nip = 6;
    else
       nip = 3;
    end
    
    ndim  = 2;
    nvert = 3;

    [x_ip, ipw] = ip_triangle(nip);                                        % IP_X [nip,2], IP_w [1,nip] with nip=6; 
    [N, dNdu]   = shp_deriv_triangle(x_ip, nnodel);                        % N [nip] cell array, with [nnodel,1] size per cell. Shape function evaluated at integration points 
    Nm = reshape(cell2mat(N),nnodel,nip);                                  % [nnodel,nip]
    
    [~, dN3ds]   = sf_dsf_tri367(x_ip',3,'cell');                          % dNdu [nip] cell array, with [nnodel,2] size per cell. Derivative od shape functions at integration points
    dN3ds        = dN3ds{1}';                                              % [nnodel3,2]
  
    if nel == nip
        error("integrate2D: input dimensions not identifiable")
    end
    if isequal(size(val),[nel,nip])
        ipval = val';                                                      % [nip,nel]
    else
        if ~isequal(size(val),size(ELEM2NODE))                             % [nnodel,nel] values defined individually for each element nnodel
            if length(val) ~= nnod
                error("integrate2D:: val must be defined in all nodes or discontinuously at [nnodel,nel]")
            else
                valn = val;
                val = val(ELEM2NODE); 
            end
        end
        ipval = Nm' * val;                                                 % [nip,nel]
    end
    
    nel = sum(mask);
    ELEM2NODE = ELEM2NODE(:,mask);
    ipval = ipval(:,mask);

    nelblk = nel;
    ECOORD_x = reshape( GCOORD(1,ELEM2NODE(1:nvert,:)), nvert, nelblk);    % REAL [nnodel, nelblk]
    ECOORD_y = reshape( GCOORD(2,ELEM2NODE(1:nvert,:)), nvert, nelblk);    % REAL [nnodel, nelblk]
    Jx       = ECOORD_x' * dN3ds;                                          % [nelblk, 2] each row : [\frac{\partial x}{\partial \xi} \frac{\partial x}{\partial \eta}]
    Jy       = ECOORD_y' * dN3ds;                                          % [nelblk, 2] each row : [\frac{\partial y}{\partial \xi} \frac{\partial y}{\partial \eta}]
    detJ = reshape(Jx(:,1).*Jy(:,2) - Jy(:,1).*Jx(:,2),1,nelblk);          % [1,nel] \frac{\partial x}{\partial \xi} * \frac{\partial y}{\partial \eta} - \frac{\partial y}{\partial ji} * \frac{\partial x}{\partial eta}
    if any(detJ < 0)
        if exist("valn","var") && plot_error
            elboo = detJ <= 0;
            find(elboo)
            figure(); plotN(valn, GCOORD, ELEM2NODE, [], "", 0);
            hold on; plot_meshF(ELEM2NODE, GCOORD, [], 'white');
            hold on; plot_meshF(ELEM2NODE, GCOORD, elboo, 'red',3);
            minmax(GCOORD(:,unique(ELEM2NODE(:,elboo))))
        end
        error('negative |J|');
    end
    
    Ie = sum(ipval .* repmat(ipw(:),1,nel)) .* detJ;                       % [1,nel] element integrals
    Id = sum(Ie);

    if getmean
        Ie = sum(ipval .* repmat(ipw(:),1,nel)) .* 2; 
        Id = Id / sum(detJ) * 2;
    end
end % function integrate2D

% example: not run
% [Tempbar, Tempel] = integrate2D(GCOORD, ELEM2NODE, Temp,[],true,true);
