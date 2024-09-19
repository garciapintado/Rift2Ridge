function grad_u_xy = sampleGrad2D(GCOORD, ELEM2NODE, u_n, xy)
    % +++ purpose +++
    % sample the gradient of a 2D triangular mesh field at requested locations
    %
    % INPUT
    % xy :: [2,p] matrix of sample coordinates
    % 
    % OUTPUT
    % grad_u_xy :: [2,p] gradient of field u sampled at xy locations
    %
    % Examples
    % 1] transect
    % figure(); plot_tF(Temp,GCOORD, ELEM2NODE)
    % xy = [linspace(min(GCOORD(1,:)),max(GCOORD(1,:)),500); repelem(-6500,1,500)];
    % hold on; plot(xy(1,:)/1000,xy(2,:)/1000,'x')
    % tic; grad_temp_xy = sampleGrad2D(GCOORD, ELEM2NODE(1:3,:), Temp, xy); toc % by direct sampling from the field
    % if 1 > 2
    %     tic; grad_temp = getGradient(ELEM2NODE(1:3,:),GCOORD, Temp);            % by first obtaining the gradient field
    %          grad_temp_xy_x = sample2D(GCOORD, ELEM2NODE(1:3,:), grad_temp(1,:), xy);
    %          grad_temp_xy_y = sample2D(GCOORD, ELEM2NODE(1:3,:), grad_temp(2,:), xy); toc
    % end
    % figure(); plot_tF(grad_temp(2,:),GCOORD, ELEM2NODE)
    % xy_chain = [0 cumsum(diff2D(xy))];
    % figure(); plot(xy_chain, grad_temp_xy(1,:)); hold on; plot(xy_chain, grad_temp_xy_x, 'color', 'red')
    % figure(); plot(xy_chain, grad_temp_xy(2,:)); hold on; plot(xy_chain, grad_temp_xy_y, 'color', 'red')
    %
    % Javier Garcia-Pintado, MARUM 2021

    [nnodel,~] = size(ELEM2NODE);
    nnod = max(ELEM2NODE,[],'all');
    
    if ~ismember(length(u_n), [nnod,size(GCOORD,2)])
        error("sampleGrad2D:: u_n size not compliant with GCOORD")
    end

    if length(GCOORD) > nnod
        GCOORD = GCOORD(:,1:nnod); 
    end
    if length(u_n) > nnod
        u_n = u_n(1:nnod); 
    end

    xyel = tsearch3(GCOORD, uint32(ELEM2NODE(1:3,:)), xy);                 % [1,p] parent element figure(); plot_meshF(EL2NOD, GCOORD)
    if any(xyel == 0)                                                      %                      hold on; plot(xy(1,:)/1000, xy(2,:)/1000,'x','color','blue')
        error("sample2D:: locations out of mesh")                          %                      hold on; plot(xy(1,xyel==0)/1000, xy(2,xyel==0)/1000,'o','color','red')
    end
    
    xyc = getLocalCoo2D(GCOORD,ELEM2NODE(1:3,xyel), xy);                   % [2,p] local coordinates 
    
    [~, dNds]   = shp_deriv_triangle(xyc', nnodel); 
    dNdsm = cell2mat(dNds');      % [nnodel,p*2]
    dNdsm_x = dNdsm(:,1:2:end-1); % [nnodel,p] already transposed (i.e. N_eta')
    dNdsm_z= dNdsm(:,2:2:end);   % [nnodel,p]
    
    [~,dN3ds]   = sf_dsf_tri367(xyc',3,'cell');
    dN3ds       = dN3ds{1}';      

    % get Jacobian
    ndim = 2;
    p = size(xy,2);
    ECOORD_x   = reshape( GCOORD(1,ELEM2NODE(1:3,xyel)), 3, p );  % [nvert=3, p]
    ECOORD_z   = reshape( GCOORD(2,ELEM2NODE(1:3,xyel)), 3, p );  % [nvert=3, p]
    Jx         = ECOORD_x' * dN3ds;                                            % [p, 2] each row : [\frac{\partial x}{\partial \xi} \frac{\partial x}{\partial \eta}]
    Jz         = ECOORD_z' * dN3ds;                                            % [p, 2] each row : [\frac{\partial y}{\partial \xi} \frac{\partial y}{\partial \eta}]
    detJ       = Jx(:,1).*Jz(:,2) - Jz(:,1).*Jx(:,2);                          % [p, 1] 
    if any(detJ<0)
        error('negative |J|')
    end
    invdetJ    = 1./detJ;
    invJx      = zeros(p, ndim);
    invJz      = zeros(p, ndim);
    invJx(:,1) = +Jz(:,2).*invdetJ;
    invJx(:,2) = -Jz(:,1).*invdetJ;
    invJz(:,1) = -Jx(:,2).*invdetJ;
    invJz(:,2) = +Jx(:,1).*invdetJ;

    dNdx =  dNdsm_x .* invJx(:,1)' ...                                     % [nnodel,p]
          + dNdsm_z .* invJx(:,2)';
    dNdz =  dNdsm_x .* invJz(:,1)' ...                                     % [nnodel,p]
          + dNdsm_z .* invJz(:,2)';
    
    grad_u_xy = [sum(dNdx .* reshape(u_n(ELEM2NODE(:,xyel)),nnodel,p)); ...
                 sum(dNdz .* reshape(u_n(ELEM2NODE(:,xyel)),nnodel,p))];
end
