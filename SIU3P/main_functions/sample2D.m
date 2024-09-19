function val = sample2D(GCOORD, ELEM2NODE, valn, xy)
    % +++ purpose +++
    % sample a 2D triangular mesh
    %
    % INPUT
    % xy :: [2,p] matrix of sample coordinates
    % 
    % OUTPUT
    % val :: [1,p] samples at xy locations
    %
    % Examples
    % 1] transect
    % xy = [repelem(0,1,500);linspace(-150000,0,500)];
    % xytemp = sample2D(GCOORD, ELEM2NODE, Temp, xy);
    % figure(); plot(xytemp, xy(2,:))
    % 
    % 2] approximate marginal integral over x:
    % X = approxXmarginalIntegrate2D(GCOORD, ELEM2NODE, Temp, GEOn);
    % figure(); plot(X.seq, X.int);
    % 
    % Javier Garcia-Pintado, MARUM 2021

    [nnodel,~] = size(ELEM2NODE);
    nnod = max(ELEM2NODE,[],'all');
    
    if length(GCOORD) > nnod
        GCOORD = GCOORD(:,1:nnod);
    end
    
    xyel = tsearch3(GCOORD, uint32(ELEM2NODE(1:3,:)), xy);                 % parent element
    if any(xyel == 0)
        error("sample2D:: locations out of mesh")
    end
    xyc = getLocalCoo2D(GCOORD,ELEM2NODE(1:3,xyel), xy);                   % [2,p] local coordinates 
    
    N = shp_triangle(xyc', nnodel);                                        % [nnodel,p]
    
    if numel(valn) ~= nnod
        error("sample2D:: valn size not compliant with GCOORD")
    end

    val = sum(N .* reshape(valn(ELEM2NODE(:,xyel)),nnodel,[]));
 
end
