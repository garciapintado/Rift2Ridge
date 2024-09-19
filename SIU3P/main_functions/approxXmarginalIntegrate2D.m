function X = approxXmarginalIntegrate2D(GCOORD, ELEM2NODE, valn, GEOn, nx, ny, dd)
    % function intx = approxXmarginalIntegrate2D(GCOORD, ELEM2NODE, valn, GEOn, nx, ny, dd)
    % +++ purpose +++ 
    % approximate marginal integral in 2D Euclidean space over x [horizontal] dimension
    % The function just estimates the integral (X.int) over a series of vertical
    % transects placed at specific x locations (X.seq). The result is such
    % integration of the {X.seq,Xint} values by the trapezoidal rule
    % exactly matches the 2D integration of the input variable over the domain
    % 
    %
    % INPUT
    % GCOORD    :: REAL [2,nnod]
    % ELEM2NODE :: INTEGER [nnodel,nel]
    % valn      :: REAL [1,nnod]
    % GEOn      :: struct, mesh geometry
    % nx        :: INTEGER, OPTIONAL, (default nx=500) number of points in the x direction
    % ny        :: INTEGER, OPTIONAL, (default nx=200) number of points in the y direction
    % dd        :: REAL, OPTIONAL. default to 1.0E-03. margin to add to border domains to avoid tsearch3() sampling errors at the bundaries
    %
    % OUTPUT
    % X.
    %  .seq     :: REAL [1,nx] marginal integral over the horizontal dimension
    %  .int     :: REAL [1,nx]     
    %      
    % Javier Garcia-Pintado, MARUM, 2021
    
    X = [];
    
    if nargin < 5
        nx = 500;
    end
    if nargin < 6
        ny = 200;
    end
    if nargin < 7
        dd = 1.0E-03;
    end
    
    

    X.seq = linspace(GEOn(1).coo(1,1)+dd,GEOn(1).coo(1,end)-dd,nx);         % [1,nx]
    ytop = interp1(GEOn(end-1).coo(1,:),GEOn(end-1).coo(2,:),X.seq)-dd;     % [1,nx] topography minus dd perturbation to avoid out-of-domain samples
    ybot = GEOn(1).coo(2,1)+dd;
    x = repmat(X.seq,ny,1);                                                 % [ny,nx]
    ystd = linspace(1,0,ny);                                                % [1,ny] standarised y sequence
    y = ystd' * (ytop-ybot) + ybot;
    xyval = reshape(sample2D(GCOORD, ELEM2NODE, valn, [x(:)';y(:)']),ny,nx); % [ny,nx] figure(); plotGEO(GEOn,1:10,''); hold on; plot(x/1000,y/1000,'+')
    X.int = sum(0.5 * (xyval(1:end-1,:) + xyval(2:end,:)) .* (- diff(y)));   % [nx,1] integral over x 
    
    Xmint = sum(diff(X.seq) .* (X.int(1:end-1) + X.int(2:end)) / 2);
    if Xmint ~= 0
        X.int = X.int * (integrate2D(GCOORD,ELEM2NODE,valn)/Xmint);
    end
end % function approxXmarginalIntegrate2D

  
