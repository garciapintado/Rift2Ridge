function xyc = getLocalGoo2D(GCOORD,ELEM2NODE, xy)
    % xyc = local_coords_2d(GCOORD,EL2NOD,els,gX_PT)
    % +++ purpose +++
    % return canonical [local] coordinates global xy locations
    % 
    % INPUT
    %   GCOORD    :: REAL [2,nnod]    : coordinates of all nodes in mesh
    %   ELEM2NODE :: [nnodel,p] : connectivity matrix in the element parent for each xy
    %   xy        :: REAL [2,p], global point coordinates
    %
    % OUTPUT
    %   xyc     :: REAL [2,p], local coordinates of points within the corresponding parent element
    %
    % Javier Garcia-Pintado, MARUM, 2021


    [nnodel,p] = size(ELEM2NODE);
    if p ~= size(xy,2)
        error("getLocalCoo2D:: ELEM2NODE and xy non-compliant sizes")
    end

    x = reshape(GCOORD(1,ELEM2NODE(1:3,:)),3,[]);  % [3,p]
    y = reshape(GCOORD(2,ELEM2NODE(1:3,:)),3,[]);  % [3,p]

    xp  = xy(1,:);
    yp  = xy(2,:);

    xyc      = nan(2,p);
    denom    = -x(1,:).*y(3,:)+x(1,:).*y(2,:)-x(2,:).*y(1,:)+x(2,:).*y(3,:)+x(3,:).*y(1,:)-x(3,:).*y(2,:);
    xyc(1,:) = -(-x(1,:).*yp+x(1,:).*y(3,:)-x(3,:).*y(1,:)+xp.*y(1,:)+x(3,:).*yp-xp.*y(3,:)) ./ denom;
    xyc(2,:) = (x(1,:).*y(2,:)-x(1,:).*yp-x(2,:).*y(1,:)+x(2,:).*yp+xp.*y(1,:)-xp.*y(2,:)) ./ denom;

end % function getLocalCoo2D
