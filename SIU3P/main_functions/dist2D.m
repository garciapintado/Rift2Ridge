function ans = dist2D(X,Y)
% +++ purpose +++
% parallel Euclidean distance in 2D for two sets of points with identical size
%
% INPUT
% X :: [p,2]   
% Y :: [p,2]  

    if ~isequal(size(X),size(Y))
        error("dist2D:: matrices with different sizes")
    end
    ans = sqrt((X(:,1)-Y(:,1)).^2 + (X(:,2)-Y(:,2)).^2);

end
