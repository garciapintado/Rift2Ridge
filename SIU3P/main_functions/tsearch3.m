function xyel = tsearch3(GCOORD, ELEM2NODE, xy, p)
    % +++ purpose +++
    % search elements on which xy samples lie.
    % For relatively big datasets a python-based tsearch is executed.
    % For big datasests, python is way faster than mutils::tsearch2() if no triangulation is given as tsearch2() input
    %
    % INPUT
    % GCOORD    :: REAL [2,nnod]
    % ELEM2NODE :: INTEGER [3,nel]
    % xy        :: REAL [2,ncoo]
    % p         :: REAL, perturbation added to xy to re-attempt allocation of non-allocated parent elements
    %
    % OUTPUT
    % xyel      :: [1,ncoo]
    %
    % Details:
    % Often, tsearchPy and/or tsearch2, both called from this function, fail in locating parent triangles at points
    % falling exactly at vertex or edge location. A perturbation (p>0.0) serves to re-try allocation of these failed searchs.
    %
    %
    % Author: Javier Garcia-Pintado, MARUM, 2021
    
    if nargin < 4
        p = 0.001;             % 1 mm
    end

    if size(xy,2) <= 35000     % size for which similar matlab vs python performance
        xyel = tsearch2(GCOORD, uint32(ELEM2NODE(1:3,:)), xy);
    else
        try                    % big dataset
            xyel = tsearchPy(GCOORD, uint32(ELEM2NODE(1:3,:)), xy, true);
        catch
            try
                disp('tsearch3: tsearchPy() error: using unverified .mat')
                xyel = tsearchPy(GCOORD, uint32(ELEM2NODE(1:3,:)), xy, false); % hold on; plot(GCOORD(1,xyel==0)/1000,GCOORD(2,xyel==0)/1000,'x','markersize',15,'color','red')
            catch    
                disp('tsearch3: tsearchPy() error: using matlab tsearch2()')
                xyel = tsearch2(GCOORD, uint32(ELEM2NODE(1:3,:)), xy);
            end
        end
    end

    if p ~= 0.0
        xyouboo = xyel == 0;
        if any(xyouboo)
            xmin = min(GCOORD(1,:)) + p;
            xmax = max(GCOORD(1,:)) - p;
            xy(1,:) = max(xmin,xy(1,:));
            xy(1,:) = min(xmax,xy(1,:));
            xy(2,:) = xy(2,:) - p;
            xyel(xyouboo) = tsearch2(GCOORD, uint32(ELEM2NODE(1:3,:)), xy(:,xyouboo));
        end
        xyouboo = xyel == 0;
        if any(xyouboo)
            xy(2,:) = xy(2,:) + 2*p;
            xyel(xyouboo) = tsearch2(GCOORD, uint32(ELEM2NODE(1:3,:)), xy(:,xyouboo));
        end
    end
end
