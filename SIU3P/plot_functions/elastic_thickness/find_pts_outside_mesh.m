function [iout,gX_PT] = find_pts_outside_mesh(MESH,gX_PT)

FigNo   = 0;
verbose = 0;

nPT     = size(gX_PT,2);
iout    = true(1,nPT);

inod_top  = ismember(MESH.PointID,MESH.PointID_top);
inod_bot  = ismember(MESH.PointID,MESH.PointID_bot);
inod_rght = ismember(MESH.PointID,MESH.PointID_rght);
inod_left = ismember(MESH.PointID,MESH.PointID_left);

z_top     = MESH.GCOORD(2,inod_top);
z_bot     = MESH.GCOORD(2,inod_bot);
x_rght    = MESH.GCOORD(1,inod_rght);
x_left    = MESH.GCOORD(1,inod_left);

% Fit a rectangular box into the mesh and mark all points within this 
% region as beeing "in"
ztol      = 0.01;
xtol      = 0.01;
zmin_top  = min(z_top);
zmax_bot  = max(z_bot);
xmin_rght = min(x_rght);
xmax_left = max(x_left);
iin       = gX_PT(1,:)>=xmax_left+xtol & gX_PT(1,:)<=xmin_rght-xtol & ...
            gX_PT(2,:)<=zmin_top -ztol & gX_PT(2,:)>=zmax_bot +ztol;
iout(iin) = 0;

% Check remaining points
ind       = find(iout);
if ~isempty(ind)
    % CHECK BOUNDARIES BY CREATING A CLOSED POLYGON ALONG THE DOMAIN
    % BOUNDARIES; THEN CALL "inpolygon"
    
    % sort bottom boundary nodes such that x increases
    [x_bot,p]  = sort(MESH.GCOORD(1,inod_bot),'ascend');
    z_bot      = z_bot(p);
    
    % sort rght boundary nodes such that z increases
    [z_rght,p] = sort(MESH.GCOORD(2,inod_rght),'ascend');
    x_rght     = x_rght(p);

    % sort top boundary nodes such that x decreases
    [x_top,p]  = sort(MESH.GCOORD(1,inod_top),'descend');
    z_top      = z_top(p);
    
    % sort left boundary nodes such that z decreases
    [z_left,p] = sort(MESH.GCOORD(2,inod_left),'descend');
    x_left     = x_left(p);

    % Form polygon around domain
    xz_bnd        = [x_bot x_rght x_top x_left;
                     z_bot z_rght z_top z_left];
    xz_bnd        = unique(xz_bnd','rows','stable')';
    xz_bnd(:,end) = xz_bnd(:,1);
    
    if FigNo
        figure(FigNo);clf
        plot(xz_bnd(1,:),xz_bnd(2,:),'k.-');
        hold on
    end
    
%     [iin,ion] = inpolygon(gX_PT(1,ind),gX_PT(2,ind),x_bnd,z_bnd);
%     iin       = ind(iin);
%     iout(iin) = 0;
end

% Deal with points that are truly outside the mesh
ind   = find(iout);
if ~isempty(ind) && nargout>1
    ztol = 0.01;
    xtol = 0.01;
    
    for it=1:4
        % CHECK TOP BOUNDARY
        zbnd_top = interp1(x_top,z_top,gX_PT(1,ind),'linear','extrap');
        iout     = gX_PT(2,ind)>zbnd_top-ztol;
        if any(isnan(zbnd_top))
            1;
        end
        if FigNo; scatter(gX_PT(1,ind(iout)),gX_PT(2,ind(iout)),10,'r'); end
        gX_PT(2,ind(iout)) = zbnd_top(iout) - ztol;
        if FigNo; scatter(gX_PT(1,ind(iout)),gX_PT(2,ind(iout)),10,'g'); end


        % CHECK BOTTOM BOUNDARY
        zbnd_bot = interp1(x_bot,z_bot,gX_PT(1,ind),'linear','extrap');
        iout     = gX_PT(2,ind)<zbnd_bot+ztol;
        if any(isnan(zbnd_bot))
            1;
        end
        if FigNo; scatter(gX_PT(1,ind(iout)),gX_PT(2,ind(iout)),10,'r'); end
        gX_PT(2,ind(iout)) = zbnd_bot(iout) + ztol;
        if FigNo; scatter(gX_PT(1,ind(iout)),gX_PT(2,ind(iout)),10,'g'); end


        % CHECK LEFT BOUNDARY
        xbnd_left = interp1(z_left,x_left,gX_PT(2,ind),'linear','extrap');
        iout      = gX_PT(1,ind)<xbnd_left+xtol;
        if any(isnan(xbnd_left))
            1;
        end
        if FigNo; scatter(gX_PT(1,ind(iout)),gX_PT(2,ind(iout)),10,'r'); end
        gX_PT(1,ind(iout)) = xbnd_left(iout) + xtol;
        if FigNo; scatter(gX_PT(1,ind(iout)),gX_PT(2,ind(iout)),10,'g'); end


        % CHECK RIGHT BOUNDARY
        xbnd_rght = interp1(z_rght,x_rght,gX_PT(2,ind),'linear','extrap');
        iout      = gX_PT(1,ind)>xbnd_rght-xtol;
        if any(isnan(xbnd_rght))
            1;
        end
        if FigNo; scatter(gX_PT(1,ind(iout)),gX_PT(2,ind(iout)),10,'r'); end
        gX_PT(1,ind(iout)) = xbnd_rght(iout) - xtol;
        if FigNo; scatter(gX_PT(1,ind(iout)),gX_PT(2,ind(iout)),10,'g'); end


        % CHECK THAT CORRECTION IS OK
        iin = inpolygon(gX_PT(1,ind),gX_PT(2,ind),xz_bnd(1,:),xz_bnd(2,:));
        if any(~iin)
            if FigNo
                scatter(gX_PT(1,ind(~iin)),gX_PT(2,ind(~iin)),50,'r');
            end
        else
            break
        end
    end
end

if verbose
    nout = length(find(iout));
    fprintf(' %1i (%.2f%%) points outside of domain',...
        nout,100*nout/nPT);
    if nargout==2
        fprintf('; shifted onto boundary.\n');
    else
        fprintf('\n');
    end
end

end