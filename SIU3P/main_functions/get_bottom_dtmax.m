function dtmax = get_bottom_dtmax(GCOORD,ELEM2NODE, DISPL, dt, ybot)
% dtmax = get_bottom_dtmax(GCOORD, ELEM2NODE, DISPL, dt, ybot) 
% 
% +++ purpose +++
% get maximum allowed dtmax so that for a given velocity (DISPL) 
% no triangular element goes fully below ybot

%--------------------------------------------------------------------------
% Author: Javier García-Pintado 2020-01-21
%--------------------------------------------------------------------------
    
    GCOORD_a = GCOORD + DISPL*dt;                                                 % GCOORD updated by mechanical solution - done later in main
    nod_belowbot = GCOORD_a(2,:) <= ybot;
    ele_belowbot = sum(reshape(nod_belowbot(ELEM2NODE(1:3,:)),3,[])) == 3;         % elements which have been moved completely below Ylim(1)
    if any(ele_belowbot)
        [nodmaxa,nolid] = max(reshape(GCOORD_a(2,ELEM2NODE(1:3,ele_belowbot)),3,[]));  % nolid: in [1,3] is local vertex index for the top vertices in ele_belowbot
        nogid = ELEM2NODE(sub2ind(size(ELEM2NODE), nolid, find(ele_belowbot)));        % size(nogid) == sum(ele_belowbot)
        dtmax = 4./5. * min((GCOORD(2,nogid) - ybot) ./ ( -DISPL(2,nogid)));        % maximum allowed dt, 4/5 is a security factor to allow for not too flat triangles
    else
        dtmax = dt;
    end
    % foo_b = GCOORD(:,ELEM2NODE(1:3,ele_belowbot));
    % foo_a = GCOORD_a(:,ELEM2NODE(1:3,ele_belowbot));
    % figure(); plot_meshF(ELEM2NODE, GCOORDa); plot_geomF(GEOMETRY, Geo_id)
    % hold on; plot(foo_b(1,:)/1000,foo_b(2,:)/1000,'o','color',[.0 .0 .2]);
    %          plot(foo_a(1,:)/1000,foo_a(2,:)/1000,'o','color','red');
    %          plot(GCOORD(1,nogid)/1000,GCOORD(2,nogid)/1000,'x','color',[.0 .0 .2])
    %          plot(GCOORD_a(1,nogid)/1000,GCOORD_a(2,nogid)/1000,'x','color','red')
    
end % function    
