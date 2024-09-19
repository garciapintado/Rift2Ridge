function xo = findZeroes(x, y)
    signs = sign(y(:));
    ix_zeros   = find(signs == 0);                                         % exact crossings
    ix_crosses = find(abs(diff(signs)) == 2);                              % mid-segment crossings
    ids = [ix_zeros; ix_crosses];
    xo = zeros(length(ids),1);
    for k=1:length(ids)
       i = ids(k);
       xo(k) = interp1(y([i i+1]), x([i i+1]),0);
    end                                                                    % plot(x,y); yline(0,'color',[.3 .3 .3])
                                                                           % xline(xo)
    if 1 > 2 % example call                                   
       xo = findZeroes(Topography(1,:), Topography(2,:) - sea_level); 
       figure(); plot(Topography(1,:), Topography(2,:))
       yline(sea_level);
       hold on; scatter(xo,repelem(sea_level,1,length(xo)),'o');
    end
    
end