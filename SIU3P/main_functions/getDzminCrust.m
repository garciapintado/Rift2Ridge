function dzmincrust = getDzminCrust(GEO)
    % +++ purpose +++
    % give an unprecise but fast global estimate of minimum crust thickness
    % this is intended as approximate measure to initiate hdrothermal circulation 
    % There is no need to be more precise as global criterion for this purpose

    hids = find([GEO.horizontal]);
    dzmincrust = min(GEO(hids(end)).coo(2,:)) -  max(GEO(hids(2)).coo(2,:));

end
