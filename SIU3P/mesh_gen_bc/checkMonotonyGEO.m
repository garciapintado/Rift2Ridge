function monotonic = checkMonotonyGEO(GEO, tolerance)
  % +++ purpose +++
  % check that x coordinates in quasi-horizontal GEO interfaces increase
  % monotonically
  
  % Javier Garcia_Pintado, MAURM, 2020
  
  monotonic = true;
  
  for i=1:length(GEO)
    if GEO(i).horizontal
        l2r = diff(GEO(i).coo(1,[1,end])) > 0. ;
        if ~l2r
            error("checkMonotonyGEO: non left to right GEO horizontal interface at GEO("+i+")")
        end
        diffx = diff(GEO(i).coo(1,:));
        backids = find(diffx <= 0.);
        if any(backids)
            monotonic = false;                                                                                      % figure(); plotGEO(GEO);
            disp("WARNING: checkMonotonyGEO: x-non monotonic horizontal interfaces at GEO("+i+").coo(1,["+num2str(backids)+"])") % hold on; plot(GEO(i).coo(1,backids)/1000, GEO(i).coo(2,backids)/1000,'x','markersize',15) 
            if i == max(find([GEO.horizontal])) || any(diffx < -tolerance)
                error("x non-monotonicity allowance exceeded")
            end
        end
        clear backids;
    end
  end
end % function

