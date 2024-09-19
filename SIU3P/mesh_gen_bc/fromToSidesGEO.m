function GEO = fromToSidesGEO(GEO)
    % add GEO(i).gidbounds for sides interfaces, identifying the
    % neighbouring interfaces bounding GEO(i)
    
    imax = length(GEO);
    rightids = 2:3:imax; % domain right side interfaces
    leftids = 4:3:imax;  % domain left side interfaces
    
    for i=1:length(GEO)
        GEO(i).class = "layer";                                      % in preparation for fault interfaces class in {"lside","rside","layer","fault"} 
        if ismember(GEO(i).gid, leftids)
            GEO(i).class = "lside";
            GEO(i).boundgids = [max(i-4,1) i-1];
            GEO(i).boundcoo = [GEO(max(i-4,1)).coo(:,1) GEO(i-1).coo(:,1)];
        end
        if ismember(GEO(i).gid, rightids)
            GEO(i).class = "rside";
            GEO(i).boundgids = [max(i-2,1) i+1];
            GEO(i).boundcoo = [GEO(max(i-2,1)).coo(:,end) GEO(i+1).coo(:,end)];
        end
    end   
end % end function fromToSidesGEO
