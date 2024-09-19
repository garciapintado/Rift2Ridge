function GEO = linkGEO(GEO)
    % generate GEO.to:: internal GEO(i) -> GEO(j) pointers based on
    % geometric coordinates
    for i=1:length(GEO)
        GEO(i).n = size(GEO(i).coo,2);                                     % in case this was not updated
        GEO(i).to = zeros(2,GEO(i).n);
        for k=1:length(GEO)
            if k <= i                                                      % only higher order interfaces as possible parents
                continue;
            end
            
            [iink_x,kids_x] = ismember(GEO(i).coo(1,:), GEO(k).coo(1,:));
            [iink_y,kids_y] = ismember(GEO(i).coo(2,:), GEO(k).coo(2,:));
            iink = iink_x & iink_y & (kids_x == kids_y);
            if sum(iink) == 0
                continue;
            end
            kids = uint32(kids_x);
            kids(~iink) = 0;
            GEO(i).to(1,iink) = kids(iink);
            GEO(i).to(2,iink) = k;
        end  
    end
end % end function linkGEO()
