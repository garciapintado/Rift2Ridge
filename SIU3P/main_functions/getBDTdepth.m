function [BDTdepth, BDTat] = getBDTdepth(ErP, E2all, Temp, GCOORD, ELEM2NODE, GEOn, eids, PHY) 
    % estimate location and depth of the brittle-ductile transition (BDT) from the sea bottom
    %
    % Details:
    % a] define "isplastic" as the domain where ErP > ErV
    % b] locate coordinate of the maximum ErP
    % c] define "inspreading" as an horizontal 20km-wide domain around the coordinates of ErP_max
    % d] obtain the maximum ErV within both domains above [the "subdomain", hereafter] 
    % e] get a normalized ErV field ("obj1" \in [0,1]) [which equals 0.0 out of the subdomain]
    % f] locate deepest coordinate of this subdomain
    % g] get a normalized vertical-distance field ("obj2" \in [0,1]) [which equals 0.0 out of the subdomain] 
    % h]
    %
    % Javier Garcia-Pintado, 2021

    sphwin = 15e03; % [m] half-window over which the BDT is evaluated around the maximum ErP
    
    ErV = E2all - ErP;
    
    [nel,nip] = size(ErP);    
    [GIPx,GIPy] = ip_coord(GCOORD, ELEM2NODE, nel, nip);                         % both in [nel,nip]
    isplastic = (2*ErP > E2all);
    
    alt = 1;
    switch alt
        case 1 % slower but less jumpy
            ErPn = ipval2nodval(ELEM2NODE(1:6,:), ErP, eids, true); 
            ErPxint = approxXmarginalIntegrate2D(GCOORD, ELEM2NODE(1:6,:),ErPn, GEOn);     % figure(); plot(ErPxint.seq/1000,ErPxint.int) 
            [~,maxid] =  max(ErPxint.int);
            ErPxcent = ErPxint.seq(maxid);
            %E2n = ipval2nodval(ELEM2NODE(1:6,:), E2all, eids, true);
            %Ernxint = approxXmarginalIntegrate2D(GCOORD, ELEM2NODE(1:6,:),E2n, GEOn);      % hold on;  plot(Ernxint.seq/1000,Ernxint.int, 'color','red')
        case 2
            [rid1,cid1] = find(ErP == max(ErP,[],'all'));                                  % figure(); plot_ip_valF(ErP, GCOORD, ELEM2NODE, GEOn(9).coo, [], "Plastic strain");
            ErPxcent = GIPx(rid1,cid1);                                                  % hold on; plot(GIPx(rid1,cid1)/1000, GIPy(rid1,cid1)/1000,'s','markersize',10)
    end % switch   
    inspreading = abs(GIPx - ErPxcent) < sphwin;                            % xline((GIPxErPmax - sphwin)/1000); xline((GIPxErPmax + sphwin)/1000)
    % obj1: ErV within the brittle subdomain
    obj1 = ErV / max(ErV(isplastic & inspreading)) .* isplastic .* inspreading;    % normalized maximum [max(obj1,[],'all')==1] | figure(); plot_ip_valF(obj1, GCOORD, ELEM2NODE, GEOn(9).coo, [], "obj1");

    % obj2: ErP within the plastic subdomain
    obj2 = ErP / max(ErP(~isplastic & inspreading)) .* ~isplastic .* inspreading;  % figure(); plot_ip_valF(obj4, GCOORD, ELEM2NODE, GEOn(9).coo, [], "obj4"); 
    
    % obj3: horizontal distance from spreading center
    obj3 = exp(-0.5*((ErPxcent - GIPx)/sphwin).^2) .* inspreading;
       
    % obj4: the deeper the better within the britle "suddomain"
    GIPyminPlastic = min(GIPy(isplastic & inspreading));
    obj4 = abs(GIPy - (GIPyminPlastic)) .* isplastic .* inspreading;      
    obj4 = (1 - obj4/max(obj4(:))) .* isplastic .* inspreading;                    % figure(); plot_ip_valF(obj2, GCOORD, ELEM2NODE, GEOn(9).coo, [], "obj2");
    
    % obj5: Temperature [BDT likelihood decaying away from the range (300,600)]
    TempQP = nodval2ipval(ELEM2NODE(1:6,:), Temp, 6);                             % [nel,nip]
    obj5 = 1 - variomodels(max(0.0,abs(TempQP-PHY.BDT_Tcent)-PHY.BDT_Trang), 0, 1, PHY.BDT_Trang, 'gau') .* inspreading; % figure(); plot(TempQP,obj5)

    %[rid4,cid4] = find(obj4 == max(obj4,[],'all'));
    %BDTp = [GIPx(rid4,cid4), GIPy(rid4,cid4)];                                     % hold on; plot(BDTp(1)/1000, BDTp(2)/1000,'o','markersize',10,'color',[.0 .5 .5])                     % 

    w = 1/5;
    obj = w * (obj1 + obj2 + obj3 + obj4 + obj5);                                   % figure(); plot_ip_valF(obj, GCOORD, ELEM2NODE, GEOn(9).coo, [], "obj");
    [rid,cid] = find(obj == max(obj,[],'all'));
    
    BDTat = [GIPx(rid,cid), GIPy(rid,cid)];
    BDTdepth = interp1(GEOn(end-1).coo(1,:),GEOn(end-1).coo(2,:),BDTat(1)) - BDTat(2);
    
    % figure(); plot_ip_valF(ErV, GCOORD, ELEM2NODE, GEOn(9).coo, [], "Viscous strain");
    % figure(); plot_ip_valF(isplastic*1.0,    GCOORD, ELEM2NODE, GEOn(9).coo, [],          "isplastic");
    % hold on; plot(BDTat(1)/1000, BDTat(2)/1000,'o','markersize',10,'color',[.5 .0 .5])
end
