function [ocrust, ocrust_melt, ocrust_serp, isocrustip] = getOceanicCrust(GCOORD, ELEM2NODE, GEO, MELTING, Dserp, PHY, ext_rate, ntimes)
    % [ocrust, ocrust_melt, ocrust_serp, isocrustip] = getOceanicCrust(GCOORD, ELEM2NODE, GEO, MELTING, Dserp, PHY, ext_rate, ntimes)
    %
    % gets an estimation of the oceanic crust under the following assumptions:
    % a) crust thickness at each x-location is approximated by the integral along the
    % vertical of the amount serpentinized volume plus the amount of allocated melt
    % b) the generated melt plus the generated total serpentinization
    % should equal the total estimated volume of oceanic crust (area in 2D)
    
    % INPUT
    % MELTING.
    %        .Crust_thickness: REAL [1×5000 double]
    %        .Trackp: [4×108162 double]
    % Dserp   :: REAL [1,nnod] degree of serpentinization in [0,1]
    % PHY     :: struct with parameters to estimate the oceanic crustal thickness
    %    .serp_dx   [m] 
    %    .melt_dx   [m]
    %
    % OUTPUT
    % ocrust      :: REAL [2,:], [m] polyline representing the bottom of the oceanic crust
    % isocrustip  :: LOGICAL [nel,nip] true for integration points within the oceanic crust 
    % 
    % ext_rate & ntimes: only for plotting purposes
    %
    % Javier GP, MARUM, 2021

    nip = 6;
    nel = size(ELEM2NODE,2);
    
    moho = GEO(3).coo;                                                     % [2,nmoho]                  
    ocrust = moho;                                                         % total oceanic crust
    ocrust_melt = moho;                                                    % only melting product component
    ocrust_serp = moho;                                                    % only serpentinization component

    % estimate serpentinization thickness                                  a] prediction
    if any(Dserp > 0.0)
        dx = PHY.serp_dx;                                                      % [m] approximate sampling spacing [discretization] of the horizontal estimation for serpentinization 
        incx = diff(minmax(GCOORD(1,:)));                                      % [m] horizontal size of domain
        nx =  ceil(incx/dx);                                                   % number of horizontal discretization nodes                                      
        dserpx = approxXmarginalIntegrate2D(GCOORD, ELEM2NODE, Dserp, GEO, nx); % figure(); plot(dserpx.seq/1000, dserpx.int/1000)                                                              
        dserpxm = max(0,interp1(dserpx.seq,dserpx.int,moho(1,:),'makima'));    % [1,nmoho] serpentinization depth at x-moho
                                                                               % b] correction: make conservative
        serp_int  = sum(diff(moho(1,:)) .* (dserpxm(2:end) + dserpxm(1:end-1)) / 2);      
        serp_in2D = integrate2D(GCOORD, ELEM2NODE, Dserp);                     % figure(); plot_tF(Dserp, GCOORD, ELEM2NODE); 
        if serp_int ~= 0.0
            dserpxm = dserpxm * serp_in2D/serp_int;                            % hold on; plot(moho(1,:)/1000,moho(2,:)/1000,'color','white','linewidth',1)
        end                                                                    % hold on; plot(moho(1,:)/1000,(moho(2,:) - dserpxm)/1000,'color',[.5 .2 .0],'linewidth',2)
        ocrust_serp(2,:) = ocrust_serp(2,:) - dserpxm;
        ocrust(2,:) = ocrust(2,:) - dserpxm;
    end                                                                       % hold on; scatter(MELTING.Trackp(1,:)/1000,MELTING.Trackp(2,:)/1000,5,'green')
    
    % estimate thicknes of magmatic bodies                                 % hold on; plot_dikesF(MELTING, ext_rate, ntimes); 
    if max(MELTING.Crust_thickness) > 0.0
        nxlags = 10;                                                           % 1 for no lag
        xlags = 0:1/nxlags:1;
        xlags = xlags(1:end-1); % normalised xlag
        xlags = xlags * PHY.melt_dx;
        
        trk_area = (pi*MELTING.Trackp(3,:).^2)';                               % [m2]; [ntrk,1] area represented by each tracking point
        melt_in2D = sum(trk_area);                                             % [m2] integral of melting areas 
        meltxms = [];
        for im=1:nxlags                                                        % number of regionalized samples
            dx = PHY.melt_dx;                                                      % [m] horizontal discretization block
            x = unique([moho(1,1) (xlags(im)+moho(1,1)):dx:moho(1,end) moho(1,end)]);
            nx = length(x);
            blkint = zeros(1,nx-1);                                                % [m]; [nx-1] init block heights
            meltxid  = interp1(x,1:nx,MELTING.Trackp(1,:),'previous')';            %       [ntrk,1] discretization block for each tracking point [i.e. block where each tracking point lies]
            meltxint = accumarray(meltxid,trk_area);
            blkint(unique(meltxid)) = meltxint(unique(meltxid));                   % horizontal block-wise integrals
            dx = diff(x);
            meltxm = max(0,interp1((x(1:end-1)+x(2:end))/2,blkint./dx,moho(1,:),'makima')); % hold on; plot(moho(1,:)/1000,(moho(2,:) - meltxm)/1000,'color','red','linewidth',2)
            melt_int  = sum(diff(moho(1,:)) .* (meltxm(2:end) + meltxm(1:end-1)) / 2);
            meltxms(im).dz = meltxm * melt_in2D/melt_int;                     % b] correction: make conservative
        end
        meltxmbar = mean(reshape([meltxms.dz],[],nxlags),2)';                  % [1,nmoho]
        melt_int  = sum(diff(moho(1,:)) .* (meltxmbar(2:end) + meltxmbar(1:end-1)) / 2);
        if melt_int ~= 0
            meltxmbar = meltxmbar * melt_in2D/melt_int; 
        end
        ocrust_melt(2,:) = ocrust_melt(2,:) - meltxmbar;                       % hold on; plot(ocrust(1,:)/1000,ocrust(2,:)/1000,'color','brown','linewidth',2)
        ocrust(2,:) = ocrust(2,:) - meltxmbar;
    end

    isocrustip = false(nel,nip);                                           % [nel,nip]
    [gipx,gipy] = ip_coord(GCOORD, ELEM2NODE, nel, nip);                   % [nel,nip]
    isocrustip(gipy > interp1(moho(1,:),ocrust(2,:),gipx) & gipy < interp1(moho(1,:),moho(2,:),gipx)) = true;  % scatter(ipx(isocrustip)/1000,ipy(isocrustip)/1000,5,'cyan')
      
  end % function
  