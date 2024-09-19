function [ISOCHRONS, NBASE] = resample_isoc(GCOORD, Point_id, ELEM2NODE, ...
    ISOCHRONS, Basement, tp_isoc, res)
% [ISOCHRONS, NBASE] = RESAMPLE_ISOC(GCOORD,POINT_ID,ELEM2NODE,ISOCHRONS,
% BASEMENT,TP_ISOC,RES) spatial horizontal resample of ISOCHRONS add points to ISOCHRONS where
% points are separated by a distance larger than RES and it also resamples
% the BASEMENT in the same way.

%--------------------------------------------------------------------------
% Function written by Yanfang Xin, PhD student at University of
% Bremen, 04-03-2018. 
% Simplified rewritting by Javier García-Pintado, MARUM, 2020. 
%--------------------------------------------------------------------------

    % Calculate topography
    [Topography,~] = find_topo(GCOORD, ELEM2NODE, Point_id);
    % Length of the topography
    ltopo = size(Topography,2);

    ISOC = [];
    isoids = unique(ISOCHRONS(3,:));
    for i=1:length(isoids)
    ISOC(i).iso = ISOCHRONS(:,ISOCHRONS(3,:)==isoids(i));
    end
    
    ics = [1 2 4]; % assumes only 4 columns and that the third one is 'istep'
    asNA = -999;
    % Time lines resampling
    % ---------------------
    if tp_isoc % isochrons are tracked
        ISOCHRONS = [];
        for i = 1:length(ISOC)
            it = ISOC(i).iso(3,1);
            CISOC = ISOC(i).iso;                                          % [:,ncoo], where ncoo is number of coordinates in this isochron
            ncoo = size(CISOC,2);
            dist_p = sqrt(diff(CISOC(1,:)).^2+diff(CISOC(2,:)).^2);       % distance between points
            larseg = dist_p > res;                                        % LOGICAL [1,ncoo] too large segments
            NISOC = repelem(asNA, size(CISOC,1), size(CISOC,2) + sum(larseg)); % initialise new isochron
            add_np = cumsum([0 larseg]);                                 % [1,ncoo]
            NISOC(:,(1:ncoo)+add_np) = CISOC;                            % fill in existing data
            is0 = find(all(NISOC == asNA));
            NISOC(ics,is0) = (NISOC(ics,is0-1) + NISOC(ics,is0+1)) / 2;
            NISOC(3,is0) = it;
            ISOCHRONS = [ISOCHRONS NISOC];
        end
    end % if tp_isoc

    % Basement
    % --------
    ics = [1 2];
    CISOC = Basement;
    ncoo = size(CISOC,2);
    dist_p  = sqrt(diff(CISOC(1,:)).^2+diff(CISOC(2,:)).^2);
    larseg  = dist_p > res;
    NBASE = repelem(asNA,size(CISOC,1),size(CISOC,2) + sum(larseg)); % initialise new isochron
    add_np = cumsum([0 larseg]);                              % [1,ncoo]
    NBASE(:,(1:ncoo)+add_np) = CISOC; 
    is0 = find(all(NBASE == asNA));
    NBASE(ics,is0) = (NBASE(ics,is0-1) + NBASE(ics,is0+1)) / 2;
end % function