function TRACKP = trackPoints(GCOORD, ELEM2NODE, GEOn, ...
    TRACKP, Vel, dt)
    
    % trackPoints(GCOORD, ELEM2NODE, TRACKP, Vel, dt)
    % +++purpose +++
    % Calculates new positions for tracked points TRACKP after deformation of the current
    % time step DT, using the velocities Vel, defined in the nodes of a mesh
    % GCOORD, ELEM2NODE.
    %
    % IO
    % TRACKP :: [2,ntrk]
    %
    % Javier Garcia-Pintado: 2020-03 - supersedes various tracking point functions by MA Martinez
    %                                  merged and optimised code [python option for tsearch() and parallel outer location search] 
    %
    % TRACKP = TRACKP(:,1:50000);
    
    TRACKP(1,:) = max(TRACKP(1,:),min(GCOORD(1,:)));                       % truncate coordinates to the left of domain
    TRACKP(1,:) = min(TRACKP(1,:),max(GCOORD(1,:)));                       % truncate coordinates to the right of domain
    
    TRACKP_ybot = interp1(GEOn(1).coo(1,:),GEOn(1).coo(2,:),...            % [2,nzeros] : bottom level at outer TRACK x coordinates 
                        TRACKP(1,:));                                      % plot(Surface_sort(1,:),Surface_sort(2,:),'o-')
    TRACKP_ytop = interp1(GEOn(end-1).coo(1,:),GEOn(end-1).coo(2,:),...    % [2,nzeros] : surface level at outer TRACK x coordinates 
                        TRACKP(1,:));                                      % plot(Surface_sort(1,:),Surface_sort(2,:),'o-')
    isbelow = TRACKP_ybot - TRACKP(2,:) >= -1.0E-03;                       % LOGICAL [1,ntrack]
    isabove = TRACKP(2,:) - TRACKP_ytop >= -1.0E-03;                       % LOGICAL [1,ntrack]
    xy = TRACKP(1:2,:);
    xy(2,isbelow) = TRACKP_ybot(isbelow) + 1.0E-03;
    xy(2,isabove) = TRACKP_ytop(isabove) - 1.0E-03;                        % figure(); plotGEO(GEOn); hold on; plot(xy(1,:)/1000, xy(2,:)/1000,'x','color','red')

    tic;
    xyel = tsearch3(GCOORD, uint32(ELEM2NODE(1:3,:)), xy);
    disp(['trackPoints:: tsearch3() computing time: ',num2str(toc)])

    iserr = xyel==0;                                                       % [1,ncoo]
    if any(iserr)
       error("trackPoints:: points with unlocated parent element") 
    end

    Vx_tp = remesh_val(xyel, GCOORD, xy, Vel(1,:), ELEM2NODE);             % Interpolate velocities at the trackpoints
    Vy_tp = remesh_val(xyel, GCOORD, xy, Vel(2,:), ELEM2NODE);
    TRACKP(1:2,:) = TRACKP(1:2,:) + [Vx_tp; Vy_tp] .* dt;                  % update track point positions
    
end % function