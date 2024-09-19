function La = resampleLine(Lb, dmax, method)
    % La = resampleLine(Lb, s, method) 
    % +++ purpose +++
    % add nodes to Lb so that the distance between succesive nodes does not exceed a given maximum spacing 
    % Note: the function differs from resample_line() in that here, existing nodes are preserve 
    % 
    % Lb     :: [2,nb] input line
    % dmax   :: maximum spacing [m]
    % method :: argument passed to matlab interp1() function
    %
    % Author: Javier Garcia-Pintado, MARUM, 2021

    if nargin < 3
        method = "makima";
    end

    nb = size(Lb,2);
    nseg = nb - 1;

    dchain = diff2D(Lb);                                       % space between nodes
    nnn = floor(abs(dchain)/dmax);                             % [1,nseg] number of new nodes for each segment
    
    idsb = 1:nb;                                               % [1,nb]
    idsbaug = repelem(idsb, [nnn+1 1]);                        % [1,na] expanded indices
    idsbnn = find(nnn ~= 0);                                   % node indices in Lb previous to insertion of new nodes
    chainb = [0 cumsum(dchain)]; % chainage
    chaina = chainb(idsbaug);    
    
    nk = numel(idsbnn);                                        % number of segments to be auugmented 
    for k=1:nk                                                 % fill block of new coordinates
        iseg = idsbnn(k); % == inode
        chnew = linspace(chainb(iseg),chainb(iseg+1),nnn(iseg)+2);
        idsa  = find(idsbaug==iseg);
        chaina(:,idsa(2:end)) = chnew(2:end-1);
    end
    La = [interp1(chainb, Lb(1,:), chaina, method);
          interp1(chainb, Lb(2,:), chaina, method)];
    
    % examples - not run
    if 1 > 2 
        % example 1
        x = rand(1,8);
        y = rand(1,8);
        Lb = [x;y];
        dmax = 0.1;
        La = resampleLine(Lb, dmax);
       
        figure(); plot(Lb(1,:),Lb(2,:),'.-')
        hold on;  plot(La(1,:),La(2,:),'o-')
        
        % example 2
        dmax = 0.2;
        x = linspace(-pi,3*pi,10);
        y = sin(x);
        Lb = (rotation([x' y'], 80*pi/180))';
        La = resampleLine(Lb, dmax);
        Ls = resample_line(Lb,dmax,"makima"); 
        figure(); plot(Lb(1,:),Lb(2,:),'.-')
        hold on;  plot(La(1,:),La(2,:),'o')
        hold on;  plot(Ls(1,:),Ls(2,:),'s')
    end
end % function
