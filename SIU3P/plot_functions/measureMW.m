%% MEASURE MARGIN WIDTH SCRIPT
% This script loops through a number of mains named as file_name followed
% by a number and finding name of saved files, then it loads the file and
% measure the margin width from the highest part of the shoulder, to the
% break of the crust.

clear;

%==========================================================================
%% INPUT
%==========================================================================
% Directory of saved files
directory_p     = '/Volumes/GeoPhysik/miguel/data/SIA3/';
% Common name of the main files
file_name       = 'main';
% Main file numbers to loop through
%ID              = 1110:1147; % GR
ID              = 1148:1185; % WQ
%ID              = 1186:1223; % AN

% Geometric parameters
% --------------------
% Thickness of the crust to be considered breaked [m]
thick_bk    = 250;
% Minimum distance between shoulders [m]
min_d_shou  = 30000;

%==========================================================================
%% INITIALISATION OF VARIABLES
%==========================================================================
h = figure(1);
MW = zeros(20,length(ID));
% Variables to load
VAR = {'GCOORD','ELEM2NODE','Point_id','Corner_id','Cornin_id','GEOMETRY','nip','nel', ...
    ...
    'surf_proc_s','tp_isoc','dt','ma','istep', ...
    ...
    'I2','ISOCHRONS','SP','Basement'};

%==========================================================================
%% LOOP THROUGH TESTS
%==========================================================================
for nnp = 1:length(ID)
    % Load file
    % ---------
    % Find name of the saved file in main
    newFileData = fileread([file_name,num2str(ID(nnp)),'.m']);
    ii = strfind(newFileData,'directory = ');
    jj = strfind(newFileData,';'); jj  = jj(jj>ii); jj = jj(1);
    row_s = newFileData(ii:jj);
    k = strfind(row_s,'/'); k = k(end)+1;
    l = jj-ii-2;
    name_t = row_s(k:l);
    
    % Load variables
    dir_p = [directory_p,'/',name_t,'/'];
    step = num2str(lastest(dir_p));
    try
        load([dir_p,'_',step],VAR{:})
        file_e = 1;
    catch ME
        file_e = 0;
    end
    
    % Plot
    % ----
    clf
    if file_e==1
        MESH.EL2NOD = ELEM2NODE;
        MESH.GCOORD = GCOORD/1000;
        plot_val(MESH,I2.p,size(ELEM2NODE,2),6)
        hold on
        plot_box
        plot(ISOCHRONS(1,ISOCHRONS(3,:)==0)/1000, ...
            ISOCHRONS(2,ISOCHRONS(3,:)==0)/1000,'--')
        plot_sea
        
        [Topography,~] = find_topo(GCOORD,ELEM2NODE,Point_id);
        moho = GCOORD(:,Point_id==3);
        mohoS = interp1(moho(1,:),moho(2,:),Topography(1,:));
        Cthick = Topography(2,:)-mohoS;
        % 30 km
        uu = Cthick<30*1000;
        C30 = find([0 uu(1:end-2)+uu(2:end-1)+uu(3:end) 0]==2);
        for nn = 1:length(C30)
            plot([Topography(1,C30(nn)) Topography(1,C30(nn))]/1000, ...
                [Topography(2,C30(nn)) mohoS(1,C30(nn))]/1000,'r')
        end
        % 35 km
        moho = GCOORD(:,Point_id==3);
        mohoS = interp1(moho(1,:),moho(2,:),Topography(1,:));
        Cthick = Topography(2,:)-mohoS;
        % 30 km
        uu = Cthick<35*1000;
        C35 = find([0 uu(1:end-2)+uu(2:end-1)+uu(3:end) 0]==2);
        for nn = 1:length(C35)
            plot([Topography(1,C35(nn)) Topography(1,C35(nn))]/1000, ...
                [Topography(2,C35(nn)) mohoS(1,C35(nn))]/1000,'g')
        end
        
        if exist('axp','var')
            axis(axp)
        end
        
        disp(['TEST: ',name_t])
        if surf_proc_s
            plot_isoc
            % Calculate maximum sediment thickness
            ST = -ISOCHRONS(2,ISOCHRONS(3,:)==0)+Topography(2,:);
            [~,sST] = sort(ST,'descend');
            maxST = sST(1);
            mST = ST(maxST);
            cmST = [ISOCHRONS(1:2,maxST) Topography(1:2,maxST)];
            plot(cmST(1,:)/1000,cmST(2,:)/1000,'r','LineWidth',2)
            disp(['Max. sed. thickness:  ', ...
                num2str(mST/1000),' km       Breakup:  ', ...
                num2str(dt*istep/ma),' ma'])
        else
            mST = 0;
            cmST = zeros(2,2);
            disp(['Max. sed. thickness:  ', ...
                'NO SED.       Breakup:  ', ...
                num2str(dt*istep/ma),' ma'])
        end
        title([name_t,'      ',num2str(dt*istep/ma),' ma'],'Interpreter', ...
            'none')
        caxis([0 5])
        
        % Prepare geometric variables
        % ---------------------------
        % Save topography or basal estimated topography
        if surf_proc_s
            Topo = ISOCHRONS(2,ISOCHRONS(3,:)==0);
            Base = [min([ISOCy(1:end,:); Topography(2,:)])];
            if SP.ksea==0
                SP.kdecay = 0;
            end
        else
            Topo = Topography(2,:);
            Base = Topography(2,:);
            SP.c = 0;
            SP.ksea = 0;
            SP.kdecay = 0;
        end
        
        rerun = 1;
        while rerun
            
            % Identify shoulders
            % ------------------
            % Short topography to find maximums
            [~,inds] = sort(Topo,'descend');
            % Find first maximum that would be consider as shoulder or area
            % of highest uplift
            SHOULDERS = [Topography(1,inds(1)); Topo(inds(1))];
            % Find second maximum that would be consider as shoulder or 
            % area of highest uplift considering the minimum distance 
            % between shoulders (min_d_shou)
            shoui2 = inds(find(abs(Topography(1,inds(2:end))-SHOULDERS(1,:)) ...
                >=min_d_shou,1,'first')+1);
            SHOULDERS = [SHOULDERS [Topography(1,shoui2); Topo(shoui2)]];
            [SHOULDERS(1,:),indS] = sort(SHOULDERS(1,:));
            SHOULDERS(2,:) = SHOULDERS(2,indS);
            % Plot shoulders for QC
            pointsS = plot(SHOULDERS(1,:)/1000,SHOULDERS(2,:)/1000,'.', ...
                'MarkerSize',50);
            
            % Identify breakup points
            % -----------------------
            % Crustal thickness in for top nodes
            Point_id(Cornin_id([1 2])) = 3;
            Moho = GCOORD(:,Point_id==3);
            [Moho(1,:),iM] = sort(Moho(1,:));
            Moho(2,:) = Moho(2,iM);
            MohoTop = interp1(Moho(1,:),Moho(2,:),Topography(1,:));
            cthick = Base-MohoTop;
            % Finding breakup points
            indB = [find(cthick<=thick_bk,1,'first') ...
                find(cthick<=thick_bk,1,'last')];
            BREAKP = [Topography(1,indB); Base(indB)];
            % Plot shoulders for QC
            pointsB = plot(BREAKP(1,:)/1000,BREAKP(2,:)/1000,'.', ...
                'MarkerSize',50);
            
            % Modify parameters during the run using m option
            drawnow
            inp = input('(1) save or (2) modify: ');
            if isempty(inp)
                inp = 1;
            end
            switch inp
                case 2
                    keyboard
                    delete(pointsS)
                    delete(pointsB)
                case 1
                    rerun = 0;
            end
        end
    end
    
    % Store information
    % -----------------
    % MW LEFT   MW RIGHT   MAX SED THICK   
    if file_e==1
        if ~isempty(BREAKP)
            MW(1:2,nnp)     = abs(SHOULDERS(1,:)-BREAKP(1,:));
            MW(3,nnp)       = mST;
            MW(4,nnp)       = dt*istep/ma;
            MW(5,nnp)       = SP.c;
            MW(6,nnp)       = SP.ksea;
            MW(7,nnp)       = SP.kdecay;
            MW(9:12,nnp)    = BREAKP(:);
            MW(13:16,nnp)   = SHOULDERS(:);
            MW(17:20,nnp)   = cmST(:);
        else
            MW(:,nnp)   = 0;
        end
    else
        MW(:,nnp)   = 0;
        disp(['TEST: ',name_t])
        disp(['NO FILE FOUND!'])
        inp = input('(1) continue or (2) modify: ');
        if isempty(inp)
            inp = 1;
        end
        switch inp
            case 2
                keyboard
            case 1
                rerun = 0;
        end
    end
    
    % Get current axis
    uuu = gca;
    axp = [uuu.XLim uuu.YLim];
    
end