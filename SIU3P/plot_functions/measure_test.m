%% INPUT
directory_p = '/Volumes/GeoPhysik/miguel/data/SIA3/';
file_name = 'main';
ID = 1121:1223;

%% LOOP THROUGH TESTS
for n = 1:length(ID)
    % Find name in main
    newFileData = fileread([file_name,num2str(ID(n)),'.m']);
    ii = strfind(newFileData,'directory = ');
    jj = strfind(newFileData,';'); jj  = jj(jj>ii); jj = jj(1);
    row_s = newFileData(ii:jj);
    k = strfind(row_s,'/'); k = k(end)+1;
    l = jj-ii-2;
    
    name_t = row_s(k:l);
    
    % Load variables
    dir_p = [directory_p,'/',name_t,'/'];
    step = num2str(lastest(dir_p));
    load([dir_p,'_',step])
    
    % Plot
    h = figure(1);
    set(h,'Visible','off')
    
    clf
    MESH.EL2NOD = ELEM2NODE;
    MESH.GCOORD = GCOORD/1000;
    plot_val(MESH,I2.p,size(ELEM2NODE,2),6)
    hold on
    plot_box
    plot(ISOCHRONS(1,ISOCHRONS(3,:)==0)/1000,ISOCHRONS(2,ISOCHRONS(3,:)==0)/1000,'--')
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
    
    disp(['TEST: ',name_t])
    if surf_proc_s
        plot_isoc
        disp(max(-ISOCHRONS(2,ISOCHRONS(3,:)==0)/1000+Topography(2,:)/1000))
    end
    caxis([0 5])
    tstep
    measure
    input('Next')
end