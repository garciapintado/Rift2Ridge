% TESTS CONSORTIUM________________________________________________________________

load('/Users/elena/WORK/TESTS_SS_2016/TESTS_CONSORTIUM/FIGURES-TESTS_CONSORTIUM_noHR/MG35_3.mat')


% load('/Users/elena/WORK/TESTS_SS_2016/TESTS_CONSORTIUM/FIGURES-TESTS_CONSORTIUM_noHR/MG40_3.mat')


% load('/Users/elena/WORK/TESTS_SS_2016/TESTS_CONSORTIUM/FIGURES-TESTS_CONSORTIUM_noHR/WQ35_3.mat')

 
% load('/Users/elena/WORK/TESTS_SS_2016/TESTS_CONSORTIUM/FIGURES-TESTS_CONSORTIUM_noHR/WQ40_3.mat')folder=1;
%function  plot_series_serp_melt_underplating_crust(folder,directory1)
%directory1 = '/Users/elena/WORK/TESTS_DIKE/MG_WQ_v3_v5_M35_M40_with_dike_study2_8_models/TESTS_DIKE/RESULTS_/1_v3_M35_MG_dike';

%% Main

load([directory1,'/','_',num2str(1),'.mat'])

for ii = 1:10:length(save_step)%length(files)-3  %3861:10:7001
    
    load([directory1,'/','_',num2str(save_step(ii)),'.mat'])
    
    d = figure(3); clf
    %d = figure('Visible','off');clf;
    %% MAPPING FROM INTEGRATION POINTS TO NODES - IP CANNOT BE PLOTTED
    
    EL2N      = zeros(nel,3); % new connectivity matrix
    GCOORD_N  = zeros(2,nel*3);
    TAU_xxn   = zeros(3,nel);
    TAU_yyn   = zeros(3,nel);
    Pres_n    = zeros(3,nel);
    Mu_n      = zeros(3,nel);
    E2_n      = zeros(3,nel);
    T_n      = zeros(3,nel);
    Dpl_n      = zeros(3,nel);
    Dserp_n   = zeros(3,nel);
    dF_n   = zeros(3,nel);
    dFserp_n = zeros(3,nel);
    
    EL2N(1,:) = 1:3;
    nip = 6;
    nnodel = 6;
    [IP_X, IP_w]    = ip_triangle(nip);
    [   Nbig]    = shp_triangle(IP_X, nnodel);
    for i=1:nel
        is         = (i-1)*3+1; ie = (i-1)*3+3;
        GCOORD_N(:,is:ie) = GCOORD(:,ELEM2NODE([1 2 3],i));
        EL2N(i,:) = is:ie;
        T_n(is:ie) = Temp(ELEM2NODE([1 2 3],i));
        Dpl_n(is:ie) = Dpl(ELEM2NODE([1 2 3],i));
        Dummy      = Nbig'\TAU_xx(i,:)';
        TAU_xxn(:,i)= Dummy(1:3);
        Dummy      = Nbig'\TAU_yy(i,:)';
        TAU_yyn(:,i)= Dummy(1:3);
        Dummy      = Nbig'\Mu_all(i,:)';
        Mu_n(:,i)= Dummy(1:3);
        Dummy      = Nbig'\PRES_IP(i,:)';
        Pres_n(:,i)= Dummy(1:3);
        Dummy      = Nbig'\E2all(i,:)';
        E2_n(:,i)= Dummy(1:3);
        Dserp_n(is:ie) = Dserp(ELEM2NODE([1 2 3],i));
        dF_n(is:ie) = dF(ELEM2NODE([1 2 3],i));
        dFserp_n(is:ie) = dFserp(ELEM2NODE([1 2 3],i));
    end
    
    % Calculates the density field
    DENSITY = Rho(Phases)';
    DENSITY = repmat(DENSITY,3,1);
    
    %% SERPENTINIZATION
    
    subplot(3,1,1)
    %title(['Serpentinization ',num2str((istep*dt/ma)),' Ma'])
    title(['Accumulated serpentinization ',num2str((istep*dt/ma)),' Ma'])
    xlabel('Distance [km]')
    ylabel('Depth [km]')
    %patch('faces',EL2N,'vertices',GCOORD_N'/1000,'facevertexcdata',dFserp_n(:),'FaceColor','flat')
    patch('faces',EL2N,'vertices',GCOORD_N'/1000,'facevertexcdata',Dserp_n(:),'FaceColor','flat')
    xmin=-200;
    xmax= 200;
    ymin= -40;
    ymax= 2.5;
    boxplotx= [xmin xmin xmax xmax xmin ];
    boxploty = [ymax ymin ymin ymax ymax] ;
    axis ([xmin xmax ymin ymax])
    
    shading interp
    colorbar
    hold on
    %plot(Boxx,Boxy,'k')
    hold on;
    hc = colorbar('location','eastoutside');
    xlabel(hc,'%Serp')
    hold on;
    cmax = max(max(dFserp_n));
    display(cmax)
    
    caxis([0 1])
    %set(gca, 'visible', 'off')
    %set(findall(gca, 'type', 'text'), 'visible', 'on')
    hold on;
    %% Plots interface
    % Plot box and interfaces
    % -----------------------
    plot_box
    drawnow
    hold off
    hold on
    drawnow
    
    
    hold on
%       %% UNDERPLATING
%     hold on;
%     %if exist('TP_xmelt')~=0
%     if tp_melt==1
%         if Crust_thickness(istep)~=0
%             %FOR DIKE_SILL_______________________________
%             k = find(ismember(TRACKP_melt(3,:),istep-1),1,'last'); %the same as max(find(ismember(TRACKP_melt(3,:),istep-1))); and more efficient 
%             %all the columns where there is answer 1 are from the previous time step, so I have to plot till there, till the maximum.
%             plot(TRACKP_melt(1,1:k)/1000,TRACKP_melt(2,1:k)/1000,'.','Color',[0, 0.5, 0],'MarkerSize',10)
%             %____________________________________________
%             drawnow
%         end
%     end
%     hold on;
%     if(Crust_thickness(istep)~=0)
%         hold on
%         %FOR DIKE____________________________________
%         %                fill(xdike/1000,ydike/1000,'m')
%         %                hold on;
%         %                plot(xdike/1000,ydike/1000,'.-m')
%         %FOR DIKE_SILL_______________________________
%         fill(x_ign_body/1000,y_ign_body/1000,'m')
%         hold on;
%         plot(x_ign_body/1000,y_ign_body/1000,'.-m')
%         %____________________________________________
%     end
    hold on;
    plot(boxplotx,boxploty, '-r')
    hold on;
%     if tp_melt==1
%         if Crust_thickness(istep)~=0
%             for k=1:istep-1
%                 plot(TRACKP_melt{k}(1,:)/1000,TRACKP_melt{k}(2,:)/1000,'.','Color',[0, 0.5, 0],'MarkerSize',7)
%             end
%             %____________________________________________
%             drawnow
%         end
%     end
%     hold on;
%     if(Crust_thickness(istep)~=0)
%         hold on
%         fill(x_ign_body/1000,y_ign_body/1000,'m')
%         hold on;
%         plot(x_ign_body/1000,y_ign_body/1000,'.-m')
%     end
%     hold on;
    %% CONTOURS TEMPERATURE
    
    % x_=GCOORD_N(1,:)';
    % y_=GCOORD_N(2,:)';
    % z_=T_n(:);
    % xi=linspace(min(x_),max(x_),100);
    % yi=linspace(min(y_),max(y_),100);
    % [XI, YI]=meshgrid(xi,yi);
    % ZI = griddata(x_,y_,z_,XI,YI);
    
    %figure(12);
    %cs=contour(XI/1000,YI/1000,ZI, 'LineWidth',2,'Color','k')
    % [cs,h]=contour(XI/1000,YI/1000,ZI,[200,350], 'LineWidth',1,'Color','r');
    % hold on
    % cc=clabel(cs,h,'Color','r');
    % drawnow
    hold on;
    
    %% MELTING
    subplot(3,1,2)
    patch('faces',EL2N,'vertices',GCOORD_N'/1000,'facevertexcdata',Dpl_n(:),'FaceColor','flat')
    shading interp
    axis([-250 250 -400 5])
    colorbar
    hold on
    
    title(['Area of melt production ',num2str((istep*dt/ma)),' Ma'])
    hold on;
    hc = colorbar('location','eastoutside');
    xlabel(hc,'%Melt')
    hold on;
    cmin =0 ;
    cmax = max(max(Dpl_n));
    display(cmax)
    %caxis([0 1e-3])
    caxis([0 0.1])
    hold on;
    
    plot_box
    axis([-250 250 -400 5])
    drawnow
    hold off
    hold on
    drawnow
    
    
         %% UNDERPLATING
    hold on;
    %if exist('TP_xmelt')~=0
    if tp_melt==1
        if Crust_thickness(istep)~=0
            %FOR DIKE_SILL_______________________________
            k = find(ismember(TRACKP_melt(3,:),istep-1),1,'last'); %the same as max(find(ismember(TRACKP_melt(3,:),istep-1))); and more efficient 
            %all the columns where there is answer 1 are from the previous time step, so I have to plot till there, till the maximum.
            plot(TRACKP_melt(1,1:k)/1000,TRACKP_melt(2,1:k)/1000,'.','Color',[0, 0.5, 0],'MarkerSize',10)
            %____________________________________________
            drawnow
        end
    end
    hold on;
    if(Crust_thickness(istep)~=0)
        hold on
        %FOR DIKE____________________________________
        %                fill(xdike/1000,ydike/1000,'m')
        %                hold on;
        %                plot(xdike/1000,ydike/1000,'.-m')
        %FOR DIKE_SILL_______________________________
        fill(x_ign_body/1000,y_ign_body/1000,'m')
        hold on;
        plot(x_ign_body/1000,y_ign_body/1000,'.-m')
        %____________________________________________
    end
    hold on;
    plot(boxplotx,boxploty, '-r')
    hold on;
    %% CRUSTAL THICKNESS
    subplot(3,1,3)
    time =[dt/ma:dt/ma:dt*istep/ma];
    plot(time(1:istep),Crust_thickness(1:istep)/1000,'o','MarkerSize',3) %,'linewidth',30)
    xlabel('Time [Ma]','FontSize',9)
    ylabel('Crustal thickness [km]','FontSize',9)
    title(['Magmatic crustal thickness ',num2str((istep*dt/ma)),' Ma'])
    axis([0 45 -0.1 7])
    grid on
    %set(gca, 'visible', 'off')
    drawnow;
    %% SAVE IMAGE
    %cd ([num2str(folder),'_study_SH_serp_melt_underplating_location_fixed'])
    %print(d,'-dpng','-r150',[num2str(folder) '_dike_serp_melt_crust_' num2str(istep*dt/ma) '.png']);
    %cd ..
    display(pwd)
    display(num2str(istep*dt/ma))
    %movefile([num2str(folder) '_dike_serp_melt_crust_' num2str(istep*dt/ma) '.png'],'\Users\elena\WORK\TESTS_DIKE\TEST_DIKE_SHEAR_HEATING\Plots\6_study2\'[num2str(folder) '_dike_serp_melt_crust_' num2str(istep*dt/ma) '.png']')
    
    
    %print(d,'-dpng','-r150',['..\',num2str(j_folders(ij)),'_study2','\',num2str(folder) '_dike_serp_melt_crust_' num2str(istep*dt/ma) '.png']);
    %webFigureHandle = webfigure(d);
    
    
    
    %   print(d,'-dpng','-r150',['TEST_craton_MG40_Tserp_noHR_serp_melt' num2str(istep*dt/ma) '.png']);
    %   print(d,'-dpng','-r150',['TEST_craton_MG40_Tserp_HR_serp_melt' num2str(istep*dt/ma) '.png']);
    %   print(d,'-dpng','-r150',['TEST_craton_MG40_noTserp_HR_serp_melt' num2str(istep*dt/ma) '.png']);
    
    
    %   print(d,'-dpng','-r150',['TEST_foldbelt_MG35_Tserp_noHR_serp_melt' num2str(istep*dt/ma) '.png']);
    %   print(d,'-dpng','-r150',['TEST_foldbelt_MG35_Tserp_HR_serp_melt' num2str(istep*dt/ma) '.png']);
    %   print(d,'-dpng','-r150',['TEST_foldbelt_MG40_Tserp_noHR_serp_melt' num2str(istep*dt/ma) '.png']);
    %   print(d,'-dpng','-r150',['TEST_foldbelt_MG40_Tserp_HR_serp_melt' num2str(istep*dt/ma) '.png']);
    %
    dir_o = '/Users/marta/WORK/PROGS_RIFT/DYN_PROGS/CRATON_NOCRATON_8-2015/PROGS_AUGUST_2015_MIGUEL_VERS_3PHASES_TEST4/FIGURES/WQ40_Tserp_NOHR_LC20_1300_ss_LAB120_CONTpress_newCe_remeshold_3lyr_DensDpl_RH';
    
    mkdir(dir_o);
    print(d,'-dpng','-r150',[dir_o,'/', 'MELT_SERP', num2str(istep*dt/ma) '.png']);
    %saveas(f,[dir_o,'/','MELT_SERP',num2str((istep*dt/ma)),'.png'],'png')
    
    %   print(d,'-dpng','-r150',['TEST_foldbelt_WQ35_Tserp_HR_serp_melt' num2str(istep*dt/ma) '.png']);
    %   print(d,'-dpng','-r150',['TEST_foldbelt_WQ40_Tserp_noHR_serp_melt' num2str(istep*dt/ma) '.png']);
    %   print(d,'-dpng','-r150',['TEST_foldbelt_WQ40_Tserp_HR_serp_melt' num2str(istep*dt/ma) '.png']);
    
    %close(d);
end
%end