function plot_p_axes2(GCOORD, ELEM2NODE, F_xx, F_xy, F_yx, F_yy, r, sratr, edgecolors)
    % Raghu's plot for deviatoric stresses
    nip         = 6;
    nnodel1      = 6;
    [IP_X, IP_w] = ip_triangle(nip);
    [Nbig]       = shp_triangle(IP_X, nnodel1);
    for i=1:nel
        is         = (i-1)*nnodel1+1; ie = (i-1)*nnodel1+nnodel1;
        GCOORD_N(:,is:ie) = GCOORD(:,ELEM2NODE(1:nnodel1,i));
        EL2N(i,:) = is:ie;
        Dummy      = Nbig'\TAU_xx(i,:)';
        TAU_xxn(:,i)= Dummy(1:nnodel1);
        Dummy      = Nbig'\TAU_yy(i,:)';
        TAU_yyn(:,i)= Dummy(1:nnodel1);
        Dummy      = Nbig'\TAU_xy(i,:)';
        TAU_xyn(:,i)= Dummy(1:nnodel1);
    end
    % Calculate principal stresses
    TAU_1 = ((TAU_xxn + TAU_yyn)./2) + sqrt(((TAU_xxn - TAU_yyn)./2).^2 + TAU_xyn.^2); %%Sigma_1
    TAU_3 = ((TAU_xxn + TAU_yyn)./2) - sqrt(((TAU_xxn - TAU_yyn)./2).^2 + TAU_xyn.^2); %% Singma_3
    angn = rad2deg(atan((2.*TAU_xyn)./(TAU_xxn - TAU_yyn)))/2.;
    TAU_p1 = (((TAU_xxn + TAU_yyn))/2.) + ((((TAU_xxn - TAU_yyn))/2.).*cosd(2* angn)) + TAU_xyn.* sind(2*angn);
    thetap = angn;
    thetap((abs(TAU_1(:,:) -  TAU_p1(:,:)) > 10)) = thetap(((abs(TAU_1(:,:) -  TAU_p1(:,:)) > 10))) +90; %%Correct for the angle thetap
    TAU_p2 = (((TAU_xxn + TAU_yyn))/2.) + ((((TAU_xxn - TAU_yyn))/2.).*cosd(2* thetap)) + TAU_xyn.* sind(2*thetap); %Sigma_1 = TAU_1 above
    TAU_p3 = (((TAU_xxn + TAU_yyn))/2.) - ((((TAU_xxn - TAU_yyn))/2.).*cosd(2* (thetap))) - TAU_xyn.* sind(2*(thetap)); %Sigma_3

    figure(1);
    eids = "645"; % rift2ridge2D
    thetapn = ipval2nodval(ELEM2NODE, thetap',eids,false);
      plot_tF(ipval2nodval(ELEM2NODE, thetap',eids,true), GCOORD, ELEM2NODE);
    colorbar

    hold on;
    Xcenter = GCOORD(1,ELEM2NODE(7,:));
    Ycenter = GCOORD(2,ELEM2NODE(7,:));
    TAU_p2ip = ipval2nodval(ELEM2NODE, TAU_p2',eids,false);
    TAU_p3ip = ipval2nodval(ELEM2NODE, TAU_p3',eids,false);
    
    %T =(TAU_p2ip(7,:)/1e6)./max(TAU_p2ip(7,:)/1e6);
    Mag = sqrt((TAU_p2ip(:,:).^2 + TAU_p3ip(:,:).^2 ));
    F_1= TAU_p2ip(7,:)'./Mag(7,:)';
    F_2= TAU_p2ip(7,:)'./Mag(7,:)';
    F_3=TAU_p3ip(7,:)'./Mag(7,:)';
    F_4=TAU_p3ip(7,:)'./Mag(7,:)';
    
    %%Sigma - 1  %%to plot sigma_1 umcomment next two lines
    %q = quiver(Xcenter/1e3,Ycenter/1e3,F_1'.*cosd(thetapn(7,:)),F_2'.*sind(thetapn(7,:)),'-b');
    %q.ShowArrowHead = 'off';
    %hold on;
    %%Sigma - 3 %%to plot sigma_3 umcomment next two lines
    %q1 = quiver(Xcenter/1e3,Ycenter/1e3,F_3'.*cosd(90+thetapn(7,:)),F_4'.*sind(90+thetapn(7,:)),'-k');
    
    %q1.ShowArrowHead = 'off';
end
