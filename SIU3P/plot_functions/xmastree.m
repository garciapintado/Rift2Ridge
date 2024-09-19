function [STR] = xmastree(file,E,xloc)
% XMASTREE(FILE,E,XLOC) loads a step of a model given by the directory FILE
% and plots the differential stresses (black) for a determined strain rate 
% E and the actual differential stresses calculated from the stresses given 
% by mechanical solver (blue), for a depth profile specified by the user 
% using a GUI graph showing the actual differential stresses of the model.
% The horizontal location can also be selected by specifying it into XLOC.
% If E is leaved empty [] the function will take the boundary condition 
% velocities and calculate the strain rate. Press a click outside of the
% temperature plot to finish the script and plot the differential stress 
% profile in a new figure.

%--------------------------------------------------------------------------
% Function written by Miguel Andres-Martinez, PhD student at Royal Holloway
% University of London, 20-11-2014. Email: mandresmartinez87@gmail.com
% Modified by Elena Ros a little bit :P
%--------------------------------------------------------------------------

% 21/04/2015 MA
    % FIXED a bug where the x-coordinate read from the manual input was
    % not multiplied by 1000 and, therefore, the interpolation of
    % temperatures and phases was not calculated in the correct position
    % FIXED calculated pressures. Now pressures are calculated from the
    % weight of the whole rock column (including topography)
    % ADDED plot_box to the temperature field
    
% 29/05/2015 ER
    % Added strength (TAU_II) field plot and actual strength profile
    
% 17/06/2015 MA
    % Improved interpolations of the strength profiles by using shape
    % functions
    
% 20/06/2018 MA
    % Improvements in brittle and viscous strength calculations. For
    % brittle, cos(Phi) is now multiplying the cohesion and lithostatic
    % pressures are calculated by integrating densities in the column and
    % not by using a constant density of 2700 kg/m^3. For viscous changes
    % are: 1) factors for transforming from principal stresses to second
    % invariant formulation are now included in the stress calculation, 2)
    % multiple rheologies per layer are now included in the calculation,
    % and 3) activation volume is now included in the flow law    
    
% 05/09/2018 MA
    % Added the possibility to specify the horizontal location by
    % specifying the variable XLOC

%==========================================================================
% SET UP
%==========================================================================
figure(11)
load(file)

res_p = 1; % [Km]
mpa = 1e6; % MPa
km = 1000; % Km

% Strain calculation
if isempty(E)
    switch bc_v
        case 'ext_rate'
            E = ext_rate/max(GCOORD(1,:));
        case 'ext_erate'
            E = ext_erate;
    end
end

%==========================================================================
% PLOT TEMPERATURES
%==========================================================================
% Map values of the ips into the nodes
    nip1         = 6;
    nnodel1      = 6;
    [IP_X, IP_w] = ip_triangle(nip1);
    [Nbig]       = shp_triangle(IP_X, nnodel1);
for i=1:nel
    is         = (i-1)*nnodel1+1; ie = (i-1)*nnodel1+nnodel1;
    GCOORD_N(:,is:ie) = GCOORD(:,ELEM2NODE(1:nnodel1,i));
    EL2N(i,:) = is:ie;
    T_n(is:ie) = Temp(ELEM2NODE(1:nnodel1,i));    
    Dummy      = Nbig'\TAU_xx(i,:)';
    TAU_xxn(:,i)= Dummy(1:nnodel1);
    Dummy      = Nbig'\TAU_yy(i,:)';
    TAU_yyn(:,i)= Dummy(1:nnodel1);
    Dummy      = Nbig'\TAU_xy(i,:)';
    TAU_xyn(:,i)= Dummy(1:nnodel1);
    Dummy      = Nbig'\RHO(i,:)';
    Rho_n(:,i)= Dummy(1:nnodel1);
    for n = 1:length(RHEOL.var)
        Dummy  = Nbig'\RHEOL.var{n}(i,:)';
        Var_n{n}(:,i) = Dummy(1:nnodel1);
    end
end
% Calculate principal stresses
% TAU_1 = ((TAU_xxn + TAU_yyn)./2) + sqrt(((TAU_xxn - TAU_yyn)./2).^2 + TAU_xyn.^2);
% TAU_3 = ((TAU_xxn + TAU_yyn)./2) - sqrt(((TAU_xxn - TAU_yyn)./2).^2 + TAU_xyn.^2);
% Strength = (TAU_1 - TAU_3)/mpa; % in MPa
% Strength = sqrt((1/2*(TAU_xxn-TAU_yyn)).^2+TAU_xyn.^2)/mpa;
Strength = sqrt(1/2*(TAU_xxn.^2+TAU_yyn.^2)+TAU_xyn.^2)/mpa;
  
% % Plot temperatures
%subplot(1,2,1)
% patch('faces',EL2N,'vertices',GCOORD_N'/1000,'facevertexcdata',T_n(:),'FaceColor','flat')
% shading interp
% axis tight
% colormap(jet)
% colorbar
% hold on
% plot_box

% % IP
% [   N, dNdu]    = shp_deriv_triangle(IP_X, size(ELEM2NODE,1));
% ECOORD_x = reshape(GCOORD(1,ELEM2NODE), size(ELEM2NODE,1), nel);
% ECOORD_y = reshape(GCOORD(2,ELEM2NODE), size(ELEM2NODE,1), nel);
% for ip=1:nip
%     Ni      =        N{ip};
%     GIP_x   = Ni'*ECOORD_x;
%     GIP_y   = Ni'*ECOORD_y;
%     GIP_xF(:,ip)   = GIP_x;
%     GIP_yF(:,ip)   = GIP_y;
% end
% T_1 = ((TAU_xx + TAU_yy)./2) + sqrt(((TAU_xx - TAU_yy)./2).^2 + TAU_xy.^2); 
% T_3 = ((TAU_xx + TAU_yy)./2) - sqrt(((TAU_xx - TAU_yy)./2).^2 + TAU_xy.^2);
% ST = (T_1-T_3)/mpa;

% Plot principal stresses
subplot(1,2,1)
patch('faces',EL2N(:,1:3),'vertices',GCOORD_N'/1000,'facevertexcdata',Strength(:),'FaceColor','flat')
shading interp
%axis tight
axis([min(GCOORD_N(1,:)/1000) max(GCOORD_N(1,:)/1000)  -120 max(GCOORD_N(2,:)/1000)] )
colormap(jet)
%colormap(parula)
colorbar
hold on
plot_box

% hold on
% indx_plot = GIP_xF >= 90000 & GIP_xF <= 120000 & GIP_yF >= -60000 & GIP_yF <= -30000;
% scatter(GIP_xF(indx_plot)/1000,GIP_yF(indx_plot)/1000,100,'w','fill')
% scatter(GIP_xF(indx_plot)/1000,GIP_yF(indx_plot)/1000,50,ST(indx_plot),'fill')

%%xtick = roundn(min(GCOORD_N(1,:)/1000),1):20:roundn(max(GCOORD_N(1,:)/1000),1);
%xticklabels = {'1y', '10y', '20y', '30y', '40y', '50y', '55y', '60y'};
%%set(gca, 'XTick', xtick);
%%rotateXLabels(gca, 90)
%set(gca, 'XTickLabel', xticklabels);
%%ytick = roundn(min(GCOORD_N(2,:)/1000),1):10: roundn(max(GCOORD_N(2,:)/1000),1);
%yticklabels = {'Jan 2005', 'Jun 2005', 'Jan 2006', 'Jun 2006', ...etc };
%%set(gca, 'YTick', ytick);
%set(gca, 'YTickLabel', yticklabels);

title(['Strength profile ',num2str((istep*dt/ma)),' Ma'])
xlabel('Distance [km]')
ylabel('Depth [km]')
hold on

%==========================================================================
% GUI AND LOCATION OF THE PROFILE
%==========================================================================
% Initial envelope
% ----------------
% Vertical line
Y = min(GCOORD(2,:)):res_p:max(GCOORD(2,:));
% Initial X
X = min(GCOORD(1,:))*ones(size(Y));
% Plot
plot(X([1 end])/1000,Y([1 end])/1000,'w')
% Get line for later remove
last_line = get(gca, 'children');

if nargin==3
    xs = true;
else
    xs = false;
end

keep_plot = 1;
% Plotting loop
while keep_plot
    if xs
        keep_plot = 0;
        x = xloc/1000;
        y = (max(GCOORD(2,:))+min(GCOORD(2,:)))/2000;
    else
        % Take pointer input
        [x,y] = ginput(1);
    end
    % If inside the area for plotting plot new envelope
    if (x*1000 >= min(GCOORD(1,:)) && x*1000 <= max(GCOORD(1,:))) && ...
            (y*1000 >= min(GCOORD(2,:)) && y*1000 <= max(GCOORD(2,:)))
        subplot(1,2,1)
        hold on
        % Delete previous line in the subplot 1
        delete(last_line(1));
        % Recalculate vertical line (otherwise it may untake the points
        % from the topography of previous lines starting deeper)
        Y = min(GCOORD(2,:)):res_p:max(GCOORD(2,:));
        % Calculate new vector of X
        X = x*1000*ones(size(Y));
        % Plot new line in subplot 1
        plot(X([1 end])/1000,Y([1 end])/1000,'w')
        % Mark last line for later remove
        last_line = get(gca, 'children');
        hold off
        % Remove points above topography
        topo = interp1(GCOORD(1,Point_id==(max(Point_id)-1)), ...
                GCOORD(2,Point_id==(max(Point_id)-1)),x*1000);
        X(Y>topo) = [];
        Y(Y>topo) = [];
            
        % Calculating envelope (from principal stresses)
        %Strength_line = griddata(GCOORD_N(1,:)',GCOORD_N(2,:)',Strength(:),X,Y);
        Tris = tsearch2(GCOORD,uint32(ELEM2NODE(1:3,:)),[X; Y]);

        Ind = find(Tris==0);
        if(~isempty(Ind))
            for i=1:length(Ind)
                [~, Tris(Ind(i))] = min(sqrt((GCOORD(1,ELEM2NODE(7,:)) - ...
                    X(Ind(i))).^2 + (GCOORD(2,ELEM2NODE(7,:)) - ...
                    Y(Ind(i))).^2));
            end
        end
        
        if(any(isnan(Tris)))
            error('remeshing failed in move_contours');
        end
        
        xx = x;
        yy = y;
        
        % Calculates the local coordinates of the new ip
        xp = X; % New ip
        yp = Y; % New ip
        
        x = reshape(GCOORD(1,ELEM2NODE(1:3,Tris)),3,size(xp,2)); % Reshape old coord
        y = reshape(GCOORD(2,ELEM2NODE(1:3,Tris)),3,size(yp,2)); % Reshape old coord
        
        xi = -(-x(1,:).*yp+x(1,:).*y(3,:)-x(3,:).*y(1,:)+xp.*y(1,:)+x(3,:) ...
            .*yp-xp.*y(3,:))./(-x(1,:).*y(3,:)+x(1,:).*y(2,:)-x(2,:).*y(1,:) ...
            +x(2,:).*y(3,:)+x(3,:).*y(1,:)-x(3,:).*y(2,:)); % Local x cood for ip
        yi = (x(1,:).*y(2,:)-x(1,:).*yp-x(2,:).*y(1,:)+x(2,:).*yp+xp ...
            .*y(1,:)-xp.*y(2,:))./(-x(1,:).*y(3,:)+x(1,:).*y(2,:)-x(2,:).*y(1,:)...
            +x(2,:).*y(3,:)+x(3,:).*y(1,:)-x(3,:).*y(2,:)); % Local y coord for ip
        
        % Load shape functions
        
        [NN] = shp_triangle([xi' yi'], nip);
        
        % Calculate the value of F for the new ips
        
        Strength_line = sum(NN.*Strength(1:6,Tris));

        % Calculating envelope (from temperatures)
        Temp_el = Temp(ELEM2NODE(1:6,:));
        T = sum(NN.*Temp_el(1:6,Tris));
%         subplot(1,2,2)
%         plot(T,Y/km,'k')
        x = xx;
        y = yy;

%==========================================================================
% CALCULATE PHASES
%==========================================================================
        % Find interfaces ids
        Interf_id = [1 3:3:max(Point_id)-1];
        
        % Find points where interfaces intersect with the profile
        Interfaces = zeros(1,length(Interf_id));
        for n = 1:length(Interf_id)
            Interfaces(n) = interp1(GCOORD(1,Point_id==Interf_id(n)), ...
                GCOORD(2,Point_id==Interf_id(n)),x*1000);
        end
        
        % Assign phases to Y
        Profile_pha = zeros(length(Y),1);
        for n = 1:length(Interf_id)-1
            Profile_pha(Y>=Interfaces(n) & Y<=Interfaces(n+1)) = n;
        end
%==========================================================================
% BRITTLE
%==========================================================================
        % Calculate lithostatic pressures at the profile
        RHOi = sum(NN.*Rho_n(1:6,Tris));
        Plit = cumsum(RHOi*G(2)*res_p);
        Plit = Plit-Plit(1);
        Plit = -Plit(end:-1:1);
        
        % Brittle shear stress (Mohr-Coulomb)
        T_b = sin(SS.Phi(1))*Plit + cos(SS.Phi(1))*SS.C(1);
        T_bw = sin(SS.Phi(2))*Plit + cos(SS.Phi(2))*SS.C(1);
%         % Plot
%         subplot(1,2,2)
%         hold on
%         plot(T_b/mpa,Y/km,'r')

%==========================================================================
% CREEP BEHAVIOUR
%==========================================================================
% Hirth and Kohlstedt (2003)

% Initialise variables for loop
T_dis = zeros(size(Y));
T_dif = zeros(size(Y));

% Loop for rheologic variation in the same phase
for n = 1:size(RHEOL.var,2)
    % Get rheologic factors for multy-rheologies layers mapped into the profile
    Var = sum(NN.*Var_n{n}(1:6,Tris));
    
    if exist('RHEOL','var')
        % Factors for scaling triaxial and uniaxial experiment
        % parameters (GERYA 2010)
        Sc_dis = ...
            1./(2.^((RHEOL.Ndis(Profile_pha,n)-1)./RHEOL.Ndis(Profile_pha,n)).* ...
            3.^((RHEOL.Ndis(Profile_pha,n)+1)./(2*RHEOL.Ndis(Profile_pha,n))));
        Sc_dif = 1/3;
        
        % Dislocation shear stress
        T_dis = T_dis + Var.*Sc_dis'.* ...
            (E./RHEOL.Adis(Profile_pha,n)').^ ...
            (1./RHEOL.Ndis(Profile_pha,n)').* ...
            exp((RHEOL.Qdis(Profile_pha,n)' + ...
            Plit.*RHEOL.Vdis(Profile_pha,n)')./ ...
            (R*(T+273).*RHEOL.Ndis(Profile_pha,n)'));
        
        % Diffusion shear stress
        T_dif = T_dif + Var.*Sc_dif'.* ...
            (E./RHEOL.Adif(Profile_pha,n)').^ ...
            (1./RHEOL.Ndif(Profile_pha,n)').* ...
            exp((RHEOL.Qdif(Profile_pha,n)' + ...
            Plit.*RHEOL.Vdif(Profile_pha,n)')./ ...
            (R*(T+273).*RHEOL.Ndif(Profile_pha,n)'));
    else
        % Factors for scaling triaxial and uniaxial experiment
        % parameters (GERYA 2010)
        Sc_dis = ...
            1./(2.^((Ndis(Profile_pha,n)-1)./Ndis(Profile_pha,n)).* ...
            3.^((Ndis(Profile_pha,n)+1)./(2*Ndis(Profile_pha,n))));
        Sc_dif = 1/3;
        
        % Dislocation shear stress
        T_dis = T_dis + Var.*Sc_dis'.* ...
            (E./Adis(Profile_pha,n)').^ ...
            (1./Ndis(Profile_pha,n)').* ...
            exp((Qdis(Profile_pha,n)' + ...
            Plit.*Vdis(Profile_pha,n)')./ ...
            (R*(T+273).*Ndis(Profile_pha,n)'));
        
        % Diffusion shear stress
        T_dif = T_dif + Var.*Sc_dif'.* ...
            (E./Adif(Profile_pha,n)').^ ...
            (1./Ndif(Profile_pha,n)').* ...
            exp((Qdif(Profile_pha,n)' + ...
            Plit.*Vdif(Profile_pha,n)')./ ...
            (R*(T+273).*Ndif(Profile_pha,n)'));
    end
end
%         % Plot
%         subplot(1,2,2)
%         hold on
%         plot(T_dis/mpa,Y/km,'.r')
        
%         % Plot
%         subplot(1,2,2)
%         hold on
%         plot(T_dif/mpa,Y/km,'.g')

%==========================================================================
% PLOT
%==========================================================================
        % Minimum shear stress
        T_min = min([T_b; T_dis; T_dif]);
        T_minw = min([T_bw; T_dis; T_dif]);
        % Differential stress
        
        disp(['Force: ',num2str(sum(Strength_line*1e3)),' N'])
        
        % Plot
        subplot(1,2,2)
        plot(Strength_line,Y/km,'b','LineWidth',2)
        hold on
        plot(T_min/mpa,Y/km,'k')
        plot(T_minw/mpa,Y/km,'r')
        hold off
        %area(T_min/mpa,Y/1000)
        
%%        xtick =  roundn(min(Strength_line),1):10: roundn(max(Strength_line),1);
%         %xticklabels = {'Jan 2005', 'Jun 2005', 'Jan 2006', 'Jun 2006', ...etc };
%%        set(gca, 'XTick', xtick);
%%        rotateXLabels(gca, 90)
        
%%        ytick =  roundn(min(Y/1000),1):10: roundn(max(Y/1000),1);
        %yticklabels = {'Jan 2005', 'Jun 2005', 'Jan 2006', 'Jun 2006', ...etc };
%%        set(gca, 'YTick', ytick);
        %set(gca, 'YTickLabel', yticklabels)
        %axis([roundn(min(T_min/mpa),1)  roundn(max(T_min/mpa),1) roundn(min(Y/1000),1) roundn(max(Y/1000),1)])
        %axis([roundn(min(T_min/mpa),1)   roundn(max(T_min/mpa),1) round(max(max(GCOORD_N(2,:)/1000))) roundn(max(Y/1000),1)])
        %axis([roundn(min(Strength_line),1) roundn(max(Strength_line),1) roundn(min(Y/1000),1) max(max(GCOORD_N(2,:)/1000))])
%%        axis([roundn(min(Strength_line),1) roundn(max(Strength_line),1)   -120   max(max(GCOORD_N(2,:)/1000))])
        
        set(gca,'XMinorGrid','on')
        title('Differential stress profile')
        xlabel('Diff. stress [MPa]')
        ylabel('Depth [Km]')

%==========================================================================
% END
%==========================================================================
        
        % If outside the area of plotting, finish plotting
    else
        hold off
        keep_plot = 0;
    end
end

figure(12)
subplot(1,2,1)
plot(T_min/mpa,Y/km,'k','linewidth',2)
title('Differential stress profile')
xlabel('Diff. stress [MPa]')
ylabel('Depth [km]')
axis([0 max(T_min/mpa)+50 -100 max(Y/km)])
daspect([10 1 2])

STR.yield   = T_min/mpa;
STR.yield_s = T_minw/mpa;
STR.y       = Y/km;
STR.E       = E;
STR.T       = T;
STR.str     = Strength_line;

subplot(1,2,2)
plot(Strength_line,Y/km,'k','linewidth',2)
title('Differential stress profile')
xlabel('Diff. stress [MPa]')
ylabel('Depth [km]')
axis([0 max(T_min/mpa)+50 -100 max(Y/km)])
daspect([10 1 2])