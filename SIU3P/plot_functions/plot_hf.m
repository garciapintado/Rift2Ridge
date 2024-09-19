% VERTICAL HEAT FLOW PLOT
%   Optional input variables:
%
%               USE                     OPTIONS                 DEFAULTS
%
% color_int     Color of the box        #interfaces x 3 vector  Black
%               and interfaces          with values [0 1]
%
% line_width    Width of the box        #interfaces x 1 vector  1
%               and interfaces          with width values
%
% plot_gt       Plot geotherm in an     0 don't plot geotherm   1
%               additional subplot      1 plot geotherm

%--------------------------------------------------------------------------
% Function written Miguel Andres-Martinez, 
% 17-02-2015. Email: mandresmartinez87@gmail.com
%--------------------------------------------------------------------------

%==========================================================================
% INPUT
%==========================================================================
% Resolution of the spatial derivative [m]
dz = 10;

% Resolution of the regular grid
res_grid_x = 1000;
res_grid_y = 1000;
y_down_lim = min(GCOORD(2,:));

Fkp = TriScatteredInterp(GCOORD(1,ELEM2NODE(7,:))', ...
    GCOORD(2,ELEM2NODE(7,:))',K(Phases),'natural');
Ftp = TriScatteredInterp(GCOORD(1,:)',GCOORD(2,:)',Temp','natural');

%==========================================================================
% PLOT VERTICAL HEAT FLOW AT THE TOP OF THE MODEL
%==========================================================================
h1 = subplot(2,1,1);
[Topography,topo2nodes] = find_topo(GCOORD,ELEM2NODE,Point_id);
if surf_proc_s
    xp = Topography(1,:);
    yp = Topography(2,:);
    Tbasement = Ftp(xp,yp);
else
    xp = Topography(1,:);
    yp = Topography(2,:);
    Tbasement = Temp(topo2nodes);
end
ylower = yp-dz;

Tlower = Ftp(xp,ylower);
Klower = Fkp(xp,ylower);
Klower(isnan(Klower)) = K(max(Phases));

gradTl = -Klower.*(Tbasement-Tlower)/dz;

plot(xp/1000,gradTl)
title('Basement heat flow')
xlabel('Distance [km]')
ylabel('Heat flow [W/m^2]')

%==========================================================================
% PLOT TEMPERATURE
%==========================================================================
h2 = subplot(2,1,2);
plot_gt = 0;
plot_t

if surf_proc_s
    color_base = [1 0 0];
    line_widthb = 2;
    plot_basement
end

% %==========================================================================
% % CALCULATE VERTICAL HEAT FLOW
% %==========================================================================
% % Building grid
% Xp = unique([min(GCOORD(1,:)):res_grid_x:max(GCOORD(1,:)) max(GCOORD(1,:))]);
% Yp = unique([y_down_lim:res_grid_y:max(GCOORD(2,:)) max(GCOORD(2,:))]);
% [XX,YY] = meshgrid(Xp,Yp);
% 
% % Temperature at the regular grid
% Tr = Ftp(XX,YY);
% % Diffusivity at the regular grid
% Kr = Fkp(XX,YY);
% Kr(isnan(Kr)) = K(max(Phases));
% 
% % Calculate gradient
% gradT = -Kr.*[diff(Tr)/res_grid_y; Xp*NaN];
% 
% %==========================================================================
% % PLOT VERTICAL HEAT FLOW FIELD
% %==========================================================================
% h2 = subplot(2,1,2);
% contourf(XX/1000,YY/1000,gradT,100,'LineStyle','none')
% axis_p = [h2.XLim h2.YLim];
% colormap(jet)
% colorbar
% hold on
% patch([Topography(1,:) Topography(1,end) Topography(1,1) Topography(1,1)]/1000, ...
%     [Topography(2,:) 20000 20000 Topography(2,1)]/1000,'w','LineStyle','none')
% plot_box
% set(h2,'Xlim',axis_p(1:2))
% set(h2,'Ylim',axis_p(3:4))

linkaxes([h2 h1],'x')