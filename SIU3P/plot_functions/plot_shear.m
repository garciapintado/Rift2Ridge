% Resolution of the vertical line for the shear stress envelope
res_p = 1; % [km]

TAU_shear = 2*Mu_all.*E2all;

% Shape functions
[IP_X,~] = ip_triangle(nip);
[N,~] = shp_deriv_triangle(IP_X, nnodel);

% GIPS
GIPxp = zeros(size(TAU_shear));
GIPyp = zeros(size(TAU_shear));
for ip = 1:nip
    Ni = N{ip};
    ECOORD_xp = reshape(GCOORD(1,ELEM2NODE(1:6,:)),nnodel,nel);
    ECOORD_yp = reshape(GCOORD(2,ELEM2NODE(1:6,:)),nnodel,nel);
    GIPxp(:,ip) = Ni'*ECOORD_xp;
    GIPyp(:,ip) = Ni'*ECOORD_yp;
end

%MAPPING FROM INTEGRATION POINTS TO NODES - IP CANNOT BE PLOTTED
EL2N      = zeros(nel,3); % new connectivity matrix
GCOORD_N  = zeros(2,nel*3);
TAU_s     = zeros(3,nel);
EL2N(1,:) = 1:3;
nip = 6;
nnodel = 6;
[IP_X, IP_w]    = ip_triangle(nip);

for i=1:nel
    is         = (i-1)*3+1; ie = (i-1)*3+3;
    GCOORD_N(:,is:ie) = GCOORD(:,ELEM2NODE([1 2 3],i));
    EL2N(i,:) = is:ie;
    Dummy     = Nbig'\TAU_shear(i,:)';
    TAU_s(:,i)= Dummy(1:3);
end

subplot(1,2,1)
title('Shear stres map')
xlabel('Distance [km]')
ylabel('Depth [km]')
patch('faces',EL2N,'vertices',GCOORD_N'/1000,'facevertexcdata',TAU_s(:),'FaceColor','flat')
axis tight
shading interp
colorbar
hold on

Interf_member = ismember(ELEM2NODE,find(Point_id>0));
Interf_el = find(sum(Interf_member)==3);
Interf_nod = find(Interf_member(1:3,Interf_el));
Indx = ELEM2NODE(1:3,Interf_el);
Indx = Indx(Interf_nod);

for n = 1:size(Interf_el,2)
    plot(GCOORD(1,[Indx(n*2-1) Indx(n*2)])/km, ...
        GCOORD(2,[Indx(n*2-1) Indx(n*2)])/km,'k','LineWidth',2)
end

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

% Calculating envelope
Envelope = griddata(GIPxp(:),GIPyp(:),TAU_shear(:),X,Y);

subplot(1,2,2)

% GUI
keep_plot = 1;
% Plotting loop
while keep_plot
    % Take pointer input
    [x,y] = ginput(1);
    % If inside the area for plotting plot new envelope
    if (x*1000 >= min(GCOORD(1,:)) && x*1000 <= max(GCOORD(1,:))) && ...
            (y*1000 >= min(GCOORD(2,:)) && y*1000 <= max(GCOORD(2,:)))
        subplot(1,2,1)
        hold on
        % Delete previous line in the subplot 1
        delete(last_line(1));
        % Calculate new vector of X
        X = x*1000*ones(size(Y));
        % Plot new line in subplot 1
        plot(X([1 end])/1000,Y([1 end])/1000,'w')
        % Mark last line for later remove
        last_line = get(gca, 'children');
        hold off
        
        % Calculating envelope
        Envelope = griddata(GIPxp(:),GIPyp(:),TAU_shear(:),X,Y);
        subplot(1,2,2)
        plot(Envelope,Y/km,'k')        
        
    % If outside the area of plotting, finish plotting
    else
        hold off
        keep_plot = 0;
    end
end

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    