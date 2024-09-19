% Plots an element at the beginning and at the ending of a time step,
% together with the integrations points at the beginning and ending of the
% time step, and the backwards path of the ending time-step integration
% points in order to see if the final integration points backtrack to the
% initial ones.

%% Run until UPDATE OF COORDINATES (DISPL) in the main file

choosen_element = 116; % Choose an element to analize

plot(GCOORD(1,[ELEM2NODE(1:3,choosen_element)' ELEM2NODE(1,choosen_element)']),...
    GCOORD(2,[ELEM2NODE(1:3,choosen_element)' ELEM2NODE(1,choosen_element)']))
hold on
displ       = reshape(Vel,[ndof,nnod]);
GCOORD_new  = GCOORD + displ*dt;
plot(GCOORD_new(1,[ELEM2NODE(1:3,choosen_element)' ELEM2NODE(1,choosen_element)']),...
    GCOORD_new(2,[ELEM2NODE(1:3,choosen_element)' ELEM2NODE(1,choosen_element)']),'k')
hold on

ecoord_x        = reshape( GCOORD(1,ELEM2NODE(:,choosen_element)), 7, 1);
ecoord_y        = reshape( GCOORD(2,ELEM2NODE(:,choosen_element)), 7, 1);

[IP_X, IP_w]    = ip_triangle(nip);
[   N, dNdu]    = shp_deriv_triangle(IP_X, 7);

for ip = 1:nip
    Ni      =        N{ip};
    
    gip_x   = Ni'*ecoord_x; %x and y coordinate of the integrations point
    gip_y   = Ni'*ecoord_y;
    plot(gip_x,gip_y,'xb')
    
   hold on
end



ecoord_x_new        = reshape( GCOORD_new(1,ELEM2NODE(:,choosen_element)), 7, 1);
ecoord_y_new        = reshape( GCOORD_new(2,ELEM2NODE(:,choosen_element)), 7, 1);

for ip = 1:nip
    Ni      =        N{ip};
    
    gip_x   = Ni'*ecoord_x_new; %x and y coordinate of the integrations point
    gip_y   = Ni'*ecoord_y_new;
    plot(gip_x,gip_y,'ok')
    
   hold on
end

displ_x        = reshape( DISPL(1,ELEM2NODE(:,choosen_element)), 7, 1);
displ_y        = reshape( DISPL(2,ELEM2NODE(:,choosen_element)), 7, 1);

for ip = 1:nip
    Ni      =        N{ip};
    
    gip_x   = Ni'*ecoord_x; %x and y coordinate of the integrations point
    gip_y   = Ni'*ecoord_y;
    vel_x   = Ni'*displ_x; %x and y coordinate of the integrations point
    vel_y   = Ni'*displ_y;
    plot([gip_x; vel_x*dt+gip_x],[gip_y; vel_y*dt+gip_y],'r')
    
   hold on
end