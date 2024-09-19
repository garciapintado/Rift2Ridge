% Plot mechanical boundary conditions

ivx = Bc_ind_fs(floor((Bc_ind_fs)/2)~=(Bc_ind_fs)/2);
ivy = Bc_ind_fs(floor((Bc_ind_fs)/2)==(Bc_ind_fs)/2);

icx = (ivx+1)/2;
icy = ivy/2;
[ic,inda,indb] = unique([icx icy]);

bcvx = zeros(1,length(ic));
bcvy = zeros(1,length(ic));

valx = Bc_val_fs(floor((Bc_ind_fs)/2)~=(Bc_ind_fs)/2);
valy = Bc_val_fs(floor((Bc_ind_fs)/2)==(Bc_ind_fs)/2);

bcvx(indb(1:length(icx))) = valx;
bcvy(indb(length(icx)+1:end)) = valy;

quiver(GCOORD(1,ic)/km, GCOORD(2,ic)/km, ...
    bcvx, bcvy)%,scale_v,'Color',c_ar);