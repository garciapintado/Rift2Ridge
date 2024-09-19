function plot_nodes(MESH,Val,nel)

GCOORD = MESH.GCOORD;
ELEM2NODE = MESH.EL2NOD;

EL2N      = zeros(nel,3); % new connectivity matrix
GCOORD_N  = zeros(2,nel*3);
VV_n      = zeros(3,nel);
EL2N(1,:) = 1:3;

for i=1:nel
    is         = (i-1)*3+1; ie = (i-1)*3+3;
    GCOORD_N(:,is:ie) = GCOORD(:,ELEM2NODE([1 2 3],i));
    EL2N(i,:) = is:ie;
    VV_n(:,i)= Val(ELEM2NODE(1:3,i));
end
patch('faces',EL2N,'vertices',GCOORD_N','facevertexcdata',VV_n(:),'FaceColor','flat')
colorbar

shading interp
end