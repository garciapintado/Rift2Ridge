function [value_new] = remesh_val(tris,GCOORD,GCOORD_O,value,ELEM2NODE)


nnodel = size(ELEM2NODE,1);
nnod   = length(tris);
x      = reshape(GCOORD(1,ELEM2NODE(1:3,tris)),3,nnod);
y      = reshape(GCOORD(2,ELEM2NODE(1:3,tris)),3,nnod);

xp     = GCOORD_O(1,:);
yp     = GCOORD_O(2,:);


xi     = -(-x(1,:).*yp+x(1,:).*y(3,:)-x(3,:).*y(1,:)+xp.*y(1,:)+x(3,:).*yp-xp.*y(3,:))./(-x(1,:).*y(3,:)+x(1,:).*y(2,:)-x(2,:).*y(1,:)+x(2,:).*y(3,:)+x(3,:).*y(1,:)-x(3,:).*y(2,:));
yi     = (x(1,:).*y(2,:)-x(1,:).*yp-x(2,:).*y(1,:)+x(2,:).*yp+xp.*y(1,:)-xp.*y(2,:))./(-x(1,:).*y(3,:)+x(1,:).*y(2,:)-x(2,:).*y(1,:)+x(2,:).*y(3,:)+x(3,:).*y(1,:)-x(3,:).*y(2,:));


[NN] = shp_triangle([xi' yi'], nnodel);                                    % [nnodel,length(xp)]

if length(tris)>1
    V = value(ELEM2NODE(:,tris));
else
    V = value(ELEM2NODE(:,tris))';
end
value_new = sum(NN .* V);

