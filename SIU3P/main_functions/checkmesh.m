function [qn,amin,amax] = checkmesh(gcoord,nodes)

%check for remeshing
l12 = sqrt((gcoord(1,nodes(1,:))-gcoord(1,nodes(2,:))).^2+(gcoord(2,nodes(1,:))-gcoord(2,nodes(2,:))).^2);
l23 = sqrt((gcoord(1,nodes(2,:))-gcoord(1,nodes(3,:))).^2+(gcoord(2,nodes(2,:))-gcoord(2,nodes(3,:))).^2);
l31 = sqrt((gcoord(1,nodes(3,:))-gcoord(1,nodes(1,:))).^2+(gcoord(2,nodes(3,:))-gcoord(2,nodes(1,:))).^2);

% qn  = min([2.*l12.*l41./(l11.^2+l41.^2) 2.*l11.*l23./(l11.^2+l23.^2) 2.*l23.*l34./(l23.^2+l34.^2) 2.*l34.*l41./(l34.^2+l41.^2)]);

% qn = 0;
v13 = [gcoord(1,nodes(3,:))-gcoord(1,nodes(1,:));gcoord(2,nodes(3,:))-gcoord(2,nodes(1,:))];
v12 = [gcoord(1,nodes(2,:))-gcoord(1,nodes(1,:));gcoord(2,nodes(2,:))-gcoord(2,nodes(1,:))];
v23 = [gcoord(1,nodes(3,:))-gcoord(1,nodes(2,:));gcoord(2,nodes(3,:))-gcoord(2,nodes(2,:))];

aaa = v13.*v12;aa = aaa(1:2:end-1) + aaa(2:2:end);
aaa = sqrt(sum(v13.^2,1)) .*  sqrt(sum(v12.^2,1));
a1  = acos(aa./aaa).*360/(2*pi);
aaa = -v12.*v23;aa = aaa(1:2:end-1) + aaa(2:2:end);
aaa = sqrt(sum(v12.^2,1)) .*  sqrt(sum(v23.^2,1));
a2  = acos(aa./aaa).*360/(2*pi);
a3 = 180 - a1 - a2;
amin=min([a1 a2 a3]);amax=max([a1 a2 a3]);

area = l12.*l31.*sin(a1.*pi./180)./2;
% area1  = l12.*l23.*sin(a2.*pi./180)./2;
qn = min((4*sqrt(3)*area)./(l12.^2 + l23.^2 +l31.^2));
