% TODO document this function

inter_g = 10*km;

minx = min(GCOORD(1,:));
maxx = max(GCOORD(1,:));
miny = min(GCOORD(2,:));
maxy = max(GCOORD(2,:));

X = minx:inter_g:maxx;
Y = miny:inter_g:maxy;

[XX,YY] = meshgrid(X,Y);

TriFx = TriScatteredInterp(GCOORD(1,:)',GCOORD(2,:)',DISPL(1,:)');
TriFy = TriScatteredInterp(GCOORD(1,:)',GCOORD(2,:)',DISPL(2,:)');

VXp = TriFx(XX,YY);
VYp = TriFy(XX,YY);

sts = streamslice(XX/km,YY/km,VXp,VYp);
set(sts,'Color','black','LineWidth',4)
sts = streamslice(XX/km,YY/km,VXp,VYp);
set(sts,'Color','white')