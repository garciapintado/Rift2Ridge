function elboo = selectElements(GCOORD, EL2NOD, xlims, zlims)
% select those elements whose center is within the queried bounds

[x_ip,~] = ip_triangle(1);
N = shp_triangle(x_ip,3);
GCOx = reshape(GCOORD(1,EL2NOD(1:3,:)),3,[]);
GCOz = reshape(GCOORD(2,EL2NOD(1:3,:)),3,[]);
xc = N' * GCOx;        % element centers - x
zc = N' * GCOz;        % element centers - z 
elboo =   xc >= xlims(1) & xc <= xlims(2) ...
        & zc >= zlims(1) & zc <= zlims(2);
end