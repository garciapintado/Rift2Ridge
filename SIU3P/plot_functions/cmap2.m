% Generate a color map based in 2 given colors

c1 = [0 0.8 0.5];
c2 = [1 0.2 0.2];
cint = 100;
colormap([(c1(1):((c2(1)-c1(1))/cint):c2(1))' ...
    (c1(2):(c2(2)-c1(2))/cint:c2(2))' (c1(3):(c2(3)-c1(3))/cint:c2(3))'])