function [] = plot_ip_val(Val, xip, zip, gc, e2n)

[IP_X,~]    = ip_triangle_m2tri(6);
[NP,~]      = sf_dsf_tri367_N(IP_X,3,'cell');

NL          = [NP{4} NP{5} NP{6}];
Valn          = NL'\Val(:,4:6)';

plot3(xip/1000,zip/1000,Val,'.','MarkerSize',30)
hold on

for n = 1:size(e2n,2)
    patch('Faces',[1 2 3],'Vertices',[gc(:,e2n(1:3,n))'/1000 ...
        Valn(1:3,n)],'EdgeColor','white')
end