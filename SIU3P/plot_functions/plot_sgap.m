% Plots the geometry nodes of Sgap
hold on

Sgap_id = ones(1,sum(Geo_id==1)-1);

for ip = 1:((max(Geo_id)-1)/3)
    u1 = 3*ip-1;
    u2 = 3*ip;
    u3 = 3*ip+1;
    Sgap_id = [Sgap_id u1*ones(1,sum(Geo_id==u1)+1)];
    Sgap_id = [Sgap_id u2*ones(1,sum(Geo_id==u2)-1)];
    Sgap_id = [Sgap_id u3*ones(1,sum(Geo_id==u3)+1)];
end

Sgap_g = Sgap_id(Sgap)/3;
Sgap_n = find(Sgap);
inc_sgap = 1:-1:(-(max(Geo_id)-1)/3+2);

second_point_gap = Sgap_n + inc_sgap(Sgap_g);
first_point_gap = second_point_gap-1;

for jp = 1:sum(Sgap)
    plot(GEOMETRY(1,[first_point_gap(jp),second_point_gap(jp)])/1000, ...
        GEOMETRY(2,[first_point_gap(jp),second_point_gap(jp)])/1000, ...
        '-or','LineWidth',3)
end