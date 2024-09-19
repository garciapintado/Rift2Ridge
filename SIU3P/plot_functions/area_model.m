function [area_mod] = area_model(GEOMETRY,Geo_id)

x = [GEOMETRY(1,Geo_id==1) GEOMETRY(1,Geo_id==2) GEOMETRY(1,Geo_id==5) GEOMETRY(1,Geo_id==8) GEOMETRY(1,Geo_id==9) GEOMETRY(1,Geo_id==10) GEOMETRY(1,Geo_id==7) GEOMETRY(1,Geo_id==4)];
y = [GEOMETRY(2,Geo_id==1) GEOMETRY(2,Geo_id==2) GEOMETRY(2,Geo_id==5) GEOMETRY(2,Geo_id==8) GEOMETRY(2,Geo_id==9) GEOMETRY(2,Geo_id==10) GEOMETRY(2,Geo_id==7) GEOMETRY(2,Geo_id==4)];
plot(x,y,'-','LineWidth',2)

area_mod = polyarea(x,y);