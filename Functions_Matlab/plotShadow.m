function plotsd = plotShadow(year,month,day,hour,minute,second)

l_st = 30;
l_loc = 21.78997222;
lat = 38.29113889;
r_real = 149597870691; % earth-sun distance (=1AU)

% import PV and GRH points
PV_points = readtable("G:\Το Drive μου\Scripts MatLab\Shadows\scripts\data\PV_points.csv");
GRH_points = readtable("G:\Το Drive μου\Scripts MatLab\Shadows\scripts\data\GRH_const_1_points.csv");

[~,Zenith,Azimuth] = SunPos(year,month,day,hour,minute,second,l_st,l_loc,lat);

[x_sun,y_sun,z_sun] = SunPosCart(r_real,Zenith,Azimuth);

[x_shad,y_shad] = PVShadowPoint(PV_points.x,PV_points.y,PV_points.z,x_sun,y_sun,z_sun);

shadow_point = table(x_shad,y_shad,'VariableNames',["x","y"]);
shadow_point.Variables = round(shadow_point.Variables,2,"decimals");

% greenhouse
Greenhouse_CU_1 = boundary(GRH_points.x,GRH_points.y); 
Greenhouse_CU_1_value = polyarea(GRH_points.x(Greenhouse_CU_1),GRH_points.y(Greenhouse_CU_1));

% Shadow
Shadow_plane_CU_1 = boundary(shadow_point.x(1:4),shadow_point.y(1:4)); 
Shadow_plane_CU_1_value = polyarea(shadow_point.x(Shadow_plane_CU_1),shadow_point.y(Shadow_plane_CU_1));

%
figure1 = figure('WindowState','maximized','Color','w');
p1 = plot(GRH_points.x,GRH_points.y,'o','Color','k','LineWidth',1.5);
hold on
p2 = plot(GRH_points.x(Greenhouse_CU_1),GRH_points.y(Greenhouse_CU_1),'Color','k','LineWidth',1.5);
hold on
p3 = plot(shadow_point.x(1:4),shadow_point.y(1:4),'o','Color','r','LineWidth',1.5);
hold on
p4 = plot(shadow_point.x(Shadow_plane_CU_1),shadow_point.y(Shadow_plane_CU_1),'Color','r','LineWidth',1.5);
legend('Greenhouse Points','Greenhouse Boundaries','Shadow Points','Shadow Boundaries','Location','northwest')
xlim([-5 20])
ylim([-5 20])

end