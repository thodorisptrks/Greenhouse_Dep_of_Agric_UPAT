function [x_sd,y_sd] = PVShadowPoint(x_pv,y_pv,z_pv,x_sun,y_sun,z_sun)
a = (x_sun - x_pv)./(y_sun - y_pv);
b = (z_sun - z_pv)./(y_sun - y_pv);
c = (x_sun - x_pv)./(z_sun - z_pv);

c_z = c.*z_pv;
b_z = z_pv./b;

x_sd = x_pv - c_z;
y_sd = y_pv - b_z;
end

