% Shade coordinates calcuations
function [x_sd,y_sd] = PVShadowPoint_xyz(x_pv,y_pv,z_pv,x_sun,y_sun,z_sun,z_sd)
    b = (z_sun - z_pv)./(y_sun - y_pv);
    c = (x_sun - x_pv)./(z_sun - z_pv);
    x_sd = x_pv - c.*(z_pv - z_sd);
    y_sd = y_pv - (1./b).*(z_pv - z_sd);
end