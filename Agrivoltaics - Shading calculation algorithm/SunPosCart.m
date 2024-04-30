% Sun position (spherical to cartesian)
function [x_sun,y_sun,z_sun] = SunPosCart(r,Zenith,Azimuth)
    x_sun = r.*sind(Zenith).*cosd(Azimuth);
    y_sun = r.*sind(Zenith).*sind(Azimuth);
    z_sun = r.*cosd(Zenith);
    x_sun((90 - Zenith) < 0) = NaN;
    y_sun((90 - Zenith) < 0) = NaN;
    z_sun((90 - Zenith) < 0) = NaN; 
end