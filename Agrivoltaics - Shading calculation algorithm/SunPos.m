% Soalr Zenith Angle and Solar Azimuth Angle calculations
function [Zenith,Azimuth] = SunPos(yr,mth,d,hr,mnt,sc,lon_loc,lat,timezone)
    d = datetime([yr mth d hr mnt sc],"Format","dd-MMM-uuuu HH:mm:ss"); 
    % fractional year - gamma  (in radians)
    if leapyear(yr) == 0
        gamma = (2*pi/365).*(day(d,"dayofyear") - 1 + (((hr - 12) + mnt/60)./24));
    else
        gamma = (2*pi/366).*(day(d,"dayofyear") - 1 + (((hr - 12) + mnt/60)./24));
    end
    % equation of time - eq_time (in minutes) 
    eq_time = 229.18*(0.000075 + 0.001868.*cos(gamma) - ...
	    0.032077.*sin(gamma) - 0.014615.*cos(2.*gamma) - ...
	    0.040849.*sin(2.*gamma));
    % solar declination - decl (in radians)
    decl = 0.006918 - ...
	    0.399912.*cos(gamma) + 0.070257.*sin(gamma) -...
	    0.006758.*cos(2.*gamma) + 0.000907.*sin(2.*gamma) -...
	    0.002697.*cos(3.*gamma) + 0.00148.*sin(3.*gamma);
    % hour angle - ha (in degrees) 
    ha = ((((hr.*60) + mnt) + (4*((360 - 15*timezone) - (360 - lon_loc)) + eq_time))/4) - 180;
    % solar zenith angle - theta_z (in degrees)
    Zenith = acosd((cosd(lat)*cos(decl).*cosd(ha)) + (sind(lat)*sin(decl)));
    % solar azimuth angle - gamma_s (in degrees)
    gamma_s_1 = abs(acosd((cosd(Zenith)*sind(lat) - sin(decl))./(sind(Zenith)*cosd(lat))));
    for k = 1:length(gamma_s_1)
	    if ha(k) > 0, gamma_s_1(k) = gamma_s_1(k);
		    elseif ha(k) < 0, gamma_s_1(k) = -gamma_s_1(k); end
    end
    Azimuth = gamma_s_1 + 180;
end