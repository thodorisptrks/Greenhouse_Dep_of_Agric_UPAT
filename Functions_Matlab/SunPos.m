function [Altitude,Zenith,Azimuth] = SunPos(yr,mth,d,hr,mnt,sc,lon_loc,lat,timezone)

if sc~=0
    msg = "Seconds must be zero.";
    error(msg)
end

X = [yr mth d hr mnt sc];
d = datetime(X);
d.Format = "dd-MMM-uuuu HH:mm:ss"; 

% fractional year - gamma  (in radians) 
gamma = (2*pi/365).*(day(d,"dayofyear") - 1 + ((hr - 12)./24));

% equation of time - eq_time (in minutes) 
eq_time = 229.18*(0.000075 + 0.001868.*cos(gamma) - ...
    0.032077.*sin(gamma) - 0.014615.*cos(2.*gamma) - ...
    0.040849.*sin(2.*gamma));

% solar declination - decl (in radians)
decl = 0.006918 - ...
    0.399912.*cos(gamma) + 0.070257.*sin(gamma) -...
    0.006758.*cos(2.*gamma) + 0.000907.*sin(2.*gamma) -...
    0.002697.*cos(3.*gamma) + 0.00148.*sin(3.*gamma);

% time offset - t_offset (in minutes)
L_st = 360 - 15*timezone; % in degrees
L_loc = 360 - lon_loc; % in degrees

time_offset = 4*(L_st - L_loc) + eq_time;

% standard time - st_time (in minutes)
st_time = (hr.*60) + mnt;

% true solar time - tst (in minutes)
tst = st_time + time_offset;

% hour angle - ha (in degrees)
ha = (tst/4) - 180;

% solar zenith angle - theta_z (in degrees)
theta_z = acosd((cosd(lat)*cos(decl).*cosd(ha)) +...
    (sind(lat)*sin(decl)));
Zenith = theta_z;

% solar altitude angle - alpha_s (in degrees)

alpha_s = 90 - theta_z;

% negative values are turned to zero
alpha_s(alpha_s<0) = 0;

Altitude = alpha_s;

% solar azimuth angle - gamma_s (in degrees)

gamma_s_1 = abs(acosd((cosd(theta_z)*sind(lat) -...
    sin(decl))./(sind(theta_z)*cosd(lat))));

for i = 1:length(gamma_s_1)
    if ha(i) > 0
        gamma_s_1(i) = gamma_s_1(i);
    elseif ha(i) < 0
        gamma_s_1(i) = -gamma_s_1(i);
    end
end

gamma_s = gamma_s_1 + 180;
Azimuth = gamma_s;

% display(Altitude);
% display(Zenith);
% display(Azimuth);
end