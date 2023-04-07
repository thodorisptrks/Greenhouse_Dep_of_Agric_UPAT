function r = EarthSunDist(yr,mth,d,hr,mnt,sc)

if sc~=0
    msg = "Seconds must be zero.";
    error(msg)
end

% set daytime
X1 = datetime(yr,1,1,0,0,0);
X1.Format = "dd-MMM-uuuu HH:mm:ss"; 
X2 = datetime(yr,12,eomday(yr,12),23,50,0);
X2.Format = "dd-MMM-uuuu HH:mm:ss"; 

X3 = (X1:minutes(10):X2)';

ten_min_idx = (1:1:find(X3 == X2))';

X_inp = datetime(yr,mth,d,hr,mnt,sc);
X_inp.Format = "dd-MMM-uuuu HH:mm:ss"; 

ten_min = ten_min_idx(X_inp == X3);

% elleiptical orbit

AU = 149597870700; % astonomical unit to meters

aphelion = 1.0167*AU; % Aphelion in meters
perihelion = 0.98329*AU; % Perihelion in meters

ecc = 0.0167086; % eccentricity

alpha = aphelion/(1+ecc); % elliptical parameter alpha
gamma = alpha*ecc; % elliptical parameter beta
beta = sqrt(alpha.^2 - gamma.^2); % elliptical parameter gamma

% coordinates of sun
x = alpha.*cosd(ten_min/144);
y = beta.*sind(ten_min/144);

% Earth - Sun Distance
r = sqrt(x.^2 + y.^2);
end
