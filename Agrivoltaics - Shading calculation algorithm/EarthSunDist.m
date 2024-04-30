% Earth-Sun Distance calculations
function r = EarthSunDist(yr,mth,d,hr,mnt,sc)
    % set daytime
    X_year = (datetime(yr,1,1,0,0,0):minutes(1):datetime(yr,12,31,23,59,0))';
    one_min_idx = (1:1:find(X_year == datetime(yr,12,31,23,59,0)))';
    X_inp = datetime(yr,mth,d,hr,mnt,sc,"Format","dd-MMM-uuuu HH:mm:ss");
    one_min = one_min_idx(X_inp == X_year);
    % elleiptical orbit 
    perihelion = 0.983*149597870700; % Perihelion in meters
    alpha = perihelion/(1-0.01671); % elliptical parameter alpha
    gamma = alpha*0.01671; % elliptical parameter gamma
    beta = sqrt(alpha.^2 - gamma.^2); % elliptical parameter beta
    % Earth - Sun Distance
    r = sqrt((alpha.*cosd(0.986/1440.*one_min) - gamma).^2 + (beta.*sind(0.986/1440.*one_min)).^2);
end