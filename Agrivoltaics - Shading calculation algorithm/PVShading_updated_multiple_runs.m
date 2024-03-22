function [surface_shade_in_CU_sum_final,surface_shade_in_CU_perc,...
    surface_shade_sc_in_CU_sum_final,surface_shade_sc_in_CU_perc,...
    reply_in_PV_array_final,reply_in_sc_array_final,sc_line_sh,PV_sh,CU_sh] = ...
    PVShading_updated_multiple_runs(yr,mth,d,hr,mnt,sc,timezone,...
    grh_dir,numCU_wdwise,numCU_lnwise,grh_width,grh_length,gut_hght,rigde_hght,...
    l_loc,lat,...
    PV_width,PV_length,...
    sc_length,...
    pos_PV_CU_1,pos_PV_CU_2,...
    z_shad,...
    varargin)
    
    %% datetime limitations
    if d < 1 || d > eomday(yr,mth), error("This day does not exist."); end
    if hr < 0 || hr > 23, error("Hours must be equal to a value between 0 and 23."); end
    if mnt ~= [0:1:59], error("Minutes must be equal to a value between 0 and 59."); end
    if sc ~= 0, error("Seconds must be zero."); end
    if l_loc < 0 || l_loc >360, error("Longitude must be a value between 0 and 360 degrees."); end
    if lat < -90 || lat > 90, error("Longitude must be a value between -90 and 90 degrees."); end

    %% 1. Zenith, Azimuth
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
    
    %% 2. Earth-Sun Distance
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

    %% 3. Sun position (spherical to cartesian)
    function [x_sun,y_sun,z_sun] = SunPosCart(r,Zenith,Azimuth)
	    x_sun = r.*sind(Zenith).*cosd(Azimuth);
	    y_sun = r.*sind(Zenith).*sind(Azimuth);
	    z_sun = r.*cosd(Zenith);
	    x_sun((90 - Zenith) < 0) = NaN;
	    y_sun((90 - Zenith) < 0) = NaN;
	    z_sun((90 - Zenith) < 0) = NaN; 
    end

    %% 4. Find the greenhouse coordinates
    function GRH_const_points = GreenhouseCoord(grh_dir,numCU_wdwise,numCU_lnwise,grh_width,grh_length,gut_hght,rigde_hght)
	    
	    % limitation 1
	    assert(any(strcmp(grh_dir, {'NS','EW'})), 'Greenhouse direction must be either North-South or East-West')
	    
	    % limitation 2
	    function validate_positive_integer_of_CU(direction, numCU)
		    if numCU < 0 || abs(numCU - round(numCU)) > 0.0001, error('%s must be a positive integer', direction); end
	    end
	    validate_positive_integer_of_CU('Number of construction units', numCU_wdwise);
	    validate_positive_integer_of_CU('Number of construction units', numCU_lnwise);
	    if numCU_wdwise == 0 && numCU_lnwise == 0, error('Must define at least one construction unit.'); end
    
	    % x,y,z
	    function GRH_const_points = GrhXYZ(grh_dir,k,l,grh_width,grh_length,gut_hght,rigde_hght)
		    GRH_z_coord = [0,0,0,0,...
			    gut_hght,gut_hght,gut_hght,gut_hght,...
			    rigde_hght,rigde_hght]';
		    % Case 1
		    if grh_dir == 'EW'
			    GRH_x_coord = [grh_width*k,grh_width*(k+1),grh_width*k,grh_width*(k+1),...
				    grh_width*k,grh_width*(k+1),grh_width*k,grh_width*(k+1),...
				    (grh_width/2)*(2*k+1),(grh_width/2)*(2*k+1)]';
				    
			    GRH_y_coord = [grh_length*l,grh_length*l,grh_length*(l+1),grh_length*(l+1),...
				    grh_length*l,grh_length*l,grh_length*(l+1),grh_length*(l+1),...
				    grh_length*l,grh_length*(l+1)]';
				    
			    GRH_const_points = table(GRH_x_coord,GRH_y_coord,GRH_z_coord,'VariableNames',["x","y","z"],...
				    'RowNames',["Point 1","Point 2","Point 3","Point 4","Point 5","Point 6","Point 7","Point 8","Point 9","Point 10"]);
		    % Case 2
		    elseif grh_dir == 'NS'
			    GRH_x_coord = [grh_length*l,grh_length*l,grh_length*(l+1),grh_length*(l+1),...
				    grh_length*l,grh_length*l,grh_length*(l+1),grh_length*(l+1),...
				    grh_length*l,grh_length*(l+1)]';
				    
			    GRH_y_coord = [grh_width*k,grh_width*(k+1),grh_width*k,grh_width*(k+1),...
				    grh_width*k,grh_width*(k+1),grh_width*k,grh_width*(k+1),...
				    (grh_width/2)*(2*k+1),(grh_width/2)*(2*k+1)]';
    
			    GRH_const_points = table(GRH_x_coord,GRH_y_coord,GRH_z_coord,'VariableNames',["x","y","z"],...
				    'RowNames',["Point 1","Point 2","Point 3","Point 4","Point 5","Point 6","Point 7","Point 8","Point 9","Point 10"]);
		    end
	    end
    
	    % Sub-Case 1
	    if numCU_wdwise > 0 && (numCU_lnwise == 0 || numCU_lnwise == 1)
		    for k=0:numCU_wdwise-1
			    for l=0
				    GRH_const_points{l+1,k+1} = GrhXYZ(grh_dir,k,l,grh_width,grh_length,gut_hght,rigde_hght);
			    end
		    end
	    % Sub-Case 2
	    elseif numCU_lnwise > 0 && (numCU_wdwise == 0 || numCU_wdwise == 1)
		    for k=0
			    for l=0:numCU_lnwise-1
				    GRH_const_points{l+1,k+1} = GrhXYZ(grh_dir,k,l,grh_width,grh_length,gut_hght,rigde_hght);
			    end
		    end
	    % Sub-Case 3
	    elseif numCU_wdwise > 1 && numCU_lnwise > 1
		    for k=0:numCU_wdwise-1
			    for l=0:numCU_lnwise-1
				    GRH_const_points{l+1,k+1} = GrhXYZ(grh_dir,k,l,grh_width,grh_length,gut_hght,rigde_hght);
			    end
		    end
	    end
    end

    %% 5. find the shade coordinates
    function [x_sd,y_sd] = PVShadowPoint_xyz(x_pv,y_pv,z_pv,x_sun,y_sun,z_sun,z_sd)
	    b = (z_sun - z_pv)./(y_sun - y_pv);
	    c = (x_sun - x_pv)./(z_sun - z_pv);
	    x_sd = x_pv - c.*(z_pv - z_sd);
	    y_sd = y_pv - (1./b).*(z_pv - z_sd);
    end

    %% Earth - Sun
    % Zenith, Azimuth
    [Zenith,Azimuth] = SunPos(yr,mth,d,hr,mnt,sc,l_loc,lat,timezone);

    % earth-sun distance in meters
    r = EarthSunDist(yr,mth,d,hr,mnt,sc);    

    % Sun position (spherical to cartesian)
    [x_sun,y_sun,z_sun] = SunPosCart(r,Zenith,Azimuth);

    %% Define Greenhouse and PV dimensions
    % Find Greenhouse coordinates (2 tables -> x,y,z for 10 points)
    % Define North and South Construction Unit
        
    % Greenhouse Coordinates
    CU_coord = GreenhouseCoord(grh_dir,numCU_wdwise,numCU_lnwise,grh_width,grh_length,gut_hght,rigde_hght);
    CU_bounds_2d = [1;2;4;3;1];
    % mono gia ta 4 prwta shmeia (auta sto edafos)
    CU_bounds_3d = [6,2,4;1,2,5;1,3,2;1,5,3;3,4,2;7,4,3;7,3,5;6,5,2;4,8,6;7,8,4;5,6,9;5,9,7;10,9,6;10,7,9;8,10,6;7,10,8];
    
    % disp(class(CU_coord{1,1}))
    CU_tables = cell(size(CU_coord));
    CU_arrays = cell(size(CU_coord));
    
    for i=1:size(CU_coord,2)
        for j=1:size(CU_coord,1)
            CU_tables{j,i} = CU_coord{j,i};
            CU_arrays{j,i} = table2array(CU_coord{j,i});
        end
    end

    %% Binary resentation of the Greenhouse
    bounds_greenhouse = [1;2;3;4;1];
    num_CUs = (1:1:numel(CU_arrays))';

    for i=1:size(CU_tables,2)
        for m=1:size(CU_tables,1)
            Greenhouse_corners{m,i} =...
                [(CU_arrays{m,i}(CU_bounds_2d,1)*100) + 1,(CU_arrays{m,i}(CU_bounds_2d,2)*100) + 1]; % the dimensions in cm plus 1
        end
    end

    %% define PV positions
    pos_PV = {pos_PV_CU_1,pos_PV_CU_2};
    
    %% for solar cell lines
    num_sc_lines = 4;
    fi = atand(2*(rigde_hght - gut_hght)/(grh_width));
    d_c = sc_length; % width of each solar cell in m
    bounds_shade = [1;2;4;3;1];
    
    
    %% x and z sc coordinates
    x_z = [6*0.01 , (6*0.01)+d_c , (13.5*0.01)+d_c , (13.5*0.01)+(2*d_c) , (24*0.01)+(2*d_c) ,...
        (24*0.01)+(3*d_c) , (31.5*0.01)+(3*d_c) , (31.5*0.01)+(4*d_c)]';

    y_z = x_z;
    %% if it's night i dont have anything, if it's not continue
    if isnan(x_sun) || isnan(y_sun) || isnan(z_sun)
        surface_shade_in_CU_sum_final = NaN;
        surface_shade_in_CU_perc = NaN;
        surface_shade_sc_in_CU_sum_final = NaN;
        surface_shade_sc_in_CU_perc = NaN;
        reply_in_PV_array_final = NaN;
        reply_in_sc_array = NaN;
        reply_in_sc_array_final = NaN;
        sc_line_sh = NaN;
        PV_sh = NaN;
        CU_sh = NaN;

        % fig_2D_all = figure('WindowState','maximized','Color','w');
        % close(fig_2D_all)
        % 
        % fig_2D_in_grh = figure('WindowState','maximized','Color','w');
        % close(fig_2D_in_grh)
        % 
        % fig_2D_in_CU_north = figure('WindowState','maximized','Color','w');
        % close(fig_2D_in_CU_north)
        % 
        % fig_3D = figure('WindowState','maximized','Color','w');
        % warning("Shading in nighttime does not exist.")
        % close(fig_3D)
    
    elseif ~isnan(x_sun) && ~isnan(y_sun) && ~isnan(z_sun)
        if size(CU_tables) ~= size(pos_PV)
            error("Greenhouse units are not identical to photovoltaic arrays.");
        else
            if PV_width > (sqrt((rigde_hght - gut_hght).^2 + (grh_width/2).^2)) || PV_width <= 0
                error("The width of the PV Unit must be a positive integer and not exceed" + ...
                    " the length of the inclined plane of the roof");
            elseif PV_width <= (sqrt((rigde_hght - gut_hght).^2 + (grh_width/2).^2)) && PV_width > 0
                for n=1:size(pos_PV,2)
                    if grh_length < numel(pos_PV{n})*PV_length
                        error("The PV units exceed the limits of the greenhouse roof");
                    else
                        if z_shad > rigde_hght
                            error("The shadow cannot be calculated over the greenhouse.")
                        else 
                            if grh_dir == 'EW'
                                for i=1:size(CU_tables,2)
                                    for m=1:size(CU_tables,1)
                                        for j=pos_PV{m,i}
                                            
                                            % EAST-WEST ORIENTATION
                                            % X-coordinate for PV for
                                            PV_x_coord(:,j,i,m) = [(cosd(atand(2*(rigde_hght-gut_hght)/grh_width)) *...
                                                (sqrt(((rigde_hght - gut_hght).^2) + ((grh_width.^2)/4)) - PV_width)) + (i-1)*grh_width,...
                                                CU_tables{m,i}.x("Point 9"),...
                                                (cosd(atand(2*(rigde_hght-gut_hght)/grh_width)) *...
                                                (sqrt(((rigde_hght - gut_hght).^2) + ((grh_width.^2)/4)) - PV_width)) + (i-1)*grh_width,...
                                                CU_tables{m,i}.x("Point 9")]';

                                            x_sc(:,j,i,m) = PV_x_coord(1,j,i,m) + x_z .* cosd(fi);

                                            for t=1:num_sc_lines
                                                x_sc_coord(:,j,i,m,t) =...
                                                    vertcat(x_sc((2*t-1),j,i),x_sc((2*t),j,i),x_sc((2*t-1),j,i),x_sc((2*t),j,i));
                                            end

                                            % Y-coordinate
                                            PV_y_coord(:,j,i,m) = [PV_length*(j-1) + (m-1)*grh_length,...
                                                PV_length*(j-1) + (m-1)*grh_length,...
                                                PV_length*j + (m-1)*grh_length,...
                                                PV_length*j + (m-1)*grh_length]';

                                            for t=1:num_sc_lines
                                                y_sc_coord(:,j,i,m,t) = [min(PV_y_coord(:,j,i,m)) + 8.5*0.01,...
                                                    min(PV_y_coord(:,j,i,m)) + 8.5*0.01,...
                                                    max(PV_y_coord(:,j,i,m)) - 8.5*0.01,...                  
                                                    max(PV_y_coord(:,j,i,m)) - 8.5*0.01];
                                            end

                                            % Z-coordinate
                                            PV_z_coord(:,j,i,m) = [(sind(atand(2*(rigde_hght-gut_hght)/grh_width)) *...
                                                (sqrt(((rigde_hght - gut_hght).^2) + ((grh_width.^2)/4)) - PV_width)) + gut_hght,...
                                                CU_tables{m,i}.z("Point 9"),...
                                                (sind(atand(2*(rigde_hght-gut_hght)/grh_width)) *...
                                                (sqrt(((rigde_hght - gut_hght).^2) + ((grh_width.^2)/4)) - PV_width)) + gut_hght,...
                                                CU_tables{m,i}.z("Point 9")]';

                                            z_sc(:,j,i,m) = PV_z_coord(1,j,i,m) + x_z .* sind(fi);

                                            for t=1:num_sc_lines
                                                z_sc_coord(:,j,i,m,t) =...
                                                    vertcat(z_sc((2*t-1),j,i),z_sc((2*t),j,i),z_sc((2*t-1),j,i),z_sc((2*t),j,i));
                                            end

                                            % FOR PLOTTING
                                            PV_x_coord_1{:,j,i,m} = PV_x_coord(:,j,i,m);
                                            PV_y_coord_1{:,j,i,m} = PV_y_coord(:,j,i,m);
                                            PV_z_coord_1{:,j,i,m} = PV_z_coord(:,j,i,m);

                                            for t=1:num_sc_lines
                                                x_sc_coord_1{:,j,i,m,t} = x_sc_coord(:,j,i,m,t);
                                                y_sc_coord_1{:,j,i,m,t} = y_sc_coord(:,j,i,m,t);
                                                z_sc_coord_1{:,j,i,m,t} = z_sc_coord(:,j,i,m,t);
                                            end

                                            X1{m,i} = PV_x_coord_1(:,:,i,m);
                                            Y1{m,i} = PV_y_coord_1(:,:,i,m);
                                            Z1{m,i} = PV_z_coord_1(:,:,i,m);

                                            for t=1:num_sc_lines
                                                X1_sc{m,i}{j}(:,t) = x_sc_coord_1{:,j,i,m,t};
                                                Y1_sc{m,i}{j}(:,t) = y_sc_coord_1{:,j,i,m,t};
                                                Z1_sc{m,i}{j}(:,t) = z_sc_coord_1{:,j,i,m,t};
                                            end

                                            Xi{m,i}{j} = linspace(min(X1{m,i}{:,j}),max(X1{m,i}{:,j}),2);
                                            Yi{m,i}{j} = linspace(min(Y1{m,i}{:,j}),max(Y1{m,i}{:,j}),2);
                                            ZiX{m,i}{j} = linspace(min(Z1{m,i}{j}),max(Z1{m,i}{j}),2);
                                            ZiY{m,i}{j} = linspace(min(Z1{m,i}{j}),max(Z1{m,i}{j}),2)';

                                            for t=1:num_sc_lines
                                                Xi_sc{m,i}{j}{t,:} = linspace(min(X1_sc{m,i}{j}(:,t)),max(X1_sc{m,i}{j}(:,t)),2);
                                                Yi_sc{m,i}{j}{t,:} = linspace(min(Y1_sc{m,i}{j}(:,t)),max(Y1_sc{m,i}{j}(:,t)),2);
                                                Zi_sc_X{m,i}{j}{t,:} = linspace(min(Z1_sc{m,i}{j}(:,t)),max(Z1_sc{m,i}{j}(:,t)),2);
                                                Zi_sc_Y{m,i}{j}{:,t} = linspace(min(Z1_sc{m,i}{j}(:,t)),max(Z1_sc{m,i}{j}(:,t)),2);
                                            end

                                            [Xii{m,i}{j},Yii{m,i}{j}] = meshgrid(Xi{m,i}{j},Yi{m,i}{j});
                                            Zii{m,i}{j} = meshgrid(ZiX{m,i}{j},ZiY{m,i}{j});

                                            for t=1:num_sc_lines
                                                [Xii_sc{m,i}{j}{t,t},Yii_sc{m,i}{j}{t,t}] =...
                                                    meshgrid(Xi_sc{m,i}{j}{t,:},Yi_sc{m,i}{j}{t,:});
                                                Zii_sc{m,i}{j}{t,t} = meshgrid(Zi_sc_X{m,i}{j}{t,:},Zi_sc_Y{m,i}{j}{:,t});
                                            end

                                            Xii_sc{m,i}{j} =...
                                                {Xii_sc{m,i}{j}{1,1};Xii_sc{m,i}{j}{2,2};Xii_sc{m,i}{j}{3,3};Xii_sc{m,i}{j}{4,4}};
                                            Yii_sc{m,i}{j} =...
                                                {Yii_sc{m,i}{j}{1,1};Yii_sc{m,i}{j}{2,2};Yii_sc{m,i}{j}{3,3};Yii_sc{m,i}{j}{4,4}};
                                            Zii_sc{m,i}{j} =...
                                                {Zii_sc{m,i}{j}{1,1};Zii_sc{m,i}{j}{2,2};Zii_sc{m,i}{j}{3,3};Zii_sc{m,i}{j}{4,4}};

                                            % FOR SHADING
                                            [x_shad{m,i}{j},y_shad{m,i}{j}] =...
                                                PVShadowPoint_xyz(X1{m,i}{j},Y1{m,i}{j},Z1{m,i}{j},x_sun,y_sun,z_sun,z_shad);

                                            Xi_shad{m,i}{j} = [x_shad{m,i}{j}(1,1),x_shad{m,i}{j}(2,1);x_shad{m,i}{j}(3,1),x_shad{m,i}{j}(4,1)];
                                            Yi_shad{m,i}{j} = [y_shad{m,i}{j}(1,1),y_shad{m,i}{j}(2,1);y_shad{m,i}{j}(3,1),y_shad{m,i}{j}(4,1)];
                                            Zi_shad{m,i}{j} = zeros(size(Xi_shad{m,i}{j}));

                                            if Zi_shad{m,i}{j}(:) == 0
                                                Zi_shad{m,i}{j}(:) = z_shad;
                                            end

                                            for t=1:num_sc_lines
                                                [x_shad_sc{m,i}{j}(:,t),y_shad_sc{m,i}{j}(:,t)] =...
                                                    PVShadowPoint_xyz(X1_sc{m,i}{j}(:,t),Y1_sc{m,i}{j}(:,t),Z1_sc{m,i}{j}(:,t),x_sun,y_sun,z_sun,z_shad);

                                                X1_shad_sc{m,i}{j}{t,:} = x_shad_sc{m,i}{j}(:,t);
                                                Y1_shad_sc{m,i}{j}{t,:} = y_shad_sc{m,i}{j}(:,t);
                                                
                                                Xi_shad_sc{m,i}{j}{t} = [X1_shad_sc{m,i}{j}{t}(1,1),X1_shad_sc{m,i}{j}{t}(2,1);...
                                                    X1_shad_sc{m,i}{j}{t}(3,1),X1_shad_sc{m,i}{j}{t}(4,1)];
                                                Yi_shad_sc{m,i}{j}{t} = [Y1_shad_sc{m,i}{j}{t}(1,1),Y1_shad_sc{m,i}{j}{t}(2,1);...
                                                    Y1_shad_sc{m,i}{j}{t}(3,1),Y1_shad_sc{m,i}{j}{t}(4,1)];
                                                Zi_shad_sc{m,i}{j}{t} = zeros(size(Xi_shad_sc{m,i}{j}{t}));

                                                if Zi_shad_sc{m,i}{j}{t}(:) == 0
                                                    Zi_shad_sc{m,i}{j}{t}(:) = z_shad;
                                                end
                                            end                                           
                                        end
                                    end
                                end

                                % Boundaries of the whole PVs Shading
                                Xi_shad_cm = Xi_shad;
                                Xi_shad_cm_in_all = Xi_shad_cm;

                                Yi_shad_cm = Yi_shad;
                                Yi_shad_cm_in_all = Yi_shad_cm;
                                
                                Xi_shad_sc_cm = Xi_shad_sc;
                                Xi_shad_sc_cm_in_all = Xi_shad_sc_cm;

                                Yi_shad_sc_cm = Yi_shad_sc;
                                Yi_shad_sc_cm_in_all = Yi_shad_sc_cm;

                                % convert all the dimensions of the shade formed by the solar cells into cm
                                for i=1:size(CU_tables,2)
                                    for m=1:size(CU_tables,1)
                                        for j=pos_PV{m,i}
                                            if ~isempty(Xi_shad{m,i}{j})
                                                Xi_shad_cm{m,i}{j}(:) = (round(Xi_shad{m,i}{j}(:),2,"decimals") * 100) + 1;
                                                Xi_shad_cm_in_all{m,i}{j}(:) = (round(Xi_shad{m,i}{j}(:),2,"decimals") * 100) + 1;
                                            end

                                            Xi_shad_cm_in_all{m,i}{j}(Xi_shad_cm_in_all{m,i}{j}(:) < 0) = 1;
                                            Xi_shad_cm_in_all{m,i}{j}(Xi_shad_cm_in_all{m,i}{j}(:) > int16((size(CU_tables,2)*grh_width*100)+1)) =...
                                                int16((size(CU_tables,2)*grh_width*100)+1);

                                            if ~isempty(Yi_shad{m,i}{j})
                                                Yi_shad_cm{m,i}{j}(:) = (round(Yi_shad{m,i}{j}(:),2,"decimals") * 100) + 1;
                                                Yi_shad_cm_in_all{m,i}{j}(:) = (round(Yi_shad{m,i}{j}(:),2,"decimals") * 100) + 1;
                                            end

                                            Yi_shad_cm_in_all{m,i}{j}(Yi_shad_cm_in_all{m,i}{j}(:) < 0) = 1;
                                            Yi_shad_cm_in_all{m,i}{j}(Yi_shad_cm_in_all{m,i}{j}(:) > int16((grh_length*100)+1)) =...
                                                int16((grh_length*100)+1);

                                            for t=1:num_sc_lines
                                                if ~isempty(Xi_shad_sc{m,i}{j})
                                                    Xi_shad_sc_cm{m,i}{j}{1,t}(:) = (round(Xi_shad_sc{m,i}{j}{1,t}(:),2,"decimals") * 100) + 1;
                                                    Xi_shad_sc_cm_in_all{m,i}{j}{1,t}(:) = (round(Xi_shad_sc{m,i}{j}{1,t}(:),2,"decimals") * 100) + 1;
                                                end

                                                Xi_shad_sc_cm_in_all{m,i}{j}{1,t}(Xi_shad_sc_cm_in_all{m,i}{j}{1,t}(:) < 0) = 1;
                                                Xi_shad_sc_cm_in_all{m,i}{j}{1,t}(Xi_shad_sc_cm_in_all{m,i}{j}{1,t}(:) > int16((size(CU_tables,2)*grh_width*100)+1)) =...
                                                    int16((size(CU_tables,2)*grh_width*100)+1);

                                                if ~isempty(Yi_shad_sc{m,i}{j})
                                                    Yi_shad_sc_cm{m,i}{j}{1,t}(:) = (round(Yi_shad_sc{m,i}{j}{1,t}(:),2,"decimals") * 100) + 1;
                                                    Yi_shad_sc_cm_in_all{m,i}{j}{1,t}(:) = (round(Yi_shad_sc{m,i}{j}{1,t}(:),2,"decimals") * 100) + 1;
                                                end
                                                
                                                Yi_shad_sc_cm_in_all{m,i}{j}{1,t}(Yi_shad_sc_cm_in_all{m,i}{j}{1,t}(:) < 0) = 1;
                                                Yi_shad_sc_cm_in_all{m,i}{j}{1,t}(Yi_shad_sc_cm_in_all{m,i}{j}{1,t}(:) > int16((grh_length*100)+1)) =...
                                                    int16((grh_length*100)+1);
                                            end
                                        end
                                    end
                                end

                                % make all integers
                                for i=1:size(CU_tables,2)
                                    for m=1:size(CU_tables,1)
                                        for j=pos_PV{m,i}
                                            for t=1:num_sc_lines
                                                Xi_shad_sc_cm_in_all{m,i}{j}{1,t}(:) = int16(Xi_shad_sc_cm_in_all{m,i}{j}{1,t}(:));
                                                Yi_shad_sc_cm_in_all{m,i}{j}{1,t}(:) = int16(Yi_shad_sc_cm_in_all{m,i}{j}{1,t}(:));
                                            end
                                        end
                                    end
                                end

                                Xi_shad_cm_in_CU = Xi_shad_cm_in_all;
                                Yi_shad_cm_in_CU = Yi_shad_cm_in_all;
                                Xi_shad_sc_cm_in_CU = Xi_shad_sc_cm_in_all;
                                Yi_shad_sc_cm_in_CU = Yi_shad_sc_cm_in_all;

                                for i=1:size(CU_tables,2)
                                    for m=1:size(CU_tables,1)
                                        for j=pos_PV{m,i}
                                            for t=1:num_sc_lines
                                                if ~isempty(Xi_shad_cm_in_CU{m,i}{j})
                                                    Xi_shad_cm_in_CU{m,i}{j}(Xi_shad_cm_in_CU{m,i}{j}(:) < min(Greenhouse_corners{1,2}(:,1))) = min(Greenhouse_corners{1,2}(:,1));
					                                Yi_shad_cm_in_CU{m,i}{j}(Xi_shad_cm_in_CU{m,i}{j}(:) < min(Greenhouse_corners{1,2}(:,2))) = min(Greenhouse_corners{1,2}(:,2));
                                
					                                Xi_shad_cm_in_CU{m,i}{j}(Xi_shad_cm_in_CU{m,i}{j}(:) > max(Greenhouse_corners{1,2}(:,1))) = max(Greenhouse_corners{1,2}(:,1));
					                                Yi_shad_cm_in_CU{m,i}{j}(Yi_shad_cm_in_CU{m,i}{j}(:) > max(Greenhouse_corners{1,2}(:,2))) = max(Greenhouse_corners{1,2}(:,2));
                                    
					                                Xi_shad_sc_cm_in_CU{m,i}{j}{:,t}(Xi_shad_sc_cm_in_CU{m,i}{j}{:,t}(:) < min(Greenhouse_corners{1,2}(:,1))) = min(Greenhouse_corners{1,2}(:,1));
					                                Yi_shad_sc_cm_in_CU{m,i}{j}{:,t}(Yi_shad_sc_cm_in_CU{m,i}{j}{:,t}(:) < min(Greenhouse_corners{1,2}(:,2))) = min(Greenhouse_corners{1,2}(:,2));
                                
					                                Xi_shad_sc_cm_in_CU{m,i}{j}{:,t}(Xi_shad_sc_cm_in_CU{m,i}{j}{:,t}(:) > max(Greenhouse_corners{1,2}(:,1))) = max(Greenhouse_corners{1,2}(:,1));
					                                Yi_shad_sc_cm_in_CU{m,i}{j}{:,t}(Yi_shad_sc_cm_in_CU{m,i}{j}{:,t}(:) > max(Greenhouse_corners{1,2}(:,2))) = max(Greenhouse_corners{1,2}(:,2));
                                                else
					                                Xi_shad_cm_in_CU{m,i}{j} = Xi_shad_cm_in_CU{m,i}{j};
                                                    Yi_shad_cm_in_CU{m,i}{j} = Yi_shad_cm_in_CU{m,i}{j};
					                                
					                                Xi_shad_sc_cm_in_CU{m,i}{j}{:,t} = Xi_shad_sc_cm_in_CU{m,i}{j};
                                                    Yi_shad_sc_cm_in_CU{m,i}{j}{:,t} = Yi_shad_sc_cm_in_CU{m,i}{j};
                                                end
                                            end
                                        end
                                    end
                                end

                                % find the surfaces 
                                for i=1:size(CU_tables,2)
                                    for m=1:size(CU_tables,1)
                                        for j=pos_PV{m,i}
                                            for t=1:num_sc_lines
                                                if ~isempty(Xi_shad_cm_in_CU{m,i}{j})
                                                    [~,surface_grh(m,i)] = boundary(Greenhouse_corners{m,i});
                                                    [~,surface_PVs(i,j,m)] = boundary(PV_x_coord(:,j,i,m),PV_y_coord(:,j,i,m));
                                                    [~,surface_scs(t,j,i,m)] = boundary(x_sc_coord(:,j,i,m,t),y_sc_coord(:,j,i,m,t));
                                    
                                                    [~,surface_shade_in_CU{m,i}(j)] = boundary(Xi_shad_cm_in_CU{m,i}{j}(:),Yi_shad_cm_in_CU{m,i}{j}(:));
                                                    [~,surface_shade_sc_in_CU{m,i}(t,j)] =...
                                                        boundary(Xi_shad_sc_cm_in_CU{m,i}{j}{1,t}(:),Yi_shad_sc_cm_in_CU{m,i}{j}{1,t}(:));
                                                else
                                                    surface_shade_in_CU{m,i}(j) = NaN;
                                                    surface_shade_sc_in_CU{m,i}(t,j) = NaN;
                                                end
                                                
                                                surface_grh(m,i) = surface_grh(m,i).*0.0001; % se m2
                                                surface_PVs(i,j) = surface_PVs(i,j)/cosd(fi); % se m2
                                                surface_scs(t,j,i,m) = surface_scs(t,j,i,m)/cosd(fi); % se m2
                                
                                                surface_shade_in_CU{m,i}(1,j) = surface_shade_in_CU{m,i}(1,j)*0.0001;
                                                surface_shade_sc_in_CU{m,i}(t,j) = surface_shade_sc_in_CU{m,i}(t,j)*0.0001;
                                            end
                                        end
                                
                                        surface_shade_in_CU_sum(m,i) = sum(surface_shade_in_CU{m,i},'all','omitnan');
                                        surface_shade_sc_in_CU_sum(m,i) = sum(surface_shade_sc_in_CU{m,i},'all','omitnan');
                                    end
                                end
                                       
                                surface_shade_in_CU_sum_final = sum(surface_shade_in_CU_sum,'all'); % se m2
                                surface_shade_sc_in_CU_sum_final = sum(surface_shade_sc_in_CU_sum,'all'); % se m2
                                
                                surface_shade_in_CU_perc = (surface_shade_in_CU_sum_final/surface_grh(1,2))*100;
                                surface_shade_sc_in_CU_perc = (surface_shade_sc_in_CU_sum_final/surface_grh(1,2))*100;

                                % check if the pyrantometer is in or outside of the shaded surface
                                if isempty(varargin)
                                    reply_in_PV_array_final = NaN;
                                    reply_in_sc_array_final = NaN;
                                    sc_line_sh = NaN;
                                    PV_sh = NaN;
                                    CU_sh = NaN;
                                    return;
                                elseif ~isempty(varargin) && mod(length(varargin),2) ~= 0
                                    error("To define the station's coordinates x and y coordinates are needed.")
                                elseif ~isempty(varargin) && mod(length(varargin),2) == 0
                                    for p = 1:2:numel(varargin)
                                        x_st = varargin{p};
                                        y_st = varargin{p+1};
                        
                                        x_st_1 = repelem(x_st,2,2);
                                        y_st_1 = repelem(y_st,2,2);

                                        % point inside the whole PV array
                                        s_ref = cell(size(CU_tables));
                                        s_t_total = cell(size(CU_tables));

                                        for i = 1:size(CU_tables,2)
                                            for m = 1:size(CU_tables,1)
                                                for j = pos_PV{m,i}
                                                    s_t_total{m,i}{j} = cell(1,1);

                                                    [~,s_ref{m,i}{j}] = boundary(reshape(Xi_shad_cm_in_all{m,i}{j},[],1),...
                                                        reshape(Yi_shad_cm_in_all{m,i}{j},[],1));

                                                    [~,s_t_1{m,i}{j}] = boundary(reshape(horzcat(Xi_shad_cm_in_all{m,i}{j}(:,1),x_st_1(:,1)),[],1),...
                                                        reshape(horzcat(Yi_shad_cm_in_all{m,i}{j}(:,1),y_st_1(:,1)),[],1));

                                                    [~,s_t_2{m,i}{j}] = boundary(reshape(horzcat(Xi_shad_cm_in_all{m,i}{j}(1,:),x_st_1(1,:)),[],1),...
                                                        reshape(horzcat(Yi_shad_cm_in_all{m,i}{j}(1,:),y_st_1(1,:)),[],1));

                                                    [~,s_t_3{m,i}{j}] = boundary(reshape(horzcat(Xi_shad_cm_in_all{m,i}{j}(:,2),x_st_1(:,2)),[],1),...
                                                        reshape(horzcat(Yi_shad_cm_in_all{m,i}{j}(:,2),y_st_1(:,2)),[],1));

                                                    [~,s_t_4{m,i}{j}] = boundary(reshape(horzcat(Xi_shad_cm_in_all{m,i}{j}(2,:),x_st_1(2,:)),[],1),...
                                                        reshape(horzcat(Yi_shad_cm_in_all{m,i}{j}(2,:),y_st_1(2,:)),[],1));

                                                    s_t_total{m,i}{j} = s_t_1{m,i}{j} + s_t_2{m,i}{j} + s_t_3{m,i}{j} + s_t_4{m,i}{j};

                                                    s_ref{m,i}{j}(:) = round(s_ref{m,i}{j}(:),1);
                                                    s_t_total{m,i}{j}(:) = round(s_t_total{m,i}{j}(:),1);

                                                    s_ref{m,i}{j}(:) = int16(s_ref{m,i}{j}(:));
                                                    s_t_total{m,i}{j}(:) = int16(s_t_total{m,i}{j}(:));

                                                    if s_ref{m,i}{j} ~= 0 && s_ref{m,i}{j} == s_t_total{m,i}{j}
                                                        reply_in_PV{m,i}{j} = 1;
                                                    elseif s_ref{m,i}{j} == 0 || s_ref{m,i}{j} ~= s_t_total{m,i}{j}
                                                        reply_in_PV{m,i}{j} = 0;
                                                    end
                                                end
                                            end
                                        end


                                        for i = 1:size(CU_tables,2)
                                            for m = 1:size(CU_tables,1)
                                                for j = pos_PV{m,i}
                                                    reply_in_PV_array(:,j,m,i) = reply_in_PV{m,i}{j};
                                                end
                                            end
                                        end

                                        if any(reply_in_PV_array(:) == 1)
                                            reply_in_PV_array_final = 1;
                                        else
                                            reply_in_PV_array_final = 0;
                                        end

                                        % point inside the solar cell
                                        s_sc_ref = cell(size(CU_tables));
                                        s_sc_t_total = cell(size(CU_tables));
                        
                                        for i = 1:size(CU_tables,2)
                                            for m = 1:size(CU_tables,1)
                                                for j = pos_PV{m,i}
                                                    s_sc_t_total{m,i}{j} = cell(1,num_sc_lines);

                                                    for t = 1:num_sc_lines
                                                        [~,s_sc_ref{m,i}{j}{t}] = boundary(reshape(Xi_shad_sc_cm_in_all{m,i}{j}{t},[],1),...
                                                            reshape(Yi_shad_sc_cm_in_all{m,i}{j}{t},[],1));
                        
                                                        [~,s_sc_t_1{m,i}{j}{t}] = boundary(reshape(horzcat(Xi_shad_sc_cm_in_all{m,i}{j}{t}(:,1),x_st_1(:,1)),[],1),...
                                                            reshape(horzcat(Yi_shad_sc_cm_in_all{m,i}{j}{t}(:,1),y_st_1(:,1)),[],1));
                        
                                                        [~,s_sc_t_2{m,i}{j}{t}] = boundary(reshape(horzcat(Xi_shad_sc_cm_in_all{m,i}{j}{t}(1,:),x_st_1(1,:)),[],1),...
                                                            reshape(horzcat(Yi_shad_sc_cm_in_all{m,i}{j}{t}(1,:),y_st_1(1,:)),[],1));
                        
                                                        [~,s_sc_t_3{m,i}{j}{t}] = boundary(reshape(horzcat(Xi_shad_sc_cm_in_all{m,i}{j}{t}(:,2),x_st_1(:,2)),[],1),...
                                                            reshape(horzcat(Yi_shad_sc_cm_in_all{m,i}{j}{t}(:,2),y_st_1(:,2)),[],1));
                        
                                                        [~,s_sc_t_4{m,i}{j}{t}] = boundary(reshape(horzcat(Xi_shad_sc_cm_in_all{m,i}{j}{t}(2,:),x_st_1(2,:)),[],1),...
                                                            reshape(horzcat(Yi_shad_sc_cm_in_all{m,i}{j}{t}(2,:),y_st_1(2,:)),[],1));
                        
                                                        s_sc_t_total{m,i}{j}{t} = s_sc_t_1{m,i}{j}{t} + s_sc_t_2{m,i}{j}{t} + s_sc_t_3{m,i}{j}{t} + s_sc_t_4{m,i}{j}{t};

                                                        s_sc_ref{m,i}{j}{1,t}(:) = round(s_sc_ref{m,i}{j}{1,t}(:),1);
                                                        s_sc_t_total{m,i}{j}{1,t}(:) = round(s_sc_t_total{m,i}{j}{1,t}(:),1);

                                                        s_sc_ref{m,i}{j}{1,t}(:) = int16(s_sc_ref{m,i}{j}{1,t}(:));
                                                        s_sc_t_total{m,i}{j}{1,t}(:) = int16(s_sc_t_total{m,i}{j}{1,t}(:));
                        
                                                        if s_sc_ref{m,i}{j}{t} ~= 0 && s_sc_ref{m,i}{j}{t} == s_sc_t_total{m,i}{j}{t}
                                                            reply_in_sc{m,i}{j}{t} = 1;
                                                        elseif s_sc_ref{m,i}{j}{t} == 0 || s_sc_ref{m,i}{j}{t} ~= s_sc_t_total{m,i}{j}{t}
                                                            reply_in_sc{m,i}{j}{t} = 0;
                                                        end
                                                    end                        
                                                end                    
                                            end
                                        end
                                    end

                                    for i = 1:size(CU_tables,2)
                                        for m = 1:size(CU_tables,1)
                                            for j = pos_PV{m,i}
                                                for t = 1:num_sc_lines
                                                    reply_in_sc_array(t,j,m,i) = reply_in_sc{m,i}{j}{t};                            
                                                end
                                            end
                                        end
                                    end
                        
                                    if any(reply_in_sc_array(:) == 1)
                                        reply_in_sc_array_final = 1;
                                        [sc_line_sh,PV_sh,CU_sh] = ind2sub(size(reply_in_sc_array),find(reply_in_sc_array(:) == 1));
                                    else
                                        reply_in_sc_array_final = 0;
                                        sc_line_sh = NaN;
                                        PV_sh = NaN;
                                        CU_sh = NaN;
                                    end
                                end
                                                            
                            elseif grh_dir == 'NS'
                                for i=1:size(CU_tables,2)
                                    for m=1:size(CU_tables,1)
                                        for j=pos_PV{m,i}
                                            
                                            % X-coordinate for PV for NORTH-SOUTH ORIENTATION
                                            PV_x_coord(:,j,i,m) = [PV_length*(j-1) + (m-1)*grh_length,...
                                                PV_length*(j-1) + (m-1)*grh_length,...
                                                PV_length*j + (m-1)*grh_length,...
                                                PV_length*j + (m-1)*grh_length]';

                                            for t=1:num_sc_lines
                                                x_sc_coord(:,j,i,m,t) = [min(PV_x_coord(:,j,i,m)) + 8.5*0.01,...
                                                    min(PV_x_coord(:,j,i,m)) + 8.5*0.01,...
                                                    max(PV_x_coord(:,j,i,m)) - 8.5*0.01,...
                                                    max(PV_x_coord(:,j,i,m)) - 8.5*0.01];
                                            end

                                            % Y-coordinate
                                            PV_y_coord(:,j,i,m) = [(cosd(atand(2*(rigde_hght-gut_hght)/grh_width)) *...
                                                (sqrt(((rigde_hght - gut_hght).^2) + ((grh_width.^2)/4)) - PV_width)) + (i-1)*grh_width,...
                                                CU_tables{m,i}.y("Point 9"),...
                                                (cosd(atand(2*(rigde_hght-gut_hght)/grh_width)) *...
                                                (sqrt(((rigde_hght - gut_hght).^2) + ((grh_width.^2)/4)) - PV_width)) + (i-1)*grh_width,...
                                                CU_tables{m,i}.y("Point 9")]';

                                            y_sc(:,j,i,m) = PV_y_coord(1,j,i,m) + y_z .* cosd(fi);

                                            for t=1:num_sc_lines
                                                y_sc_coord(:,j,i,m,t) =...
                                                    vertcat(y_sc((2*t-1),j,i),y_sc((2*t),j,i),y_sc((2*t-1),j,i),y_sc((2*t),j,i));
                                            end

                                            % Z-coordinate
                                            PV_z_coord(:,j,i,m) = [(sind(atand(2*(rigde_hght-gut_hght)/grh_width)) *...
                                                (sqrt(((rigde_hght - gut_hght).^2) + ((grh_width.^2)/4)) - PV_width)) + gut_hght,...
                                                CU_tables{m,i}.z("Point 9"),...
                                                (sind(atand(2*(rigde_hght-gut_hght)/grh_width)) *...
                                                (sqrt(((rigde_hght - gut_hght).^2) + ((grh_width.^2)/4)) - PV_width)) + gut_hght,...
                                                CU_tables{m,i}.z("Point 9")]';

                                            z_sc(:,j,i,m) = PV_z_coord(1,j,i,m) + y_z .* sind(fi);

                                            for t=1:num_sc_lines
                                                z_sc_coord(:,j,i,m,t) =...
                                                    vertcat(z_sc((2*t-1),j,i),z_sc((2*t),j,i),z_sc((2*t-1),j,i),z_sc((2*t),j,i));
                                            end

                                            % FOR PLOTTING
                                            PV_x_coord_1{:,j,i,m} = PV_x_coord(:,j,i,m);
                                            PV_y_coord_1{:,j,i,m} = PV_y_coord(:,j,i,m);
                                            PV_z_coord_1{:,j,i,m} = PV_z_coord(:,j,i,m);

                                            for t=1:num_sc_lines
                                                x_sc_coord_1{:,j,i,m,t} = x_sc_coord(:,j,i,m,t);
                                                y_sc_coord_1{:,j,i,m,t} = y_sc_coord(:,j,i,m,t);
                                                z_sc_coord_1{:,j,i,m,t} = z_sc_coord(:,j,i,m,t);
                                            end

                                            X1{m,i} = PV_x_coord_1(:,:,i,m);
                                            Y1{m,i} = PV_y_coord_1(:,:,i,m);
                                            Z1{m,i} = PV_z_coord_1(:,:,i,m);

                                            for t=1:num_sc_lines
                                                X1_sc{m,i}{j}(:,t) = x_sc_coord_1{:,j,i,m,t};
                                                Y1_sc{m,i}{j}(:,t) = y_sc_coord_1{:,j,i,m,t};
                                                Z1_sc{m,i}{j}(:,t) = z_sc_coord_1{:,j,i,m,t};
                                            end

                                            Xi{m,i}{j} = linspace(min(X1{m,i}{:,j}),max(X1{m,i}{:,j}),2);
                                            Yi{m,i}{j} = linspace(min(Y1{m,i}{:,j}),max(Y1{m,i}{:,j}),2);
                                            ZiX{m,i}{j} = linspace(min(Z1{m,i}{j}),max(Z1{m,i}{j}),2);
                                            ZiY{m,i}{j} = linspace(min(Z1{m,i}{j}),max(Z1{m,i}{j}),2)';

                                            for t=1:num_sc_lines
                                                Xi_sc{m,i}{j}{t,:} = linspace(min(X1_sc{m,i}{j}(:,t)),max(X1_sc{m,i}{j}(:,t)),2);
                                                Yi_sc{m,i}{j}{t,:} = linspace(min(Y1_sc{m,i}{j}(:,t)),max(Y1_sc{m,i}{j}(:,t)),2);
                                                Zi_sc_X{m,i}{j}{t,:} = linspace(min(Z1_sc{m,i}{j}(:,t)),max(Z1_sc{m,i}{j}(:,t)),2);
                                                Zi_sc_Y{m,i}{j}{:,t} = linspace(min(Z1_sc{m,i}{j}(:,t)),max(Z1_sc{m,i}{j}(:,t)),2);
                                            end

                                            [Xii{m,i}{j},Yii{m,i}{j}] = meshgrid(Xi{m,i}{j},Yi{m,i}{j});
                                            Zii{m,i}{j} = meshgrid(ZiX{m,i}{j},ZiY{m,i}{j});

                                            Xii{m,i}{j} = Xii{m,i}{j}';
                                            Yii{m,i}{j} = Yii{m,i}{j}';

                                            for t=1:num_sc_lines
                                                [Xii_sc{m,i}{j}{t,t},Yii_sc{m,i}{j}{t,t}] =...
                                                    meshgrid(Xi_sc{m,i}{j}{t,:},Yi_sc{m,i}{j}{t,:});
                                                Zii_sc{m,i}{j}{t,t} = meshgrid(Zi_sc_X{m,i}{j}{t,:},Zi_sc_Y{m,i}{j}{:,t});
                                            end

                                            Xii_sc{m,i}{j} =...
                                                {Xii_sc{m,i}{j}{1,1}';Xii_sc{m,i}{j}{2,2}';Xii_sc{m,i}{j}{3,3}';Xii_sc{m,i}{j}{4,4}'};
                                            Yii_sc{m,i}{j} =...
                                                {Yii_sc{m,i}{j}{1,1}';Yii_sc{m,i}{j}{2,2}';Yii_sc{m,i}{j}{3,3}';Yii_sc{m,i}{j}{4,4}'};
                                            Zii_sc{m,i}{j} =...
                                                {Zii_sc{m,i}{j}{1,1};Zii_sc{m,i}{j}{2,2};Zii_sc{m,i}{j}{3,3};Zii_sc{m,i}{j}{4,4}};

                                            % FOR SHADING
                                            [x_shad{m,i}{j},y_shad{m,i}{j}] =...
                                                PVShadowPoint_xyz(X1{m,i}{j},Y1{m,i}{j},Z1{m,i}{j},x_sun,y_sun,z_sun,z_shad);

                                			Xi_shad{m,i}{j} = [x_shad{m,i}{j}(1,1),x_shad{m,i}{j}(2,1);x_shad{m,i}{j}(3,1),x_shad{m,i}{j}(4,1)];
                                            Yi_shad{m,i}{j} = [y_shad{m,i}{j}(1,1),y_shad{m,i}{j}(2,1);y_shad{m,i}{j}(3,1),y_shad{m,i}{j}(4,1)];
                                            Zi_shad{m,i}{j} = zeros(size(Xi_shad{m,i}{j}));

                                            if Zi_shad{m,i}{j}(:) == 0
                                                Zi_shad{m,i}{j}(:) = z_shad;
                                            end

                                            for t=1:num_sc_lines
                                                [x_shad_sc{m,i}{j}(:,t),y_shad_sc{m,i}{j}(:,t)] =...
                                                    PVShadowPoint_xyz(X1_sc{m,i}{j}(:,t),Y1_sc{m,i}{j}(:,t),Z1_sc{m,i}{j}(:,t),x_sun,y_sun,z_sun,z_shad);

                                                X1_shad_sc{m,i}{j}{t,:} = x_shad_sc{m,i}{j}(:,t);
                                                Y1_shad_sc{m,i}{j}{t,:} = y_shad_sc{m,i}{j}(:,t);

                                                Xi_shad_sc{m,i}{j}{t} = [X1_shad_sc{m,i}{j}{t}(1,1),X1_shad_sc{m,i}{j}{t}(2,1);...
                                                    X1_shad_sc{m,i}{j}{t}(3,1),X1_shad_sc{m,i}{j}{t}(4,1)];
                                                Yi_shad_sc{m,i}{j}{t} = [Y1_shad_sc{m,i}{j}{t}(1,1),Y1_shad_sc{m,i}{j}{t}(2,1);...
                                                    Y1_shad_sc{m,i}{j}{t}(3,1),Y1_shad_sc{m,i}{j}{t}(4,1)];
                                                Zi_shad_sc{m,i}{j}{t} = zeros(size(Xi_shad_sc{m,i}{j}{t}));

                                                if Zi_shad_sc{m,i}{j}{t}(:) == 0
                                                    Zi_shad_sc{m,i}{j}{t}(:) = z_shad;
                                                end                                                                                                                                                
                                            end
                                        end
                                    end
                                end

                                % Boundaries of the whole PVs Shading 
                                Xi_shad_cm = Xi_shad;
                                Xi_shad_cm_in_all = Xi_shad_cm;

                                Yi_shad_cm = Yi_shad;
                                Yi_shad_cm_in_all = Yi_shad_cm;
                                
                                Xi_shad_sc_cm = Xi_shad_sc;
                                Xi_shad_sc_cm_in_all = Xi_shad_sc_cm;
                                
                                Yi_shad_sc_cm = Yi_shad_sc;
                                Yi_shad_sc_cm_in_all = Yi_shad_sc_cm;

                                % convert all the dimensions of the shade formed by the solar cells into cm
                                for i=1:size(CU_tables,2)
                                    for m=1:size(CU_tables,1)
                                        for j=pos_PV{m,i}
                                            if ~isempty(Xi_shad{m,i}{j})
                                                Xi_shad_cm{m,i}{j}(:) = (round(Xi_shad{m,i}{j}(:),2,"decimals") * 100) + 1;
                                                Xi_shad_cm_in_all{m,i}{j}(:) = (round(Xi_shad{m,i}{j}(:),2,"decimals") * 100) + 1;
                                            end
                                            Xi_shad_cm_in_all{m,i}{j}(Xi_shad_cm_in_all{m,i}{j}(:) < 0) = 1;
                                            Xi_shad_cm_in_all{m,i}{j}(Xi_shad_cm_in_all{m,i}{j}(:) > int16((size(CU_tables,2)*grh_length*100)+1)) =...
                                                int16((size(CU_tables,2)*grh_length*100)+1);

                                            if ~isempty(Yi_shad{m,i}{j})
                                                Yi_shad_cm{m,i}{j}(:) = (round(Yi_shad{m,i}{j}(:),2,"decimals") * 100) + 1;
                                                Yi_shad_cm_in_all{m,i}{j}(:) = (round(Yi_shad{m,i}{j}(:),2,"decimals") * 100) + 1;
                                            end
                                            Yi_shad_cm_in_all{m,i}{j}(Yi_shad_cm_in_all{m,i}{j}(:) < 0) = 1;
                                            Yi_shad_cm_in_all{m,i}{j}(Yi_shad_cm_in_all{m,i}{j}(:) > int16((grh_width*100)+1)) =...
                                                int16((grh_width*100)+1);            

                                            for t=1:num_sc_lines
                                                if ~isempty(Xi_shad_sc{m,i}{j})
                                                    Xi_shad_sc_cm{m,i}{j}{1,t}(:) = (round(Xi_shad_sc{m,i}{j}{1,t}(:),2,"decimals") * 100) + 1;
                                                    Xi_shad_sc_cm_in_all{m,i}{j}{1,t}(:) = (round(Xi_shad_sc{m,i}{j}{1,t}(:),2,"decimals") * 100) + 1;
                                                end
                                                Xi_shad_sc_cm_in_all{m,i}{j}{1,t}(Xi_shad_sc_cm_in_all{m,i}{j}{1,t}(:) < 0) = 1;
                                
                                                Xi_shad_sc_cm_in_all{m,i}{j}{1,t}(Xi_shad_sc_cm_in_all{m,i}{j}{1,t}(:) > int16((size(CU_tables,2)*grh_length*100)+1)) =...
                                                    int16((size(CU_tables,2)*grh_length*100)+1);
                                
                                                if ~isempty(Yi_shad_sc{m,i}{j})
                                                    Yi_shad_sc_cm{m,i}{j}{1,t}(:) = (round(Yi_shad_sc{m,i}{j}{1,t}(:),2,"decimals") * 100) + 1;
                                                    Yi_shad_sc_cm_in_all{m,i}{j}{1,t}(:) = (round(Yi_shad_sc{m,i}{j}{1,t}(:),2,"decimals") * 100) + 1;
                                                end
                                                Yi_shad_sc_cm_in_all{m,i}{j}{1,t}(Yi_shad_sc_cm_in_all{m,i}{j}{1,t}(:) < 0) = 1;
                                
                                                Yi_shad_sc_cm_in_all{m,i}{j}{1,t}(Yi_shad_sc_cm_in_all{m,i}{j}{1,t}(:) > int16((grh_width*100)+1)) =...
                                                    int16((grh_width*100)+1);
                                            end
                                        end
                                    end
                                end

                                % make all integers
                                for i=1:size(CU_tables,2)
                                    for m=1:size(CU_tables,1)
                                        for j=pos_PV{m,i}
                                            for t=1:num_sc_lines
                                                Xi_shad_sc_cm_in_all{m,i}{j}{1,t}(:) = int16(Xi_shad_sc_cm_in_all{m,i}{j}{1,t}(:));
                                                Yi_shad_sc_cm_in_all{m,i}{j}{1,t}(:) = int16(Yi_shad_sc_cm_in_all{m,i}{j}{1,t}(:));
                                            end
                                        end
                                    end
                                end

                                Xi_shad_cm_in_CU = Xi_shad_cm_in_all;
                                Yi_shad_cm_in_CU = Yi_shad_cm_in_all;
                                Xi_shad_sc_cm_in_CU = Xi_shad_sc_cm_in_all;
                                Yi_shad_sc_cm_in_CU = Yi_shad_sc_cm_in_all;
                                
                                for i=1:size(CU_tables,2)
                                    for m=1:size(CU_tables,1)
                                        for j=pos_PV{m,i}
                                            for t=1:num_sc_lines
                                                if ~isempty(Xi_shad_cm_in_CU{m,i}{j})
                                                    Xi_shad_cm_in_CU{m,i}{j}(Xi_shad_cm_in_CU{m,i}{j}(:) < min(Greenhouse_corners{1,1}(:,1))) = min(Greenhouse_corners{1,1}(:,1));
				                                    Yi_shad_cm_in_CU{m,i}{j}(Xi_shad_cm_in_CU{m,i}{j}(:) < min(Greenhouse_corners{1,1}(:,2))) = min(Greenhouse_corners{1,1}(:,2));
                                
				                                    Xi_shad_cm_in_CU{m,i}{j}(Xi_shad_cm_in_CU{m,i}{j}(:) > max(Greenhouse_corners{1,1}(:,1))) = max(Greenhouse_corners{1,1}(:,1));
				                                    Yi_shad_cm_in_CU{m,i}{j}(Yi_shad_cm_in_CU{m,i}{j}(:) > max(Greenhouse_corners{1,1}(:,2))) = max(Greenhouse_corners{1,1}(:,2));
                                    
				                                    Xi_shad_sc_cm_in_CU{m,i}{j}{:,t}(Xi_shad_sc_cm_in_CU{m,i}{j}{:,t}(:) < min(Greenhouse_corners{1,1}(:,1))) = min(Greenhouse_corners{1,1}(:,1));
				                                    Yi_shad_sc_cm_in_CU{m,i}{j}{:,t}(Yi_shad_sc_cm_in_CU{m,i}{j}{:,t}(:) < min(Greenhouse_corners{1,1}(:,2))) = min(Greenhouse_corners{1,1}(:,2));
                                
				                                    Xi_shad_sc_cm_in_CU{m,i}{j}{:,t}(Xi_shad_sc_cm_in_CU{m,i}{j}{:,t}(:) > max(Greenhouse_corners{1,1}(:,1))) = max(Greenhouse_corners{1,1}(:,1));
				                                    Yi_shad_sc_cm_in_CU{m,i}{j}{:,t}(Yi_shad_sc_cm_in_CU{m,i}{j}{:,t}(:) > max(Greenhouse_corners{1,1}(:,2))) = max(Greenhouse_corners{1,1}(:,2));
                                                else
				                                    Xi_shad_cm_in_CU{m,i}{j} = Xi_shad_cm_in_CU{m,i}{j};
                                                    Yi_shad_cm_in_CU{m,i}{j} = Yi_shad_cm_in_CU{m,i}{j};
				                                    
				                                    Xi_shad_sc_cm_in_CU{m,i}{j}{:,t} = Xi_shad_sc_cm_in_CU{m,i}{j};
                                                    Yi_shad_sc_cm_in_CU{m,i}{j}{:,t} = Yi_shad_sc_cm_in_CU{m,i}{j};
                                                end
                                            end
                                        end
                                    end
                                end

                                % find the surfaces
                                for i=1:size(CU_tables,2)
                                    for m=1:size(CU_tables,1)
                                        for j=pos_PV{m,i}
                                            for t=1:num_sc_lines
                                                if ~isempty(Xi_shad_cm_in_CU{m,i}{j})
                                                    [~,surface_grh(m,i)] = boundary(Greenhouse_corners{m,i});
                                                    [~,surface_PVs(i,j,m)] = boundary(PV_x_coord(:,j,i,m),PV_y_coord(:,j,i,m));
                                                    [~,surface_scs(t,j,i,m)] = boundary(x_sc_coord(:,j,i,m,t),y_sc_coord(:,j,i,m,t));
                                    
                                                    [~,surface_shade_in_CU{m,i}(j)] = boundary(Xi_shad_cm_in_CU{m,i}{j}(:),Yi_shad_cm_in_CU{m,i}{j}(:));
                                                    [~,surface_shade_sc_in_CU{m,i}(t,j)] =...
                                                        boundary(Xi_shad_sc_cm_in_CU{m,i}{j}{1,t}(:),Yi_shad_sc_cm_in_CU{m,i}{j}{1,t}(:));
                                                else
                                                    surface_shade_in_CU{m,i}(j) = NaN;
                                                    surface_shade_sc_in_CU{m,i}(t,j) = NaN;
                                                end
                                                
                                                surface_grh(m,i) = surface_grh(m,i).*0.0001; % se m2
                                                surface_PVs(i,j) = surface_PVs(i,j)/cosd(fi); % se m2
                                                surface_scs(t,j,i,m) = surface_scs(t,j,i,m)/cosd(fi); % se m2
                                
                                                surface_shade_in_CU{m,i}(1,j) = surface_shade_in_CU{m,i}(1,j)*0.0001;
                                                surface_shade_sc_in_CU{m,i}(t,j) = surface_shade_sc_in_CU{m,i}(t,j)*0.0001;
                                            end
                                        end
                                
                                        surface_shade_in_CU_sum(m,i) = sum(surface_shade_in_CU{m,i},'all','omitnan');
                                        surface_shade_sc_in_CU_sum(m,i) = sum(surface_shade_sc_in_CU{m,i},'all','omitnan');
                                    end
                                end
                                       
                                surface_shade_in_CU_sum_final = sum(surface_shade_in_CU_sum,'all'); % se m2
                                surface_shade_sc_in_CU_sum_final = sum(surface_shade_sc_in_CU_sum,'all'); % se m2
                                
                                surface_shade_in_CU_perc = (surface_shade_in_CU_sum_final/surface_grh(1,2))*100;
                                surface_shade_sc_in_CU_perc = (surface_shade_sc_in_CU_sum_final/surface_grh(1,2))*100;

                                % check if the pyranometer is in or outside of the shaded surface

                                if isempty(varargin)
                                    reply_in_PV_array_final = NaN;
                                    reply_in_sc_array_final = NaN;
                                    sc_line_sh = NaN;
                                    PV_sh = NaN;
                                    CU_sh = NaN;
                                    return;
                                elseif ~isempty(varargin) && mod(length(varargin),2) ~= 0
                                    error("To define the station's coordinates x and y coordinates are needed.")
                                elseif ~isempty(varargin) && mod(length(varargin),2) == 0
                                    for p = 1:2:numel(varargin)
                                        x_st = varargin{p};
                                        y_st = varargin{p+1};
                        
                                        x_st_1 = repelem(x_st,2,2);
                                        y_st_1 = repelem(y_st,2,2);

                                        % point inside the whole PV array
                                        s_ref = cell(size(CU_tables));
                                        s_t_total = cell(size(CU_tables));

                                        for i = 1:size(CU_tables,2)
                                            for m = 1:size(CU_tables,1)
                                                for j = pos_PV{m,i}
                                                    s_t_total{m,i}{j} = cell(1,1);

                                                    [~,s_ref{m,i}{j}] = boundary(reshape(Xi_shad_cm_in_all{m,i}{j},[],1),...
                                                        reshape(Yi_shad_cm_in_all{m,i}{j},[],1));

                                                    [~,s_t_1{m,i}{j}] = boundary(reshape(horzcat(Xi_shad_cm_in_all{m,i}{j}(:,1),x_st_1(:,1)),[],1),...
                                                        reshape(horzcat(Yi_shad_cm_in_all{m,i}{j}(:,1),y_st_1(:,1)),[],1));

                                                    [~,s_t_2{m,i}{j}] = boundary(reshape(horzcat(Xi_shad_cm_in_all{m,i}{j}(1,:),x_st_1(1,:)),[],1),...
                                                        reshape(horzcat(Yi_shad_cm_in_all{m,i}{j}(1,:),y_st_1(1,:)),[],1));

                                                    [~,s_t_3{m,i}{j}] = boundary(reshape(horzcat(Xi_shad_cm_in_all{m,i}{j}(:,2),x_st_1(:,2)),[],1),...
                                                        reshape(horzcat(Yi_shad_cm_in_all{m,i}{j}(:,2),y_st_1(:,2)),[],1));

                                                    [~,s_t_4{m,i}{j}] = boundary(reshape(horzcat(Xi_shad_cm_in_all{m,i}{j}(2,:),x_st_1(2,:)),[],1),...
                                                        reshape(horzcat(Yi_shad_cm_in_all{m,i}{j}(2,:),y_st_1(2,:)),[],1));

                                                    s_t_total{m,i}{j} = s_t_1{m,i}{j} + s_t_2{m,i}{j} + s_t_3{m,i}{j} + s_t_4{m,i}{j};

                                                    s_ref{m,i}{j}(:) = round(s_ref{m,i}{j}(:),1);
                                                    s_t_total{m,i}{j}(:) = round(s_t_total{m,i}{j}(:),1);                                                    

                                                    s_ref{m,i}{j}(:) = int16(s_ref{m,i}{j}(:));
                                                    s_t_total{m,i}{j}(:) = int16(s_t_total{m,i}{j}(:));

                                                    if s_ref{m,i}{j} ~= 0 && s_ref{m,i}{j} == s_t_total{m,i}{j}
                                                        reply_in_PV{m,i}{j} = 1;
                                                    elseif s_ref{m,i}{j} == 0 || s_ref{m,i}{j} ~= s_t_total{m,i}{j}
                                                        reply_in_PV{m,i}{j} = 0;
                                                    end
                                                end
                                            end
                                        end


                                        for i = 1:size(CU_tables,2)
                                            for m = 1:size(CU_tables,1)
                                                for j = pos_PV{m,i}
                                                    reply_in_PV_array(:,j,m,i) = reply_in_PV{m,i}{j};
                                                end
                                            end
                                        end

                                        if any(reply_in_PV_array(:) == 1)
                                            reply_in_PV_array_final = 1;
                                        else
                                            reply_in_PV_array_final = 0;
                                        end
                                        
                                        % point inside the solar cell 
                                        s_sc_ref = cell(size(CU_tables));
                                        s_sc_t_total = cell(size(CU_tables));
                        
                                        for i = 1:size(CU_tables,2)
                                            for m = 1:size(CU_tables,1)
                                                for j = pos_PV{m,i}
                                                    s_sc_t_total{m,i}{j} = cell(1,num_sc_lines);

                                                    for t = 1:num_sc_lines
                                                        [~,s_sc_ref{m,i}{j}{t}] = boundary(reshape(Xi_shad_sc_cm_in_all{m,i}{j}{t},[],1),...
                                                            reshape(Yi_shad_sc_cm_in_all{m,i}{j}{t},[],1));
                        
                                                        [~,s_sc_t_1{m,i}{j}{t}] = boundary(reshape(horzcat(Xi_shad_sc_cm_in_all{m,i}{j}{t}(:,1),x_st_1(:,1)),[],1),...
                                                            reshape(horzcat(Yi_shad_sc_cm_in_all{m,i}{j}{t}(:,1),y_st_1(:,1)),[],1));
                        
                                                        [~,s_sc_t_2{m,i}{j}{t}] = boundary(reshape(horzcat(Xi_shad_sc_cm_in_all{m,i}{j}{t}(1,:),x_st_1(1,:)),[],1),...
                                                            reshape(horzcat(Yi_shad_sc_cm_in_all{m,i}{j}{t}(1,:),y_st_1(1,:)),[],1));
                        
                                                        [~,s_sc_t_3{m,i}{j}{t}] = boundary(reshape(horzcat(Xi_shad_sc_cm_in_all{m,i}{j}{t}(:,2),x_st_1(:,2)),[],1),...
                                                            reshape(horzcat(Yi_shad_sc_cm_in_all{m,i}{j}{t}(:,2),y_st_1(:,2)),[],1));
                        
                                                        [~,s_sc_t_4{m,i}{j}{t}] = boundary(reshape(horzcat(Xi_shad_sc_cm_in_all{m,i}{j}{t}(2,:),x_st_1(2,:)),[],1),...
                                                            reshape(horzcat(Yi_shad_sc_cm_in_all{m,i}{j}{t}(2,:),y_st_1(2,:)),[],1));
                        
                                                        s_sc_t_total{m,i}{j}{t} = s_sc_t_1{m,i}{j}{t} + s_sc_t_2{m,i}{j}{t} + s_sc_t_3{m,i}{j}{t} + s_sc_t_4{m,i}{j}{t};

                                                        s_sc_ref{m,i}{j}{1,t}(:) = round(s_sc_ref{m,i}{j}{1,t}(:),1);
                                                        s_sc_t_total{m,i}{j}{1,t}(:) = round(s_sc_t_total{m,i}{j}{1,t}(:),1);

                                                        s_sc_ref{m,i}{j}{1,t}(:) = int16(s_sc_ref{m,i}{j}{1,t}(:));
                                                        s_sc_t_total{m,i}{j}{1,t}(:) = int16(s_sc_t_total{m,i}{j}{1,t}(:));
                                                       
                                                        if s_sc_ref{m,i}{j}{t} ~= 0 && s_sc_ref{m,i}{j}{t} == s_sc_t_total{m,i}{j}{t}
                                                            reply_in_sc{m,i}{j}{t} = 1;
                                                        elseif s_sc_ref{m,i}{j}{t} == 0 || s_sc_ref{m,i}{j}{t} ~= s_sc_t_total{m,i}{j}{t}
                                                            reply_in_sc{m,i}{j}{t} = 0;
                                                        end
                                                    end                        
                                                end                    
                                            end
                                        end
                                    end

                                    for i = 1:size(CU_tables,2)
                                        for m = 1:size(CU_tables,1)
                                            for j = pos_PV{m,i}
                                                for t = 1:num_sc_lines
                                                    reply_in_sc_array(t,j,m,i) = reply_in_sc{m,i}{j}{t};                            
                                                end
                                            end
                                        end
                                    end
                        
                                    if any(reply_in_sc_array(:) == 1)
                                        reply_in_sc_array_final = 1;
                                        [sc_line_sh,PV_sh,CU_sh] = ind2sub(size(reply_in_sc_array),find(reply_in_sc_array(:) == 1));
                                    else
                                        reply_in_sc_array_final = 0;
                                        sc_line_sh = NaN;
                                        PV_sh = NaN;
                                        CU_sh = NaN;
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end
