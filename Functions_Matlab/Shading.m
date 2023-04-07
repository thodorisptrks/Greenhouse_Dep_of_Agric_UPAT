function [greenhouse_north_shaded_surface,greenhouse_north_shaded_surface_perc] = Shading(yr,mth,d,hr,mnt,scd,timezone,l_loc,lat)

% Altitude, Zenith, Azimuth
[~,Zenith,Azimuth] = SunPos(yr,mth,d,hr,mnt,scd,l_loc,lat,timezone);

% earth-sun distance in meters
r = EarthSunDist(yr,mth,d,hr,mnt,scd);

% Sun position (spherical to cartesian)
[x_sun,y_sun,z_sun] = SunPosCart(r,Zenith,Azimuth);

if isnan(x_sun) || isnan(y_sun) || isnan(z_sun)
    greenhouse_north_shaded_surface = NaN;
    greenhouse_north_shaded_surface_perc = NaN;

elseif ~isnan(x_sun) && ~isnan(y_sun) && ~isnan(z_sun)   
    % Greenhouse and PV dimensions
    grh_width = 3.46;
    grh_length = 16.45;
    gut_hgth = 3.15;
    PV_width = 1.033;
    PV_length = 2.089;
    r_slope = 31.5;

    % Greenhouse Coordinates
    grh_coord = GreenhouseCoord(2,0,grh_width,grh_length,gut_hgth,r_slope);
    
    % Construction unit south
    CU_south = grh_coord{1,1};
    CU_south_array = table2array(grh_coord{1,1});
    CU_south_bounds_2d = boundary(CU_south_array(1:4,1),CU_south_array(1:4,2));
    CU_south_bounds_3d = boundary(CU_south_array);

    % Construction unit north
    CU_north = grh_coord{1,2};
    CU_north_array = table2array(grh_coord{1,2});
    CU_north_bounds_2d = boundary(CU_north_array(1:4,1),CU_north_array(1:4,2)); % mono gia ta 4 prwta shmeia (auta sto edafos)
    CU_north_bounds_3d = boundary(CU_north_array);

    % PV Coordinates for construction unit South
    num_PVs_CU_south = [2,3,6,7];

    % for construction unit south
    for i=num_PVs_CU_south
        PV_x = [(grh_width/2)-PV_width*cosd(r_slope),CU_south.x("Point 9"),...
            (grh_width/2)-PV_width*cosd(r_slope),CU_south.x("Point 9")]';
    
        PV_y(:,i) = [PV_length*i,PV_length*i,PV_length*(i+1),PV_length*(i+1)]';
    
        PV_z = [PV_width*tand(r_slope)*cosd(r_slope) + gut_hgth,CU_south.z("Point 9"),...
            PV_width*tand(r_slope)*cosd(r_slope) + gut_hgth,CU_south.z("Point 9")]';
        
        PV_points_CU_south_tables{i,:} = table(PV_x,PV_y(:,i),PV_z,'VariableNames',["x","y","z"]);
        PV_points_CU_south_arrays{i,:} = [PV_x,PV_y(:,i),PV_z];
    end

    % for the plot (construction unit south)
    for i=num_PVs_CU_south
        X1_south{i} = (PV_points_CU_south_arrays{i}(:,1))';
        Y1_south{i} = (PV_points_CU_south_arrays{i}(:,2))';
        Z1_south{i} = (PV_points_CU_south_arrays{i}(:,3))';
        N = 2;
    
        Xi_south{i} = linspace(min(X1_south{i}),max(X1_south{i}),N);
        Yi_south{i} = linspace(min(Y1_south{i}),max(Y1_south{i}),N);
        ZiX_south{i} = linspace(min(Z1_south{i}),max(Z1_south{i}),N);
        ZiY_south{i} = linspace(min(Z1_south{i}),max(Z1_south{i}),N)';
    
        [Xii_south{i},Yii_south{i}] = meshgrid(Xi_south{i},Yi_south{i});
        Zii_south{i} = meshgrid(ZiX_south{i},ZiY_south{i}) ;
    end

    % for the shadows (construction unit south)
    for i=num_PVs_CU_south
        [x_shad_CU_south{i},y_shad_CU_south{i}] = PVShadowPoint(PV_points_CU_south_arrays{i}(:,1),...
            PV_points_CU_south_arrays{i}(:,2),PV_points_CU_south_arrays{i}(:,3),x_sun,z_sun,z_sun);
    
        X1_shad_CU_south{i} = x_shad_CU_south{i}';
        Y1_shad_CU_south{i} = y_shad_CU_south{i}';
        N = 2;
    
        Xi_shad_CU_south{i} = linspace(min(X1_shad_CU_south{i}),max(X1_shad_CU_south{i}),N);
        Yi_shad_CU_south{i} = linspace(min(Y1_shad_CU_south{i}),max(Y1_shad_CU_south{i}),N);
    
        [Xii_shad_CU_south{i},Yii_shad_CU_south{i}] = meshgrid(Xi_shad_CU_south{i},Yi_shad_CU_south{i});
        Zii_shad_CU_south{i} = zeros(size(Xii_shad_CU_south{i}));
    end

    % PV Coordinates for construction unit North
    num_PVs_CU_north = (0:1:7);
    
    % for construction unit north
    for i=num_PVs_CU_north
        PV_x = [((grh_width/2)-PV_width*cosd(r_slope)) + grh_width,CU_north.x("Point 9"),...
            ((grh_width/2)-PV_width*cosd(r_slope)) + grh_width,CU_north.x("Point 9")]';
    
        PV_y(:,i+1) = [PV_length*i,PV_length*i,PV_length*(i+1),PV_length*(i+1)]';
    
        PV_z = [PV_width*tand(r_slope)*cosd(r_slope) + gut_hgth,CU_north.z("Point 9"),...
            PV_width*tand(r_slope)*cosd(r_slope) + gut_hgth,CU_north.z("Point 9")]';
        
        PV_points_CU_north_tables{i+1,:} = table(PV_x,PV_y(:,i+1),PV_z,'VariableNames',["x","y","z"]);
        PV_points_CU_north_arrays{i+1,:} = [PV_x,PV_y(:,i+1),PV_z];
    end

    % for the plot (construction unit north)
    for i=num_PVs_CU_north+1
        X1_north{i} = (PV_points_CU_north_arrays{i}(:,1))';
        Y1_north{i} = (PV_points_CU_north_arrays{i}(:,2))';
        Z1_north{i} = (PV_points_CU_north_arrays{i}(:,3))';
        N = 2;
    
        Xi_north{i} = linspace(min(X1_north{i}),max(X1_north{i}),N);
        Yi_north{i} = linspace(min(Y1_north{i}),max(Y1_north{i}),N);
        ZiX_north{i} = linspace(min(Z1_north{i}),max(Z1_north{i}),N);
        ZiY_north{i} = linspace(min(Z1_north{i}),max(Z1_north{i}),N)';
    
        [Xii_north{i},Yii_north{i}] = meshgrid(Xi_north{i},Yi_north{i});
        Zii_north{i} = meshgrid(ZiX_north{i},ZiY_north{i}) ;
    end
    
    % for the shadows (construction unit north)
    for i=num_PVs_CU_north+1
        [x_shad_CU_north{i},y_shad_CU_north{i}] = PVShadowPoint(PV_points_CU_north_arrays{i}(:,1),...
            PV_points_CU_north_arrays{i}(:,2),PV_points_CU_north_arrays{i}(:,3),x_sun,z_sun,z_sun);
    
        X1_shad_CU_north{i} = x_shad_CU_north{i}';
        Y1_shad_CU_north{i} = y_shad_CU_north{i}';
        N = 2;
    
        Xi_shad_CU_north{i} = linspace(min(X1_shad_CU_north{i}),max(X1_shad_CU_north{i}),N);
        Yi_shad_CU_north{i} = linspace(min(Y1_shad_CU_north{i}),max(Y1_shad_CU_north{i}),N);
    
        [Xii_shad_CU_north{i},Yii_shad_CU_north{i}] = meshgrid(Xi_shad_CU_north{i},Yi_shad_CU_north{i});
        Zii_shad_CU_north{i} = zeros(size(Xii_shad_CU_north{i}));
    end

    % Greenhouse boundaries and area
    Greenhouse_CU_south_area = polyarea(CU_south_array(CU_south_bounds_2d,1),CU_south_array(CU_south_bounds_2d,2));
    Greenhouse_CU_north_area = polyarea(CU_north_array(CU_north_bounds_2d,1),CU_south_array(CU_north_bounds_2d,2));

    % binary presentation greenhouse
    camvas_greenhouse_south = zeros(25*100);
    camvas_greenhouse_north = zeros(25*100);

    % Greenhouse CU South
    Greenhouse_CU_south_bounds_vals = [CU_south_array(CU_south_bounds_2d,1),CU_south_array(CU_south_bounds_2d,2)]*100; % 
    % oi pragmatikes times twn oriwn toy thermokipiou (oi gwnies) 
    % epi 100 gia na exw akeraies times (metra se ekatosta) 
    % kathws kathe stili grammi antistoixei se ena 1 ekatosto
    Greenhouse_CU_south_bounds_vals(:) = Greenhouse_CU_south_bounds_vals(:) + 500; % 
    % prosthetw 500 gia na afisw swro gurw gurw toso ston x oso kai ston y aksona
    
    camvas_greenhouse_south(Greenhouse_CU_south_bounds_vals(1,2):Greenhouse_CU_south_bounds_vals(3,2),...
        Greenhouse_CU_south_bounds_vals(1,1):Greenhouse_CU_south_bounds_vals(3,1)) = 1; %
    % ta midenika ginontai assoi ekei pou exw thermokipio
    
    CU_south_bounds_2d_new = boundary(Greenhouse_CU_south_bounds_vals(:,1),Greenhouse_CU_south_bounds_vals(:,2)); % ta nea oria
    
    % Greenhouse CU North
    Greenhouse_CU_north_bounds_vals = [CU_north_array(CU_north_bounds_2d,1),CU_north_array(CU_north_bounds_2d,2)]*100; % 
    % oi pragmatikes times twn oriwn toy thermokipiou (oi gwnies) 
    % epi 100 gia na exw akeraies times (metra se ekatosta)
    % kathws kathe stili grammi antistoixei se ena 1 ekatosto
    Greenhouse_CU_north_bounds_vals(:) = Greenhouse_CU_north_bounds_vals(:) + 500; %
    % prosthetw 500 gia na afisw swro gurw gurw toso ston x oso kai ston y aksona
    
    camvas_greenhouse_north(Greenhouse_CU_north_bounds_vals(1,2):Greenhouse_CU_north_bounds_vals(3,2),...
        Greenhouse_CU_north_bounds_vals(1,1):Greenhouse_CU_north_bounds_vals(3,1)) = 1; %
    % ta midenika ginontai assoi ekei pou exw thermokipio
    
    CU_north_bounds_2d_new = boundary(Greenhouse_CU_north_bounds_vals(:,1),Greenhouse_CU_north_bounds_vals(:,2)); % ta nea oria

    % Shadows boundaries and area
    % PVs CU North
    X_whole_shadow_north = [Xii_shad_CU_north{1,1}(1,1),Xii_shad_CU_north{1,1}(1,2),...
        Xii_shad_CU_north{1,8}(2,1),Xii_shad_CU_north{1,8}(2,2)]';
    Y_whole_shadow_north = [Yii_shad_CU_north{1,1}(1,1),Yii_shad_CU_north{1,1}(1,2),...
        Yii_shad_CU_north{1,8}(2,1),Yii_shad_CU_north{1,8}(2,2)]';
    shadow_north_corners = [X_whole_shadow_north,Y_whole_shadow_north];
    
    shadow_boundaries_north = boundary(shadow_north_corners(:,1),shadow_north_corners(:,2)); % north shadow boundaries
    shadow_area_north = polyarea(shadow_north_corners(shadow_boundaries_north,1),...
        shadow_north_corners(shadow_boundaries_north,2)); % north shadow surface
    
    % PVs CU South - West Side
    X_shadow_south_west = [Xii_shad_CU_south{1,2}(1,1),Xii_shad_CU_south{1,2}(1,2),...
        Xii_shad_CU_south{1,3}(2,1),Xii_shad_CU_south{1,3}(2,2)]';
    Y_shadow_south_west = [Yii_shad_CU_south{1,2}(1,1),Yii_shad_CU_south{1,2}(1,2),...
        Yii_shad_CU_south{1,3}(2,1),Yii_shad_CU_south{1,3}(2,2)]';
    shadow_south_west_corners = [X_shadow_south_west,Y_shadow_south_west];
    
    shadow_boundaries_south_west = boundary(shadow_south_west_corners(:,1),shadow_south_west_corners(:,2)); % south-west shadow boundaries
    shadow_area_south_west = polyarea(shadow_south_west_corners(shadow_boundaries_south_west,1),...
        shadow_south_west_corners(shadow_boundaries_south_west,2)); % south-west shadow surface
    
    % PVs CU South - East Side
    X_shadow_south_east = [Xii_shad_CU_south{1,6}(1,1),Xii_shad_CU_south{1,6}(1,2),...
        Xii_shad_CU_south{1,7}(2,1),Xii_shad_CU_south{1,7}(2,2)]';
    Y_shadow_south_east = [Yii_shad_CU_south{1,6}(1,1),Yii_shad_CU_south{1,6}(1,2),...
        Yii_shad_CU_south{1,7}(2,1),Yii_shad_CU_south{1,7}(2,2)]';
    shadow_south_east_corners = [X_shadow_south_east,Y_shadow_south_east];
    
    shadow_boundaries_south_east = boundary(shadow_south_east_corners(:,1),shadow_south_east_corners(:,2)); % south-east shadow boundaries
    shadow_area_south_east = polyarea(shadow_south_east_corners(shadow_boundaries_south_east,1),...
        shadow_south_east_corners(shadow_boundaries_south_east,2)); % south-east shadow surface

    % binary presentation shadows
    camvas_shadow = zeros(25*100);

    % PVs CU North
    shadow_north_bounds_vals = ...
        [round(shadow_north_corners(shadow_boundaries_north,1),2,"decimals"),...
        round(shadow_north_corners(shadow_boundaries_north,2),2,"decimals")]*100; % 
        % oi pragmatikes times twn oriwn toy thermokipiou (oi gwnies)
        % epi 100 gia na exw akeraies times (metatropi metra se ekatosta)
        % kathws kathe stili grammi antistoixei se ena 1 ekatosto
    shadow_north_bounds_vals(:) = shadow_north_bounds_vals(:) + 500; %
    % prosthetw 500 gia na afisw swro gurw gurw toso ston x oso kai ston y aksona

    for i=1:size(shadow_north_bounds_vals,1)
        for j = 1: size(shadow_north_bounds_vals,2)
            if shadow_north_bounds_vals(i,j) <= 0, shadow_north_bounds_vals(i,j) = 0; end
            if shadow_north_bounds_vals(i,j) > 2500, shadow_north_bounds_vals(i,j) = 2500; end
        end
    end

    if shadow_north_bounds_vals(1,2) < shadow_north_bounds_vals(3,2) && shadow_north_bounds_vals(1,1) < shadow_north_bounds_vals(3,1)
        camvas_shadow(shadow_north_bounds_vals(1,2):shadow_north_bounds_vals(3,2),...
            shadow_north_bounds_vals(1,1):shadow_north_bounds_vals(3,1)) = 1;
    
    elseif shadow_north_bounds_vals(1,2) > shadow_north_bounds_vals(3,2) && shadow_north_bounds_vals(1,1) > shadow_north_bounds_vals(3,1)
        camvas_shadow(shadow_north_bounds_vals(3,2):shadow_north_bounds_vals(1,2),...
            shadow_north_bounds_vals(3,1):shadow_north_bounds_vals(1,1)) = 1;
    
    elseif shadow_north_bounds_vals(1,2) < shadow_north_bounds_vals(3,2) && shadow_north_bounds_vals(1,1) > shadow_north_bounds_vals(3,1)
        camvas_shadow(shadow_north_bounds_vals(1,2):shadow_north_bounds_vals(3,2),...
            shadow_north_bounds_vals(3,1):shadow_north_bounds_vals(1,1)) = 1;
    
    elseif shadow_north_bounds_vals(1,2) > shadow_north_bounds_vals(3,2) && shadow_north_bounds_vals(1,1) < shadow_north_bounds_vals(3,1)
        camvas_shadow(shadow_north_bounds_vals(3,2):shadow_north_bounds_vals(1,2),...
            shadow_north_bounds_vals(1,1):shadow_north_bounds_vals(3,1)) = 1;
    end
    % ta midenika ginontai assoi ekei pou exw thermokipio

    shadow_north_boundaries_new = ...
        boundary(shadow_north_bounds_vals(:,1),shadow_north_bounds_vals(:,2)); % ta nea oria

    % PVs CU South West
    shadow_south_west_bounds_vals = ...
        [round(shadow_south_west_corners(shadow_boundaries_south_west,1),2,"decimals"),...
        round(shadow_south_west_corners(shadow_boundaries_south_west,2),2,"decimals")]*100; %
        % oi pragmatikes times twn oriwn toy thermokipiou (oi gwnies)
        % epi 100 gia na exw akeraies times (metra se ekatosta) 
        % kathws kathe stili grammi antistoixei se ena 1 ekatosto
    shadow_south_west_bounds_vals(:) = shadow_south_west_bounds_vals(:) + 500; %
    % prosthetw 500 gia na afisw swro gurw gurw toso ston x oso kai ston y aksona

    for i=1:size(shadow_south_west_bounds_vals,1)
        for j = 1: size(shadow_south_west_bounds_vals,2)
            if shadow_south_west_bounds_vals(i,j) <= 0, shadow_south_west_bounds_vals(i,j) = 0; end
            if shadow_south_west_bounds_vals(i,j) > 2500, shadow_south_west_bounds_vals(i,j) = 2500; end
        end
    end

    if shadow_south_west_bounds_vals(1,2) < shadow_south_west_bounds_vals(3,2) && shadow_south_west_bounds_vals(1,1) < shadow_south_west_bounds_vals(3,1)
        camvas_shadow(shadow_south_west_bounds_vals(1,2):shadow_south_west_bounds_vals(3,2),...
            shadow_south_west_bounds_vals(1,1):shadow_south_west_bounds_vals(3,1)) = 1;
    
    elseif shadow_south_west_bounds_vals(1,2) > shadow_south_west_bounds_vals(3,2) && shadow_south_west_bounds_vals(1,1) > shadow_south_west_bounds_vals(3,1)
        camvas_shadow(shadow_south_west_bounds_vals(3,2):shadow_south_west_bounds_vals(1,2),...
            shadow_south_west_bounds_vals(3,1):shadow_south_west_bounds_vals(1,1)) = 1;
    
    elseif shadow_south_west_bounds_vals(1,2) < shadow_south_west_bounds_vals(3,2) && shadow_south_west_bounds_vals(1,1) > shadow_south_west_bounds_vals(3,1)
        camvas_shadow(shadow_south_west_bounds_vals(1,2):shadow_south_west_bounds_vals(3,2),...
            shadow_south_west_bounds_vals(3,1):shadow_south_west_bounds_vals(1,1)) = 1;
    
    elseif shadow_south_west_bounds_vals(1,2) > shadow_south_west_bounds_vals(3,2) && shadow_south_west_bounds_vals(1,1) < shadow_south_west_bounds_vals(3,1)
        camvas_shadow(shadow_south_west_bounds_vals(3,2):shadow_south_west_bounds_vals(1,2),...
            shadow_south_west_bounds_vals(1,1):shadow_south_west_bounds_vals(3,1)) = 1;
    end
    % ta midenika ginontai assoi ekei pou exw thermokipio

    shadow_south_west_boundaries_new = ...
        boundary(shadow_south_west_bounds_vals(:,1),shadow_south_west_bounds_vals(:,2)); % ta nea oria

        % PVs CU South East
    shadow_south_east_bounds_vals = ...
        [round(shadow_south_east_corners(shadow_boundaries_south_east,1),2,"decimals"),...
        round(shadow_south_east_corners(shadow_boundaries_south_east,2),2,"decimals")]*100; %
        % oi pragmatikes times twn oriwn toy thermokipiou (oi gwnies)
        % epi 100 gia na exw akeraies times (metra se ekatosta) 
        % kathws kathe stili grammi antistoixei se ena 1 ekatosto
    shadow_south_east_bounds_vals(:) = shadow_south_east_bounds_vals(:) + 500; %
    % prosthetw 100 gia na afisw swro gurw gurw toso ston x oso kai ston y aksona

    for i=1:size(shadow_south_east_bounds_vals,1)
        for j = 1: size(shadow_south_east_bounds_vals,2)
            if shadow_south_east_bounds_vals(i,j) <= 0, shadow_south_east_bounds_vals(i,j) = 0; end
            if shadow_south_east_bounds_vals(i,j) > 2500, shadow_south_east_bounds_vals(i,j) = 2500; end
        end
    end

    if shadow_south_east_bounds_vals(1,2) < shadow_south_east_bounds_vals(3,2) && shadow_south_east_bounds_vals(1,1) < shadow_south_east_bounds_vals(3,1)
        camvas_shadow(shadow_south_east_bounds_vals(1,2):shadow_south_east_bounds_vals(3,2),...
            shadow_south_east_bounds_vals(1,1):shadow_south_east_bounds_vals(3,1)) = 1;
    
    elseif shadow_south_east_bounds_vals(1,2) > shadow_south_east_bounds_vals(3,2) && shadow_south_east_bounds_vals(1,1) > shadow_south_east_bounds_vals(3,1)
        camvas_shadow(shadow_south_east_bounds_vals(3,2):shadow_south_east_bounds_vals(1,2),...
            shadow_south_east_bounds_vals(3,1):shadow_south_east_bounds_vals(1,1)) = 1;
    
    elseif shadow_south_east_bounds_vals(1,2) < shadow_south_east_bounds_vals(3,2) && shadow_south_east_bounds_vals(1,1) > shadow_south_east_bounds_vals(3,1)
        camvas_shadow(shadow_south_east_bounds_vals(1,2):shadow_south_east_bounds_vals(3,2),...
            shadow_south_east_bounds_vals(3,1):shadow_south_east_bounds_vals(1,1)) = 1;
    
    elseif shadow_south_east_bounds_vals(1,2) > shadow_south_east_bounds_vals(3,2) && shadow_south_east_bounds_vals(1,1) < shadow_south_east_bounds_vals(3,1)
        camvas_shadow(shadow_south_east_bounds_vals(3,2):shadow_south_east_bounds_vals(1,2),...
            shadow_south_east_bounds_vals(1,1):shadow_south_east_bounds_vals(3,1)) = 1;
    end
    % ta midenika ginontai assoi ekei pou exw thermokipio
    
    shadow_south_east_boundaries_new = ...
        boundary(shadow_south_east_bounds_vals(:,1),shadow_south_east_bounds_vals(:,2)); % ta nea oria

    % Find the shadowing percentage
    canvas_length = length(find(camvas_greenhouse_north(:,:) == 1 | camvas_greenhouse_north(:,:) == 0));
    greenhouse_north_surface = length(find(camvas_greenhouse_north(:,:) == 1));
    total_shading_surface = length(find(camvas_shadow(:,:) == 1));
    
    greenhouse_north_shaded_surface = length(find(camvas_greenhouse_north(:,:) == 1 & camvas_shadow(:,:) == 1));
    greenhouse_north_shaded_surface_perc = (greenhouse_north_shaded_surface/greenhouse_north_surface)*100;
end
end













