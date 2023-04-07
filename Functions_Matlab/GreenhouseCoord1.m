function [GRH_const_points,PV_points] = ...
    GreenhouseCoord1(numCU_wdwise,numCU_lnwise,grh_width,grh_length,gut_hght,r_slope,PV_roof_side,PV_width,PV_length)

if numCU_wdwise < 0 || (numCU_wdwise - round(numCU_wdwise,0,"decimals")) ~= 0
    msg_1 = "Number of construction units must be positive integer";
    error(msg_1)
elseif numCU_lnwise < 0 || (numCU_lnwise - round(numCU_lnwise,0,"decimals")) ~= 0
    error(msg_1)
elseif numCU_wdwise == 0 && numCU_lnwise ==0
    msg_2 = "Must be defined at least one construction unit.";
    error(msg_2)
end

if numCU_wdwise > 0 && (numCU_lnwise == 0 || numCU_lnwise == 1)
    for i=0:numCU_wdwise-1

        GRH_const_x(:,i+1) = [grh_width*i,grh_width*(i+1),grh_width*i,grh_width*(i+1),...
            grh_width*i,grh_width*(i+1),grh_width*i,grh_width*(i+1),...
            (grh_width/2)*(2*i+1),(grh_width/2)*(2*i+1)]';
        
        GRH_const_y = [0,0,grh_length,grh_length,...
            0,0,grh_length,grh_length,...
            0,grh_length]';
        
        GRH_const_z = [0,0,0,0,...
            gut_hght,gut_hght,gut_hght,gut_hght,...
            (gut_hght + (grh_width*tand(r_slope)/2)),(gut_hght + (grh_width*tand(r_slope)/2))]';
        
        GRH_const_points{:,i+1} = table(GRH_const_x(:,i+1),GRH_const_y,GRH_const_z,'VariableNames',["x","y","z"],...
            'RowNames',["Point 1","Point 2","Point 3","Point 4","Point 5","Point 6","Point 7","Point 8","Point 9","Point 10"]);
    end

elseif numCU_lnwise > 0 && (numCU_wdwise == 0 || numCU_wdwise == 1)
    for i=0:numCU_lnwise-1

        GRH_const_x = [0,grh_width,0,grh_width,...
            0,grh_width,0,grh_width,...
            grh_width/2,grh_width/2]';
        
        GRH_const_y(:,i+1) = [grh_length*i,grh_length*i,grh_length*(i+1),grh_length*(i+1),...
            grh_length*i,grh_length*i,grh_length*(i+1),grh_length*(i+1),...
            grh_length*i,grh_length*(i+1)]';
        
        GRH_const_z = [0,0,0,0,...
            gut_hght,gut_hght,gut_hght,gut_hght,...
            (gut_hght + (grh_width*tand(r_slope)/2)),(gut_hght + (grh_width*tand(r_slope)/2))]';
        
        GRH_const_points{i+1,:} = table(GRH_const_x,GRH_const_y(:,i+1),GRH_const_z,'VariableNames',["x","y","z"],...
            'RowNames',["Point 1","Point 2","Point 3","Point 4","Point 5","Point 6","Point 7","Point 8","Point 9","Point 10"]);
    end

elseif numCU_wdwise > 1 && numCU_lnwise > 1
    for i=0:numCU_wdwise-1
        for j=0:numCU_lnwise-1
            
            GRH_const_x(:,i+1) = [grh_width*i,grh_width*(i+1),grh_width*i,grh_width*(i+1),...
                grh_width*i,grh_width*(i+1),grh_width*i,grh_width*(i+1),...
                (grh_width/2)*(2*i+1),(grh_width/2)*(2*i+1)]';
            
            GRH_const_y(:,j+1) = [grh_length*j,grh_length*j,grh_length*(j+1),grh_length*(j+1),...
                grh_length*j,grh_length*j,grh_length*(j+1),grh_length*(j+1),...
                grh_length*j,grh_length*(j+1)]';
            
            GRH_const_z = [0,0,0,0,...
                gut_hght,gut_hght,gut_hght,gut_hght,...
                (gut_hght + (grh_width*tand(r_slope)/2)),(gut_hght + (grh_width*tand(r_slope)/2))]';
            
            GRH_const_points{j+1,i+1} = table(GRH_const_x(:,i+1),GRH_const_y(:,j+1),GRH_const_z,'VariableNames',["x","y","z"],...
                'RowNames',["Point 1","Point 2","Point 3","Point 4","Point 5","Point 6","Point 7","Point 8","Point 9","Point 10"]);
        end
    end
    
    if roof_side == "South"
            PV_points_x = [grh_width-PV_width*cosd(r_slope), GRH_const_points.x("Point 9"),grh_width-PV_width*cosd(r_slope),GRH_const_points.x("Point 9")];
    
            PV_points_y = [0,GRH_const_points.y("Point 9"),...
                PV_length,PV_length];
    
            PV_points_z = [PV_width*tand(r_slope)*cosd(r_slope) + gut_hght, GRH_const_points.z("Point 9"),...
                PV_width*tand(r_slope)*cosd(r_slope) + gut_hght, GRH_const_points.z("Point 9")];
            
            PV_points = table(PV_points_x,PV_points_y,PV_points_z,'VariableNames',["x","y","z"],...
                'RowNames',["Point 1","Point 2","Point 3","Point 4"]);
        end
    end
end


