% Greenhouse coordinates calculations
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