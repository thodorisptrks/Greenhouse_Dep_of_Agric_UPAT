function [days_array] = day_of_year(year)
    
    if leapyear(year) == 0
        days_array = (1:1:365)';
    else
        days_array = (1:1:366)';
    end

end