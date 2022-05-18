function [data_mins,data_maxes,data_norm] = normlzn(inputs,rmax,rmin)
    
    data_mins = min(inputs,[],1);
    data_maxes = max(inputs,[],1);
    
    data_norm = ((rmax - rmin)*((inputs - data_mins)./(data_maxes - data_mins))) + rmin;

end
