%function that would take profile likelihood data as input and checks for
%practical identifiability of the parameters


function practical_identifiability_status = practical_identifiability_checker(pop_model_indx,profile_para_data,practical_identifiability_status,r2)
    
    
    for para_indx = (unique(profile_para_data(:,end)))'
    
    x_rate = profile_para_data(profile_para_data(:,end) == para_indx, 1:end-2);
    y = profile_para_data(profile_para_data(:,end) == para_indx, end-1);
        switch pop_model_indx
        
            case {1,2,5,8,9}
        
                if(para_indx <=4) % i.e. only for pop models with dimensional parameters
                    x_time = 1./x_rate;
                else
                    x_time = x_rate;
                end
        
            case{3,4,6,7,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27}
       
                x_time = x_rate; 
        
        end
        % based on the choosen alpha = 0.954, define the delta chi sqr avlue
        
        delta_chi_sq_theta = 4;
        
        if(isempty(r2))
            switch pop_model_indx
                case{6} % i.e. population model with 1 parameter
                    delta_chi_sq_all = 4;
                case {1,2,10,11,12,20} % i.e. population model with 4 parameters
                    delta_chi_sq_all = 9.70;
                case {3,7,18} % i.e. population model with 3 parameters
                    delta_chi_sq_all = 8.02;
                case {4,5,13,14,15,16,17,19,29} % i.e. population model with 5 parameters
                    delta_chi_sq_all = 11.3;
                case {8,9,22,23,24} % i.e. population model with 6 parameters
                    delta_chi_sq_all = 12.8;
                case {21} % i.e. population model with 7 parameters
                    delta_chi_sq_all = 14.067;
            end
        else
            switch pop_model_indx
                case{6} % i.e. population model with 1 parameter
                    delta_chi_sq_all = 4;
                case{3} % i.e. population model with 2 parameter
                    delta_chi_sq_all = 6.17;
                case {4,5,13,14,15,16,17,19,29} % i.e. population model with 4 parameters
                    delta_chi_sq_all = 9.70;
                case {1,2,10,11,12,20} % i.e. population model with 3 parameters
                    delta_chi_sq_all = 8.02;
                case {8,9,22,23,24} % i.e. population model with 5 parameters
                    delta_chi_sq_all = 11.3;
                case {21} % i.e. population model with 6 parameters
                    delta_chi_sq_all = 12.8;
            end
        end
    
%         y = y*(num_rep*length(timept)*num_init_cond); % rescaling chi-sqr levels
        [minima,posn] = min(y); %finds the minimum in y and the corresponding index
        minima_position = x_time(posn);
        
        chi_sq_threshold_theta = minima + delta_chi_sq_theta;% the chi_sq threshold for single parameter
        chi_sq_threshold_all = minima + delta_chi_sq_all;%the chi_sq threshold for all params
        
        %finding the points at which goodness 'y' intersects the
        %chi_sq-threshold
        intersections = find((y(1:end-1) - chi_sq_threshold_all).*(y(2:end) - chi_sq_threshold_all) < 0);
        
        if (length(intersections) == 2)
            practical_identifiability_status(para_indx) = 1;
        end

    
    end

    if(~isempty(r2))
    practical_identifiability_status(1) = 1;
    end
end
