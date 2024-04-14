function cost = obj_fun(transform_para,data_type, data,pop_model_indx, num_rep,timept,interpolated_timept,init_cond, interpolation)


ode_sim_data = ODE_simulation(pop_model_indx, num_rep,timept,interpolated_timept,interpolation,data, init_cond, transform_para);

[~, msgid] = lastwarn;

if strcmp(msgid,'MATLAB:ode45:IntegrationTolNotMet')
    tol_issue = true;
else
    tol_issue = false;
end

lastwarn('');

switch pop_model_indx

    case {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15}

        if( (pop_model_indx == 1 && data_type == 2)) % pop model 1 simulated cell numbers instead of cell fraction; so if we comparing it with Exp data then we should change ODE soultion to fraction;
            ode_sim_data = ode_sim_data./sum(ode_sim_data,2);
        end

        if(length(data(:,2)) ~= length(ode_sim_data(:,1)) || tol_issue)
            cost = zeros(size(data(:,2)));
        else
            cost = (data(:,2) - ode_sim_data(:,1))./data(:,4); % cost to fit Mesenchymal fraction dynamics across all replicate and time points
        end

    case {16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 29,30,31, 32, 33}

        sim_data = ode_sim_data(:,3)./(sum(ode_sim_data(:,3:4),2)); % to find the fraction of Venus Label EpCAM -ve cells in the total fraction of Venus labelled cells
        
        if(length(data(:,2)) ~= length(sim_data) || tol_issue)
            cost = zeros(size(data(:,2)));
        else
            cost = (data(:,2) - sim_data)./data(:,4); % cost to fit Mesenchymal fraction dynamics across all replicate and time points
        end
end
end