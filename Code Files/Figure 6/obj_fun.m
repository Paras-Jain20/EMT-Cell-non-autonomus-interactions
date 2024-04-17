function cost = obj_fun(transform_para,data_type, data,pop_model_indx, num_rep,timept,interpolated_timept,interpolation,original_num_timepts)

temp_ode_data = ODE_simulation(pop_model_indx, num_rep,sort(unique(timept)),interpolated_timept,interpolation,data,original_num_timepts,transform_para);

switch pop_model_indx

    case {1,2,3,4,5,6,7,8,9,10,11,12,13,14}
        ode_sim_data = ODE_simulation(pop_model_indx, num_rep,timept,interpolated_timept,interpolation,data,transform_para);

        if(length(timept)>5)
            temp_ode_data = ODE_simulation(pop_model_indx, num_rep,timept,interpolated_timept,interpolation,data,transform_para);
            ode_sim_data = [ode_sim_data; temp_ode_data];
        end

        if( (pop_model_indx == 1 && data_type == 2)) % pop model 1 simulated cell numbers instead of cell fraction; so if we comparing it with Exp data then we should change ODE soultion to fraction;
            ode_sim_data = ode_sim_data./sum(ode_sim_data,2);
        end

        cost = (data(:,2) - ode_sim_data(:,1))./data(:,4); % cost to fit Mesenchymal fraction dynamics across all replicate and time points

    case {15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27}
                
        % to separate the dynamics of cell number to calculate Venus labelled EpCAM low fraction from the
        % ODE simulations
        temp_ode_data_1 = temp_ode_data;
        num_additional_timepts = length(timept) - original_num_timepts;
        while(num_additional_timepts > 0)
            if(sum(timept(original_num_timepts + num_additional_timepts) ~= timept(1:original_num_timepts)) == original_num_timepts)
                indx_1 = temp_ode_data_1(:,1) ~= timept(original_num_timepts + num_additional_timepts);
                temp_ode_data_1 = temp_ode_data_1(indx_1,:);
            end
            num_additional_timepts = num_additional_timepts - 1;
        end
        temp_ode_data_1 = temp_ode_data_1(:,2:end); % removing time column

        % to separate the dynamics of cell number to calculate M cell population from the
        % ODE simulations
        temp_ode_data_2 = temp_ode_data;
        num_additional_timepts = original_num_timepts;
        while(num_additional_timepts > 0)
            if(sum(timept(num_additional_timepts) ~= timept(original_num_timepts+1:end)) == (length(timept) - original_num_timepts))
            indx_1 = temp_ode_data_2(:,1) ~= timept(num_additional_timepts);
            temp_ode_data_2 = temp_ode_data_2(indx_1,:);
            end
            num_additional_timepts = num_additional_timepts - 1;
        end
        temp_ode_data_2 = temp_ode_data_2(:,2:end); % removing time column

        % ode_sim_data = [temp_ode_data_1(:,3)./(sum(temp_ode_data_1(:,3:4),2)); sum(temp_ode_data_2(:,[1 3]),2)]; % concatinating two output variable - Venus labelled EpCAM low fraction and total EpCAM low population

        ode_sim_data = [temp_ode_data_1(:,3)./(sum(temp_ode_data_1(:,3:4),2)); (sum(temp_ode_data_2(:,1:4),2))]; % concatinating two output variable - Venus labelled EpCAM low fraction and total cell number

        % sim_data = ode_sim_data(:,3)./(sum(ode_sim_data(:,3:4),2)); % to find the fraction of Venus Label EpCAM -ve cells in the total fraction of Venus labelled cells
        cost = (data(:,2) - ode_sim_data)./data(:,4); % cost to fit Mesenchymal fraction dynamics across all replicate and time points
end
end