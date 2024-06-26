function joint_ode_sim_data = ODE_simulation(pop_model_indx, num_rep,timept, interpolated_timept,interpolation, data,original_num_timepts,varargin)

joint_ode_sim_data = [];

handles = feval(@model_ODE);

if(pop_model_indx > 14 && pop_model_indx < 28)
    num_ini_dist = 5; % 5 different initial population distribution in the Yamamoto study
else
    num_ini_dist = 2; % 2 different initial population distribution in the Bhatia study
end

for ini_dist_indx = 1:num_ini_dist
    for rep_indx = 1:num_rep
        switch pop_model_indx

            case 1 % pop_model with cell number dynamics

                if(ini_dist_indx == 1)
                    initial_pop = [100000, 100];
                else
                    initial_pop = [100, 100000];
                end

            case {2,3,4,5,6,7,8,9,10,11,12,13,14}

                if(ini_dist_indx == 1)
                    initial_pop = [0.999, 0.001];
                else
                    initial_pop = [0.001, 0.999];
                end

            case {15, 16,17,18, 19, 20, 21, 22, 23, 24, 25, 26, 27} % yamamoto Venus labelled - non-labelled cells; Here, the fractional distribution at the zeroth time point is given in the data
                if(ini_dist_indx == 1)

                    Ep_neg_vl_frac = 1*data((ini_dist_indx-1)*length(timept(1:original_num_timepts))*num_rep +(rep_indx-1)*length(timept(1:original_num_timepts))+1,2);
                    Ep_pos_vl_frac = 1*data((ini_dist_indx-1)*length(timept(1:original_num_timepts))*num_rep +(rep_indx-1)*length(timept(1:original_num_timepts))+1,3);
                    initial_pop = [0, 0, Ep_neg_vl_frac , Ep_pos_vl_frac]; % Venus -ve EpCAM -ve, Venus -ve EpCAM +ve, Venus +ve EpCAM -ve, Venus +ve EpCAM +ve

                elseif(ini_dist_indx == 2)

                    Ep_neg_vl_frac = 0.9*data((ini_dist_indx-1)*length(timept(1:original_num_timepts))*num_rep +(rep_indx-1)*length(timept(1:original_num_timepts))+1,2);
                    Ep_pos_vl_frac = 0.9*data((ini_dist_indx-1)*length(timept(1:original_num_timepts))*num_rep +(rep_indx-1)*length(timept(1:original_num_timepts))+1,3);
                    initial_pop = [0.1, 0, Ep_neg_vl_frac , Ep_pos_vl_frac]; % Venus -ve EpCAM -ve, Venus -ve EpCAM +ve, Venus +ve EpCAM -ve, Venus +ve EpCAM +ve

                elseif(ini_dist_indx == 3)

                    Ep_neg_vl_frac = 0.1*data((ini_dist_indx-1)*length(timept(1:original_num_timepts))*num_rep +(rep_indx-1)*length(timept(1:original_num_timepts))+1,2);
                    Ep_pos_vl_frac = 0.1*data((ini_dist_indx-1)*length(timept(1:original_num_timepts))*num_rep +(rep_indx-1)*length(timept(1:original_num_timepts))+1,3);
                    initial_pop = [0.9, 0, Ep_neg_vl_frac , Ep_pos_vl_frac]; % Venus -ve EpCAM -ve, Venus -ve EpCAM +ve, Venus +ve EpCAM -ve, Venus +ve EpCAM +ve

                elseif(ini_dist_indx == 4)

                    Ep_neg_vl_frac = 0.1*data((ini_dist_indx-1)*length(timept(1:original_num_timepts))*num_rep +(rep_indx-1)*length(timept(1:original_num_timepts))+1,2);
                    Ep_pos_vl_frac = 0.1*data((ini_dist_indx-1)*length(timept(1:original_num_timepts))*num_rep +(rep_indx-1)*length(timept(1:original_num_timepts))+1,3);
                    initial_pop = [0, 0.9, Ep_neg_vl_frac , Ep_pos_vl_frac]; % Venus -ve EpCAM -ve, Venus -ve EpCAM +ve, Venus +ve EpCAM -ve, Venus +ve EpCAM +ve

                else

                    Ep_neg_vl_frac = 0.1*data((ini_dist_indx-1)*length(timept(1:original_num_timepts))*num_rep +(rep_indx-1)*length(timept(1:original_num_timepts))+1,2);
                    Ep_pos_vl_frac = 0.1*data((ini_dist_indx-1)*length(timept(1:original_num_timepts))*num_rep +(rep_indx-1)*length(timept(1:original_num_timepts))+1,3);
                    initial_pop = [0.9, 0, Ep_neg_vl_frac , Ep_pos_vl_frac]; % Venus -ve EpCAM -ve, Venus -ve EpCAM +ve, Venus +ve EpCAM -ve, Venus +ve EpCAM +ve

                end
        end


        switch pop_model_indx

            case {1,2,3,4,5,6,7,8,9,10,11,12,13,14}
                if(interpolation == 0)
                    [t,ode_sim_data] = ode45(@(t,y) handles{pop_model_indx}(t,y,varargin{1}), [0 timept],initial_pop); % Zeroth time point is not present in the Bhatia et al data
                else
                    [t,ode_sim_data] = ode45(@(t,y) handles{pop_model_indx}(t,y,varargin{1}), [0 interpolated_timept],initial_pop); % Zeroth time point is not present in the Bhatia et al data
                end
                joint_ode_sim_data = [joint_ode_sim_data;ode_sim_data(2:end,:)]; % we are excluding the zeroth time point solution value to compare the filtered solution to the experimental data
            case {15, 16 ,17,18}
                if(interpolation == 0)
                    [t,ode_sim_data] = ode45(@(t,y) handles{pop_model_indx}(t,y,varargin{1}), [timept],initial_pop); % Zeroth time point is present in the Yamamoto et al data therefore we donot include it explicitly
                else
                    [t,ode_sim_data] = ode45(@(t,y) handles{pop_model_indx}(t,y,varargin{1}), [interpolated_timept],initial_pop);
                end
                joint_ode_sim_data = [joint_ode_sim_data;ode_sim_data(1:end,:)];
            case{19, 20, 21, 22, 23, 24, 25, 26, 27}
                initial_pop = 1000000 * initial_pop;
                    if(interpolation == 0)
                        [t,ode_sim_data] = ode45(@(t,y) handles{pop_model_indx}(t,y,varargin{1}), [timept],initial_pop); % Zeroth time point is present in the Yamamoto et al data therefore we donot include it explicitly
                    else
                        [t,ode_sim_data] = ode45(@(t,y) handles{pop_model_indx}(t,y,varargin{1}), [interpolated_timept],initial_pop); % we start with 1 million cells in the culture.
                    end
                    joint_ode_sim_data = [joint_ode_sim_data;t ode_sim_data(1:end,:)];
        end

    end
end

end