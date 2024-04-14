%%%%%%%%% main_function %%%%%%%%%
close all

% whether we are using the test or experimental data: test_data = 1; exp_data = 2
data_type = 2;

test_data_indx = 1; % necessary dummy variable to just to fill in variable list in profile likelihood function
noise_factor = 1; % necessary dummy variable to just to fill in variable list in profile likelihood function
r1 = 1/54; % M cells division rates; scaling factor for non dimensional population model
inverse_r2 =  35;% r_e
r2 = 1/(inverse_r2*r1); % normalized E cell division rate wrt M cell; if r2 is empty set then it is a variable parameter

default_total_itr = 3;% number of iterative parameter search to be performed

pop_model = 15; % population_model; select a population growth model
for pop_model_indx = 1:max(pop_model)
    if(~isempty(find(pop_model_indx == pop_model,1)))
        for g_t_ratio = 250 % growth_transition_ratio; max ratio between growth and transition
            close all

            % type of optimum parameters search to be performed:
            % use existing parameters = choose 0
            % optimizer based : choose 1
            % uniform sampling based : choose 2; this method is not updated so might not work
            search_method = 1; % parameter search method

            if(pop_model_indx < 15) % for models capturing dynamics of Bhatia et al. data

                cd('C:\Users\user\OneDrive - Indian Institute of Science\GitHub\EMT-State-Transition-2.0\Parameters sets from experimental data')

                if(isfile(['para_Bhatia_data_pd_mod_' num2str(pop_model_indx) '_s_m_' num2str(search_method) '_b_t_ratio_' num2str(g_t_ratio) '_inverse_r1_' num2str(1/r1) '.csv']))

                    para_data = table2array(readtable(['para_Bhatia_data_pd_mod_' num2str(pop_model_indx) '_s_m_' num2str(search_method) '_b_t_ratio_' num2str(g_t_ratio) '_inverse_r1_' num2str(1/r1) '.csv']));
                    num_para_set = size(para_data,1);
                    if(num_para_set >= default_total_itr)
                        select_para = (para_data(:,1:end-1));
                        goodness = (para_data(:,end));
                        para_opt = false; % to perform parameter search using optimizer or not
                        total_itr = 0;
                    else

                        para_opt = true;
                        total_itr = default_total_itr-num_para_set; % max 3 iterations have to be performed
                    end
                else
                    para_opt = true; % parameter_optimization; if the data file is not present then perform the parameter search
                    total_itr = default_total_itr;
                end

                initial_dist = {'high', 'low'; 'low' , 'high'}; % for Bhatia et al data; Mes:Epi

                data = [];

                cd('C:\Users\user\OneDrive - Indian Institute of Science\Projects\State transition\Experimental data');

                % importing Bhatia et al data
                for ini_dist_indx = 1:2

                    temp_data = (readmatrix((['Bhatia et al 2019\M_' initial_dist{ini_dist_indx,1} '_E_' initial_dist{ini_dist_indx,2} '_dynamics.xlsx'])));

                    temp_data = [temp_data(:,1) temp_data(:,2:end)/100]; % changing population percentage into fractions
                    data = [data;temp_data];
                end
                % data charactersitics
                num_rep = 3;% number_of_replicates
                num_init_cond = 2; % number of different initial condition in the data
                init_cond = 1:num_init_cond;% initial_condition_index
                timept = (2:2:8)*(7*24)*r1; % time is unitless as it is scaled by m cells division rate (r1)
                interpolated_timept = r1:r1:timept(end); % interpolating time to get more data points to plot ODE solution trajectories

            else
                cd('C:\Users\user\OneDrive - Indian Institute of Science\GitHub\EMT-State-Transition-2.0\Parameters sets from experimental data')
                if(isfile(['para_Yamamoto_data_pd_mod_' num2str(pop_model_indx) '_s_m_' num2str(search_method) '_b_t_ratio_' num2str(g_t_ratio) '_inverse_r1_' num2str(1/r1) '_inverse_r2_' num2str(inverse_r2) '.csv']))

                    para_data = table2array(readtable(['para_Yamamoto_data_pd_mod_' num2str(pop_model_indx) '_s_m_' num2str(search_method) '_b_t_ratio_' num2str(g_t_ratio) '_inverse_r1_' num2str(1/r1) '_inverse_r2_' num2str(inverse_r2) '.csv']));

                    num_para_set = size(para_data,1);
                    if(num_para_set >= default_total_itr)
                        select_para = (para_data(:,1:end-1));
                        goodness = (para_data(:,end));
                        para_opt = false; % to perform parameter search using optimizer
                        total_itr = 0;
                    else
                        para_opt = true;
                        total_itr = default_total_itr-num_para_set;
                    end

                else
                    para_opt = true; % if the data file is not present then perform the para search
                    total_itr = default_total_itr;

                end

                %               Importing Yamamoto et al. 2017 data; models having
                %               pop_model_indx >= 15
                cd('C:\Users\user\OneDrive - Indian Institute of Science\Projects\State transition\Experimental data\Yamamoto et al. 2017')
                data = (readmatrix('Yamamoto et al 2017 Data 1 cell fraction.xlsx'));
                num_rep = 3;% number_of_replicates
                num_init_cond = 5; % number of different initial condition in the data
                init_cond = 1:num_init_cond; % initial_condition_index
                timept = (0:3:12)*24*r1; % time is unitless as it is scaled by m cells division rate (r1)

                %                 cd('C:\Users\user\OneDrive - Indian Institute of Science\GitHub\EMT-State-Transition-2.0\Codes files ram\data post sampling')
                %                 data = (readmatrix('data_post_iterative_sampling_Yamamoto_Exp_data_pd_mod_21_b_t_ratio_250_inverse_r1_54_inverse_r2_35'));
                %                 num_rep = 3;
                %                 num_init_cond = 5; % number of different initial condition in the data
                %                 init_cond = 1:num_init_cond;
                %                 timept = data(1:(size(data,1)/(num_rep*num_init_cond)),1)'; % time is unitless as it is scaled by m cells division rate (r1)

                interpolated_timept = 0:r1:timept(end);  % interpolating time to get more data points to plot ODE solution trajectories
            end

            num_time_pts = length(timept); % number_of_time_points
            %             cd('C:\Users\user\OneDrive - Indian Institute of Science\GitHub\EMT-State-Transition-2.0\Codes files')
            cd('C:\Users\user\OneDrive - Indian Institute of Science\GitHub\EMT-State-Transition-2.0\Manuscript documents\Final Codes files')
            %                             plot_func(select_para,goodness,data_type, data,pop_model_indx, num_rep,timept, interpolated_timept, interpolation,r1, init_cond, g_t_ratio);

            if(para_opt)

                interpolation = 0; % to add more time points for the ODE simulated data; must be zero if we using objective function to minimize error with the data

                % parameter search function
                [select_para,goodness] = para_search(search_method,pop_model_indx,total_itr, data_type, data, num_rep, init_cond, timept, interpolated_timept, interpolation, num_time_pts,g_t_ratio,r2); % select_para are filtered for re > rm

                interpolation = 1; % to add more time points for the ODE simulated data; must be zero if we using objective function to minimize error with the data
                % plot data and simulated trajectories for the best fit parameter set
                plot_func(select_para,goodness,data_type, data,pop_model_indx, num_rep,timept, interpolated_timept, interpolation,r1,r2,inverse_r2, init_cond, g_t_ratio);


                cd('C:\Users\user\OneDrive - Indian Institute of Science\GitHub\EMT-State-Transition-2.0\Parameters sets from experimental data')

                if(pop_model_indx <15 )
                    writematrix(['para_Bhatia_data_pd_mod_' num2str(pop_model_indx) '_s_m_' num2str(search_method) '_b_t_ratio_' num2str(g_t_ratio) '_inverse_r1_' num2str(1/r1) '.csv'],[select_para,goodness],'-append');
                else
                    writematrix([select_para,goodness],['para_Yamamoto_data_pd_mod_' num2str(pop_model_indx) '_s_m_' num2str(search_method) '_b_t_ratio_' num2str(g_t_ratio) '_inverse_r1_' num2str(1/r1) '_inverse_r2_' num2str(inverse_r2) '.csv']);
                end
            end

            % Estimation of confidence interval using profile likelihood
            interpolation = 0; % to add more time points for the ODE simulated data; must be zero if we using objective function to minimize error with the data
            profile_sample_size = 2000; % total number of samples between lower and upper bound of the parameter

            cd('C:\Users\user\OneDrive - Indian Institute of Science\GitHub\EMT-State-Transition-2.0\Parameters sets from experimental data')

            %             if((pop_model_indx >=19 && (~isfile(['profile_lik_para_Yamamoto_Exp_data_pd_mod_' num2str(pop_model_indx) '_b_t_ratio_' num2str(g_t_ratio) '_inverse_r1_' num2str(1/r1) '_para_steps_' num2str(profile_sample_size) '.csv'])))...
            %                     || (pop_model_indx <= 14 && ~isfile(['profile_lik_para_Bhatia_Exp_data_pd_mod_' num2str(pop_model_indx) '_b_t_ratio_' num2str(g_t_ratio) '_inverse_r1_' num2str(1/r1)  '_para_steps_' num2str(profile_sample_size) '.csv'])) )
            if(true)
                cd('C:\Users\user\OneDrive - Indian Institute of Science\GitHub\EMT-State-Transition-2.0\Manuscript documents\Final Codes files')

                % params: profile likelihood of the model parameters;
                % chi_sqr_err: chi-square values
                [params, chi_sqr_err] = profile_likelihood_est(search_method,pop_model_indx,total_itr,data_type, data, num_rep, init_cond,timept, interpolated_timept, interpolation, num_time_pts,g_t_ratio,profile_sample_size, noise_factor, test_data_indx, r1, r2);

                % saving the parameters and chi-sqr values obtained from the
                % analysis in a single matrix

                para_data = [];

                for i = 1:length(params) % since params is the cell vector
                    if(~isempty(params{i}))
                        para_data =  [para_data; params{i}(:,1:end-1) chi_sqr_err{i}' params{i}(:,end)];
                    end
                end

                %                 cd('C:\Users\user\OneDrive - Indian Institute of Science\GitHub\EMT-State-Transition-2.0\Parameters sets from experimental data')

                if(isempty(r2))
                    if(pop_model_indx <15 && pop_model_indx == 4 )
                        writematrix(para_data, ['profile_lik_para_Bhatia_Exp_data_pd_mod_' num2str(pop_model_indx) '_b_t_ratio_' num2str(g_t_ratio) '_inverse_r1_' num2str(1/r1)  '_para_steps_' num2str(profile_sample_size) '_negative_comp.csv']);

                    elseif(pop_model_indx <15 && pop_model_indx ~= 4)
                        writematrix(para_data, ['profile_lik_para_Bhatia_Exp_data_pd_mod_' num2str(pop_model_indx) '_b_t_ratio_' num2str(g_t_ratio) '_inverse_r1_' num2str(1/r1)  '_para_steps_' num2str(profile_sample_size) '.csv']);

                    elseif(pop_model_indx > 15)
                        writematrix(para_data,['profile_lik_para_Yamamoto_Exp_data_pd_mod_' num2str(pop_model_indx) '_b_t_ratio_' num2str(g_t_ratio) '_inverse_r1_' num2str(1/r1) '_para_steps_' num2str(profile_sample_size)  '.csv']);

                    end
                else
                    if(pop_model_indx <15 )
                        writematrix(para_data,['profile_lik_para_Bhatia_Exp_data_pd_mod_' num2str(pop_model_indx) '_b_t_ratio_' num2str(g_t_ratio) '_inverse_r1_' num2str(1/r1) '_inverse_r2_' num2str(inverse_r2)  '_para_steps_' num2str(profile_sample_size) '.csv']);
                    else
                        writematrix(para_data,['profile_lik_para_Yamamoto_Exp_data_pd_mod_' num2str(pop_model_indx) '_b_t_ratio_' num2str(g_t_ratio) '_inverse_r1_' num2str(1/r1) '_inverse_r2_' num2str(inverse_r2) '_para_steps_' num2str(profile_sample_size)  '.csv']);
                        %                         writematrix(para_data,['profile_lik_para_post_itr_samp_pd_mod_' num2str(pop_model_indx) '_b_t_ratio_' num2str(g_t_ratio) '_inverse_r1_' num2str(1/r1) '_inverse_r2_' num2str(inverse_r2) '_para_steps_' num2str(profile_sample_size)  '.csv']);
                    end
                end
            end
        end
    end
end