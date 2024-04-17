%the iterative sampling code
clear all
close all

%some relevant parameters that go into the profile likeli function
profile_sample_size = 2000;
data_type = 2;
test_data_indx = 1; % necessary dummy variable to just to fill in variable list in profile likelihood function
noise_factor = 1; % necessary dummy variable to just to fill in variable list in profile likelihood function
r1 = 1/54; % m cells division rates; scaling factor for non dimensional population model
inverse_r2 = 35;
r2 = 1/(inverse_r2*r1);
num_rep = 3;
search_method = 1;

num_para_steps = 100;% argument in the iterative improvement para bounds function

%relevant directories

% code_file_path = 'D:\Ramanarayanan\Documents\BS-MS\3rd yr\SEM2\mohit jolly\optimisation\Results in Aug and Sep 23\code files ram';
% para_set_from_exp_path = 'D:\Ramanarayanan\Documents\BS-MS\3rd yr\SEM2\mohit jolly\optimisation\Results in Aug and Sep 23\Parameter sets from experimental data';
% experimental_data_path = 'D:\Ramanarayanan\Documents\BS-MS\3rd yr\SEM2\mohit jolly\optimisation\Results in July 23\Experimental data';
% profile_liklihood_path = 'D:\Ramanarayanan\Documents\BS-MS\3rd yr\SEM2\mohit jolly\optimisation\Results in Aug and Sep 23\Parameter sets from experimental data';
% iterative_sampling_results_path = 'D:\Ramanarayanan\Documents\BS-MS\3rd yr\SEM2\mohit jolly\optimisation\Results in Aug and Sep 23\Results';

code_file_path = 'C:\Users\user\OneDrive - Indian Institute of Science\GitHub\EMT-State-Transition-2.0\code files ram';
para_set_from_exp_path = 'C:\Users\user\OneDrive - Indian Institute of Science\GitHub\EMT-State-Transition-2.0\Parameters sets from experimental data';
experimental_data_path = 'C:\Users\user\OneDrive - Indian Institute of Science\Projects\State transition\Experimental data';
profile_liklihood_path = 'C:\Users\user\OneDrive - Indian Institute of Science\GitHub\EMT-State-Transition-2.0\Parameters sets from experimental data';
% iterative_sampling_results_path = 'C:\Users\user\OneDrive - Indian Institute of Science\GitHub\EMT-State-Transition-2.0\code files ram\data post sampling';
iterative_sampling_results_path = 'C:\Users\user\OneDrive - Indian Institute of Science\GitHub\EMT-State-Transition-2.0\Manuscript documents\Revision\Parameter identifiability improvement';

% select a population growth model
for pop_model_indx = 21%[20 21 26 27 25 23]%[1 2 3 10 11 12 13]
    for birth_tran_ratio = 250%[5 10 20 50 100]%[50 100 150 200 250 300] %200 500] %[5 10 20 50 100] % max ratio between birth and transition
        cd(profile_liklihood_path)
        %         if((pop_model_indx < 15 && ~isfile(['data_post_iterative_sampling_Bhatia_Exp_data_pd_mod_' num2str(pop_model_indx) '_b_t_ratio_' num2str(birth_tran_ratio) '_inverse_r1_'   '.csv'])) ...
        %                 || (pop_model_indx > 15 &&  ~isfile(['data_post_iterative_sampling_Yamamoto_Exp_data_pd_mod_' num2str(pop_model_indx) '_b_t_ratio_' num2str(birth_tran_ratio) '_inverse_r1_' num2str(1/r1) '.csv']))...
        %                 || (pop_model_indx == 4 && ~isfile(['data_post_iterative_sampling_Bhatia_Exp_data_pd_mod_' num2str(pop_model_indx) '_b_t_ratio_' num2str(birth_tran_ratio) '_inverse_r1_' num2str(1/r1) '_negative_comp.csv'])))
        if(true)
            %fetching the relevant experimental data
            if pop_model_indx < 15

                %getting the relevant profile likelihood data
                cd(profile_liklihood_path)
                if(pop_model_indx == 4)
                    profile_likli_data = readmatrix(['profile_lik_para_Bhatia_Exp_data_pd_mod_' num2str(pop_model_indx) '_b_t_ratio_' num2str(birth_tran_ratio) '_inverse_r1_' num2str(1/r1)  '_para_steps_' num2str(profile_sample_size) '_negative_comp.csv']);

                    %                     profile_likli_data = readmatrix(['profile_lik_para_Bhatia_Exp_data_pd_mod_' num2str(pop_model_indx) '_b_t_ratio_' num2str(birth_tran_ratio) '_inverse_r1_' num2str(1/r1) '_dynamic_para_steps_with_default_value_' num2str(profile_sample_size) '_negative_comp.csv']);
                else
                    %                     profile_likli_data = readmatrix(['profile_lik_para_pd_mod_' num2str(pop_model_indx) '_b_t_ratio_' num2str(birth_tran_ratio) '_inverse_r1_' num2str(1/r1) '_inverse_r2_' num2str(inverse_r2) '_dynamic_para_steps_with_default_value_' num2str(profile_sample_size) '.csv']);
                    profile_likli_data = readmatrix(['profile_lik_para_Bhatia_Exp_data_pd_mod_' num2str(pop_model_indx) '_b_t_ratio_' num2str(birth_tran_ratio) '_inverse_r1_' num2str(1/r1) '_para_steps_' num2str(profile_sample_size) '.csv']);
                end

                initial_dist = {'high', 'low'; 'low' , 'high'}; % for Bhatia et al data; Mes:Epi

                data = [];

                cd(experimental_data_path);

                % importing Bhatia et al data
                for ini_dist_indx = 1:2

                    temp_data = readmatrix(['Bhatia et al 2019\M_' initial_dist{ini_dist_indx,1} '_E_' initial_dist{ini_dist_indx,2} '_dynamics.xlsx']);

                    temp_data = [temp_data(:,1)*(7*24)*r1 temp_data(:,2:end)/100]; % changing population percentage into fractions
                    data = [data;temp_data];
                end
                % data charactersitics
                num_rep = 3;
                num_init_cond = 2; % number of different initial condition in the data
                timept = (2:2:8)*(7*24)*r1; % time is unitless as it is scaled by m cells division rate (r1)
                interpolated_timept = (0.5:0.5:8)*(7*24)*r1; %r1:r1:timept(end); % interpolating time to get more data points to plot ODE solution trajectories

            else

                cd(profile_liklihood_path)
                if(isempty(r2))
                    profile_likli_data = readmatrix(['profile_lik_para_Yamamoto_Exp_data_pd_mod_' num2str(pop_model_indx) '_b_t_ratio_' num2str(birth_tran_ratio) '_inverse_r1_' num2str(1/r1)  '_para_steps_' num2str(profile_sample_size) '.csv']);
                else
                    profile_likli_data = readmatrix(['profile_lik_para_Yamamoto_Exp_data_pd_mod_' num2str(pop_model_indx) '_b_t_ratio_' num2str(birth_tran_ratio) '_inverse_r1_' num2str(1/r1) '_inverse_r2_' num2str(inverse_r2) '_para_steps_' num2str(profile_sample_size) '.csv']);

                end
                %                 profile_likli_data = readmatrix(['profile_lik_para_Yamamoto_Exp_data_pd_mod_' num2str(pop_model_indx) '_b_t_ratio_' num2str(birth_tran_ratio) '_inverse_r1_' num2str(1/r1) '_dynamic_para_steps_with_default_value_' num2str(profile_sample_size) '.csv']);
                %                 profile_likli_data = readmatrix(['profile_lik_para_Yamamoto_Exp_data_pd_mod_' num2str(pop_model_indx) '_b_t_ratio_' num2str(birth_tran_ratio) '_inverse_r1_' num2str(1/r1) '_para_steps_' num2str(1000) '.csv']);

                % Importing Yamamoto et al. 2017 data; models having
                % pop_model_indx >= 15
                cd([experimental_data_path '\Yamamoto et al. 2017'])
                data = table2array(readtable('Yamamoto et al 2017 Data 1 cell fraction.xlsx'));
                data(:,1) = data(:,1)*24*r1;
                num_rep = 3;
                num_init_cond = 5; % number of different initial condition in the data
                timept = (0:3:12)*24*r1; % time is unitless as it is scaled by m cells division rate (r1)
                original_num_timepts = length(timept); % number of the time points in the original data; a constant
                interpolated_timept = timept; %0:r1:timept(end);  % interpolating time to get more data points to plot ODE solution trajectories
            end

            current_profile_data = profile_likli_data;
            current_timept = timept;

            cd(code_file_path)
if(isempty(r2))
            initial_practical_identifiability_status = practical_identifiability_checker_modified(pop_model_indx,current_profile_data,zeros(1,length(unique(current_profile_data(:,end)))),r2);
else
            initial_practical_identifiability_status = practical_identifiability_checker_modified(pop_model_indx,current_profile_data,zeros(1,length(unique(current_profile_data(:,end)))+1),r2);
end
            initial_practical_identifiability_score = sum(initial_practical_identifiability_status)/length(initial_practical_identifiability_status);

            %the switch for the while loop
            practically_identifiable = false;
            sampling_iteration_threshold = 4;

            %initialisation of variables before the while loop
            sampling_iteration_indx = 1;
            current_profile_data = profile_likli_data;
            current_timept = timept;
            current_data = data;
            current_practical_identifiability_status = initial_practical_identifiability_status;
            current_practical_identifiability_score = initial_practical_identifiability_score;

            disp(['Pop model ' num2str(pop_model_indx) '; Prac Identifiability Status with original data ' num2str(current_practical_identifiability_status) ])

            while(~practically_identifiable && (sampling_iteration_indx <= sampling_iteration_threshold))

                %current_practical_identifiability_status = practical_identifiability_checker(pop_model_indx,current_profile_data,profile_sample_size,r1,num_rep,current_timept,num_init_cond);

                % A function that takes in the vector of identifiability
                %status and generates psuedo data for each of the parameters
                %and then returns the best psuedo data and reduced profile
                %likelihood corresponding to this data.


                [updated_data,updated_timept,updated_profile_liklihood,updated_practical_identifiability_status,terminate] = practical_identifiability_updater(pop_model_indx,current_practical_identifiability_status,current_profile_data,profile_sample_size, num_para_steps,current_data,num_rep,num_init_cond,current_timept,interpolated_timept,r1,r2,birth_tran_ratio,original_num_timepts);

                updated_practical_identifiability_status_score = sum(updated_practical_identifiability_status)/length(updated_practical_identifiability_status);
                %if all the parameters are identifiable, stop iteration.
                %else update the variables in the loop with the outputs
                %from the above function
                if(updated_practical_identifiability_status_score == 1)
                    practically_identifiable = true;
                end

                current_profile_data = updated_profile_liklihood;
                current_timept = updated_timept;
                current_data = updated_data;

                previous_iteration_practical_identifiability_status = current_practical_identifiability_status;
                current_practical_identifiability_status = updated_practical_identifiability_status;


                disp(['Pop model ' num2str(pop_model_indx) '; After ' num2str(sampling_iteration_indx) ' iteration Practical Identifiability Status is ' num2str(current_practical_identifiability_status) ])
                sampling_iteration_indx = sampling_iteration_indx+1;

                if terminate == 1
                    %                         practically_identifiable = true
                    disp(['Unable to add new data, stopping the simulation'])
                    break
                end
            end

            %generating the profile likelihood for the updated data that
            %came out of the iterative sampling

            interpolation = 0;
            num_time_pts = length(updated_timept);
            unidentifiable_para_only = 1;% 1 for generating profile liklihood of identifiable parameters which we left out in the last iteration above
            identifiable_para_prior_to_last_iteration = ~(previous_iteration_practical_identifiability_status);

            [params, chi_sqr_err] = profile_likelihood_est(pop_model_indx,data_type, updated_data, num_rep,num_init_cond,updated_timept, interpolated_timept, interpolation,birth_tran_ratio,profile_sample_size, noise_factor, test_data_indx, r1,unidentifiable_para_only,identifiable_para_prior_to_last_iteration,r2,original_num_timepts);


            % saving the parameters and chi-sqr values obtained from the
            % analysis in a single matrix
            para_data = [];

            for i = 1:length(params) % since params is the cell vector
                if(~isempty(params{i}))
                    para_data =  [para_data; params{i}(:,1:end-1) chi_sqr_err{i}' params{i}(:,end)];
                end
            end

            current_profile_data = [current_profile_data; para_data];


            %checking out the final status of all parameters post sampling
            practical_identifiability_status_post_iterative_sampling = practical_identifiability_checker(pop_model_indx,current_profile_data,current_practical_identifiability_status, r2);
            disp(['Final Practical Identifiability Status is ' num2str(practical_identifiability_status_post_iterative_sampling) ])

            %saving the updated data and updated profile likelihood data
            cd(iterative_sampling_results_path)
            if(isempty(r2))
                if(pop_model_indx  < 15)
                    if(pop_model_indx ~= 4)
                        writematrix(updated_data,['data_post_iterative_sampling_Bhatia_Exp_data_pd_mod_' num2str(pop_model_indx) '_b_t_ratio_' num2str(birth_tran_ratio) '_inverse_r1_' num2str(1/r1) '.csv']);
                        %                     writematrix(current_profile_data,['profile_lik_para_post_itr_samp_pd_mod_' num2str(pop_model_indx) '_b_t_ratio_' num2str(birth_tran_ratio) '_inverse_r1_' num2str(1/r1)  '_dynamic_para_steps_default_value_' num2str(profile_sample_size) '.csv']);
                        writematrix(current_profile_data,['profile_lik_para_post_itr_samp_pd_mod_' num2str(pop_model_indx) '_b_t_ratio_' num2str(birth_tran_ratio) '_inverse_r1_' num2str(1/r1)  '_para_steps_' num2str(profile_sample_size) '.csv']);

                    else
                        writematrix(updated_data,['data_post_iterative_sampling_Bhatia_Exp_data_pd_mod_' num2str(pop_model_indx) '_b_t_ratio_' num2str(birth_tran_ratio) '_inverse_r1_' num2str(1/r1) '_negative_comp.csv']);
                        writematrix(current_profile_data,['profile_lik_para_post_itr_samp_pd_mod_' num2str(pop_model_indx) '_b_t_ratio_' num2str(birth_tran_ratio) '_inverse_r1_' num2str(1/r1)  '_para_steps_' num2str(profile_sample_size) '_negative_comp.csv']);

                    end
                else
                    writematrix(updated_data,['data_post_iterative_sampling_Yamamoto_Exp_data_pd_mod_' num2str(pop_model_indx) '_b_t_ratio_' num2str(birth_tran_ratio) '_inverse_r1_' num2str(1/r1)  '.csv']);
                    %                 writematrix(current_profile_data,['profile_lik_para_post_itr_samp_pd_mod_' num2str(pop_model_indx) '_b_t_ratio_' num2str(birth_tran_ratio) '_inverse_r1_' num2str(1/r1)  '_dynamic_para_steps_default_value_' num2str(profile_sample_size) '.csv']);
                    writematrix(current_profile_data,['profile_lik_para_post_itr_samp_pd_mod_' num2str(pop_model_indx) '_b_t_ratio_' num2str(birth_tran_ratio) '_inverse_r1_' num2str(1/r1)  '_para_steps_' num2str(profile_sample_size) '.csv']);
                end
            else
                if(pop_model_indx  < 15)
                    if(pop_model_indx ~= 4)
                        writematrix(updated_data,['data_post_iterative_sampling_Bhatia_Exp_data_pd_mod_' num2str(pop_model_indx) '_b_t_ratio_' num2str(birth_tran_ratio) '_inverse_r1_' num2str(1/r1) '_inverse_r2_' num2str(inverse_r2) '.csv']);
                        %                     writematrix(current_profile_data,['profile_lik_para_post_itr_samp_pd_mod_' num2str(pop_model_indx) '_b_t_ratio_' num2str(birth_tran_ratio) '_inverse_r1_' num2str(1/r1)  '_dynamic_para_steps_default_value_' num2str(profile_sample_size) '.csv']);
                        writematrix(current_profile_data,['profile_lik_para_post_itr_samp_pd_mod_' num2str(pop_model_indx) '_b_t_ratio_' num2str(birth_tran_ratio) '_inverse_r1_' num2str(1/r1)  '_para_steps_' num2str(profile_sample_size) '_inverse_r2_' num2str(inverse_r2) '.csv']);

                    else
                        writematrix(updated_data,['data_post_iterative_sampling_Bhatia_Exp_data_pd_mod_' num2str(pop_model_indx) '_b_t_ratio_' num2str(birth_tran_ratio) '_inverse_r1_' num2str(1/r1) '_inverse_r2_' num2str(inverse_r2) '_negative_comp.csv']);
                        writematrix(current_profile_data,['profile_lik_para_post_itr_samp_pd_mod_' num2str(pop_model_indx) '_b_t_ratio_' num2str(birth_tran_ratio) '_inverse_r1_' num2str(1/r1) '_inverse_r2_' num2str(inverse_r2) '_para_steps_' num2str(profile_sample_size) '_negative_comp.csv']);

                    end
                else
                    writematrix(updated_data,['data_post_iterative_sampling_Yamamoto_Exp_data_pd_mod_' num2str(pop_model_indx) '_b_t_ratio_' num2str(birth_tran_ratio) '_inverse_r1_' num2str(1/r1) '_inverse_r2_' num2str(inverse_r2) '_output_2_total_cell_population.csv']);
                    %                 writematrix(current_profile_data,['profile_lik_para_post_itr_samp_pd_mod_' num2str(pop_model_indx) '_b_t_ratio_' num2str(birth_tran_ratio) '_inverse_r1_' num2str(1/r1)  '_dynamic_para_steps_default_value_' num2str(profile_sample_size) '.csv']);
                    writematrix(current_profile_data,['profile_lik_para_post_itr_samp_pd_mod_' num2str(pop_model_indx) '_b_t_ratio_' num2str(birth_tran_ratio) '_inverse_r1_' num2str(1/r1) '_inverse_r2_' num2str(inverse_r2) '_para_steps_' num2str(profile_sample_size) '_output_2_total_cell_population.csv']);
                end
            end

        else
            disp('Sample iteration is already present for the given parametric condition')
        end


    end
end


