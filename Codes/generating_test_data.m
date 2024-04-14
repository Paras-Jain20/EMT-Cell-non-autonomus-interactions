% data fetched below is not being used in generating test data, it is just
% to fill in the variable while call a function

total_test_data = 1000;

cd('C:\Users\user\OneDrive - Indian Institute of Science\Projects\State transition\Experimental data');
data = 1; % dummy variable of no use other than initializing para_search function; can take any value.

cd('C:\Users\user\OneDrive - Indian Institute of Science\GitHub\EMT-State-Transition-2.0\Code files')
% charactersitics of test data
data_type = 1; % test_data = 1; exp_data = 2
interpolation = 0;

r1 = 1;
inverse_r2 = 35;
r2 = []; %1/(r1*inverse_r2);

pop_model = [27 25 32 33];% [20 26 27 25 21 23 32 33];
% pop_model = [3 11 12 14 10 13 4];


search_method = 2; % parameter search method
parfor pop_model_indx = 1:length(pop_model)%[3 11 12 14 10 13]% % population growth model
    % if(find(pop_model_indx == pop_model,1))
        if(pop_model(pop_model_indx) <= 15) % for models capturing dynamics of Bhatia et al. data

            cd('C:\Users\user\OneDrive - Indian Institute of Science\GitHub\EMT-State-Transition-2.0\Parameters sets from experimental data')

            initial_dist = {'high', 'low'; 'low' , 'high'}; % for Bhatia et al data; Mes:Epi

            data = [];

            cd('C:\Users\user\OneDrive - Indian Institute of Science\Projects\State transition\Experimental data');

            % importing Bhatia et al data
            for ini_dist_indx = 1:2

                temp_data = readmatrix(['Bhatia et al 2019\M_' initial_dist{ini_dist_indx,1} '_E_' initial_dist{ini_dist_indx,2} '_dynamics.xlsx']);

                temp_data = [temp_data(:,1) temp_data(:,2:end)/100]; % changing population percentage into fractions
                data = [data;temp_data];
            end
            % data charactersitics
            num_rep = 3;
            num_init_cond = 2; % number of different initial condition in the data
            init_cond = [1 2];
            timept = 1:8; % time is unitless as it is scaled by m cells division rate (r1)
            interpolated_timept = 1; % interpolating time to get more data points to plot ODE solution trajectories
            num_readout = 2; % number of states captured in the experimental data

        elseif(pop_model(pop_model_indx) >=19)
            % Importing Yamamoto et al. 2017 data; models having
            % pop_model(pop_model_indx) >= 15
            cd('C:\Users\user\OneDrive - Indian Institute of Science\Projects\State transition\Experimental data\Yamamoto et al. 2017')
            data = table2array(readtable('Yamamoto et al 2017 Data 1 cell fraction.xlsx'));
            num_rep = 3;
            init_cond = [1 2 3 4 5];
            num_init_cond = 5; % number of different initial condition in the data
            timept = (0:3:12); % time is unitless as it is scaled by m cells division rate (r1)
            interpolated_timept = 1;  % interpolating time to get more data points to plot ODE solution trajectories
            num_readout = 2; % number of states captured in the experimental dat

        end

        num_time_pts = length(timept); % number of time pts per replicate
        
%         combined_test_data = zeros(total_test_data*num_time_pts*num_init_cond*num_rep,num_readout+2+1)

        for noise_factor = [5 25 50]
            total_itr = 1; % total number of parameters set to be checked/sampled

            birth_tran_ratio = 10; % max ratio between birth and transition
            test_data_matrix = zeros(total_test_data*(num_init_cond*length(timept)*num_rep),num_readout*2 + 2); % the first column is for time and last column denotes the nth test data for the given model 
            for test_data_indx  = 1:total_test_data
                disp(['Generating data for pop model ' num2str(pop_model(pop_model_indx)) ' test data number ' num2str(test_data_indx) ' noise factor ' num2str(noise_factor)])
                cd('C:\Users\user\OneDrive - Indian Institute of Science\GitHub\EMT-State-Transition-2.0\Parameter sets and generated pseudo data')
%                 if(~isfile(['test_data_' num2str(test_data_indx) '_pop_mod_' num2str(pop_model(pop_model_indx)) '_b_t_ratio_' num2str(birth_tran_ratio) '_noise_factor_' num2str(noise_factor) '.csv']))
                if(true)

                    cd('C:\Users\user\OneDrive - Indian Institute of Science\GitHub\EMT-State-Transition-2.0\Code files')
                    [select_para,goodness] = para_search(search_method,pop_model(pop_model_indx),total_itr,data_type, data, num_rep,init_cond, timept, interpolated_timept, interpolation, num_time_pts,birth_tran_ratio, r2); % select_para are filtered for re > rm

                    ode_sim_data = ODE_simulation(pop_model(pop_model_indx), num_rep,timept,interpolated_timept,interpolation, data,init_cond, select_para);

                    % changing cell number data into fractions
                    if(pop_model(pop_model_indx)>= 19)
                        ode_sim_data = [ode_sim_data(:,3)./(sum(ode_sim_data(:,3:4),2)) ode_sim_data(:,4)./(sum(ode_sim_data(:,3:4),2))];
                    end

                    % 2 belongs to number of distinct initial conditions
                    test_data = zeros(length(timept)*num_rep*num_init_cond,4);
                    if(pop_model(pop_model_indx) == 1)
                        std_level = (sqrt(ode_sim_data))/noise_factor; % for models with cell numbers as output
                    else
                        std_level = (-(ode_sim_data - 1/2).^2 +1/4)/noise_factor; % for models with cell fraction as output or models which have 4 distinct subpopulations
                    end

                    % adding noise to the data to get multiple replicates for a para set
                    while(true)
                        test_data(:,1:num_readout+1) = [repmat(timept',num_rep*num_init_cond ,1) ode_sim_data + std_level.*randn(size(ode_sim_data))];
                        if(pop_model(pop_model_indx) == 1)
                            break;
                        else
                            if(isempty(find(test_data(:,2:num_readout+1)<0,1)))
                                if(isempty(find(test_data(:,2:num_readout+1)>1,1)))
                                    break;
                                end
                            end
                        end
                    end

                    for readout_indx = 1:num_readout
                        for i = 1:length(timept)
                            for j = 1:num_init_cond
                                indx = (i+(j-1)*num_rep*length(timept):length(timept):(j-1)*num_rep*length(timept)+(num_rep-1)*length(timept)+i);
                                test_data(indx,1+num_readout+readout_indx) = std(test_data(indx,readout_indx+1),0,'all');
                            end
                        end
                    end

                    test_data_matrix(1+(test_data_indx-1)*size(test_data,1):test_data_indx*size(test_data,1),:) = [test_data ones(size(test_data,1),1)*test_data_indx];

%                     test_data = splitvars(table(test_data));
%                     parameters = splitvars(table(select_para));

                end
            end

            test_data_matrix = splitvars(table(test_data_matrix));
            test_data_matrix.Properties.VariableNames = ["Time", "Out_1", "Out_2", "Std_out_1", "Std_out_2", "Test_data"];
            directory = pwd;
            cd('C:\Users\user\OneDrive - Indian Institute of Science\GitHub\EMT-State-Transition-2.0\Parameter sets and generated pseudo data')
            writetable(test_data_matrix,['test_data_pop_mod_' num2str(pop_model(pop_model_indx)) '_b_t_ratio_' num2str(birth_tran_ratio) '_noise_factor_' num2str(noise_factor) '.csv']);
            %                     writetable(parameters,['common_test_para_' num2str(test_data_indx) '_pop_mod_' num2str(pop_model(pop_model_indx)) '_b_t_ratio_' num2str(birth_tran_ratio) '_noise_factor_' num2str(noise_factor) '.csv']);
            cd(directory)

        end
    % else
    %     continue;
    % end
end