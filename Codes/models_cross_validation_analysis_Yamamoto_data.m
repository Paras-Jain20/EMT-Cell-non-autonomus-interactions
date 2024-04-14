%%%%%%%%% main_function %%%%%%%%%
close all

% whether we are using the test or experimental data: test_data = 1; exp_data = 2
data_type = 2;

% set of population models
pop_model =  [20 26 27 25 21 23]; %[3 4 11 12 14 10 13]; %
birth_tran_ratio = 250;

table_variable_names =     {["Init Cond 1"    "Init Cond 2"    "Init Cond 3"    "Init Cond 4"    "Init Cond 5"    "Pop mod"    "Itr indx"];
    ["Init Cond 1"    "Init Cond 2"    "Pop mod"    "Itr indx"]};

total_itr =10; % total number of parameters set to be checked/sampled

% experimental data characteristics
if(~isempty(find(pop_model > 15,1)))
    num_rep = 3;
    num_init_cond = 5; % number of different initial condition in the data
    init_cond = [1 2 3 4 5];
    % Importing experimental data
    cd('C:\Users\user\OneDrive - Indian Institute of Science\Projects\State transition\Experimental data\Yamamoto et al. 2017')
    data = table2array(readtable('Yamamoto et al 2017 Data 1 cell fraction.xlsx'));
elseif (~isempty(find(pop_model < 15,1)))
    num_rep = 3;
    num_init_cond = 2; % number of different initial condition in the data
    init_cond = [1 2];
    % Importing experimental data
    cd('C:\Users\user\OneDrive - Indian Institute of Science\Projects\State transition\Experimental data\')
    initial_dist = {'high', 'low'; 'low' , 'high'}; % for Bhatia et al data; Mes:Epi
    data = [];
    for ini_dist_indx = 1:2

        temp_data = xlsread(['Bhatia et al 2019\M_' initial_dist{ini_dist_indx,1} '_E_' initial_dist{ini_dist_indx,2} '_dynamics.xlsx']);

        temp_data = [temp_data(:,1) temp_data(:,2:end)/100]; % changing population percentage into fractions
        data = [data;temp_data];
    end
end


test_data_indx = 1; % necessary dummy variable to just to fill in variable list in profile likelihood function
noise_factor = 1; % necessary dummy variable to just to fill in variable list in profile likelihood function
r1 = 1/54; % m cells division rates; scaling factor for non dimensional population model
inverse_r2 = 35;
r2 = 1/(inverse_r2*r1);

%%%%%%%%%%%%%%%%%%%%%%% cross 1 validation analysis %%%%%%%%%%%%%%%%%%%%%
cd('C:\Users\user\OneDrive - Indian Institute of Science\GitHub\EMT-State-Transition-2.0\Parameter sets and generated pseudo data')

if(isempty(r2))
if(~isempty(find(pop_model > 15,1)))
    if(isfile(['cross_1_validation_chi_sqaure_values_Yamamoto_Exp_data_b_t_ratio_' num2str(birth_tran_ratio) '_inverse_r1_' num2str(1/r1)  '.csv']))
        file_present = true;
        cost = table2array(readtable(['cross_1_validation_chi_sqaure_values_Yamamoto_Exp_data_b_t_ratio_' num2str(birth_tran_ratio) '_inverse_r1_' num2str(1/r1)  '.csv']));
        pop_mod_itr_indx = size(cost,1)+1;
    else
        cost = zeros(length(pop_model)*total_itr,length(init_cond)+2);
        file_present = false;
        pop_mod_itr_indx = 1;
    end
elseif (~isempty(find(pop_model < 15,1)))
    if(isfile(['cross_1_validation_chi_sqaure_values_Bhatia_Exp_data_b_t_ratio_' num2str(birth_tran_ratio) '_inverse_r1_' num2str(1/r1)  '.csv']))
        file_present = true;
        cost = table2array(readtable(['cross_1_validation_chi_sqaure_values_Bhatia_Exp_data_b_t_ratio_' num2str(birth_tran_ratio) '_inverse_r1_' num2str(1/r1)  '.csv']));
        pop_mod_itr_indx = size(cost,1)+1;
    else
        cost = zeros(length(pop_model)*total_itr,length(init_cond)+2);
        file_present = false;
        pop_mod_itr_indx = 1;
    end
end
else
if(~isempty(find(pop_model > 15,1)))
    if(isfile(['cross_1_validation_chi_sqaure_values_Yamamoto_Exp_data_b_t_ratio_' num2str(birth_tran_ratio) '_inverse_r1_' num2str(1/r1)  '.csv']))
        file_present = true;
        cost = table2array(readtable(['cross_1_validation_chi_sqaure_values_Yamamoto_Exp_data_b_t_ratio_' num2str(birth_tran_ratio) '_inverse_r1_' num2str(1/r1)  '_inverse_r2_' num2str(inverse_r2)  '.csv']));
        pop_mod_itr_indx = size(cost,1)+1;
    else
        cost = zeros(length(pop_model)*total_itr,length(init_cond)+2);
        file_present = false;
        pop_mod_itr_indx = 1;
    end
elseif (~isempty(find(pop_model < 15,1)))
    if(isfile(['cross_1_validation_chi_sqaure_values_Bhatia_Exp_data_b_t_ratio_' num2str(birth_tran_ratio) '_inverse_r1_' num2str(1/r1)  '.csv']))
        file_present = true;
        cost = table2array(readtable(['cross_1_validation_chi_sqaure_values_Bhatia_Exp_data_b_t_ratio_' num2str(birth_tran_ratio) '_inverse_r1_' num2str(1/r1) '_inverse_r2_' num2str(inverse_r2) '.csv']));
        pop_mod_itr_indx = size(cost,1)+1;
    else
        cost = zeros(length(pop_model)*total_itr,length(init_cond)+2);
        file_present = false;
        pop_mod_itr_indx = 1;
    end
end
end


% select a population growth model
pop_model_seq = 1;

for pop_model_indx = pop_model%[3 11 12 14 10 13]%[20 26 27 25 21 23]%[1 2 3 10 11 12 13]


        close all

        % type of optimum parameters search to be performed:
        % use existing parameters = choose 0
        % optimizer based : choose 1
        % uniform sampling based : choose 2; this method is not updated so might not work
        search_method = 1; % parameter search method

        % Importing Yamamoto et al. 2017 data; models having
        % data's temporal characteristics
        if(pop_model_indx >15)
            timept = (0:3:12)*24*r1; % time is unitless as it is scaled by m cells division rate (r1)
            interpolated_timept = 0:r1:timept(end);  % interpolating time to get more data points to plot ODE solution trajectories

        elseif(pop_model_indx)
            timept = (2:2:8)*7*24*r1; % time is unitless as it is scaled by m cells division rate (r1)
            interpolated_timept = r1:r1:timept(end);  % interpolating time to get more data points to plot ODE solution trajectories
        end


        num_time_pts = length(timept);

        for itr_indx = 1:total_itr
            if(isempty(find(sum(cost(:,end-1:end) == [pop_model(pop_model_seq) itr_indx],2) == 2,1)))

                for init_cond_indx = init_cond


                    test_data_indx = 1+(init_cond_indx-1)*num_rep*num_time_pts:(init_cond_indx)*num_rep*num_time_pts;
                    test_data = data(1+(init_cond_indx-1)*num_rep*num_time_pts:(init_cond_indx)*num_rep*num_time_pts,:);

                    train_data_indx = ones(size(data,1),1);
                    train_data_indx(test_data_indx) = 0;
                    train_data = data(logical(train_data_indx),:);

                    cd('C:\Users\user\OneDrive - Indian Institute of Science\GitHub\EMT-State-Transition-2.0\Code files')

                    %                 if(true)
                    interpolation = 0; % to add more time points for the ODE simulated data; must be zero if we using objective function to minimize error with the data
                    % parameter search function
                    init_cond_dyn = init_cond(init_cond ~= init_cond(init_cond_indx)); % change set of initial condition for cross validation analysis
                    [select_para,goodness] = para_search(search_method,pop_model_indx,1, data_type, train_data, num_rep,init_cond_dyn, timept, interpolated_timept, interpolation, num_time_pts,birth_tran_ratio,r2); % select_para are filtered for re > rm

                    temp_cost = obj_fun(select_para,data_type, test_data,pop_model_indx, num_rep,timept,interpolated_timept,init_cond_indx,interpolation);
                    cost(pop_mod_itr_indx,init_cond_indx) = sum(temp_cost.^2);

                end
                cost(pop_mod_itr_indx,end-1) = pop_model(pop_model_seq);
                cost(pop_mod_itr_indx,end) = itr_indx;
                pop_mod_itr_indx = pop_mod_itr_indx + 1;

            end
        end

        pop_model_seq = pop_model_seq + 1;
end
cost = splitvars(table(cost));

cd('C:\Users\user\OneDrive - Indian Institute of Science\GitHub\EMT-State-Transition-2.0\Parameter sets and generated pseudo data')

if(isempty(r2))
if(~isempty(find(pop_model > 15,1)))
    cost.Properties.VariableNames = table_variable_names{1};
    writetable(cost, ['cross_1_validation_chi_sqaure_values_Yamamoto_Exp_data_b_t_ratio_' num2str(birth_tran_ratio) '_inverse_r1_' num2str(1/r1)  '.csv']);
elseif (~isempty(find(pop_model < 15,1)))
    cost.Properties.VariableNames = table_variable_names{2};
    writetable(cost, ['cross_1_validation_chi_sqaure_values_Bhatia_Exp_data_b_t_ratio_' num2str(birth_tran_ratio) '_inverse_r1_' num2str(1/r1)  '.csv']);
end
else
    if(~isempty(find(pop_model > 15,1)))
        cost.Properties.VariableNames = table_variable_names{1};
        writetable(cost, ['cross_1_validation_chi_sqaure_values_Yamamoto_Exp_data_b_t_ratio_' num2str(birth_tran_ratio) '_inverse_r1_' num2str(1/r1) '_inverse_r2_' num2str(inverse_r2) '.csv']);
    elseif (~isempty(find(pop_model < 15,1)))
        cost.Properties.VariableNames = table_variable_names{2};
        writetable(cost, ['cross_1_validation_chi_sqaure_values_Bhatia_Exp_data_b_t_ratio_' num2str(birth_tran_ratio) '_inverse_r1_' num2str(1/r1) '_inverse_r2_' num2str(inverse_r2) '.csv']);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%% control for cross 1 validation analysis %%%%%%%%%%%%%%%%%%%%%

cd('C:\Users\user\OneDrive - Indian Institute of Science\GitHub\EMT-State-Transition-2.0\Parameter sets and generated pseudo data')

if(isempty(r2))
    if(~isempty(find(pop_model > 15,1)))
        if(isfile(['Control_cross_1_validation_chi_sqaure_values_Yamamoto_Exp_data_b_t_ratio_' num2str(birth_tran_ratio) '_inverse_r1_' num2str(1/r1)  '.csv']))
            file_present = true;
            control_cost = table2array(readtable(['Control_cross_1_validation_chi_sqaure_values_Yamamoto_Exp_data_b_t_ratio_' num2str(birth_tran_ratio) '_inverse_r1_' num2str(1/r1) '_inverse_r2_' num2str(inverse_r2)   '.csv']));
            pop_mod_itr_indx = size(control_cost,1)+1;
        else
            control_cost = zeros(length(pop_model)*total_itr,length(init_cond)+2);
            file_present = false;
            pop_mod_itr_indx = 1;
        end
    elseif (~isempty(find(pop_model < 15,1)))
        if(isfile(['Control_cross_1_validation_chi_sqaure_values_Bhatia_Exp_data_b_t_ratio_' num2str(birth_tran_ratio) '_inverse_r1_' num2str(1/r1)  '.csv']))
            file_present = true;
            control_cost = table2array(readtable(['Control_cross_1_validation_chi_sqaure_values_Bhatia_Exp_data_b_t_ratio_' num2str(birth_tran_ratio) '_inverse_r1_' num2str(1/r1) '_inverse_r2_' num2str(inverse_r2)  '.csv']));
            pop_mod_itr_indx = size(control_cost,1)+1;
        else
            control_cost = zeros(length(pop_model)*total_itr,length(init_cond)+2);
            file_present = false;
            pop_mod_itr_indx = 1;
        end
    end
else
    if(~isempty(find(pop_model > 15,1)))
        if(isfile(['Control_cross_1_validation_chi_sqaure_values_Yamamoto_Exp_data_b_t_ratio_' num2str(birth_tran_ratio) '_inverse_r1_' num2str(1/r1)  '.csv']))
            file_present = true;
            control_cost = table2array(readtable(['Control_cross_1_validation_chi_sqaure_values_Yamamoto_Exp_data_b_t_ratio_' num2str(birth_tran_ratio) '_inverse_r1_' num2str(1/r1)  '.csv']));
            pop_mod_itr_indx = size(control_cost,1)+1;
        else
            control_cost = zeros(length(pop_model)*total_itr,length(init_cond)+2);
            file_present = false;
            pop_mod_itr_indx = 1;
        end
    elseif (~isempty(find(pop_model < 15,1)))
        if(isfile(['Control_cross_1_validation_chi_sqaure_values_Bhatia_Exp_data_b_t_ratio_' num2str(birth_tran_ratio) '_inverse_r1_' num2str(1/r1)  '.csv']))
            file_present = true;
            control_cost = table2array(readtable(['Control_cross_1_validation_chi_sqaure_values_Bhatia_Exp_data_b_t_ratio_' num2str(birth_tran_ratio) '_inverse_r1_' num2str(1/r1)  '.csv']));
            pop_mod_itr_indx = size(control_cost,1)+1;
        else
            control_cost = zeros(length(pop_model)*total_itr,length(init_cond)+2);
            file_present = false;
            pop_mod_itr_indx = 1;
        end
    end
end

% select a population growth model
pop_model_seq = 1;

for pop_model_indx = pop_model%[3 11 12 14 10 13]%[20 26 27 25 21 23]%[1 2 3 10 11 12 13]
     

        close all

        % type of optimum parameters search to be performed:
        % use existing parameters = choose 0
        % optimizer based : choose 1
        % uniform sampling based : choose 2; this method is not updated so might not work
        search_method = 1; % parameter search method

        % Importing Yamamoto et al. 2017 data; models having
        % data's temporal characteristics
        if(pop_model_indx >15)
            timept = (0:3:12)*24*r1; % time is unitless as it is scaled by m cells division rate (r1)
            interpolated_timept = 0:r1:timept(end);  % interpolating time to get more data points to plot ODE solution trajectories

        elseif(pop_model_indx)
            timept = (2:2:8)*7*24*r1; % time is unitless as it is scaled by m cells division rate (r1)
            interpolated_timept = r1:r1:timept(end);  % interpolating time to get more data points to plot ODE solution trajectories
        end


        num_time_pts = length(timept);

        for itr_indx = 1:total_itr
            if(isempty(find(sum(control_cost(:,end-1:end) == [pop_model(pop_model_seq) itr_indx],2) == 2,1)))

                for init_cond_indx = init_cond


                    test_data_indx = 1+(init_cond_indx-1)*num_rep*num_time_pts:(init_cond_indx)*num_rep*num_time_pts;
                    test_data = data(1+(init_cond_indx-1)*num_rep*num_time_pts:(init_cond_indx)*num_rep*num_time_pts,:);

                    train_data = data;

                    cd('C:\Users\user\OneDrive - Indian Institute of Science\GitHub\EMT-State-Transition-2.0\Code files')

                    %                 if(true)
                    interpolation = 0; % to add more time points for the ODE simulated data; must be zero if we using objective function to minimize error with the data
                    % parameter search function
%                     init_cond_dyn = init_cond(init_cond ~= init_cond(init_cond_indx)); % change set of initial condition for cross validation analysis
                    [select_para,goodness] = para_search(search_method,pop_model_indx,1, data_type, train_data, num_rep,init_cond, timept, interpolated_timept, interpolation, num_time_pts,birth_tran_ratio,r2); % select_para are filtered for re > rm

                    temp_control_cost = obj_fun(select_para,data_type, test_data,pop_model_indx, num_rep,timept,interpolated_timept,init_cond_indx,interpolation);
                    control_cost(pop_mod_itr_indx,init_cond_indx) = sum(temp_control_cost.^2);

                end
                control_cost(pop_mod_itr_indx,end-1) = pop_model(pop_model_seq);
                control_cost(pop_mod_itr_indx,end) = itr_indx;
                pop_mod_itr_indx = pop_mod_itr_indx + 1;

            end
        end

        pop_model_seq = pop_model_seq + 1;
end
control_cost = splitvars(table(control_cost));

cd('C:\Users\user\OneDrive - Indian Institute of Science\GitHub\EMT-State-Transition-2.0\Parameter sets and generated pseudo data')
if(isempty(r2))
    if(~isempty(find(pop_model > 15,1)))
        cost.Properties.VariableNames = table_variable_names{1};
        writetable(cost, ['Control_cross_1_validation_chi_sqaure_values_Yamamoto_Exp_data_b_t_ratio_' num2str(birth_tran_ratio) '_inverse_r1_' num2str(1/r1)  '.csv']);
    elseif (~isempty(find(pop_model < 15,1)))
        cost.Properties.VariableNames = table_variable_names{2};
        writetable(cost, ['Control_cross_1_validation_chi_sqaure_values_Bhatia_Exp_data_b_t_ratio_' num2str(birth_tran_ratio) '_inverse_r1_' num2str(1/r1)  '.csv']);
    end
else
    if(~isempty(find(pop_model > 15,1)))
        cost.Properties.VariableNames = table_variable_names{1};
        writetable(cost, ['Control_cross_1_validation_chi_sqaure_values_Yamamoto_Exp_data_b_t_ratio_' num2str(birth_tran_ratio) '_inverse_r1_' num2str(1/r1) '_inverse_r2_' num2str(inverse_r2) '.csv']);
    elseif (~isempty(find(pop_model < 15,1)))
        cost.Properties.VariableNames = table_variable_names{2};
        writetable(cost, ['Control_cross_1_validation_chi_sqaure_values_Bhatia_Exp_data_b_t_ratio_' num2str(birth_tran_ratio) '_inverse_r1_' num2str(1/r1) '_inverse_r2_' num2str(inverse_r2) '.csv']);
    end
end


if(length(init_cond)>2)
    %%%%%%%%%%%%%%%%%%%%%%%%% cross 2 validation analysis %%%%%%%%%%%%%%%%%%%%%
    init_cond_comb = [];
    for i = 1:length(init_cond)
        for j = i+1:length(init_cond)
            init_cond_comb  = [init_cond_comb; i j];
        end
    end

    cost = zeros(length(pop_model)*total_itr,length(init_cond_comb));

    for init_cond_indx = 1:size(init_cond_comb,1)
        % select a population growth model
        pop_model_seq = 1;

        for pop_model_indx = pop_model%[3 11 12 14 10 13]%[20 26 27 25 21 23]%[1 2 3 10 11 12 13]
            for birth_tran_ratio = 150%[5 10 20 50 100]%[50 100 150 200 250 300] %200 500] %[5 10 20 50 100] % max ratio between birth and transition
                for r1_indx = 1:length(r1_set)
                    close all

                    % type of optimum parameters search to be performed:
                    % use existing parameters = choose 0
                    % optimizer based : choose 1
                    % uniform sampling based : choose 2; this method is not updated so might not work
                    search_method = 1; % parameter search method

                    switch pop_model_indx

                        case {1,2,5} % models with time dimension
                            r1 = 1;
                        case {3,6,10,11,12,13,14,15,16,17,18,19,20,21,22, 23, 24, 25, 26, 27} % models without time dimension
                            r1 = r1_set(r1_indx);
                    end

                    % Importing Yamamoto et al. 2017 data; models having
                    % data's temporal characteristics
                    timept = (0:3:12)*24*r1; % time is unitless as it is scaled by m cells division rate (r1)
                    interpolated_timept = 0:r1:timept(end);  % interpolating time to get more data points to plot ODE solution trajectories

                    num_time_pts = length(timept);

                    test_data_indx = [1+(init_cond_comb(init_cond_indx,1)-1)*num_rep*num_time_pts:init_cond_comb(init_cond_indx,1)*num_rep*num_time_pts 1+(init_cond_comb(init_cond_indx,2)-1)*num_rep*num_time_pts:init_cond_comb(init_cond_indx,2)*num_rep*num_time_pts];
                    test_data = data(test_data_indx,:);

                    train_data_indx = ones(size(data,1),1);
                    train_data_indx(test_data_indx) = 0;
                    train_data = data(logical(train_data_indx),:);

                    cd('C:\Users\user\OneDrive - Indian Institute of Science\GitHub\EMT-State-Transition-2.0\Code files')

                    %                 if(true)
                    interpolation = 0; % to add more time points for the ODE simulated data; must be zero if we using objective function to minimize error with the data
                    % parameter search function
                    init_cond_dyn = init_cond(logical((init_cond ~= init_cond_comb(init_cond_indx,1)) .* (init_cond ~= init_cond_comb(init_cond_indx,2)))); % change set of initial condition for cross validation analysis
                    [select_para,goodness] = para_search(search_method,pop_model_indx,total_itr, data_type, train_data, num_rep,init_cond_dyn, timept, interpolated_timept, interpolation, num_time_pts,birth_tran_ratio); % select_para are filtered for re > rm

                    for itr_indx = 1:total_itr
                        temp_cost = obj_fun(select_para(itr_indx,:),data_type, test_data,pop_model_indx, num_rep,timept,interpolated_timept,init_cond_comb(init_cond_indx,:),interpolation);

                        cost(pop_model_seq+(length(pop_model)*(itr_indx-1)),init_cond_indx) = sum(temp_cost.^2);
                    end
                    %                 end
                end
            end
            pop_model_seq = pop_model_seq + 1;
        end
    end

    cd('C:\Users\user\OneDrive - Indian Institute of Science\GitHub\EMT-State-Transition-2.0\Parameter sets and generated pseudo data')

    dlmwrite(['cross_2_validation_chi_sqaure_values_Yamamoto_Exp_data_b_t_ratio_' num2str(birth_tran_ratio) '_inverse_r1_' num2str(1/r1)  '.csv'],cost, '-append');
end

if(length(init_cond)>2)
    %%%%%%%%%%%%%%%%%%%%%%%%% cross 2 validation analysis %%%%%%%%%%%%%%%%%%%%%
    init_cond_comb = [];
    for i = 1:length(init_cond)
        for j = i+1:length(init_cond)
            init_cond_comb  = [init_cond_comb; i j];
        end
    end

    cost = zeros(length(pop_model)*total_itr,length(init_cond_comb));

    for init_cond_indx = 1:size(init_cond_comb,1)
        % select a population growth model
        pop_model_seq = 1;

        for pop_model_indx = pop_model%[3 11 12 14 10 13]%[20 26 27 25 21 23]%[1 2 3 10 11 12 13]
            for birth_tran_ratio = 150%[5 10 20 50 100]%[50 100 150 200 250 300] %200 500] %[5 10 20 50 100] % max ratio between birth and transition
                for r1_indx = 1:length(r1_set)
                    close all

                    % type of optimum parameters search to be performed:
                    % use existing parameters = choose 0
                    % optimizer based : choose 1
                    % uniform sampling based : choose 2; this method is not updated so might not work
                    search_method = 1; % parameter search method

                    switch pop_model_indx

                        case {1,2,5} % models with time dimension
                            r1 = 1;
                        case {3,6,10,11,12,13,14,15,16,17,18,19,20,21,22, 23, 24, 25, 26, 27} % models without time dimension
                            r1 = r1_set(r1_indx);
                    end

                    % Importing Yamamoto et al. 2017 data; models having
                    % data's temporal characteristics
                    timept = (0:3:12)*24*r1; % time is unitless as it is scaled by m cells division rate (r1)
                    interpolated_timept = 0:r1:timept(end);  % interpolating time to get more data points to plot ODE solution trajectories

                    num_time_pts = length(timept);

                    test_data_indx = [1+(init_cond_comb(init_cond_indx,1)-1)*num_rep*num_time_pts:init_cond_comb(init_cond_indx,1)*num_rep*num_time_pts 1+(init_cond_comb(init_cond_indx,2)-1)*num_rep*num_time_pts:init_cond_comb(init_cond_indx,2)*num_rep*num_time_pts];
                    test_data = data(test_data_indx,:);

                    train_data = data;


                    %                 if(true)
                    interpolation = 0; % to add more time points for the ODE simulated data; must be zero if we using objective function to minimize error with the data
                    % parameter search function
                    cd('C:\Users\user\OneDrive - Indian Institute of Science\GitHub\EMT-State-Transition-2.0\Parameters sets from experimental data\')

                    if(isfile(['para_Yamamoto_data_pd_mod_' num2str(pop_model_indx) '_s_m_' num2str(1) '_b_t_ratio_' num2str(birth_tran_ratio) '_inverse_r1_' num2str(1/r1) '.csv']))


                        para_data = table2array(readtable(['para_Yamamoto_data_pd_mod_' num2str(pop_model_indx) '_s_m_' num2str(1) '_b_t_ratio_' num2str(birth_tran_ratio) '_inverse_r1_' num2str(1/r1) '.csv']));

                        select_para = para_data(:,1:end-1);

                    else
                        cd('C:\Users\user\OneDrive - Indian Institute of Science\GitHub\EMT-State-Transition-2.0\Code files')

                        [select_para,goodness] = para_search(search_method,pop_model_indx,total_itr, data_type, train_data, num_rep,init_cond, timept, interpolated_timept, interpolation, num_time_pts,birth_tran_ratio); % select_para are filtered for re > rm

                    end

                    cd('C:\Users\user\OneDrive - Indian Institute of Science\GitHub\EMT-State-Transition-2.0\Code files')


                    for itr_indx = 1:total_itr
                        temp_cost = obj_fun(select_para(itr_indx,:),data_type, test_data,pop_model_indx, num_rep,timept,interpolated_timept,init_cond_comb(init_cond_indx,:),interpolation);

                        cost(pop_model_seq+(length(pop_model)*(itr_indx-1)),init_cond_indx) = sum(temp_cost.^2);
                    end
                    %                 end
                end
            end
            pop_model_seq = pop_model_seq + 1;
        end
    end

    cd('C:\Users\user\OneDrive - Indian Institute of Science\GitHub\EMT-State-Transition-2.0\Parameter sets and generated pseudo data')

    dlmwrite(['Control_cross_2_validation_chi_sqaure_values_Yamamoto_Exp_data_b_t_ratio_' num2str(birth_tran_ratio) '_inverse_r1_' num2str(1/r1)  '.csv'],cost, '-append');
end