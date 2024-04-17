function [select_para,goodness] = para_search(search_method,pop_model_indx,total_itr, data_type, data, num_rep,timept,interpolated_timept,interpolation, num_time_pts, birth_tran_ratio)

rng('shuffle')
goodness = zeros(total_itr,1);

switch pop_model_indx
    case{6} % i.e. population model with 1 parameter
        num_model_para = 1;
    case {1,2,10,11,12,20} % i.e. population model with 4 parameters
        num_model_para = 4;
    case {3,7,18} % i.e. population model with 3 parameters
        num_model_para = 3;
    case {4,5,13,14,15,16,17,19, 26,27} % i.e. population model with 5 parameters
        num_model_para = 5;
    case {8,9,22,23,24,25} % i.e. population model with 6 parameters
        num_model_para = 6;
    case {21} % i.e. population model with 7 parameters
        num_model_para = 7;
end

select_para = zeros(total_itr,num_model_para);

parfor itr_num = 1:total_itr
    switch pop_model_indx

        case {1,2} % model: exponential grwoth with transition - cell numbers and cell fraction dynamics
            parameters = zeros(4,1);
            parameters(1,1) = 1/(20 + 40*rand); % grwoth rate of m cells
            parameters(2,1) = 1/(20 + 40*rand); % growth rate of e cells
            parameters(3,1) = parameters(1,1)/(1+birth_tran_ratio*rand); % M E transition rate
            parameters(4,1) = parameters(2,1)/(1+birth_tran_ratio*rand); % E M transition rate

            lower_bound = [1/60 1/60 1/(60*birth_tran_ratio) 1/(60*birth_tran_ratio)];
            upper_bound = [1/20 1/20 1/20 1/20];

        case {3, 18} % model: exponential grwoth with transition - cell fraction dynamics non dimensional
            parameters = zeros(3,1);
            parameters(1,1) = (1 + 2*rand); % growth rate of e cells wrt growth rate of m cells
            parameters(2,1) = 1/(1+birth_tran_ratio*rand); % M E transition rate wrt growth rate of m cells
            parameters(3,1) = 1/(1+birth_tran_ratio*rand); % E M transition rate wrt growth rate of m cells

            lower_bound = [1 1/(birth_tran_ratio) 1/(birth_tran_ratio)];
            upper_bound = [3 1 1];

        case {4}
            parameters = zeros(5,1);
            parameters(1,1) = (1 + 2*rand); % growth rate of e cells
            parameters(2,1) = 1/(1+birth_tran_ratio*rand); % M E basal transition rate
            parameters(3,1) = 1/(1+birth_tran_ratio*rand); % E M basal transition rate
            parameters(4,1) = 1 * rand; % Retention constant of E or M population
            parameters(5,1) = 1 * rand; % Influence constant of M population

            lower_bound = [1 1/(birth_tran_ratio)  1/(birth_tran_ratio)  -1 -1];
            upper_bound = [3 1 1 1 1];

        case {5}

            parameters = zeros(2,1);
            parameters(1,1) = 1/(20 + 40*rand); % grwoth rate of m cells
            parameters(2,1) = 1/(20 + 40*rand); % growth rate of e cells

            lower_bound = [1/60 1/60];
            upper_bound = [1/20 1/20];

        case {6}
            parameters = zeros(1,1);
            parameters(1,1) = (1 + 2*rand); % grwoth rate of m cells

            lower_bound = [1];
            upper_bound = [3];

        case {7}
            parameters = zeros(3,1);
            parameters(1,1) = (1 + 2*rand); % growth rate of e cells
            parameters(2,1) = 1 * rand; % Retention constant of E or M population
            parameters(3,1) = 1 * rand; % Influence constant of M population

            lower_bound = [1 0 0];
            upper_bound = [3 1 1];



        case {8,9}
            parameters = zeros(6,1);
            parameters(1,1) = 1/(20 + 40*rand); % grwoth rate of m cells
            parameters(2,1) = 1/(20 + 40*rand); % growth rate of e cells
            parameters(3,1) = parameters(1,1)/(1+birth_tran_ratio*rand); % M E basal transition rate
            parameters(4,1) = parameters(2,1)/(1+birth_tran_ratio*rand); % E M basal transition rate
            parameters(5,1) = 1 * rand; % Retention constant of E or M population
            parameters(6,1) = 4 * rand;% influence constant of M population

            lower_bound = [1/60 1/60 1/(60*birth_tran_ratio)  1/(60*birth_tran_ratio)  0 0];
            upper_bound = [1/20 1/20 1/20 1/20 1 4];

        case {10}
            parameters = zeros(4,1);
            parameters(1,1) = (1 + 2*rand); % growth rate of e cells wrt growth rate of m cells
            parameters(2,1) = 1/(1+birth_tran_ratio*rand); % M E basal transition rate wrt growth rate of m cells
            parameters(3,1) = 1/(1+birth_tran_ratio*rand); % E M basal transition rate wrt growth rate of m cells
            parameters(4,1) = 4 * rand; % Influence constant of M pop

            lower_bound = [1 1/(birth_tran_ratio)  1/(birth_tran_ratio)  0];
            upper_bound = [3 1 1 4];


        case {11,12}
            parameters = zeros(4,1);
            parameters(1,1) = (1 + 2*rand); % growth rate of e cells
            parameters(2,1) = 1/(1+birth_tran_ratio*rand); % M E basal transition rate
            parameters(3,1) = 1/(1+birth_tran_ratio*rand); % E M basal transition rate
            parameters(4,1) = 1 * rand; % Retention constant of E or M pop

            lower_bound = [1 1/(birth_tran_ratio)  1/(birth_tran_ratio)  0];
            upper_bound = [3 1 1 1];


        case {13}
            parameters = zeros(5,1);
            parameters(1,1) = (1 + 2*rand); % growth rate of e cells
            parameters(2,1) = 1/(1+birth_tran_ratio*rand); % M E basal transition rate
            parameters(3,1) = 1/(1+birth_tran_ratio*rand); % E M basal transition rate
            parameters(4,1) = 1 * rand; % Retention constant of E or M population
            parameters(5,1) = 4 * rand; % Influence constant of M population

            lower_bound = [1 1/(birth_tran_ratio)  1/(birth_tran_ratio)  0 0];
            upper_bound = [3 1 1 1 4];

        case {14}
            parameters = zeros(5,1);
            parameters(1,1) = (1 + 2*rand); % growth rate of e cells
            parameters(2,1) = 1/(1+birth_tran_ratio*rand); % M E basal transition rate
            parameters(3,1) = 1/(1+birth_tran_ratio*rand); % E M basal transition rate
            parameters(4,1) = 1 * rand; % Retention constant of E or M population
            parameters(5,1) = 1 * rand; % Influence constant of M population

            lower_bound = [1 1/(birth_tran_ratio)  1/(birth_tran_ratio)  0 0];
            upper_bound = [3 1 1 1 1];


        case {15}
            parameters = zeros(5,1);
            parameters(1,1) = (1 + 2*rand); % growth rate of e cells
            parameters(2,1) = 1/(1+birth_tran_ratio*rand); % M E basal transition rate
            parameters(3,1) = 1/(1+birth_tran_ratio*rand); % E M basal transition rate
            parameters(4,1) = 4 * rand; % Influence constant of E population
            parameters(5,1) = 4 * rand; % Influence constant of M population

            lower_bound = [1 1/(birth_tran_ratio)  1/(birth_tran_ratio)  0 0];
            upper_bound = [3 1 1 10 10];


        case {16,17}
            parameters = zeros(5,1);
            parameters(1,1) = (1 + 2*rand); % growth rate of e cells
            parameters(2,1) = 1/(1+birth_tran_ratio*rand); % M E basal transition rate
            parameters(3,1) = 1/(1+birth_tran_ratio*rand); % E M basal transition rate
            parameters(4,1) = 1 * rand; % Influence constant of E population
            parameters(5,1) = 4 * rand; % Influence constant of M populatio


            %             lower_bound = [1 1/(birth_tran_ratio)  1/(birth_tran_ratio)  0 0];
            %             upper_bound = [3 1 1 1 10];
            lower_bound = [0  1/(birth_tran_ratio) 1/(birth_tran_ratio)  0 0];
            upper_bound = [3 1 1 1 10];

        case {19}
            parameters = zeros(2,1);
            parameters(1,1) = (1 + 2*rand); % growth rate of e cells
            parameters(2,1) = 10000000; % carrying capacity of M cells

            lower_bound = [1 1000000];
            upper_bound = [5 20000000];

        case {20}
            parameters = zeros(4,1);
            parameters(1,1) = (1 + 2*rand); % growth rate of e cells
            parameters(2,1) = 1/(1+birth_tran_ratio*rand); % non-dimensional rescaled M to E transiton rates
            parameters(3,1) = 1/(1+birth_tran_ratio*rand); % non-dimensional rescaled E to M transiton rates
            parameters(4,1) = 10000000; % carrying capacity of M cells

            lower_bound = [1 1/birth_tran_ratio 1/birth_tran_ratio 1000000];
            upper_bound = [5 1 1 20000000];

        case {21}
            parameters = zeros(7,1);
            parameters(1,1) = (1 + 2*rand); % growth rate of e cells
            parameters(2,1) = 1/(1+birth_tran_ratio*rand); % non-dimensional rescaled M to E transiton rates
            parameters(3,1) = 1/(1+birth_tran_ratio*rand); % non-dimensional rescaled E to M transiton rates
            parameters(4,1) = rand; % E competition to M cells
            parameters(5,1) = rand; % M competition to E cells
            parameters(6,1) = 10000000; % carrying capacity of M cells
            parameters(7,1) = 10000000; % carrying capacity of E cells

            lower_bound = [1 1/birth_tran_ratio 1/birth_tran_ratio 0 0  1000000 1000000 ];
            upper_bound = [5 1 1 10 10 20000000 20000000];

        case {22}
            parameters = zeros(6,1);
            parameters(1,1) = (1 + 2*rand); % growth rate of e cells
            parameters(2,1) = 1/(1+birth_tran_ratio*rand); % non-dimensional rescaled M to E transiton rates
            parameters(3,1) = 1/(1+birth_tran_ratio*rand); % non-dimensional rescaled E to M transiton rates
            parameters(4,1) = 10000000; % carrying capacity of M cells
            parameters(5,1) = 4 * rand; % Influence constant of E population
            parameters(6,1) = 4 * rand; % Influence constant of M population

            lower_bound = [1 1/birth_tran_ratio 1/birth_tran_ratio 1000000 0 0];
            upper_bound = [5 1 1 20000000 10 10];

        case {23}
            parameters = zeros(6,1);
            parameters(1,1) = (1 + 2*rand); % growth rate of e cells
            parameters(2,1) = 1/(1+birth_tran_ratio*rand); % non-dimensional rescaled M to E transiton rates
            parameters(3,1) = 1/(1+birth_tran_ratio*rand); % non-dimensional rescaled E to M transiton rates
            parameters(4,1) = 10000000; % carrying capacity of M cells
            parameters(5,1) = 1 * rand; % Retetion constant of E population
            parameters(6,1) = 4 * rand; % Influence constant of M population

            lower_bound = [1 1/birth_tran_ratio 1/birth_tran_ratio 1000000 0 0];
            upper_bound = [5 1 1 20000000 1 10];

        case {24}
            parameters = zeros(6,1);
            parameters(1,1) = (1 + 2*rand); % growth rate of e cells
            parameters(2,1) = 1/(1+birth_tran_ratio*rand); % non-dimensional rescaled M to E transiton rates
            parameters(3,1) = 1/(1+birth_tran_ratio*rand); % non-dimensional rescaled E to M transiton rates
            parameters(4,1) = 10000000; % carrying capacity of M cells

            parameters(5,1) = 1 * rand; % retention constant of M population
            parameters(6,1) = 4 * rand; % Influence constant of M population

            lower_bound = [1  1/birth_tran_ratio 1/birth_tran_ratio 1000000 0 0];
            upper_bound = [5 1 1 20000000 1 10];

        case {25}
            parameters = zeros(6,1);
            parameters(1,1) = (1 + 2*rand); % growth rate of e cells
            parameters(2,1) = 1/(1+birth_tran_ratio*rand); % non-dimensional rescaled M to E transiton rates
            parameters(3,1) = 1/(1+birth_tran_ratio*rand); % non-dimensional rescaled E to M transiton rates
            parameters(4,1) = 10000000; % carrying capacity of M cells

            parameters(5,1) = 1 * rand; % retention constant of M population
            parameters(6,1) = 1 * rand; % Influence constant of M population

            lower_bound = [1  1/birth_tran_ratio 1/birth_tran_ratio 1000000 0 0];
            upper_bound = [5 1 1 20000000 1 1];

        case {26,27}
            parameters = zeros(5,1);
            parameters(1,1) = (1 + 2*rand); % growth rate of e cells
            parameters(2,1) = 1/(1+birth_tran_ratio*rand); % non-dimensional rescaled M to E transiton rates
            parameters(3,1) = 1/(1+birth_tran_ratio*rand); % non-dimensional rescaled E to M transiton rates
            parameters(4,1) = 10000000; % carrying capacity of M cells

            parameters(5,1) = 1 * rand; % retention constant of M population

            lower_bound = [1  1/birth_tran_ratio 1/birth_tran_ratio 1000000 0];
            upper_bound = [5 1 1 20000000 1];

        case {28}
            parameters = zeros(13,1);
            parameters(1,1) = rand; % growth rate of e cells
            parameters(2,1) = rand; % non-dimensional rescaled M to E transiton rates
            parameters(3,1) = rand; % non-dimensional rescaled E to M transiton rates
            parameters(4,1) = rand; % carrying capacity of M cells
            parameters(5,1) = rand; % retention constant of M population
            parameters(6,1) = rand; % retention constant of M population
            parameters(7,1) = rand; % growth rate of e cells
            parameters(8,1) = rand; % non-dimensional rescaled M to E transiton rates
            parameters(9,1) = rand; % non-dimensional rescaled E to M transiton rates
            parameters(10,1) = rand; % carrying capacity of M cells
            parameters(11,1) = rand; % retention constant of M population
            parameters(12,1) = rand; % retention constant of M population
            parameters(13,1) = rand; % retention constant of M population

            lower_bound = [0 0 0 0 0 0 0 0 0 0 0 0 0];
            upper_bound = [1 1 1 1 1 1 1 1 1 1 1 1 1];

        case {29}
            parameters = zeros(5,1);
            parameters(1,1) = (1 + 2*rand); % growth rate of e cells
            parameters(2,1) = rand; % E competition to M cells
            parameters(3,1) = rand; % M competition to E cells
            parameters(4,1) = 10000000; % carrying capacity of M cells
            parameters(5,1) = 10000000; % carrying capacity of E cells


            lower_bound = [1 0 0  1000000 1000000 ];
            upper_bound = [5 10 10 20000000 20000000];

    end
    
    
    if(search_method == 1)
        cost = @(parameters) obj_fun(parameters, data_type,data,pop_model_indx, num_rep,timept,interpolated_timept, interpolation);
        options = optimoptions('lsqnonlin','Display','off');
        
        [temp_select_para, temp_goodness] = lsqnonlin(cost,parameters,lower_bound, upper_bound, options); % lsqnonlin is a non-linear optimizer; sim_traj: simulate trajectory
    else
        temp_select_para = parameters;
        %                 ode_sim_data = ODE_simulation(pop_model_indx, num_rep,timept,parameters);
        %                 temp_goodness = sum((data(:,2) - ode_sim_data(:,1)).^2)/(2*num_rep*num_time_pts);
        temp_goodness = 0;
    end
    
    % if(temp_select_para(1)<=temp_select_para(2) || temp_select_para(1)>=1) % To ascertain that epithelial cells divide faster than mesenchymal cell but this condition can be relaxed;
    % temp_select_para(1)<=temp_select_para(2) for models not scaled by r1
    % temp_select_para(1) >=1 for models scaled by r1
    %     switch pop_model_indx
    %         case {1,2,5}
    %             if(temp_select_para(1)<=temp_select_para(2)) % To ascertain that epithelial cells divide faster than mesenchymal cell but this condition can be relaxed;
    %                 select_para(itr_num,:) = temp_select_para';
    %                 goodness(itr_num) = temp_goodness;
    % %                 disp(['iteration number ' num2str(itr_num)]);
    %             end
    %         otherwise
    select_para(itr_num,:) = temp_select_para';
    goodness(itr_num) = temp_goodness;
    disp(['iteration number ' num2str(itr_num)]);
    %     end
    
end
end

