
% Code file to compare minimum chi-square values for population models
% across b-t-ratio

close all

% whether we are using the test or experimental data: test_data = 1; exp_data = 2
data_type = 2;

% type of optimum parameters search to be performed:
% use existing parameters = choose 0
% optimizer based : choose 1
% uniform sampling based : choose 2; this method is not updated so might not work
search_method = 1; % parameter search method

test_data_indx = 1; % necessary dummy variable to just to fill in variable list in profile likelihood function
noise_factor = 1; % necessary dummy variable to just to fill in variable list in profile likelihood function
r1 = [1/54]; % m cells division rates; scaling factor for non dimensional population model
inverse_r2 =  35;% r_e
r2 =  1/(inverse_r2*r1(1)); % normalised r_e; is set to 1.5 which is normalised with r1; if r2 is empty set then it is a variable parameter

default_total_itr =3;

birth_tran_ratio_set = 150;

Pop_model_names = ["GC_{s}&T" , "GC_{a}&T", "GC_{s}&T-Er", "GC_{s}&T-Mr" , "GC_{s}&T-EMr",  "GC_{s}&T-Mi-Er", "GC_{s}I&T" ];
% Pop_model_names = ["GC_{s}&T" , "GC_{a}&T", "GC_{s}&T-Mi-Er",  "GC_{s}&T-EMr", "GC_{s}&T-Er", "GC_{s}&T-Mr" ,   "GC_{s}I&T" ];

goodness_matrix = zeros(7,3);

pop_model = [34 20 20.1 35 21 21.1 38 26 26.1 39 27 27.1 37 25 25.1 36 23 23.1 40 33 33.1];%

i = 1; % to track population model
% select a population growth model
for pop_model_indx = pop_model%1:max(pop_model)%[20 26 27 25 23];%[11 14]%[3 11 12 14 10 13]%[1 2 3 10 11 12 13]
    temp_pop_model_indx = pop_model_indx; 
    disp(temp_pop_model_indx)
    j = 1; % to track b t ratio model
    %     if(~isempty(find(temp_pop_model_indx == pop_model,1)))

    switch temp_pop_model_indx
        case {34 35 36 37 38 39 40}
            r1 = 1; % m cells division rates; scaling factor for non dimensional population model
            r2 =  [];% 1/(inverse_r2*r1(1)); % normalised r_e; is set to 1.5 which is normalised with r1; if r2 is empty set then it is a variable parameter
        case {20 21 23 25 26 27 33}
            r1 = 1/54; % m cells division rates; scaling factor for non dimensional population model
            r2 =  [];% 1/(inverse_r2*r1(1)); % normalised r_e; is set to 1.5 which is normalised with r1; if r2 is empty set then it is a variable parameter
        case {20.1 21.1 23.1 25.1 26.1 27.1 33.1}
            temp_pop_model_indx = floor(temp_pop_model_indx);
            r1 = 1/54; % m cells division rates; scaling factor for non dimensional population model
            r2 = 1/(inverse_r2*r1(1)); % normalised r_e; is set to 1.5 which is normalised with r1; if r2 is empty set then it is a variable parameter
    end

    for birth_tran_ratio = birth_tran_ratio_set % max ratio between birth and transition
        close all

        if(temp_pop_model_indx < 15) % for models capturing dynamics of Bhatia et al. data

            cd('C:\Users\user\OneDrive - Indian Institute of Science\GitHub\EMT-State-Transition-2.0\Parameters sets from experimental data')

            if(isfile(['para_Bhatia_data_pd_mod_' num2str(temp_pop_model_indx) '_s_m_' num2str(search_method) '_b_t_ratio_' num2str(birth_tran_ratio) '_inverse_r1_' num2str(1/r1) '.csv']))

                para_data = table2array(readtable(['para_Bhatia_data_pd_mod_' num2str(temp_pop_model_indx) '_s_m_' num2str(search_method) '_b_t_ratio_' num2str(birth_tran_ratio) '_inverse_r1_' num2str(1/r1) '.csv']));
                num_para_set = size(para_data,1);
                if(num_para_set >= default_total_itr)
                    select_para = (para_data(:,1:end-1));
                    goodness = (para_data(:,end));
                    para_opt = false; % to perform parameter search using optimizer or not
                    total_itr = 0;
                else

                    para_opt = true;
                    total_itr = default_total_itr-num_para_set; % max 100 iterations have to be performed
                end
            else
                para_opt = true; % if the data file is not present then perform the para search
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
            num_rep = 3;
            num_init_cond = 2; % number of different initial condition in the data
            init_cond = 1:num_init_cond;
            timept = (2:2:8)*(7*24)*r1; % time is unitless as it is scaled by m cells division rate (r1)
            %                 interpolated_timept = r1:r1:timept(end); % interpolating time to get more data points to plot ODE solution trajectories
            interpolated_timept = r1:r1:1000; % interpolating time to get more data points to plot ODE solution trajectories

        else
            cd('C:\Users\user\OneDrive - Indian Institute of Science\GitHub\EMT-State-Transition-2.0\Parameters sets from experimental data')
            if(~isempty(r2) && r1 ~= 1 && isfile(['para_Yamamoto_data_pd_mod_' num2str(temp_pop_model_indx) '_s_m_' num2str(search_method) '_b_t_ratio_' num2str(birth_tran_ratio) '_inverse_r1_' num2str(1/r1) '_inverse_r2_' num2str(inverse_r2) '.csv'])...
                    || isempty(r2) && r1 ~= 1 && isfile(['para_Yamamoto_data_pd_mod_' num2str(temp_pop_model_indx) '_s_m_' num2str(search_method) '_b_t_ratio_' num2str(birth_tran_ratio) '_inverse_r1_' num2str(1/r1) '.csv']) ...
                    || isempty(r2) && r1 == 1 && isfile(['para_Yamamoto_data_pd_mod_' num2str(temp_pop_model_indx) '_s_m_' num2str(search_method) '_b_t_ratio_' num2str(birth_tran_ratio) '.csv']))

                if(~isempty(r2) && r1 ~= 1)
                    para_data = table2array(readtable(['para_Yamamoto_data_pd_mod_' num2str(temp_pop_model_indx) '_s_m_' num2str(search_method) '_b_t_ratio_' num2str(birth_tran_ratio) '_inverse_r1_' num2str(1/r1) '_inverse_r2_' num2str(inverse_r2) '.csv']));
                elseif(isempty(r2) && r1 ~= 1)
                    para_data = table2array(readtable(['para_Yamamoto_data_pd_mod_' num2str(temp_pop_model_indx) '_s_m_' num2str(search_method) '_b_t_ratio_' num2str(birth_tran_ratio) '_inverse_r1_' num2str(1/r1) '.csv']));
                elseif(isempty(r2) && r1 == 1)
                    para_data = table2array(readtable(['para_Yamamoto_data_pd_mod_' num2str(temp_pop_model_indx) '_s_m_' num2str(search_method) '_b_t_ratio_' num2str(birth_tran_ratio)  '.csv']));
                end

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

            %                 Importing Yamamoto et al. 2017 data; models having
            %                 pop_model_indx >= 15
            cd('C:\Users\user\OneDrive - Indian Institute of Science\Projects\State transition\Experimental data\Yamamoto et al. 2017')
            data = (readmatrix('Yamamoto et al 2017 Data 1 cell fraction.xlsx'));
            num_rep = 3;
            num_init_cond = 5; % number of different initial condition in the data
            init_cond = 1:num_init_cond;
            timept = (0:3:12)*24*r1; % time is unitless as it is scaled by m cells division rate (r1)

            interpolated_timept = 0:r1:timept(end);  % interpolating time to get more data points to plot ODE solution trajectories
        end

        num_time_pts = length(timept);

        cd('C:\Users\user\OneDrive - Indian Institute of Science\GitHub\EMT-State-Transition-2.0\Code files')

        if(para_opt)
            interpolation = 0; % to add more time points for the ODE simulated data; must be zero if we using objective function to minimize error with the data

            % parameter search function
            [select_para,goodness] = para_search(search_method,temp_pop_model_indx,total_itr, data_type, data, num_rep, init_cond, timept, interpolated_timept, interpolation, num_time_pts,birth_tran_ratio,r2); % select_para are filtered for re > rm
            row_indx = floor((i-1)/3)+1;
            column_indx = mod(i-1,3)+1;
            goodness_matrix(row_indx,column_indx) = min(goodness);

            cd('C:\Users\user\OneDrive - Indian Institute of Science\GitHub\EMT-State-Transition-2.0\Parameters sets from experimental data')
           
            %                 if(pop_model_indx <15 )
            %                     writematrix(['para_mod_' num2str(pop_model_indx) '_s_m_' num2str(search_method) '_b_t_ratio_' num2str(birth_tran_ratio) '_inverse_r1_' num2str(1/r1) '.csv'],[select_para,goodness],'-append');
            %                 else
            %                     writematrix([select_para,goodness],['para_pd_mod_' num2str(pop_model_indx) '_s_m_' num2str(search_method) '_b_t_ratio_' num2str(birth_tran_ratio) '_inverse_r1_' num2str(1/r1) '_inverse_r2_' num2str(inverse_r2) '.csv']);
            %                 end

            % if(~isempty(r2) && r1 ~= 1)
            %     para_data = writematrix(['para_Yamamoto_data_pd_mod_' num2str(temp_pop_model_indx) '_s_m_' num2str(search_method) '_b_t_ratio_' num2str(birth_tran_ratio) '_inverse_r1_' num2str(1/r1) '_inverse_r2_' num2str(inverse_r2) '.csv'],[select_para,goodness],'-append');
            % elseif(isempty(r2) && r1 ~= 1)
            %     para_data = writematrix(['para_Yamamoto_data_pd_mod_' num2str(temp_pop_model_indx) '_s_m_' num2str(search_method) '_b_t_ratio_' num2str(birth_tran_ratio) '_inverse_r1_' num2str(1/r1) '.csv'],[select_para,goodness],'-append');
            % elseif(isempty(r2) && r1 == 1)
            %     para_data = writematrix(['para_Yamamoto_data_pd_mod_' num2str(temp_pop_model_indx) '_s_m_' num2str(search_method) '_b_t_ratio_' num2str(birth_tran_ratio)  '.csv'],[select_para,goodness],'-append');
            % end
        end
        j = j + 1;
    end
    i = i + 1;
    %     end
end


figure_position = [0.4630 0.2639 0.4791 0.4648];

norm_goodness_matrix = goodness_matrix/(length(timept)*num_rep*num_init_cond);
% to plot the minimum chi-sqr value against pop model for various bt ratios
figure('Units','normalized','Position',figure_position)
X = categorical(Pop_model_names);
X = reordercats(X,Pop_model_names);
heatmap(X,{'none', 'r1', 'r1 r2'},norm_goodness_matrix',"ColorScaling","log",Title="Normalised Minimum Chi-square");
ax = gca;
grid on
xlabel('Population model')
ylabel('Scenario')
ax.FontSize = 15;
colormap parula
colorbar off
saveas(gcf,['Yamamoto_data_Normalised_Min_chi_sqr_vs_pop_model_with_without_const_r_s_m_' num2str(search_method) '.png'])


% determination of number of parameter for different models
pop_model = [20 21 26 27 25 23 33];
num_model_parameters = ones(size(goodness_matrix));
i = 1; % to track population model

for pop_model_indx = pop_model
    % below is the number of parameters when both r1 and r2 are defined and
    % constant
    for j = 1:3
    switch pop_model_indx
        case {4,5,13,14,15,16,17,19,29,26,27} % i.e. population model with 4 parameters
            num_para = 4;
        case {1,2,10,11,12,20} % i.e. population model with 3 parameters
            num_para = 3;
        case {8,9,22,23,24,25,33} % i.e. population model with 5 parameters
            num_para = 5;
        case {21} % i.e. population model with 6 parameters
            num_para = 6;
    end
    num_model_parameters(i, j) = num_para + 3-j;
    end
    i = i+1;
end


% calculation of p-value of the chi-square statistics - estimate of
% goodness of fit of the model
p_value = chi2cdf(goodness_matrix,(length(timept)*num_rep*num_init_cond) - num_model_parameters,'upper');

figure('Units','normalized','Position',figure_position)
X = categorical(Pop_model_names);
X = reordercats(X,Pop_model_names);
heatmap(X,{'none', 'r1', 'r1 r2'},log10(p_value'),"ColorScaling","log",Title="log_{10}(p-value chi-sqr statistic)");
ax = gca;
grid on
xlabel('Population model')
ylabel('Scenario')
ax.FontSize = 15;
colormap parula
colorbar off
saveas(gcf,['Yamamoto_data_p_statistic_vs_pop_model_with_without_const_r_s_m_' num2str(search_method)  '.png'])

[aic,bic] = aicbic(reshape(-1/2*goodness_matrix,1,[]),reshape(num_model_parameters,1,[]),reshape(ones(size(goodness_matrix)).*(num_time_pts*num_rep*num_init_cond),1,[]));
aicc = aic + (2*reshape(num_model_parameters,1,[]).*(reshape(num_model_parameters,1,[])+1))./((num_time_pts*num_rep*num_init_cond) - reshape(num_model_parameters,1,[]) - 1);

aicc  = reshape(aicc,size(goodness_matrix));
figure('Units','normalized','Position',figure_position)
X = categorical(Pop_model_names);
X = reordercats(X,Pop_model_names);
heatmap(X,{'none', 'r1', 'r1 r2'},aicc',"ColorScaling","log",Title="AICc");
ax = gca;
grid on
xlabel('Population model')
ylabel('Scenario')
ax.FontSize = 15;
colormap parula
colorbar off
saveas(gcf,['Yamamoto_data_AICc_vs_pop_model_with_without_const_r_s_m_' num2str(search_method) '.png'])


% to plot the minimum chi-sqr value against pop model for various bt ratios
figure('Units','normalized','Position',figure_position)
X = categorical(Pop_model_names);
X = reordercats(X,Pop_model_names);
heatmap(X,{'none', 'r1', 'r1 r2'},goodness_matrix',"ColorScaling","log",Title="Minimum Chi-square");
ax = gca;
grid on
xlabel('Population model')
ylabel('Scenario')
ax.FontSize = 15;
colormap parula
colorbar off
saveas(gcf,['Yamamoto_data_Min_chi_sqr_vs_pop_model_with_without_const_r_s_m_' num2str(search_method) '.png'])