
% Code file to compare minimum chi-square values for population models
% across b-t-ratio

close all

% whether we are using the test or experimental data: test_data = 1; exp_data = 2
data_type = 2;

test_data_indx = 1; % necessary dummy variable to just to fill in variable list in profile likelihood function
noise_factor = 1; % necessary dummy variable to just to fill in variable list in profile likelihood function
r1 = [1/54]; % m cells division rates; scaling factor for non dimensional population model
inverse_r2 =  35;% r_e
r2 =  1/(inverse_r2*r1); % normalised r_e; is set to 1.5 which is normalised with r1; if r2 is empty set then it is a variable parameter

default_total_itr = 10;
pop_model = [20 57 55 56 54];%
% pop_model = [3 28 52 51 53];
birth_tran_ratio = 200;%[50 100 150 200 250 300];
% Pop_model_names = ["GC_{s}&T" , "GC_{a}&T", "GC_{s}&T-Er", "GC_{s}&T-Mr" , "GC_{s}&T-EMr",  "GC_{s}&T-Mi-Er", "GC_{s}I&T" ];
Pop_model_names = ["T2" , "T3-HM", "T3-L-HM", "T3-EH" , "T3-L-EH"];
% Pop_model_names = ["T2" , "T3-HM", "T3-L-HM", "T3-EH" , "T3-L-EH"];

goodness_matrix = zeros(length(pop_model),1);

% select a population growth model
for pop_model_indx = 1:length(pop_model)
            disp(['Running the analysis for ' num2str(pop_model(pop_model_indx)) ' ' num2str(birth_tran_ratio)])
            close all

            % type of optimum parameters search to be performed:
            % use existing parameters = choose 0
            % optimizer based : choose 1
            % uniform sampling based : choose 2; this method is not updated so might not work
            search_method = 1; % parameter search method

            if(pop_model(pop_model_indx) < 15 || pop_model(pop_model_indx) == 28 || pop_model(pop_model_indx) == 51 || pop_model(pop_model_indx) == 52 || pop_model(pop_model_indx) == 53) % for models capturing dynamics of Bhatia et al. data

                cd('C:\Users\Asus\OneDrive - Indian Institute of Science\GitHub\EMT-State-Transition-2.0\Parameters sets from experimental data')

                if(isfile(['para_pd_mod_' num2str(pop_model(pop_model_indx)) '_s_m_' num2str(search_method) '_b_t_ratio_' num2str(birth_tran_ratio) '_inverse_r1_' num2str(1/r1) '.csv']))

                    para_data = table2array(readtable(['para_pd_mod_' num2str(pop_model(pop_model_indx)) '_s_m_' num2str(search_method) '_b_t_ratio_' num2str(birth_tran_ratio) '_inverse_r1_' num2str(1/r1) '.csv']));
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

                cd('C:\Users\Asus\OneDrive - Indian Institute of Science\Projects\State transition\Experimental data');

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
                cd('C:\Users\Asus\OneDrive - Indian Institute of Science\GitHub\EMT-State-Transition-2.0\Parameters sets from experimental data')
                if(isfile(['para_pd_mod_' num2str(pop_model(pop_model_indx)) '_s_m_' num2str(search_method) '_b_t_ratio_' num2str(birth_tran_ratio) '_inverse_r1_' num2str(1/r1) '_inverse_r2_' num2str(inverse_r2) '.csv']))

                    para_data = table2array(readtable(['para_pd_mod_' num2str(pop_model(pop_model_indx)) '_s_m_' num2str(search_method) '_b_t_ratio_' num2str(birth_tran_ratio) '_inverse_r1_' num2str(1/r1) '_inverse_r2_' num2str(inverse_r2) '.csv']));
    
                    select_para = (para_data(:,1:end-1));
                    goodness = (para_data(:,end));

                    num_para_set = size(para_data,1);

                    if(num_para_set >= default_total_itr)
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
                %                 pop_model(pop_model_indx) >= 15
                cd('C:\Users\Asus\OneDrive - Indian Institute of Science\Projects\State transition\Experimental data\Yamamoto et al. 2017')
                data = (readmatrix('Yamamoto et al 2017 Data 1 cell fraction.xlsx'));
                num_rep = 3;
                num_init_cond = 5; % number of different initial condition in the data
                init_cond = 1:num_init_cond;
                timept = (0:3:12)*24*r1; % time is unitless as it is scaled by m cells division rate (r1)
                interpolated_timept = 0:r1:timept(end);  % interpolating time to get more data points to plot ODE solution trajectories
            end

            num_time_pts = length(timept);

            cd('C:\Users\Asus\OneDrive - Indian Institute of Science\GitHub\EMT-State-Transition-2.0\Code files')

            if(para_opt)
                interpolation = 0; % to add more time points for the ODE simulated data; must be zero if we using objective function to minimize error with the data

                % parameter search function
                [temp_select_para,temp_goodness] = para_search(search_method,pop_model(pop_model_indx),total_itr, data_type, data, num_rep, init_cond, timept, interpolated_timept, interpolation, num_time_pts,birth_tran_ratio,r2); % select_para are filtered for re > rm

                cd('C:\Users\Asus\OneDrive - Indian Institute of Science\GitHub\EMT-State-Transition-2.0\Parameters sets from experimental data')

                if(pop_model(pop_model_indx) <15 || pop_model(pop_model_indx) == 28 || pop_model(pop_model_indx) == 51 || pop_model(pop_model_indx) == 52 || pop_model(pop_model_indx) == 53 )
                    writematrix([temp_select_para,temp_goodness],['para_pd_mod_' num2str(pop_model(pop_model_indx)) '_s_m_' num2str(search_method) '_b_t_ratio_' num2str(birth_tran_ratio) '_inverse_r1_' num2str(1/r1) '.csv'],'WriteMode', 'append');
                else
                    writematrix([temp_select_para,temp_goodness],['para_pd_mod_' num2str(pop_model(pop_model_indx)) '_s_m_' num2str(search_method) '_b_t_ratio_' num2str(birth_tran_ratio) '_inverse_r1_' num2str(1/r1) '_inverse_r2_' num2str(inverse_r2) '.csv'],'WriteMode','append');
                end

            else
                temp_goodness = [];
            end
            
            goodness = [goodness;temp_goodness];

            goodness_matrix(pop_model_indx) = min(goodness);

end

cd('C:\Users\Asus\OneDrive - Indian Institute of Science\GitHub\EMT-State-Transition-2.0\Manuscript documents\Revision\Figures')

figure_position = [0.4630 0.2639 0.4791 0.4648];

% to plot the minimum chi-sqr value against pop model for various bt ratios
figure('Units','normalized','Position',figure_position)
X = categorical(Pop_model_names);
X = reordercats(X,Pop_model_names);
% heatmap(X,birth_tran_ratio_set,norm_goodness_matrix',"ColorScaling","log",Title="Normalised Minimum Chi-square");
% bar(X,log10(goodness_matrix))
bar(X,(goodness_matrix))
ax = gca;
axis square
grid on
xlabel('Population model')
ylabel('Chi-square')
ax.FontSize = 16;
% colormap parula
% colorbar off
% saveas(gcf,['Yamamoto_data_Normalised_Min_chi_sqr_vs_pop_model_various_b_t_ratio_s_m_' num2str(search_method)  '_inverse_r1_' num2str(1/r1) '_inverse_r2_' num2str(inverse_r2) '.png'])
% saveas(gcf,['Bhatia_data_Min_chi_sqr_vs_pop_model_G&T_two_three_states_s_m_' num2str(search_method)  '_inverse_r1_' num2str(1/r1) '.png'])
saveas(gcf,['Yamamoto_data_Min_chi_sqr_vs_pop_model_G&T_two_three_states_s_m_' num2str(search_method)  '_inverse_r1_' num2str(1/r1) '_inverse_r2_' num2str(inverse_r2) '.png'])


figure_position = [0.4630 0.2639 0.4791 0.4648];

norm_goodness_matrix = goodness_matrix/(length(timept)*num_rep*num_init_cond);
% to plot the minimum chi-sqr value against pop model for various bt ratios
figure('Units','normalized','Position',figure_position)
X = categorical(Pop_model_names);
X = reordercats(X,Pop_model_names);
% heatmap(X,birth_tran_ratio_set,norm_goodness_matrix',"ColorScaling","log",Title="Normalised Minimum Chi-square");
% bar(X,log10(goodness_matrix))
bar(X,(norm_goodness_matrix))
ax = gca;
axis square
grid on
xlabel('Population model')
ylabel('Normalised Chi-square')
ax.FontSize = 16;
% colormap parula
% colorbar off
% saveas(gcf,['Yamamoto_data_Normalised_Min_chi_sqr_vs_pop_model_various_b_t_ratio_s_m_' num2str(search_method)  '_inverse_r1_' num2str(1/r1) '_inverse_r2_' num2str(inverse_r2) '.png'])
% saveas(gcf,['Bhatia_data_Normalised_Min_chi_sqr_vs_pop_model_G&T_two_three_states_s_m_' num2str(search_method)  '_inverse_r1_' num2str(1/r1) '.png'])
saveas(gcf,['Yamamoto_data_Normalised_Min_chi_sqr_vs_pop_model_G&T_two_three_states_s_m_' num2str(search_method)  '_inverse_r1_' num2str(1/r1) '_inverse_r2_' num2str(inverse_r2) '.png'])


% determination of number of parameter for different models
num_model_parameters = ones(size(goodness_matrix));
i = 1; % to track population model
for pop_model_indx = pop_model
    if(isempty(r2))
        switch pop_model_indx
            case{6} % i.e. population model with 1 parameter
                delta_chi_sq_all = 4;
                num_para = 1;
            case {1,2,10,11,12,20} % i.e. population model with 4 parameters
                delta_chi_sq_all = 9.70;
                num_para = 4;
            case {3,18} % i.e. population model with 3 parameters
                delta_chi_sq_all = 8.02;
                num_para = 3;
            case {4,5,7,13,14,15,16,17,19, 26,27,4.1} % i.e. population model with 5 parameters
                delta_chi_sq_all = 11.3;
                num_para = 5;
            case {8,9,22,23,24,25,33,52,53} % i.e. population model with 6 parameters
                delta_chi_sq_all = 12.8;
                num_para = 6;
            case {21} % i.e. population model with 7 parameters
                delta_chi_sq_all = 14.067;
                num_para = 7;
            case{28,51}
                num_para = 8;
        end
    else

        switch pop_model_indx
            case{6} % i.e. population model with 1 parameter
                delta_chi_sq_all = 4;
            case{3} % i.e. population model with 2 parameter
                delta_chi_sq_all = 6.17;
                num_para = 2;
            case {4,5,13,14,15,16,17,19,29,26,27} % i.e. population model with 4 parameters
                delta_chi_sq_all = 9.70;
                num_para = 4;
            case {1,2,10,11,12,20} % i.e. population model with 3 parameters
                delta_chi_sq_all = 8.02;
                num_para = 3;
            case {8,9,22,23,24,25,33} % i.e. population model with 5 parameters
                delta_chi_sq_all = 11.3;
                num_para = 5;
            case {21} % i.e. population model with 6 parameters
                delta_chi_sq_all = 12.8;
                num_para = 6;
            case {54,55}
                num_para = 6;
            case {56,57}
                num_para = 8;
        end

    end
    num_model_parameters(i,:) = num_para.*num_model_parameters(i,:);
    i = i+1;
end

% calculation of p-value of the chi-square statistics - estimate of
% goodness of fit of the model
p_value = chi2cdf(goodness_matrix,(length(timept)*num_rep*num_init_cond) - num_model_parameters,'upper');

figure('Units','normalized','Position',figure_position)
X = categorical(Pop_model_names);
X = reordercats(X,Pop_model_names);
% heatmap(X,birth_tran_ratio_set,log10(p_value'),"ColorScaling","log",Title="log_{10}(p-value chi-sqr statistic)");
bar(X,-log10(p_value'))
ax = gca;
grid on
xlabel('Population model')
ylabel('-log_{10}(p-value)')
ax.FontSize = 16;
axis square
colormap parula
colorbar off
% saveas(gcf,['Yamamoto_data_p_statistic_vs_pop_model_various_b_t_ratio_s_m_' num2str(search_method)  '_inverse_r1_' num2str(1/r1) '_inverse_r2_' num2str(inverse_r2) '.png'])
% saveas(gcf,['Bhatia_data_p_statistic_vs_pop_model_G&T_two_three_states_s_m_' num2str(search_method)  '_inverse_r1_' num2str(1/r1) '.png'])
saveas(gcf,['Yamamoto_data_p_statistic_vs_pop_model_G&T_two_three_states_s_m_' num2str(search_method)  '_inverse_r1_' num2str(1/r1) '_inverse_r2_' num2str(inverse_r2) '.png'])

[aic,bic] = aicbic(-1/2*goodness_matrix,num_model_parameters.*ones(size(goodness_matrix)),(num_time_pts*num_rep*num_init_cond).*ones(size(goodness_matrix)));
aicc = aic + 2*num_model_parameters.*(num_model_parameters + 1)./((num_time_pts*num_rep*num_init_cond) - num_model_parameters - 1);

figure('Units','normalized','Position',figure_position)
X = categorical(Pop_model_names);
X = reordercats(X,Pop_model_names);
% heatmap(X,birth_tran_ratio_set,log10(p_value'),"ColorScaling","log",Title="log_{10}(p-value chi-sqr statistic)");
bar(X,aicc)
ax = gca;
grid on
xlabel('Population model')
ylabel('AICc')
axis square
ax.FontSize = 16;
% saveas(gcf,['Yamamoto_data_p_statistic_vs_pop_model_various_b_t_ratio_s_m_' num2str(search_method)  '_inverse_r1_' num2str(1/r1) '_inverse_r2_' num2str(inverse_r2) '.png'])
% saveas(gcf,['Bhatia_data_AICc_vs_pop_model_G&T_two_three_states_s_m_' num2str(search_method)  '_inverse_r1_' num2str(1/r1) '.png'])
saveas(gcf,['Yamamoto_data_AICc_vs_pop_model_G&T_two_three_states_s_m_' num2str(search_method)  '_inverse_r1_' num2str(1/r1) '_inverse_r2_' num2str(inverse_r2) '.png'])