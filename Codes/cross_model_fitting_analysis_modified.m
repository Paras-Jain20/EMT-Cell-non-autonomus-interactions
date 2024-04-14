%%%%%%%%% main_function structural identifiability %%%%%%%%%
close all

% whether we are using the test or experimental data: test_data = 1; exp_data = 2
data_type = 1;

% type of optimum parameters search to be performed:
% optimizer based : choose 1
% uniform sampling based : choose 2
search_method = 1; % parameter search method
total_test_data = 1000;
total_itr = 1;
interpolation = 0;

r1 = 1;
% pop_model = [20 21 26 27 25 23];
% pop_model_names = ["G&T" , "GC&T", "G&T-Er" , "G&T-Mr", "G&T-EMr", "G&T-Mi-Er"];

pop_model= [3 4 11 12 14 10 13];
pop_model_names = ["G&T" , "GC&T", "G&T-Mr", "G&T-Er", "G&T-EMr", "G&T-Mi", "G&T-Mi-Er"];


bic_matrix = zeros(length(pop_model),length(pop_model),total_test_data); % for AIC BIC analysis
aic_matrix = zeros(length(pop_model),length(pop_model),total_test_data);
aicc_matrix = zeros(length(pop_model),length(pop_model),total_test_data);

birth_tran_ratio = 10;

for noise_factor = 5%[25 50]

    % select a population growth model
    for cross_pop_model_indx = 1:length(pop_model)


        switch pop_model(cross_pop_model_indx)
            case{6} % i.e. population model with 1 parameter
                num_model_para = 1;
            case{19} % i.e. population model with 1 parameter
                num_model_para = 2;
            case {1,2,10,11,12,20} % i.e. population model with 4 parameters
                num_model_para = 4;
            case {3,18} % i.e. population model with 3 parameters
                num_model_para = 3;
            case {4,5,7,13,14,15,16,17,29, 26,27} % i.e. population model with 5 parameters
                num_model_para = 5;
            case {8,9,22,23,24,25} % i.e. population model with 6 parameters
                num_model_para = 6;
            case {21} % i.e. population model with 7 parameters
                num_model_para = 7;
        end

        % data charactersitics
        if(pop_model(cross_pop_model_indx) < 15) % for models capturing dynamics of Bhatia et al. data
            num_rep = 3;
            init_cond = [1 2]; % initial condition vector
            num_init_cond = length(init_cond);
            r1 = 1;
            timept = 1:8; % time is unitless as it is scaled by m cells division rate (r1)
            interpolated_timept = 1; % interpolating time to get more data points to plot ODE solution trajectories
        else
            num_rep = 3;
            init_cond = [1 2 3 4 5]; % initial condition vector
            num_init_cond = length(init_cond);
            r1 = 1;
            timept = (0:3:12); % time is unitless as it is scaled by m cells division rate (r1)
            interpolated_timept = 1;  % interpolating time to get more data points to plot ODE solution trajectories
        end

        num_time_pts = length(timept); % number of time pts per replicate

        temp_bic_matrix = zeros(length(pop_model),total_test_data); % for AIC BIC analysis
        temp_aic_matrix = zeros(length(pop_model),total_test_data);
        temp_aicc_matrix = zeros(length(pop_model),total_test_data);

        combined_cross_model_data = zeros(total_test_data*length(pop_model),num_model_para+3); % for a cross pop model

        cd('C:\Users\user\OneDrive - Indian Institute of Science\GitHub\EMT-State-Transition-2.0\Parameter sets and generated pseudo data')

        if(pop_model(cross_pop_model_indx) > 15 && isfile(['Cross_model_fitting_data_Yamamoto_models_cross_model_' num2str(pop_model(cross_pop_model_indx)) '_b_t_ratio_' num2str(birth_tran_ratio) '_noise_factor ' num2str(noise_factor) '.xlsx'])...
                || pop_model(cross_pop_model_indx) < 15 && isfile(['Cross_model_fitting_data_Bhatia_models_cross_model_' num2str(pop_model(cross_pop_model_indx)) '_b_t_ratio_' num2str(birth_tran_ratio) '_noise_factor ' num2str(noise_factor) '.xlsx']))
            file_present = true;
            if(pop_model(cross_pop_model_indx) > 15)
                data = readmatrix(['Cross_model_fitting_data_Yamamoto_models_cross_model_' num2str(pop_model(cross_pop_model_indx)) '_b_t_ratio_' num2str(birth_tran_ratio) '_noise_factor ' num2str(noise_factor) '.xlsx']);

            else
                data = readmatrix(['Cross_model_fitting_data_Bhatia_models_cross_model_' num2str(pop_model(cross_pop_model_indx)) '_b_t_ratio_' num2str(birth_tran_ratio) '_noise_factor ' num2str(noise_factor) '.xlsx']);
            end
            combined_cross_model_data(1:size(data,1),1:size(data,2)) = data;
            empty_row_indx = size(data,1) + 1;
        else
            file_present = false;
            empty_row_indx = 1;
        end

        for pop_model_indx = 1:length(pop_model)
            for test_data_indx = 1:total_test_data
                disp(['Running cross fitting analysis for pop_model ' num2str(pop_model(pop_model_indx)) ' cross_pop_model ' num2str(pop_model(cross_pop_model_indx)) ' test data ' num2str(test_data_indx) ' noise factor ' num2str(noise_factor)])

                if(file_present)

                    cross_fitting_row_indx = find(sum(data(:,end-1:end) == [test_data_indx pop_model(pop_model_indx)],2) == 2,1);

                    if(isempty(cross_fitting_row_indx))

                        cd('C:\Users\user\OneDrive - Indian Institute of Science\GitHub\EMT-State-Transition-2.0\Parameter sets and generated pseudo data')
                        test_data = readmatrix(['test_data_' num2str(test_data_indx) '_pop_mod_' num2str(pop_model(pop_model_indx)) '_b_t_ratio_' num2str(birth_tran_ratio) '_noise_factor_' num2str(noise_factor) '.csv']);

                        cd('C:\Users\user\OneDrive - Indian Institute of Science\GitHub\EMT-State-Transition-2.0\Code files')
                        [select_para,goodness] = para_search(search_method,pop_model(cross_pop_model_indx),total_itr,data_type, test_data, num_rep,init_cond, timept, interpolated_timept, interpolation, num_time_pts,birth_tran_ratio); % select_para are filtered for re > rm

                        cross_model_data  = [select_para,goodness];
                        combined_cross_model_data(empty_row_indx,1:end-2) = cross_model_data(1,:);
                        combined_cross_model_data(empty_row_indx,end-1) = test_data_indx;
                        combined_cross_model_data(empty_row_indx,end) = pop_model(pop_model_indx);

                        empty_row_indx = empty_row_indx + 1;

                    else
                        cross_model_data =  data(cross_fitting_row_indx,:);
                    end
                else
                    cd('C:\Users\user\OneDrive - Indian Institute of Science\GitHub\EMT-State-Transition-2.0\Parameter sets and generated pseudo data')
                    test_data = readmatrix(['test_data_' num2str(test_data_indx) '_pop_mod_' num2str(pop_model(pop_model_indx)) '_b_t_ratio_' num2str(birth_tran_ratio) '_noise_factor_' num2str(noise_factor) '.csv']);

                    cd('C:\Users\user\OneDrive - Indian Institute of Science\GitHub\EMT-State-Transition-2.0\Code files')
                    [select_para,goodness] = para_search(search_method,pop_model(cross_pop_model_indx),total_itr,data_type, test_data, num_rep,init_cond, timept, interpolated_timept, interpolation, num_time_pts,birth_tran_ratio); % select_para are filtered for re > rm

                    cross_model_data  = [select_para,goodness];
                    combined_cross_model_data(empty_row_indx,1:end-2) = cross_model_data(1,:);
                    combined_cross_model_data(empty_row_indx,end-1) = test_data_indx;
                    combined_cross_model_data(empty_row_indx,end) = pop_model(pop_model_indx);

                    empty_row_indx = empty_row_indx + 1;

                end


                cross_model_goodness = (-1/2)*min(cross_model_data(:,end-2));


                [aic,bic] = aicbic(cross_model_goodness,num_model_para,(num_time_pts*num_rep*num_init_cond));
                temp_aic_matrix( pop_model_indx, test_data_indx) = aic;
                temp_bic_matrix( pop_model_indx, test_data_indx) = bic;
                temp_aicc_matrix( pop_model_indx, test_data_indx) = aic + (2*num_model_para*(num_model_para + 1))/((num_time_pts*num_rep*num_init_cond) - num_model_para - 1);


            end
        end

        aic_matrix(:,cross_pop_model_indx,:) = temp_aic_matrix;
        bic_matrix(:,cross_pop_model_indx, :) = temp_bic_matrix;
        aicc_matrix(:,cross_pop_model_indx, :) = temp_aicc_matrix;

        if((file_present == true && size(data,1) < total_test_data*length(pop_model)) || file_present == false)
            cd('C:\Users\user\OneDrive - Indian Institute of Science\GitHub\EMT-State-Transition-2.0\Parameter sets and generated pseudo data')

            if(~isempty(find(pop_model > 15,1)))
                writematrix(combined_cross_model_data, ['Cross_model_fitting_data_Yamamoto_models_cross_model_' num2str(pop_model(cross_pop_model_indx)) '_b_t_ratio_' num2str(birth_tran_ratio) '_noise_factor ' num2str(noise_factor) '.xlsx']);
            elseif (~isempty(find(pop_model < 15,1)))
                writematrix(combined_cross_model_data, ['Cross_model_fitting_data_Bhatia_models_cross_model_' num2str(pop_model(cross_pop_model_indx)) '_b_t_ratio_' num2str(birth_tran_ratio) '_noise_factor ' num2str(noise_factor) '.xlsx']);
            end
        end
    end


    % cd('C:\Users\user\OneDrive - Indian Institute of Science\GitHub\EMT-State-Transition-2.0\Parameter sets and generated pseudo data')
    % if(~isempty(find(pop_model > 15,1)))
    %     save(['Cross model fitting Yamamoto models information criteria matrix noise factor ' num2str(noise_factor) ],'aic_matrix', 'bic_matrix','aicc_matrix');
    % elseif (~isempty(find(pop_model< 15,1)))
    %     save(['Cross model fitting Bhatia models information criteria matrix noise factor ' num2str(noise_factor)],'aic_matrix', 'bic_matrix','aicc_matrix');
    % end

    close all
    
    data_model_accuracy_aicc = zeros(length(pop_model),length(pop_model));
    data_model_accuracy_bic = zeros(length(pop_model),length(pop_model));
    data_model_accuracy_aic = zeros(length(pop_model),length(pop_model));

    for test_data_indx = 1:size(aicc_matrix,3)
        for row_indx = 1:size(aicc_matrix,1)
    
            [~,indx] = min(aicc_matrix(row_indx,:,test_data_indx));
            data_model_accuracy_aicc(row_indx, indx) = data_model_accuracy_aicc(row_indx, indx) + 1;
    
            [~,indx] = min(bic_matrix(row_indx,:,test_data_indx));
            data_model_accuracy_bic(row_indx, indx) = data_model_accuracy_bic(row_indx, indx) + 1;
    
            [~,indx] = min(aic_matrix(row_indx,:,test_data_indx));
            data_model_accuracy_aic(row_indx, indx) = data_model_accuracy_aic(row_indx, indx) + 1;
        end
    end
    
%     % rearranging rows and columns required for Yamamoto models analysis
%     data_model_accuracy_aic = data_model_accuracy_aic(:,[1 5 2 3 4 6]);
%     data_model_accuracy_aic = data_model_accuracy_aic([1 5 2 3 4 6],:);
%     
%     data_model_accuracy_aicc = data_model_accuracy_aicc(:,[1 5 2 3 4 6]);
%     data_model_accuracy_aicc = data_model_accuracy_aicc([1 5 2 3 4 6],:);
%     
%     data_model_accuracy_bic = data_model_accuracy_bic(:,[1 5 2 3 4 6]);
%     data_model_accuracy_bic = data_model_accuracy_bic([1 5 2 3 4 6],:);
    
    data_model_accuracy_aic  = 100*data_model_accuracy_aic./total_test_data;
    data_model_accuracy_bic  = 100*data_model_accuracy_bic./total_test_data;
    
    data_model_accuracy_aicc  = 100*data_model_accuracy_aicc./total_test_data;
    
    cd('C:\Users\user\OneDrive - Indian Institute of Science\GitHub\EMT-State-Transition-2.0\Manuscript figures\Figure 2 and related SI')
    
    figure(1)
    h1 = heatmap(categorical(pop_model_names),categorical(pop_model_names),data_model_accuracy_aicc);
    h1.XLabel = 'cross model';
    h1.YLabel = 'true model';
    colorbar off
    h1.Title = 'Percentage preference as per AICc';
    h1.ColorLimits = [0 100];
    h1.FontSize = 14;
    
    if(~isempty(find(pop_model> 15,1)))
        print(gcf, ['Cross Yamamoto model fitting analysis AICc model selection criteria noise factor ' num2str(noise_factor) '.png'],'-dpng','-r300');
    elseif (~isempty(find(pop_model < 15,1)))
        print(gcf, ['Cross Bhatia model fitting analysis AICc model selection criteria noise factor ' num2str(noise_factor) '.png'],'-dpng','-r300');
    end
    
    
    figure(2)
    h1 = heatmap(categorical(pop_model_names),categorical(pop_model_names),data_model_accuracy_bic);
    h1.XLabel = 'cross model';
    h1.YLabel = 'true model';
    colorbar off
    h1.Title = 'Percentage preference as per BIC';
    h1.ColorLimits = [0 100];
    h1.FontSize = 14;
    
    if(~isempty(find(pop_model > 15,1)))
        print(gcf, ['Cross Yamamoto model fitting analysis BIC model selection criteria noise factor ' num2str(noise_factor) '.png'],'-dpng','-r300');
    elseif (~isempty(find(pop_model < 15,1)))
        print(gcf, ['Cross Bhatia model fitting analysis BIC model selection criteria noise factor ' num2str(noise_factor) '.png'],'-dpng','-r300');
    end
    
    
    figure(3)
    h1 = heatmap(categorical(pop_model_names),categorical(pop_model_names),data_model_accuracy_aic);
    h1.XLabel = 'cross model';
    h1.YLabel = 'true model';
    colorbar off
    h1.Title = 'Percentage preference as per AIC';
    h1.ColorLimits = [0 100];
    h1.FontSize = 14;
    
    if(~isempty(find(pop_model > 15,1)))
        print(gcf, ['Cross Yamamoto model fitting analysis AIC model selection criteria noise factor ' num2str(noise_factor) '.png'],'-dpng','-r300');
    elseif (~isempty(find(pop_model < 15,1)))
        print(gcf, ['Cross Bhatia model fitting analysis AIC model selection criteria noise factor ' num2str(noise_factor) '.png'],'-dpng','-r300');
    end
end