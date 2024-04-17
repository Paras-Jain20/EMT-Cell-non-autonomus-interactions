
% for Yamamoto data
pop_model = [20 21 26 27 25 23 33];
pop_model_names = ["GC_{s}&T" , "GC_{a}&T", "GC_{s}&T-Mr", "GC_{s}&T-Er" , "GC_{s}&T-EMr",  "GC_{s}&T-Mi-Er", "GC_{s}I&T" ];
num_rep = 3;
num_init = 2; % data generated before 15th July considered value 2 for this parameter 
num_timept = 5;

% for Bhatia data
% pop_model= [3 4 11 12 14 10 13];
% pop_model_names = ["G&T" , "GC&T", "G&T-Mr", "G&T-Er","G&T-EMr", "G&T-Mi", "G&T-Mi-Er"];
% num_rep = 3;
% num_init = 2;
% num_timept = 4;

mymap = [[0.8500 0.3250 0.0980]
        [0.9290 0.6940 0.1250]
        [0 0.4470 0.7410]];

num_pop_model = length(pop_model);
r1 = 1/54;
inverse_r2 = 35;
r2 = 1/(r1*inverse_r2);
data_src = {'Bhatia', 'Yamamoto'};
profile_sample_size = 2000;

% max_para = 5; % for Bhatia models
max_para = 6;% = 6 for r2 and r1 normalised; = 7 for r1 noramlised


for data_src_indx = 2%1:length(data_src)
    for birth_tran_ratio = 250%[5 10 20 50 100]%[10 50 100 500 1000]
        iden_mat = zeros(num_pop_model,max_para);
        actual_lower_bound = zeros(num_pop_model,max_para);
        actual_upper_bound = zeros(num_pop_model,max_para);
        struct_non_iden_mat = zeros(num_pop_model,max_para);
        struc_iden_mat = zeros(num_pop_model,max_para);
        prac_iden_mat = zeros(num_pop_model,max_para);
        goodness_mat = zeros(num_pop_model,max_para);
        close all
        for pop_model_indx = 1:length(pop_model)
                
            [lower_bound,upper_bound] = parameter_bounds(pop_model(pop_model_indx),birth_tran_ratio,r2);

                cd('C:\Users\Asus\OneDrive - Indian Institute of Science\GitHub\EMT-State-Transition-2.0\Parameters sets from experimental data')
%                 cd('C:\Users\Asus\OneDrive - Indian Institute of Science\GitHub\EMT-State-Transition-2.0\code files ram\data post sampling')

                
%                 para_data = table2array(readtable(['profile_lik_para_' data_src{data_src_indx} '_Exp_data_pd_mod_' num2str(pop_model(pop_model_indx)) '_b_t_ratio_' num2str(birth_tran_ratio) '_inverse_r1_' num2str(1/r1) '_inverse_r2_' num2str(inverse_r2) '_K_20_mill' '_para_steps_' num2str(profile_sample_size) '.csv']));
                  para_data = table2array(readtable(['profile_lik_para_' data_src{data_src_indx} '_Exp_data_pd_mod_' num2str(pop_model(pop_model_indx)) '_b_t_ratio_' num2str(birth_tran_ratio) '_inverse_r1_' num2str(1/r1) '_inverse_r2_' num2str(inverse_r2) '_para_steps_' num2str(profile_sample_size) '.csv']));
%                     para_data = table2array(readtable(['profile_lik_para_post_itr_samp_pd_mod_' num2str(pop_model(pop_model_indx)) '_b_t_ratio_' num2str(birth_tran_ratio) '_inverse_r1_' num2str(54) '_inverse_r2_' num2str(35) '_para_steps_' num2str(profile_sample_size) '.csv']));
                  
                    difference_indx = 1;

%                   if(pop_model(pop_model_indx) ~= 4)
% 
% %                         para_data = table2array(readtable(['profile_lik_para_' data_src{data_src_indx} '_Exp_data_pd_mod_' num2str(pop_model(pop_model_indx)) '_b_t_ratio_' num2str(birth_tran_ratio) '_inverse_r1_' num2str(1/r1) '_para_steps_' num2str(profile_sample_size) '.csv']));
% %                           para_data = readmatrix(['profile_lik_para_post_itr_samp_pd_mod_' num2str(pop_model(pop_model_indx)) '_b_t_ratio_' num2str(birth_tran_ratio) '_inverse_r1_' num2str(1/r1)  '_dynamic_para_steps_default_value_' num2str(profile_sample_size) '.csv']);
% 
%                   else
% %                         para_data = readmatrix(['profile_lik_para_post_itr_samp_pd_mod_' num2str(pop_model(pop_model_indx)) '_b_t_ratio_' num2str(birth_tran_ratio) '_inverse_r1_' num2str(1/r1)  '_dynamic_para_steps_default_value_' num2str(profile_sample_size) '_negative_comp.csv']);
% 
% %                         para_data = table2array(readtable(['profile_lik_para_' data_src{data_src_indx} '_Exp_data_pd_mod_' num2str(pop_model(pop_model_indx)) '_b_t_ratio_' num2str(birth_tran_ratio) '_inverse_r1_' num2str(1/r1) '_para_steps_' num2str(profile_sample_size) '_negative_comp.csv']));
%                     para_data = table2array(readtable(['profile_lik_para_post_itr_samp_pd_mod_' num2str(pop_model(pop_model_indx)) '_b_t_ratio_' num2str(birth_tran_ratio) '_inverse_r1_' num2str(1/54) '_inverse_r2_' num2str(35) '_para_steps_' num2str(profile_sample_size) '.csv']));
% 
%                   end
%                   difference_indx = 0;
                for para_indx = (1+difference_indx):(size(para_data,2)-2)
                    
                x_rate_all = para_data(para_data(:,end) == para_indx,1:end-2);
                x_rate = para_data(para_data(:,end) == para_indx,para_indx);

                y = para_data(para_data(:,end) == para_indx,end-1);                    

                    % based on the choosen alpha = 0.954, define the delta chi sqr avlue
                                        
                    if(isempty(r2))
                        switch pop_model(pop_model_indx)
                            case{6} % i.e. population model with 1 parameter
                                delta_chi_sq_all = 4;
                            case {1,2,10,11,12,20,32} % i.e. population model with 4 parameters
                                delta_chi_sq_all = 9.70;
                            case {3,7,18} % i.e. population model with 3 parameters
                                delta_chi_sq_all = 8.02;
                            case {4,5,13,14,15,16,17,19,29} % i.e. population model with 5 parameters
                                delta_chi_sq_all = 11.3;
                            case {8,9,22,23,24,33} % i.e. population model with 6 parameters
                                delta_chi_sq_all = 12.8;
                            case {21} % i.e. population model with 7 parameters
                                delta_chi_sq_all = 14.067;
                        end
                    else
                        switch pop_model(pop_model_indx)
                            case{6} % i.e. population model with 1 parameter
                                delta_chi_sq_all = 4;
                            case{3} % i.e. population model with 2 parameter
                                delta_chi_sq_all = 6.17;
                            case {4,5,13,14,15,16,17,19,29} % i.e. population model with 4 parameters
                                delta_chi_sq_all = 9.70;
                            case {1,2,10,11,12,20,32} % i.e. population model with 3 parameters
                                delta_chi_sq_all = 8.02;
                            case {8,9,22,23,24,33} % i.e. population model with 5 parameters
                                delta_chi_sq_all = 11.3;
                            case {21} % i.e. population model with 6 parameters
                                delta_chi_sq_all = 12.8;
                        end
                    end
                    
                    [minima,posn] = min(y); %finds the minimum in y and the corresponding index
                    
                    chi_sq_threshold_all = minima + delta_chi_sq_all;%the chi_sq threshold for all params
                    
                    %finding the confidence range
                    
                    left_indx = find((y(1:posn) - chi_sq_threshold_all > 0)); % here, left indx is for rate para
                    
                    right_indx = find((y(posn:end) - chi_sq_threshold_all > 0),1); % here, right indx is for rate para
                    
                    if(~isempty(left_indx))
                    actual_lower_bound(pop_model_indx,para_indx-difference_indx) = x_rate(left_indx(end));
                    else
                    actual_lower_bound(pop_model_indx,para_indx-difference_indx) = lower_bound(para_indx);
                    end

                    if(~isempty(right_indx))
                    actual_upper_bound(pop_model_indx,para_indx-difference_indx) = x_rate(posn + right_indx -1 );
                    else
                    actual_upper_bound(pop_model_indx,para_indx-difference_indx) = upper_bound(para_indx);
                    end

                    if(~(~isempty(left_indx) || ~isempty(right_indx)))

                            iden_mat(pop_model_indx,para_indx-difference_indx) = -1;

                    elseif(xor(~isempty(left_indx), ~isempty(right_indx)))

                            iden_mat(pop_model_indx,para_indx-difference_indx) = 0;
                                        
                    elseif(~isempty(left_indx) && ~isempty(right_indx))

                            iden_mat(pop_model_indx,para_indx-difference_indx) = 1;
                    end
            
                    goodness_mat(pop_model_indx,para_indx+1) = min(y);
                    
                end

                iden_mat(pop_model_indx,[para_indx+1-difference_indx:end]) = NaN;
actual_lower_bound(pop_model_indx,[para_indx+1-difference_indx:end]) = NaN;
actual_upper_bound(pop_model_indx,[para_indx+1-difference_indx:end]) = NaN;
            
        end
% to reaarnage the parameter sequence for the GCI&T model
        iden_mat(end,:) = iden_mat(end,[4 5 1 2 3 6]);
        actual_lower_bound(end,:) = actual_lower_bound(end,[4 5 1 2 3 6]);
        actual_upper_bound(end,:) = actual_upper_bound(end,[4 5 1 2 3 6]);

%         cd('C:\Users\Asus\OneDrive - Indian Institute of Science\GitHub\EMT-State-Transition-2.0\Experimental Data Results')
        
        figure(1)
%         yvalues = string(pop_model);
        xvalues = string(1:max_para);%while doing analysis for yamamoto models execpt 21
%         xvalues = {'r_e','t_{me}','t_{em}','\alpha','\beta'};
        yvalues = pop_model_names;
        h = heatmap(xvalues,yvalues,iden_mat);
        h.Title = 'Identifiability comaparison';
        h.XLabel = 'Parameters';
        h.YLabel = 'Pop model';
        h.FontSize = 14;
        h.CellLabelColor = 'none';
        h.Units = 'normalized';
        colormap(mymap);
        h.ColorLimits = [-1 1];
        colorbar off;
        
                figure(2)
%         yvalues = string(pop_model);
        xvalues = string(1:max_para);%while doing analysis for yamamoto models execpt 21
%         xvalues = {'r_e','t_{me}','t_{em}','\alpha','\beta'};
        yvalues = pop_model_names;
        h = heatmap(xvalues,yvalues,actual_lower_bound);
        h.Title = 'Lower bound';
        h.XLabel = 'Parameters';
        h.YLabel = 'Pop model';
        h.FontSize = 14;
%         h.CellLabelColor = 'none';
%         h.Units = 'normalized';
%         colormap(mymap);
%         h.ColorLimits = [-1 1];
        colormap white
        colorbar off;
        
        figure(3)
%         yvalues = string(pop_model);
        xvalues = string(1:max_para);%while doing analysis for yamamoto models execpt 21
%         xvalues = {'r_e','t_{me}','t_{em}','\alpha','\beta'};
        yvalues = pop_model_names;
        h = heatmap(xvalues,yvalues,actual_upper_bound);
        h.Title = 'Upper bound';
        h.XLabel = 'Parameters';
        h.YLabel = 'Pop model';
        h.FontSize = 14;
%         h.CellLabelColor = 'none';
%         h.Units = 'normalized';
%         colormap(mymap);
%         h.ColorLimits = [-1 1];
        colormap white
        colorbar off;
        

%         h.Position = [0.1300 0.1100 0.7179 0.179];
%         cd('C:\Users\Asus\OneDrive - Indian Institute of Science\GitHub\EMT-State-Transition-2.0\Manuscript figures\Figure 3 and related SI')
        if(difference_indx == 0)
            saveas(figure(1),['Identifiability ' data_src{data_src_indx} ' data b_t_ratio ' num2str(birth_tran_ratio) ' inverse r1 ' num2str(1/r1) '.png'])
                       saveas(figure(2),['Lower bound ' data_src{data_src_indx} ' data b_t_ratio ' num2str(birth_tran_ratio) ' inverse r1 ' num2str(1/r1)  '.png'])
            saveas(figure(3),['Upper bound ' data_src{data_src_indx} ' data b_t_ratio ' num2str(birth_tran_ratio) ' inverse r1 ' num2str(1/r1)  '.png'])

        else
            saveas(figure(1),['Identifiability ' data_src{data_src_indx} ' data b_t_ratio ' num2str(birth_tran_ratio) ' inverse r1 ' num2str(1/r1) '_inverse_r2_' num2str(inverse_r2) '.png'])
            saveas(figure(2),['Lower bound ' data_src{data_src_indx} ' data b_t_ratio ' num2str(birth_tran_ratio) ' inverse r1 ' num2str(1/r1) '_inverse_r2_' num2str(inverse_r2) '.png'])
            saveas(figure(3),['Upper bound ' data_src{data_src_indx} ' data b_t_ratio ' num2str(birth_tran_ratio) ' inverse r1 ' num2str(1/r1) '_inverse_r2_' num2str(inverse_r2) '.png'])


        end
        
%         figure
%         yvalues = string(pop_model);
% %         xvalues = string(1:8);
%          xvalues = string(1:7);%while doing analysis for yamamoto models execpt 21
% %        xvalues = {'r_m','r_e','t_{me}','t_{em}','\alpha','\beta'};
%         h = heatmap(xvalues,yvalues,goodness_mat,'CellLabelColor','none');
% %        h = heatmap(goodness_mat,'CellLabelColor','none')
%         
%         h.Title = 'Goodness comaparison';
%         h.XLabel = 'Parameters';
%         h.YLabel = 'Pop model';
%         colormap('parula')
% %         colorbar off
%         saveas(gcf,['Goodness ' data_src{data_src_indx} ' data b_t_ratio ' num2str(birth_tran_ratio) ' inverse r1 ' num2str(1/r1) '_inverse_r2_' num2str(inverse_r2) '_K_20_mill' '.png'])
    end
end
