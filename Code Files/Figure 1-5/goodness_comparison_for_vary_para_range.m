clear all
close all
search_method = 1;
% for Bhatia et al Data
% pop_model_indx_set = [3 11 12 14 10 13 4 4.1];
profile_sample_size = 1000;  % number of parameter sets present in each data file
r1 = 1/50;
% inverse_r2 = 35;
r2= []; %1/(r1*inverse_r2);
pop_model_indx_set = [3 4 11 12 14 10 13];
% Pop_model_names = ["G&T" , "GC&T"];
% 
Pop_model_names = ["G&T" , "GI&T", "G&T-Mr","G&T-Er", "G&T-EMr", "G&T-Mi", "G&T-Mi-Mr"];
b_t_ratio_set = [5 10 20 50 100];
% for Yamamoto et al. data
% pop_model_indx_set = [20 21 26 27 25 23 33]; %[20 22 23 24 25 19 21];
% pop_model_indx_set = 27%[20 21 26 27 25 23]; %[20 22 23 24 25 19 21];
% Pop_model_names = ["GC_{s}&T" , "GC_{a}&T", "GC_{s}&T-Er", "GC_{s}&T-Mr" , "GC_{s}&T-EMr",  "GC_{s}&T-Mi-Er", "GC_{s}I&T" ];

% pop_model_indx_set = [21]; %[20 22 23 24 25 19 21];
% Pop_model_names = ["GC&T"];
% b_t_ratio_set = [50 100 150 200 250 300];
% profile_sample_size = 2000;  % number of parameter sets present in each data file
% r1 = 1/54;
% inverse_r2 = 35;
% r2= 1/(r1*inverse_r2);
max_para =3; % compare Epi growth and transition rates which are common to all models; influence parameters have different roles depending on the model
cd('C:\Users\Asus\OneDrive - Indian Institute of Science\GitHub\EMT-State-Transition-2.0\Parameters sets from experimental data')

goodness_matrix = zeros(length(pop_model_indx_set),length(b_t_ratio_set));
aic= zeros(length(pop_model_indx_set),length(b_t_ratio_set));
bic = zeros(length(pop_model_indx_set),length(b_t_ratio_set));
aicc = zeros(length(pop_model_indx_set),length(b_t_ratio_set));


actual_b_t_ratio = zeros(profile_sample_size,2,length(pop_model_indx_set),length(b_t_ratio_set));

combined_data = zeros(length(pop_model_indx_set)*length(b_t_ratio_set)*max_para*profile_sample_size,5);


itr_num = 1; % to be incremented to analysis second version of model 4

i = 1; % variable to track pop_model_indx in the literation
for pop_model_indx = pop_model_indx_set
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
            case {8,9,22,23,24,25} % i.e. population model with 6 parameters
                delta_chi_sq_all = 12.8;
                num_para = 6;
            case {21} % i.e. population model with 7 parameters
                delta_chi_sq_all = 14.067;
                num_para = 7;
        end
    else

        switch pop_model_indx
            case{6} % i.e. population model with 1 parameter
                delta_chi_sq_all = 4;
            case{3} % i.e. population model with 2 parameter
                delta_chi_sq_all = 6.17;
            case {4,5,13,14,15,16,17,19,29,26,27} % i.e. population model with 4 parameters
                delta_chi_sq_all = 9.70;
                num_para = 5;
            case {1,2,10,11,12,20} % i.e. population model with 3 parameters
                delta_chi_sq_all = 8.02;
                num_para = 4;
            case {8,9,22,23,24,25} % i.e. population model with 5 parameters
                delta_chi_sq_all = 11.3;
                num_para = 6;
            case {21} % i.e. population model with 6 parameters
                delta_chi_sq_all = 12.8;
                num_para = 7;
        end

    end
    j = 1; % variable to track b_t_ratio in the literation
    for b_t_ratio = b_t_ratio_set
if(isempty(r2))
        if(pop_model_indx <15 )
            %             parameters = table2array(readtable(['para_Bhatia_data_pd_mod_' num2str(pop_model_indx) '_s_m_' num2str(search_method) '_b_t_ratio_' num2str(b_t_ratio) '_inverse_r1_' num2str(1/r1) '.csv']));
            if( pop_model_indx ~=4.1)
                profile_para = table2array(readtable(['profile_lik_para_Bhatia_Exp_data_pd_mod_' num2str(pop_model_indx) '_b_t_ratio_' num2str(b_t_ratio) '_inverse_r1_' num2str(1/r1)  '_para_steps_' num2str(profile_sample_size) '.csv']));
            else
                profile_para = table2array(readtable(['profile_lik_para_Bhatia_Exp_data_pd_mod_' num2str(4) '_b_t_ratio_' num2str(b_t_ratio) '_inverse_r1_' num2str(1/r1)  '_para_steps_' num2str(profile_sample_size) '_negative_comp.csv']));
            end
            num_rep = 3;
            num_init_cond = 2;
            num_time_pts = 4;
        else
            %             parameters = table2array(readtable(['para_Bhatia_data_pd_mod_' num2str(pop_model_indx) '_s_m_' num2str(search_method) '_b_t_ratio_' num2str(b_t_ratio) '_inverse_r1_' num2str(1/r1) '.csv']));
            profile_para = table2array(readtable(['profile_lik_para_Yamamoto_Exp_data_pd_mod_' num2str(pop_model_indx) '_b_t_ratio_' num2str(b_t_ratio) '_inverse_r1_' num2str(1/r1)  '_para_steps_' num2str(profile_sample_size) '.csv']));
            num_rep = 3;
            num_init_cond = 5;
            num_time_pts = 5;
        end
else
        if(pop_model_indx <15 )
            %             parameters = table2array(readtable(['para_Bhatia_data_pd_mod_' num2str(pop_model_indx) '_s_m_' num2str(search_method) '_b_t_ratio_' num2str(b_t_ratio) '_inverse_r1_' num2str(1/r1) '.csv']));
            if( pop_model_indx ~=4.1)
                profile_para = table2array(readtable(['profile_lik_para_Bhatia_Exp_data_pd_mod_' num2str(pop_model_indx) '_b_t_ratio_' num2str(b_t_ratio) '_inverse_r1_' num2str(1/r1)  '_para_steps_' num2str(profile_sample_size) '.csv']));
            else
                profile_para = table2array(readtable(['profile_lik_para_Bhatia_Exp_data_pd_mod_' num2str(4) '_b_t_ratio_' num2str(b_t_ratio) '_inverse_r1_' num2str(1/r1)  '_para_steps_' num2str(profile_sample_size) '_negative_comp.csv']));
            end
            num_rep = 3;
            num_init_cond = 2;
            num_time_pts = 4;
        else
            %             parameters = table2array(readtable(['para_Bhatia_data_pd_mod_' num2str(pop_model_indx) '_s_m_' num2str(search_method) '_b_t_ratio_' num2str(b_t_ratio) '_inverse_r1_' num2str(1/r1) '.csv']));
            profile_para = table2array(readtable(['profile_lik_para_Yamamoto_Exp_data_pd_mod_' num2str(pop_model_indx) '_b_t_ratio_' num2str(b_t_ratio) '_inverse_r1_' num2str(1/r1) '_inverse_r2_' num2str(inverse_r2)  '_para_steps_' num2str(profile_sample_size) '.csv']));
            num_rep = 3;
            num_init_cond = 5;
            num_time_pts = 5;
        end
end


%         if(pop_model_indx ~= 4 && pop_model_indx ~= 4.1)
%             profile_para(:,end) = profile_para(:,end) .* (num_rep*num_init_cond*num_time_pts);
%         end
        %             goodness_matrix(i,j) = parameters(1,end);

        [goodness_matrix(i,j),min_error_indx] = min(profile_para(:,end-1));
        [aic(i,j),bic(i,j)] = aicbic(-1/2*goodness_matrix(i,j),num_para,(num_time_pts*num_rep*num_init_cond));
        aicc( i, j) = aic(i,j) + (2*num_para*(num_para + 1))/((num_time_pts*num_rep*num_init_cond) - num_para - 1);

        %         goodness_matrix(i,j) = parameters(min_error_indx,end);
        if(pop_model_indx~= 33)
        for para_indx = 2:size(profile_para,2)-2 % for parameters 2 (t1) and 3(t3) in the nondimensionalised models
                        
            combined_data((i-1)*profile_sample_size*length(b_t_ratio_set)*max_para + (j-1)*max_para*profile_sample_size + (para_indx-1)*profile_sample_size + 1 : (i-1)*profile_sample_size*length(b_t_ratio_set)*max_para + (j-1)*max_para*profile_sample_size + (para_indx-1)*profile_sample_size + profile_sample_size,1) = (pop_model_indx)*ones(profile_sample_size,1);
            combined_data((i-1)*profile_sample_size*length(b_t_ratio_set)*max_para + (j-1)*max_para*profile_sample_size + (para_indx-1)*profile_sample_size + 1 : (i-1)*profile_sample_size*length(b_t_ratio_set)*max_para + (j-1)*max_para*profile_sample_size + (para_indx-1)*profile_sample_size + profile_sample_size,2) = b_t_ratio*ones(profile_sample_size,1);
            combined_data((i-1)*profile_sample_size*length(b_t_ratio_set)*max_para + (j-1)*max_para*profile_sample_size + (para_indx-1)*profile_sample_size + 1 : (i-1)*profile_sample_size*length(b_t_ratio_set)*max_para + (j-1)*max_para*profile_sample_size + (para_indx-1)*profile_sample_size + profile_sample_size,3) = para_indx*ones(profile_sample_size,1);
%             combined_data((i-1)*profile_sample_size*length(b_t_ratio_set)*max_para + (j-1)*max_para*profile_sample_size + (para_indx-1)*profile_sample_size + 1 : (i-1)*profile_sample_size*length(b_t_ratio_set)*max_para + (j-1)*max_para*profile_sample_size + (para_indx-1)*profile_sample_size + profile_sample_size,4) = profile_para((para_indx-1)*profile_sample_size+1:para_indx*profile_sample_size,para_indx);
            combined_data((i-1)*profile_sample_size*length(b_t_ratio_set)*max_para + (j-1)*max_para*profile_sample_size + (para_indx-1)*profile_sample_size + 1 : (i-1)*profile_sample_size*length(b_t_ratio_set)*max_para + (j-1)*max_para*profile_sample_size + (para_indx-1)*profile_sample_size + profile_sample_size,4) = profile_para(profile_para(:,end) == para_indx,para_indx);

            combined_data((i-1)*profile_sample_size*length(b_t_ratio_set)*max_para + (j-1)*max_para*profile_sample_size + (para_indx-1)*profile_sample_size + 1 : (i-1)*profile_sample_size*length(b_t_ratio_set)*max_para + (j-1)*max_para*profile_sample_size + (para_indx-1)*profile_sample_size + profile_sample_size,5) = (delta_chi_sq_all*ones(profile_sample_size,1)+min(profile_para(profile_para(:,end) == para_indx,end-1))) - profile_para(profile_para(:,end) == para_indx,end-1);
        end
        else
        for para_indx = 5:size(profile_para,2)-2 % for parameters 2 (t1) and 3(t2) in the nondimensionalised models
                        
            combined_data((i-1)*profile_sample_size*length(b_t_ratio_set)*max_para + (j-1)*max_para*profile_sample_size + (para_indx-1)*profile_sample_size + 1 : (i-1)*profile_sample_size*length(b_t_ratio_set)*max_para + (j-1)*max_para*profile_sample_size + (para_indx-1)*profile_sample_size + profile_sample_size,1) = (pop_model_indx)*ones(profile_sample_size,1);
            combined_data((i-1)*profile_sample_size*length(b_t_ratio_set)*max_para + (j-1)*max_para*profile_sample_size + (para_indx-1)*profile_sample_size + 1 : (i-1)*profile_sample_size*length(b_t_ratio_set)*max_para + (j-1)*max_para*profile_sample_size + (para_indx-1)*profile_sample_size + profile_sample_size,2) = b_t_ratio*ones(profile_sample_size,1);
            combined_data((i-1)*profile_sample_size*length(b_t_ratio_set)*max_para + (j-1)*max_para*profile_sample_size + (para_indx-1)*profile_sample_size + 1 : (i-1)*profile_sample_size*length(b_t_ratio_set)*max_para + (j-1)*max_para*profile_sample_size + (para_indx-1)*profile_sample_size + profile_sample_size,3) = (para_indx-3)*ones(profile_sample_size,1);
%             combined_data((i-1)*profile_sample_size*length(b_t_ratio_set)*max_para + (j-1)*max_para*profile_sample_size + (para_indx-1)*profile_sample_size + 1 : (i-1)*profile_sample_size*length(b_t_ratio_set)*max_para + (j-1)*max_para*profile_sample_size + (para_indx-1)*profile_sample_size + profile_sample_size,4) = profile_para((para_indx-1)*profile_sample_size+1:para_indx*profile_sample_size,para_indx);
            combined_data((i-1)*profile_sample_size*length(b_t_ratio_set)*max_para + (j-1)*max_para*profile_sample_size + (para_indx-1)*profile_sample_size + 1 : (i-1)*profile_sample_size*length(b_t_ratio_set)*max_para + (j-1)*max_para*profile_sample_size + (para_indx-1)*profile_sample_size + profile_sample_size,4) = profile_para(profile_para(:,end) == para_indx,para_indx);

            combined_data((i-1)*profile_sample_size*length(b_t_ratio_set)*max_para + (j-1)*max_para*profile_sample_size + (para_indx-1)*profile_sample_size + 1 : (i-1)*profile_sample_size*length(b_t_ratio_set)*max_para + (j-1)*max_para*profile_sample_size + (para_indx-1)*profile_sample_size + profile_sample_size,5) = (delta_chi_sq_all*ones(profile_sample_size,1)+min(profile_para(profile_para(:,end) == para_indx,end-1))) - profile_para(profile_para(:,end) == para_indx,end-1);
        end
        end
        if(false)
%         if(pop_model_indx == 4.1 || pop_model_indx == 21)
            % scatter plot 
            figure('Units','normalized','Position',[0.3177 0.4722 0.2901 0.3750])
            ax =gca;
            scatter_data = [];
            for para_indx = [3 4 5]

                confidence_sign = (delta_chi_sq_all*ones(profile_sample_size,1)+min(profile_para(profile_para(:,end) == para_indx,end-1))) - profile_para(profile_para(:,end) == para_indx,end-1);
                filter_profile_data = profile_para(profile_para(:,end) == para_indx,:);            
                filter_profile_data = filter_profile_data(confidence_sign>=0,:);
                scatter_data = [scatter_data; filter_profile_data];
            end
            [~,order_indices] = sort(scatter_data(:,3),1);
            scatter_data = scatter_data(order_indices,:);
            s = scatter(scatter_data(:,4),scatter_data(:,5),[],scatter_data(:,3),"filled");
            cb = colorbar;
%             cb.Label = 'normalized E growth rate';
            axis square
            ylabel('Influence parameter \beta')
            xlabel('Influence parameter \alpha')
%             AlphaData = ((scatter_data(:,1)-min(scatter_data(:,1)))/max(scatter_data(:,1)));
%             AlphaData(AlphaData < 0.5,:) = 0.7;
%             alpha(AlphaData);
%             s.AlphaData = AlphaData;
%             s.MarkerFaceAlpha = 'flat';
            grid on
            ax.FontSize = 15;
            % axis square
            cwd = cd;
            cd('C:\Users\Asus\OneDrive - Indian Institute of Science\GitHub\EMT-State-Transition-2.0\Manuscript figures\Figure 2 and related SI')
            if(pop_model_indx == 4.1)
                saveas(gcf,['Bhatia_data_scatter_plot_influence_parameter_EM_transition_color_label_pop_model_' num2str(4)  '_b_t_ratio_ ' num2str(b_t_ratio) '_inverse_r1_' num2str(1/r1) '.png'])
            else
                saveas(gcf,['Yamamoto_data_scatter_plot_influence_parameter_EM_transition_color_label_pop_model_' num2str(pop_model_indx) '_b_t_ratio_ ' num2str(b_t_ratio) '_inverse_r1_' num2str(1/r1) '.png'])
            end

            cd(cwd);
        end


        j = j+1;

    end
    i = i+1;
end

% % to generate data frame
combined_data_table = splitvars(table(combined_data));
combined_data_table.Properties.VariableNames = {'Pop_model' 'b_t_ratio' 'para_indx','para_value','confidence_sign'}; % if confidence sign is positive then parameter is in 95 percent confidence interval otherwisw not
% cd('C:\Users\Asus\OneDrive - Indian Institute of Science\GitHub\EMT-State-Transition-2.0\Parameters sets from experimental data')
% writetable(combined_data_table,['Yamamoto_data_actual_b_t_ratio_combined_profile_like_data_r1_' num2str(1/r1) '.xlsx']);

figure_position = [0.6036 0.4028 0.3385 0.3259];
close all
% 
% cd('C:\Users\Asus\OneDrive - Indian Institute of Science\GitHub\EMT-State-Transition-2.0\Manuscript figures\Figure 2 and related SI')

% pop_model_names_seq = linspace(1.5,1.5*length(pop_model_indx_set),length(pop_model_indx_set));

% % to plot M grow vs M E transition rates box plots
% figure('Units','normalized','Position',figure_position)
% 
% ax =gca;
% filter_table = combined_data_table(combined_data_table.para_indx == 1,:);
% filter_table = filter_table(filter_table.confidence_sign>=0,:);
% % pop_model = categorical(filter_table.Pop_model,pop_model_indx_set,string(pop_model_indx_set));
% pop_model = zeros(height(filter_table),1);
% for j = 1:length(pop_model_indx_set)
%     pop_model(filter_table.Pop_model == pop_model_indx_set(j)) = pop_model_names_seq(j);
% end
% boxchart(pop_model,(filter_table.para_value),'GroupByColor',filter_table.b_t_ratio);
% xlim([-1+min(pop_model_names_seq) 1+max(pop_model_names_seq)])
% ax.XTick = pop_model_names_seq;
% ax.XTickLabels = Pop_model_names;
% lg = legend(string(b_t_ratio_set),'Location','northeastoutside');
% title(lg,'g t ratio')
% 
% ylabel('Normalised E growth rate')
% xlabel('Population model')
% grid on
% % axis square
% 
% ax.FontSize = 14;
% saveas(gcf,['Bhatia_data_E_grow_rate_norm_vs_pop_model_various_b_t_ratio_s_m_' num2str(search_method)  '_inverse_r1_' num2str(1/r1) '.png'])


% to plot M grow vs M E transition rates box plots
figure('Units','normalized','Position',figure_position)
ax =gca;
filter_table = combined_data_table(combined_data_table.b_t_ratio == 250,:);
filter_table = filter_table(filter_table.confidence_sign>=0,:);
filter_table = filter_table(logical((filter_table.para_indx == 2) + (filter_table.para_indx == 3) ),:);

% pop_model = zeros(height(filter_table),1);
% for j = 1:length(pop_model_indx_set)
%     pop_model(filter_table.Pop_model == pop_model_indx_set(j)) = pop_model_names_seq(j);
% end
boxchart(categorical(filter_table.para_indx),(filter_table.para_value),'GroupByColor',filter_table.Pop_model);
% xlim([-1+min(pop_model_names_seq) 1+max(pop_model_names_seq)])
% ax.XTick = pop_model_names_seq;
ax.XTickLabels = ["t_{me}" , "t_{em}"];
lg = legend(Pop_model_names,'Location','northeastoutside');
title(lg,'Pop model')
grid on
ylabel('Values')
xlabel('Transition rates')
grid on
% axis square
ax.FontSize = 14;
saveas(gcf,['Yamamoto_data_M_E_trans_rate_norm_vs_pop_model_various_b_t_ratio_s_m_' num2str(search_method)  '_inverse_r1_' num2str(1/r1) '_inverse_r2_' num2str(inverse_r2) '.png'])




% to plot M grow vs M E transition rates box plots
figure('Units','normalized','Position',figure_position)
ax =gca;
filter_table = combined_data_table(combined_data_table.para_indx == 2,:);
filter_table = filter_table(filter_table.confidence_sign>=0,:);
pop_model = zeros(height(filter_table),1);
for j = 1:length(pop_model_indx_set)
    pop_model(filter_table.Pop_model == pop_model_indx_set(j)) = pop_model_names_seq(j);
end
boxchart(pop_model,(filter_table.para_value),'GroupByColor',filter_table.b_t_ratio);
xlim([-1+min(pop_model_names_seq) 1+max(pop_model_names_seq)])
ax.XTick = pop_model_names_seq;
ax.XTickLabels = Pop_model_names;
lg = legend(string(b_t_ratio_set),'Location','northeastoutside');
title(lg,'g t ratio')

ylabel('Normalised M E trans rate')
xlabel('Population model')
grid on
% axis square
ax.FontSize = 14;
saveas(gcf,['Yamamoto_data_M_E_trans_rate_norm_vs_pop_model_various_b_t_ratio_s_m_' num2str(search_method)  '_inverse_r1_' num2str(1/r1) '_inverse_r2_' num2str(inverse_r2) '.png'])

% to plot M grow vs E M transition rates box plots
figure('Units','normalized','Position',figure_position)
ax = gca;
filter_table = combined_data_table(combined_data_table.para_indx == 3,:);
filter_table = filter_table(filter_table.confidence_sign>=0,:);
pop_model = zeros(height(filter_table),1);
for j = 1:length(pop_model_indx_set)
    pop_model(filter_table.Pop_model == pop_model_indx_set(j)) = pop_model_names_seq(j);
end
boxchart(pop_model,(filter_table.para_value),'GroupByColor',filter_table.b_t_ratio);
xlim([-1+min(pop_model_names_seq) 1+max(pop_model_names_seq)])
ax.XTick = pop_model_names_seq;
ax.XTickLabels = Pop_model_names;
lg = legend(string(b_t_ratio_set),'Location','northeastoutside');
title(lg,'g t ratio')

xlabel('Population model')
% ylabel('log_{10}(Minimum Chi-sqr)')
ylabel('Normalised E M trans rate')
grid on
% axis square

ax.FontSize = 14;
saveas(gcf,['Yamamoto_data_E_M_trans_rate_norm_vs_pop_model_various_b_t_ratio_s_m_' num2str(search_method)  '_inverse_r1_' num2str(1/r1) '_inverse_r2_' num2str(inverse_r2) '.png'])

figure('Units','normalized','Position',figure_position)
ax =gca;
filter_table = combined_data_table(combined_data_table.para_indx == 2,:);
filter_table = filter_table(filter_table.confidence_sign>=0,:);
pop_model = zeros(height(filter_table),1);
for j = 1:length(pop_model_indx_set)
    pop_model(filter_table.Pop_model == pop_model_indx_set(j)) = pop_model_names_seq(j);
end
boxchart(pop_model,1./(filter_table.para_value),'GroupByColor',filter_table.b_t_ratio);
xlim([-1+min(pop_model_names_seq) 1+max(pop_model_names_seq)])
ax.XTick = pop_model_names_seq;
ax.XTickLabels = Pop_model_names;
lg = legend(string(b_t_ratio_set),'Location','northeastoutside');
title(lg,'g t ratio')

ylabel('M grow to M E trans rate ratio')
xlabel('Population model')
grid on
ax.FontSize = 14;
% axis square

saveas(gcf,['Yamamoto_data_actual_M_grow_to_M_E_trans_rate_ratio_vs_pop_model_various_b_t_ratio_s_m_' num2str(search_method)  '_inverse_r1_' num2str(1/r1) '_inverse_r2_' num2str(inverse_r2) '.png'])

% to plot M grow vs E M transition rates box plots
figure('Units','normalized','Position',figure_position)
ax = gca;
filter_table = combined_data_table(combined_data_table.para_indx == 3,:);
filter_table = filter_table(filter_table.confidence_sign>=0,:);
pop_model = zeros(height(filter_table),1);
for j = 1:length(pop_model_indx_set)
    pop_model(filter_table.Pop_model == pop_model_indx_set(j)) = pop_model_names_seq(j);
end
boxchart(pop_model,1./(filter_table.para_value),'GroupByColor',filter_table.b_t_ratio);
xlim([-1+min(pop_model_names_seq) 1+max(pop_model_names_seq)])
ax.XTick = pop_model_names_seq;
ax.XTickLabels = Pop_model_names;
lg = legend(string(b_t_ratio_set),'Location','northeastoutside');
title(lg,'g t ratio')
xlabel('Population model')
% ylabel('log_{10}(Minimum Chi-sqr)')
ylabel('M grow to E M trans rate ratio')
grid on
ax.FontSize = 14;
% axis square

saveas(gcf,['Yamamoto_data_actual_M_grow_to_E_M_trans_rate_ratio_vs_pop_model_various_b_t_ratio_s_m_' num2str(search_method)  '_inverse_r1_' num2str(1/r1) '_inverse_r2_' num2str(inverse_r2) '.png'])
% 

% % % transition rate influence parameter 
% figure('Units','normalized','Position',[0.3177 0.4722 0.2901 0.3750])
% ax =gca;
% filter_table = combined_data_table(combined_data_table.para_indx == 4,:);
% filter_table = filter_table(filter_table.confidence_sign>=0,:);
% pop_model = zeros(height(filter_table),1);
% for j = 1:length(pop_model_indx_set)
%     pop_model(filter_table.Pop_model == pop_model_indx_set(j)) = pop_model_names_seq(j);
% end
% boxchart(pop_model,(filter_table.para_value),'GroupByColor',filter_table.b_t_ratio);
% xlim([-1+min(pop_model_names_seq) 1+max(pop_model_names_seq)])
% ax.XTick = pop_model_names_seq;
% ax.XTickLabels = Pop_model_names;
% lg = legend(string(b_t_ratio_set),'Location','northeastoutside');
% title(lg,'g t ratio')
% 
% ylabel('Influence parameter \alpha')
% xlabel('Population model')
% grid on
% ax.FontSize = 15;
% % axis square
% ylim([-1 1])
% ax.YTick = -1:0.5:1;
% saveas(gcf,['Yamamoto_data_influence_parameter_alpha_vs_pop_model_various_b_t_ratio_s_m_' num2str(search_method)  '_inverse_r1_' num2str(1/r1) '.png'])
% 
% figure('Units','normalized','Position',[0.3177 0.4722 0.2901 0.3750])
% ax =gca;
% filter_table = combined_data_table(combined_data_table.para_indx == 5,:);
% filter_table = filter_table(filter_table.confidence_sign>=0,:);
% pop_model = zeros(height(filter_table),1);
% for j = 1:length(pop_model_indx_set)
%     pop_model(filter_table.Pop_model == pop_model_indx_set(j)) = pop_model_names_seq(j);
% end
% boxchart(pop_model,(filter_table.para_value),'GroupByColor',filter_table.b_t_ratio);
% xlim([-1+min(pop_model_names_seq) 1+max(pop_model_names_seq)])
% ax.XTick = pop_model_names_seq;
% ax.XTickLabels = Pop_model_names;
% lg = legend(string(b_t_ratio_set),'Location','northeastoutside');
% title(lg,'g t ratio')
% ylim([-1 1])
% ax.YTick = -1:0.5:1;
% ylabel('Influence parameter \beta')
% xlabel('Population model')
% grid on
% ax.FontSize = 15;
% % axis square
% 
% saveas(gcf,['Yamamoto_data_influence_parameter_beta_vs_pop_model_various_b_t_ratio_s_m_' num2str(search_method)  '_inverse_r1_' num2str(1/r1) '.png'])


% % % transition rate influence parameter 
% figure('Units','normalized','Position',[0.3177 0.4722 0.2901 0.3750])
% ax =gca;
% filter_table = combined_data_table(combined_data_table.para_indx == 6,:);
% filter_table = filter_table(filter_table.confidence_sign>=0,:);
% pop_model = zeros(height(filter_table),1);
% for j = 1:length(pop_model_indx_set)
%     pop_model(filter_table.Pop_model == pop_model_indx_set(j)) = pop_model_names_seq(j);
% end
% boxchart(pop_model,(filter_table.para_value),'GroupByColor',filter_table.b_t_ratio);
% xlim([-1+min(pop_model_names_seq) 1+max(pop_model_names_seq)])
% ax.XTick = pop_model_names_seq;
% ax.XTickLabels = Pop_model_names;
% lg = legend(string(b_t_ratio_set),'Location','northeastoutside');
% title(lg,'g t ratio')
% 
% ylabel('K_m')
% xlabel('Population model')
% grid on
% ax.FontSize = 15;
% % axis square
% ylim([0 2*10^7])
% % ax.YTick = -1:0.5:1;
% saveas(gcf,['Yamamoto_data_M_carrying_capacity_vs_pop_model_various_b_t_ratio_s_m_' num2str(search_method)  '_inverse_r1_' num2str(1/r1) '.png'])
% 
% figure('Units','normalized','Position',[0.3177 0.4722 0.2901 0.3750])
% ax =gca;
% filter_table = combined_data_table(combined_data_table.para_indx == 7,:);
% filter_table = filter_table(filter_table.confidence_sign>=0,:);
% pop_model = zeros(height(filter_table),1);
% for j = 1:length(pop_model_indx_set)
%     pop_model(filter_table.Pop_model == pop_model_indx_set(j)) = pop_model_names_seq(j);
% end
% boxchart(pop_model,(filter_table.para_value),'GroupByColor',filter_table.b_t_ratio);
% xlim([-1+min(pop_model_names_seq) 1+max(pop_model_names_seq)])
% ax.XTick = pop_model_names_seq;
% ax.XTickLabels = Pop_model_names;
% lg = legend(string(b_t_ratio_set),'Location','northeastoutside');
% title(lg,'g t ratio')
% ylim([0 2*10^7])
% % ax.YTick = -1:0.5:1;
% ylabel('K_e')
% xlabel('Population model')
% grid on
% ax.FontSize = 15;
% % axis square
% 
% saveas(gcf,['Yamamoto_data_E_carrying_capacity_vs_pop_model_various_b_t_ratio_s_m_' num2str(search_method)  '_inverse_r1_' num2str(1/r1) '.png'])



% to plot the minimum chi-sqr value against pop model for various bt ratios
figure('Units','normalized','Position',figure_position)
X = categorical(Pop_model_names);
X = reordercats(X,Pop_model_names);
% heatmap(X,b_t_ratio_set,goodness_matrix',"ColorScaling","log",Title="Chi-square");
bar(X,log10(goodness_matrix))
% bar(X,(goodness_matrix))
ax = gca;
% grid on
% colormap parula
xlabel('Population model')
ylabel('log_{10}(Minimum Chi-sqr)')
% ylabel('g t ratio')
% ylabel('Minimum Chi-sqr')

lg = legend(string(b_t_ratio_set),'Location','northeastoutside');
title(lg,'g t ratio')
% axis square
ax.FontSize = 14;
colorbar off
% saveas(gcf,['Yamamoto_data_Min_chi_sqr_vs_pop_model_various_b_t_ratio_s_m_' num2str(search_method)  '_inverse_r1_' num2str(1/r1)  '_heartmap.png'])
saveas(gcf,['Yamamoto_data_Min_chi_sqr_vs_pop_model_various_b_t_ratio_s_m_' num2str(search_method)  '_inverse_r1_' num2str(1/r1) '_inverse_r2_' num2str(inverse_r2) '.png'])


% to plot the minimum chi-sqr value against pop model for various bt ratios
figure('Units','normalized','Position',[0.3177 0.5296 0.3141 0.3176])
X = categorical(Pop_model_names);
X = reordercats(X,Pop_model_names);
heatmap(X,b_t_ratio_set,aicc',"ColorScaling","log",Title="AICc");
% bar(X,log10(aicc))
ax = gca;
% grid on
colormap parula
xlabel('Population model')
% ylabel('log_{10}(AICc)')
ylabel('g t ratio')
% ylabel('AICc')

% lg = legend(string(b_t_ratio_set),'Location','northeastoutside');
% title(lg,'g t ratio')
% axis square
ax.FontSize = 14;
colorbar off

saveas(gcf,['Yamamoto_data_AICc_vs_pop_model_various_b_t_ratio_s_m_' num2str(search_method)  '_inverse_r1_' num2str(1/r1) '_heartmap.png'])

% to plot the minimum chi-sqr value against pop model for various bt ratios
figure('Units','normalized','Position',figure_position)
X = categorical(Pop_model_names);
X = reordercats(X,Pop_model_names);
bar(X,log10(aic))
% bar(X,(goodness_matrix))
ax = gca;
grid on
xlabel('Population model')
ylabel('log_{10}(AIC)')
% ylabel('Minimum Chi-sqr')

lg = legend(string(b_t_ratio_set),'Location','northeastoutside');
title(lg,'g t ratio')
% axis square
ax.FontSize = 15;
saveas(gcf,['Yamamoto_data_AIC_vs_pop_model_various_b_t_ratio_s_m_' num2str(search_method)  '_inverse_r1_' num2str(1/r1) '_heartmap.png'])


% to plot the minimum chi-sqr value against pop model for various bt ratios
% figure('Units','normalized','Position',figure_position)
X = categorical(Pop_model_names);
X = reordercats(X,Pop_model_names);
bar(X,log10(bic))
% bar(X,(goodness_matrix))
ax = gca;
grid on
xlabel('Population model')
ylabel('log_{10}(BIC)')
% ylabel('Minimum Chi-sqr')

lg = legend(string(b_t_ratio_set),'Location','northeastoutside');
title(lg,'g t ratio')
% axis square
ax.FontSize = 15;
saveas(gcf,['Yamamoto_data_BIC_vs_pop_model_various_b_t_ratio_s_m_' num2str(search_method)  '_inverse_r1_' num2str(1/r1) '.png'])


% % bar(X,1./actual_b_t_ratio_M_grow_to_M_E_trans)
% ax = gca;
% grid on
% xlabel('Population model')
% ylabel('actual g t ratio M grow to M E trans')
% legend(string(b_t_ratio_set),'Location','northeastoutside')
% ax.FontSize = 14;
% saveas(gcf,['Bhatia_data_actual_b_t_ratio_M_grow_to_M_E_trans_vs_pop_model_various_b_t_ratio_s_m_' num2str(search_method)  '_inverse_r1_' num2str(1/r1) '.png'])
% 
% bar(X,1./actual_b_t_ratio_M_grow_to_E_M_trans)
% ax = gca;
% grid on
% xlabel('Population model')
% ylabel('actual g t ratio M grow to E M trans')
% legend(string(b_t_ratio_set),'Location','northeastoutside')
% ax.FontSize = 14;
% saveas(gcf,['Bhatia_data_actual_b_t_ratio_M_grow_to_E_M_trans_vs_pop_model_various_b_t_ratio_s_m_' num2str(search_method)  '_inverse_r1_' num2str(1/r1) '.png'])

