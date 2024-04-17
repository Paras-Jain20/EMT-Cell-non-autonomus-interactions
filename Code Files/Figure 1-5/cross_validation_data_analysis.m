% cross validation data analysis %
close all
search_method = 1;

total_itr= 10;
pop_model = [20 21 26 27 25 23 33];
% pop_model_names = ["G&T" , "GC&T", "G&T-Er" , "G&T-Mr", "G&T-EMr", "G&T-Mi-Er"];
pop_model_names = ["GC_{s}&T" , "GC_{a}&T", "GC_{s}&T-Mr", "GC_{s}&T-Er" , "GC_{s}&T-EMr", "GC_{s}&T-Mi-Er","GC_{s}I&T" ];
birth_tran_ratio = 250;
init_cond = [1 2 3 4 5];
r1 = 1/54;
inverse_r2 = 35;
r2 = 1/(r1*inverse_r2);


% pop_model= [3 4 11 12 14 10 13];
% pop_model_names = ["G&T" , "GI&T", "G&T-Mr", "G&T-Er" , "G&T-EMr", "G&T-Mi", "G&T-Mi-Mr"];
% birth_tran_ratio = 50;
% init_cond = [1 2];
% r1 = 1/50;
% r2 = [];


%%%% control cross 1 validation analysis %%%%%%
total_itr= 10;

cd('C:\Users\Asus\OneDrive - Indian Institute of Science\GitHub\EMT-State-Transition-2.0\Parameter sets and generated pseudo data')
if(~isempty(find(pop_model > 16,1)))
    temp_cost = table2array(readtable(['Control_cross_1_validation_chi_sqaure_values_Yamamoto_Exp_data_b_t_ratio_' num2str(birth_tran_ratio) '_inverse_r1_' num2str(1/r1) '_inverse_r2_' num2str(inverse_r2)  '.csv']));
elseif (~isempty(find(pop_model < 16,1)))
    temp_cost = table2array(readtable(['Control_cross_1_validation_chi_sqaure_values_Bhatia_Exp_data_b_t_ratio_' num2str(birth_tran_ratio) '_inverse_r1_' num2str(1/r1)  '.csv']));
end
control_cost = zeros(length(pop_model),length(init_cond),total_itr);
% rearranging cost matrix

for pop_model_seq = 1:length(pop_model)
    for itr_indx = 1:total_itr
        cost_indx = find(sum(temp_cost(:,end-1:end) == [pop_model(pop_model_seq) itr_indx],2) == 2,1);
        control_cost(pop_model_seq,:,itr_indx) = temp_cost(cost_indx,1:length(init_cond));
    end
end

min_control_cost = min(control_cost,[],3);

total_itr= 10;

%%% cross 1 validation data analysis %%%

cd('C:\Users\Asus\OneDrive - Indian Institute of Science\GitHub\EMT-State-Transition-2.0\Parameter sets and generated pseudo data')
if(~isempty(find(pop_model > 16,1)))
    temp_cost = table2array(readtable(['cross_1_validation_chi_sqaure_values_Yamamoto_Exp_data_b_t_ratio_' num2str(birth_tran_ratio) '_inverse_r1_' num2str(1/r1) '_inverse_r2_' num2str(inverse_r2) '.csv']));
elseif (~isempty(find(pop_model < 16,1)))
    temp_cost = table2array(readtable(['cross_1_validation_chi_sqaure_values_Bhatia_Exp_data_b_t_ratio_' num2str(birth_tran_ratio) '_inverse_r1_' num2str(1/r1)  '.csv']));
end
cost = zeros(length(pop_model),length(init_cond),total_itr);
% rearranging cost matrix
for pop_model_seq = 1:length(pop_model)
    for itr_indx = 1:total_itr
        cost_indx = find(sum(temp_cost(:,end-1:end) == [pop_model(pop_model_seq) itr_indx],2) == 2,1);
        cost(pop_model_seq,:,itr_indx) = temp_cost(cost_indx,1:length(init_cond));
    end
end

min_cost = min(cost,[],3);
norm_min_cost = min_cost./min_control_cost;
% norm_min_cost = min_cost;

figure
h = heatmap(init_cond, categorical(pop_model_names),norm_min_cost,'CellLabelColor','none');
colormap parula
xlabel('Held out initial conditions');
ylabel('Population model');
ax = gca;
ax.FontSize = 14;
c = colorbar;
c.Label.String = 'log_{10}(chi-square)';

min_cost = min(cost,[],3);
% norm_min_cost = min_cost./min_control_cost;
norm_min_cost = min_control_cost;

figure
h = heatmap(init_cond, categorical(pop_model_names),(min_control_cost));
colormap white
xlabel('Held out initial conditions');
ylabel('Population model');
ax = gca;
ax.FontSize = 14;
c = colorbar;
c.Label.String = 'log_{10}(chi-square)';

for itr_indx = 1:length(pop_model)
    subplot(2,5,itr_indx)
    ax = gca;
    b = bar(log10(norm_min_cost(itr_indx,:)));
    b.BarWidth = 0.25;
    xlim([0 length(init_cond)+1]);
    ylim([min(log10(norm_min_cost),[],'all')-0.5 max(log10(norm_min_cost),[],'all')+0.5])
    ax.XTickLabels = (init_cond);
    xlabel('Held out initial conditions');
%     ylabel('log_{10}(fractional change in chi-square)');
    ylabel('log_{10}(chi-square)');

    title(['pop model ' pop_model_names(itr_indx)])
    axis square
    grid on
    ax.FontSize = 12;
end

% cd('C:\Users\Asus\OneDrive - Indian Institute of Science\GitHub\EMT-State-Transition-2.0\Manuscript figures\Figure 2 and related SI')
if(~isempty(find(pop_model > 16,1)))
    print(gcf, ['Cross 1 model validation Yamamoto data b_t_ratio ' num2str(birth_tran_ratio) ' inverse r1 ' num2str(1/r1) '_inverse_r2_' num2str(inverse_r2) '.png'],'-dpng','-r300');
elseif (~isempty(find(pop_model < 16,1)))
    print(gcf, ['Cross 1 model validation Bhatia data b_t_ratio ' num2str(birth_tran_ratio) ' inverse r1 ' num2str(1/r1) '.png'],'-dpng','-r300');
end

figure
ax = gca;
cat_pop_model_names = reordercats(categorical(pop_model_names),pop_model_names);
b = bar(cat_pop_model_names, log10(sum(norm_min_cost,2)));
b.BarWidth = 0.25;
% xlim([0 length(pop_model)+1]);
ylim([min(log10(sum(norm_min_cost,2)),[],'all')-0.5 max(log10(sum(norm_min_cost,2)),[],'all')+0.5])
% ax.XTickLabels = (pop_model);
xlabel('Pop model');
% ylabel('log_{10}(fractional change in chi-square sum)');
ylabel('log_{10}(chi-square sum)');
axis square
grid on
ax.FontSize = 12;

% cd('C:\Users\Asus\OneDrive - Indian Institute of Science\GitHub\EMT-State-Transition-2.0\Manuscript figures\Figure 2 and related SI')
if(~isempty(find(pop_model > 16,1)))
    print(gcf, ['Cross 1 model validation summed Chi square Yamamoto data b_t_ratio ' num2str(birth_tran_ratio) ' inverse r1 ' num2str(1/r1) '_inverse_r2_' num2str(inverse_r2) '.png'],'-dpng','-r300');
elseif (~isempty(find(pop_model < 16,1)))
    print(gcf, ['Cross 1 model validation summed Chi square Bhatia data b_t_ratio ' num2str(birth_tran_ratio) ' inverse r1 ' num2str(1/r1) '.png'],'-dpng','-r300');
end


if(length(init_cond)>2)
    % Control cross 2 validation data analysis %
    init_cond_string = {'1,2' , '1,3' , '1,4' , '1,5' , '2,3' , '2,4' , '2,5' , '3,4' , '3,5' , '4,5'};
    
    init_cond_comb = [];
    for itr_indx = 1:length(init_cond)
        for j = itr_indx+1:length(init_cond)
            init_cond_comb  = [init_cond_comb; itr_indx j];
        end
    end
    
    cd('C:\Users\Asus\OneDrive - Indian Institute of Science\GitHub\EMT-State-Transition-2.0\Parameter sets and generated pseudo data')
    temp_cost = xlsread(['Control_cross_2_validation_chi_sqaure_values_Yamamoto_Exp_data_b_t_ratio_' num2str(birth_tran_ratio) '_inverse_r1_' num2str(1/r1)  '.csv']);
    control_cost = zeros(length(pop_model),size(init_cond_comb,1),total_itr);
    % rearranging cost matrix
    for itr_indx = 1:total_itr
        control_cost(:,:,itr_indx) = temp_cost(1+(itr_indx-1)*length(pop_model):itr_indx*length(pop_model),:);
    end
    
    % cross 2 validation data analysis %
    
    cd('C:\Users\Asus\OneDrive - Indian Institute of Science\GitHub\EMT-State-Transition-2.0\Parameter sets and generated pseudo data')
    temp_cost = xlsread(['cross_2_validation_chi_sqaure_values_Yamamoto_Exp_data_b_t_ratio_' num2str(birth_tran_ratio) '_inverse_r1_' num2str(1/r1)  '.csv']);
    cost = zeros(length(pop_model),size(init_cond_comb,1),total_itr);
    % rearranging cost matrix
    for itr_indx = 1:total_itr
        cost(:,:,itr_indx) = temp_cost(1+(itr_indx-1)*length(pop_model):itr_indx*length(pop_model),:);
    end
    
    mean_norm_cost = mean(cost./control_cost,3);
    
    figure('units','normalized','Position', [0 0 1 1]);
    for itr_indx = 1:length(pop_model)
        subplot(2,3,itr_indx)
        ax = gca;
        bar(log10(mean_norm_cost(itr_indx,:)));
        ax.XTickLabels = categorical(init_cond_string);
        xlabel('test initial conditions');
        ylabel('log_{10}(fractional change in chi-square)');
        title(['pop model ' num2str(pop_model(itr_indx))])
        ylim([min(log10(mean_norm_cost),[],'all')-0.5 max(log10(mean_norm_cost),[],'all')+0.5])

        axis square
        grid on
        ax.FontSize = 12;
        
    end
    
    print(gcf, ['Cross 2 model validation Yamamoto data b_t_ratio ' num2str(birth_tran_ratio) ' inverse r1 ' num2str(1/r1) '.png'],'-dpng','-r300');


    figure
    ax = gca;
    b = bar(log10(sum(mean_norm_cost,2)));
    b.BarWidth = 0.25;
    xlim([0 length(pop_model)+1]);
    ylim([min(log10(sum(mean_norm_cost,2)),[],'all')-0.5 max(log10(sum(mean_norm_cost,2)),[],'all')+0.5])
    ax.XTickLabels = (pop_model);
    xlabel('test initial conditions');
    ylabel('log10(fractional change in chi-square sum)');
    axis square
    grid on
    ax.FontSize = 12;
    
    print(gcf, ['Cross 2 model validation summed Chi square Yamamoto data b_t_ratio ' num2str(birth_tran_ratio) ' inverse r1 ' num2str(1/r1) '.png'],'-dpng','-r300');

end