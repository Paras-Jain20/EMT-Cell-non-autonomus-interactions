clear all

% for Yamamoto data
pop_model = [33];%[20 26 27 25 21 23];%[6 3 11 12 14 10 13];
num_rep = 3;
num_init = 5; % data generated before 15th July considered value 2 for this parameter
num_timept = 5;

% for Bhatia data
% pop_model = 53 ; %[6 3 4 11 12 14 10 13];
% num_rep = 3;
% num_init = 2;
% num_timept = 4;

% para_name = {'r_e', 't_{me}' , 't_{em}', '\alpha', '\beta', 'K_m', 'K_e'};
% para_name = {'r_e', 't_{me}' , 't_{em}', '\alpha', '\beta', 'K_m', 'K_e'};
% para_name = {'r_e', 't_{me}' , 't_{em}', 'K', 'K_e'};

para_name = {'r_e','K', '\sigma', '\mu','t_{me}' , 't_{em}' };

% para_name = {'r_{e}(r_{h})','t_{mh}', 't_{hm}', 't_{he}', 't_{eh}', 'K' ,'H_{E}'};

num_pop_model = length(pop_model);
r1 = 1/54;
inverse_r2 = 35;
r2 = 1/(inverse_r2*r1);
data_src = {'Bhatia', 'Yamamoto'};
profile_sample_size = 2000;

post_itr_samp = 0;

profile_likelihood_data = 'C:\Users\user\OneDrive - Indian Institute of Science\GitHub\EMT-State-Transition-2.0\Parameters sets from experimental data';
results = 'C:\Users\user\OneDrive - Indian Institute of Science\GitHub\EMT-State-Transition-2.0\Manuscript documents\Revision\Figures';

for data_src_indx = 2%1:length(data_src)
    for birth_tran_ratio = 250%[5 10 20 50 100]%[10 50 100 500 1000]

        close all

        for pop_model_indx = 1:length(pop_model)

            if(isempty(r2))
                if(post_itr_samp == 0)
                    cd(profile_likelihood_data)
%                     para_data = table2array(readtable(['profile_lik_para_' data_src{data_src_indx} '_Exp_data_pd_mod_' num2str(pop_model(pop_model_indx)) '_b_t_ratio_' num2str(birth_tran_ratio) '_inverse_r1_' num2str(1/r1) '_para_steps_' num2str(profile_sample_size) '.csv']));
                    para_data = table2array(readtable(['profile_lik_para_data_pd_mod_' num2str(pop_model(pop_model_indx)) '_b_t_ratio_' num2str(birth_tran_ratio) '_inverse_r1_' num2str(1/r1) '_para_steps_' num2str(profile_sample_size) '.csv']));

                else
                    cd('C:\Users\user\OneDrive - Indian Institute of Science\GitHub\EMT-State-Transition-2.0\code files ram\data post sampling')
                    para_data = table2array(readtable(['profile_lik_para_post_itr_samp_pd_mod_' num2str(pop_model(pop_model_indx)) '_b_t_ratio_' num2str(birth_tran_ratio) '_inverse_r1_' num2str(1/r1) '_para_steps_' num2str(profile_sample_size) '.csv']));
                end
            else
                if(post_itr_samp == 0)
                    cd(profile_likelihood_data)
                    para_data = table2array(readtable(['profile_lik_para_' data_src{data_src_indx} '_Exp_data_pd_mod_' num2str(pop_model(pop_model_indx)) '_b_t_ratio_' num2str(birth_tran_ratio) '_inverse_r1_' num2str(1/r1) '_inverse_r2_' num2str(inverse_r2) '_para_steps_' num2str(profile_sample_size) '.csv']));
                    % para_data = table2array(readtable(['profile_lik_para_data_pd_mod_' num2str(pop_model(pop_model_indx)) '_b_t_ratio_' num2str(birth_tran_ratio) '_inverse_r1_' num2str(1/r1) '_inverse_r2_' num2str(inverse_r2) '_para_steps_' num2str(profile_sample_size) '.csv']));

                else
                    % cd('C:\Users\user\OneDrive - Indian Institute of Science\GitHub\EMT-State-Transition-2.0\code files ram\data post sampling')
                    % para_data = table2array(readtable(['profile_lik_para_post_itr_samp_pd_mod_' num2str(pop_model(pop_model_indx)) '_b_t_ratio_' num2str(birth_tran_ratio) '_inverse_r1_' num2str(1/r1) '_inverse_r2_' num2str(inverse_r2) '_para_steps_' num2str(profile_sample_size) '.csv']));
                    cd('C:\Users\user\OneDrive - Indian Institute of Science\GitHub\EMT-State-Transition-2.0\Manuscript documents\Revision\Parameter identifiability improvement')
                    para_data = table2array(readtable(['profile_lik_post_samp_pd_mod_' num2str(pop_model(pop_model_indx)) '_b_t_' num2str(birth_tran_ratio) '_r1_' num2str(1/r1) '_r2_' num2str(inverse_r2) '_para_steps_' num2str(profile_sample_size) '_output_total_cell_population.csv']));

                end
            end
            num_para = size(para_data,2)-2;

            %                 para_data = table2array(readtable(['profile_lik_para_' data_src{data_src_indx} '_Exp_data_pd_mod_' num2str(pop_model(pop_model_indx)) '_b_t_ratio_' num2str(birth_tran_ratio) '_inverse_r1_' num2str(1/r1) '_inverse_r2_' num2str(1/r2)  '_para_steps_' num2str(profile_sample_size) '.csv']));

            %                 para_data = table2array(readtable(['profile_lik_para_' data_src{data_src_indx}  '_Exp_data_pd_mod_' num2str(pop_model(pop_model_indx)) '_b_t_ratio_' num2str(birth_tran_ratio) '_inverse_r1_' num2str(1/r1)  '_dynamic_para_steps_with_default_value_' num2str(profile_sample_size) '.csv']));

            fig1 = figure('Units','normalized','Position',[0 0 1 1]);
            fig2 = figure('Units','normalized','Position',[0 0 1 1]);

            for para_indx = 2:(num_para)
                x_rate_all = para_data(para_data(:,end) == para_indx,1:end-2);
                x_rate = para_data(para_data(:,end) == para_indx,para_indx);

                y = para_data(para_data(:,end) == para_indx,end-1);

                % based on the choosen alpha = 0.954, define the delta chi sqr avlue

                delta_chi_sq_theta = 4;
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
                        case {8,9,22,23,24,33,53,54} % i.e. population model with 6 parameters
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
                        case {21,54} % i.e. population model with 6 parameters
                            delta_chi_sq_all = 12.8;
                    end
                end


                [minima,posn] = min(y); %finds the minimum in y and the corresponding index
                minima_position = x_rate(posn);

                chi_sq_threshold_theta = minima + delta_chi_sq_theta;% the chi_sq threshold for single parameter
                chi_sq_threshold_all = minima + delta_chi_sq_all;%the chi_sq threshold for all params

                %finding the confidence range

                left_indx = find((y(1:posn) - chi_sq_threshold_all > 0)); % here, left indx is for rate para

                right_indx = find((y(posn:end) - chi_sq_threshold_all > 0),1); % here, right indx is for rate para

                if(~isempty(left_indx))

                    left_indx = left_indx(end)+0;

                else
                    left_indx = 1;

                end

                if(~isempty(right_indx))
                    right_indx = posn+right_indx-1;

                else
                    right_indx = length(y);


                end
                figure(fig1);
                para_indx_seq = 1;
                for i = 1:num_para
                    subplot((num_para),(num_para), ((num_para) - para_indx) *((num_para)) + para_indx_seq)
                    %                         plot((x_rate(left_indx+1:right_indx)),(x_rate_all(left_indx+1:right_indx,i) - x_rate_all(left_indx:right_indx-1,i))./x_rate_all(left_indx:right_indx-1,i),'LineWidth',1.5);
                    plot(x_rate(left_indx:right_indx),x_rate_all(left_indx:right_indx,i),'LineWidth',2);

                    %                         xlabel(['Para ' num2str(para_indx) ])
                    %                         ylabel(['Para ' num2str(i)])
                    xlabel([para_name{para_indx} ])
                    ylabel([para_name{i}])
                    axis square
                    grid on
                    para_indx_seq = para_indx_seq + 1;
                    xlim([min(x_rate(left_indx:right_indx))-0.01*min(x_rate(left_indx:right_indx))  max(x_rate(left_indx:right_indx))+0.01*max(x_rate(left_indx:right_indx))])
                    ylim([min(x_rate_all(left_indx:right_indx,i))-0.01*min(x_rate_all(left_indx:right_indx,i))  max(x_rate_all(left_indx:right_indx,i))+0.01*max(x_rate_all(left_indx:right_indx,i))])
                end

                figure(fig2);
                subplot(2,4,para_indx)
                ax= gca;
                plot(x_rate(left_indx:right_indx),y(left_indx:right_indx),'LineWidth',4)
                hold on
                plot(x_rate(left_indx:right_indx),chi_sq_threshold_all*ones(size(x_rate(left_indx:right_indx))),'--k','LineWidth',4)
                xlabel(para_name{para_indx} )
                ylabel(['chi-sqr value'])
                ylim([min(y(left_indx:right_indx))-0.01*min(y(left_indx:right_indx))  max(y(left_indx:right_indx))+0.01*max(y(left_indx:right_indx))])
                ax.FontSize = 24.5;
                axis square
                grid on

            end


            cd(results)
            % saveas(fig1,['Parameter Covariabilty pop model ' num2str(pop_model(pop_model_indx)) ' b_t_ratio ' num2str(birth_tran_ratio) ' inverse r1 ' num2str(1/r1) ' inverse r2 ' num2str(1/r2) '.png']);
            % saveas(fig2,['Profile Likelihood pop model ' num2str(pop_model(pop_model_indx)) ' b_t_ratio ' num2str(birth_tran_ratio) ' inverse r1 ' num2str(1/r1) ' inverse r2 ' num2str(1/r2) '.png']);
            if(isempty(r2))
                if(post_itr_samp == 0)
                    print(fig1,['Parameter Covariabilty pop model ' num2str(pop_model(pop_model_indx)) ' b_t_ratio ' num2str(birth_tran_ratio) ' inverse r1 ' num2str(1/r1) '_para_steps_' num2str(profile_sample_size)  '.png'],'-dpng','-r300');
                    print(fig2,['Profile Likelihood pop model ' num2str(pop_model(pop_model_indx)) ' b_t_ratio ' num2str(birth_tran_ratio) ' inverse r1 ' num2str(1/r1) '_para_steps_' num2str(profile_sample_size) '.png'],'-dpng','-r300');
                else
                    print(fig1,['Parameter Covariabilty post itr samp  pop model ' num2str(pop_model(pop_model_indx)) ' b_t_ratio ' num2str(birth_tran_ratio) ' inverse r1 ' num2str(1/r1) '_para_steps_' num2str(profile_sample_size)  '.png'],'-dpng','-r300');
                    print(fig2,['Profile Likelihood post itr samp pop model ' num2str(pop_model(pop_model_indx)) ' b_t_ratio ' num2str(birth_tran_ratio) ' inverse r1 ' num2str(1/r1) '_para_steps_' num2str(profile_sample_size) '.png'],'-dpng','-r300');
                end
            else
                if(post_itr_samp == 0)
                    print(fig1,['Parameter Covariabilty pop model ' num2str(pop_model(pop_model_indx)) ' b_t_ratio ' num2str(birth_tran_ratio) ' inverse r1 ' num2str(1/r1) ' inverse r2 ' num2str(inverse_r2) '_para_steps_' num2str(profile_sample_size)  '.png'],'-dpng','-r300');
                    print(fig2,['Profile Likelihood pop model ' num2str(pop_model(pop_model_indx)) ' b_t_ratio ' num2str(birth_tran_ratio) ' inverse r1 ' num2str(1/r1) ' inverse r2 ' num2str(inverse_r2) '_para_steps_' num2str(profile_sample_size) '.png'],'-dpng','-r300');
                else
                    print(fig1,['Parameter Covariabilty post itr samp  pop model ' num2str(pop_model(pop_model_indx)) ' b_t_ratio ' num2str(birth_tran_ratio) ' inverse r1 ' num2str(1/r1) ' inverse r2 ' num2str(inverse_r2) '_para_steps_' num2str(profile_sample_size)  '.png'],'-dpng','-r300');
                    print(fig2,['Profile Likelihood post itr samp pop model ' num2str(pop_model(pop_model_indx)) ' b_t_ratio ' num2str(birth_tran_ratio) ' inverse r1 ' num2str(1/r1) ' inverse r2 ' num2str(inverse_r2) '_para_steps_' num2str(profile_sample_size) '.png'],'-dpng','-r300');
                end
            end

            % saveas(fig1,['Parameter Covariabilty pop model ' num2str(pop_model(pop_model_indx)) ' b_t_ratio ' num2str(birth_tran_ratio) ' inverse r1 ' num2str(1/r1) '.png']);
            % saveas(fig2,['Profile Likelihood pop model ' num2str(pop_model(pop_model_indx)) ' b_t_ratio ' num2str(birth_tran_ratio) ' inverse r1 ' num2str(1/r1) '.png']);

        end

    end
end
