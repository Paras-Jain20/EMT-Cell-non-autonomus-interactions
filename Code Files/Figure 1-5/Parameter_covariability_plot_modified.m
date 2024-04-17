clear all

% for Yamamoto data
% pop_model = 3%[20 26 27 25 21 23];%[6 3 11 12 14 10 13];
% num_rep = 3;
% num_init = 5; % data generated before 15th July considered value 2 for this parameter
% num_timept = 5;

% for Bhatia data
pop_model = 3 ; %[6 3 4 11 12 14 10 13];
num_rep = 3;
num_init = 2;
num_timept = 4;


num_pop_model = length(pop_model);
r1 = 1/50;
r2 = [];%1/33;
data_src = {'Bhatia', 'Yamamoto'};
profile_sample_size = 1000;
for data_src_indx = 1%1:length(data_src)
    for birth_tran_ratio = 100%[5 10 20 50 100]%[10 50 100 500 1000]

        close all

        for pop_model_indx = 1:length(pop_model)

                cd('C:\Users\Asus\OneDrive - Indian Institute of Science\GitHub\EMT-State-Transition-2.0\Parameters sets from experimental data')
                                profile_para = table2array(readtable(['profile_lik_para_' data_src{data_src_indx} '_Exp_data_pd_mod_' num2str(pop_model(pop_model_indx)) '_b_t_ratio_' num2str(birth_tran_ratio) '_inverse_r1_' num2str(1/r1) '_para_steps_' num2str(profile_sample_size) '.csv']));
%                 profile_para = table2array(readtable(['profile_lik_para_' data_src{data_src_indx} '_Exp_data_pd_mod_' num2str(pop_model(pop_model_indx)) '_b_t_ratio_' num2str(birth_tran_ratio) '_inverse_r1_' num2str(1/r1) '_inverse_r2_' num2str(1/r2)  '_para_steps_' num2str(profile_sample_size) '.csv']));
                num_para = size(profile_para,2)-2;
                %                 para_data = table2array(readtable(['profile_lik_para_' data_src{data_src_indx}  '_Exp_data_pd_mod_' num2str(pop_model(pop_model_indx)) '_b_t_ratio_' num2str(birth_tran_ratio) '_inverse_r1_' num2str(1/r1)  '_dynamic_para_steps_with_default_value_' num2str(profile_sample_size) '.csv']));

                fig1 = figure('Units','normalized','Position',[0 0 1 1]);
                fig2 = figure('Units','normalized','Position',[0 0 1 1]);

                                    % based on the choosen alpha = 0.954, define the delta chi sqr avlue

                    delta_chi_sq_theta = 4;
                    if(isempty(r2))
                        switch pop_model(pop_model_indx)
                            case{6} % i.e. population model with 1 parameter
                                delta_chi_sq_all = 4;
                            case {1,2,10,11,12,20} % i.e. population model with 4 parameters
                                delta_chi_sq_all = 9.70;
                            case {3,7,18} % i.e. population model with 3 parameters
                                delta_chi_sq_all = 8.02;
                            case {4,5,13,14,15,16,17,19,29} % i.e. population model with 5 parameters
                                delta_chi_sq_all = 11.3;
                            case {8,9,22,23,24} % i.e. population model with 6 parameters
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
                            case {1,2,10,11,12,20} % i.e. population model with 3 parameters
                                delta_chi_sq_all = 8.02;
                            case {8,9,22,23,24} % i.e. population model with 5 parameters
                                delta_chi_sq_all = 11.3;
                            case {21} % i.e. population model with 6 parameters
                                delta_chi_sq_all = 12.8;
                        end
                    end

                filter_profile_para = [];

                for para_indx = 1:num_para


                    %finding the confidence range

                    confidence_sign = (delta_chi_sq_all*ones(profile_sample_size,1)+min(profile_para(profile_para(:,end) == para_indx,end-1))) - profile_para(profile_para(:,end) == para_indx,end-1);
                    temp_profile_para = profile_para(profile_para(:,end) == para_indx,:);
                    temp_profile_para = temp_profile_para(confidence_sign>=0,:);
                    filter_profile_para = [filter_profile_para; temp_profile_para];

%                     figure(fig1);
%                     para_indx_seq = 1;
%                     for i = 1:(size(para_data,2)-1)
%                         subplot((size(para_data,2)-1),(size(para_data,2)-1), ((size(para_data,2)-1) - para_indx) *((size(para_data,2)-1)) + para_indx_seq)
%                         %                         plot((x_rate(left_indx+1:right_indx)),(x_rate_all(left_indx+1:right_indx,i) - x_rate_all(left_indx:right_indx-1,i))./x_rate_all(left_indx:right_indx-1,i),'LineWidth',1.5);
%                         plot(x_rate(left_indx:right_indx),x_rate_all(left_indx:right_indx,i),'LineWidth',2);
% 
%                         xlabel(['Para ' num2str(para_indx) ])
%                         ylabel(['Para ' num2str(i)])
%                         axis square
%                         grid on
%                         para_indx_seq = para_indx_seq + 1;
%                         xlim([min(x_rate(left_indx:right_indx))-0.01*min(x_rate(left_indx:right_indx))  max(x_rate(left_indx:right_indx))+0.01*max(x_rate(left_indx:right_indx))])
%                         ylim([min(x_rate_all(left_indx:right_indx,i))-0.01*min(x_rate_all(left_indx:right_indx,i))  max(x_rate_all(left_indx:right_indx,i))+0.01*max(x_rate_all(left_indx:right_indx,i))])
%                     end
% 
%                     figure(fig2);
%                     subplot(2,4,para_indx)
%                     ax= gca;
%                     plot(x_rate(left_indx:right_indx),y(left_indx:right_indx),'LineWidth',2)
%                     hold on
%                     plot(x_rate(left_indx:right_indx),chi_sq_threshold_all*ones(size(x_rate(left_indx:right_indx))),'--k','LineWidth',2)
%                     xlabel(['Para ' num2str(para_indx) ])
%                     ylabel(['chi-sqr value'])
%                     ylim([min(y(left_indx:right_indx))-0.01*min(y(left_indx:right_indx))  max(y(left_indx:right_indx))+0.01*max(y(left_indx:right_indx))])
%                     ax.FontSize = 14;
%                     axis square
%                     grid on

                end

                    figure(fig1);
                    for para_indx = 1:num_para
                        for cross_para_indx = para_indx:num_para

                            subplot(num_para,num_para, (num_para - para_indx) *(num_para) + cross_para_indx)
                            %                         plot((x_rate(left_indx+1:right_indx)),(x_rate_all(left_indx+1:right_indx,i) - x_rate_all(left_indx:right_indx-1,i))./x_rate_all(left_indx:right_indx-1,i),'LineWidth',1.5);
                            scatter(filter_profile_para(:,para_indx),filter_profile_para(:,cross_para_indx),"filled");

                            xlabel(['Para ' num2str(para_indx) ])
                            ylabel(['Para ' num2str(cross_para_indx)])
                            axis square
                            grid on
                            xlim([min(filter_profile_para(:,para_indx))-0.01*min(filter_profile_para(:,para_indx))  max(filter_profile_para(:,para_indx))+0.01*max(filter_profile_para(:,para_indx))])
                            ylim([min(filter_profile_para(:,cross_para_indx))-0.01*min(filter_profile_para(:,cross_para_indx))  max(filter_profile_para(:,cross_para_indx))+0.01*max(filter_profile_para(:,cross_para_indx))])

                        end
                    end

%                     figure(fig2);
%                     subplot(2,4,para_indx)
%                     ax= gca;
%                     plot(x_rate(left_indx:right_indx),y(left_indx:right_indx),'LineWidth',2)
%                     hold on
%                     plot(x_rate(left_indx:right_indx),chi_sq_threshold_all*ones(size(x_rate(left_indx:right_indx))),'--k','LineWidth',2)
%                     xlabel(['Para ' num2str(para_indx) ])
%                     ylabel(['chi-sqr value'])
%                     ylim([min(y(left_indx:right_indx))-0.01*min(y(left_indx:right_indx))  max(y(left_indx:right_indx))+0.01*max(y(left_indx:right_indx))])
%                     ax.FontSize = 14;
%                     axis square
%                     grid on

            
            cd('C:\Users\Asus\OneDrive - Indian Institute of Science\GitHub\EMT-State-Transition-2.0\Experimental Data Results')
%             saveas(fig1,['Parameter Covariabilty pop model ' num2str(pop_model(pop_model_indx)) ' b_t_ratio ' num2str(birth_tran_ratio) ' inverse r1 ' num2str(1/r1) ' inverse r2 ' num2str(1/r2) '.png']);
%             saveas(fig2,['Profile Likelihood pop model ' num2str(pop_model(pop_model_indx)) ' b_t_ratio ' num2str(birth_tran_ratio) ' inverse r1 ' num2str(1/r1) ' inverse r2 ' num2str(1/r2) '.png']);

            %             saveas(fig1,['Parameter Covariabilty pop model ' num2str(pop_model(pop_model_indx)) ' b_t_ratio ' num2str(birth_tran_ratio) ' inverse r1 ' num2str(1/r1) '.png']);
            %             saveas(fig2,['Profile Likelihood pop model ' num2str(pop_model(pop_model_indx)) ' b_t_ratio ' num2str(birth_tran_ratio) ' inverse r1 ' num2str(1/r1) '.png']);

        end

    end
end
