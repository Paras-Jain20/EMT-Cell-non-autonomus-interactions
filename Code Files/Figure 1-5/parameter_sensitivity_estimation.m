close all
profile_sample_size =1000;
interpolation = 1;
% color_code = flip([20 0 0; 40 0 0; 60 0 0; 80 0 0; 100 0 0;  120 0 0; 140 0 0; 160 0 0; 180 0 0; 200 0 0;  240 0 0;...
%     240 20 20; 240 40 40; 240 60 60; 240 80 80; 240 100 100; 240 120 120; 240 140 140; 240 160 160; 240 180 180; 240 200 200; 240 220 220; 240 240 240]./255,1);
data_src = {'Bhatia', 'Yamamoto'};
num_para_plots = 100;
pre_samp = 1;
r1 = 1/50;
inverse_r2 = 35;
r2 = []; %1/(inverse_r2*r1);
for data_src_indx = 1%:length(data_src)
    for birth_tran_ratio = 50%[10 50 100 500 1000]
        for pop_model_indx = 4%[20 23 25 26 27] %[4 10 11 12 13 14]

            if(isempty(r2))
                switch pop_model_indx
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
                switch pop_model_indx
                    case{6} % i.e. population model with 1 parameter
                        delta_chi_sq_all = 4;
                    case{3} % i.e. population model with 2 parameter
                        delta_chi_sq_all = 6.17;
                    case {4,5,13,14,15,16,17,19,29} % i.e. population model with 4 parameters
                        delta_chi_sq_all = 9.70;
                    case {1,2,10,11,12,20} % i.e. population model with 3 parameters
                        delta_chi_sq_all = 8.02;
                    case {8,9,22,23,24,33} % i.e. population model with 5 parameters
                        delta_chi_sq_all = 11.3;
                    case {21} % i.e. population model with 6 parameters
                        delta_chi_sq_all = 12.8;
                end
            end



            if(pop_model_indx < 15)
                % importing bhatia et al. 2019 data
                initial_dist = {'high', 'low'; 'low' , 'high'}; % for Bhatia et al data; Mes:Epi
                % High Low fraction: set ini_dist_indx = 1
                % Low High fraction: set ini_dist_indx = 2
                data = [];

                if(pre_samp == 1)
                    cd('C:\Users\Asus\OneDrive - Indian Institute of Science\Projects\State transition\Experimental data');

                    for ini_dist_indx = 1:2

                        temp_data = readmatrix(['Bhatia et al 2019\M_' initial_dist{ini_dist_indx,1} '_E_' initial_dist{ini_dist_indx,2} '_dynamics.xlsx']);

                        temp_data = [temp_data(:,1) temp_data(:,2:end)/100]; % changing population percentage into fractions
                        data = [data;temp_data];


                    end

                    num_rep = 3;
                    num_init_cond = 2; % number of different initial condition in the data
                    init_cond = 1:num_init_cond;
                    timept = (data(1:(size(data,1)/(num_rep*num_init_cond)),1) * 7* 24*r1)'; % time is unitless as it is scaled by m cells division rate (r1)
                    interpolated_timept = 7*12*r1:7*12*r1:timept(end); % interpolating time to get more data points to plot ODE solution trajectories


                else

                    cd('C:\Users\Asus\OneDrive - Indian Institute of Science\GitHub\EMT-State-Transition-2.0\code files ram\data post sampling')
                    data = table2array(readtable('data_post_iterative_sampling_Bhatia_Exp_data_pd_mod_4_b_t_ratio_50_inverse_r1_50_negative_comp'));
                    num_rep = 3;
                    num_init_cond = 2; % number of different initial condition in the data
                    init_cond = 1:num_init_cond;
                    timept = (data(1:(size(data,1)/(num_rep*num_init_cond)),1))'; % time is unitless as it is scaled by m cells division rate (r1)
                    interpolated_timept = 7*12*r1:7*12*r1:timept(end); % interpolating time to get more data points to plot ODE solution trajectories
                    % ODE data charactersitics; Since, we are simulating
                    % ODE trajectories for each parameter set in the 95 confidence we don't need to have multiple replicates
                end



            else

                % Importing Yamamoto et al. 2017 data
                if(pre_samp == 1)
                    cd('C:\Users\Asus\OneDrive - Indian Institute of Science\Projects\State transition\Experimental data\Yamamoto et al. 2017')
                    data = table2array(readtable('Yamamoto et al 2017 Data 1 cell fraction.xlsx'));
                    num_rep = 3;
                    num_init_cond = 5; % number of different initial condition in the data
                    init_cond = 1:num_init_cond;
                    timept = (data(1:(size(data,1)/(num_rep*num_init_cond)),1) * 24*r1)'; % time is unitless as it is scaled by m cells division rate (r1)
                    interpolated_timept =0:24*r1:timept(end); % interpolating time to get more data points to plot ODE solution trajectories
                else
                    cd('C:\Users\Asus\OneDrive - Indian Institute of Science\GitHub\EMT-State-Transition-2.0\code files ram\data post sampling')
                    data = table2array(readtable(['data_post_iterative_sampling_Yamamoto_Exp_data_pd_mod_21_b_t_ratio_250_inverse_r1_54_inverse_r2_35']));
                    num_rep = 3;
                    num_init_cond = 5; % number of different initial condition in the data
                    init_cond = 1:num_init_cond;
                    timept = (data(1:(size(data,1)/(num_rep*num_init_cond)),1))'; % time is unitless as it is scaled by m cells division rate (r1)
                    interpolated_timept =0:24*r1:timept(end)+24*r1; % interpolating time to get more data points to plot ODE solution trajectories
                end



            end
            if(isempty(r2))
                para_start_indx = 1;

                if(pre_samp == 1)
                    cd('C:\Users\Asus\OneDrive - Indian Institute of Science\GitHub\EMT-State-Transition-2.0\Parameters sets from experimental data')
                    if(pop_model_indx ~= 4)
                        para_data = table2array(readtable(['profile_lik_para_' data_src{data_src_indx} '_Exp_data_pd_mod_' num2str(pop_model_indx) '_b_t_ratio_' num2str(birth_tran_ratio) '_inverse_r1_' num2str(1/r1) '_para_steps_' num2str(profile_sample_size) '.csv']));

                    else
                        para_data = table2array(readtable(['profile_lik_para_' data_src{data_src_indx} '_Exp_data_pd_mod_' num2str(pop_model_indx) '_b_t_ratio_' num2str(birth_tran_ratio) '_inverse_r1_' num2str(1/r1) '_para_steps_' num2str(profile_sample_size) '_negative_comp.csv']));
                    end
                else
                    cd('C:\Users\Asus\OneDrive - Indian Institute of Science\GitHub\EMT-State-Transition-2.0\code files ram\data post sampling')

                    if(pop_model_indx ~= 4)
                        para_data = table2array(readtable(['profile_lik_para_post_itr_samp_pd_mod_' num2str(pop_model_indx) '_b_t_ratio_' num2str(birth_tran_ratio) '_inverse_r1_' num2str(1/r1) '_para_steps_' num2str(profile_sample_size) '.csv']));

                    else
                        para_data = table2array(readtable(['profile_lik_para_post_itr_samp_pd_mod_' num2str(pop_model_indx) '_b_t_ratio_' num2str(birth_tran_ratio) '_inverse_r1_' num2str(1/r1) '_para_steps_' num2str(profile_sample_size) '_negative_comp.csv']));

                    end

                end
            else
                para_start_indx = 2;

                if(pre_samp == 1)
                    cd('C:\Users\Asus\OneDrive - Indian Institute of Science\GitHub\EMT-State-Transition-2.0\Parameters sets from experimental data')
                    if(pop_model_indx ~= 4)
                        para_data = table2array(readtable(['profile_lik_para_' data_src{data_src_indx} '_Exp_data_pd_mod_' num2str(pop_model_indx) '_b_t_ratio_' num2str(birth_tran_ratio) '_inverse_r1_' num2str(1/r1) '_inverse_r2_' num2str(inverse_r2)  '_para_steps_' num2str(profile_sample_size) '.csv']));

                    else
                        para_data = table2array(readtable(['profile_lik_para_' data_src{data_src_indx} '_Exp_data_pd_mod_' num2str(pop_model_indx) '_b_t_ratio_' num2str(birth_tran_ratio) '_inverse_r1_' num2str(1/r1) '_inverse_r2_' num2str(inverse_r2) '_para_steps_' num2str(profile_sample_size) '_negative_comp.csv']));
                    end
                else
                    cd('C:\Users\Asus\OneDrive - Indian Institute of Science\GitHub\EMT-State-Transition-2.0\code files ram\data post sampling')

                    if(pop_model_indx ~= 4)
                        para_data = table2array(readtable(['profile_lik_para_post_itr_samp_pd_mod_' num2str(pop_model_indx) '_b_t_ratio_' num2str(birth_tran_ratio) '_inverse_r1_' num2str(1/r1) '_inverse_r2_' num2str(inverse_r2) '_para_steps_' num2str(profile_sample_size) '.csv']));

                    else
                        para_data = table2array(readtable(['profile_lik_para_post_itr_samp_pd_mod_' num2str(pop_model_indx) '_b_t_ratio_' num2str(birth_tran_ratio) '_inverse_r1_' num2str(1/r1) '_inverse_r2_' num2str(inverse_r2) '_para_steps_' num2str(profile_sample_size) '_negative_comp.csv']));

                    end

                end
            end
            num_para = size(para_data,2)-2;

            %             fig1 = figure('units','normalized','Position',[0.2052 0.3120 0.5047 0.3611]);
            fig1 = figure('units','normalized','Position',[0 0 1 1]);


            for para_indx = 1:(size(para_data,2)-2)

                cd('C:\Users\Asus\OneDrive - Indian Institute of Science\GitHub\EMT-State-Transition-2.0\Code files')

                x_rate_all = para_data(para_data(:,end) == para_indx,1:end-2);

                x_rate = para_data(para_data(:,end) == para_indx,para_indx);

                y = para_data(para_data(:,end) == para_indx,end-1);

                [minima,posn] = min(y); %finds the minimum in y and the corresponding index

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


                para_indx_step = ceil((right_indx-left_indx)/num_para_plots);

                if(pop_model_indx < 15)
                    ode_sim_data = zeros(length(interpolated_timept)+1,length(left_indx:para_indx_step:right_indx)*num_init_cond*num_rep); % +1 to include additional experimental time point at day 0
                else
                    ode_sim_data = zeros(length(interpolated_timept),length(left_indx:para_indx_step:right_indx)*num_init_cond*num_rep); % +1 to include additional experimental time point at day 0

                end


                if(para_indx_step ~=0)

                    itr_indx = 1; % as j deosn't starts at one.

                    for j = left_indx:para_indx_step:right_indx

                        temp_ode_sim_data = ODE_simulation(pop_model_indx,num_rep,timept, interpolated_timept,interpolation,data,init_cond, x_rate_all(j,:)); % 1 is for num_rep of ODE simulated data

                        if(pop_model_indx > 15)
                            temp_ode_sim_data = temp_ode_sim_data(:,3)./sum(temp_ode_sim_data(:,3:4),2);
                        else
                            temp_ode_sim_data = temp_ode_sim_data(:,1)./sum(temp_ode_sim_data(:,1:2),2);
                        end
                        % concatenating the ODE solutions together for each increment in the para value for which profile likelihood is presently cosnidered

                        % reshaping to timept row wise
                        if(pop_model_indx < 15)
                            ode_sim_data(:,itr_indx:itr_indx + (num_rep*num_init_cond) -1) = [0.999 0.999 0.999 0.001 0.001 0.001; reshape(temp_ode_sim_data,length(interpolated_timept),num_init_cond*num_rep)];
                            %                         ode_sim_data(:,itr_indx:itr_indx + (num_rep*num_init_cond) -1) = [reshape(temp_ode_sim_data,length(interpolated_timept),num_init_cond*num_rep)];

                        else
                            ode_sim_data(:,itr_indx:itr_indx + (num_rep*num_init_cond) -1) = reshape(temp_ode_sim_data,length(interpolated_timept),num_init_cond*num_rep);
                        end

                        itr_indx = itr_indx + num_init_cond*num_rep;

                        %                         if(pop_model_indx < 15)
                        %
                        %                             %                             ode_sim_data(:,color_code_indx) = temp_ode_sim_data(num_rep*length(interpolated_timept)*(ini_cond_indx-1)+(1-1)*length(interpolated_timept)+1: num_rep*length(interpolated_timept)*(ini_cond_indx-1) + 1*length(interpolated_timept),1);
                        %                             plot(interpolated_timept/r1,ode_sim_data(num_rep*length(interpolated_timept)*(ini_cond_indx-1)+(1-1)*length(interpolated_timept)+1: num_rep*length(interpolated_timept)*(ini_cond_indx-1) + 1*length(interpolated_timept),1),'b-','LineWidth',2) % plotting M fraction dynamics
                        %                             hold on
                        %                         else
                        %                             %                             ode_sim_data(:,color_code_indx) = temp_ode_sim_data(num_rep*length(interpolated_timept)*(ini_cond_indx-1)+(1-1)*length(interpolated_timept)+1: num_rep*length(interpolated_timept)*(ini_cond_indx-1) + 1*length(interpolated_timept),3)./(sum(temp_ode_sim_data(num_rep*length(interpolated_timept)*(ini_cond_indx-1)+(1-1)*length(interpolated_timept)+1: num_rep*length(interpolated_timept)*(ini_cond_indx-1) + 1*length(interpolated_timept),3:4),2));
                        %                             plot(interpolated_timept/r1,ode_sim_data(num_rep*length(interpolated_timept)*(ini_cond_indx-1)+(1-1)*length(interpolated_timept)+1: num_rep*length(interpolated_timept)*(ini_cond_indx-1) + 1*length(interpolated_timept),3)./(sum(ode_sim_data(num_rep*length(interpolated_timept)*(ini_cond_indx-1)+(1-1)*length(interpolated_timept)+1: num_rep*length(interpolated_timept)*(ini_cond_indx-1) + 1*length(interpolated_timept),3:4),2)),'b-','LineWidth',2) % plotting M fraction dynamics
                        %                             hold on
                        %                         end

                    end

                    ini_dist_indx =1;
                    for itr_indx = 1:size(ode_sim_data,2)
                        gcf = fig1;
                        subplot(num_init_cond,num_para, para_indx + num_para *(ini_dist_indx-1))
                        if(pop_model_indx < 15)
                            plot([0 interpolated_timept],ode_sim_data(:,itr_indx),'Color',[0.8500 0.3250 0.0980],'LineStyle','-');
                        else
                            plot([interpolated_timept],ode_sim_data(:,itr_indx),'Color',[0.8500 0.3250 0.0980],'LineStyle','-');

                        end
                        hold on
                        if(mod(itr_indx,num_rep) == 0 && mod(itr_indx,num_rep*num_init_cond) ~= 0)
                            ini_dist_indx = ini_dist_indx + 1;
                        elseif(mod(itr_indx,num_rep) == 0 && mod(itr_indx,num_rep*num_init_cond) == 0)
                            ini_dist_indx = 1;
                        end
                    end

                    for ini_dist_indx = 1:num_init_cond
                        gcf = fig1;
                        subplot(num_init_cond,num_para, para_indx + num_para *(ini_dist_indx-1))

                        temp_timept_indx = zeros(length(timept),1);
                        for j = 1:length(timept)
                            if(~isempty(find(abs(timept(j)/(24*r1) - [0 3 6 9 12]) < 0.001,1)))
                                temp_timept_indx(j) = 1;
                            end
                        end

                        for i = 1:num_rep
                            % plotting experimental data one replicate at a time
                            temp_timept = timept(logical(temp_timept_indx));
                            temp_data = data(num_rep*length(timept)*(ini_dist_indx-1)+(i-1)*length(timept)+1: num_rep*length(timept)*(ini_dist_indx-1) + i*length(timept),2);
                            temp_data = temp_data(logical(temp_timept_indx));
                            plot(temp_timept,temp_data,'ko','MarkerSize',7);

                            %                                 plot(timept,data(num_rep*length(timept)*(ini_dist_indx-1)+(i-1)*length(timept)+1: num_rep*length(timept)*(ini_dist_indx-1) + i*length(timept),2),'ko','MarkerSize',7);

                            temp_timept = timept(~logical(temp_timept_indx));
                            temp_data = data(num_rep*length(timept)*(ini_dist_indx-1)+(i-1)*length(timept)+1: num_rep*length(timept)*(ini_dist_indx-1) + i*length(timept),2);
                            temp_data = temp_data(~logical(temp_timept_indx));
                            plot(temp_timept,temp_data,'mo','MarkerSize',7);


                            %                     plot(timept/r1,data(num_rep*length(timept)*(ini_cond_indx-1)+(i-1)*length(timept)+1: num_rep*length(timept)*(ini_cond_indx-1) + i*length(timept),3),'bo','MarkerSize',5);
                        end
                        if(pop_model_indx < 15)
                            ax = gca;
                            ax.FontSize = 15;
                            grid on
                            ax.GridAlpha = 0.5;
                            % ax.YTick = 0:0.2:1;
                            ax.XTick = (0:2:8)*7*24*r1;
                            ax.XTickLabels = (0:2:8);
                            %                             xlabel('time (weeks)')
                            %                             ylabel([{'EpCAM^{low} fraction'}])
                            axis square

                        else
                            ax = gca;
                            ax.FontSize = 15;
                            grid on
                            ax.GridAlpha = 0.5;
                            ax.XTick = (0:3:12)*24*r1;
                            ax.XTickLabels = (0:3:12);
                            %                             xlabel('time (days)')
                            %                             ylabel([{'Venus cells'}, {'EpCAM^{low} fraction'}])
                            axis square

                        end


                    end
                    %                     mean_sim_data = mean(ode_sim_data,2);
                    %                     std_sim_data = std(ode_sim_data,0,2);
                    %                     xconf = [interpolated_timept/r1 interpolated_timept(end:-1:1)/r1];
                    %                     yconf = [(mean_sim_data+std_sim_data)' flip((mean_sim_data-std_sim_data)') ];
                    %                     p = fill(xconf,yconf,'red');
                    %                     p.FaceColor = [1 0.8 0.8];
                    %                     p.EdgeColor = 'none';
                    %                     hold on
                    %                     plot(interpolated_timept/r1,mean_sim_data,'r-')


                    %                     ax = gca;
                    %                     axis square
                    %                     ax.FontSize = 12;
                    %                     grid on
                    %                     ax.GridAlpha = 0.5;
                    %                     % ax.YTick = 0:0.2:1;
                    %                     ylim([0 1])
                    %                     xlabel('time (hrs)')
                    %                     ylabel('Population fraction')
                    %                     title(['Para sense Initial Cond ' num2str(ini_cond_indx)])
                    %                     directory = pwd;
                    %                     cd('C:\Users\Asus\OneDrive - Indian Institute of Science\GitHub\EMT-State-Transition-2.0\Experimental Data Results')
                    %                     saveas(gcf,['Para ' num2str(para_indx) ' sense pop_model ' num2str(pop_model_indx) ' b_t_ratio ' num2str(birth_tran_ratio) ' inverse r1 ' num2str(1/r1) ' ini_cond ' num2str(ini_cond_indx) '.png'])
                    %                     cd(directory)
                end
            end

%             cd('C:\Users\Asus\OneDrive - Indian Institute of Science\GitHub\EMT-State-Transition-2.0\Manuscript figures\Figure 3 and related SI')
cd('C:\Users\Asus\OneDrive - Indian Institute of Science\GitHub\EMT-State-Transition-2.0\Manuscript documents\Revision\Figures')
            if(isempty(r2))
                if(pre_samp == 1)
                    saveas(gcf,['Parameter Senstivity pre sampling pop model ' num2str(pop_model_indx) ' b t ratio ' num2str(birth_tran_ratio) ' inverse r1 ' num2str(1/r1) '_para_steps_' num2str(profile_sample_size) '.png']);
                else
                    saveas(gcf,['Parameter Senstivity post sampling pop model ' num2str(pop_model_indx) ' b t ratio ' num2str(birth_tran_ratio) ' inverse r1 ' num2str(1/r1) '_para_steps_' num2str(profile_sample_size) '.png']);

                end
            else
                if(pre_samp == 1)
                    saveas(gcf,['Parameter Senstivity pre sampling pop model ' num2str(pop_model_indx) ' b t ratio ' num2str(birth_tran_ratio) ' inverse r1 ' num2str(1/r1) '_inverse_r2_' num2str(inverse_r2) '_para_steps_' num2str(profile_sample_size) '.png']);
                else
                    saveas(gcf,['Parameter Senstivity post sampling pop model ' num2str(pop_model_indx) ' b t ratio ' num2str(birth_tran_ratio) ' inverse r1 ' num2str(1/r1) '_inverse_r2_' num2str(inverse_r2) '_para_steps_' num2str(profile_sample_size) '.png']);

                end
            end




        end
    end
end
