
close all
% whether we are using the test or experimental data: test_data = 1; exp_data = 2
data_type = 2;

test_data_indx = 1; % necessary dummy variable to just to fill in variable list in profile likelihood function
noise_factor = 1; % necessary dummy variable to just to fill in variable list in profile likelihood function
r1 = 1/54; % m cells division rates; scaling factor for non dimensional population model
inverse_r2 = 33 ;% r_e is set to 1.5 which is normalised with r1
r2 = 1/(inverse_r2*r1); % normalised r_e; is set to 1.5 which is normalised with r1; if r2 is empty set then it is a variable parameter
default_total_itr = 10;
% r1 = 1/50;
% r2 = [];
interpolation = 0;

% fig1 = figure('units','normalized','Position',[0.5 0.5 0.25 0.25]);
% fig2 = figure('units','normalized','Position',[0.5 0.5 0.25 0.25]);
% fig3 = figure('units','normalized','Position',[0.5 0.5 0.25 0.25]);
% fig4 = figure('units','normalized','Position',[0.5 0.5 0.25 0.25]);
% fig5 = figure('units','normalized','Position',[0.5 0.5 0.25 0.25]);

fig1 = figure('units','normalized','Position',[0.4797 0.4019 0.2703 0.3481]);
fig2 = figure('units','normalized','Position',[0.4797 0.4019 0.2703 0.3481]);
fig3 = figure('units','normalized','Position',[0.4797 0.4019 0.2703 0.3481]);
fig4 = figure('units','normalized','Position',[0.4797 0.4019 0.2703 0.3481]);
fig5 = figure('units','normalized','Position',[0.4797 0.4019 0.2703 0.3481]);


% line_plot_style = {'-' ,'--', '-.' ':'};
line_plot_style = {'-' ,'-', '-' '-'};

Linecolor = {[0.8500 0.3250 0.0980] ,[0 0.4470 0.7410], [0.4940 0.1840 0.5560], [0.4660 0.6740 0.1880], [0.6350 0.0780 0.1840], [0.3010 0.7450 0.9330]};

% orange, blue, purple, green, violet, light blue
% select a population growth model
pop_mod_seq = 1;
for pop_model_indx = [48 49 50]%[3 4 11 12 14 10 13]%[19 29 32 20 21 33]
    for birth_tran_ratio = 200%[5 10 20 100]%[50 100 150 200 250 300] %200 500] %[5 10 20 50 100] % max ratio between birth and transition

        % type of optimum parameters search to be performed:
        % use existing parameters = choose 0
        % optimizer based : choose 1
        % uniform sampling based : choose 2; this method is not updated so might not work
        search_method = 1; % parameter search method

        if(pop_model_indx < 15) % for models capturing dynamics of Bhatia et al. data

            initial_dist = {'high', 'low'; 'low' , 'high'}; % for Bhatia et al data; Mes:Epi

            data = [];

            cd('C:\Users\Asus\OneDrive - Indian Institute of Science\Projects\State transition\Experimental data');

            % importing Bhatia et al data
            for ini_dist_indx = 1:2

                temp_data = table2array(readtable((['Bhatia et al 2019\M_' initial_dist{ini_dist_indx,1} '_E_' initial_dist{ini_dist_indx,2} '_dynamics.xlsx'])));

                temp_data = [temp_data(:,1) temp_data(:,2:end)/100]; % changing population percentage into fractions
                data = [data;temp_data];
            end
            % data charactersitics
            num_rep = 3;
            num_init_cond = 2; % number of different initial condition in the data
            init_cond = 1:num_init_cond;
            timept = (2:2:8)*(7*24)*r1; % time is unitless as it is scaled by m cells division rate (r1)
            interpolated_timept = r1:r1:timept(end); % interpolating time to get more data points to plot ODE solution trajectories

            cd('C:\Users\Asus\OneDrive - Indian Institute of Science\GitHub\EMT-State-Transition-2.0\Parameters sets from experimental data')

            if(isfile(['para_Bhatia_data_pd_mod_' num2str(pop_model_indx) '_s_m_' num2str(search_method) '_b_t_ratio_' num2str(birth_tran_ratio) '_inverse_r1_' num2str(1/r1) '.csv']))
               
                para_data = table2array(readtable(['para_Bhatia_data_pd_mod_' num2str(pop_model_indx) '_s_m_' num2str(search_method) '_b_t_ratio_' num2str(birth_tran_ratio) '_inverse_r1_' num2str(1/r1) '.csv']));

                select_para = (para_data(:,1:end-1));
                goodness = (para_data(:,end));
            else
                cd('C:\Users\Asus\OneDrive - Indian Institute of Science\GitHub\EMT-State-Transition-2.0\Code files')

                [select_para,goodness] = para_search(1,pop_model_indx,5, data_type, data, num_rep, init_cond, timept, interpolated_timept, interpolation, num_time_pts,birth_tran_ratio,r2); % select_para are filtered for re > rm
            end

        else
            %             cd('C:\Users\Asus\OneDrive - Indian Institute of Science\GitHub\EMT-State-Transition-2.0\Parameters sets from experimental data')
            %             para_data = table2array(readtable(['para_Yamamoto_data_pd_mod_' num2str(pop_model_indx) '_s_m_' num2str(search_method) '_b_t_ratio_' num2str(birth_tran_ratio) '_inverse_r1_' num2str(1/r1)  '.csv']));



            % Importing Yamamoto et al. 2017 data; models having
            % pop_model_indx >= 15
            cd('C:\Users\Asus\OneDrive - Indian Institute of Science\Projects\State transition\Experimental data\Yamamoto et al. 2017')
            data = table2array(readtable('Yamamoto et al 2017 Data 1 cell fraction.xlsx'));
            num_rep = 3;
            num_init_cond = 5; % number of different initial condition in the data
            init_cond = 1:num_init_cond;
            timept = (0:3:12)*24*r1; % time is unitless as it is scaled by m cells division rate (r1)
            interpolated_timept = 0:r1:timept(end);  % interpolating time to get more data points to plot ODE solution trajectories
            interpolation = 0;
            cd('C:\Users\Asus\OneDrive - Indian Institute of Science\GitHub\EMT-State-Transition-2.0\Code files')

            [select_para,goodness] = para_search(1,pop_model_indx,5, data_type, data, num_rep, init_cond, timept, interpolated_timept, interpolation, num_time_pts,birth_tran_ratio,r2); % select_para are filtered for re > rm
            %             select_para = (para_data(:,1:end-1));
            %
            %             goodness = (para_data(:,end));
        end

        interpolation = 1;
        num_time_pts = length(timept);

        [~,order] = sort(goodness);
        select_para = select_para(order,:);

        cd('C:\Users\Asus\OneDrive - Indian Institute of Science\GitHub\EMT-State-Transition-2.0\Code files')
        num_rep = 3; % to limit the number of line plots
        ode_sim_data = ODE_simulation(pop_model_indx,num_rep,timept, interpolated_timept, interpolation,data, init_cond, select_para(1,:));

        if((pop_model_indx == 1 && data_type == 2))
            ode_sim_data = ode_sim_data./sum(ode_sim_data,2);
        end

        for ini_cond_indx = init_cond

            % Fractional Dynamics for each initial condition is plotted separately
            switch ini_cond_indx
                case 1
                    figure(fig1)
                case 2
                    figure(fig2)
                case 3
                    figure(fig3)
                case 4
                    figure(fig4)
                case 5
                    figure(fig5)
            end

            for i = 1:num_rep
                if (interpolation == 0) % use experimental time points to plot ODE solutions
                    if(pop_model_indx <= 15)

                        plot(timept/r1,ode_sim_data(num_rep*length(timept)*(ini_cond_indx-1)+(i-1)*length(timept)+1: num_rep*length(timept)*(ini_cond_indx-1) + i*length(timept),1),['r'  Linecolor{pop_mod_seq} line_plot_style{pop_mod_seq}],'LineWidth',2,'MarkerSize',8) % plotting M fraction dynamics
                        hold on
                        %                         plot(timept/r1,ode_sim_data(num_rep*length(timept)*(ini_cond_indx-1)+(i-1)*length(timept)+1: num_rep*length(timept)*(ini_cond_indx-1) + i*length(timept),2),'b','LineStyle', line_plot_style(ini_cond_indx),'LineWidth',2) % plotting E fraction dynamics
                    elseif(pop_model_indx > 15 )
                        plot(timept/r1,ode_sim_data(num_rep*length(timept)*(ini_cond_indx-1)+(i-1)*length(timept)+1: num_rep*length(timept)*(ini_cond_indx-1) + i*length(timept),3)./(sum(ode_sim_data(num_rep*length(timept)*(ini_cond_indx-1)+(i-1)*length(timept)+1: num_rep*length(timept)*(ini_cond_indx-1) + i*length(timept),3:4),2)),['r'  Linecolor{pop_mod_seq} line_plot_style{pop_mod_seq}],'LineWidth',2,'MarkerSize',8) % plotting M fraction dynamics
                        hold on
                        %                         plot(timept/r1,ode_sim_data(num_rep*length(timept)*(ini_cond_indx-1)+(i-1)*length(timept)+1: num_rep*length(timept)*(ini_cond_indx-1) + i*length(timept),4)./(sum(ode_sim_data(num_rep*length(timept)*(ini_cond_indx-1)+(i-1)*length(timept)+1: num_rep*length(timept)*(ini_cond_indx-1) + i*length(timept),3:4),2)),'b','LineStyle', line_plot_style(ini_cond_indx),'LineWidth',2) % plotting E fraction dynamics
                    end
                else % use interpolated time points to plot ODE solutions
                    if(pop_model_indx < 15)

                        plot(interpolated_timept/r1,ode_sim_data(num_rep*length(interpolated_timept)*(ini_cond_indx-1)+(i-1)*length(interpolated_timept)+1: num_rep*length(interpolated_timept)*(ini_cond_indx-1) + i*length(interpolated_timept),1),'Color',Linecolor{pop_mod_seq},'LineWidth',2) % plotting M fraction dynamics
                        hold on
                        %                         plot(interpolated_timept/r1,ode_sim_data(num_rep*length(interpolated_timept)*(ini_cond_indx-1)+(i-1)*length(interpolated_timept)+1: num_rep*length(interpolated_timept)*(ini_cond_indx-1) + i*length(interpolated_timept),2),'Color',Linecolor{pop_mod_seq},'LineWidth',2) % plotting E fraction dynamics
                    elseif(pop_model_indx >= 15 && i == 1) % to plot ODE solution for only one replicate
                        % plotting Venus Labelled EpCAM -ve cells fraction in the total Venus Labelled fraction
                        plot(interpolated_timept/r1,ode_sim_data(num_rep*length(interpolated_timept)*(ini_cond_indx-1)+(i-1)*length(interpolated_timept)+1: num_rep*length(interpolated_timept)*(ini_cond_indx-1) + i*length(interpolated_timept),3)./(sum(ode_sim_data(num_rep*length(interpolated_timept)*(ini_cond_indx-1)+(i-1)*length(interpolated_timept)+1: num_rep*length(interpolated_timept)*(ini_cond_indx-1) + i*length(interpolated_timept),3:4),2)),'Color',Linecolor{pop_mod_seq},'LineWidth',2)
                        hold on
                        %                         plot(interpolated_timept/r1,ode_sim_data(num_rep*length(interpolated_timept)*(ini_cond_indx-1)+(i-1)*length(interpolated_timept)+1: num_rep*length(interpolated_timept)*(ini_cond_indx-1) + i*length(interpolated_timept),4)./(sum(ode_sim_data(num_rep*length(interpolated_timept)*(ini_cond_indx-1)+(i-1)*length(interpolated_timept)+1: num_rep*length(interpolated_timept)*(ini_cond_indx-1) + i*length(interpolated_timept),3:4),2)),'Color',Linecolor{pop_mod_seq},'LineWidth',2) % plotting E fraction dynamics
                    end
                end
                % plotting experimental data one replicate at a time
                plot(timept/r1,data(num_rep*length(timept)*(ini_cond_indx-1)+(i-1)*length(timept)+1: num_rep*length(timept)*(ini_cond_indx-1) + i*length(timept),2),'ko','MarkerSize',10);
                %                 plot(timept/r1,data(num_rep*length(timept)*(ini_cond_indx-1)+(i-1)*length(timept)+1: num_rep*length(timept)*(ini_cond_indx-1) + i*length(timept),3),'bo','MarkerSize',7);
            end

            %     directory = pwd;
            %     cd('C:\Users\Asus\OneDrive - Indian Institute of Science\GitHub\EMT-State-Transition-2.0\Experimental Data Results')
            %     if(pop_model_indx < 15)
            %         saveas(gcf,['Bhatia data Model Fit Ini Cond ' num2str(ini_cond_indx) ' pop mod ' num2str(pop_model_indx) ' b_t_ratio ' num2str(b_t_ratio)  ' r1 ' num2str(1/r1) '.png']);
            %     else
            %         saveas(gcf,['Yamamoto data Model Fit Ini Cond ' num2str(ini_cond_indx) ' pop mod ' num2str(pop_model_indx) ' b_t_ratio ' num2str(b_t_ratio)  ' r1 ' num2str(1/r1) '.png']);
            %     end
            %     cd(directory)
        end
    end
    pop_mod_seq = pop_mod_seq + 1;
end

cd('C:\Users\Asus\OneDrive - Indian Institute of Science\GitHub\EMT-State-Transition-2.0\Experimental Data Results')
figure(fig1)
ax = gca;
axis square
ax.FontSize = 14;
grid on
ax.GridAlpha = 0.5;
% ax.YTick = 0:0.2:1;
% ax.XTick = (0:2:8)*7*24;
% ax.XTickLabels = (0:2:8);
% xlabel('time (weeks)')

ax.XTick = (0:3:12)*24;
ax.XTickLabels = (0:3:12);
xlabel('time (days)')

ylabel([{'Venus cells'}, {'EpCAM^{low} fraction'}])
saveas(gcf,['Growth and transition Influence Yamamoto data Model Fit Ini Cond ' num2str(1)  ' b_t_ratio ' num2str(birth_tran_ratio)  ' r1 ' num2str(1/r1) '.png']);

% title(['Data and System dynamics Initial Cond ' num2str(1)])

figure(fig2)
ax = gca;
axis square
ax.FontSize = 14;
grid on
ax.GridAlpha = 0.5;
% ax.YTick = 0:0.2:1;
% ax.XTick = (0:2:8)*7*24;
% ax.XTickLabels = (0:2:8);
% xlabel('time (weeks)')

ax.XTick = (0:3:12)*24;
ax.XTickLabels = (0:3:12);
xlabel('time (days)')
ylabel([{'Venus cells'}, {'EpCAM^{low} fraction'}])
saveas(gcf,['Growth and transition Influence Yamamoto data Model Fit Ini Cond ' num2str(2)  ' b_t_ratio ' num2str(birth_tran_ratio)  ' r1 ' num2str(1/r1) '.png']);

% title(['Data and System dynamics Initial Cond ' num2str(1)])

figure(fig3)
ax = gca;
axis square
ax.FontSize = 14;
grid on
ax.GridAlpha = 0.5;
% ax.YTick = 0:0.2:1;
% ax.XTick = (0:2:8)*7*24;
% ax.XTickLabels = (0:2:8);
% xlabel('time (weeks)')
ax.XTick = (0:3:12)*24;
ax.XTickLabels = (0:3:12);
xlabel('time (days)')
ylabel([{'Venus cells'}, {'EpCAM^{low} fraction'}])
saveas(gcf,['Growth and transition Influence Yamamoto data Model Fit Ini Cond ' num2str(3)  ' b_t_ratio ' num2str(birth_tran_ratio)  ' r1 ' num2str(1/r1) '.png']);


figure(fig4)
ax = gca;
axis square
ax.FontSize = 14;
grid on
ax.GridAlpha = 0.5;
% ax.YTick = 0:0.2:1;
% ax.XTick = (0:2:8)*7*24;
% ax.XTickLabels = (0:2:8);
% xlabel('time (weeks)')
ax.XTick = (0:3:12)*24;
ax.XTickLabels = (0:3:12);
xlabel('time (days)')
ylabel([{'Venus cells'}, {'EpCAM^{low} fraction'}])
saveas(gcf,['Growth and transition Influence Yamamoto data Model Fit Ini Cond ' num2str(4)  ' b_t_ratio ' num2str(birth_tran_ratio)  ' r1 ' num2str(1/r1) '.png']);

figure(fig5)
ax = gca;
axis square
ax.FontSize = 14;
grid on
ax.GridAlpha = 0.5;
% ax.YTick = 0:0.2:1;
% ax.XTick = (0:2:8)*7*24;
% ax.XTickLabels = (0:2:8);
% xlabel('time (weeks)')
ax.XTick = (0:3:12)*24;
ax.XTickLabels = (0:3:12);
xlabel('time (days)')
ylabel([{'Venus cells'}, {'EpCAM^{low} fraction'}])
saveas(gcf,['Growth and transition Influence Yamamoto data Model Fit Ini Cond ' num2str(5)  ' b_t_ratio ' num2str(birth_tran_ratio)  ' r1 ' num2str(1/r1) '.png']);

% title(['Data and System dynamics Initial Cond ' num2str(2)])
