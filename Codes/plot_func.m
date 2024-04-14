
function plot_func(select_para,goodness,data_type, data,pop_model_indx,num_rep,timept, interpolated_timept, interpolation,r1,r2,inverse_r2, init_cond,b_t_ratio)

[~,order] = sort(goodness);
select_para = select_para(order,:);
for order_indx = 1:length(order)
    ode_sim_data = ODE_simulation(pop_model_indx,num_rep,timept, interpolated_timept, interpolation,data, init_cond, select_para(order_indx,:));

    if((pop_model_indx == 1 && data_type == 2))
        ode_sim_data = ode_sim_data./sum(ode_sim_data,2);
    end

    Linecolor = {[0.8500 0.3250 0.0980] ,[0 0.4470 0.7410], [0.4940 0.1840 0.5560], [0.4660 0.6740 0.1880]};
    for ini_cond_indx = init_cond
        % Fractional Dynamics for each initial condition is plotted separately
        figure('units','normalized','Position',[0.5000 0.5000 0.2089 0.2500]);
        for i = 1:num_rep
            if (interpolation == 0) % use experimental time points to plot ODE solutions
                if(pop_model_indx < 15)

                    plot(timept/r1,ode_sim_data(num_rep*length(timept)*(ini_cond_indx-1)+(i-1)*length(timept)+1: num_rep*length(timept)*(ini_cond_indx-1) + i*length(timept),1),'r-','LineWidth',2) % plotting M fraction dynamics
                    hold on
                    plot(timept/r1,ode_sim_data(num_rep*length(timept)*(ini_cond_indx-1)+(i-1)*length(timept)+1: num_rep*length(timept)*(ini_cond_indx-1) + i*length(timept),2),'b-','LineWidth',2) % plotting E fraction dynamics
                elseif(pop_model_indx >= 15 && pop_model_indx < 28)
                    plot(timept/r1,ode_sim_data(num_rep*length(timept)*(ini_cond_indx-1)+(i-1)*length(timept)+1: num_rep*length(timept)*(ini_cond_indx-1) + i*length(timept),3)./(sum(ode_sim_data(num_rep*length(timept)*(ini_cond_indx-1)+(i-1)*length(timept)+1: num_rep*length(timept)*(ini_cond_indx-1) + i*length(timept),3:4),2)),'r-','LineWidth',2) % plotting M fraction dynamics
                    hold on
                    plot(timept/r1,ode_sim_data(num_rep*length(timept)*(ini_cond_indx-1)+(i-1)*length(timept)+1: num_rep*length(timept)*(ini_cond_indx-1) + i*length(timept),4)./(sum(ode_sim_data(num_rep*length(timept)*(ini_cond_indx-1)+(i-1)*length(timept)+1: num_rep*length(timept)*(ini_cond_indx-1) + i*length(timept),3:4),2)),'b-','LineWidth',2) % plotting E fraction dynamics
                else
                    plot(timept/r1,ode_sim_data(num_rep*length(timept)*(ini_cond_indx-1)+(i-1)*length(timept)+1: num_rep*length(timept)*(ini_cond_indx-1) + i*length(timept),1)./sum(ode_sim_data(num_rep*length(timept)*(ini_cond_indx-1)+(i-1)*length(timept)+1: num_rep*length(timept)*(ini_cond_indx-1) + i*length(timept),:),2),'r-','LineWidth',2) % plotting M fraction dynamics
                    hold on
                    plot(timept/r1,ode_sim_data(num_rep*length(timept)*(ini_cond_indx-1)+(i-1)*length(timept)+1: num_rep*length(timept)*(ini_cond_indx-1) + i*length(timept),2)./sum(ode_sim_data(num_rep*length(timept)*(ini_cond_indx-1)+(i-1)*length(timept)+1: num_rep*length(timept)*(ini_cond_indx-1) + i*length(timept),:),2),'b-','LineWidth',2) % plotting E fraction dynamics
                end
            else % use interpolated time points to plot ODE solutions
                if(pop_model_indx < 15)

                    plot(interpolated_timept/r1,ode_sim_data(num_rep*length(interpolated_timept)*(ini_cond_indx-1)+(i-1)*length(interpolated_timept)+1: num_rep*length(interpolated_timept)*(ini_cond_indx-1) + i*length(interpolated_timept),1),'r-','LineWidth',2) % plotting M fraction dynamics
                    hold on
                    %                 plot(interpolated_timept/r1,ode_sim_data(num_rep*length(interpolated_timept)*(ini_cond_indx-1)+(i-1)*length(interpolated_timept)+1: num_rep*length(interpolated_timept)*(ini_cond_indx-1) + i*length(interpolated_timept),2),'b-','LineWidth',2) % plotting E fraction dynamics
                elseif(pop_model_indx >= 15 && pop_model_indx ~=28)
%                     % plotting total cell number with time
%                     %                 plot(interpolated_timept/r1,sum(ode_sim_data(num_rep*length(interpolated_timept)*(ini_cond_indx-1)+(i-1)*length(interpolated_timept)+1: num_rep*length(interpolated_timept)*(ini_cond_indx-1) + i*length(interpolated_timept),[1 3]),2)./(sum(ode_sim_data(num_rep*length(interpolated_timept)*(ini_cond_indx-1)+(i-1)*length(interpolated_timept)+1: num_rep*length(interpolated_timept)*(ini_cond_indx-1) + i*length(interpolated_timept),1:4),2)),'r-','LineWidth',2)
%                     plot(interpolated_timept/r1,sum(ode_sim_data(num_rep*length(interpolated_timept)*(ini_cond_indx-1)+(i-1)*length(interpolated_timept)+1: num_rep*length(interpolated_timept)*(ini_cond_indx-1) + i*length(interpolated_timept),:),2),'-','Color', Linecolor{1},'LineWidth',3)
%                     hold on
%                     % plotting EpCAM low cell number with time
%                     plot(interpolated_timept/r1,sum(ode_sim_data(num_rep*length(interpolated_timept)*(ini_cond_indx-1)+(i-1)*length(interpolated_timept)+1: num_rep*length(interpolated_timept)*(ini_cond_indx-1) + i*length(interpolated_timept),[1 3]),2),'-','Color', Linecolor{2},'LineWidth',3)
%                     % plotting EpCAM high cell number with time
%                     plot(interpolated_timept/r1,sum(ode_sim_data(num_rep*length(interpolated_timept)*(ini_cond_indx-1)+(i-1)*length(interpolated_timept)+1: num_rep*length(interpolated_timept)*(ini_cond_indx-1) + i*length(interpolated_timept),[2 4]),2),'-','Color', Linecolor{3},'LineWidth',3)
%                     hold on

                    % plotting total cell number with time
                    %                 plot(interpolated_timept/r1,sum(ode_sim_data(num_rep*length(interpolated_timept)*(ini_cond_indx-1)+(i-1)*length(interpolated_timept)+1: num_rep*length(interpolated_timept)*(ini_cond_indx-1) + i*length(interpolated_timept),[1 3]),2)./(sum(ode_sim_data(num_rep*length(interpolated_timept)*(ini_cond_indx-1)+(i-1)*length(interpolated_timept)+1: num_rep*length(interpolated_timept)*(ini_cond_indx-1) + i*length(interpolated_timept),1:4),2)),'r-','LineWidth',2)
%                     plot(interpolated_timept/r1,sum(ode_sim_data(num_rep*length(interpolated_timept)*(ini_cond_indx-1)+(i-1)*length(interpolated_timept)+1: num_rep*length(interpolated_timept)*(ini_cond_indx-1) + i*length(interpolated_timept),[1 3]),2)./(sum(ode_sim_data(num_rep*length(interpolated_timept)*(ini_cond_indx-1)+(i-1)*length(interpolated_timept)+1: num_rep*length(interpolated_timept)*(ini_cond_indx-1) + i*length(interpolated_timept),1:4),2)),'r-','LineWidth',2)
%                     hold on
%                     plotting EpCAM low cell fraction with time
                    plot(interpolated_timept/r1,sum(ode_sim_data(num_rep*length(interpolated_timept)*(ini_cond_indx-1)+(i-1)*length(interpolated_timept)+1: num_rep*length(interpolated_timept)*(ini_cond_indx-1) + i*length(interpolated_timept),[1 3]),2)./(sum(ode_sim_data(num_rep*length(interpolated_timept)*(ini_cond_indx-1)+(i-1)*length(interpolated_timept)+1: num_rep*length(interpolated_timept)*(ini_cond_indx-1) + i*length(interpolated_timept),1:4),2)),'-','Color', Linecolor{2},'LineWidth',3)
%                     plotting EpCAM high cell fraction with time
                    plot(interpolated_timept/r1,sum(ode_sim_data(num_rep*length(interpolated_timept)*(ini_cond_indx-1)+(i-1)*length(interpolated_timept)+1: num_rep*length(interpolated_timept)*(ini_cond_indx-1) + i*length(interpolated_timept),[2 4]),2)./(sum(ode_sim_data(num_rep*length(interpolated_timept)*(ini_cond_indx-1)+(i-1)*length(interpolated_timept)+1: num_rep*length(interpolated_timept)*(ini_cond_indx-1) + i*length(interpolated_timept),1:4),2)),'-','Color', Linecolor{3},'LineWidth',3)
                    hold on

                    % plot(interpolated_timept/r1,ode_sim_data(num_rep*length(interpolated_timept)*(ini_cond_indx-1)+(i-1)*length(interpolated_timept)+1: num_rep*length(interpolated_timept)*(ini_cond_indx-1) + i*length(interpolated_timept),4)./(sum(ode_sim_data(num_rep*length(interpolated_timept)*(ini_cond_indx-1)+(i-1)*length(interpolated_timept)+1: num_rep*length(interpolated_timept)*(ini_cond_indx-1) + i*length(interpolated_timept),3:4),2)),'b-','LineWidth',2) % plotting E fraction dynamics
                else
                    plot(interpolated_timept/r1,ode_sim_data(num_rep*length(interpolated_timept)*(ini_cond_indx-1)+(i-1)*length(interpolated_timept)+1: num_rep*length(interpolated_timept)*(ini_cond_indx-1) + i*length(interpolated_timept),1)./sum(ode_sim_data(num_rep*length(interpolated_timept)*(ini_cond_indx-1)+(i-1)*length(interpolated_timept)+1: num_rep*length(interpolated_timept)*(ini_cond_indx-1) + i*length(interpolated_timept),:),2),'r-','LineWidth',2) % plotting M fraction dynamics
                    hold on
                    %                 plot(interpolated_timept/r1,ode_sim_data(num_rep*length(interpolated_timept)*(ini_cond_indx-1)+(i-1)*length(interpolated_timept)+1: num_rep*length(interpolated_timept)*(ini_cond_indx-1) + i*length(interpolated_timept),2)./sum(ode_sim_data(num_rep*length(interpolated_timept)*(ini_cond_indx-1)+(i-1)*length(interpolated_timept)+1: num_rep*length(interpolated_timept)*(ini_cond_indx-1) + i*length(interpolated_timept),:),2),'b-','LineWidth',2) % plotting E fraction dynamics
                end
            end
            % plotting experimental data one replicate at a time
%             plot(timept/r1,data(num_rep*length(timept)*(ini_cond_indx-1)+(i-1)*length(timept)+1: num_rep*length(timept)*(ini_cond_indx-1) + i*length(timept),2),'ko','MarkerSize',7);
            %         plot(timept/r1,data(num_rep*length(timept)*(ini_cond_indx-1)+(i-1)*length(timept)+1: num_rep*length(timept)*(ini_cond_indx-1) + i*length(timept),3),'bo','MarkerSize',7);
        end



        if(pop_model_indx < 15)
            ax = gca;
            axis square
            ax.FontSize = 16;
            grid on
            ax.GridAlpha = 0.5;
            if(ini_cond_indx == 1)
                ax.YTick = 0:0.2:1;
            else
                ax.YTick = 0:0.1:0.2;
            end
            ax.XTick = (0:2:8)*7*24;
            ax.XTickLabels = (0:2:8);
            xlabel('time (weeks)')
            ylabel([{'EpCAM^{low} fraction'}])
        else
            ax = gca;
%             axis square
            ax.FontSize = 16;
            grid on
            ax.GridAlpha = 0.5;
            ax.XTick = (0:3:24)*24;
            ax.XTickLabels = (0:3:24);
            xlabel('time (days)')
            ax.YTick = 0:0.2:1;
            ylabel(['Cell fraction']);

        end
        %     title(['Data and System dynamics Initial Cond ' num2str(ini_cond_indx)])

        directory = pwd;
        cd('C:\Users\user\OneDrive - Indian Institute of Science\GitHub\EMT-State-Transition-2.0\Manuscript figures\Figure 3 and related SI')
        if(isempty(r2))
            if(pop_model_indx < 15)
                saveas(gcf,['Bhatia data Model Fit Ini Cond ' num2str(ini_cond_indx) ' pop mod ' num2str(pop_model_indx) ' b_t_ratio ' num2str(b_t_ratio)  ' r1 ' num2str(1/r1) ' para indx ' num2str(order_indx) '.png']);
            else
                saveas(gcf,['Cell fraction Yamamoto data Model Fit Ini Cond ' num2str(ini_cond_indx) ' pop mod ' num2str(pop_model_indx) ' b_t_ratio ' num2str(b_t_ratio)  ' r1 ' num2str(1/r1) ' para indx ' num2str(order_indx) '.png']);
            end
        else
            if(pop_model_indx < 15)
                saveas(gcf,['Bhatia data Model Fit Ini Cond ' num2str(ini_cond_indx) ' pop mod ' num2str(pop_model_indx) ' b_t_ratio ' num2str(b_t_ratio)  ' r1 ' num2str(1/r1) ' r2 ' num2str(inverse_r2)  ' para indx ' num2str(order_indx) '.png']);
            else
                saveas(gcf,['Cell fraction Yamamoto data Model Fit Ini Cond ' num2str(ini_cond_indx) ' pop mod ' num2str(pop_model_indx) ' b_t_ratio ' num2str(b_t_ratio)  ' r1 ' num2str(1/r1) ' r2 ' num2str(inverse_r2)  ' para indx ' num2str(order_indx) '.png']);
            end
        end
        cd(directory)
    end
end
end