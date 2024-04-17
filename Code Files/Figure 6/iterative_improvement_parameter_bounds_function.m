
% profile_sample_size =1000;
% num_para_steps = 100; max number of parameter sets to be considered
% between lower and upper bound of a profile likelihood parameter
% data : experimental data; Note: code adds day zero timepoint in Bhatia el al data for caluculation of optimal sample time:
% num_rep = 1;
% num_init_cond = 5; % number of different initial condition in the data
% init_cond = [1 2] for bhatia el al. or [1 2 3 4 5] for Yamamoto;
% timept = experimental data points; this will change at each iteration on adding new timepoint to the existing exp timepoints; % time is unitless as it is scaled by m cells division rate (considered to be 1/36 hrs)
% interpolated_timept = this won't change with iteration; value should be set to 0:24*r1:timept(end);
% profile_para_data = profile data with will update
% after each iteration
% interpolation = 1;


function [updated_data, updated_data_timept, terminate ] = iterative_improvement_parameter_bounds_function(pop_model_indx, data, profile_para_data,profile_sample_size, num_para_steps, unidentifiable_para_indx, num_rep,num_init_cond,init_cond,timept,interpolated_timept,r2,original_num_timepts)

interpolation = 1;

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
        case {8,9,22,23,24} % i.e. population model with 5 parameters
            delta_chi_sq_all = 11.3;
        case {21} % i.e. population model with 6 parameters
            delta_chi_sq_all = 12.8;
    end
end

% below the std_dev_init_cond_para_time matrix stores the standard
% deviation of the ODE trajectories of profile likelihood parameters for a
% given initial condition (init_cond), parameter (para) and time point (time)
if(pop_model_indx < 15)
    std_dev_init_cond_para_time = zeros(num_init_cond, length(unidentifiable_para_indx),length(interpolated_timept)+1); % it contains the std deviation within the simulated ODE trajectories from a profile likelihood of a para;  +1 to include additional experimental time point at day 0

else
    std_dev_init_cond_para_time = zeros(num_init_cond, length(unidentifiable_para_indx),length(interpolated_timept)); % to exclude zeroth day timepoint in the Yamamoto data
end

% below the ode_sim_data_across_para_init_cond_time cell array stores the
% ODE trajectories of profile likelihood parameters for a
% given initial condition (init_cond), parameter (para) and time point (time)
ode_sim_data_across_para_init_cond_time = cell(length(unidentifiable_para_indx),1);

para_indx_seq = 1;

for para_indx = unidentifiable_para_indx

    x_rate = profile_para_data(profile_para_data(:,end) == para_indx, 1:end-2);
    y = profile_para_data(profile_para_data(:,end) == para_indx, end-1);

    [minima,posn] = min(y); %finds the minimum in y and the corresponding index

    chi_sq_threshold_all = minima + delta_chi_sq_all;%the chi_sq threshold for all params

    %finding the confidence range

    left_indx = find((y(1:posn) - chi_sq_threshold_all > 0)); % here, left indx is for rate para

    right_indx = find((y(posn:end) - chi_sq_threshold_all > 0),1); % here, right indx is for rate para

    if(~isempty(left_indx))
        left_indx = left_indx(end)+1;
    else
        left_indx = 1;
    end

    if(~isempty(right_indx))
        right_indx = posn+right_indx-2;
    else
        right_indx = length(y);
    end

    para_indx_step = ceil((right_indx-left_indx)/num_para_steps);

    % below is to skip considering those parameter in finding out optimal time point whose profile likelihood has less than three points below the chi-sqr threshold
    if(para_indx_step  >= 3)
        if(pop_model_indx < 15)
            ode_sim_data = zeros(length(interpolated_timept)+1,length(left_indx:para_indx_step:right_indx)*num_init_cond*num_rep); % +1 to include additional experimental time point at day 0
        else
            ode_sim_data = zeros(length(interpolated_timept),length(left_indx:para_indx_step:right_indx)*num_init_cond*num_rep);
        end
    else
        if(pop_model_indx < 15)
            ode_sim_data = zeros(length(interpolated_timept)+1,3*num_init_cond*num_rep); % +1 to include additional experimental time point at day 0
        else
            ode_sim_data = zeros(length(interpolated_timept),3*num_init_cond*num_rep);
        end
    end

    itr_indx = 1; % as j deosn't starts at one.
    % the condition below takes care of both at least three point below chi-sqr threshold and structural-non identifiable parameter (those whose variation doesnot effect other parameters)
    if((size(unique(x_rate(left_indx:para_indx_step:right_indx,[1:para_indx-1 para_indx+1:size(x_rate,2)]),'rows'),1)>=3))
        for j = left_indx:para_indx_step:right_indx

            temp_ode_sim_data = ODE_simulation(pop_model_indx,num_rep,timept, interpolated_timept,interpolation,data,original_num_timepts, x_rate(j,:));
            
            temp_ode_sim_data = temp_ode_sim_data(:,2:end); % removing the time column from the ode simulation data

            % concatenating the ODE solutions together for each increment in the para value for which profile likelihood is presently cosnidered
            if(pop_model_indx < 15)
                ode_sim_data(:,itr_indx:itr_indx + (num_rep*num_init_cond) -1) = [0.999 0.999 0.999 0.001 0.001 0.001; reshape(temp_ode_sim_data(:,1)./sum(temp_ode_sim_data,2),length(interpolated_timept),num_init_cond*num_rep)];
            else
                % capturing total cell numbers of the population
                ode_sim_data(:,itr_indx:itr_indx + (num_rep*num_init_cond) -1) = reshape(sum(temp_ode_sim_data(:,1:4),2),length(interpolated_timept),num_init_cond*num_rep);
               
                % % capturing M cells fraction (irrespective of Venus label) of the total population
                % ode_sim_data(:,itr_indx:itr_indx + (num_rep*num_init_cond) -1) = reshape(sum(temp_ode_sim_data(:,[1,3]),2)./sum(temp_ode_sim_data(:,1:4),2),length(interpolated_timept),num_init_cond*num_rep);
                
                % % capturing M cells numbers (irrespective of Venus label) of the total population
                % ode_sim_data(:,itr_indx:itr_indx + (num_rep*num_init_cond) -1) = reshape(sum(temp_ode_sim_data(:,[1,3]),2),length(interpolated_timept),num_init_cond*num_rep);

            end

            itr_indx = itr_indx + num_init_cond*num_rep;

        end

        ode_sim_data_across_para_init_cond_time{para_indx_seq} =  ode_sim_data; % saving the ode data for the current parameter variation
        
        % saving the std of ODE solutions at each time point
        % for a given set of initial condition and parameter
        % for which profile likelihood is considered in the
        % loop
        for init_cond_indx = 1:num_init_cond
            
            % the rep_matrix (replicate matrix) below concatenates the ODE
            % trajectories (multiple time points along the  row and three column for replicates) for different parameters in a cloumn wise fashion 
            rep_matrix = zeros(size(ode_sim_data,1),size(ode_sim_data,2)/num_init_cond);

            for rep_indx = 1:num_rep

                rep_matrix(:,1+(rep_indx-1)*(size(ode_sim_data,2)/(num_rep*num_init_cond)):rep_indx*(size(ode_sim_data,2)/(num_rep*num_init_cond))) = ode_sim_data(:,num_rep * (init_cond_indx-1) + rep_indx:num_init_cond*num_rep:end);

            end
            
            % calculating the standard deviation for a given timepoint for
            % all variation of the parameter (and across replicates)

            std_dev_init_cond_para_time(init_cond_indx, para_indx_seq,:) =  std(rep_matrix,0,2);

        end

    else
        ode_sim_data_across_para_init_cond_time{para_indx_seq} =  ode_sim_data;
    end
    para_indx_seq = para_indx_seq + 1;
end

% to check for maximum standard deviation across initial condition and parameter index

if(pop_model_indx < 15)
    % sum of the std deviation across initial
    % conditions
    std_dev_para_vs_time = reshape(sum(std_dev_init_cond_para_time,1),length(unidentifiable_para_indx),length(interpolated_timept)+1);
else
    %  std_prod_time = reshape(prod(std_dev_time,[1 2]),1,length(interpolated_timept));
    std_dev_para_vs_time = reshape(sum(std_dev_init_cond_para_time,1),length(unidentifiable_para_indx),length(interpolated_timept));
end


terminate = 0;

% here, descend_order_para_vs_time_indices stores the linear (1D) index of
% parameter and time point which gives maximum standard deviation in
% std_dev_para_vs_time matrix in the descending order

[~, descend_order_para_vs_time_indices] = sort(reshape(std_dev_para_vs_time,[],1),'descend');  

for descend_order_search = descend_order_para_vs_time_indices' % (transpose is required as descend_order_para_vs_time_indices is a column vector)

    if(rem(descend_order_search,(length(unidentifiable_para_indx))) == 0)
        opt_sample_para_indx =  length(unidentifiable_para_indx); % opt_sample_para_indx: optimal sample parameter index
        opt_sample_time_indx = floor(descend_order_search/(length(unidentifiable_para_indx))); % opt_sample_time_indx: optimal sample time index
    else
        opt_sample_para_indx = rem(descend_order_search,(length(unidentifiable_para_indx)));
        opt_sample_time_indx = floor(descend_order_search/(length(unidentifiable_para_indx))) + 1;
    end

    %the lines of codes till now gets us the opt_para_indx and the opt time
    %indx

    if(pop_model_indx < 15)
        temp_interpolated_timept = [0 interpolated_timept];
        opt_sample_time = temp_interpolated_timept(opt_sample_time_indx); % optimal sampling time
    else
        opt_sample_time = interpolated_timept(opt_sample_time_indx);
    end

    % here since we are measuring another output variable, it could be
    % measured at the already given experimental timepts; that's why
    % commenting the lines below
    if(length(timept) == original_num_timepts) 
        if(opt_sample_time  == 0)% to exclude the sampling of 0th time point
            if(descend_order_search == descend_order_para_vs_time_indices(end))
                terminate = 1;
            end
            continue;
        else
            break
        end
    else
        if(opt_sample_time  == 0 || (sum(timept(original_num_timepts+1:end) == opt_sample_time) ~= 0)) % to exclude the sampling of 0th time point as well as a point time point that is already sampled before
            if(descend_order_search == descend_order_para_vs_time_indices(end))
                terminate = 1;
            end
            continue;
        else
            break
        end
    end

end

% sampling new data at the new time point and parameter index
ode_data = ode_sim_data_across_para_init_cond_time{opt_sample_para_indx}(opt_sample_time_indx,:);
% initial condition wise ODE data
ode_data = reshape(ode_data,num_init_cond*num_rep,length(ode_data)/(num_init_cond*num_rep));

ode_data_init_cond_wise = zeros(num_init_cond,size(ode_data,2)*num_rep);
for i = 1:num_init_cond
    ode_data_init_cond_wise(i,:) = reshape(ode_data((i-1)*num_rep+1:i*num_rep,:),[],1);
end

total_perm = 100;
sampled_data = zeros(num_init_cond,num_rep,total_perm);
sampled_data_std = zeros(num_init_cond,total_perm);

for permutation_indx = 1:total_perm
    permutation = randperm(size(ode_data_init_cond_wise,2));
    sampled_data(:,:,permutation_indx) = ode_data_init_cond_wise(:,permutation(1:num_rep));
    sampled_data_std(:,permutation_indx) = std(sampled_data(:,:,permutation_indx),0,2);
end

reshaped_sampled_data_std = reshape(prod(sampled_data_std,1),[],total_perm); 
[~,max_std_dev_permutation_indx] = max(reshaped_sampled_data_std,[],'all');

if terminate == 1
    updated_data = data;
    updated_data_timept = timept;
else

updated_data = zeros(size(data,1) + num_rep*num_init_cond,size(data,2));

%%% following lines of code concatenates the new data of another output
%%% below the experimental data of Venus labelled EpCAM low cell fraction
updated_data(1:size(data(1:original_num_timepts*num_rep*num_init_cond,:),1),:) = data(1:original_num_timepts*num_rep*num_init_cond,:);

if(length(timept) == original_num_timepts)
    interpolation_indx = 1;
    updated_data_timept = [timept(1:original_num_timepts) opt_sample_time];

else
    uni_exist_timept = unique(timept(original_num_timepts+1:end)); % unique exisiting time points
    interpolation_indx = find(sort([uni_exist_timept opt_sample_time]) == opt_sample_time,1);

    updated_data_timept = [timept(1:original_num_timepts) sort([uni_exist_timept opt_sample_time])];
end

for i = 1:num_init_cond
    for j = 1:num_rep
        updated_data(original_num_timepts*num_rep*num_init_cond + 1 + (j-1)*(length(timept)-original_num_timepts+1) + (i-1)*num_rep*(length(timept)-original_num_timepts+1): original_num_timepts*num_rep*num_init_cond + (length(timept)-original_num_timepts+1) +(j-1)*(length(timept)-original_num_timepts+1) + (i-1)*num_rep*(length(timept)-original_num_timepts+1),:) =  ...
            [data(original_num_timepts*num_rep*num_init_cond + 1 + (j-1)*(length(timept)-original_num_timepts) + (i-1)*num_rep*(length(timept)-original_num_timepts):original_num_timepts*num_rep*num_init_cond + interpolation_indx-1 + (j-1)*(length(timept)-original_num_timepts) + (i-1)*num_rep*(length(timept)-original_num_timepts),:);...
            opt_sample_time sampled_data(i,j,max_std_dev_permutation_indx)  1-sampled_data(i,j,max_std_dev_permutation_indx)  sampled_data_std(i,max_std_dev_permutation_indx) j i 2;...
            data(original_num_timepts*num_rep*num_init_cond + interpolation_indx + (j-1)*(length(timept)-original_num_timepts) + (i-1)*num_rep*(length(timept)-original_num_timepts):original_num_timepts*num_rep*num_init_cond + (length(timept)-original_num_timepts) +(j-1)*(length(timept)-original_num_timepts) + (i-1)*num_rep*(length(timept)-original_num_timepts),:)];
    end
end

%%% following lines of code adds the new data of same output
%%% as the experimental data of Venus labelled EpCAM low cell
%%% fraction according to the ascending timepts

% for i = 1:num_init_cond
%     for j = 1:num_rep
%         interpolation_indx = find(data(1+(j-1)*length(timept) + (i-1)*num_rep*length(timept): length(timept) +(j-1)*length(timept) + (i-1)*num_rep*length(timept),1) >= opt_sample_time,1);
%         updated_data(1+(j-1)*(length(timept)+1) + (i-1)*num_rep*(length(timept)+1): (length(timept)+1) +(j-1)*(length(timept)+1) + (i-1)*num_rep*(length(timept)+1),:) =  ...
%             [data(1+(j-1)*length(timept) + (i-1)*num_rep*length(timept):interpolation_indx-1 + (j-1)*length(timept) + (i-1)*num_rep*length(timept),:); opt_sample_time sampled_data(i,j,max_std_indx)  1-sampled_data(i,j,max_std_indx)  sampled_data_std(i,max_std_indx); data(interpolation_indx + (j-1)*length(timept) + (i-1)*num_rep*length(timept):length(timept) +(j-1)*length(timept) + (i-1)*num_rep*length(timept),:)];
%     end
% end
%     interpolation_indx = find(timept > opt_sample_time,1);
%     updated_data_timept = [timept(1:interpolation_indx-1) opt_sample_time timept(interpolation_indx:end)];

disp('Done with adding data at new time point')
end


end