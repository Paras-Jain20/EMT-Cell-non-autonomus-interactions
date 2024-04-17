function [param_val,chi_sq] = prof_likelihood_algo(pop_model_indx,unidentifiable_para_only,practical_identifiability_status,birth_tran_ratio,data_type ,data,num_rep,num_init_cond,timept,interpolated_timept, interpolation,lower_bound,upper_bound,profile_sample_size, noise_factor,test_data_indx, r1,r2,original_num_timepts)

unidentifiable_para_indx = find(practical_identifiability_status == 0);

param_val = cell(length(lower_bound),1);
chi_sq = cell(length(lower_bound),1); % here 1000 is the length of theta vector defined below


switch pop_model_indx % i.e. model whose parameters have dimension
    case {1,2,5,7,8,9} % these models have dimensional parameters

        % bound in units of time
        upper_bound_time = 1./lower_bound(1:4); % first four parameters (r1, r2, t1, t2) are inversed to get there units in time
        lower_bound_time = 1./upper_bound(1:4);

        if(pop_model_indx ~=2 && pop_model_indx ~=1)
            upper_bound_time = [upper_bound_time lower_bound(5:end)];
            lower_bound_time  = [lower_bound_time  upper_bound(5:end)];
        end

        for para_indx = 1:length(lower_bound) % to run the analyze for the parameters sequentially

            disp(['Running profile likelihood analysis for para indx ' num2str(para_indx)])

            theta = linspace(upper_bound_time(para_indx),lower_bound_time(para_indx),profile_sample_size);

            if(para_indx <=4)
                theta = 1./theta;
            end

            for i  = 1:size(theta,2)

                new_lower_bound = [lower_bound(1:para_indx-1) theta(i) lower_bound(para_indx+1:end)];
                new_upper_bound = [upper_bound(1:para_indx-1) theta(i) upper_bound(para_indx+1:end)];

                cost = @(params) obj_fun(params,data_type,data,pop_model_indx,num_rep,timept, interpolated_timept, interpolation);
                % note: I have not ramdomly selected the initial guess for each
                % increment of theta values. This may affect the results if there
                % are local minima
                if(i == 1)
                    initial_guess = new_lower_bound;
                else
                    initial_guess = param_val(i-1,:,para_indx);
                    initial_guess(para_indx) =  theta(i);
                end
                options = optimoptions('lsqnonlin','Display','off');
                [param_val(i,:,para_indx),chi_sq(i,para_indx)] = lsqnonlin(cost,initial_guess,new_lower_bound,new_upper_bound, options); % here, optimized para from last iteration serves as a initial guess
                %                 chi_sq(i,para_indx) = chi_sq(i,para_indx)/(num_init_cond*num_rep*length(timept)); % here, num_init_cond represents the number of initial population conditions
                chi_sq(i,para_indx) = chi_sq(i,para_indx); % here, num_init_cond represents the number of initial population conditions

            end
            %             profile_likelihood_plot(param_val(:,:,para_indx),chi_sq(:,para_indx),para_indx,pop_model_indx,birth_tran_ratio, noise_factor, data_type,test_data_indx, r1, num_rep, timept, num_init_cond)
        end

    case {3,4,6,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28}% for models whose parameters are dimensionless

        if(isempty(r2))
            if (unidentifiable_para_only == 0)
                para_indx_vals = 1:length(lower_bound);
            else
                para_indx_vals = unidentifiable_para_indx;
            end
        else
            if (unidentifiable_para_only == 0)
                para_indx_vals = 2:length(lower_bound);
            else
%                 para_indx_vals = unidentifiable_para_indx;
                para_indx_vals = unidentifiable_para_indx(unidentifiable_para_indx ~= 1);

            end
        end

        parfor para_indx = 1:max(para_indx_vals)%1:length(lower_bound) % to run the analyze for the parameters sequentially

            if(~isempty(find(para_indx == para_indx_vals,1)))
                disp(['Running profile likelihood analysis for pop model ' num2str(pop_model_indx) ' para indx ' num2str(para_indx)])

                theta = linspace(lower_bound(para_indx),upper_bound(para_indx),profile_sample_size); % no inversion of parameters is required as they are dimensionless

                for i  = 1:size(theta,2)
                                        
                    disp(['Running profile likelihood analysis for para indx ' num2str(para_indx) ' Current para iteration ' num2str(i)  ' Current para value ' num2str(theta(i)) ' Upper bound ' num2str(upper_bound(para_indx))])

                    new_lower_bound = [lower_bound(1:para_indx-1) theta(i) lower_bound(para_indx+1:end)];
                    new_upper_bound = [upper_bound(1:para_indx-1) theta(i) upper_bound(para_indx+1:end)];

                    cost = @(params) obj_fun(params, data_type,data,pop_model_indx,num_rep,timept, interpolated_timept, interpolation,original_num_timepts);
                    % note: I have not ramdomly selected the initial guess for each
                    % increment of theta values. This may affect the results if there
                    % are local minima
                    if(i == 1)
                        initial_guess = new_lower_bound;
                    else
                        initial_guess = param_val{para_indx}(i-1,:);
                        initial_guess(para_indx) =  theta(i);
                    end

                    options = optimoptions('lsqnonlin','Display','off');

                    [param_val{para_indx}(i,:),chi_sq{para_indx}(i)] = lsqnonlin(cost,initial_guess,new_lower_bound,new_upper_bound, options); % here, optimized para from last iteration serves as a initial guess
                    %                 chi_sq(i,para_indx) = chi_sq(i,para_indx)/(num_init_cond*2*length(timept)); % here, num_init_cond represents the number of initial population conditions; num_init_cond value was 2 even for Yamamoto data until date July 15, 2023
                    %                 chi_sq(i,para_indx) = chi_sq(i,para_indx); % here, num_init_cond represents the number of initial population conditions; num_init_cond value was 2 even for Yamamoto data until date July 15, 2023

                end
                %             profile_likelihood_plot(param_val(:,:,para_indx),chi_sq(:,para_indx),para_indx,pop_model_indx,birth_tran_ratio, noise_factor, data_type,test_data_indx, r1, num_rep, timept, num_init_cond)
                param_val{para_indx}(:,length(lower_bound)+1) = para_indx * ones(length(param_val{para_indx}(:,1)),1);
                disp(['Ran profile likelihood analysis for para indx ' num2str(para_indx)])

            end
        end

end
end