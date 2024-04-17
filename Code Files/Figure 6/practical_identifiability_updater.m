%function for updating each iteration of the sampling procedure

function [updated_data,updated_timept,updated_profile_liklihood,updated_practical_identifiability_status,terminate] = practical_identifiability_updater(pop_model_indx,practical_identifiability_status,current_profile_data,profile_sample_size, num_para_steps,current_data,num_rep,num_init_cond,current_timept,interpolated_timept,r1,r2,birth_tran_ratio,original_num_timepts)

unidentifiable_para_indx = find(practical_identifiability_status == 0);%gets the indices of the unidentifiable params

% getting the optimum timept to sample data
[updated_data, updated_timept,terminate] = iterative_improvement_parameter_bounds_function(pop_model_indx,current_data, current_profile_data,profile_sample_size, num_para_steps, unidentifiable_para_indx,num_rep,num_init_cond,1:num_init_cond,current_timept,interpolated_timept,r2,original_num_timepts);

if terminate == 0
    %sampled_data_profile_data_store = {};%stores the profile liklihood generated from data post sampling for each unidentifiable parameter
    
    %generating the profile liklihood from the sampled data
    %some input variables of profile likelihood function that are
    %constant as far as the whole sampling analysis is concerned
    data_type = 2;
    interpolation = 0;
    noise_factor = 1;
    test_data_indx = 1;
    
    unidentifiable_para_only = 1;% 1 for generating profile liklihood of unidentifiable params only and 0 for generating prof liklihood of all params
    
    [params, chi_sqr_err] = profile_likelihood_est(pop_model_indx,data_type, updated_data, num_rep,num_init_cond,updated_timept, interpolated_timept, interpolation,birth_tran_ratio,profile_sample_size, noise_factor, test_data_indx, r1,unidentifiable_para_only,practical_identifiability_status,r2,original_num_timepts);
    % saving the parameters and chi-sqr values obtained from the
    % analysis in a single matrix
    para_data = [];
    
    for i = 1:length(params) % since params is the cell vector
        if(~isempty(params{i}))
            para_data =  [para_data; params{i}(:,1:end-1) chi_sqr_err{i}' params{i}(:,end)];
        end
    end
    
    updated_profile_liklihood = para_data;
    %checking the practical identifiability of the profile likelihood
    %generated from the data post sampling
    updated_practical_identifiability_status = practical_identifiability_checker_modified(pop_model_indx,para_data,practical_identifiability_status,r2);
    %sampled_data_profile_data_store{end+1} = para_data;
    
    
    %finding the index of the maximum practical identifiability score
    %[maximum_score,max_score_indx] = max(practical_identifiability_score_store);
    
    % updated_data = sampled_data_store{max_score_indx};
    % updated_timept = updated_timept_store{max_score_indx};
    % updated_profile_liklihood = sampled_data_profile_data_store{max_score_indx};
    % updated_practical_identifiability_status_score = maximum_score;

else
    updated_profile_liklihood = current_profile_data;
    updated_practical_identifiability_status = practical_identifiability_status;
    
end
end
