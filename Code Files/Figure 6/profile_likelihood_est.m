

function [params, chi_sqr_err] = profile_likelihood_est(pop_model_indx,data_type, data, num_rep,num_init_cond,timept,interpolated_timept, interpolation,birth_tran_ratio,profile_sample_size, noise_factor, test_data_indx, r1,unidentifiable_para_only,practical_identifiability_status,r2,original_num_timepts)
% rng('default')

if(isempty(r2))
switch pop_model_indx

    case {1,2}
        lower_bound = [1/60 1/60 1/(60*birth_tran_ratio) 1/(60*birth_tran_ratio)];
        upper_bound = [1/20 1/20 1/20 1/20];

    case {3,18}
        lower_bound = [1 1/(birth_tran_ratio) 1/(birth_tran_ratio)];
        upper_bound = [3 1 1];

    case {4}

        lower_bound = [1 1/(birth_tran_ratio)  1/(birth_tran_ratio)  -1 -1];
        upper_bound = [3 1 1 1 1];

    case {5}

        lower_bound = [1/60 1/60 1/(60*birth_tran_ratio)  1/(60*birth_tran_ratio)  0];
        upper_bound = [1/20 1/20 1/20 1/20 4];

    case {6}

        lower_bound = [1];
        upper_bound = [5];
    case {7}

        lower_bound = [1 0 0];
        upper_bound = [3 1 1];

    case {8,9}

        lower_bound = [1/60 1/60 1/(60*birth_tran_ratio)  1/(60*birth_tran_ratio)  0 0];
        upper_bound = [1/20 1/20 1/20 1/20 1 4];
    case {10}

        lower_bound = [1 1/(birth_tran_ratio)  1/(birth_tran_ratio)  0];
        upper_bound = [3 1 1 4];

    case {11,12}

        lower_bound = [1 1/(birth_tran_ratio)  1/(birth_tran_ratio)  0];
        upper_bound = [3 1 1 1];

    case {13}

        lower_bound = [1 1/(birth_tran_ratio)  1/(birth_tran_ratio)  0 0];
        upper_bound = [3 1 1 1 4];

    case {14}

        lower_bound = [1 1/(birth_tran_ratio)  1/(birth_tran_ratio)  0 0];
        upper_bound = [3 1 1 1 1];

    case 15
        lower_bound = [0 1/(birth_tran_ratio)  1/(birth_tran_ratio)  0 0];
        upper_bound = [3 1 1 10 10];

    case {16,17}
        lower_bound = [0 1/(birth_tran_ratio)  1/(birth_tran_ratio)  0 0];
        upper_bound = [3 1 1 1 10];

    case {19}
        lower_bound = [1 1000000];
        upper_bound = [3 20000000];
    case {20}
        lower_bound = [1 1/birth_tran_ratio 1/birth_tran_ratio 1000000];
        upper_bound = [3 1 1 20000000];
    case {21}
        lower_bound = [1 1/birth_tran_ratio 1/birth_tran_ratio 0 0  1000000 1000000 ];
        upper_bound = [3 1 1 10 10 20000000 20000000];
    case {22}
        lower_bound = [1 1/birth_tran_ratio 1/birth_tran_ratio 1000000 0 0];
        upper_bound = [3 1 1 20000000 10 10];
    case {23,24}
        lower_bound = [1 1/birth_tran_ratio 1/birth_tran_ratio 1000000 0 0];
        upper_bound = [3 1 1 20000000 1 10];
    case {25}
        lower_bound = [1 1/birth_tran_ratio 1/birth_tran_ratio 1000000 0 0];
        upper_bound = [3 1 1 20000000 1 1];
    case {26,27}
        lower_bound = [1 1/birth_tran_ratio 1/birth_tran_ratio 1000000 0];
        upper_bound = [3 1 1 20000000 1];
    case {28}
        lower_bound = [0 0 0 0 0 0 0 0 0 0 0 0 0];
        upper_bound = [1 1 1 1 1 1 1 1 1 1 1 1 1];
    case {29}

        lower_bound = [1 0 0  1000000 1000000 ];
        upper_bound = [3 10 10 20000000 20000000];

end
else
switch pop_model_indx

    case {1,2}
        lower_bound = [1/60 1/60 1/(60*birth_tran_ratio) 1/(60*birth_tran_ratio)];
        upper_bound = [1/20 1/20 1/20 1/20];

    case {3,18}
        lower_bound = [r2 1/(birth_tran_ratio) 1/(birth_tran_ratio)];
        upper_bound = [r2 1 1];

    case {4}

        lower_bound = [r2 1/(birth_tran_ratio)  1/(birth_tran_ratio)  -1 -1];
        upper_bound = [r2 1 1 1 1];

    case {5}

        lower_bound = [1/60 1/60 1/(60*birth_tran_ratio)  1/(60*birth_tran_ratio)  0];
        upper_bound = [1/20 1/20 1/20 1/20 4];

    case {6}

        lower_bound = [1];
        upper_bound = [5];
    case {7}

        lower_bound = [1 0 0];
        upper_bound = [3 1 1];

    case {8,9}

        lower_bound = [1/60 1/60 1/(60*birth_tran_ratio)  1/(60*birth_tran_ratio)  0 0];
        upper_bound = [1/20 1/20 1/20 1/20 1 4];
    case {10}

        lower_bound = [r2 1/(birth_tran_ratio)  1/(birth_tran_ratio)  0];
        upper_bound = [r2 1 1 4];

    case {11,12}

        lower_bound = [r2 1/(birth_tran_ratio)  1/(birth_tran_ratio)  0];
        upper_bound = [r2 1 1 1];

    case {13}

        lower_bound = [r2 1/(birth_tran_ratio)  1/(birth_tran_ratio)  0 0];
        upper_bound = [r2 1 1 1 4];

    case {14}

        lower_bound = [r2 1/(birth_tran_ratio)  1/(birth_tran_ratio)  0 0];
        upper_bound = [r2 1 1 1 1];

    case 15
        lower_bound = [0 1/(birth_tran_ratio)  1/(birth_tran_ratio)  0 0];
        upper_bound = [3 1 1 10 10];

    case {16,17}
        lower_bound = [0 1/(birth_tran_ratio)  1/(birth_tran_ratio)  0 0];
        upper_bound = [3 1 1 1 10];

    case {19}
        lower_bound = [r2 1000000];
        upper_bound = [r2 20000000];
    case {20}
        lower_bound = [r2 1/birth_tran_ratio 1/birth_tran_ratio 1000000];
        upper_bound = [r2 1 1 20000000];
    case {21}
        lower_bound = [r2 1/birth_tran_ratio 1/birth_tran_ratio 0 0  100000 1000000 ];
        upper_bound = [r2 1 1 10 10 20000000 20000000];
    case {22}
        lower_bound = [r2 1/birth_tran_ratio 1/birth_tran_ratio 1000000 0 0];
        upper_bound = [r2 1 1 20000000 10 10];
    case {23,24}
        lower_bound = [r2 1/birth_tran_ratio 1/birth_tran_ratio 1000000 0 0];
        upper_bound = [r2 1 1 20000000 1 10];
    case {25}
        lower_bound = [r2 1/birth_tran_ratio 1/birth_tran_ratio 1000000 0 0];
        upper_bound = [r2 1 1 20000000 1 1];
    case {26,27}
        lower_bound = [r2 1/birth_tran_ratio 1/birth_tran_ratio 1000000 0];
        upper_bound = [r2 1 1 20000000 1];
    case {28}
        lower_bound = [0 0 0 0 0 0 0 0 0 0 0 0 0];
        upper_bound = [1 1 1 1 1 1 1 1 1 1 1 1 1];
    case {29}

        lower_bound = [r2 0 0  1000000 1000000 ];
        upper_bound = [r2 10 10 20000000 20000000];

end
end

% [params, chi_sqr_err] = prof_likelihood_algo_modified(pop_model_indx,unidentifiable_para_only,practical_identifiability_status,birth_tran_ratio, data_type,data,num_rep,num_init_cond,timept,interpolated_timept, interpolation,lower_bound,upper_bound,profile_sample_size, noise_factor, test_data_indx, r1,r2);
[params, chi_sqr_err] = prof_likelihood_algo(pop_model_indx,unidentifiable_para_only,practical_identifiability_status,birth_tran_ratio, data_type,data,num_rep,num_init_cond,timept,interpolated_timept, interpolation,lower_bound,upper_bound,profile_sample_size, noise_factor, test_data_indx, r1,r2,original_num_timepts);

end