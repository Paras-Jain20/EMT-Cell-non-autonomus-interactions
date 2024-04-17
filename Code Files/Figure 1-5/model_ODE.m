function out = model_ODE

out{1} = @exp_grow_with_trans_cell_numbers; % cell number dynamics
out{2} = @exp_grow_with_trans_cell_fractions; % cell fraction dynamics
out{3} = @exp_grow_with_trans_cell_fractions_non_dim; % cell fraction dynamics non dim
out{4} = @exp_grow_competition_with_trans_cell_fraction_non_dim; % cell fraction dynamics non dim
out{5} = @exp_grow_cell_fraction;
out{6} = @exp_grow_cell_fraction_non_dim;
out{7} = @exp_grow_competition_cell_fraction_non_dim;
out{8} = @exp_grow_trans_ret_of_M_inf_of_M_on_E_cell_fractions;
out{9} = @exp_grow_trans_ret_of_E_inf_of_M_on_E_cell_fractions;
out{10} = @exp_grow_trans_influence_of_M_on_E_cell_fractions_non_dim;
out{11} = @exp_grow_trans_retention_of_M_cell_fractions_non_dim;
out{12} = @exp_grow_trans_retention_of_E_cell_fractions_non_dim;
out{13} = @exp_grow_trans_ret_of_M_inf_of_M_on_E_cell_fractions_non_dim;
out{14} = @exp_grow_trans_ret_of_E_ret_of_M_cell_fractions_non_dim;
out{15} = @exp_grow_trans_ret_of_E_inf_of_M_on_E_cell_fractions_non_dim;

% models for Yamamoto et al. data
out{16} = @exp_grow_trans_influence_of_M_reten_of_E_cell_frac_ND_VL_NL; % ND: non dimensional; VL: venus labelled; NL: Venus non labelled cells
out{17} = @exp_grow_trans_influence_of_M_reten_of_M_cell_frac_ND_VL_NL;
out{18} = @exp_grow_trans_ND_VL_NL;
out{19} = @logistic_grow_ND_VL_NL;
out{20} = @logistic_grow_trans_ND_VL_NL;
out{21} = @logistic_grow_comp_trans_ND_VL_NL;
out{22} = @logistic_grow_trans_influence_of_both_M_E_ND_VL_NL;
out{23} = @logistic_grow_trans_influence_of_M_reten_of_E_ND_VL_NL;
out{24} = @logistic_grow_trans_influence_of_M_reten_of_M_ND_VL_NL;
out{25} = @logistic_grow_trans_reten_of_E_reten_of_M_ND_VL_NL;
out{26} = @logistic_grow_trans_reten_of_E_ND_VL_NL;
out{27} = @logistic_grow_trans_reten_of_M_ND_VL_NL;
out{28} = @three_state_transition_model_non_dim_H_M_combined;
out{29} = @logistic_grow_comp_ND_VL_NL;
out{30} = @logistic_grow_comp_trans_independent_parameters_ND_VL_NL;
out{31} = @logistic_grow_comp_trans_asymmetric_growth_ND_VL_NL;
out{32} = @logistic_grow_influence_ND_VL_NL;
out{33} = @logistic_grow_influence_trans_ND_VL_NL;
out{34} = @logistic_grow_trans_VL_NL;
out{35} = @logistic_grow_comp_trans_VL_NL;
out{36} = @logistic_grow_trans_influence_of_M_reten_of_E_VL_NL; 
out{37} = @logistic_grow_trans_reten_of_E_reten_of_M_VL_NL;
out{38} = @logistic_grow_trans_reten_of_E_VL_NL;
out{39} = @logistic_grow_trans_reten_of_M_VL_NL;
out{40} = @logistic_grow_influence_trans_VL_NL;
out{41} = @logistic_grow_VL_NL;
out{42} = @logistic_grow_comp_VL_NL;
out{43} = @logistic_grow_influence_VL_NL;
out{44} = @logistic_grow_trans_influence_of_M_reten_of_M_VL_NL;

out{45} = @logistic_grow_comp_trans_reten_of_E_ND_VL_NL;
out{46} = @logistic_grow_comp_trans_reten_of_M_ND_VL_NL;
out{47} = @logistic_grow_comp_trans_influence_of_M_ND_VL_NL;
out{48} = @logistic_grow_influence_trans_retention_of_E_ND_VL_NL;
out{49} = @logistic_grow_influence_trans_retention_of_M_ND_VL_NL;
out{50} = @logistic_grow_influence_trans_influence_of_M_ND_VL_NL;

out{51} = @three_state_transition_model_non_dim_E_H_combined;

out{52} = @linear_three_state_transition_model_non_dim_H_M_combined;
out{53} = @linear_three_state_transition_model_non_dim_E_H_combined;

out{54} = @linear_three_state_logistic_grow_trans_ND_VL_NL_E_H_combined;
out{55} = @linear_three_state_logistic_grow_trans_ND_VL_NL_H_M_combined;
out{56} = @three_state_logistic_grow_trans_ND_VL_NL_E_H_combined;
out{57} = @three_state_logistic_grow_trans_ND_VL_NL_H_M_combined;
end


function dydt = exp_grow_with_trans_cell_numbers(t,y,para)
dydt(1,1) = y(1)*(para(1)-para(3))+y(2)*para(4);
dydt(2,1) = y(2)*(para(2)-para(4))+y(1)*para(3);
end


function dydt = exp_grow_competition_with_trans_cell_fraction_non_dim(t,y,para)
    dydt(1,1) = y(1)*(1-para(4)*y(2))  - para(2)*y(1) + para(3)*y(2) - (y(1)*(1-para(4)*y(2)) + para(1)*y(2)*(1-para(5)*y(1)))* y(1);
    dydt(2,1) = para(1)*y(2)*(1-para(5)*y(1))  + para(2)*y(1) - para(3)*y(2) - (y(1)*(1-para(4)*y(2)) + para(1)*y(2)*(1-para(5)*y(1)))* y(2);
end

function dydt = exp_grow_competition_cell_fraction_non_dim(t,y,para)
    dydt(1,1) = y(1)*(1-para(2)*y(2))  - (y(1)*(1-para(2)*y(2)) + para(1)*y(2)*(1-para(3)*y(1)))* y(1);
    dydt(2,1) = para(1)*y(2)*(1-para(3)*y(1))  - (y(1)*(1-para(2)*y(2)) + para(1)*y(2)*(1-para(3)*y(1)))* y(2);
end

function dydt = exp_grow_with_trans_cell_fractions(t,y,para)
dydt(1,1) = y(1)*(para(1)-para(3))+y(2)*(para(4)-para(2)*y(1)) - (y(1)^2)*para(1);
dydt(2,1) = y(2)*(para(2)-para(4))+y(1)*(para(3)-para(1)*y(2)) - (y(2)^2)*para(2);
end

function dydt = exp_grow_with_trans_cell_fractions_non_dim(t,y,para)
dydt(1,1) = y(1)*(1-para(2))+y(2)*(para(3)-para(1)*y(1)) - (y(1)^2);
dydt(2,1) = y(2)*(para(1)-para(3))+y(1)*(para(2)-y(2)) - (y(2)^2)*para(1);
end

function dydt = exp_grow_cell_fraction(t,y,para)
    dydt(1,1) = para(1)*y(1) - (para(1)*y(1) + para(2)*y(2))* y(1);
    dydt(2,1) = y(2)*para(2)  - (para(1)*y(1) + para(2)*y(2))* y(2);
end

function dydt = exp_grow_cell_fraction_non_dim(t,y,para)
    dydt(1,1) = y(1) - (y(1) + para(1)*y(2))* y(1);
    dydt(2,1) = y(2)*para(1)  - (y(1) + para(1)*y(2))* y(2);
end


function dydt = exp_grow_trans_retention_of_E_cell_fractions(t,y,para)
    dydt(1,1) = y(1)*para(1)  -para(3)*y(1) +y(2)*para(4)*(1- para(5)*y(2)/sum(y)) - (para(1)*y(1) + para(2)*y(2))* y(1);
    dydt(2,1) = y(2)*para(2)  -para(4)*(1- para(5)*y(2)/sum(y))*y(2) +y(1)*para(3) - (para(1)*y(1) + para(2)*y(2))* y(2);
end

function dydt = exp_grow_trans_influence_of_M_on_E_cell_fractions(t,y,para)
    dydt(1,1) = y(1)*para(1)  -para(3)*y(1) + y(2)*para(4)*(1+ para(5)*y(1)/sum(y)) - (para(1)*y(1) + para(2)*y(2))* y(1);
    dydt(2,1) = y(2)*para(2)  -para(4)*(1+ para(5)*y(1)/sum(y))*y(2) +y(1)*para(3) - (para(1)*y(1) + para(2)*y(2))* y(2);
end

function dydt = exp_grow_trans_ret_of_M_inf_of_M_on_E_cell_fractions(t,y,para)
    dydt(1,1) = y(1)*para(1)  -para(3)*y(1)*(1-para(5)*y(1)/sum(y)) + y(2)*para(4)*(1+ para(6)*y(1)/sum(y)) - (para(1)*y(1) + para(2)*y(2))* y(1);
    dydt(2,1) = y(2)*para(2)  -para(4)*(1+ para(6)*y(1)/sum(y))*y(2) + y(1)*para(3)*(1-para(5)*y(1)/sum(y)) - (para(1)*y(1) + para(2)*y(2))* y(2);
end

function dydt = exp_grow_trans_ret_of_E_inf_of_M_on_E_cell_fractions(t,y,para)
    dydt(1,1) = y(1)*para(1)  -para(3)*y(1) + y(2)*para(4)*(1+ para(6)*y(1)/sum(y) - para(5)*y(2)/sum(y)) - (para(1)*y(1) + para(2)*y(2))* y(1);
    dydt(2,1) = y(2)*para(2)  -para(4)*(1+ para(6)*y(1)/sum(y) - para(5)*y(2)/sum(y))*y(2) + y(1)*para(3) - (para(1)*y(1) + para(2)*y(2))* y(2);
end

function dydt = exp_grow_trans_influence_of_M_on_E_cell_fractions_non_dim(t,y,para)
    dydt(1,1) = y(1)  -para(2)*y(1) + y(2)*para(3)*(1+ para(4)*y(1)) - (y(1) + para(1)*y(2))* y(1);
    dydt(2,1) = y(2)*para(1)  -para(3)*(1+ para(4)*y(1))*y(2) +y(1)*para(2) - (y(1) + para(1)*y(2))* y(2);
end

function dydt = exp_grow_trans_retention_of_M_cell_fractions_non_dim(t,y,para)
    dydt(1,1) = y(1)  - para(2)*(1- para(4)*y(1))*y(1) +y(2)*para(3) - (y(1) + para(1)*y(2))* y(1);
    dydt(2,1) = y(2)*para(1)  - para(3)*y(2) +y(1)*para(2)*(1- para(4)*y(1)) - (y(1) + para(1)*y(2))* y(2);
end

function dydt = exp_grow_trans_retention_of_E_cell_fractions_non_dim(t,y,para)
    dydt(1,1) = y(1)  -para(2)*y(1) +y(2)*para(3)*(1- para(4)*y(2)) - (y(1) + para(1)*y(2))* y(1);
    dydt(2,1) = y(2)*para(1)  -para(3)*(1- para(4)*y(2))*y(2) +y(1)*para(2) - (y(1) + para(1)*y(2))* y(2);
end

function dydt = exp_grow_trans_ret_of_M_inf_of_M_on_E_cell_fractions_non_dim(t,y,para)
    dydt(1,1) = y(1)  -para(2)*y(1)*(1-para(4)*y(1)) + y(2)*para(3)*(1+ para(5)*y(1)) - (y(1) + para(1)*y(2))* y(1);
    dydt(2,1) = y(2)*para(1)  -para(3)*(1+ para(5)*y(1))*y(2) + y(1)*para(2)*(1-para(4)*y(1)) - (y(1) + para(1)*y(2))* y(2);
end

function dydt = exp_grow_trans_ret_of_E_ret_of_M_cell_fractions_non_dim(t,y,para)
    dydt(1,1) = y(1)  -para(2)*(1-para(5)*y(1))*y(1) + y(2)*para(3)*(1 - para(4)*y(2)) - (y(1) + para(1)*y(2))* y(1);
    dydt(2,1) = y(2)*para(1)  -para(3)*(1 - para(4)*y(2))*y(2) + y(1)*(1-para(5)*y(1))*para(2) - (y(1) + para(1)*y(2))* y(2);
end

function dydt = exp_grow_trans_ret_of_E_inf_of_M_on_E_cell_fractions_non_dim(t,y,para)
    dydt(1,1) = y(1)  -para(2)*y(1) + y(2)*para(3)*(1 - para(4)*y(2) + para(5)*y(1)) - (y(1) + para(1)*y(2))* y(1);
    dydt(2,1) = y(2)*para(1)  -para(3)*(1 - para(4)*y(2) + para(5)*y(1))*y(2) + y(1)*para(2) - (y(1) + para(1)*y(2))* y(2);
end

function dydt = exp_grow_trans_influence_of_both_M_E_cell_frac_ND_VL_NL(t,y,para)    
    dydt(1,1) = y(1)  -para(2)*(1+para(4)*(y(2)+y(4)))*y(1) + y(2)*para(3)*(1+ para(5)*(y(1)+y(3))) - (y(1) + y(3) + para(1)*(y(2)+y(4)))* y(1);
    dydt(2,1) = y(2)*para(1)  - y(2)*para(3)*(1+ para(5)*(y(1)+y(3))) + para(2)*(1+para(4)*(y(2)+y(4)))*y(1) - (y(1) + y(3) + para(1)*(y(2)+y(4)))* y(2);
    dydt(3,1) = y(3)  -para(2)*(1+para(4)*(y(2)+y(4)))*y(3) + y(4)*para(3)*(1+ para(5)*(y(1)+y(3))) - (y(1) + y(3) + para(1)*(y(2)+y(4)))* y(3);
    dydt(4,1) = y(4)*para(1)  - y(4)*para(3)*(1+ para(5)*(y(1)+y(3))) + para(2)*(1+para(4)*(y(2)+y(4)))*y(3) - (y(1) + y(3) + para(1)*(y(2)+y(4)))* y(4);
end

function dydt = exp_grow_trans_influence_of_M_reten_of_E_cell_frac_ND_VL_NL(t,y,para)
    % para(4) - retention parameter; para(5): influence parameter
    dydt(1,1) = y(1)  -para(2)*y(1) + y(2)*para(3)*(1+ para(5)*(y(1)+y(3)) - para(4)*(y(2)+y(4))) - (y(1) + y(3) + para(1)*(y(2)+y(4)))* y(1);
    dydt(2,1) = y(2)*para(1)  - y(2)*para(3)*(1+ para(5)*(y(1)+y(3)) - para(4)*(y(2)+y(4))) + para(2)*y(1) - (y(1) + y(3) + para(1)*(y(2)+y(4)))* y(2);
    dydt(3,1) = y(3)  -para(2)*y(3) + y(4)*para(3)*(1+ para(5)*(y(1)+y(3)) - para(4)*(y(2)+y(4))) - (y(1) + y(3) + para(1)*(y(2)+y(4)))* y(3);
    dydt(4,1) = y(4)*para(1)  - y(4)*para(3)*(1+ para(5)*(y(1)+y(3)) - para(4)*(y(2)+y(4))) + para(2)*y(3) - (y(1) + y(3) + para(1)*(y(2)+y(4)))* y(4);
end

function dydt = exp_grow_trans_influence_of_M_reten_of_M_cell_frac_ND_VL_NL(t,y,para)  
    % para(4) - retention parameter; para(5): influence parameter
    dydt(1,1) = y(1)  -para(2)*(1 - para(4)*(y(1)+y(3)))*y(1) + y(2)*para(3)*(1+ para(5)*(y(1)+y(3))) - (y(1) + y(3) + para(1)*(y(2)+y(4)))* y(1);
    dydt(2,1) = y(2)*para(1)  - y(2)*para(3)*(1+ para(5)*(y(1)+y(3))) + para(2)*(1 - para(4)*(y(1)+y(3)))*y(1) - (y(1) + y(3) + para(1)*(y(2)+y(4)))* y(2);
    dydt(3,1) = y(3)  -para(2)*(1 - para(4)*(y(1)+y(3)))*y(3) + y(4)*para(3)*(1+ para(5)*(y(1)+y(3))) - (y(1) + y(3) + para(1)*(y(2)+y(4)))* y(3);
    dydt(4,1) = y(4)*para(1)  - y(4)*para(3)*(1+ para(5)*(y(1)+y(3))) + para(2)*(1 - para(4)*(y(1)+y(3)))*y(3) - (y(1) + y(3) + para(1)*(y(2)+y(4)))* y(4);
end

function dydt = exp_grow_trans_ND_VL_NL(t,y,para)    
    dydt(1,1) = y(1)  -para(2)*y(1) + y(2)*para(3) - (y(1) + y(3) + para(1)*(y(2)+y(4)))* y(1);
    dydt(2,1) = y(2)*para(1)  - y(2)*para(3) + para(2)*y(1) - (y(1) + y(3) + para(1)*(y(2)+y(4)))* y(2);
    dydt(3,1) = y(3)  -para(2)*y(3) + y(4)*para(3) - (y(1) + y(3) + para(1)*(y(2)+y(4)))* y(3);
    dydt(4,1) = y(4)*para(1)  - y(4)*para(3) + para(2)*y(3) - (y(1) + y(3) + para(1)*(y(2)+y(4)))* y(4);
end

function dydt = logistic_grow_ND_VL_NL(t,y,para)
    
    dydt(1,1) = y(1)*(1-((y(1)+y(3))+(y(2)+y(4)))/para(2));
    dydt(2,1) = para(1)*y(2)*(1-((y(2)+y(4))+(y(1)+y(3)))/para(2));
    dydt(3,1) = y(3)*(1-((y(1)+y(3))+(y(2)+y(4)))/para(2));
    dydt(4,1) = para(1)*y(4)*(1-((y(2)+y(4))+(y(1)+y(3)))/para(2));
    
end

function dydt = logistic_grow_trans_ND_VL_NL(t,y,para)
    
    dydt(1,1) = y(1)*(1-sum(y)/para(4)) - para(2)*y(1) + para(3)*y(2);
    dydt(2,1) = para(1)*y(2)*(1-sum(y)/para(4)) + para(2)*y(1) - para(3)*y(2);
    dydt(3,1) = y(3)*(1-sum(y)/para(4)) - para(2)*y(3) + para(3)*y(4);
    dydt(4,1) = para(1)*y(4)*(1-sum(y)/para(4)) + para(2)*y(3) - para(3)*y(4);
    
end

function dydt = logistic_grow_comp_trans_ND_VL_NL(t,y,para)
    
    dydt(1,1) = y(1)*(1-((y(1)+y(3))+(y(2)+y(4))*para(4))/para(6)) - para(2)*y(1) + para(3)*y(2);
    dydt(2,1) = para(1)*y(2)*(1-((y(2)+y(4))+(y(1)+y(3))*para(5))/para(7)) + para(2)*y(1) - para(3)*y(2);
    dydt(3,1) = y(3)*(1-((y(1)+y(3))+(y(2)+y(4))*para(4))/para(6)) - para(2)*y(3) + para(3)*y(4);
    dydt(4,1) = para(1)*y(4)*(1-((y(2)+y(4))+(y(1)+y(3))*para(5))/para(7)) + para(2)*y(3) - para(3)*y(4);
    
end

function dydt = logistic_grow_trans_influence_of_both_M_E_ND_VL_NL(t,y,para)
    
    dydt(1,1) = y(1)*(1-sum(y)/para(4)) - para(2)*(1+para(5)*(y(2)+y(4))/sum(y))*y(1) + para(3)*(1+para(6)*(y(1)+y(3))/sum(y))*y(2);
    dydt(2,1) = para(1)*y(2)*(1-sum(y)/para(4)) + para(2)*(1+para(5)*(y(2)+y(4))/sum(y))*y(1) - para(3)*(1+para(6)*(y(1)+y(3))/sum(y))*y(2);
    dydt(3,1) = y(3)*(1-sum(y)/para(4)) - para(2)*(1+para(5)*(y(2)+y(4))/sum(y))*y(3) + para(3)*(1+para(6)*(y(1)+y(3))/sum(y))*y(4);
    dydt(4,1) = para(1)*y(4)*(1-sum(y)/para(4)) + para(2)*(1+para(5)*(y(2)+y(4))/sum(y))*y(3) - para(3)*(1+para(6)*(y(1)+y(3))/sum(y))*y(4);
    
end

function dydt = logistic_grow_trans_influence_of_M_reten_of_E_ND_VL_NL(t,y,para)
    
    dydt(1,1) = y(1)*(1-sum(y)/para(4)) - para(2)*y(1) + para(3)*(1+(para(6)*(y(1)+y(3))/sum(y)) - (para(5)*(y(2)+y(4))/sum(y)))*y(2);
    dydt(2,1) = para(1)*y(2)*(1-sum(y)/para(4)) + para(2)*y(1) - para(3)*(1+(para(6)*(y(1)+y(3))/sum(y)) - (para(5)*(y(2)+y(4))/sum(y)))*y(2);
    dydt(3,1) = y(3)*(1-sum(y)/para(4)) - para(2)*y(3) + para(3)*(1+(para(6)*(y(1)+y(3))/sum(y)) - (para(5)*(y(2)+y(4))/sum(y)))*y(4);
    dydt(4,1) = para(1)*y(4)*(1-sum(y)/para(4)) + para(2)*y(3) - para(3)*(1+(para(6)*(y(1)+y(3))/sum(y)) - (para(5)*(y(2)+y(4))/sum(y)))*y(4);
    
end

function dydt = logistic_grow_trans_influence_of_M_reten_of_M_ND_VL_NL(t,y,para)
    
    dydt(1,1) = y(1)*(1-sum(y)/para(4)) - para(2)*(1-para(5)*(y(1)+y(3))/sum(y))*y(1) + para(3)*(1+para(6)*(y(1)+y(3))/sum(y))*y(2);
    dydt(2,1) = para(1)*y(2)*(1-sum(y)/para(4)) + para(2)*(1-para(5)*(y(1)+y(3))/sum(y))*y(1) - para(3)*(1+para(6)*(y(1)+y(3))/sum(y))*y(2);
    dydt(3,1) = y(3)*(1-sum(y)/para(4)) - para(2)*(1-para(5)*(y(1)+y(3))/sum(y))*y(3) + para(3)*(1+para(6)*(y(1)+y(3))/sum(y))*y(4);
    dydt(4,1) = para(1)*y(4)*(1-sum(y)/para(4)) + para(2)*(1-para(5)*(y(1)+y(3))/sum(y))*y(3) - para(3)*(1+para(6)*(y(1)+y(3))/sum(y))*y(4);
    
end

function dydt = logistic_grow_trans_reten_of_E_reten_of_M_ND_VL_NL(t,y,para)
    
    dydt(1,1) = y(1)*(1-sum(y)/para(4)) - para(2)*(1-para(5)*(y(1)+y(3))/sum(y))*y(1) + para(3)*(1-para(6)*(y(2)+y(4))/sum(y))*y(2);
    dydt(2,1) = para(1)*y(2)*(1-sum(y)/para(4)) + para(2)*(1-para(5)*(y(1)+y(3))/sum(y))*y(1) - para(3)*(1-para(6)*(y(2)+y(4))/sum(y))*y(2);
    dydt(3,1) = y(3)*(1-sum(y)/para(4)) - para(2)*(1-para(5)*(y(1)+y(3))/sum(y))*y(3) + para(3)*(1-para(6)*(y(2)+y(4))/sum(y))*y(4);
    dydt(4,1) = para(1)*y(4)*(1-sum(y)/para(4)) + para(2)*(1-para(5)*(y(1)+y(3))/sum(y))*y(3) - para(3)*(1-para(6)*(y(2)+y(4))/sum(y))*y(4);
    
end

function dydt = logistic_grow_trans_reten_of_E_ND_VL_NL(t,y,para)
    
    dydt(1,1) = y(1)*(1-sum(y)/para(4)) - para(2)*y(1) + para(3)*(1-para(5)*(y(2)+y(4))/sum(y))*y(2);
    dydt(2,1) = para(1)*y(2)*(1-sum(y)/para(4)) + para(2)*y(1) - para(3)*(1-para(5)*(y(2)+y(4))/sum(y))*y(2);
    dydt(3,1) = y(3)*(1-sum(y)/para(4)) - para(2)*y(3) + para(3)*(1-para(5)*(y(2)+y(4))/sum(y))*y(4);
    dydt(4,1) = para(1)*y(4)*(1-sum(y)/para(4)) + para(2)*y(3) - para(3)*(1-para(5)*(y(2)+y(4))/sum(y))*y(4);
    
end

function dydt = logistic_grow_trans_reten_of_M_ND_VL_NL(t,y,para)
    
    dydt(1,1) = y(1)*(1-sum(y)/para(4)) - para(2)*(1-para(5)*(y(1)+y(3))/sum(y))*y(1) + para(3)*y(2);
    dydt(2,1) = para(1)*y(2)*(1-sum(y)/para(4)) + para(2)*(1-para(5)*(y(1)+y(3))/sum(y))*y(1) - para(3)*y(2);
    dydt(3,1) = y(3)*(1-sum(y)/para(4)) - para(2)*(1-para(5)*(y(1)+y(3))/sum(y))*y(3) + para(3)*y(4);
    dydt(4,1) = para(1)*y(4)*(1-sum(y)/para(4)) + para(2)*(1-para(5)*(y(1)+y(3))/sum(y))*y(3) - para(3)*y(4);
    
end

function dydt = three_state_transition_model_non_dim_H_M_combined(t,y,para)
    
    dydt(1,1) =         y(1)*(1 - para(2) - para(4)) + y(2)*para(3) + y(3)*para(5) - (y(1) + y(2) + para(1)*y(3))*y(1); 
    dydt(2,1) =         y(2)*(1 - para(3) - para(6)) + y(3)*para(7) + y(1)*para(2) - (y(1) + y(2) + para(1)*y(3))*y(2);
    dydt(3,1) = y(3)*(para(1) - para(5) - para(7)) + y(1)*para(4) + y(2)*para(6) - (y(1) + y(2) + para(1)*y(3))*y(3);
    
end

function dydt = three_state_transition_model_non_dim_E_H_combined(t,y,para)
    
    dydt(1,1) =         y(1)*(1 - para(2) - para(4)) + y(2)*para(3) + y(3)*para(5) - (y(1) + para(1)*y(2) + para(1)*y(3))*y(1); 
    dydt(2,1) = y(2)*(para(1) - para(3) - para(6)) + y(3)*para(7) + y(1)*para(2) - (y(1) + para(1)*y(2) + para(1)*y(3))*y(2);
    dydt(3,1) = y(3)*(para(1) - para(5) - para(7)) + y(1)*para(4) + y(2)*para(6) - (y(1) + para(1)*y(2) + para(1)*y(3))*y(3);
    
end

function dydt = logistic_grow_comp_ND_VL_NL(t,y,para)
    
    dydt(1,1) = y(1)*(1-((y(1)+y(3))+(y(2)+y(4))*para(2))/para(4));
    dydt(2,1) = para(1)*y(2)*(1-((y(2)+y(4))+(y(1)+y(3))*para(3))/para(5));
    dydt(3,1) = y(3)*(1-((y(1)+y(3))+(y(2)+y(4))*para(2))/para(4));
    dydt(4,1) = para(1)*y(4)*(1-((y(2)+y(4))+(y(1)+y(3))*para(3))/para(5));
    
end

function dydt = logistic_grow_comp_trans_independent_parameters_ND_VL_NL(t,y,para)
    
r_e = para(1);
r_mv = para(2);
r_ev = para(3);
t_me = para(4);
t_em = para(5);
t_mev = para(4);
t_emv = para(5);
c_e = para(8);
c_m = para(9);
c_ev = para(10);
c_mv = para(11);
K_e = para(12);
K_m = para(13);
K_ev = para(14);
K_mv = para(15);

    dydt(1,1) = y(1)*(1-(y(1)+y(2)*c_e+y(3)*c_mv+y(4)*c_ev)/K_m) - t_me*y(1) + t_em*y(2);
    
    dydt(2,1) = r_e*y(2)*(1-(y(1)*c_m+y(2)+y(3)*c_mv+y(4)*c_ev)/K_e) + t_me*y(1) - t_em*y(2);
    
    dydt(3,1) = r_mv*y(3)*(1-(y(1)*c_m+y(2)*c_e+y(3)+y(4)*c_ev)/K_mv) - t_mev*y(3) + t_emv*y(4);
    
    dydt(4,1) = r_ev*y(4)*(1-(y(1)*c_m+y(2)*c_e+y(3)*c_mv+y(4))/K_ev) + t_mev*y(3) - t_emv*y(4);
    
end

% function dydt = logistic_grow_comp_trans_independent_parameters_except_transition_rates_ND_VL_NL(t,y,para)
%     
% r_e = para(1);
% r_mv = para(2);
% r_ev = para(3);
% t_me = para(4);
% t_em = para(5);
% t_mev = para(4);
% t_emv = para(5);
% c_e = para(8);
% c_m = para(9);
% c_ev = para(10);
% c_mv = para(11);
% K_e = para(12);
% K_m = para(13);
% K_ev = para(14);
% K_mv = para(15);
% 
%     dydt(1,1) = y(1)*(1-(y(1)+y(2)*c_e+y(3)*c_mv+y(4)*c_ev)/K_m) - t_me*y(1) + t_em*y(2);
%     
%     dydt(2,1) = r_e*y(2)*(1-(y(1)*c_m+y(2)+y(3)*c_mv+y(4)*c_ev)/K_e) + t_me*y(1) - t_em*y(2);
%     
%     dydt(3,1) = r_mv*y(3)*(1-(y(1)*c_m+y(2)*c_e+y(3)+y(4)*c_ev)/K_mv) - t_mev*y(3) + t_emv*y(4);
%     
%     dydt(4,1) = r_ev*y(4)*(1-(y(1)*c_m+y(2)*c_e+y(3)*c_mv+y(4))/K_ev) + t_mev*y(3) - t_emv*y(4);
%     
% end

function dydt = logistic_grow_influence_ND_VL_NL(t,y,para)

    dydt(1,1) = y(1)*(1-((y(1)+y(3))+(y(2)+y(4)))/para(2) - para(3)*(y(2)+y(4))/sum(y));
    dydt(2,1) = para(1)*y(2)*(1-((y(2)+y(4))+(y(1)+y(3)))/para(2) - para(4)*(y(1)+y(3))/sum(y));
    dydt(3,1) = y(3)*(1-((y(1)+y(3))+(y(2)+y(4)))/para(2) - para(3)*(y(2)+y(4))/sum(y));
    dydt(4,1) = para(1)*y(4)*(1-((y(2)+y(4))+(y(1)+y(3)))/para(2) - para(4)*(y(1)+y(3))/sum(y));

end

function dydt = logistic_grow_influence_trans_ND_VL_NL(t,y,para)

    dydt(1,1) = y(1)*(1-((y(1)+y(3))+(y(2)+y(4)))/para(2) - para(3)*(y(2)+y(4))/sum(y)) - para(5)* y(1) + para(6)* y(2);
    dydt(2,1) = para(1)*y(2)*(1-((y(2)+y(4))+(y(1)+y(3)))/para(2) - para(4)*(y(1)+y(3))/sum(y)) + para(5)* y(1) - para(6)* y(2);
    dydt(3,1) = y(3)*(1-((y(1)+y(3))+(y(2)+y(4)))/para(2) - para(3)*(y(2)+y(4))/sum(y)) - para(5)* y(3) + para(6)* y(4);
    dydt(4,1) = para(1)*y(4)*(1-((y(2)+y(4))+(y(1)+y(3)))/para(2) - para(4)*(y(1)+y(3))/sum(y)) + para(5)* y(3) - para(6)* y(4);

end

function dydt = logistic_grow_comp_trans_asymmetric_growth_ND_VL_NL(t,y,para)
    
r_e = para(1);
r_mv = para(2);
r_ev = para(3);
t_me = para(4);
t_em = para(5);
t_mev = para(4);
t_emv = para(5);
c_e = para(6);
c_m = para(7);
c_ev = para(6);
c_mv = para(7);
K_e = para(8);
K_m = para(9);
K_ev = para(8);
K_mv = para(9);

    dydt(1,1) = y(1)*(1-(y(1)+y(2)*c_e+y(3)*c_mv+y(4)*c_ev)/K_m) - t_me*y(1) + t_em*y(2);
    
    dydt(2,1) = r_e*y(2)*(1-(y(1)*c_m+y(2)+y(3)*c_mv+y(4)*c_ev)/K_e) + t_me*y(1) - t_em*y(2);
    
    dydt(3,1) = r_mv*y(3)*(1-(y(1)*c_m+y(2)*c_e+y(3)+y(4)*c_ev)/K_mv) - t_mev*y(3) + t_emv*y(4);
    
    dydt(4,1) = r_ev*y(4)*(1-(y(1)*c_m+y(2)*c_e+y(3)*c_mv+y(4))/K_ev) + t_mev*y(3) - t_emv*y(4);
    
end

function dydt = logistic_grow_trans_infl_of_M_reten_of_E_inde_growth_ND_VL_NL(t,y,para)
    
    dydt(1,1) = y(1)*(1-sum(y)/para(4)) - para(2)*y(1) + para(3)*(1+(para(6)*(y(1)+y(3))/sum(y)) - (para(5)*(y(2)+y(4))/sum(y)))*y(2);
    dydt(2,1) = para(1)*y(2)*(1-sum(y)/para(4)) + para(2)*y(1) - para(3)*(1+(para(6)*(y(1)+y(3))/sum(y)) - (para(5)*(y(2)+y(4))/sum(y)))*y(2);
    dydt(3,1) = y(3)*(1-sum(y)/para(4)) - para(2)*y(3) + para(3)*(1+(para(6)*(y(1)+y(3))/sum(y)) - (para(5)*(y(2)+y(4))/sum(y)))*y(4);
    dydt(4,1) = para(1)*y(4)*(1-sum(y)/para(4)) + para(2)*y(3) - para(3)*(1+(para(6)*(y(1)+y(3))/sum(y)) - (para(5)*(y(2)+y(4))/sum(y)))*y(4);
    
end

function dydt = logistic_grow_comp_trans_VL_NL(t,y,para)
    
    dydt(1,1) = para(1)*y(1)*(1-((y(1)+y(3))+(y(2)+y(4))*para(5))/para(7)) - para(3)*y(1) + para(4)*y(2);
    dydt(2,1) = para(2)*y(2)*(1-((y(2)+y(4))+(y(1)+y(3))*para(6))/para(8)) + para(3)*y(1) - para(4)*y(2);
    dydt(3,1) = para(1)*y(3)*(1-((y(1)+y(3))+(y(2)+y(4))*para(5))/para(7)) - para(3)*y(3) + para(4)*y(4);
    dydt(4,1) = para(2)*y(4)*(1-((y(2)+y(4))+(y(1)+y(3))*para(6))/para(8)) + para(3)*y(3) - para(4)*y(4);
    
end

function dydt = logistic_grow_trans_VL_NL(t,y,para)
    
    dydt(1,1) = para(1)*y(1)*(1-sum(y)/para(5)) - para(3)*y(1) + para(4)*y(2);
    dydt(2,1) = para(2)*y(2)*(1-sum(y)/para(5)) + para(3)*y(1) - para(4)*y(2);
    dydt(3,1) = para(1)*y(3)*(1-sum(y)/para(5)) - para(3)*y(3) + para(4)*y(4);
    dydt(4,1) = para(2)*y(4)*(1-sum(y)/para(5)) + para(3)*y(3) - para(4)*y(4);
    
end

function dydt = logistic_grow_trans_influence_of_M_reten_of_E_VL_NL(t,y,para)
    
    dydt(1,1) = para(1)*y(1)*(1-sum(y)/para(5)) - para(3)*y(1) + para(4)*(1+(para(7)*(y(1)+y(3))/sum(y)) - (para(6)*(y(2)+y(4))/sum(y)))*y(2);
    dydt(2,1) = para(2)*y(2)*(1-sum(y)/para(5)) + para(3)*y(1) - para(4)*(1+(para(7)*(y(1)+y(3))/sum(y)) - (para(6)*(y(2)+y(4))/sum(y)))*y(2);
    dydt(3,1) = para(1)*y(3)*(1-sum(y)/para(5)) - para(3)*y(3) + para(4)*(1+(para(7)*(y(1)+y(3))/sum(y)) - (para(6)*(y(2)+y(4))/sum(y)))*y(4);
    dydt(4,1) = para(2)*y(4)*(1-sum(y)/para(5)) + para(3)*y(3) - para(4)*(1+(para(7)*(y(1)+y(3))/sum(y)) - (para(6)*(y(2)+y(4))/sum(y)))*y(4);
    
end

function dydt = logistic_grow_trans_reten_of_E_reten_of_M_VL_NL(t,y,para)
    
    dydt(1,1) = para(1)*y(1)*(1-sum(y)/para(5)) - para(3)*(1-para(6)*(y(1)+y(3))/sum(y))*y(1) + para(4)*(1-para(7)*(y(2)+y(4))/sum(y))*y(2);
    dydt(2,1) = para(2)*y(2)*(1-sum(y)/para(5)) + para(3)*(1-para(6)*(y(1)+y(3))/sum(y))*y(1) - para(4)*(1-para(7)*(y(2)+y(4))/sum(y))*y(2);
    dydt(3,1) = para(1)*y(3)*(1-sum(y)/para(5)) - para(3)*(1-para(6)*(y(1)+y(3))/sum(y))*y(3) + para(4)*(1-para(7)*(y(2)+y(4))/sum(y))*y(4);
    dydt(4,1) = para(2)*y(4)*(1-sum(y)/para(5)) + para(3)*(1-para(6)*(y(1)+y(3))/sum(y))*y(3) - para(4)*(1-para(7)*(y(2)+y(4))/sum(y))*y(4);
    
end

function dydt = logistic_grow_trans_reten_of_E_VL_NL(t,y,para)
    
    dydt(1,1) = para(1)*y(1)*(1-sum(y)/para(5)) - para(3)*y(1) + para(4)*(1-para(6)*(y(2)+y(4))/sum(y))*y(2);
    dydt(2,1) = para(2)*y(2)*(1-sum(y)/para(5)) + para(3)*y(1) - para(4)*(1-para(6)*(y(2)+y(4))/sum(y))*y(2);
    dydt(3,1) = para(1)*y(3)*(1-sum(y)/para(5)) - para(3)*y(3) + para(4)*(1-para(6)*(y(2)+y(4))/sum(y))*y(4);
    dydt(4,1) = para(2)*y(4)*(1-sum(y)/para(5)) + para(3)*y(3) - para(4)*(1-para(6)*(y(2)+y(4))/sum(y))*y(4);
    
end

function dydt = logistic_grow_trans_reten_of_M_VL_NL(t,y,para)
    
    dydt(1,1) = para(1)*y(1)*(1-sum(y)/para(5)) - para(3)*(1-para(6)*(y(1)+y(3))/sum(y))*y(1) + para(4)*y(2);
    dydt(2,1) = para(2)*y(2)*(1-sum(y)/para(5)) + para(3)*(1-para(6)*(y(1)+y(3))/sum(y))*y(1) - para(4)*y(2);
    dydt(3,1) = para(1)*y(3)*(1-sum(y)/para(5)) - para(3)*(1-para(6)*(y(1)+y(3))/sum(y))*y(3) + para(4)*y(4);
    dydt(4,1) = para(2)*y(4)*(1-sum(y)/para(5)) + para(3)*(1-para(6)*(y(1)+y(3))/sum(y))*y(3) - para(4)*y(4);
    
end

function dydt = logistic_grow_influence_trans_VL_NL(t,y,para)

    dydt(1,1) = para(1)*y(1)*(1-((y(1)+y(3))+(y(2)+y(4)))/para(3) - para(4)*(y(2)+y(4))/sum(y)) - para(6)* y(1) + para(7)* y(2);
    dydt(2,1) = para(2)*y(2)*(1-((y(2)+y(4))+(y(1)+y(3)))/para(3) - para(5)*(y(1)+y(3))/sum(y)) + para(6)* y(1) - para(7)* y(2);
    dydt(3,1) = para(1)*y(3)*(1-((y(1)+y(3))+(y(2)+y(4)))/para(3) - para(4)*(y(2)+y(4))/sum(y)) - para(6)* y(3) + para(7)* y(4);
    dydt(4,1) = para(2)*y(4)*(1-((y(2)+y(4))+(y(1)+y(3)))/para(3) - para(5)*(y(1)+y(3))/sum(y)) + para(6)* y(3) - para(7)* y(4);

end

function dydt = logistic_grow_VL_NL(t,y,para)
    
    dydt(1,1) = para(1)*y(1)*(1-((y(1)+y(3))+(y(2)+y(4)))/para(3));
    dydt(2,1) = para(2)*y(2)*(1-((y(2)+y(4))+(y(1)+y(3)))/para(3));
    dydt(3,1) = para(1)*y(3)*(1-((y(1)+y(3))+(y(2)+y(4)))/para(3));
    dydt(4,1) = para(2)*y(4)*(1-((y(2)+y(4))+(y(1)+y(3)))/para(3));
    
end

function dydt = logistic_grow_comp_VL_NL(t,y,para)
    
    dydt(1,1) = para(1)*y(1)*(1-((y(1)+y(3))+(y(2)+y(4))*para(3))/para(5));
    dydt(2,1) = para(2)*y(2)*(1-((y(2)+y(4))+(y(1)+y(3))*para(4))/para(6));
    dydt(3,1) = para(1)*y(3)*(1-((y(1)+y(3))+(y(2)+y(4))*para(3))/para(5));
    dydt(4,1) = para(2)*y(4)*(1-((y(2)+y(4))+(y(1)+y(3))*para(4))/para(6));
    
end

function dydt = logistic_grow_influence_VL_NL(t,y,para)

    dydt(1,1) = para(1)*y(1)*(1-((y(1)+y(3))+(y(2)+y(4)))/para(3) - para(4)*(y(2)+y(4))/sum(y));
    dydt(2,1) = para(2)*y(2)*(1-((y(2)+y(4))+(y(1)+y(3)))/para(3) - para(5)*(y(1)+y(3))/sum(y));
    dydt(3,1) = para(1)*y(3)*(1-((y(1)+y(3))+(y(2)+y(4)))/para(3) - para(4)*(y(2)+y(4))/sum(y));
    dydt(4,1) = para(2)*y(4)*(1-((y(2)+y(4))+(y(1)+y(3)))/para(3) - para(5)*(y(1)+y(3))/sum(y));

end

function dydt = logistic_grow_trans_influence_of_M_reten_of_M_VL_NL(t,y,para)
    
    dydt(1,1) = para(1)*y(1)*(1-sum(y)/para(5)) - para(3)*(1-para(6)*(y(1)+y(3))/sum(y))*y(1) + para(4)*(1+para(7)*(y(1)+y(3))/sum(y))*y(2);
    dydt(2,1) = para(2)*y(2)*(1-sum(y)/para(5)) + para(3)*(1-para(6)*(y(1)+y(3))/sum(y))*y(1) - para(4)*(1+para(7)*(y(1)+y(3))/sum(y))*y(2);
    dydt(3,1) = para(1)*y(3)*(1-sum(y)/para(5)) - para(3)*(1-para(6)*(y(1)+y(3))/sum(y))*y(3) + para(4)*(1+para(7)*(y(1)+y(3))/sum(y))*y(4);
    dydt(4,1) = para(2)*y(4)*(1-sum(y)/para(5)) + para(3)*(1-para(6)*(y(1)+y(3))/sum(y))*y(3) - para(4)*(1+para(7)*(y(1)+y(3))/sum(y))*y(4);
    
end



function dydt = logistic_grow_comp_trans_reten_of_E_ND_VL_NL(t,y,para)
    
    dydt(1,1) = y(1)*(1-((y(1)+y(3))+(y(2)+y(4))*para(4))/para(6)) - para(2)*y(1) + para(3)*(1-para(8)*(y(2)+y(4))/sum(y))*y(2);
    dydt(2,1) = para(1)*y(2)*(1-((y(2)+y(4))+(y(1)+y(3))*para(5))/para(7)) + para(2)*y(1) - para(3)*(1-para(8)*(y(2)+y(4))/sum(y))*y(2);
    dydt(3,1) = y(3)*(1-((y(1)+y(3))+(y(2)+y(4))*para(4))/para(6)) - para(2)*y(3) + para(3)*(1-para(8)*(y(2)+y(4))/sum(y))*y(4);
    dydt(4,1) = para(1)*y(4)*(1-((y(2)+y(4))+(y(1)+y(3))*para(5))/para(7)) + para(2)*y(3) - para(3)*(1-para(8)*(y(2)+y(4))/sum(y))*y(4);
    
end

function dydt = logistic_grow_comp_trans_reten_of_M_ND_VL_NL(t,y,para)
    
    dydt(1,1) = y(1)*(1-((y(1)+y(3))+(y(2)+y(4))*para(4))/para(6)) - para(2)*(1-para(8)*(y(1)+y(3))/sum(y))*y(1) + para(3)*y(2);
    dydt(2,1) = para(1)*y(2)*(1-((y(2)+y(4))+(y(1)+y(3))*para(5))/para(7)) + para(2)*(1-para(8)*(y(1)+y(3))/sum(y))*y(1) - para(3)*y(2);
    dydt(3,1) = y(3)*(1-((y(1)+y(3))+(y(2)+y(4))*para(4))/para(6)) - para(2)*(1-para(8)*(y(1)+y(3))/sum(y))*y(3) + para(3)*y(4);
    dydt(4,1) = para(1)*y(4)*(1-((y(2)+y(4))+(y(1)+y(3))*para(5))/para(7)) + para(2)*(1-para(8)*(y(1)+y(3))/sum(y))*y(3) - para(3)*y(4);
    
end

function dydt = logistic_grow_comp_trans_influence_of_M_ND_VL_NL(t,y,para)
    
    dydt(1,1) = y(1)*(1-((y(1)+y(3))+(y(2)+y(4))*para(4))/para(6)) - para(2)*y(1) + para(3)*(1+para(8)*(y(1)+y(3))/sum(y))*y(2);
    dydt(2,1) = para(1)*y(2)*(1-((y(2)+y(4))+(y(1)+y(3))*para(5))/para(7)) + para(2)*y(1) - para(3)*(1+para(8)*(y(1)+y(3))/sum(y))*y(2);
    dydt(3,1) = y(3)*(1-((y(1)+y(3))+(y(2)+y(4))*para(4))/para(6)) - para(2)*y(3) + para(3)*(1+para(8)*(y(1)+y(3))/sum(y))*y(4);
    dydt(4,1) = para(1)*y(4)*(1-((y(2)+y(4))+(y(1)+y(3))*para(5))/para(7)) + para(2)*y(3) - para(3)*(1+para(8)*(y(1)+y(3))/sum(y))*y(4);
    
end

function dydt = logistic_grow_influence_trans_retention_of_E_ND_VL_NL(t,y,para)

    dydt(1,1) = y(1)*(1-((y(1)+y(3))+(y(2)+y(4)))/para(2) - para(3)*(y(2)+y(4))/sum(y)) - para(5)* y(1) + para(6)*(1-para(7)*(y(2)+y(4))/sum(y))* y(2);
    dydt(2,1) = para(1)*y(2)*(1-((y(2)+y(4))+(y(1)+y(3)))/para(2) - para(4)*(y(1)+y(3))/sum(y)) + para(5)* y(1) - para(6)*(1-para(7)*(y(2)+y(4))/sum(y))* y(2);
    dydt(3,1) = y(3)*(1-((y(1)+y(3))+(y(2)+y(4)))/para(2) - para(3)*(y(2)+y(4))/sum(y)) - para(5)* y(3) + para(6)*(1-para(7)*(y(2)+y(4))/sum(y))* y(4);
    dydt(4,1) = para(1)*y(4)*(1-((y(2)+y(4))+(y(1)+y(3)))/para(2) - para(4)*(y(1)+y(3))/sum(y)) + para(5)* y(3) - para(6)*(1-para(7)*(y(2)+y(4))/sum(y))* y(4);

end


function dydt = logistic_grow_influence_trans_retention_of_M_ND_VL_NL(t,y,para)

    dydt(1,1) = y(1)*(1-((y(1)+y(3))+(y(2)+y(4)))/para(2) - para(3)*(y(2)+y(4))/sum(y)) - para(5)*(1-para(7)*(y(1)+y(3))/sum(y))* y(1) + para(6)* y(2);
    dydt(2,1) = para(1)*y(2)*(1-((y(2)+y(4))+(y(1)+y(3)))/para(2) - para(4)*(y(1)+y(3))/sum(y)) + para(5)*(1-para(7)*(y(1)+y(3))/sum(y))* y(1) - para(6)* y(2);
    dydt(3,1) = y(3)*(1-((y(1)+y(3))+(y(2)+y(4)))/para(2) - para(3)*(y(2)+y(4))/sum(y)) - para(5)*(1-para(7)*(y(1)+y(3))/sum(y))* y(3) + para(6)* y(4);
    dydt(4,1) = para(1)*y(4)*(1-((y(2)+y(4))+(y(1)+y(3)))/para(2) - para(4)*(y(1)+y(3))/sum(y)) + para(5)*(1-para(7)*(y(1)+y(3))/sum(y))* y(3) - para(6)* y(4);

end

function dydt = logistic_grow_influence_trans_influence_of_M_ND_VL_NL(t,y,para)

    dydt(1,1) = y(1)*(1-((y(1)+y(3))+(y(2)+y(4)))/para(2) - para(3)*(y(2)+y(4))/sum(y)) - para(5)* y(1) + para(6)*(1+para(7)*(y(1)+y(3))/sum(y))* y(2);
    dydt(2,1) = para(1)*y(2)*(1-((y(2)+y(4))+(y(1)+y(3)))/para(2) - para(4)*(y(1)+y(3))/sum(y)) + para(5)* y(1) - para(6)*(1+para(7)*(y(1)+y(3))/sum(y))* y(2);
    dydt(3,1) = y(3)*(1-((y(1)+y(3))+(y(2)+y(4)))/para(2) - para(3)*(y(2)+y(4))/sum(y)) - para(5)* y(3) + para(6)*(1+para(7)*(y(1)+y(3))/sum(y))* y(4);
    dydt(4,1) = para(1)*y(4)*(1-((y(2)+y(4))+(y(1)+y(3)))/para(2) - para(4)*(y(1)+y(3))/sum(y)) + para(5)* y(3) - para(6)*(1+para(7)*(y(1)+y(3))/sum(y))* y(4);

end

function dydt = linear_three_state_transition_model_non_dim_H_M_combined(t,y,para)
    
    dydt(1,1) =         y(1)*(1 - para(2)) + y(2)*para(3) - (y(1) + y(2) + para(1)*y(3))*y(1); 
    dydt(2,1) =         y(2)*(1 - para(3) - para(4)) + y(3)*para(5) + y(1)*para(2) - (y(1) + y(2) + para(1)*y(3))*y(2);
    dydt(3,1) =         y(3)*(para(1) - para(5)) + y(2)*para(4) - (y(1) + y(2) + para(1)*y(3))*y(3);
    
end

function dydt = linear_three_state_transition_model_non_dim_E_H_combined(t,y,para)
    
    dydt(1,1) =         y(1)*(1 - para(2)) + y(2)*para(3) - (y(1) + para(1)*y(2) + para(1)*y(3))*y(1); 
    dydt(2,1) =         y(2)*(para(1) - para(3) - para(4)) + y(3)*para(5) + y(1)*para(2) - (y(1) + para(1)*y(2) + para(1)*y(3))*y(2);
    dydt(3,1) =         y(3)*(para(1) - para(5)) + y(2)*para(4) - (y(1) + para(1)*y(2) + para(1)*y(3))*y(3);
    
end

function dydt = linear_three_state_logistic_grow_trans_ND_VL_NL_E_H_combined(t,y,para)
    
    dydt(1,1) = y(1)*(1-sum(y)/para(6)) - para(2)*y(1) + para(3)*y(2);
    dydt(2,1) = para(1)*y(2)*(1-sum(y)/para(6)) + para(2)*y(1) + para(5)*y(3) - (para(3)+para(4))*y(2);
    dydt(3,1) = para(1)*y(3)*(1-sum(y)/para(6)) + para(4)*y(2) - para(5)*y(3);

    dydt(4,1) = y(4)*(1-sum(y)/para(6)) - para(2)*y(4) + para(3)*y(5);
    dydt(5,1) = para(1)*y(5)*(1-sum(y)/para(6)) + para(2)*y(4) + para(5)*y(6) - (para(3)+para(4))*y(5);
    dydt(6,1) = para(1)*y(6)*(1-sum(y)/para(6)) + para(4)*y(5) - para(5)*y(6);
    
end

function dydt = linear_three_state_logistic_grow_trans_ND_VL_NL_H_M_combined(t,y,para)
    
    dydt(1,1) = y(1)*(1-sum(y)/para(6)) - para(2)*y(1) + para(3)*y(2);
    dydt(2,1) = y(2)*(1-sum(y)/para(6)) + para(2)*y(1) + para(5)*y(3) - (para(3)+para(4))*y(2);
    dydt(3,1) = para(1)*y(3)*(1-sum(y)/para(6)) + para(4)*y(2) - para(5)*y(3);

    dydt(4,1) = y(4)*(1-sum(y)/para(6)) - para(2)*y(4) + para(3)*y(5);
    dydt(5,1) = y(5)*(1-sum(y)/para(6)) + para(2)*y(4) + para(5)*y(6) - (para(3)+para(4))*y(5);
    dydt(6,1) = para(1)*y(6)*(1-sum(y)/para(6)) + para(4)*y(5) - para(5)*y(6);
    
end

function dydt = three_state_logistic_grow_trans_ND_VL_NL_E_H_combined(t,y,para)
    
    dydt(1,1) = y(1)*(1-sum(y)/para(8)) - para(2)*y(1) + para(3)*y(2) - para(6)*y(1) + para(7)*y(3);
    dydt(2,1) = para(1)*y(2)*(1-sum(y)/para(8)) + para(2)*y(1) + para(5)*y(3) - (para(3)+para(4))*y(2);
    dydt(3,1) = para(1)*y(3)*(1-sum(y)/para(8)) + para(4)*y(2) - para(5)*y(3) + para(6)*y(1) - para(7)*y(3);

    dydt(4,1) = y(4)*(1-sum(y)/para(8)) - para(2)*y(4) + para(3)*y(5) - para(6)*y(4) + para(7)*y(6);
    dydt(5,1) = para(1)*y(5)*(1-sum(y)/para(8)) + para(2)*y(4) + para(5)*y(6) - (para(3)+para(4))*y(5);
    dydt(6,1) = para(1)*y(6)*(1-sum(y)/para(8)) + para(4)*y(5) - para(5)*y(6) + para(6)*y(4) - para(7)*y(6);
    
end

function dydt = three_state_logistic_grow_trans_ND_VL_NL_H_M_combined(t,y,para)
    
    dydt(1,1) = y(1)*(1-sum(y)/para(8)) - para(2)*y(1) + para(3)*y(2) - para(6)*y(1) + para(7)*y(3);
    dydt(2,1) = y(2)*(1-sum(y)/para(8)) + para(2)*y(1) + para(5)*y(3) - (para(3)+para(4))*y(2);
    dydt(3,1) = para(1)*y(3)*(1-sum(y)/para(8)) + para(4)*y(2) - para(5)*y(3) + para(6)*y(1) - para(7)*y(3);

    dydt(4,1) = y(4)*(1-sum(y)/para(8)) - para(2)*y(4) + para(3)*y(5) - para(6)*y(4) + para(7)*y(6);
    dydt(5,1) = y(5)*(1-sum(y)/para(8)) + para(2)*y(4) + para(5)*y(6) - (para(3)+para(4))*y(5);
    dydt(6,1) = para(1)*y(6)*(1-sum(y)/para(8)) + para(4)*y(5) - para(5)*y(6) + para(6)*y(4) - para(7)*y(6);
    
end

% function dydt = three_state_transition_model(t,y,para)
%     
%     dydt(1,1) = - y(1)*(para(1) + para(3) + para(5)) + y(2)*para(2); 
%     dydt(2,1) = para(1)*y(1) - (para(2) + para(4) + para(6))*y(2) ;
%     dydt(3,1) = y(1)*para(3) + para(6)*y(2) + para(8)*y(4) - y(3)*(para(9)+para(11)+para(7));
%     dydt(4,1) = y(2)*para(4) + para(5)*y(1) + para(7)*y(3) - y(4)*(para(10)+para(12)+para(8));
%     dydt(5,1) = y(4)*para(12) + para(9)*y(3) + para(14)*y(6) - y(5)*(para(13));
%     dydt(6,1) = y(4)*para(10) + para(11)*y(3) + para(13)*y(5) - y(6)*(para(14));
%     
% end