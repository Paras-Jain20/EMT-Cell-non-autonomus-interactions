function out = model_ODE

out{1} = [];
out{2} = []; % cell fraction dynamics
out{3} = @exp_grow_with_trans_cell_fractions_non_dim; % cell fraction dynamics non dim
out{4} = @exp_grow_influence_with_trans_cell_fraction_non_dim; % cell fraction dynamics non dim
out{5} = [];
out{6} = @exp_grow_cell_fraction_non_dim;
out{7} = @exp_grow_influence_cell_fraction_non_dim;
out{8} = [];
out{9} = [];
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
out{28} = [];
out{29} = @logistic_grow_comp_ND_VL_NL;
out{30} = [];
out{31} = [];
out{32} = @logistic_grow_influence_ND_VL_NL;
out{33} = @logistic_grow_influence_trans_ND_VL_NL;

end


function dydt = exp_grow_with_trans_cell_numbers(t,y,para)
dydt(1,1) = y(1)*(para(1)-para(3))+y(2)*para(4);
dydt(2,1) = y(2)*(para(2)-para(4))+y(1)*para(3);
end


function dydt = exp_grow_influence_with_trans_cell_fraction_non_dim(t,y,para)
    dydt(1,1) = y(1)*(1-para(4)*y(2))  - para(2)*y(1) + para(3)*y(2) - (y(1)*(1-para(4)*y(2)) + para(1)*y(2)*(1-para(5)*y(1)))* y(1);
    dydt(2,1) = para(1)*y(2)*(1-para(5)*y(1))  + para(2)*y(1) - para(3)*y(2) - (y(1)*(1-para(4)*y(2)) + para(1)*y(2)*(1-para(5)*y(1)))* y(2);
end

function dydt = exp_grow_influence_cell_fraction_non_dim(t,y,para)
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

function dydt = logistic_grow_comp_ND_VL_NL(t,y,para)
    
    dydt(1,1) = y(1)*(1-((y(1)+y(3))+(y(2)+y(4))*para(2))/para(4));
    dydt(2,1) = para(1)*y(2)*(1-((y(2)+y(4))+(y(1)+y(3))*para(3))/para(5));
    dydt(3,1) = y(3)*(1-((y(1)+y(3))+(y(2)+y(4))*para(2))/para(4));
    dydt(4,1) = para(1)*y(4)*(1-((y(2)+y(4))+(y(1)+y(3))*para(3))/para(5));
    
end

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