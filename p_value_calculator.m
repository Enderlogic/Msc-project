function [p_value] = p_value_calculator(data, data_rep, criteria)
% criteria: number_of_zero; norm_of_gradient
if strcmp(criteria, 'number_of_zero')
    cri_data = sum(diff(data > 0) ~= 0);
    cri_rep = sum(diff(data_rep > 0) ~= 0);
elseif strcmp(criteria, 'norm_of_gradient')
    cri_data = norm(gradient(data));
    [~, ygrad] = gradient(data_rep);
    cri_rep = vecnorm(ygrad);
else
    error('The criteria is invalid!')
end
p_value = 2 * min(sum(cri_rep > cri_data) / size(data_rep, 2), 1 - sum(cri_rep > cri_data) / size(data_rep, 2));
end

