clear
close all
%% Setting for data, model and PPC
option_data = 'real'; % choose ficticious data (fic) or real data (real)
% data_name = {'Periodic'};
data_name = {'3'}; % choose the real data from the 1st file to the 14th file
% model_cov = {'Matern', '3'};
ratio_train = 0.8; % set the ratio of training set among the whole set
num_rep = 50; % set the amount of samples

model_cov = {'prod', 'Periodic','LIN'}; % 'SE'; 'Periodic'; 'Matern', '3'; 'Matern', '5', 'Exp', 'RQ' for single kernel
                                        % 'SE'; 'LIN'; 'Periodic' for 'sum'
                                        % and 'prod' mixing kernel
                                        % example: 'sum', 'SE', 'LIN' =
                                        % SE + LIN kernel
sample_type = 'both'; % 'para'; 'para&obs'; 'both'
save = false; % true: save the parameters and data; false: not save

criteria = 'chi_square'; % choose one criteria for PPC or MMD ('number_of_zero'; 'norm_of_gradient'; 'chi_square'; 'mmd')
%% Load existing model otherwise generate new model
filename = strcat('model\', option_data, '_', strjoin(data_name), '_', strjoin(model_cov), '_', sample_type, '.mat');
if exist(filename, 'file') == 2
    load(filename)
else
    [data_test, data_test_rep, sample_hyp, x_train, data_train, x_test] = model_generation(option_data, data_name, ratio_train, model_cov, sample_type, num_rep, save);
end
%% Compute posterior p-value
p_value = p_value_calculator(data_test, data_test_rep, criteria, model_cov, sample_hyp, x_train, data_train, x_test);
% disp(['p-value:',num2str(p_value)])