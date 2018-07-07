clear
close all
%% Setting for data, model and PPC
option_data = 'fic'; % choose ficticious data (fic) or real data (real)

% set the hyperparameters of ficticious data set (not all of them will be used)
hyp_data.magnitude = 10; % squared magnitude
hyp_data.lengthScale = 2;
hyp_data.period = 2;
hyp_data.sigma = 0.1;
data_cov = {'Periodic'}; % choose a gaussian kernel for data ('SE'; 'Periodic'; 'Matern', 1; 'Matern', 2)
num_data = 1000; % choose sampling number
random_seed = 'on'; % choose random seed for data, if seed is 'on', the data will be the same

ratio_train = 0.6; % set the ratio of training set among the whole set

model_cov = {'Matern', 2}; %choose the kernel of model ('SE'; 'Periodic'; 'Matern', 1; 'Matern', 2)

num_rep = 100; % set the amount of samples
sample_type = 'para'; % choose the type of sampling ('para': p(y_rep|'para', x_heldout); 'para&obs': p(y_rep|'para', y_obs, x_heldout))

criteria = 'chi_square'; % choose one criteria for PPC ('number_of_mean'; 'norm_of_gradient'; 'chi_square')
%% Data generation
if strcmp(option_data, 'fic')
    [x, data] = data_generation(num_data, random_seed, hyp_data, data_cov);
elseif strcmp(option_data, 'real')
    data = xlsread('mean-daily-temperature-fisher-ri.xlsx');
    num_data = size(data, 1);
    x = (1 : num_data)';
else
    error('The source of observation is invalid')
end

num_data_train = round(num_data * ratio_train);
num_data_test = num_data - num_data_train;
data_train = data(1 : num_data_train);
x_train = x(1 : num_data_train);
data_test = data(num_data_train + 1 : num_data);
x_test = x(num_data_train + 1 : num_data);
disp('Data generation complete!')
%% Compute MLE solution using gpml
%initilise the hyperparameters for MLE
hyp_mle.magnitude = 10; %squared magnitude
hyp_mle.lengthScale = 3;
hyp_mle.period = 2;
hyp_mle.sigma = 0.1;

hyp_mle = mle_gpml(hyp_mle, x_train, data_train, model_cov);
disp('MLE solution complete!')
%% Draw samples from posterior distribution of hyperparameters
% set prior by MLE result
prior.lengthScale = prior_loggaussian('s2', hyp_mle.lengthScale, 'mu', hyp_mle.lengthScale);
% prior.lengthScale = prior_gaussian('s2', 0.1, 'mu', 10);
prior.magnSigma2 = prior_loggaussian('s2', hyp_mle.magnitude, 'mu', hyp_mle.magnitude);
% prior.magnSigma2 = prior_gaussian('s2', hyp_mle.magnitude, 'mu', hyp_mle.magnitude);
prior.sigma2 = prior_loggaussian('s2', hyp_mle.sigma, 'mu', hyp_mle.sigma);
% prior.sigma2 = prior_gaussian('s2', hyp_mle.sigma, 'mu', hyp_mle.sigma);
prior.period = prior_loggaussian('s2', hyp_mle.period, 'mu', hyp_mle.period);

sample_hyp = sample_parameter(num_rep, x_train, data_train, prior, model_cov);
disp('Parameter sampling complete!')
%% Generate replicated data
data_test_rep = sample_data(num_rep, sample_hyp, model_cov, sample_type, x_train, data_train, x_test);
disp('Data sampling complete!')
%% Compute posterior p-value
p_value = p_value_calculator(data_test, data_test_rep, criteria, sample_type, model_cov, sample_hyp, x_train, data_train, x_test);
disp(['p-value:',num2str(p_value)])