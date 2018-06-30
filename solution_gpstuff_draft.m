clear
close all
%% Generate observed data
A_before = datetime(1987, 12, 1, 12, 0, 0);
A = datetime(1988, 1, 1, 12, 0, 0);
B = datetime(1991, 12, 31, 12, 0, 0);
B_next = datetime(1992, 1, 31, 12, 0, 0);

num_data = 1461;
num_data_train = round(num_data * 0.8);
num_data_test = num_data - num_data_train;

date = linspace(A, B, num_data)';
date_num = datenum(date);
date_train_num = date_num(1 : num_data_train);
date_test_num = date_num(num_data_train + 1 : num_data);

magnitude = 10;
lengthScale = 50;
lengthScale_sexp = 500000000;
sigma = 0.001;
period = 380;

cov_sexp = @(X1, X2, magnitude, lengthScale) magnitude * exp(- (X1 * ones(1, size(X2, 1)) - ones(size(X1, 1), 1) * X2').^2/2/(lengthScale^2));
cov_exp = @(X1, X2, magnitude, lengthScale) magnitude * exp(- sqrt((X1 * ones(1, size(X2, 1)) - ones(size(X1, 1), 1) * X2').^2/(lengthScale^2)));
cov_periodic = @(X1, X2, magnitude, period, lengthScale) magnitude * exp(- 2 * sin(pi * (X1 * ones(1, size(X2, 1)) - ones(size(X1, 1), 1) * X2')/ period).^2 / (lengthScale)^2);% - (X1 * ones(1, size(X2, 1)) - ones(size(X1, 1), 1) * X2').^2 / 2 / (lengthScale_sexp)^2);
cov = cov_periodic(date_num, date_num, magnitude, period, lengthScale);
% cov = cov_sexp(date_num, date_num, magnitude, lengthScale);

data = mvnrnd(zeros(1, num_data), cov)';
data = mvnrnd(data, sigma * eye(num_data))';
% data = data - mean(data);

data_train = data(1 : num_data_train);
data_test = data(num_data_train + 1 : num_data);
%% Import observed data
% A_before = datetime(1987, 12, 1, 12, 0, 0);
% A = datetime(1988, 1, 1, 12, 0, 0);
% B = datetime(1991, 12, 31, 12, 0, 0);
% B_next = datetime(1992, 1, 31, 12, 0, 0);
% 
% data = xlsread('mean-daily-temperature-fisher-ri.xlsx');
% % data = data(14 : 1474,2);
% num_data = size(data, 1);
% num_data_train = round(num_data * 0.8);
% num_data_test = num_data - num_data_train;
% 
% data_train = data(1 : num_data_train);
% data_test = data(num_data_train + 1 : num_data);
% 
% date = linspace(A, B, num_data)';
% date_train = date(1 : num_data_train);
% date_test = date(num_data_train + 1 : num_data);
% % date_dense = linspace(A_before,B_next,num_data*48)';
% date_num = datenum(date);
% date_train_num = date_num(1 : num_data_train);
% date_test_num = date_num(num_data_train + 1 : num_data);
% % date_dense_num = datenum(date_dense);

%% Construct and initialize the model
% Initialize likelihood function, covariance function and base function for mean function.
% Set the prior for the parameters of covariance function and mean function
prior_lengthScale = prior_gaussian('s2', 1000);
prior_magnSigma2 = prior_gaussian('s2', 1000);
prior_sigma2 = prior_gaussian('s2', 1000);
prior_lengthScale_sexp = prior_t('s2', 100^2, 'nu', 3);
prior_period = prior_t('s2', 100, 'nu', 4);
lik = lik_gaussian('sigma2', 1);
% gpcf = gpcf_sexp('lengthScale', 0.5, 'magnSigma2', 1, 'lengthScale_prior', prior_lengthScale, 'magnSigma2_prior', prior_magnSigma2);
% gpcf2 = gpcf_exp('lengthScale', 0.5, 'magnSigma2', 1, 'lengthScale_prior', prior_lengthScale, 'magnSigma2_prior', prior_magnSigma2);
gpcf3 = gpcf_periodic('lengthScale', 1.3, 'magnSigma2', 5, 'lengthScale_prior', prior_lengthScale, 'magnSigma2_prior', prior_magnSigma2, 'period', 12, 'lengthScale_sexp', 10, 'decay', 0, 'lengthScale_sexp_prior', prior_lengthScale_sexp, 'period_prior', prior_period);
% gpmf = 'none';
% gpmf = gpmf_linear('prior_mean',0,'prior_cov',1);

lik = lik_gaussian(lik, 'sigma2_prior', prior_sigma2);
% Initialize gp structure
% gp = gp_set('lik', lik, 'cf', gpcf);%, 'meanf', gpmf);
% gp2 = gp_set('lik', lik, 'cf', gpcf2);%, 'meanf', gpmf);
gp3 = gp_set('lik', lik, 'cf', gpcf3);%, 'meanf', gpmf);

%% Update hyperparameters
opt=optimset('TolX',1e-4,'TolFun',1e-4); 
% gp=gp_optim(gp,date_train_num,data_train,'opt',opt);
% gp2=gp_optim(gp2,date_train_num,data_train,'opt',opt);
gp3=gp_optim(gp3,date_train_num,data_train,'opt',opt);

%% Monte Carlo sampling
fprintf('Sampling of hyperparameters with HMC.\n')
num_rep = 5;
num_thin = round(0.2 * num_rep);

% [rgp,g,opt]=gp_mc(gp, date_train_num, data_train, 'nsamples', num_rep + num_thin - 1, 'display',num_thin);
% [rgp,g,opt]=gp_mc(gp2, date_train_num, data_train, 'nsamples', num_rep + num_thin - 1, 'display',num_thin);
[rgp,g,opt]=gp_mc(gp3, date_train_num, data_train, 'nsamples', num_rep + num_thin - 1, 'display',num_thin);

% Here we remove the burn-in in the samples
rgp=thin(rgp, num_thin);

%% Generate replicated data by parameters
for n = 1 : num_rep
    magnitude = rgp.cf{1}.magnSigma2(n);
    lengthScale = rgp.cf{1}.lengthScale(n);
%     lengthScale_sexp = rgp.cf{1}.lengthScale_sexp(n);
%     lengthScale_sexp = 500000000;
    period = rgp.cf{1}.period(n);
    sigma = rgp.lik.sigma2(n);

%     cov_rep = cov_sexp(date_num, date_num, magnitude, lengthScale);
    cov_rep = cov_periodic(date_num, date_num, magnitude, period, lengthScale);
    
    data_ft_rep(:, n) = mvnrnd(zeros(1, num_data), cov_rep);
    data_yt_rep(:, n) = mvnrnd(data_ft_rep(:,  n), sigma * eye(num_data));
    
    plot(data_yt_rep(:, n), '.')
    
    hold on
    plot(data, '.')
    
    vline(num_data_train, 'g', 'boundary between train and test')
    hline(0,'--')
    
    legend('sampled manually', 'original data')
    title('p(y|parameters) (periodic kernel) (set value)')
%     title('p(y|parameters) (SE kernel) (set value)')
end

%% Calculating posterior p-value for replicated data generated by parameters
% gradient_obs = gradient(data_train);
% norm_gradient_obs = norm(gradient_obs);
% % gradient_ft_rep = zeros(num_data_test, num_rep);
% % norm_gradient_ft_rep = zeros(1, num_rep);
% gradient_yt_rep = zeros(size(data_yt_rep));
% norm_gradient_yt_rep = zeros(1, num_rep);
% 
% for i = 1 : num_rep
% %     gradient_ft_rep(:, i) = gradient(data_ft_rep(:, i));
% %     norm_gradient_ft_rep(i) = norm(gradient_ft_rep(:, i));
%     gradient_yt_rep(:, i) = gradient(data_yt_rep(:, i));
%     norm_gradient_yt_rep(i) = norm(gradient_yt_rep(:, i));
% end
% 
% figure
% histogram(gradient_yt_rep(5, :), 20);
% hold on
% plot([gradient_obs(5) gradient_obs(5)], get(gca, 'ylim'));
% 
% % p_value_ft_local = sum(gradient_ft_rep > gradient_real * ones(1, num_rep), 2) / num_rep;
% % p_value_ft_local_twoside = 2 * min(p_value_ft_local, ones(size(data, 1), 1) - p_value_ft_local);
% % p_value_ft_local_binary = (p_value_ft_local_twoside < 0.1 * ones(size(data, 1), 1));
% 
% % p_value_ft_norm = sum(norm_gradient_ft_rep > norm_gradient_real * ones(1, num_rep))/num_rep;
% % p_value_ft_norm_twoside = 2 * min(p_value_ft_norm, 1 - p_value_ft_norm);
% 
% p_value_yt_norm = sum(norm_gradient_yt_rep > norm_gradient_obs * ones(1, num_rep))/num_rep;
% p_value_yt_norm_twoside = 2 * min(p_value_yt_norm, 1 - p_value_yt_norm);

%% Generate replicated data by parameters and observations
for n = 1 : num_rep
    magnitude = rgp.cf{1}.magnSigma2(n);
    lengthScale = rgp.cf{1}.lengthScale(n);
%     lengthScale_sexp = rgp.cf{1}.lengthScale_sexp(n);
%     lengthScale_sexp = 500000000;
    period = rgp.cf{1}.period(n);
    sigma = rgp.lik.sigma2(n);

%     magnitude = 10;
%     lengthScale = 50;
%     sigma = 1;
    
    lik = lik_gaussian('sigma2', sigma);
%     gpcf = gpcf_sexp('lengthScale', lengthScale, 'magnSigma2', magnitude);
%     gpcf2 = gpcf_exp('lengthScale', lengthScale, 'magnSigma2', magnitude);
    gpcf3 = gpcf_periodic('lengthScale', lengthScale, 'magnSigma2', magnitude,'period', period, 'lengthScale_sexp', lengthScale_sexp, 'decay', 1);

    % Initialize gp structure
%     gp_temp = gp_set('lik', lik, 'cf', gpcf);
%     gp_temp2 = gp_set('lik', lik, 'cf', gpcf2);
    gp_temp3 = gp_set('lik', lik, 'cf', gpcf3);
    [data_ft_rep(:,n), data_yt_rep(:,n)] = gp_rnd(gp_temp3, date_train_num, data_train, date_num);
    
    plot(data_yt_rep(:, n),'.')
    
%     cov_rr = cov_sexp(date_num, date_num, magnitude, lengthScale);
%     cov_ro = cov_sexp(date_num, date_train_num, magnitude, lengthScale);
%     cov_oo = cov_sexp(date_train_num, date_train_num, magnitude, lengthScale);

    cov_rr = cov_periodic(date_num, date_num, magnitude, period, lengthScale, lengthScale_sexp);
    cov_ro = cov_periodic(date_num, date_train_num, magnitude, period, lengthScale, lengthScale_sexp);
    cov_oo = cov_periodic(date_train_num, date_train_num, magnitude, period, lengthScale, lengthScale_sexp);
    
    mean_rep = cov_ro * ((cov_oo + sigma * eye(num_data_train)) \ data_train);
    cov_rep = cov_rr - cov_ro * ((cov_oo + sigma * eye(num_data_train)) \ cov_ro');
    
%     [mean_rep, diag_rep] = gp_pred(gp_temp3, date_train_num, data_train, date_num);
%     cov_rep = diag(diag_rep);
%     figure
%     plot(mean_rep1)
%     hold on
%     plot(mean_rep)
%     figure
%     plot(diag(cov_rep1))
%     hold on
%     plot(diag_rep)
    
    data_ft_rep(:, n) = mvnrnd(mean_rep, cov_rep);
    data_yt_rep(:, n) = mvnrnd(data_ft_rep(:, n), sigma * eye(num_data));
    hold on
    plot(data_yt_rep(:, n), '.')
    
    hold on
    plot(data, '.')
    
    hold on
    plot(mean_rep)
    
    hold on
    plotshaded(1:num_data, [mean_rep - 2 * diag(cov_rep), mean_rep + 2 * diag(cov_rep)]', 'r')
    
    vline(num_data_train, 'g', 'boundary between train and test')
    hline(0,'--')
    
    legend('sampled by gpstuff', 'sampled manually', 'original data', 'posterior mean', 'confidence interval')
    title('p(y|parameters, observation) (periodic kernel) (set value)')
%     title('p(y|parameters, observation) (SE kernel) (set value)')
end

%% Calculating posterior p-value for replicated data generated by parameters and observation
gradient_obs = gradient(data_test);
norm_gradient_obs = norm(gradient_obs);
% gradient_ft_rep = zeros(num_data_test, num_rep);
% norm_gradient_ft_rep = zeros(1, num_rep);
gradient_yt_rep = zeros(size(data_yt_rep));
norm_gradient_yt_rep = zeros(1, num_rep);

for i = 1 : num_rep
%     gradient_ft_rep(:, i) = gradient(data_ft_rep(:, i));
%     norm_gradient_ft_rep(i) = norm(gradient_ft_rep(:, i));
    gradient_yt_rep(:, i) = gradient(data_yt_rep(:, i));
    norm_gradient_yt_rep(i) = norm(gradient_yt_rep(:, i));
end

figure
histogram(gradient_yt_rep(5, :), 20);
hold on
plot([gradient_obs(5) gradient_obs(5)], get(gca, 'ylim'));

% p_value_ft_local = sum(gradient_ft_rep > gradient_real * ones(1, num_rep), 2) / num_rep;
% p_value_ft_local_twoside = 2 * min(p_value_ft_local, ones(size(data, 1), 1) - p_value_ft_local);
% p_value_ft_local_binary = (p_value_ft_local_twoside < 0.1 * ones(size(data, 1), 1));

% p_value_ft_norm = sum(norm_gradient_ft_rep > norm_gradient_real * ones(1, num_rep))/num_rep;
% p_value_ft_norm_twoside = 2 * min(p_value_ft_norm, 1 - p_value_ft_norm);

p_value_yt_norm = sum(norm_gradient_yt_rep > norm_gradient_obs * ones(1, num_rep))/num_rep;
p_value_yt_norm_twoside = 2 * min(p_value_yt_norm, 1 - p_value_yt_norm);
