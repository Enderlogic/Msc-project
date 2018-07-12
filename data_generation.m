function [x,data] = data_generation(num_data, length, random_seed, hyp_data, data_cov)
% choose the kernel function of data from:
% Matern kernel (Matern)
% SE kernel (SE)
% Periodic kernel (Periodic)

%choose fixed random data or variational random data ('on' or 'off')

%set the length of data by 'length'
x = (0 : length/(num_data - 1) : length)';
if strcmp(data_cov{1}, 'Periodic')
    covfunc_data = {@covPeriodic};
    hyp.cov = log([hyp_data.lengthScale; hyp_data.period; sqrt(hyp_data.magnitude)]);
elseif strcmp(data_cov{1}, 'SE')
    covfunc_data = {@covSEiso};
    hyp.cov = log([hyp_data.lengthScale; sqrt(hyp_data.magnitude)]);
elseif strcmp(data_cov{1}, 'Matern')
    covfunc_data = {'covMaterniso', data_cov{2}};
    hyp.cov = log([hyp_data.lengthScale, sqrt(hyp_data.magnitude)]);
elseif strcmp(data_cov{1}, 'RQ')
    covfunc_data = {@covRQiso};
    hyp.cov = log([hyp_data.lengthScale, sqrt(hyp_data.magnitude), hyp_data.alpha]);
else
    error('The type of kernel function is invalid')
end
K = feval(covfunc_data{:}, hyp.cov, x);
K = (K + K') / 2;
mu = zeros(1,num_data);
if strcmp(random_seed, 'on')
    rng(19990806)
    data = mvnrnd(mu, K);
    rng(19990806)
    data = mvnrnd(data, hyp_data.sigma * eye(num_data))';
elseif strcmp(random_seed, 'off')
    data = mvnrnd(mu, K);
    data = mvnrnd(data, hyp_data.sigma * eye(num_data))';
else
    error('The state of random seed is invalid')
end
end

