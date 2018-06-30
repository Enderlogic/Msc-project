function [data_rep] = sample_data(num_rep, sample_hyp, model_cov, varargin)
% num_rep: the number of replicated data
% sample_hyp: the posterior samples of hyperparameters
% model_cov: the covariance function usedted in the fit model

% p(yt|'parameters'): [yt] = sample_data(num_rep, sample_hyp, model_cov, x);
% p(yt|'parameters', observation): [yt] = sample_data(num_rep, sample_hyp, model_cov, x_train, data_train, x);
if size(varargin, 2) == 1
    num_data = size(varargin{1}, 1);
    data_rep = zeros(num_data, num_rep);
    mu = zeros(1, num_data); % assume using zero mean function, more non-zero mean functions will be finished later
    if strcmp(model_cov{1}, 'Periodic')
        covfunc = {@covPeriodic};
        hyp.cov = log([sample_hyp.cf{1}.lengthScale, sample_hyp.cf{1}.period, sqrt(sample_hyp.cf{1}.magnSigma2)])';
%         cov_func = @(X, hyp) reshape(hyp.cf{1, 1}.magnSigma2, [1 1 num_rep]) .* exp(- 2 * sin(pi * repmat((X * ones(1, size(X, 1)) - ones(size(X, 1), 1) * X').^2, [1 1 num_rep]) ./ reshape(hyp.cf{1}.period, [1 1 num_rep])).^2 ./ reshape(hyp.cf{1}.lengthScale.^2, [1 1 num_rep]));
    elseif strcmp(model_cov{1}, 'SE')
        covfunc = {@covSEiso};
        hyp.cov = log([sample_hyp.cf{1}.lengthScale, sqrt(sample_hyp.cf{1}.magnSigma2)])';
%         cov_func = @(X, hyp) reshape(hyp.cf{1, 1}.magnSigma2, [1 1 num_rep]) .* exp(- repmat((X * ones(1, size(X, 1)) - ones(size(X, 1), 1) * X').^2, [1 1 num_rep]) / 2 ./ reshape(hyp.cf{1}.lengthScale.^2, [1 1 num_rep]));
    elseif strcmp(model_cov{1}, 'Matern')
        covfunc = {@covMaterniso, model_cov{2}};
        hyp.cov = log([sample_hyp.cf{1}.lengthScale, sqrt(sample_hyp.cf{1}.magnSigma2)])';
%         r = @(X, hyp) abs(repmat(X * ones(1, size(X, 1)) - ones(size(X, 1), 1) * X', [1 1 num_rep]) ./ reshape(hyp.cf{1, 1}.lengthScale, [1 1 num_rep]));
%         cov_func = @(X, hyp) reshape(hyp.cf{1, 1}.magnSigma2, [1 1 num_rep]) .* (1 + sqrt(3) * r(X, hyp)) .* exp(- sqrt(3) * r(X, hyp));
    else
        error('The type of covariance function is invalid')
    end
%     K = cov_func(varargin{1}, sample_hyp);
    for i = 1 : num_rep
        K = feval(covfunc{:}, hyp.cov(:,i), varargin{1});
        K = (K + K') / 2;
        data_rep(:, i) = mvnrnd(mu, K)';
        data_rep(:, i) = mvnrnd(data_rep(:, i), sample_hyp.lik.sigma2(i) * eye(num_data));
    end
elseif size(varargin, 2) == 3
    num_data_rep = size(varargin{3}, 1);
    data_ft_rep = zeros(num_data_rep, num_rep);
    data_rep = zeros(num_data_rep, num_rep);
    for i = 1 : num_rep
        lik = lik_gaussian('sigma2', sample_hyp.lik.sigma2(i));
        if strcmp(model_cov{1}, 'Periodic')
            gpcf = gpcf_periodic('lengthScale', sample_hyp.cf{1}.lengthScale(i), ...
                'magnSigma2', sample_hyp.cf{1}.magnSigma2(i), 'period', sample_hyp.cf{1}.period(i), 'decay', 0);
        elseif strcmp(model_cov{1}, 'SE')
            gpcf = gpcf_sexp('lengthScale', sample_hyp.cf{1}.lengthScale(i), 'magnSigma2', sample_hyp.cf{1}.magnSigma2(i));
        elseif strcmp(model_cov{1}, 'Matern')
            if model_cov{2} == 1
                gpcf = gpcf_matern32('lengthScale', sample_hyp.cf{1}.lengthScale(i), 'magnSigma2', sample_hyp.cf{1}.magnSigma2(i));
            elseif model_cov{2} == 2
                gpcf = gpcf_matern52('lengthScale', sample_hyp.cf{1}.lengthScale(i), 'magnSigma2', sample_hyp.cf{1}.magnSigma2(i));
            else
                error('Only support Matern 3/2 and Matern 5/2 covariance function!')
            end
        else
            error('The type of covariance function is invalid')
        end
        gp = gp_set('lik', lik, 'cf', gpcf);
        [data_ft_rep(:,i), data_rep(:,i)] = gp_rnd(gp, varargin{1}, varargin{2}, varargin{3});
    end
else
    error('The number of input argument is wrong!')
end
end
