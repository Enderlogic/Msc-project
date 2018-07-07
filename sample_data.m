function [data_rep] = sample_data(num_rep, sample_hyp, model_cov, sample_type, x_train, data_train, x_test)
% num_rep: the number of replicated data
% sample_hyp: the posterior samples of hyperparameters
% model_cov: the covariance function usedted in the fit model

% [yt] = sample_data(num_rep, sample_hyp, model_cov, sample_type, x_train, data_train, x_test);
num_data = size(x_test, 1);
data_rep = zeros(num_data, num_rep);
if strcmp(sample_type, 'para')
    mu = zeros(1, num_data); % assume using zero mean function, more non-zero mean functions will be finished later
    if strcmp(model_cov{1}, 'Periodic')
        covfunc = {@covPeriodic};
        hyp.cov = log([sample_hyp.cf{1}.lengthScale, sample_hyp.cf{1}.period, sqrt(sample_hyp.cf{1}.magnSigma2)])';
    elseif strcmp(model_cov{1}, 'SE')
        covfunc = {@covSEiso};
        hyp.cov = log([sample_hyp.cf{1}.lengthScale, sqrt(sample_hyp.cf{1}.magnSigma2)])';
    elseif strcmp(model_cov{1}, 'Matern')
        covfunc = {@covMaterniso, model_cov{2}};
        hyp.cov = log([sample_hyp.cf{1}.lengthScale, sqrt(sample_hyp.cf{1}.magnSigma2)])';
    else
        error('The type of covariance function is invalid')
    end
    for i = 1 : num_rep
        K = feval(covfunc{:}, hyp.cov(:,i), x_test);
        K = (K + K') / 2;
        data_rep(:, i) = mvnrnd(mu, K)';
        data_rep(:, i) = mvnrnd(data_rep(:, i), sample_hyp.lik.sigma2(i) * eye(num_data));
    end
elseif strcmp(sample_type, 'para&obs')
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
        [~, data_rep(:,i)] = gp_rnd(gp, x_train, data_train, x_test);
    end
else
    error('The type of sampling is invalid!')
end
end
