function [p_value] = p_value_calculator(data_test, data_test_rep, criteria, sample_type, model_cov, sample_hyp, x_train, data_train, x_test)
% criteria: number_of_zero; norm_of_gradient; chi_square
% [p_value] = p_value_calculator(data_test, data_test_rep, criteria, sample_type, model_cov, sample_hyp, x_train, data_train, x_test)
if strcmp(criteria, 'number_of_mean')
    cri_data = sum(diff(data_test > mean(data_test)) ~= 0);
    cri_rep = sum(diff(data_test_rep > mean(data_test_rep)) ~= 0);
elseif strcmp(criteria, 'norm_of_gradient')
    cri_data = norm(gradient(data_test));
    [~, ygrad] = gradient(data_test_rep);
    cri_rep = vecnorm(ygrad);
elseif strcmp(criteria, 'chi_square')
    cri_data = zeros(size(data_test_rep, 2), 1);
    cri_rep = zeros(size(data_test_rep, 2), 1);
    if strcmp(sample_type, 'para')
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
        for i = 1 : size(data_test_rep, 2)
            K = feval(covfunc{:}, hyp.cov(:,i), x_test);
            K = (K + K') / 2;
            K = K + sample_hyp.lik.sigma2(i) * eye(size(x_test, 1));
%             K = diag(diag(K));
            cri_data(i) = data_test' / K * data_test;
%             cri_data2(i) = data_test' * data_test;
            cri_rep(i) = data_test_rep(:, i)' / K * data_test_rep(:, i);
%             cri_rep2(i) = data_test_rep(:, i)' * data_test_rep(:, i);
        end
    elseif strcmp(sample_type, 'para&obs')
        for i = 1 : size(data_test_rep, 2)
            if strcmp(model_cov{1}, 'Periodic')
                cov_periodic = @(X1, X2, magnitude, period, lengthScale) magnitude * exp(- 2 * sin(pi * (X1 * ones(1, size(X2, 1)) - ones(size(X1, 1), 1) * X2')/ period).^2 / (lengthScale)^2);
                
                cov_rr = cov_periodic(x_test, x_test, sample_hyp.cf{1}.magnSigma2(i), sample_hyp.cf{1}.period(i), sample_hyp.cf{1}.lengthScale(i));
                cov_ro = cov_periodic(x_test, x_train, sample_hyp.cf{1}.magnSigma2(i), sample_hyp.cf{1}.period(i), sample_hyp.cf{1}.lengthScale(i));
                cov_oo = cov_periodic(x_train, x_train, sample_hyp.cf{1}.magnSigma2(i), sample_hyp.cf{1}.period(i), sample_hyp.cf{1}.lengthScale(i));
                
                mean_rep = cov_ro * ((cov_oo + sample_hyp.lik.sigma2(i) * eye(size(x_train, 1))) \ data_train);
                cov_rep = cov_rr - cov_ro * ((cov_oo + sample_hyp.lik.sigma2(i) * eye(size(x_train, 1))) \ cov_ro');
                cov_rep = cov_rep + sample_hyp.lik.sigma2(i) * eye(size(x_test, 1));
                
%                 gpcf = gpcf_periodic('lengthScale', sample_hyp.cf{1}.lengthScale(i), 'magnSigma2', sample_hyp.cf{1}.magnSigma2(i), 'period', sample_hyp.cf{1}.period(i), 'decay', 0);
            elseif strcmp(model_cov{1}, 'SE')
                cov_sexp = @(X1, X2, magnitude, lengthScale) magnitude * exp(- (X1 * ones(1, size(X2, 1)) - ones(size(X1, 1), 1) * X2').^2/2/(lengthScale^2));
                
                cov_rr = cov_sexp(x_test, x_test, sample_hyp.cf{1}.magnSigma2(i), sample_hyp.cf{1}.lengthScale(i));
                cov_ro = cov_sexp(x_test, x_train, sample_hyp.cf{1}.magnSigma2(i), sample_hyp.cf{1}.lengthScale(i));
                cov_oo = cov_sexp(x_train, x_train, sample_hyp.cf{1}.magnSigma2(i), sample_hyp.cf{1}.lengthScale(i));
                
                mean_rep = cov_ro * ((cov_oo + sample_hyp.lik.sigma2(i) * eye(size(x_train, 1))) \ data_train);
                cov_rep = cov_rr - cov_ro * ((cov_oo + sample_hyp.lik.sigma2(i) * eye(size(x_train, 1))) \ cov_ro');
                cov_rep = cov_rep + sample_hyp.lik.sigma2(i) * eye(size(x_test, 1));
                
%                 gpcf = gpcf_sexp('lengthScale', sample_hyp.cf{1}.lengthScale(i), 'magnSigma2', sample_hyp.cf{1}.magnSigma2(i));
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
%             lik = lik_gaussian('sigma2', sample_hyp.lik.sigma2(i));
%             gp = gp_set('lik', lik, 'cf', gpcf);
%             [EFT, VARFT, LPYT, EYT, VARYT] = gp_pred(gp, x_train, data_train, x_test);
            cri_data(i) = (data_test - mean_rep)' / cov_rep * (data_test - mean_rep);
            cri_rep(i) = (data_test_rep(:, i) - mean_rep)' / cov_rep * (data_test_rep(:, i) - mean_rep);
        end
    else
        error('The type of sampling is invalid!')
    end
else
    error('The criteria is invalid!')
end
p_value = 2 * min(sum(cri_rep > cri_data) / size(data_test_rep, 2), 1 - sum(cri_rep > cri_data) / size(data_test_rep, 2));
% p_value2 = 2 * min(sum(cri_rep2 > cri_data2) / size(data_test_rep, 2), 1 - sum(cri_rep2 > cri_data2) / size(data_test_rep, 2));
end

