function [data_rep] = sample_data(sample_hyp, model_cov, sample_location, sample_type, x_train, data_train, x_test)
    % num_rep: the number of replicated data
    % sample_hyp: the posterior samples of hyperparameters
    % model_cov: the covariance function usedted in the fit model

    % [yt] = sample_data(num_rep, sample_hyp, model_cov, sample_type, x_train, data_train, x_test);
    num_rep = size(sample_hyp.lik.sigma2, 1);
    if strcmp(sample_location, 'test')
        num_data = size(x_test, 1);
    elseif strcmp(sample_location, 'train')
        num_data = size(x_train, 1);
    else
        error('The type of sample location is invalid!')
    end
    if strcmp(sample_type, 'para') || strcmp(sample_type, 'both')
        data_rep.para = zeros(num_data, num_rep);
        if strcmp(model_cov{1}, 'sum') || strcmp(model_cov{1}, 'prod')
            hyp.cov = [];
            if strcmp(model_cov{1}, 'sum')
                covfunc = {'covSum', {}};
            else
                covfunc = {'covProd', {}};
            end
            for i = 1 : size(model_cov, 2) - 1
                switch model_cov{i + 1}
                    case 'SE'
                        covfunc{2} = [covfunc{2}, 'covSEiso'];
                        hyp.cov = [hyp.cov; log([sample_hyp.cf{1, 1}.cf{1, i}.lengthScale'; sqrt(sample_hyp.cf{1, 1}.cf{1, i}.magnSigma2')])];
                    case 'LIN'
                        covfunc{2} = [covfunc{2}, 'covLINiso'];
                        hyp.cov = [hyp.cov; log(sample_hyp.cf{1, 1}.cf{1, i}.coeffSigma2')];
                    case 'Periodic'
                        covfunc{2} = [covfunc{2}, 'covPeriodic'];
                        hyp.cov = [hyp.cov; log([sample_hyp.cf{1, 1}.cf{1, i}.lengthScale'; sample_hyp.cf{1, 1}.cf{1, i}.period'; sqrt(sample_hyp.cf{1, 1}.cf{1, i}.magnSigma2')])];
                    otherwise
                        error('The type of covariance funciton is invalid for composition!')
                end
            end
        elseif strcmp(model_cov{1}, 'Periodic')
            covfunc = {@covPeriodic};
            hyp.cov = log([sample_hyp.cf{1}.lengthScale, sample_hyp.cf{1}.period, sqrt(sample_hyp.cf{1}.magnSigma2)])';
        elseif strcmp(model_cov{1}, 'SE')
            covfunc = {@covSEiso};
            hyp.cov = log([sample_hyp.cf{1}.lengthScale, sqrt(sample_hyp.cf{1}.magnSigma2)])';
        elseif strcmp(model_cov{1}, 'Matern')
            covfunc = {@covMaterniso, str2double(model_cov{2})};
            hyp.cov = log([sample_hyp.cf{1}.lengthScale, sqrt(sample_hyp.cf{1}.magnSigma2)])';
        elseif strcmp(model_cov{1}, 'RQ')
            covfunc = {@covRQiso};
            hyp.cov = log([sample_hyp.cf{1}.lengthScale, sqrt(sample_hyp.cf{1}.magnSigma2), sample_hyp.cf{1}.alpha])';
        elseif strcmp(model_cov{1}, 'Exp')
            covfunc = {@covMaterniso, 1};
            hyp.cov = log([sample_hyp.cf{1}.lengthScale, sqrt(sample_hyp.cf{1}.magnSigma2)])';
        else
            error('The type of covariance function is invalid')
        end
        mu = mean(data_train) * ones(num_data, 1); % assume using zero mean function, more non-zero mean functions will be finished later
        for i = 1 : num_rep
            if strcmp(sample_location, 'test')
                K = feval(covfunc{:}, hyp.cov(:,i), x_test);
            elseif strcmp(sample_location, 'train')
                K = feval(covfunc{:}, hyp.cov(:,i), x_train);
            else
                error('The type of sample location is invalid!')
            end
            K = (K + K') / 2;
            data_rep.para(:, i) = mvnrnd(mu, K + sample_hyp.lik.sigma2(i) * eye(num_data))';
        end
    end
    if strcmp(sample_type, 'para&obs') || strcmp(sample_type, 'both')
        data_rep.para_obs = zeros(num_data, num_rep);
        for i = 1 : num_rep
            lik = lik_gaussian('sigma2', sample_hyp.lik.sigma2(i));
            if strcmp(model_cov{1}, 'sum') || strcmp(model_cov{1}, 'prod')
                gpcf_com = {};
                for j = 1 : size(model_cov, 2) - 1
                    switch model_cov{j + 1}
                        case 'SE'
                            gpcf_com{end + 1} = gpcf_sexp('lengthScale', sample_hyp.cf{1, 1}.cf{1, j}.lengthScale(i), 'magnSigma2', sample_hyp.cf{1, 1}.cf{1, j}.magnSigma2(i));
                        case 'LIN'
                            gpcf_com{end + 1} = gpcf_linear('coeffSigma2', sample_hyp.cf{1, 1}.cf{1, j}.coeffSigma2(i));
                        case 'Periodic'
                            gpcf_com{end + 1} = gpcf_periodic('lengthScale', sample_hyp.cf{1, 1}.cf{1, j}.lengthScale(i), ...
                                'magnSigma2', sample_hyp.cf{1, 1}.cf{1, j}.magnSigma2(i), 'period', sample_hyp.cf{1, 1}.cf{1, j}.period(i), 'decay', 0);
                        otherwise
                            error('The type of covariance funciton is invalid for composition!')
                    end
                end
                if strcmp(model_cov{1}, 'sum')
                    gpcf = gpcf_sum('cf', gpcf_com);
                else
                    gpcf = gpcf_prod('cf', gpcf_com);
                end
            elseif strcmp(model_cov{1}, 'Periodic')
                gpcf = gpcf_periodic('lengthScale', sample_hyp.cf{1}.lengthScale(i), ...
                    'magnSigma2', sample_hyp.cf{1}.magnSigma2(i), 'period', sample_hyp.cf{1}.period(i), 'decay', 0);
            elseif strcmp(model_cov{1}, 'SE')
                gpcf = gpcf_sexp('lengthScale', sample_hyp.cf{1}.lengthScale(i), 'magnSigma2', sample_hyp.cf{1}.magnSigma2(i));
            elseif strcmp(model_cov{1}, 'Matern')
                if strcmp(model_cov{2}, '3')
                    gpcf = gpcf_matern32('lengthScale', sample_hyp.cf{1}.lengthScale(i), 'magnSigma2', sample_hyp.cf{1}.magnSigma2(i));
                elseif strcmp(model_cov{2}, '5')
                    gpcf = gpcf_matern52('lengthScale', sample_hyp.cf{1}.lengthScale(i), 'magnSigma2', sample_hyp.cf{1}.magnSigma2(i));
                else
                    error('Only support Matern 3/2 and Matern 5/2 covariance function!')
                end
            elseif strcmp(model_cov{1}, 'RQ')
                gpcf = gpcf_rq('lengthScale', sample_hyp.cf{1}.lengthScale(i), 'magnSigma2', sample_hyp.cf{1}.magnSigma2(i), 'alpha', sample_hyp.cf{1}.alpha(i));
            elseif strcmp(model_cov{1}, 'Exp')
                gpcf = gpcf_exp('lengthScale', sample_hyp.cf{1}.lengthScale(i), 'magnSigma2', sample_hyp.cf{1}.magnSigma2(i));
            else
                error('The type of covariance function is invalid')
            end
            gp = gp_set('lik', lik, 'cf', gpcf);
            if strcmp(sample_location, 'test')
                [~, data_rep.para_obs(:,i)] = gp_rnd(gp, x_train, data_train, x_test);
            elseif strcmp(sample_location, 'train')
                [~, data_rep.para_obs(:,i)] = gp_rnd(gp, x_train, data_train, x_train);
            else
                error('The type of sample location is invalid!')
            end
        end
    end
end