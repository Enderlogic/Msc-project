function [p_value] = p_value_calculator(data_test, data_rep, criteria, model_cov, sample_location, sample_hyp, x_train, data_train, x_test)
    % criteria: number_of_zero; norm_of_gradient; chi_square; mmd
    % [p_value] = p_value_calculator(data_test, data_rep, criteria, sample_type, model_cov, sample_hyp, x_train, data_train, x_test)
    if strcmp(criteria, 'mmd')
        if isfield(data_rep, 'para_obs')
            fail = 0;
            for i = 1 : size(data_rep.para_obs, 2)
                if strcmp(sample_location, 'test')
                    fail = fail + kmd([data_test; data_rep.para_obs(:, i)], [ones(size(data_test, 1), 1); -ones(size(data_rep.para_obs(:, i), 1), 1)]);
                elseif strcmp(sample_location, 'train')
                    fail = fail + kmd([data_train; data_rep.para_obs(:, i)], [ones(size(data_train, 1), 1); -ones(size(data_rep.para_obs(:, i), 1), 1)]);
                else
                    error('The type of sample location is invalid!')
                end
            end
            p_value = 1 - fail / size(data_rep.para_obs, 2);
            disp(['p-value for MMD test is: ', num2str(p_value)]);
        else
            error('MMD test could only be applied on p(y_rep|para, y_obs). Please use para&obs or both as sample_type')
        end
        elseif strcmp(criteria, 'chi_square')
            p_value_temp = zeros(size(data_rep, 2), 1);
            if isfield(data_rep, 'para')
                cri_data.para = zeros(size(data_rep, 2), 1);
%                 cri_rep.para = zeros(size(data_rep, 2), 1);
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
                    error('The type of kernel function for MLE solution is invalid')
                end
                for i = 1 : size(data_rep.para, 2)
                    if strcmp(sample_location, 'test')
                        K = feval(covfunc{:}, hyp.cov(:,i), x_test);
                        K = (K + K') / 2;
                        K = K + sample_hyp.lik.sigma2(i) * eye(size(x_test, 1));
                        cri_data.para(i) = data_test' / K * data_test;
                        p_value_temp(i) = 1 - chi2cdf(cri_data.para(i), size(x_test, 1));
                    elseif strcmp(sample_location, 'train')
                        K = feval(covfunc{:}, hyp.cov(:,i), x_train);
                        K = (K + K') / 2;
                        K = K + sample_hyp.lik.sigma2(i) * eye(size(x_train, 1));
                        cri_data.para(i) = data_train' / K * data_train;
                        p_value_temp(i) = 1 - chi2cdf(cri_data.para(i), size(x_train, 1));
                    else
                        error('The type of sample location is invalid!')
                    end
%                     cri_rep.para(i) = data_rep.para(:, i)' / K * data_rep.para(:, i);
                end
                p_value.para = mean(p_value_temp);
                disp(['p-value for p(y^rep|para) is: ', num2str(p_value.para)])
            end
            if isfield(data_rep, 'para_obs')
                for i = 1 : size(data_rep.para_obs, 2)
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
                        [~, ~, ~, Eyt, ~, Covyt] = gp_pred(gp, x_train, data_train, x_test);
                        cri_data.para_obs(i) = (data_test - Eyt)' / Covyt * (data_test - Eyt);
                        p_value_temp(i) = 1 - chi2cdf(cri_data.para_obs(i), size(data_test, 1));
                    elseif strcmp(sample_location, 'train')
                        [~, ~, ~, Eyt, ~, Covyt] = gp_pred(gp, x_train, data_train, x_train);
                        cri_data.para_obs(i) = (data_train - Eyt)' / Covyt * (data_train - Eyt);
                        p_value_temp(i) = 1 - chi2cdf(cri_data.para_obs(i), size(data_train, 1));
                    else
                        error('The type of sample location is invalid!')
                    end
%                     cri_rep.para_obs(i) = (data_rep.para_obs(:, i) - Eyt)' / Covyt * (data_rep.para_obs(:, i) - Eyt);
                end
                p_value.para_obs = mean(p_value_temp);
                disp(['p-value for p(y^rep|para, y^obs) is: ', num2str(p_value.para_obs)])
            end
    else
        if strcmp(criteria, 'norm_of_gradient')
            % compute gradient with momentum
            beta = 0.9;
            if strcmp(sample_location, 'test')
                par1 = (1 - beta) * fliplr(beta .^ (0 : size(data_test, 1) - 1))';
                par2 = beta .^ (1 - size(data_test, 1) : 0)';
                par3 = 1 - beta .^ (1 : size(data_test, 1))';

                data_test_rescale = par1 .* data_test;
                data_test_rescale = cumsum(data_test_rescale);
                data_test_rescale = data_test_rescale .* par2;
                data_test_rescale = data_test_rescale ./ par3;

                if isfield(data_rep, 'para')
                    data_test_rep_rescale = par1 .* data_rep.para;
                    data_test_rep_rescale = cumsum(data_test_rep_rescale);
                    data_test_rep_rescale = data_test_rep_rescale .* par2;
                    data_test_rep_rescale = data_test_rep_rescale ./ par3;

                    cri_data.para = norm(gradient(data_test_rescale));

                    [~, ygrad] = gradient(data_test_rep_rescale);
                    cri_rep.para = vecnorm(ygrad);
                end
                if isfield(data_rep, 'para_obs')
                    data_test_rep_rescale = par1 .* data_rep.para_obs;
                    data_test_rep_rescale = cumsum(data_test_rep_rescale);
                    data_test_rep_rescale = data_test_rep_rescale .* par2;
                    data_test_rep_rescale = data_test_rep_rescale ./ par3;

                    cri_data.para_obs = norm(gradient(data_test_rescale));
                    [~, ygrad] = gradient(data_test_rep_rescale);

                    cri_rep.para_obs = vecnorm(ygrad);
                end
            elseif strcmp(sample_location, 'train')
                par1 = (1 - beta) * fliplr(beta .^ (0 : size(data_train, 1) - 1))';
                par2 = beta .^ (1 - size(data_train, 1) : 0)';
                par3 = 1 - beta .^ (1 : size(data_train, 1))';

                data_train_rescale = par1 .* data_train;
                data_train_rescale = cumsum(data_train_rescale);
                data_train_rescale = data_train_rescale .* par2;
                data_train_rescale = data_train_rescale ./ par3;

                if isfield(data_rep, 'para')
                    data_train_rep_rescale = par1 .* data_rep.para;
                    data_train_rep_rescale = cumsum(data_train_rep_rescale);
                    data_train_rep_rescale = data_train_rep_rescale .* par2;
                    data_train_rep_rescale = data_train_rep_rescale ./ par3;

                    cri_data.para = norm(gradient(data_train_rescale));

                    [~, ygrad] = gradient(data_train_rep_rescale);
                    cri_rep.para = vecnorm(ygrad);
                end
                if isfield(data_rep, 'para_obs')
                    data_train_rep_rescale = par1 .* data_rep.para_obs;
                    data_train_rep_rescale = cumsum(data_train_rep_rescale);
                    data_train_rep_rescale = data_train_rep_rescale .* par2;
                    data_train_rep_rescale = data_train_rep_rescale ./ par3;

                    cri_data.para_obs = norm(gradient(data_train_rescale));
                    [~, ygrad] = gradient(data_train_rep_rescale);

                    cri_rep.para_obs = vecnorm(ygrad);
                end
            else
                error('The type of sample location is invalid!')
            end
        elseif strcmp(criteria, 'number_of_zero')
            if isfield(data_rep, 'para')
                if strcmp(sample_location, 'test')
                    cri_data.para = sum(diff(data_test > mean(data_test)) ~= 0);
                    cri_rep.para = sum(diff(data_rep.para > mean(data_rep.para)) ~= 0);
                elseif strcmp(sample_location, 'train')
                    cri_data.para = sum(diff(data_train > mean(data_train)) ~= 0);
                    cri_rep.para = sum(diff(data_rep.para > mean(data_rep.para)) ~= 0);
                else
                    error('The type of sample location is invalid!')
                end
            end
            if isfield(data_rep, 'para_obs')
                if strcmp(sample_location, 'test')
                    cri_data.para_obs = sum(diff(data_test > mean(data_test)) ~= 0);
                    cri_rep.para_obs = sum(diff(data_rep.para_obs > mean(data_rep.para_obs)) ~= 0);
                elseif strcmp(sample_location, 'train')
                    cri_data.para_obs = sum(diff(data_train > mean(data_train)) ~= 0);
                    cri_rep.para_obs = sum(diff(data_rep.para_obs > mean(data_rep.para_obs)) ~= 0);
                else
                    error('The type of sample location is invalid!')
                end
            end
        elseif strcmp(criteria, 'frequency')
            fft_data_test = fft(data_test);
            tem = abs(fft_data_test / size(data_test, 1));
            fre_data_test = tem(1 : size(data_test, 1) / 2 + 1);
            fre_data_test(2 : end - 1) = 2 * fre_data_test(2 : end - 1);
            f_data_test = 1 / (x_test(2) - x_test(1)) * (0 : size(data_test, 1) / 2) / size(data_test, 1);
        else
            error('The criteria is invalid!')
        end
        if isfield(data_rep, 'para')
            p_value.para = 2 * min(sum(cri_rep.para > cri_data.para) / size(data_rep.para, 2), 1 - sum(cri_rep.para > cri_data.para) / size(data_rep.para, 2));
            disp(['p-value for p(y^rep|para) is: ', num2str(p_value.para)])
        end
        if isfield(data_rep, 'para_obs')
            p_value.para_obs = 2 * min(sum(cri_rep.para_obs > cri_data.para_obs) / size(data_rep.para_obs, 2), 1 - sum(cri_rep.para_obs > cri_data.para_obs) / size(data_rep.para_obs, 2));
            disp(['p-value for p(y^rep|para, y^obs) is: ', num2str(p_value.para_obs)])
        end
    end
end