function [p_value] = p_value_calculator(data_test, data_test_rep, criteria, sample_type, model_cov, sample_hyp, x_train, data_train, x_test)
    % criteria: number_of_zero; norm_of_gradient; chi_square
    % [p_value] = p_value_calculator(data_test, data_test_rep, criteria, sample_type, model_cov, sample_hyp, x_train, data_train, x_test)
    if strcmp(criteria, 'mmd')
        if isfield(data_test_rep, 'para_obs')
            fail = 0;
            for i = 1 : size(data_test_rep, 2)
                fail = fail + kmd([data_test; data_test_rep(:, i)], [ones(size(data_test, 1), 1); -ones(size(data_test_rep(:, i), 1), 1)]);
            end
            p_value = 1 - fail / size(data_test_rep, 2);
        else
            error('MMD test could only be applied on p(y_rep|para, y_obs). Please use para&obs or both as sample_type')
        end
    else
        if strcmp(criteria, 'number_of_zero')
            if isfield(data_test_rep, 'para')
                cri_data.para = sum(diff(data_test > mean(data_test)) ~= 0);
                cri_rep.para = sum(diff(data_test_rep.para > mean(data_test_rep.para)) ~= 0);
            end
            if isfield(data_test_rep, 'para_obs')
                cri_data.para_obs = sum(diff(data_test > mean(data_test)) ~= 0);
                cri_rep.para_obs = sum(diff(data_test_rep.para_obs > mean(data_test_rep.para_obs)) ~= 0);
            end
        elseif strcmp(criteria, 'norm_of_gradient')
            % compute gradient with momentum
            beta = 0.99;
            par1 = (1 - beta) * fliplr(beta .^ (0 : size(data_test, 1) - 1))';
            par2 = beta .^ (1 - size(data_test, 1) : 0)';
            par3 = 1 - beta .^ (1 : size(data_test, 1))';
            
            data_test_rescale = par1 .* data_test;
            data_test_rescale = cumsum(data_test_rescale);
            data_test_rescale = data_test_rescale .* par2;
            data_test_rescale = data_test_rescale ./ par3;

            if isfield(data_test_rep, 'para')
                data_test_rep_rescale = par1 .* data_test_rep.para;
                data_test_rep_rescale = cumsum(data_test_rep_rescale);
                data_test_rep_rescale = data_test_rep_rescale .* par2;
                data_test_rep_rescale = data_test_rep_rescale ./ par3;

                cri_data.para = norm(gradient(data_test_rescale(10:end)));
                [~, ygrad] = gradient(data_test_rep_rescale(10:end,:));

                cri_rep.para = vecnorm(ygrad);
            end
            if isfield(data_test_rep, 'para_obs')
                data_test_rep_rescale = par1 .* data_test_rep.para_obs;
                data_test_rep_rescale = cumsum(data_test_rep_rescale);
                data_test_rep_rescale = data_test_rep_rescale .* par2;
                data_test_rep_rescale = data_test_rep_rescale ./ par3;

                cri_data.para_obs = norm(gradient(data_test_rescale(10:end)));
                [~, ygrad] = gradient(data_test_rep_rescale(10:end,:));

                cri_rep.para_obs = vecnorm(ygrad);
            end
        elseif strcmp(criteria, 'chi_square')
            if isfield(data_test_rep, 'para')
                cri_data.para = zeros(size(data_test_rep, 2), 1);
                cri_rep.para = zeros(size(data_test_rep, 2), 1);
                if strcmp(model_cov{1}, 'sum') || strcmp(model_cov{1}, 'prod')
                    hyp.cov = [];
                    if strcmp(model_cov{1}, 'sum')
                        covfunc = {'covSum'};
                    else
                        covfunc = {'covProd'};
                    end
                    cov = strings(1, (size(model_cov, 2) - 1));
                    for i = 1 : size(model_cov, 2) - 1
                        switch model_cov{i + 1}
                            case 'SE'
                                cov(i) = 'covSEiso';
                                hyp.cov = [hyp.cov; log([sample_hyp.cf{1, 1}.cf{1, i}.lengthScale'; sqrt(sample_hyp.cf{1, 1}.cf{1, i}.magnSigma2')])];
                            case 'LIN'
                                cov(i) = 'covLINiso';
                                hyp.cov = [hyp.cov; log(sample_hyp.cf{1, 1}.cf{1, i}.coeffSigma2')];
                            case 'Periodic'
                                cov(i) = 'covPeriodic';
                                hyp.cov = [hyp.cov; log([sample_hyp.cf{1, 1}.cf{1, i}.lengthScale'; sample_hyp.cf{1, 1}.cf{1, i}.period'; sqrt(sample_hyp.cf{1, 1}.cf{1, i}.magnSigma2')])];
                            otherwise
                                error('The type of covariance funciton is invalid for composition!')
                        end
                    end
                    covfunc{end + 1} = cov;
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
                    hyp.cov = log([sample_hyp.cf{1}.lengthScale, sqrt(sample_hyp.cf{1}.magnSigma2); sample_hyp.cf{1}.alpha])';
                elseif strcmp(model_cov{1}, 'Exp')
                    covfunc = {@covMaterniso, 1};
                    hyp.cov = log([sample_hyp.cf{1}.lengthScale, sqrt(sample_hyp.cf{1}.magnSigma2)])';
                else
                    error('The type of kernel function for MLE solution is invalid')
                end
                for i = 1 : size(data_test_rep.para, 2)
                    K = feval(covfunc{:}, hyp.cov(:,i), x_test);
                    K = (K + K') / 2;
                    K = K + sample_hyp.lik.sigma2(i) * eye(size(x_test, 1));
                    cri_data.para(i) = data_test' / K * data_test;
                    cri_rep.para(i) = data_test_rep.para(:, i)' / K * data_test_rep.para(:, i);
                end
            end
            if isfield(data_test_rep, 'para_obs')
                for i = 1 : size(data_test_rep.para_obs, 2)
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
                    [Eft, Varft, Lpyt, Eyt, Varyt, Covyt] = gp_pred(gp, x_train, data_train, x_test);
                    cri_data.para_obs(i) = (data_test - Eyt)' / Covyt * (data_test - Eyt);
                    cri_rep.para_obs(i) = (data_test_rep.para_obs(:, i) - Eyt)' / Covyt * (data_test_rep.para_obs(:, i) - Eyt);
                end
            end
%             cri_data = zeros(size(data_test_rep, 2), 1);
%             cri_rep = zeros(size(data_test_rep, 2), 1);
%             if strcmp(sample_type, 'para')
%                 if strcmp(model_cov{1}, 'Periodic')
%                     covfunc = {@covPeriodic};
%                     hyp.cov = log([sample_hyp.cf{1}.lengthScale, sample_hyp.cf{1}.period, sqrt(sample_hyp.cf{1}.magnSigma2)])';
%                 elseif strcmp(model_cov{1}, 'SE')
%                     covfunc = {@covSEiso};
%                     hyp.cov = log([sample_hyp.cf{1}.lengthScale, sqrt(sample_hyp.cf{1}.magnSigma2)])';
%                 elseif strcmp(model_cov{1}, 'Matern')
%                     covfunc = {@covMaterniso, model_cov{2}};
%                     hyp.cov = log([sample_hyp.cf{1}.lengthScale, sqrt(sample_hyp.cf{1}.magnSigma2)])';
%                 elseif strcmp(model_cov{1}, 'RQ')
%                     covfunc = {@covRQiso};
%                     hyp.cov = log([sample_hyp.cf{1}.lengthScale, sqrt(sample_hyp.cf{1}.magnSigma2), sample_hyp.cf{1}.alpha])';
%                 elseif strcmp(model_cov{1}, 'Exp')
%                     covfunc = {@covMaterniso, 1};
%                     hyp.cov = log([sample_hyp.cf{1}.lengthScale, sqrt(sample_hyp.cf{1}.magnSigma2)])';
%                 else
%                     error('The type of covariance function is invalid')
%                 end
%                 for i = 1 : size(data_test_rep, 2)
%                     K = feval(covfunc{:}, hyp.cov(:,i), x_test);
%                     K = (K + K') / 2;
%                     K = K + sample_hyp.lik.sigma2(i) * eye(size(x_test, 1));
%                     cri_data(i) = data_test' / K * data_test;
%                     cri_rep(i) = data_test_rep(:, i)' / K * data_test_rep(:, i);
%                 end
%             elseif strcmp(sample_type, 'para&obs')
%                 for i = 1 : size(data_test_rep, 2)
%                     if strcmp(model_cov{1}, 'Periodic')
%                         gpcf = gpcf_periodic('lengthScale', sample_hyp.cf{1}.lengthScale(i), 'magnSigma2', sample_hyp.cf{1}.magnSigma2(i), 'period', sample_hyp.cf{1}.period(i), 'decay', 0);
%                     elseif strcmp(model_cov{1}, 'SE')
%                         gpcf = gpcf_sexp('lengthScale', sample_hyp.cf{1}.lengthScale(i), 'magnSigma2', sample_hyp.cf{1}.magnSigma2(i));
%                     elseif strcmp(model_cov{1}, 'Matern')
%                         if strcmp(model_cov{2}, '3')
%                             gpcf = gpcf_matern32('lengthScale', sample_hyp.cf{1}.lengthScale(i), 'magnSigma2', sample_hyp.cf{1}.magnSigma2(i));
%                         elseif strcmp(model_cov{2}, '5')
%                             gpcf = gpcf_matern52('lengthScale', sample_hyp.cf{1}.lengthScale(i), 'magnSigma2', sample_hyp.cf{1}.magnSigma2(i));
%                         else
%                             error('Only support Matern 3/2 and Matern 5/2 covariance function!')
%                         end
%                     elseif strcmp(model_cov{1}, 'RQ')
%                         gpcf = gpcf_rq('lengthScale', sample_hyp.cf{1}.lengthScale(i), 'magnSigma2', sample_hyp.cf{1}.magnSigma2(i), 'alpha', sample_hyp.cf{1}.alpha(i));
%                     else
%                         error('The type of covariance function is invalid')
%                     end
%                     lik = lik_gaussian('sigma2', sample_hyp.lik.sigma2(i));
%                     gp = gp_set('lik', lik, 'cf', gpcf);
%                     [Eft, Varft, Lpyt, Eyt, Varyt, Covyt] = gp_pred(gp, x_train, data_train, x_test);
%                     cri_data(i) = (data_test - Eyt)' / Covyt * (data_test - Eyt);
%                     cri_rep(i) = (data_test_rep(:, i) - Eyt)' / Covyt * (data_test_rep(:, i) - Eyt);
%                 end
%             else
%                 error('The type of sampling is invalid!')
%             end
        elseif strcmp(criteria, 'frequency')
            fft_data_test = fft(data_test);
            tem = abs(fft_data_test / size(data_test, 1));
            fre_data_test = tem(1 : size(data_test, 1) / 2 + 1);
            fre_data_test(2 : end - 1) = 2 * fre_data_test(2 : end - 1);
            f_data_test = 1 / (x_test(2) - x_test(1)) * (0 : size(data_test, 1) / 2) / size(data_test, 1);
        else
            error('The criteria is invalid!')
        end
        if isfield(data_test_rep, 'para')
            p_value.para = 2 * min(sum(cri_rep.para > cri_data.para) / size(data_test_rep.para, 2), 1 - sum(cri_rep.para > cri_data.para) / size(data_test_rep.para, 2));
            disp(['p-value for p(y^rep|para) is: ', num2str(p_value.para)])
        end
        if isfield(data_test_rep, 'para_obs')
            p_value.para_obs = 2 * min(sum(cri_rep.para_obs > cri_data.para_obs) / size(data_test_rep.para_obs, 2), 1 - sum(cri_rep.para_obs > cri_data.para_obs) / size(data_test_rep.para_obs, 2));
            disp(['p-value for p(y^rep|para, y^obs) is: ', num2str(p_value.para_obs)])
        end
    end
end