function [prior] = mle(hyp_mle, x, data, model_cov)
    % choose the kernel function (mle_covfunc) for MLE solution from:
    % Matern kernel (Matern) (only accept Matern 3/2 (d = 1) and Matern 5/2 (d = 2))
    % SE kernel (SE)
    % Periodic kernel (Periodic)
    % Rational quadratic (RQ)

    % hyp_mle: the initialisation of hyperparameters
    % x: the training input location
    % data: the training data

    likfunc = @likGauss;
    hyp.lik = log(hyp_mle.sigma2);
    hyp.cov = [];
    prior = {};
    if strcmp(model_cov{1}, 'sum') || strcmp(model_cov{1}, 'prod')
        if strcmp(model_cov{1}, 'sum')
            covfunc = {'covSum', {}};
        else
            covfunc = {'covProd', {}};
        end
        for i = 1 : size(model_cov, 2) - 1
            switch model_cov{i + 1}
                case 'SE'
                    covfunc{2} = [covfunc{2}, 'covSEiso'];
                    hyp.cov = [hyp.cov; log([hyp_mle.lengthScale; sqrt(hyp_mle.magnSigma2)])];
                case 'LIN'
                    covfunc{2} = [covfunc{2}, 'covLINiso'];
                    hyp.cov = [hyp.cov; log(hyp_mle.lengthScale)];
                case 'Periodic'
                    covfunc{2} = [covfunc{2}, 'covPeriodic'];
                    hyp.cov = [hyp.cov; log([hyp_mle.lengthScale; hyp_mle.period; sqrt(hyp_mle.magnSigma2)])];
                otherwise
                    error('The type of covariance funciton is invalid for composition!')
            end
        end
        hyp = minimize(hyp, @gp, -150, @infGaussLik, [], covfunc, likfunc, x, data);
        ord = 1;
        for i = 1 : size(model_cov, 2) - 1
            switch model_cov{i + 1}
                case 'SE'
                    hyp_mle.SE.lengthScale = exp(hyp.cov(ord));
                    hyp_mle.SE.magnSigma2 = exp(hyp.cov(ord + 1))^2;
                    ord = ord + 2;
                case 'LIN'
                    hyp_mle.LIN.lengthScale = exp(hyp.cov(ord));
                    ord = ord + 1;
                case 'Periodic'
                    hyp_mle.Periodic.lengthScale = exp(hyp.cov(ord));
                    hyp_mle.Periodic.period = exp(hyp.cov(ord + 1));
                    hyp_mle.Periodic.magnSigma2 = exp(hyp.cov(ord + 2))^2;
                    ord = ord + 3;
                otherwise
                    error('The type of covariance funciton is invalid for composition!')
            end
        end
    else
        if strcmp(model_cov{1}, 'Periodic')
            covfunc = {@covPeriodic};
            hyp.cov = log([hyp_mle.lengthScale; hyp_mle.period; sqrt(hyp_mle.magnSigma2)]);
            hyp = minimize(hyp, @gp, -100, @infGaussLik, [], covfunc, likfunc, x, data);
            hyp_mle.Periodic.lengthScale = exp(hyp.cov(1));
            hyp_mle.Periodic.period = exp(hyp.cov(2));
            hyp_mle.Periodic.magnSigma2 = exp(hyp.cov(3))^2;
        elseif strcmp(model_cov{1}, 'SE')
            covfunc = {@covSEiso};
            hyp.cov = log([hyp_mle.lengthScale; sqrt(hyp_mle.magnSigma2)]);
            hyp = minimize(hyp, @gp, -100, @infGaussLik, [], covfunc, likfunc, x, data);
            hyp_mle.SE.lengthScale = exp(hyp.cov(1));
            hyp_mle.SE.magnSigma2 = exp(hyp.cov(2))^2;
        elseif strcmp(model_cov{1}, 'Matern')
            covfunc = {@covMaterniso, str2double(model_cov{2})};
            hyp.cov = log([hyp_mle.lengthScale; sqrt(hyp_mle.magnSigma2)]);
            hyp = minimize(hyp, @gp, -100, @infGaussLik, [], covfunc, likfunc, x, data);
            hyp_mle.Matern.lengthScale = exp(hyp.cov(1));
            hyp_mle.Matern.magnSigma2 = exp(hyp.cov(2))^2;
        elseif strcmp(model_cov{1}, 'RQ')
            covfunc = {@covRQiso};
            hyp.cov = log([hyp_mle.lengthScale; sqrt(hyp_mle.magnSigma2); hyp_mle.alpha]);
            hyp = minimize(hyp, @gp, -100, @infGaussLik, [], covfunc, likfunc, x, data);
            hyp_mle.RQ.lengthScale = exp(hyp.cov(1));
            hyp_mle.RQ.magnSigma2 = exp(hyp.cov(2))^2;
            hyp_mle.RQ.alpha = exp(hyp.cov(3));
        elseif strcmp(model_cov{1}, 'Exp')
            covfunc = {@covMaterniso, 1};
            hyp.cov = log([hyp_mle.lengthScale; sqrt(hyp_mle.magnSigma2)]);
            hyp = minimize(hyp, @gp, -100, @infGaussLik, [], covfunc, likfunc, x, data);
            hyp_mle.Exp.lengthScale = exp(hyp.cov(1));
            hyp_mle.Exp.magnSigma2 = exp(hyp.cov(2))^2;
        else
            error('The type of kernel function for MLE solution is invalid')
        end
    end
    %% set prior by using MLE
    hyp_mle.sigma2 = exp(hyp.lik);
    prior.sigma2 = prior_gaussian('s2', hyp_mle.sigma2, 'mu', hyp_mle.sigma2);
    if isfield(hyp_mle, 'SE')
        prior.SE.lengthScale = prior_gaussian('s2', hyp_mle.SE.lengthScale, 'mu', hyp_mle.SE.lengthScale);
        prior.SE.magnSigma2 = prior_gaussian('s2', hyp_mle.SE.magnSigma2, 'mu', hyp_mle.SE.magnSigma2);
    end
    if isfield(hyp_mle, 'LIN')
        prior.LIN.lengthScale = prior_gaussian('s2', hyp_mle.LIN.lengthScale, 'mu', hyp_mle.LIN.lengthScale);
    end
    if isfield(hyp_mle, 'Periodic')
        prior.Periodic.lengthScale = prior_gaussian('s2', hyp_mle.Periodic.lengthScale, 'mu', hyp_mle.Periodic.lengthScale);
        prior.Periodic.period = prior_gaussian('s2', hyp_mle.Periodic.period, 'mu', hyp_mle.Periodic.period);
        prior.Periodic.magnSigma2 = prior_gaussian('s2', hyp_mle.Periodic.magnSigma2, 'mu', hyp_mle.Periodic.magnSigma2);
    end
    if isfield(hyp_mle, 'Matern')
        prior.Matern.lengthScale = prior_gaussian('s2', hyp_mle.Matern.lengthScale, 'mu', hyp_mle.Matern.lengthScale);
        prior.Matern.magnSigma2 = prior_gaussian('s2', hyp_mle.Matern.magnSigma2, 'mu', hyp_mle.Matern.magnSigma2);
    end
    if isfield(hyp_mle, 'RQ')
        prior.RQ.lengthScale = prior_gaussian('s2', hyp_mle.RQ.lengthScale, 'mu', hyp_mle.RQ.lengthScale);
        prior.RQ.magnSigma2 = prior_gaussian('s2', hyp_mle.RQ.magnSigma2, 'mu', hyp_mle.RQ.magnSigma2);
        prior.RQ.alpha = prior_gaussian('s2', hyp_mle.RQ.alpha, 'mu', hyp_mle.RQ.alpha);
    end
    if isfield(hyp_mle, 'Exp')
        prior.Exp.lengthScale = prior_gaussian('s2', hyp_mle.Exp.lengthScale, 'mu', hyp_mle.Exp.lengthScale);
        prior.Exp.magnSigma2 = prior_gaussian('s2', hyp_mle.Exp.magnSigma2, 'mu', hyp_mle.Exp.magnSigma2);
    end
end