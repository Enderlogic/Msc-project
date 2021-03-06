function [prior] = mle(hyp_mle, x, data, model_cov)
    % choose the kernel function (mle_covfunc) for MLE solution from:
    % Matern kernel (Matern) (only accept Matern 3/2 (d = 1) and Matern 5/2 (d = 2))
    % SE kernel (SE)
    % Periodic kernel (Periodic)
    % Rational quadratic (RQ)

    % hyp_mle: the initialisation of hyperparameters
    % x: the training input location
    % data: the training data

%     [x, ~, ~]  = zscore(x);

    likfunc = @likGauss;
    meanfunc = @meanConst;
    hyp.mean = mean(data);
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
                    hyp.cov = [hyp.cov; log([hyp_mle.SE.lengthScale; sqrt(hyp_mle.SE.magnSigma2)])];
                case 'LIN'
                    covfunc{2} = [covfunc{2}, 'covLINiso'];
                    hyp.cov = [hyp.cov; log(hyp_mle.LIN.lengthScale)];
                case 'Periodic'
                    covfunc{2} = [covfunc{2}, 'covPeriodic'];
                    hyp.cov = [hyp.cov; log([hyp_mle.Periodic.lengthScale; hyp_mle.Periodic.period; sqrt(hyp_mle.Periodic.magnSigma2)])];
                otherwise
                    error('The type of covariance funciton is invalid for composition!')
            end
        end
        hyp = minimize(hyp, @gp, -150, @infVB, meanfunc, covfunc, likfunc, x, data);
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
            hyp.cov = log([hyp_mle.Periodic.lengthScale; hyp_mle.Periodic.period; sqrt(hyp_mle.Periodic.magnSigma2)]);
            hyp = minimize(hyp, @gp, -100, @infGaussLik, meanfunc, covfunc, likfunc, x, data);
            hyp_mle.Periodic.lengthScale = exp(hyp.cov(1));
            hyp_mle.Periodic.period = exp(hyp.cov(2));
            hyp_mle.Periodic.magnSigma2 = exp(hyp.cov(3))^2;
        elseif strcmp(model_cov{1}, 'SE')
            covfunc = {@covSEiso};
            hyp.cov = log([hyp_mle.SE.lengthScale; sqrt(hyp_mle.SE.magnSigma2)]);
            hyp = minimize(hyp, @gp, -100, @infVB, meanfunc, covfunc, likfunc, x, data);
            hyp_mle.SE.lengthScale = exp(hyp.cov(1));
            hyp_mle.SE.magnSigma2 = exp(hyp.cov(2))^2;
        elseif strcmp(model_cov{1}, 'Matern')
            covfunc = {@covMaterniso, str2double(model_cov{2})};
%             covfunc = {@covMaterniso, 5};
%             hyp.cov = log([median(pdist(x))^2; sqrt(hyp_mle.Matern.magnSigma2)]);
            hyp.cov = log([hyp_mle.Matern.lengthScale; sqrt(hyp_mle.Matern.magnSigma2)]);
            rng default
            hyp = minimize(hyp, @gp, -100, @infVB, meanfunc, covfunc, likfunc, x, data);
            hyp_mle.Matern.lengthScale = exp(hyp.cov(1));
            hyp_mle.Matern.magnSigma2 = exp(hyp.cov(2))^2;
        elseif strcmp(model_cov{1}, 'RQ')
            covfunc = {@covRQiso};
            hyp.cov = log([hyp_mle.RQ.lengthScale; sqrt(hyp_mle.RQ.magnSigma2); hyp_mle.RQ.alpha]);
            hyp = minimize(hyp, @gp, -100, @infVB, meanfunc, covfunc, likfunc, x, data);
            hyp_mle.RQ.lengthScale = exp(hyp.cov(1));
            hyp_mle.RQ.magnSigma2 = exp(hyp.cov(2))^2;
            hyp_mle.RQ.alpha = exp(hyp.cov(3));
        elseif strcmp(model_cov{1}, 'Exp')
            covfunc = {@covMaterniso, 1};
            hyp.cov = log([hyp_mle.Exp.lengthScale; sqrt(hyp_mle.Exp.magnSigma2)]);
            hyp = minimize(hyp, @gp, -100, @infGaussLik, meanfunc, covfunc, likfunc, x, data);
            hyp_mle.Exp.lengthScale = exp(hyp.cov(1));
            hyp_mle.Exp.magnSigma2 = exp(hyp.cov(2))^2;
        else
            error('The type of kernel function for MLE solution is invalid')
        end
    end
    %% set prior by using MLE
    hyp_mle.sigma2 = exp(hyp.lik);
    prior.sigma2 = prior_gamma('sh', hyp_mle.sigma2, 'is', 1);
%     prior.sigma2 = prior_loggaussian('s2', hyp_mle.sigma2, 'mu', hyp_mle.sigma2);
    if isfield(hyp_mle, 'SE')
        prior.SE.lengthScale = prior_gamma('sh', hyp_mle.SE.lengthScale, 'is', 1);
        prior.SE.magnSigma2 = prior_gamma('sh', hyp_mle.SE.magnSigma2, 'is', 1);
%         prior.SE.lengthScale = prior_loggaussian('s2', hyp_mle.SE.lengthScale, 'mu', hyp_mle.SE.lengthScale);
%         prior.SE.magnSigma2 = prior_loggaussian('s2', hyp_mle.SE.magnSigma2, 'mu', hyp_mle.SE.magnSigma2);
    end
    if isfield(hyp_mle, 'LIN')
        prior.LIN.lengthScale = prior_gamma('sh', hyp_mle.LIN.lengthScale, 'is', 1);
    end
    if isfield(hyp_mle, 'Periodic')
        prior.Periodic.lengthScale = prior_gamma('sh', hyp_mle.Periodic.lengthScale, 'is', 1);
        prior.Periodic.period = prior_gamma('sh', hyp_mle.Periodic.period, 'is', 1);
        prior.Periodic.magnSigma2 = prior_gamma('sh', hyp_mle.Periodic.magnSigma2, 'is', 1);
    end
    if isfield(hyp_mle, 'Matern')
        prior.Matern.lengthScale = prior_gamma('sh', hyp_mle.Matern.lengthScale, 'is', 1);
        prior.Matern.magnSigma2 = prior_gamma('sh', hyp_mle.Matern.magnSigma2, 'is', 1);
    end
    if isfield(hyp_mle, 'RQ')
        prior.RQ.lengthScale = prior_gamma('sh', hyp_mle.RQ.lengthScale, 'is', 1);
        prior.RQ.magnSigma2 = prior_gamma('sh', hyp_mle.RQ.magnSigma2, 'is', 1);
        prior.RQ.alpha = prior_gamma('sh', hyp_mle.RQ.alpha, 'is', 1);
    end
    if isfield(hyp_mle, 'Exp')
        prior.Exp.lengthScale = prior_gamma('sh', hyp_mle.Exp.lengthScale, 'is', 1);
        prior.Exp.magnSigma2 = prior_gamma('sh', hyp_mle.Exp.magnSigma2, 'is', 1);
    end
end