function [hyp_mle] = mle_gpml(hyp_mle, x, data, model_cov)
% choose the kernel function (mle_covfunc) for MLE solution from:
% Matern kernel (Matern) (only accept Matern 3/2 (d = 1) and Matern 5/2 (d = 2))
% SE kernel (SE)
% Periodic kernel (Periodic)
% Rational quadratic (RQ)

% hyp_mle: the initialisation of hyperparameters
% x: the training input location
% data: the training data

likfunc = @likGauss;
hyp.lik = log(hyp_mle.sigma);
if strcmp(model_cov{1}, 'Periodic')
    covfunc = {@covPeriodic};
    hyp.cov = log([hyp_mle.lengthScale; hyp_mle.period; sqrt(hyp_mle.magnitude)]);
    hyp = minimize(hyp, @gp, -100, @infGaussLik, [], covfunc, likfunc, x, data);
    hyp_mle.lengthScale = exp(hyp.cov(1));
    hyp_mle.period = exp(hyp.cov(2));
    hyp_mle.magnitude = exp(hyp.cov(3))^2;
elseif strcmp(model_cov{1}, 'SE')
    covfunc = {@covSEiso};
    hyp.cov = log([hyp_mle.lengthScale; sqrt(hyp_mle.magnitude)]);
    hyp = minimize(hyp, @gp, -100, @infGaussLik, [], covfunc, likfunc, x, data);
    hyp_mle.lengthScale = exp(hyp.cov(1));
    hyp_mle.magnitude = exp(hyp.cov(2))^2;
elseif strcmp(model_cov{1}, 'Matern')
    covfunc = {@covMaterniso, model_cov{2}};
    hyp.cov = log([hyp_mle.lengthScale; sqrt(hyp_mle.magnitude)]);
    hyp = minimize(hyp, @gp, -100, @infGaussLik, [], covfunc, likfunc, x, data);
    hyp_mle.lengthScale = exp(hyp.cov(1));
    hyp_mle.magnitude = exp(hyp.cov(2))^2;
elseif strcmp(model_cov{1}, 'RQ')
    covfunc = {@covRQiso};
    hyp.cov = log([hyp_mle.lengthScale; sqrt(hyp_mle.magnitude); hyp_mle.alpha]);
    hyp = minimize(hyp, @gp, -100, @infGaussLik, [], covfunc, likfunc, x, data);
    hyp_mle.lengthScale = exp(hyp.cov(1));
    hyp_mle.magnitude = exp(hyp.cov(2))^2;
    hyp_mle.alpha = exp(hyp.cov(3));
else
    error('The type of kernel function for MLE solution is invalid')
end
hyp_mle.sigma = exp(hyp.lik);
end

