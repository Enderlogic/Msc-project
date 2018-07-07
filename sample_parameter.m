function [sample_hyp] = sample_parameter(num_rep, x_train, data_train, prior, model_cov)
% num_rep: the number of replicated data
% model_cov: the covariance function used for fitted model
% x_train: the training location
% data_train: the training data
% prior: the prior distribution for hyperparameters

% Construct and Initialise GP model
lik = lik_gaussian('sigma2', prior.sigma2.mu, 'sigma2_prior', prior.sigma2);
if strcmp(model_cov{1}, 'Periodic')
%     gpcf = gpcf_periodic();
    gpcf = gpcf_periodic('lengthScale', prior.lengthScale.mu, 'lengthScale_prior', prior.lengthScale, ...
        'magnSigma2', prior.magnSigma2.mu, 'magnSigma2_prior', prior.magnSigma2, ...
        'decay', 0, 'period', prior.period.mu, 'period_prior', prior.period);
elseif strcmp(model_cov{1}, 'SE')
    gpcf = gpcf_sexp('lengthScale', prior.lengthScale.mu, 'lengthScale_prior', prior.lengthScale, ...
        'magnSigma2', prior.magnSigma2.mu, 'magnSigma2_prior', prior.magnSigma2);
elseif strcmp(model_cov{1}, 'Matern')
    if model_cov{2} == 1
        gpcf = gpcf_matern32('lengthScale', prior.lengthScale.mu, 'lengthScale_prior', prior.lengthScale, ...
            'magnSigma2', prior.magnSigma2.mu, 'magnSigma2_prior', prior.magnSigma2);
%         gpcf = gpcf_matern32();
    elseif model_cov{2} == 2
        gpcf = gpcf_matern52('lengthScale', prior.lengthScale.mu, 'lengthScale_prior', prior.lengthScale, ...
            'magnSigma2', prior.magnSigma2.mu, 'magnSigma2_prior', prior.magnSigma2);
    else
        error('Only support Matern 3/2 and Matern 5/2 covariance function!')
    end
else
    error('The type of kernel function for Sampling is invalid')
end
gp = gp_set('lik', lik, 'cf', gpcf);

% Update hyperparameters
opt = optimset('TolX',1e-4,'TolFun',1e-4);
gp = gp_optim(gp, x_train, data_train, 'opt', opt);

% Monte Carlo sampling
num_thin = round(0.2 * num_rep);
[sample_hyp,g,opt]=gp_mc(gp, x_train, data_train, 'nsamples', num_rep + num_thin - 1, 'display',num_thin);
sample_hyp=thin(sample_hyp, num_thin);
end

