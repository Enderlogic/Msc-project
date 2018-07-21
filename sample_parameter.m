function [sample_hyp] = sample_parameter(num_rep, x_train, data_train, prior, model_cov)
    % num_rep: the number of replicated data
    % model_cov: the covariance function used for fitted model
    % x_train: the training location
    % data_train: the training data
    % prior: the prior distribution for hyperparameters

    % Construct and Initialise GP model
    lik = lik_gaussian('sigma2', prior.sigma2.mu, 'sigma2_prior', prior.sigma2);
    if strcmp(model_cov{1}, 'sum') || strcmp(model_cov{1}, 'prod')
        gpcf_com = {};
        for i = 1 : size(model_cov, 2) - 1
            switch model_cov{i + 1}
                case 'SE'
                    gpcf_com{end + 1} = gpcf_sexp('lengthScale', prior.SE.lengthScale.mu, 'lengthScale_prior', prior.SE.lengthScale, ...
                        'magnSigma2', prior.SE.magnSigma2.mu, 'magnSigma2_prior', prior.SE.magnSigma2);
                case 'LIN'
                    gpcf_com{end + 1} = gpcf_linear('coeffSigma2', prior.LIN.lengthScale.mu, 'coeffSigma2_prior', prior.LIN.lengthScale);
                case 'Periodic'
                    gpcf_com{end + 1} = gpcf_periodic('lengthScale', prior.Periodic.lengthScale.mu,...
                        'lengthScale_prior', prior.Periodic.lengthScale, ...
                        'magnSigma2', prior.Periodic.magnSigma2.mu, 'magnSigma2_prior', prior.Periodic.magnSigma2, ...
                        'decay', 0, 'period', prior.Periodic.period.mu, 'period_prior', prior.Periodic.period);
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
        gpcf = gpcf_periodic('lengthScale', prior.Periodic.lengthScale.mu, 'lengthScale_prior', prior.Periodic.lengthScale, ...
            'magnSigma2', prior.Periodic.magnSigma2.mu, 'magnSigma2_prior', prior.Periodic.magnSigma2, ...
            'decay', 0, 'period', prior.Periodic.period.mu, 'period_prior', prior.Periodic.period);
    elseif strcmp(model_cov{1}, 'SE')
        gpcf = gpcf_sexp('lengthScale', prior.SE.lengthScale.mu, 'lengthScale_prior', prior.SE.lengthScale, ...
            'magnSigma2', prior.SE.magnSigma2.mu, 'magnSigma2_prior', prior.SE.magnSigma2);
    elseif strcmp(model_cov{1}, 'Matern')
        if strcmp(model_cov{2}, '3')
            gpcf = gpcf_matern32('lengthScale', prior.Matern.lengthScale.mu, 'lengthScale_prior', prior.Matern.lengthScale, ...
                'magnSigma2', prior.Matern.magnSigma2.mu, 'magnSigma2_prior', prior.Matern.magnSigma2);
        elseif strcmp(model_cov{2}, '5')
            gpcf = gpcf_matern52('lengthScale', prior.Matern.lengthScale.mu, 'lengthScale_prior', prior.Matern.lengthScale, ...
                'magnSigma2', prior.Matern.magnSigma2.mu, 'magnSigma2_prior', prior.Matern.magnSigma2);
        else
            error('Only support Matern 3/2 and Matern 5/2 covariance function!')
        end
    elseif strcmp(model_cov{1}, 'RQ')
        gpcf = gpcf_rq('lengthScale', prior.RQ.lengthScale.mu, 'lengthScale_prior', prior.RQ.lengthScale, ...
            'magnSigma2', prior.RQ.magnSigma2.mu, 'magnSigma2_prior', prior.RQ.magnSigma2, ...
            'alpha', prior.RQ.alpha.mu, 'alpha_prior', prior.RQ.alpha);
    elseif strcmp(model_cov{1}, 'Exp')
        gpcf = gpcf_exp('lengthScale', prior.Exp.lengthScale.mu, 'lengthScale_prior', prior.Exp.lengthScale, ...
            'magnSigma2', prior.Exp.magnSigma2.mu, 'magnSigma2_prior', prior.Exp.magnSigma2);
    else
        error('The type of kernel function for Sampling is invalid')
    end
    gp = gp_set('lik', lik, 'cf', gpcf);

    % Update hyperparameters
    opt = optimset('TolX',1e-4,'TolFun',1e-4);
    gp = gp_optim(gp, x_train, data_train, 'opt', opt);

    % Monte Carlo sampling
    num_thin = round(0.2 * num_rep);
%     % Set some options for HMC
%     hmc_opt.nsamples=1;
%     hmc_opt.decay=0.8;
%     hmc_opt.persistence=0;
%     hmc_opt.stepadj=0.04;
%     hmc_opt.steps=10;
    [sample_hyp,g,opt]=gp_mc(gp, x_train, data_train, 'nsamples', num_rep + num_thin - 1, 'display',num_thin);
    sample_hyp=thin(sample_hyp, num_thin);
end