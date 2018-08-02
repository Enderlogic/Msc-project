function [x,data] = data_generation(hyp_data, num_data)
    % choose the kernel function of data from:
    % Matern kernel ('Matern', '3'; 'Matern', '5')
    % SE kernel (SE)
    % Periodic kernel (Periodic)
    % Rational quadratic kernel (RQ)
    % Exponential kernel (Exp)

    %choose fixed random data or variational random data ('on' or 'off')

    %set the length of data by 'length'
    x = (0 : hyp_data.length/(num_data - 1) : hyp_data.length)';
    if strcmp(hyp_data.data_name{1}, 'Periodic')
        covfunc_data = {@covPeriodic};
        hyp.cov = log([hyp_data.lengthScale, hyp_data.period, sqrt(hyp_data.magnSigma2)]);
    elseif strcmp(hyp_data.data_name{1}, 'SE')
        covfunc_data = {@covSEiso};
        hyp.cov = log([hyp_data.lengthScale, sqrt(hyp_data.magnSigma2)]);
    elseif strcmp(hyp_data.data_name{1}, 'Matern')
        covfunc_data = {'covMaterniso', str2double(hyp_data.data_name{2})};
        hyp.cov = log([hyp_data.lengthScale, sqrt(hyp_data.magnSigma2)]);
    elseif strcmp(hyp_data.data_name{1}, 'RQ')
        covfunc_data = {@covRQiso};
        hyp.cov = log([hyp_data.lengthScale, sqrt(hyp_data.magnSigma2), hyp_data.alpha]);
    elseif strcmp(hyp_data.data_name{1}, 'Exp')
        covfunc_data = {@covMaterniso, 1};
        hyp.cov = log([hyp_mle.lengthScale, sqrt(hyp_mle.magnSigma2)]);
    else
        error('The type of kernel function is invalid')
    end
    K = feval(covfunc_data{:}, hyp.cov, x);
    K = (K + K') / 2;
    mu = zeros(1,num_data);
    if hyp_data.random_seed
        rng(19990806)
        data = mvnrnd(mu, K);
        rng(19990806)
        data = mvnrnd(data, hyp_data.sigma2 * eye(num_data))';
    else
        data = mvnrnd(mu, K);
        data = mvnrnd(data, hyp_data.sigma2 * eye(num_data))';
    end
end