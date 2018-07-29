function [data_test, data_rep, sample_hyp, x_train, data_train, x_test, hyp_data, prior] = model_generation(option_data, data_name, ratio_train, model_cov, sample_location, sample_type, num_rep, save_option)
    %% Data generation
    if strcmp(option_data, 'syn')
        % set the hyperparameters of ficticious data set (not all of them will be used)
        hyp_data.magnSigma2 = 10; % squared magnitude
        hyp_data.lengthScale = 1;
        hyp_data.period = 2;
        hyp_data.sigma2 = 0.01;
        hyp_data.alpha = 5;
        hyp_data.length = 20; % choose the length of data
        hyp_data.random_seed = true; % choose random seed for data, if seed is true, the data will be the same; if seed is false, the data will be different every time
        hyp_data.data_name = data_name; % choose a kernel for generation model ('SE'; 'Periodic'; 'Matern', '3' (nu = 3/2); 'Matern', '5' (nu = 5/2); 'RQ')
        num_data = 1000; % choose sampling number

        [x, data] = data_generation(hyp_data, num_data);
    elseif strcmp(option_data, 'real')
        datfiles = dir('data\*.mat');
        baseFileNames = natsortfiles({datfiles.name});
        fullFilename = fullfile('data\', baseFileNames{str2double(data_name)});
        Data = load(fullFilename);
        x = Data.X;
        data = Data.y;
        num_data = size(data, 1);
    else
        error('The source of observation is invalid')
    end

    num_data_train = round(num_data * ratio_train);
    data_train = data(1 : num_data_train);
    x_train = x(1 : num_data_train);
    data_test = data(num_data_train + 1 : num_data);
    x_test = x(num_data_train + 1 : num_data);
    disp('Data generation complete!')
    %% Compute MLE solution using gpml
    %initilise the hyperparameters for MLE
    hyp_mle.SE.magnSigma2 = 20; %squared magnitude
    hyp_mle.SE.lengthScale = 2;
    hyp_mle.Periodic.magnSigma2 = 100; %squared magnitude
    hyp_mle.Periodic.lengthScale = 10;
    hyp_mle.Periodic.period = 1;
    hyp_mle.Matern.magnSigma2 = 10;
    hyp_mle.Matern.lengthScale = 2;
    hyp_mle.Exp.magnSigma2 = 10;
    hyp_mle.Exp.lengthScale = 2;
    hyp_mle.LIN.lengthScale = 40;
    hyp_mle.RQ.magnSigma2 = 20; %squared magnitude
    hyp_mle.RQ.lengthScale = 4;
    hyp_mle.RQ.alpha = 2;
    hyp_mle.sigma2 = 0.01;

    prior = mle(hyp_mle, x_train, data_train, model_cov);
    disp('MLE solution complete!')
    %% Draw samples of hyperparameters from posterior distribution
    sample_hyp = sample_parameter(num_rep, x_train, data_train, prior, model_cov);
    disp('Parameter sampling complete!')
    %% Generate replicated data
    data_rep = sample_data(sample_hyp, model_cov, sample_location, sample_type, x_train, data_train, x_test);
    disp('Data sampling complete!')
    %% Save the useful information
    if save_option
        filename = strcat('model\', option_data, '_', strjoin(data_name), '_', strjoin(model_cov), '_', sample_type, '.mat');
        if strcmp(option_data, 'syn')
            save(filename, 'data_test', 'data_rep', 'hyp_data', 'data_name', 'x_train', 'x_test', 'data_train', 'model_cov', 'sample_location', 'prior', 'ratio_train', 'sample_hyp', 'sample_type')
        elseif strcmp(option_data, 'real')
            save(filename, 'data_test', 'data_rep', 'data_name', 'x_train', 'x_test', 'data_train', 'model_cov', 'sample_location', 'prior', 'ratio_train', 'sample_hyp', 'sample_type')
        else
            error('The source of observation is invalid')
        end
        disp('Model saving complete!')
    end
end