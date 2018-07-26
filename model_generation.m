function [data_test, data_rep, sample_hyp, x_train, data_train, x_test] = model_generation(option_data, data_name, ratio_train, model_cov, sample_location, sample_type, num_rep, save)
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
        hyp_data.data_cov = data_name; % choose a kernel for generation model ('SE'; 'Periodic'; 'Matern', '3' (nu = 3/2); 'Matern', '5' (nu = 5/2); 'RQ')
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
    hyp_mle.magnSigma2 = 10; %squared magnitude
    hyp_mle.lengthScale = 1;
    hyp_mle.period = 10;
    hyp_mle.sigma2 = 0.01;
    hyp_mle.alpha = 2;

    prior = mle(hyp_mle, x_train, data_train, model_cov);
    disp('MLE solution complete!')
    %% Draw samples of hyperparameters from posterior distribution
    sample_hyp = sample_parameter(num_rep, x_train, data_train, prior, model_cov);
    disp('Parameter sampling complete!')
    %% Generate replicated data
    data_rep = sample_data(sample_hyp, model_cov, sample_location, sample_type, x_train, data_train, x_test);
    disp('Data sampling complete!')
    %% Save the useful information
    if save
        if strcmp(option_data, 'syn')
            filename = strcat('model\', option_data, '_', strjoin(hyp_data.data_cov), '_', strjoin(model_cov), '_', sample_type, '.mat');
            save(filename, 'data_test', 'data_rep', 'hyp_data', 'x_train', 'x_test', 'data_train', 'model_cov', 'prior', 'ratio_train', 'sample_hyp', 'sample_type')
        elseif strcmp(option_data, 'real')
            filename = strcat('model\', option_data, '_', strjoin(model_cov), '_', sample_type, '.mat');
            save(filename, 'data_test', 'data_rep', 'x_train', 'x_test', 'data_train', 'model_cov', 'prior', 'ratio_train', 'sample_hyp', 'sample_type')
        else
            error('The source of observation is invalid')
        end
        disp('Model saving complete!')
    end
end