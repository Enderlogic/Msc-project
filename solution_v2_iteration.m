clear

option_data = 'syn'; % choose synthetic data (syn) or real data (real)
data_baseset = {{'SE'}, {'Periodic'}, {'RQ'}, {'Exp'}, {'Matern', '3'}, {'Matern', '5'}};
model_baseset = {{'SE'}, {'Periodic'}, {'RQ'}, {'Exp'}, {'Matern', '3'}, {'Matern', '5'}};
ratio_train = 0.6; % set the ratio of training set among the whole set.
num_rep = 5; % set the amount of samples
sample_location = 'test'; % choose the location of replications: 'test' for test set, 'train' for training set
sample_type = 'both'; % choose the type of samples: 'para': p(y^rep|para)
                        % 'para&obs'; p(y^rep|para, y^obs)
                        % 'both': draw from both of the two distributions
save = true; % true: save the parameters and data; false: not save

for i = 1 : length(data_baseset)
    for j = 1 : length(model_baseset)
        data_name = data_baseset{i};
        model_cov = model_baseset{i};
        filename = strcat('model\', option_data, '_', strjoin(data_name), '_', strjoin(model_cov), '_', sample_location, '_', sample_type, '.mat');
        if exist(filename, 'file') == 2
            load(filename)
        else
            [data_test, data_rep, sample_hyp, x_train, data_train, x_test] = model_generation(option_data, data_name, ratio_train, model_cov, sample_location, sample_type, num_rep, save);
        end
    end
end